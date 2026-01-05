# mcpath/pymol_views.py
import os
import re
import glob
from typing import Dict, List, Optional, Tuple

_NUM_RE = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")

def _parse_chain_list(chain_str: str) -> Optional[set]:
    """
    chain_str: "" => None (all chains)
               "A" or "A,B" or "A B" => set({"A","B"})
    """
    chain_str = (chain_str or "").strip()
    if not chain_str:
        return None
    toks = [t.strip() for t in re.split(r"[,\s]+", chain_str) if t.strip()]
    return set(toks) if toks else None

def parse_ca_residues_in_order(pdb_path: str, chain_str: str = "") -> List[Dict[str, str]]:
    """
    Returns residues in the same order as CA atoms appear in the PDB:
      [{"chain":"A","resi":"123","resn":"GLY"}, ...]
    Uses CA atoms only; de-duplicates altloc.
    """
    want_chains = _parse_chain_list(chain_str)
    residues: List[Dict[str, str]] = []
    seen = set()

    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue

            chain = line[21].strip()
            if want_chains is not None and chain not in want_chains:
                continue

            resn = line[17:20].strip()
            resi = line[22:26].strip()
            icode = line[26].strip()
            resi_full = f"{resi}{icode}" if icode else resi

            key = (chain, resi_full)
            if key in seen:
                continue
            seen.add(key)

            residues.append({"chain": chain, "resi": resi_full, "resn": resn})

    return residues

def read_metric_file(metric_path: str, n_expected: Optional[int] = None) -> List[float]:
    """
    Robust reader for centrality outputs.
    Supports:
      - one value per line
      - "index value" per line
      - "index ... value" per line (takes the last number as the value)
    Returns a list aligned by index order (1..N) if indices are present, else by file order.
    """
    idx_to_val: Dict[int, float] = {}
    seq_vals: List[float] = []

    with open(metric_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            nums = _NUM_RE.findall(line)
            if not nums:
                continue

            # If looks like "idx value" (or idx ... value)
            if len(nums) >= 2:
                try:
                    idx = int(float(nums[0]))
                    val = float(nums[-1])
                    if idx >= 1:
                        idx_to_val[idx] = val
                        continue
                except Exception:
                    pass

            # Otherwise treat as single value line
            try:
                seq_vals.append(float(nums[-1]))
            except Exception:
                continue

    if idx_to_val:
        max_idx = max(idx_to_val.keys())
        N = n_expected if (n_expected and n_expected > 0) else max_idx
        out = []
        for i in range(1, N + 1):
            out.append(idx_to_val.get(i, float("nan")))
        return out

    return seq_vals

def _find_metric_file(work_dir: str, kind: str) -> Optional[str]:
    """
    Autodetect a metric file produced by mcpath.closeness / mcpath.betweenness.
    We look for common substrings and prefer newest file.
    """
    kind = kind.lower()
    if kind not in ("closeness", "betweenness"):
        raise ValueError("kind must be 'closeness' or 'betweenness'")

    patterns = []
    if kind == "closeness":
        patterns = [
            "*closeness*.txt", "*closeness*.dat", "*closeness*.out", "*closeness*",
            "*close*.txt", "*close*.dat", "*close*.out"
        ]
    else:
        patterns = [
            "*betweenness*.txt", "*betweenness*.dat", "*betweenness*.out", "*betweenness*",
            "*betw*.txt", "*betw*.dat", "*betw*.out"
        ]

    candidates = []
    for pat in patterns:
        candidates.extend(glob.glob(os.path.join(work_dir, pat)))

    # Filter obvious non-data
    candidates = [
        c for c in candidates
        if os.path.isfile(c)
        and not c.endswith((".py", ".pml", ".pse", ".png", ".jpg", ".jpeg"))
    ]
    if not candidates:
        return None

    candidates.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return candidates[0]

def _top_k(values: List[float], residues: List[Dict[str, str]], k: int) -> List[Tuple[int, float, Dict[str, str]]]:
    rows = []
    n = min(len(values), len(residues))
    for i in range(n):
        v = values[i]
        if v != v:  # NaN
            continue
        rows.append((i + 1, float(v), residues[i]))  # 1-based index
    rows.sort(key=lambda t: t[1], reverse=True)
    return rows[: max(1, int(k))]

def _write_tsv(path: str, rows: List[Tuple[int, float, Dict[str, str]]]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("ranked_index\tvalue\tchain\tresi\tresn\n")
        for idx, val, r in rows:
            f.write(f"{idx}\t{val:.6g}\t{r['chain']}\t{r['resi']}\t{r['resn']}\n")

def _pml_alter_b_factors(obj: str, values: List[float], residues: List[Dict[str, str]]) -> List[str]:
    lines = []
    n = min(len(values), len(residues))
    for i in range(n):
        v = values[i]
        if v != v:  # NaN
            continue
        ch = residues[i]["chain"]
        resi = residues[i]["resi"]
        # Set B-factor on all atoms of that residue (works well with spectrum b)
        lines.append(f"alter ({obj} and chain {ch} and resi {resi}), b={v}")
    lines.append("rebuild")
    return lines

def _pml_select_top(obj: str, sel_name: str, top_rows: List[Tuple[int, float, Dict[str, str]]]) -> List[str]:
    # Build: (chain A and resi 10) or (chain B and resi 55) ...
    parts = []
    for _, _, r in top_rows:
        parts.append(f"(chain {r['chain']} and resi {r['resi']})")
    expr = " or ".join(parts) if parts else "none"
    return [f"select {sel_name}, ({obj} and ({expr}))"]

def build_pymol_views(
    pdb_path: str,
    work_dir: str,
    closeness_path: str,
    betweenness_path: str,
    chain_id: str = "",
    top_n: int = 30,
    label_top: int = 10,
    out_prefix: Optional[str] = None,
    write_png: bool = True,
) -> str:
    """
    Writes a .pml that creates two objects (prot_close, prot_betw),
    stores 2 scenes (closeness, betweenness), highlights top residues, and (optionally) writes PNGs.
    Returns the .pml path.
    """
    os.makedirs(work_dir, exist_ok=True)

    base = out_prefix or os.path.splitext(os.path.basename(pdb_path))[0]
    pml_path = os.path.join(work_dir, f"{base}_centrality_views.pml")

    residues = parse_ca_residues_in_order(pdb_path, chain_id)
    if not residues:
        raise RuntimeError("No CA residues found (check PDB format / chain filter).")

    close_vals = read_metric_file(closeness_path, n_expected=len(residues))
    betw_vals  = read_metric_file(betweenness_path, n_expected=len(residues))

    close_top = _top_k(close_vals, residues, top_n)
    betw_top  = _top_k(betw_vals,  residues, top_n)

    # Write companion TSVs (very handy for reporting)
    _write_tsv(os.path.join(work_dir, f"{base}_top_closeness.tsv"), close_top)
    _write_tsv(os.path.join(work_dir, f"{base}_top_betweenness.tsv"), betw_top)

    # PML
    obj0 = "prot"
    objC = "prot_close"
    objB = "prot_betw"

    lines = []
    lines += [
        "reinitialize",
        f"load {os.path.abspath(pdb_path)}, {obj0}",
        "remove solvent",
        "hide everything, all",
        f"show cartoon, {obj0}",
        f"color gray80, {obj0}",
        "set cartoon_transparency, 0.15",
        "set ray_opaque_background, off",
        "bg_color white",
        "",
        f"create {objC}, {obj0}",
        f"create {objB}, {obj0}",
        f"disable {obj0}",
        "",
        "# ---- apply B-factors from centrality ----",
        f"disable {objB}",
        f"enable {objC}",
    ]
    lines += _pml_alter_b_factors(objC, close_vals, residues)
    lines += [
        f"spectrum b, rainbow, {objC}",
        "",
        f"disable {objC}",
        f"enable {objB}",
    ]
    lines += _pml_alter_b_factors(objB, betw_vals, residues)
    lines += [
        f"spectrum b, rainbow, {objB}",
        "",
        "# ---- top residue selections ----",
    ]
    lines += _pml_select_top(objC, "close_top", close_top)
    lines += _pml_select_top(objB, "betw_top", betw_top)

    # Styling: closeness view (cartoon + spheres), betweenness view (cartoon + surface + sticks)
    lines += [
        "",
        "# ---- Closeness styling ----",
        f"disable {objB}",
        f"enable {objC}",
        "show spheres, close_top and name CA",
        "set sphere_scale, 0.6, close_top",
        "show sticks, close_top",
        "set stick_radius, 0.2, close_top",
        "set sphere_quality, 2",
    ]

    # Labels (top few only)
    if label_top > 0 and close_top:
        top_for_labels = close_top[: min(label_top, len(close_top))]
        lines.append("hide labels, all")
        for _, _, r in top_for_labels:
            lines.append(
                f"label ({objC} and chain {r['chain']} and resi {r['resi']} and name CA), "
                f"\"{r['resn']}{r['resi']}:{r['chain']}\""
            )

    lines += [
        "orient",
        "zoom",
        "scene closeness, store",
    ]

    if write_png:
        lines += [
            "scene closeness, recall",
            f"png {os.path.join(os.path.abspath(work_dir), base + '_closeness.png')}, dpi=300",
        ]

    lines += [
        "",
        "# ---- Betweenness styling ----",
        f"disable {objC}",
        f"enable {objB}",
        "hide spheres, all",
        "hide sticks, all",
        "show sticks, betw_top",
        "set stick_radius, 0.25, betw_top",
        "show surface, " + objB,
        "set transparency, 0.35, " + objB,
        "orient",
        "zoom",
        "scene betweenness, store",
    ]

    if label_top > 0 and betw_top:
        top_for_labels = betw_top[: min(label_top, len(betw_top))]
        lines.append("hide labels, all")
        for _, _, r in top_for_labels:
            lines.append(
                f"label ({objB} and chain {r['chain']} and resi {r['resi']} and name CA), "
                f"\"{r['resn']}{r['resi']}:{r['chain']}\""
            )

    if write_png:
        lines += [
            "scene betweenness, recall",
            f"png {os.path.join(os.path.abspath(work_dir), base + '_betweenness.png')}, dpi=300",
        ]

    lines += [
        "",
        "# Tip: after running this script, you can save a session manually:",
        f"#   save {base}_centrality_views.pse",
    ]

    with open(pml_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    return pml_path

def build_views_autodetect(
    pdb_path: str,
    work_dir: str,
    chain_id: str = "",
    top_n: int = 30,
    label_top: int = 10,
    write_png: bool = True,
) -> str:
    """
    Convenience wrapper:
      - autodetect closeness & betweenness output files in work_dir
      - then write the PyMOL PML
    """
    close_path = _find_metric_file(work_dir, "closeness")
    betw_path  = _find_metric_file(work_dir, "betweenness")

    if not close_path:
        raise FileNotFoundError(f"Could not autodetect a closeness output file in: {work_dir}")
    if not betw_path:
        raise FileNotFoundError(f"Could not autodetect a betweenness output file in: {work_dir}")

    return build_pymol_views(
        pdb_path=pdb_path,
        work_dir=work_dir,
        closeness_path=close_path,
        betweenness_path=betw_path,
        chain_id=chain_id,
        top_n=top_n,
        label_top=label_top,
        out_prefix=os.path.splitext(os.path.basename(pdb_path))[0],
        write_png=write_png,
    )
