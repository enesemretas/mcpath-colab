# mcpath/ui.py
import os, re, requests, yaml, importlib, shutil, sys, glob, subprocess, zipfile, math
from IPython.display import display, clear_output, HTML
import ipywidgets as W
import numpy as np


# ---------- singletons (so New Job/Submit always uses the same log panel) ----------
_FORM_ROOT = None
_LOG_OUT   = None

# -------------------- Colab widgets enabling --------------------
def _enable_colab_widgets():
    """Enable widget rendering in Colab (safe no-op elsewhere)."""
    try:
        from google.colab import output
        output.enable_custom_widget_manager()
    except Exception:
        pass

# -------------------- validators & helpers --------------------

def _get_upload_bytes_and_name(pdb_upload: W.FileUpload):
    """
    Works with both:
      - ipywidgets<=7  : pdb_upload.value is dict {filename: {"content":..., "metadata": {...}}}
      - ipywidgets>=8  : pdb_upload.value is tuple/list of dicts or UploadedFile objects
    Returns (bytes, filename) or (None, None)
    """
    v = getattr(pdb_upload, "value", None)
    if not v:
        return None, None

    # ipywidgets 7 style: dict
    if isinstance(v, dict):
        fname, meta = next(iter(v.items()))
        content = meta.get("content", None)
        name = (meta.get("metadata", {}) or {}).get("name", None) or fname
        return content, name

    # ipywidgets 8 style: list/tuple
    if isinstance(v, (list, tuple)) and len(v) > 0:
        item = v[0]

        # Sometimes it's a dict
        if isinstance(item, dict):
            content = item.get("content", None)
            name = item.get("name", None) or item.get("metadata", {}).get("name", None) or "upload.pdb"
            return content, name

        # Sometimes it's an UploadedFile-like object
        content = getattr(item, "content", None)
        name = getattr(item, "name", None) or "upload.pdb"
        return content, name

    return None, None



def _is_valid_pdb_code(c):
    return bool(re.fullmatch(r"[0-9A-Za-z]{4}", (c or "").strip()))

def _is_valid_chain(ch):
    """Single-chain validator (for init/final chains etc.)."""
    return bool(re.fullmatch(r"[A-Za-z0-9]", (ch or "").strip()))

def _is_valid_chain_list(ch_str: str) -> bool:
    """
    Validate a comma-separated list of chain IDs for the global Chain ID field.

    Rules:
      - Empty string -> valid (means: all chains).
      - Otherwise: comma/whitespace-separated tokens, each must be a single [A-Za-z0-9].
    """
    ch_str = (ch_str or "").strip()
    if not ch_str:
        return True
    tokens = [t.strip() for t in re.split(r"[,\s]+", ch_str) if t.strip()]
    if not tokens:
        return True
    return all(_is_valid_chain(tok) for tok in tokens)

def _is_valid_email(s):
    return (not (s or "").strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", (s or "").strip()))

def _fetch_rcsb(code: str) -> bytes:
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.content

def _chains_from_pdb_bytes(pdb_bytes: bytes):
    """
    Return sorted unique chain IDs seen in ATOM/HETATM records.
    """
    chains = []
    seen = set()
    try:
        txt = pdb_bytes.decode("utf-8", errors="ignore")
    except Exception:
        return []

    for line in txt.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        if len(line) < 22:
            continue
        ch = line[21].strip() or " "  # blank chain possible
        if ch not in seen:
            seen.add(ch)
            chains.append(ch)
    # keep original order (more useful than alphabetical)
    return chains


def _logo_widget(branding: dict):
    url = (branding.get("logo_url") or "").strip()
    if not url:
        return None
    height = int(branding.get("logo_height", 96))
    align  = (branding.get("logo_align") or "center").lower()
    try:
        img_bytes = requests.get(url, timeout=20).content
        jc = {"center": "center", "left": "flex-start", "right": "flex-end"}.get(align, "center")
        return W.HBox([W.Image(value=img_bytes, format="png", layout=W.Layout(height=f"{height}px"))],
                      layout=W.Layout(justify_content=jc))
    except Exception:
        return None

def _progress(step: int, total: int, message: str):
    print(f"[{step}/{total}] {message}")

def _show_download(path: str, label: str = "file"):
    # user requested: no download links in output
    return

def _find_first_existing(work_dir: str, candidates):
    for name in candidates:
        p = os.path.join(work_dir, name)
        if os.path.isfile(p):
            return p
    return None

# -------------------- Peak parsing + mapping --------------------
_NUM_RE = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")

# Used only for log-style strings like "16A -> ..."
# For tabular files with separate columns, we don't rely on this for chain mapping
_RESCHAIN_RE = re.compile(r"\b(\d+)\s*([A-Za-z0-9]+)\b")

def _ensure_py3dmol():
    try:
        import py3Dmol  # noqa: F401
        return True
    except Exception:
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "py3Dmol"])
            import py3Dmol  # noqa: F401
            return True
        except Exception as e:
            print(f"Warning: py3Dmol not available and install failed: {e}")
            return False

def _parse_chain_tokens(chain_str: str):
    chain_str = (chain_str or "").strip()
    if not chain_str:
        return None
    toks = [t.strip() for t in re.split(r"[,\s]+", chain_str) if t.strip()]
    return set(toks) if toks else None

def _ca_residues_in_order(pdb_path: str, chain_str: str = ""):
    """
    CA residues in PDB order:
      [{"chain":"A","resi_raw":"123A","resi_num":123}, ...]
    """
    want = _parse_chain_tokens(chain_str)
    residues = []
    seen = set()
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            ch = line[21].strip()
            if want is not None and ch not in want:
                continue

            resi = line[22:26].strip()
            icode = line[26].strip()
            resi_raw = f"{resi}{icode}" if icode else resi

            key = (ch, resi_raw)
            if key in seen:
                continue
            seen.add(key)

            m = re.match(r"\d+", resi_raw)
            resi_num = int(m.group(0)) if m else None

            residues.append({"chain": ch, "resi_raw": resi_raw, "resi_num": resi_num})
    return residues

def _read_pure_ints(path: str):
    """Extract ONLY pure integers (avoids turning floats into 0)."""
    out = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            for m in re.finditer(r"\b\d+\b", line):
                out.append(int(m.group(0)))
    return out

def _parse_reschain_pairs(path: str):
    """
    Extract (chainToken, number) pairs from lines like:
      16A  -> ...
      1913Q -> ...
    """
    pairs = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = _RESCHAIN_RE.search(line)
            if not m:
                continue
            rn = int(m.group(1))
            ch = m.group(2)
            pairs.append((ch, rn))
    return pairs

def _parse_chain_resnum_value_lines(path: str):
    """
    Parse tabular peaks lines like:
      16   1   16.1   0.349624
    Returns list of (value, chainToken, resnum).
    NOTE: chainToken in your files is numeric ("1"), NOT 'A' => we normalize later.
    """
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if (not line) or line.startswith("#"):
                continue

            # Split by whitespace (handles tabs too)
            parts = re.split(r"\s+", line)
            if len(parts) < 2:
                continue

            # Expect: resID  chainID_numeric  encoded_node  value
            # We'll attempt robustly:
            try:
                resnum = int(float(parts[0]))
            except Exception:
                continue

            chain_tok = str(parts[1]).strip()

            # last numeric in line is peak value
            nums = _NUM_RE.findall(line)
            if not nums:
                continue
            try:
                val = float(nums[-1])
            except Exception:
                continue

            rows.append((val, chain_tok, resnum))
    return rows

def _build_chain_maps(residues_in_order):
    """
    Build per-chain ordered lists and residue-number sets for auto-detection,
    AND preserve chain appearance order (needed for numeric-chain mapping).
    """
    chain_to_ordered = {}
    chain_to_resnums = {}
    chain_order = []

    for r in residues_in_order:
        ch = r["chain"]
        if ch not in chain_to_ordered:
            chain_to_ordered[ch] = []
            chain_order.append(ch)  # first-seen order
        chain_to_ordered[ch].append((r["resi_raw"], r["resi_num"]))
        if r["resi_num"] is not None:
            chain_to_resnums.setdefault(ch, set()).add(r["resi_num"])

    return chain_to_ordered, chain_to_resnums, chain_order

def _normalize_chain_token(ch, chain_to_ordered, chain_order):
    """
    Convert numeric chain tokens ("1","2",...) into actual PDB chain IDs ("A","B",...)
    based on chain order in residues_in_order.

    If the PDB really has a chain named "1", keep it as-is.
    """
    ch = str(ch).strip()
    if ch in chain_to_ordered:
        return ch
    if ch.isdigit():
        idx = int(ch)
        if 1 <= idx <= len(chain_order):
            return chain_order[idx - 1]
    return ch

def _pairs_to_peak_keys_auto(pairs, residues_in_order):
    """
    AUTO-DETECT for each chain whether numbers are:
      - PDB residue numbers (resi_num)
      - OR 1-based indices into that chain's CA list

    Also supports peaks files where chain token is numeric index (1,2,3,...).
    Returns peak_keys = {(chain, resi_raw), ...}
    """
    if not pairs:
        return set()

    chain_to_ordered, chain_to_resnums, chain_order = _build_chain_maps(residues_in_order)

    per_chain = {}
    for ch, n in pairs:
        ch_norm = _normalize_chain_token(ch, chain_to_ordered, chain_order)
        per_chain.setdefault(ch_norm, []).append(int(n))

    peak_keys = set()

    for ch, nums in per_chain.items():
        if ch not in chain_to_ordered:
            continue

        ordered = chain_to_ordered[ch]  # list of (resi_raw, resi_num)
        resnum_set = chain_to_resnums.get(ch, set())
        L = len(ordered)

        as_resnum = sum(1 for x in nums if x in resnum_set)
        as_index  = sum(1 for x in nums if 1 <= x <= L)

        # tie-break => residue numbers (safer)
        use_index = (as_index > as_resnum)

        if use_index:
            for x in nums:
                if 1 <= x <= L:
                    resi_raw, _ = ordered[x - 1]
                    peak_keys.add((ch, resi_raw))
        else:
            want = set(nums)
            for resi_raw, resi_num in ordered:
                if resi_num in want:
                    peak_keys.add((ch, resi_raw))

    return peak_keys

def _ints_to_peak_keys_auto(ints_list, residues_in_order):
    """
    If a peaks file contains only integers, decide whether they are:
      - indices into residues_in_order (1-based), OR
      - residue numbers (across all chains)
    """
    if not ints_list:
        return set()

    N = len(residues_in_order)
    if N <= 0:
        return set()

    leN = sum(1 for x in ints_list if 1 <= x <= N)
    if leN >= max(1, int(0.7 * len(ints_list))):
        keys = set()
        for idx in ints_list:
            if 1 <= idx <= N:
                r = residues_in_order[idx - 1]
                keys.add((r["chain"], r["resi_raw"]))
        return keys

    want_nums = set(ints_list)
    keys = set()
    for r in residues_in_order:
        if r["resi_num"] in want_nums:
            keys.add((r["chain"], r["resi_raw"]))
    return keys


def _load_float2(path: str):
    """
    Read two-column numeric file like:
      node_id   value
    Returns (x, v) as 1D numpy arrays.
    """
    a = np.loadtxt(path, dtype=float)
    if a.ndim == 1 and a.size >= 2:
        a = a.reshape(1, -1)
    if a.shape[1] < 2:
        raise ValueError(f"{os.path.basename(path)} must have at least 2 columns.")
    x = a[:, 0].astype(float)
    v = a[:, 1].astype(float)
    return x, v


def _peakdet(v, delta, x=None):
    """
    Python port of Eli Billauer peakdet, matching MATLAB logic.
    Returns (maxtab, mintab) as Nx2 arrays: [x_pos, peak_value].
    """
    v = np.asarray(v, dtype=float).reshape(-1)
    if x is None:
        x = np.arange(1, len(v) + 1, dtype=float)
    else:
        x = np.asarray(x, dtype=float).reshape(-1)
        if len(x) != len(v):
            raise ValueError("Input vectors v and x must have same length")

    delta = float(delta)
    if not np.isfinite(delta) or delta <= 0:
        raise ValueError("Input argument DELTA must be positive")

    maxtab = []
    mintab = []

    mn = float("inf")
    mx = -float("inf")
    mnpos = float("nan")
    mxpos = float("nan")

    lookformax = True
    for i in range(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx - delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn + delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.array(maxtab, dtype=float), np.array(mintab, dtype=float)


def _decode_node_id(node_id: float, chain_order: list):
    """
    MATLAB encoding:
      node_id = res_num + 0.1 * chain_index
    where chain_index is 1..chain_no.
    """
    node_id = float(node_id)
    res_num = int(math.floor(node_id + 1e-6))
    frac = node_id - res_num
    chain_idx = int(round(frac * 10))

    # Robustness for tiny float errors (e.g., 0.099999999 -> 1)
    if chain_idx == 0 and len(chain_order) == 1:
        chain_idx = 1

    if 1 <= chain_idx <= len(chain_order):
        return chain_order[chain_idx - 1], res_num
    return None, None


def _chain_resnum_to_raws(residues_in_order):
    """
    Map chain -> resi_num -> set(resi_raw)
    (so we can mark all insertion-code variants, mirroring MATLAB's int64(res_num) behavior).
    """
    m = {}
    for r in residues_in_order:
        ch = r["chain"]
        rn = r["resi_num"]
        if rn is None:
            continue
        m.setdefault(ch, {}).setdefault(int(rn), set()).add(r["resi_raw"])
    return m


def _mcpath_peaks_like_matlab(work_dir: str, residues_in_order):
    """
    Reproduce MATLAB peaks.m behavior from closeness_float / betweenness_float:
      - peakdet with delta = std(values) (MATLAB std => ddof=1)
      - filter maxima > mean + std/4
      - decode node_id -> (chain, resnum) using 0.1 encoding
      - write *_peaks, *_peaks_jmol, close_min_max, between_min_max, peaks, peaks_jmol
    Returns: (close_peak_keys, betw_peak_keys)
    """
    close_path = os.path.join(work_dir, "closeness_float")
    betw_path  = os.path.join(work_dir, "betweenness_float")

    if not (os.path.isfile(close_path) and os.path.isfile(betw_path)):
        raise FileNotFoundError("closeness_float / betweenness_float not found in work_dir.")

    x_c, v_c = _load_float2(close_path)
    x_b, v_b = _load_float2(betw_path)

    # MATLAB std uses N-1 (sample std) by default
    std_c = float(np.std(v_c, ddof=1)) if len(v_c) > 1 else float(np.std(v_c))
    std_b = float(np.std(v_b, ddof=1)) if len(v_b) > 1 else float(np.std(v_b))

    close_peak, _ = _peakdet(v_c, std_c, x_c)    # maxtab
    betw_peak,  _ = _peakdet(v_b, std_b, x_b)

    thr_c = float(np.mean(v_c) + std_c / 4.0)
    thr_b = float(np.mean(v_b) + std_b / 4.0)

    filtered_close = close_peak[close_peak[:, 1] > thr_c] if close_peak.size else np.zeros((0, 2))
    filtered_betw  = betw_peak[betw_peak[:, 1] > thr_b]  if betw_peak.size  else np.zeros((0, 2))

    # Chain order (MATLAB uses chainId(i) ordering). Here we infer from PDB CA order.
    _, _, chain_order = _build_chain_maps(residues_in_order)
    resmap = _chain_resnum_to_raws(residues_in_order)

    def _write_simple_peaks(fname_base: str, filt_arr: np.ndarray, empty_msg: str):
        txt_path  = os.path.join(work_dir, fname_base)
        jmol_path = os.path.join(work_dir, f"{fname_base}_jmol")

        if filt_arr.shape[0] < 1:
            with open(txt_path, "w") as f:
                f.write(empty_msg)
            with open(jmol_path, "w") as f:
                f.write(empty_msg)
            return []

        decoded = []
        for node_id, val in filt_arr:
            ch, rn = _decode_node_id(node_id, chain_order)
            if ch is None:
                continue
            decoded.append((ch, int(rn), float(val)))

        with open(txt_path, "w") as f:
            for ch, rn, _ in decoded:
                f.write(f"{int(rn)}{ch}\n")

        with open(jmol_path, "w") as f:
            for ch, rn, _ in decoded:
                f.write(f"{int(rn)}{ch} ALA\n")

        return decoded

    betw_dec  = _write_simple_peaks("betweenness_peaks", filtered_betw,  "There is no betweenness peak found.")
    close_dec = _write_simple_peaks("closeness_peaks",   filtered_close, "There is no closeness peak found.")

    # close_min_max / between_min_max
    with open(os.path.join(work_dir, "close_min_max"), "w") as f:
        f.write(f"{float(np.max(v_c) + 0.01)}\n{float(np.min(v_c) - 0.01)}\n")
    with open(os.path.join(work_dir, "between_min_max"), "w") as f:
        f.write(f"{float(np.max(v_b) + 0.01)}\n{float(np.min(v_b) - 0.01)}\n")

    # Combined peaks file:
    # MATLAB code has a bug (abs(...)<Inf always true). A meaningful equivalent is to JOIN on the same node_id.
    # We use 1-decimal rounding because encoding is res + 0.1*chainIndex.
    close_map = {round(float(node), 1): float(val) for (node, val) in filtered_close} if filtered_close.size else {}
    betw_map  = {round(float(node), 1): float(val) for (node, val) in filtered_betw}  if filtered_betw.size  else {}
    common_nodes = sorted(set(close_map).intersection(set(betw_map)))

    peaks_path = os.path.join(work_dir, "peaks")
    peaks_jmol = os.path.join(work_dir, "peaks_jmol")

    if not common_nodes:
        with open(peaks_path, "w") as f:
            f.write("There is no closeness or betweenness peak found.")
        with open(peaks_jmol, "w") as f:
            f.write("There is no closeness or betweenness peak found.")
    else:
        with open(peaks_path, "w") as f:
            for node in common_nodes:
                ch, rn = _decode_node_id(node, chain_order)
                if ch is None:
                    continue
                f.write(f"{int(rn)}{ch}\t{close_map[node]:.6f}\t{betw_map[node]:.6f}\n")
        with open(peaks_jmol, "w") as f:
            for node in common_nodes:
                ch, rn = _decode_node_id(node, chain_order)
                if ch is None:
                    continue
                f.write(f"{int(rn)}{ch} ALA\n")

    # Build peak_keys sets for B-factor PDB writing
    close_keys = set()
    for node_id, _ in filtered_close:
        ch, rn = _decode_node_id(node_id, chain_order)
        if ch is None:
            continue
        for resi_raw in resmap.get(ch, {}).get(int(rn), set()):
            close_keys.add((ch, resi_raw))

    betw_keys = set()
    for node_id, _ in filtered_betw:
        ch, rn = _decode_node_id(node_id, chain_order)
        if ch is None:
            continue
        for resi_raw in resmap.get(ch, {}).get(int(rn), set()):
            betw_keys.add((ch, resi_raw))

    return close_keys, betw_keys



# -------------------- B-factor PDB (FULL protein; peaks encoded in B) --------------------
def _write_bfactor_peak_pdb(pdb_in: str, pdb_out: str, peak_keys: set,
                            peak_b: float = 100.0, other_b: float = 0.0):
    os.makedirs(os.path.dirname(pdb_out) or ".", exist_ok=True)

    with open(pdb_in, "r", encoding="utf-8", errors="ignore") as fin, \
         open(pdb_out, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                ch = line[21].strip()
                resi = line[22:26].strip()
                icode = line[26].strip()
                resi_raw = f"{resi}{icode}" if icode else resi

                b = peak_b if (ch, resi_raw) in peak_keys else other_b
                b_str = f"{float(b):6.2f}"

                if len(line) < 66:
                    line = line.rstrip("\n").ljust(66) + "\n"
                line = line[:60] + b_str + line[66:]
            fout.write(line)

    return pdb_out

# -------------------- Viewer (UNCHANGED behavior: red spheres from B-factor PDB) --------------------
def _extract_peak_resis_from_bfactor_pdb(pdb_path: str, b_threshold: float = 50.0):
    out = {}
    seen = set()
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            if len(line) < 66:
                continue

            ch = line[21].strip()
            resi = line[22:26].strip()
            icode = line[26].strip()
            resi_raw = f"{resi}{icode}" if icode else resi

            try:
                b = float(line[60:66])
            except Exception:
                continue
            if b < b_threshold:
                continue

            m = re.match(r"\d+", resi_raw)
            if not m:
                continue
            rn = int(m.group(0))

            key = (ch, resi_raw)
            if key in seen:
                continue
            seen.add(key)

            out.setdefault(ch, set()).add(rn)

    return {ch: sorted(list(s)) for ch, s in out.items()}

def _view_bfactor_pdb_py3dmol(pdb_path: str, title: str, sphere_radius: float = 0.9, show_sticks: bool = False):
    """
    Shows FULL protein (cartoon) + RED spheres for peak residues (B>=50).
    Call this OUTSIDE ipywidgets.Output capture.
    """
    if not _ensure_py3dmol():
        print("Warning: py3Dmol not available; skipping viewer.")
        return

    _enable_colab_widgets()
    import py3Dmol

    pdb_text = open(pdb_path, "r", encoding="utf-8", errors="ignore").read()
    peaks = _extract_peak_resis_from_bfactor_pdb(pdb_path, b_threshold=50.0)

    view = py3Dmol.view(width=850, height=520)
    view.addModel(pdb_text, "pdb")
    view.setStyle({}, {"cartoon": {"color": "lightgray"}})

    for ch, resis in peaks.items():
        if not resis:
            continue
        if show_sticks:
            view.addStyle({"chain": ch, "resi": resis}, {"stick": {"color": "red"}})
        view.addStyle({"chain": ch, "resi": resis, "atom": "CA"},
                      {"sphere": {"color": "red", "radius": float(sphere_radius)}})

    view.zoomTo()
    display(HTML(f"<b>{title}</b>"))
    view.show()

# -------------------- NEW: Peak-only PDBs for packaging (does NOT affect viewers) --------------------
def _write_peak_only_pdb(pdb_in: str, pdb_out: str, peak_keys: set, peak_b: float = 100.0):
    """
    Write a PDB that contains ONLY atoms belonging to residues in peak_keys.
    Sets B-factor=peak_b for all kept atoms.
    This is only for saving/archiving; NOT used for the viewer.
    """
    os.makedirs(os.path.dirname(pdb_out) or ".", exist_ok=True)

    wrote_any = False
    with open(pdb_in, "r", encoding="utf-8", errors="ignore") as fin, \
         open(pdb_out, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            ch = line[21].strip()
            resi = line[22:26].strip()
            icode = line[26].strip()
            resi_raw = f"{resi}{icode}" if icode else resi

            if (ch, resi_raw) not in peak_keys:
                continue

            wrote_any = True
            b_str = f"{float(peak_b):6.2f}"
            if len(line) < 66:
                line = line.rstrip("\n").ljust(66) + "\n"
            line = line[:60] + b_str + line[66:]
            fout.write(line)

        if wrote_any:
            fout.write("END\n")

    return pdb_out if wrote_any else None

def _collect_chain_plot_files(work_dir: str):
    """
    Collect betweenness_chain_plot* and closeness_chain_plot* (any extension).
    """
    pats = [
        "closeness_chain_plot*",
        "betweenness_chain_plot*",
        "*closeness*chain*plot*",
        "*betweenness*chain*plot*",
    ]
    out = []
    for pat in pats:
        out.extend(glob.glob(os.path.join(work_dir, pat)))

    uniq = []
    seen = set()
    for p in out:
        ap = os.path.abspath(p)
        if os.path.isfile(ap) and ap not in seen:
            seen.add(ap)
            uniq.append(ap)
    return uniq

def _make_rar_or_zip(work_dir: str, base_name: str, file_paths: list):
    """
    Create .rar if possible, else .zip (fallback).
    Returns archive path or None.
    """
    file_paths = [os.path.abspath(p) for p in file_paths if p and os.path.isfile(p)]
    if not file_paths:
        return None

    rar_path = os.path.join(work_dir, f"{base_name}.rar")
    zip_path = os.path.join(work_dir, f"{base_name}.zip")

    # Try RAR
    rar_exe = shutil.which("rar")
    if rar_exe is None:
        try:
            subprocess.check_call(["apt-get", "update", "-qq"])
            subprocess.check_call(["apt-get", "install", "-y", "rar"])
            rar_exe = shutil.which("rar")
        except Exception:
            rar_exe = None

    if rar_exe is not None:
        try:
            # -ep1 => store only filenames
            cmd = [rar_exe, "a", "-idq", "-ep1", rar_path] + file_paths
            subprocess.check_call(cmd)
            if os.path.isfile(rar_path):
                return rar_path
        except Exception:
            pass

    # ZIP fallback
    try:
        with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
            for p in file_paths:
                z.write(p, arcname=os.path.basename(p))
        return zip_path if os.path.isfile(zip_path) else None
    except Exception:
        return None

def _download_file(path: str):
    """
    Triggers a real browser download in Colab.
    No clickable links are displayed.
    """
    if not path or (not os.path.isfile(path)):
        print("Warning: archive not found; nothing to download.")
        return
    try:
        from google.colab import files
        files.download(path)
    except Exception:
        # no links; just a message
        print("Archive saved at:", path)

# -------------------- UI helpers --------------------
def _list_or_custom_row(label: str, options, default_value, minv, maxv, step=1):
    options = sorted({int(x) for x in options})
    default_value = int(default_value)
    if default_value not in options:
        options = sorted(options + [default_value])

    lbl = W.Label(f"{label if label.endswith(':') else label + ':'}", layout=W.Layout(width="180px"))

    toggle = W.ToggleButtons(
        options=[("List", "list"), ("Custom", "custom")],
        value="list",
        description="",
        style={'button_width': '120px'},
        layout=W.Layout(display="flex", flex_flow="row", width="260px", margin="0 12px 0 0")
    )
    dropdown = W.Dropdown(options=options, value=default_value,
                          layout=W.Layout(width="280px", min_width="220px"))
    intbox   = W.BoundedIntText(value=default_value, min=minv, max=maxv, step=step,
                                layout=W.Layout(width="280px", min_width="220px"))
    value_box = W.Box([dropdown], layout=W.Layout(align_items="center"))

    def _on_mode(ch):
        value_box.children = [dropdown] if ch["new"] == "list" else [intbox]
    toggle.observe(_on_mode, names="value")

    def get_value():
        return int(dropdown.value if toggle.value == "list" else intbox.value)

    def set_list(val):
        toggle.value = "list"
        if val not in dropdown.options:
            dropdown.options = sorted({int(x) for x in list(dropdown.options) + [int(val)]})
        dropdown.value = int(val)

    def set_custom(val):
        toggle.value = "custom"
        intbox.value = int(val)

    row = W.HBox([lbl, toggle, value_box], layout=W.Layout(align_items="center", gap="12px", width="auto"))
    return row, {'toggle': toggle, 'dropdown': dropdown, 'intbox': intbox,
                 'get': get_value, 'set_list': set_list, 'set_custom': set_custom}

# -------------------- main UI --------------------
def launch(
    defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
    show_title="MCPath-style Parameters"
):
    global _FORM_ROOT, _LOG_OUT

    _enable_colab_widgets()

    cfg = yaml.safe_load(requests.get(defaults_url, timeout=30).text)
    target_url = cfg.get("target_url", "").strip()
    FN = cfg["field_names"]

    SAVE_DIR = cfg.get("save_dir", "/content")
    os.makedirs(SAVE_DIR, exist_ok=True)

    logo = _logo_widget(cfg.get("branding", {}))

    DESC = {'description_width': '180px'}
    wide = W.Layout(width="420px", min_width="360px")

    pdb_code   = W.Text(value=str(cfg.get("pdb_code", "")), description="PDB code:", layout=wide, style=DESC)
    or_lbl     = W.HTML("<b>&nbsp;&nbsp;or&nbsp;&nbsp;</b>")
    btn_submit = W.Button(description="Submit", button_style="success", icon="paper-plane")
    btn_load   = W.Button(description="Load", button_style="", icon="refresh")
    btn_submit.disabled = True  # Load yapılmadan submit olmasın
    
    # --- Global chain selection (MULTI via checkboxes) ---
    chain_all_cb = W.Checkbox(
        value=True,
        description="All chains",
        indent=False,
        layout=W.Layout(width="200px")
    )
    
    chain_checks_box = W.VBox([])  # Load sonrası chain checkboxları burada oluşacak
    
    chain_global_box = W.VBox(
        [W.HTML("<b>Chain ID selection:</b>"), chain_all_cb, chain_checks_box],
        layout=W.Layout(border="1px solid #eee", padding="8px", width="auto")
    )

    
    view_out = W.Output()  # py3Dmol veya kısa info göstermek için

    
    btn_new_job = W.Button(description="New Job",  button_style="info", icon="plus")
    # --- Upload state (shared) ---

    
    _UPLOADED_PDB = {"content": None, "name": None}
    
    file_lbl = W.Label("No file chosen")
    
    def _is_colab():
        try:
            import google.colab  # noqa
            return True
        except Exception:
            return False
    
    if _is_colab():
        btn_upload = W.Button(description="Upload PDB", icon="upload", button_style="")
        pdb_upload = None  # Colab mode: ipywidgets FileUpload kullanmıyoruz
    
        def _on_click_upload(_):
            btn_submit.disabled = True
            file_lbl.value = "Uploading..."
            try:
                from google.colab import files
                uploaded = files.upload()
                if not uploaded:
                    file_lbl.value = "No file chosen"
                    _UPLOADED_PDB["content"] = None
                    _UPLOADED_PDB["name"] = None
                    return
    
                name, content = next(iter(uploaded.items()))
                _UPLOADED_PDB["name"] = name
                _UPLOADED_PDB["content"] = content
                file_lbl.value = name
    
            except Exception as e:
                file_lbl.value = "Upload failed"
                _UPLOADED_PDB["content"] = None
                _UPLOADED_PDB["name"] = None
                print("Upload error:", e)
    
            finally:
                btn_submit.disabled = False
    
        btn_upload.on_click(_on_click_upload)
    
    else:
        pdb_upload = W.FileUpload(accept=".pdb", multiple=False, description="Choose file")
        btn_upload = None
    
        def _on_upload_change(_):
            content, name = _get_upload_bytes_and_name(pdb_upload)
            if content is None:
                file_lbl.value = "No file chosen"
                _UPLOADED_PDB["content"] = None
                _UPLOADED_PDB["name"] = None
                return
    
            if isinstance(content, memoryview):
                content = content.tobytes()
    
            _UPLOADED_PDB["content"] = content
            _UPLOADED_PDB["name"] = name or "upload.pdb"
            file_lbl.value = _UPLOADED_PDB["name"]
    
        pdb_upload.observe(_on_upload_change, names="value")

    


    chain_id = W.Text(
        value=str(cfg.get("chain_id", "")),
        description="Chain ID:",
        placeholder="A or A,B (empty = all chains)",
        layout=wide, style=DESC
    )
    email = W.Text(value=str(cfg.get("email", "")), description="Email (opt):", layout=wide, style=DESC)

    pred_type = W.RadioButtons(
        options=[
            ("Functional Residues", "functional"),
            ("Allosteric Paths (initial residue + path length)", "paths_init_len"),
            ("Allosteric Paths (initial & final residues)", "paths_init_final"),
        ],
        value="functional",
        description="Prediction:",
        layout=W.Layout(width="auto"),
        style=DESC
    )

    big_opts = cfg.get("path_length_options", [100000, 200000, 300000, 400000, 500000])
    big_default = int(cfg.get("path_length", big_opts[0] if big_opts else 100000))
    row_big, big_ctrl = _list_or_custom_row("Path length", big_opts, big_default, 1, 10_000_000, 1000)
    get_big_len = big_ctrl['get']

    init_idx = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1,
                                description="Index of initial residue:", layout=wide, style=DESC)
    init_chain_rb = W.RadioButtons(
        options=[],
        value=None,
        description="Chain of initial residue:",
        layout=W.Layout(width="auto"),
        style=DESC
    )

    short_len_opts = [5, 8, 10, 13, 15, 20, 25, 30]
    row_short, short_ctrl = _list_or_custom_row("Length of Paths", short_len_opts, 5, 1, 10_000, 1)
    get_short_len = short_ctrl['get']

    num_paths_opts_mode2 = [1000, 2000, 3000, 5000, 10000]
    row_np2, np2_ctrl = _list_or_custom_row("Number of Paths", num_paths_opts_mode2, 1000, 1, 10_000_000, 100)
    get_num_paths_2 = np2_ctrl['get']

    final_idx = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1,
                                 description="Index of final residue:", layout=wide, style=DESC)
    final_chain_rb = W.RadioButtons(
        options=[],
        value=None,
        description="Chain of final residue:",
        layout=W.Layout(width="auto"),
        style=DESC
    )

    num_paths_opts_mode3 = [1000, 2000, 3000, 5000, 10000]
    row_np3, np3_ctrl = _list_or_custom_row("Number of Paths", num_paths_opts_mode3, 1000, 1, 10_000_000, 100)
    get_num_paths_3 = np3_ctrl['get']



    if _LOG_OUT is None:
        _LOG_OUT = W.Output()
    out = _LOG_OUT

    upload_widget = btn_upload if btn_upload is not None else pdb_upload
    pdb_row = W.HBox(
        [pdb_code, W.HTML("&nbsp;"), or_lbl, upload_widget, file_lbl, btn_load],
        layout=W.Layout(align_items="center", justify_content="flex-start", flex_flow="row wrap", gap="10px")
    )

    functional_box = W.VBox([row_big])
    mode2_box = W.VBox([init_idx, init_chain_rb, row_short, row_np2])
    mode3_box = W.VBox([init_idx, init_chain_rb, final_idx, final_chain_rb, row_np3])

    def _sync_mode(*_):
        functional_box.layout.display = ""
        mode2_box.layout.display = "none"
        mode3_box.layout.display = "none"
        if pred_type.value == "paths_init_len":
            functional_box.layout.display = "none"; mode2_box.layout.display = ""
        elif pred_type.value == "paths_init_final":
            functional_box.layout.display = "none"; mode2_box.layout.display = "none"; mode3_box.layout.display = ""
    _sync_mode()
    pred_type.observe(_sync_mode, names="value")

    body_children = []
    if logo:
        body_children.append(logo)
    body_children += [
        W.HTML(f"<h3>{show_title}</h3>"),
        pdb_row,
    
        chain_global_box,   # <-- buraya eklendi
    
        W.HTML("<b>Prediction of:</b>"),
        pred_type,
        functional_box,
        mode2_box,
        mode3_box,
    
        chain_id, email,    # chain_id backend’de kalsın
        view_out,
        W.HBox([btn_submit, btn_new_job]),
        W.HTML("<hr>"),
        out
    ]

    root = W.VBox(body_children, layout=W.Layout(width="auto"))
    _FORM_ROOT = root
    display(root)

    # -------------------- handlers --------------------

    def _get_selected_global_chains():
        """
        Returns list of selected chains.
        If All is checked OR none selected => [] meaning "all chains"
        """
        if chain_all_cb.value:
            return []
        sels = []
        for cb in chain_checks_box.children:
            if getattr(cb, "value", False):
                sels.append(cb.description)
        return sels if sels else []  # none => treat as all chains
    
    def _sync_chain_id_and_init_final_options():
        # chain_id sync (backend uses this)
        sel = _get_selected_global_chains()
        chain_id.value = ",".join(sel) if sel else ""
    
        all_chains = [cb.description for cb in chain_checks_box.children]
        allowed = sel if sel else all_chains
    
        init_chain_rb.options  = [(f"Chain {c}", c) for c in allowed]
        final_chain_rb.options = [(f"Chain {c}", c) for c in allowed]
    
        if allowed:
            if init_chain_rb.value not in allowed:
                init_chain_rb.value = allowed[0]
            if final_chain_rb.value not in allowed:
                final_chain_rb.value = allowed[0]
        else:
            init_chain_rb.value = None
            final_chain_rb.value = None
    
    def _on_all_chains_toggle(change):
        if change.get("name") != "value":
            return
        if bool(change["new"]):
            for cb in chain_checks_box.children:
                cb.value = True
        _sync_chain_id_and_init_final_options()
    
    chain_all_cb.observe(_on_all_chains_toggle, names="value")
    
    def _on_any_chain_checkbox(change):
        if change.get("name") != "value":
            return
        if chain_checks_box.children:
            all_selected = all(cb.value for cb in chain_checks_box.children)
            # keep All chains consistent, but avoid recursion loops
            chain_all_cb.unobserve(_on_all_chains_toggle, names="value")
            chain_all_cb.value = all_selected
            chain_all_cb.observe(_on_all_chains_toggle, names="value")
        _sync_chain_id_and_init_final_options()


    def _update_chain_ui_and_view(pdb_bytes: bytes, pdb_name: str):
        chains = _chains_from_pdb_bytes(pdb_bytes)
    
        # Chain yoksa bile submit aç (tek chain gibi davran)
        if not chains:
            chain_checks_box.children = ()
            chain_id.value = ""
            btn_submit.disabled = False
            with view_out:
                clear_output(wait=True)
                print("Loaded PDB but no chain IDs detected (will run as 'all chains').")
            return
    
        # Create chain checkboxes (default checked)
        cbs = []
        for ch in chains:
            cb = W.Checkbox(value=True, description=ch, indent=False, layout=W.Layout(width="90px"))
            cb.observe(_on_any_chain_checkbox, names="value")
            cbs.append(cb)
        chain_checks_box.children = tuple(cbs)
    
        # Default: All chains ON
        chain_all_cb.value = True
        _sync_chain_id_and_init_final_options()
    
        btn_submit.disabled = False
    
        # py3Dmol preview (selected chains red; if "all", highlight all)
        with view_out:
            clear_output(wait=True)
            if not _ensure_py3dmol():
                print(f"Loaded {pdb_name}. Chains: {', '.join(chains)}")
                return
    
            _enable_colab_widgets()
            import py3Dmol
    
            pdb_text = pdb_bytes.decode("utf-8", errors="ignore")
            v = py3Dmol.view(width=850, height=420)
            v.addModel(pdb_text, "pdb")
            v.setStyle({}, {"cartoon": {"color": "lightgray"}})
    
            sel = _get_selected_global_chains()
            highlight = chains if (not sel) else sel
            for ch in highlight:
                v.addStyle({"chain": ch}, {"cartoon": {"color": "red"}})
    
            v.zoomTo()
            display(HTML(f"<b>Loaded:</b> {pdb_name} &nbsp; | &nbsp; <b>Chains:</b> {', '.join(chains)}"))
            v.show()
    



    def on_load(_):
        _LOG_OUT.clear_output(wait=True)
        with _LOG_OUT:
            try:
                pdb_bytes, pdb_name = _collect_pdb_bytes()
                _progress(0, 1, f"PDB loaded into UI: {pdb_name}")
            except Exception as e:
                print("Error:", e)
                btn_submit.disabled = True
                return
    
        _update_chain_ui_and_view(pdb_bytes, pdb_name)
    
    btn_load.on_click(on_load)

        
    def on_new_job(_):
        global _FORM_ROOT, _LOG_OUT
        try:
            if _LOG_OUT:
                _LOG_OUT.clear_output(wait=True)
        except Exception:
            pass
        try:
            if _FORM_ROOT:
                _FORM_ROOT.close()
        except Exception:
            pass
        _FORM_ROOT = None
        _LOG_OUT   = None
        clear_output(wait=True)
        launch(defaults_url=defaults_url, show_title=show_title)

    def _collect_pdb_bytes():
        # Önce: eğer upload yapılmışsa onu kullan
        if _UPLOADED_PDB.get("content") is not None:
            return _UPLOADED_PDB["content"], (_UPLOADED_PDB.get("name") or "upload.pdb")
    
        # Yoksa: PDB code ile indir
        code = pdb_code.value.strip()
        if not _is_valid_pdb_code(code):
            raise ValueError("PDB code must be exactly 4 alphanumeric characters (or upload a file).")
        return _fetch_rcsb(code), f"{code.upper()}.pdb"


    


    def _try_import_readpdb():
        try:
            from .readpdb_strict import run as run_readpdb
            return run_readpdb
        except Exception:
            try:
                from mcpath.readpdb_strict import run as run_readpdb
                return run_readpdb
            except Exception:
                return None

    def _copy_unique(src_path: str, dst_basename: str, work_dir: str):
        dst = os.path.join(work_dir, dst_basename)
        if not os.path.exists(dst):
            shutil.copyfile(src_path, dst)
            return dst
        k = 1
        while True:
            cand = f"{dst}_{k}"
            if not os.path.exists(cand):
                shutil.copyfile(src_path, cand)
                return cand
            k += 1

    def on_submit(_):
        # viewers list stays for the SAME B-factor PDBs (full protein) as before
        viewers = []                # list of (pdb_path, title, radius, sticks)
        archive_to_download = None  # only the archive download

        _LOG_OUT.clear_output(wait=True)

        with _LOG_OUT:
            try:
                total_steps = 5

                # Chain seçimi (SelectMultiple) -> hiçbir şey seçili değilse "all chains"
                sel = _get_selected_global_chains()
                chain_global = ",".join(sel) if sel else ""



                if not _is_valid_email(email.value.strip()):
                    raise ValueError("Invalid email format.")

                pdb_bytes, pdb_name = _collect_pdb_bytes()
                save_path = os.path.join(SAVE_DIR, pdb_name)
                with open(save_path, "wb") as f:
                    f.write(pdb_bytes)

                input_path = os.path.join(os.path.dirname(save_path), "mcpath_input.txt")
                mode = pred_type.value

                if mode == "functional":
                    rows = ["1", pdb_name, chain_global, str(get_big_len()), (email.value.strip() or "-")]
                elif mode == "paths_init_len":
                    if not _is_valid_chain(init_chain_rb.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    rows = ["2", pdb_name, chain_global, str(get_short_len()), str(get_num_paths_2()),
                            str(int(init_idx.value)), (init_chain_rb.value or "").strip(), (email.value.strip() or "-")]

                else:
                    if not _is_valid_chain(init_chain_rb.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    if not _is_valid_chain(final_chain_rb.value or ""):
                        raise ValueError("Chain of final residue must be a single character.")
                    rows = ["3", pdb_name, chain_global,
                            str(int(init_idx.value)), (init_chain_rb.value or "").strip(),
                            str(int(final_idx.value)), (final_chain_rb.value or "").strip(),
                            (email.value.strip() or "-")]


                with open(input_path, "w") as f:
                    for r in rows:
                        f.write(str(r).strip() + "\n")

                _progress(1, total_steps, "Input validated and PDB saved.")

                run_readpdb = _try_import_readpdb()
                cor_path = None
                if run_readpdb is None:
                    print("Warning: readpdb_strict not found; skipping .cor generation.")
                else:
                    try:
                        cor_path = run_readpdb(input_path=input_path)
                    except TypeError:
                        fallback_chain = (chain_global.split(",")[0].strip() if chain_global else "A")
                        cor_path = run_readpdb(pdb_path=save_path, chain=fallback_chain)


                    if cor_path and os.path.isfile(cor_path):
                        _progress(2, total_steps, ".cor coordinate file generated.")
                    else:
                        print("Warning: .cor file was not generated or not found.")

                # -------- Step 3: atomistic LJ --------
                try:
                    if cor_path and os.path.isfile(cor_path):
                        base_no_ext = os.path.splitext(save_path)[0]
                        want_cor = f"{base_no_ext}.cor"
                        if os.path.abspath(cor_path) != os.path.abspath(want_cor):
                            shutil.copyfile(cor_path, want_cor)

                        here = os.path.dirname(os.path.abspath(__file__))
                        param_path = os.path.join(here, "vdw_cns.param")
                        top_path   = os.path.join(here, "pdb_cns.top")

                        if not (os.path.isfile(param_path) and os.path.isfile(top_path)):
                            print("Warning: vdw_cns.param or pdb_cns.top missing; skipping atomistic LJ step.")
                        else:
                            pkg_dir = os.path.dirname(os.path.abspath(__file__))
                            root_dir = os.path.dirname(pkg_dir)
                            if root_dir not in sys.path:
                                sys.path.append(root_dir)

                            atom_mod = importlib.import_module("mcpath.atomistic")
                            importlib.reload(atom_mod)

                            base_for_atom = os.path.splitext(save_path)[0]
                            atom_mod.atomistic(
                                base_for_atom,
                                param_file=param_path,
                                top_file=top_path,
                                rcut=5.0,
                                kT=1.0,
                                save_txt=True
                            )
                            _progress(3, total_steps, "Atomistic LJ calculation completed.")
                    else:
                        print("Warning: No .cor available; skipping atomistic LJ step.")
                except Exception as e:
                    print(f"Warning: atomistic run failed: {e}")

                # -------- Step 4: infinite path + closeness/betweenness --------
                work_dir = os.path.dirname(save_path)

                if mode == "functional" and cor_path and os.path.isfile(cor_path):
                    try:
                        inf_mod = importlib.import_module("mcpath.infinite")
                        importlib.reload(inf_mod)

                        old_cwd = os.getcwd()
                        try:
                            os.chdir(work_dir)
                            steps = int(get_big_len())
                            path_arr = inf_mod.infinite(
                                os.path.basename(os.path.splitext(save_path)[0]),
                                path_length=steps,
                                pottype='1'
                            )
                            pl = path_arr.shape[1]
                            path_src = f"{os.path.basename(os.path.splitext(save_path)[0])}_atomistic_{pl}steps_infinite.path"
                            path_src = os.path.join(work_dir, path_src)
                        finally:
                            os.chdir(old_cwd)

                        cor_src  = f"{os.path.splitext(save_path)[0]}.cor"
                        atom_src = f"{os.path.splitext(save_path)[0]}_atomistic.out"

                        _copy_unique(cor_src,  "coor_file", work_dir=os.path.dirname(cor_src))
                        _copy_unique(atom_src, "atom_file", work_dir=os.path.dirname(atom_src))
                        _copy_unique(path_src, "path_file", work_dir=os.path.dirname(path_src))

                        close_mod = importlib.import_module("mcpath.closeness")
                        importlib.reload(close_mod)
                        betw_mod = importlib.import_module("mcpath.betweenness")
                        importlib.reload(betw_mod)

                        old_cwd2 = os.getcwd()
                        try:
                            os.chdir(work_dir)
                            close_mod.main()
                            betw_mod.main()
                        finally:
                            os.chdir(old_cwd2)

                        # ---- Step 4b: peaks (MATLAB-like) -> B-factor PDBs (FULL protein) + peak-only PDBs for archive ----
                        base = os.path.splitext(os.path.basename(save_path))[0]
                        residues_in_order = _ca_residues_in_order(save_path, chain_str=chain_global)
                        
                        close_peak_keys, betw_peak_keys = set(), set()
                        try:
                            close_peak_keys, betw_peak_keys = _mcpath_peaks_like_matlab(work_dir, residues_in_order)
                            print(f"[peaks.m-like] closeness keys: {len(close_peak_keys)} | betweenness keys: {len(betw_peak_keys)}")
                        except Exception as e_pk:
                            print(f"Warning: MATLAB-like peak selection failed; falling back to file parsing. Reason: {e_pk}")
                            # (Optional) fallback: keep your old parsing logic here if you want.
                        
                        # -------- CLOSENESS PDBs --------
                        close_pdb_out = os.path.join(work_dir, f"{base}_CLOSENESS_peaks_bfac.pdb")
                        close_only_pdb = None
                        if close_peak_keys:
                            _write_bfactor_peak_pdb(save_path, close_pdb_out, close_peak_keys, peak_b=100.0, other_b=0.0)
                            viewers.append((close_pdb_out, "Closeness peaks (B=100 → RED spheres)", 0.9, False))
                        
                            close_only_pdb = os.path.join(work_dir, f"{base}_CLOSENESS_peaks_ONLY.pdb")
                            close_only_pdb = _write_peak_only_pdb(save_path, close_only_pdb, close_peak_keys, peak_b=100.0)
                            if not close_only_pdb:
                                print("Warning: closeness PEAK-ONLY PDB not written (no atoms matched).")
                        else:
                            print("Warning: closeness peak set is empty; skipping closeness PDB write.")
                        
                        # -------- BETWEENNESS PDBs --------
                        betw_pdb_out = os.path.join(work_dir, f"{base}_BETWEENNESS_peaks_bfac.pdb")
                        betw_only_pdb = None
                        if betw_peak_keys:
                            _write_bfactor_peak_pdb(save_path, betw_pdb_out, betw_peak_keys, peak_b=100.0, other_b=0.0)
                            viewers.append((betw_pdb_out, "Betweenness peaks (B=100 → RED spheres)", 0.9, False))
                        
                            betw_only_pdb = os.path.join(work_dir, f"{base}_BETWEENNESS_peaks_ONLY.pdb")
                            betw_only_pdb = _write_peak_only_pdb(save_path, betw_only_pdb, betw_peak_keys, peak_b=100.0)
                            if not betw_only_pdb:
                                print("Warning: betweenness PEAK-ONLY PDB not written (no atoms matched).")
                        else:
                            print("Warning: betweenness peak set is empty; skipping betweenness PDB write.")
                        
                        

                        
                        # -------- NEW: make archive (original PDB + 2 peak-only PDBs + chain plots) --------
                        plot_files = _collect_chain_plot_files(work_dir)

                        files_to_pack = [save_path]
                        if close_only_pdb: files_to_pack.append(close_only_pdb)
                        if betw_only_pdb:  files_to_pack.append(betw_only_pdb)
                        files_to_pack += plot_files

                        archive_base = f"{base}_MCPath_package"
                        archive_to_download = _make_rar_or_zip(work_dir, archive_base, files_to_pack)
                        if archive_to_download:
                            print(f"[archive] Ready: {os.path.basename(archive_to_download)}")
                        else:
                            print("Warning: archive could not be created.")

                        _progress(4, total_steps, "Path generation and centrality analysis completed.")

                    except Exception as e_inf:
                        print(f"Warning: infinite/closeness pipeline failed: {e_inf}")
                elif mode == "functional":
                    print("Warning: Skipping path/centrality analysis (missing .cor or earlier failure).")

                # -------- Step 5: remote submission (if configured) --------
                data = {"prediction_mode": mode, FN["chain_id"]: chain_global}
                if pdb_code.value.strip():
                    data[FN["pdb_code"]] = pdb_code.value.strip().upper()
                if email.value.strip():
                    data[FN["email"]] = email.value.strip()

                if mode == "functional":
                    data[FN["path_length"]] = str(get_big_len())
                elif mode == "paths_init_len":
                    data.update({
                        "length_paths": int(get_short_len()),
                        "number_paths": int(get_num_paths_2()),
                        "index_initial": int(init_idx.value),
                        "chain_initial": (init_chain.value or "").strip()
                    })
                else:
                    data.update({
                        "index_initial": int(init_idx.value),
                        "chain_initial": (init_chain_rb.value or "").strip(),
                        "index_final": int(final_idx.value),
                        "chain_final": (final_chain_rb.value or "").strip(),
                        "number_paths": int(get_num_paths_3())
                    })

                files = {FN["pdb_file"]: (pdb_name, pdb_bytes, "chemical/x-pdb")}
                if not target_url:
                    _progress(5, total_steps, "Local processing completed (no remote submission configured).")
                else:
                    try:
                        r = requests.post(target_url, data=data, files=files, timeout=180)
                        _progress(5, total_steps, f"Remote submission completed (HTTP {r.status_code}).")
                    except Exception as e_sub:
                        print(f"Warning: remote submission failed: {e_sub}")

            except Exception as e:
                print("Error:", e)

        # ---- IMPORTANT: render py3Dmol OUTSIDE the ipywidgets.Output capture ----
        if viewers:
            display(HTML("<hr><h4>3D Views</h4>"))
            for (pdb_path, title, radius, sticks) in viewers:
                _view_bfactor_pdb_py3dmol(pdb_path, title, sphere_radius=radius, show_sticks=sticks)

        # ---- Download ONLY the archive (no links printed) ----
        if archive_to_download:
            _download_file(archive_to_download)

    btn_new_job.on_click(on_new_job)
    btn_submit.on_click(on_submit)
