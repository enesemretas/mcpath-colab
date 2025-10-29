# mcpath/ui.py
import os, re, requests, shutil
from IPython.display import display, clear_output
import ipywidgets as W
import numpy as np
import yaml

# -------------------- validators --------------------
def _is_valid_pdb_code(c): return bool(re.fullmatch(r"[0-9A-Za-z]{4}", c.strip()))
def _is_valid_chain(ch):   return bool(re.fullmatch(r"[A-Za-z0-9]", ch.strip()))
def _is_valid_email(s):    return (not s.strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s.strip()))

def _fetch_rcsb(code: str) -> bytes:
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    print(f"Fetching PDB from RCSB: {url} (≤15s timeout)")
    r = requests.get(url, timeout=15); r.raise_for_status()
    return r.content

def _list_or_custom(label: str, options, default_value, minv, maxv, step=1, desc_style=None):
    options = sorted({int(x) for x in options})
    default_value = int(default_value)
    if default_value not in options:
        options = sorted(options + [default_value])

    common_layout = W.Layout(width="460px", min_width="420px")
    mode = W.ToggleButtons(
        options=[("List", "list"), ("Custom", "custom")],
        value="list",
        description="Input:",
        layout=W.Layout(width="260px"),
        style=desc_style
    )
    dd = W.Dropdown(
        options=options, value=default_value,
        description=(label if label.endswith(":") else f"{label}:"),
        layout=common_layout, style=desc_style
    )
    txt = W.BoundedIntText(
        value=default_value, min=minv, max=maxv, step=step,
        description=(label if label.endswith(":") else f"{label}:"),
        layout=common_layout, style=desc_style
    )

    container = W.VBox([dd])
    def _on_mode(change):
        container.children = [dd] if change["new"] == "list" else [txt]
    mode.observe(_on_mode, names="value")

    def get_value():
        return int(dd.value if mode.value == "list" else txt.value)

    return mode, container, dd, txt, get_value

def _logo_widget(branding: dict):
    url = (branding.get("logo_url") or "").strip()
    if not url: return None
    height = int(branding.get("logo_height", 96))
    align  = (branding.get("logo_align") or "center").lower()
    try:
        img_bytes = requests.get(url, timeout=10).content
        img = W.Image(value=img_bytes, format="png", layout=W.Layout(height=f"{height}px"))
        jc = {"center":"center", "left":"flex-start", "right":"flex-end"}.get(align, "center")
        return W.HBox([img], layout=W.Layout(justify_content=jc))
    except Exception:
        return None

# -------------------- Python port of readpdb.m --------------------
AA_CODE = {
    "GLY": 1, "ALA": 2, "VAL": 3, "ILE": 4, "LEU": 5,
    "SER": 6, "THR": 7, "ASP": 8, "ASN": 9, "GLU":10,
    "GLN":11, "LYS":12, "ARG":13, "CYS":14, "MET":15,
    "PHE":16, "TYR":17, "TRP":18, "HIS":19, "PRO":20
}
AA_CONST = {
    1:  75.0,  2:  89.1,  3: 117.2,  4: 131.2,  5: 131.2,
    6: 105.1,  7: 119.1,  8: 133.1,  9: 132.1, 10: 147.1,
   11: 146.0, 12: 146.2, 13: 174.2, 14: 121.2, 15: 149.2,
   16: 165.2, 17: 181.2, 18: 204.2, 19: 155.2, 20: 115.1,
}
AA_MIN_ATOMS = {
    1:  2,  2:  5,  3:  7,  4:  8,  5:  8,
    6:  6,  7:  6,  8:  8,  9:  8, 10:  9,
   11:  9, 12:  9, 13: 11, 14:  6, 15:  7,
   16: 11, 17: 12, 18: 14, 19: 10, 20:  7,
}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
ALTNAME_NORMALIZE = {
    "HIP":"HIS", "HIE":"HIS", "HID":"HIS", "HSD":"HIS", "HSE":"HIS", "HSP":"HIS",
    "LYN":"LYS", "GLH":"GLU", "CYM":"CYS", "CYX":"CYS", "ASH":"ASN", "TYM":"TYR"
}

def _trim_to_first_model(pdb_lines):
    has_model = any(l.startswith("MODEL") for l in pdb_lines)
    if not has_model: return pdb_lines
    new_lines = []; in_model = False
    for L in pdb_lines:
        if L.startswith("MODEL"):
            if in_model: break
            in_model = True
            continue
        if L.startswith("ENDMDL") and in_model: break
        if not in_model:
            if L.startswith("ATOM") or L.startswith("HETATM") or L.startswith("TER"):
                new_lines.append(L)
        else:
            new_lines.append(L)
    return new_lines

def _safe_sub(s, a, b):
    a = max(1, a); b = min(len(s), b)
    return s[a-1:b] if a <= b else ""

def _parse_atoms_from_path(pdb_path):
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.read().splitlines()
    lines = _trim_to_first_model(lines)
    atoms = []
    for L in lines:
        if not (L.startswith("ATOM") or L.startswith("HETATM")): continue
        if len(L) < 54: continue
        try:
            AtomSerNo = int(_safe_sub(L,7,11).strip() or "0")
            AtomName  = _safe_sub(L,13,16).strip()
            altLoc    = _safe_sub(L,17,17).strip()
            resName   = _safe_sub(L,18,20).strip()
            chainID   = _safe_sub(L,22,22).strip()
            resSeq    = int(_safe_sub(L,23,26).strip() or "0")
            iCode     = _safe_sub(L,27,27).strip()
            X         = float(_safe_sub(L,31,38))
            Y         = float(_safe_sub(L,39,46))
            Z         = float(_safe_sub(L,47,54))
            element   = _safe_sub(L,77,78).strip()
        except: 
            continue
        atoms.append({
            "AtomSerNo": AtomSerNo, "AtomName": AtomName, "altLoc": altLoc,
            "resName": resName, "chainID": chainID, "resSeq": resSeq,
            "iCode": iCode, "X": X, "Y": Y, "Z": Z, "element": element
        })
    return atoms

def _normalize_resname(nm): 
    u = (nm or "").upper(); 
    return ALTNAME_NORMALIZE.get(u, u)

def _is_hydrogen(atom):
    el = (atom["element"] or "").upper()
    nm = (atom["AtomName"] or "").upper()
    if el == "H" or nm == "H": return True
    if nm.startswith("H"): return True
    return False

def _chain_rank_map(atoms):
    ranks = {}; r = 0; last = None
    for a in atoms:
        ch = (a["chainID"] or "").upper()
        if ch != last:
            if ch not in ranks:
                r += 1; ranks[ch] = r
            last = ch
    return ranks

def _sidechain_coord(resname, atom_group):
    rn = (resname or "").upper()
    ca = next((a for a in atom_group if (a["AtomName"] or "").upper() == "CA"), None)
    if rn == "GLY":
        return (ca["X"], ca["Y"], ca["Z"]) if ca else (np.nan, np.nan, np.nan)
    cb = next((a for a in atom_group if (a["AtomName"] or "").upper() == "CB"), None)
    if cb is not None:
        return (cb["X"], cb["Y"], cb["Z"])
    pts = [(a["X"], a["Y"], a["Z"]) for a in atom_group
           if (a["AtomName"] or "").upper() not in BACKBONE and not _is_hydrogen(a)]
    if pts:
        arr = np.array(pts, dtype=float)
        c = arr.mean(axis=0)
        return (float(c[0]), float(c[1]), float(c[2]))
    return (ca["X"], ca["Y"], ca["Z"]) if ca else (np.nan, np.nan, np.nan)

def readpdb_py(pdb_path: str, chainID: str | None):
    """
    Python port of your readpdb.m:
    - Writes "<pdb_path>.cor"
    - If residues fail atom-count thresholds, writes "problem".
    """
    chainID = (chainID or "").strip()[:1].upper()
    atoms_all = _parse_atoms_from_path(pdb_path)

    # Normalize names, drop hydrogens, handle altLoc (keep first seen per chain/res/atom)
    cleaned = []
    seen_alt = set()
    for a in atoms_all:
        if _is_hydrogen(a): 
            continue
        a["resName"] = _normalize_resname(a["resName"])
        a["chainID"] = (a["chainID"] or "").upper()
        key = (a["chainID"], a["resSeq"], (a["AtomName"] or "").upper())
        if a["altLoc"]:
            if key in seen_alt: 
                continue
            seen_alt.add(key)
        cleaned.append(a)

    ranks_full = _chain_rank_map(cleaned)
    if chainID and chainID not in ranks_full:
        print(f'WARNING: requested chain "{chainID}" not found. No chain filter applied.')
        chosen = cleaned
        ranks = ranks_full
    else:
        chosen = [a for a in cleaned if (not chainID) or a["chainID"] == chainID]
        ranks = ranks_full

    # Group by (chain,resSeq) in file order
    groups = []
    last = None
    for a in chosen:
        key = (a["chainID"], a["resSeq"])
        if key != last:
            groups.append([a["chainID"], a["resSeq"], [a]])
            last = key
        else:
            groups[-1][2].append(a)

    # Build rows for standard AAs with CA
    result_rows = []
    problems = []
    for ch, resi, grp in groups:
        rn = grp[0]["resName"]
        if rn not in AA_CODE: 
            continue
        ca = next((a for a in grp if (a["AtomName"] or "").upper() == "CA"), None)
        if ca is None:
            continue
        aa_code = AA_CODE[rn]
        chain_rank = ranks.get(ch, 0)
        res_atom_count = sum(1 for a in grp if not _is_hydrogen(a))
        scx, scy, scz = _sidechain_coord(rn, grp)
        const = AA_CONST[aa_code]
        result_rows.append([
            float(resi), float(aa_code),
            float(ca["X"]), float(ca["Y"]), float(ca["Z"]),
            float(scx), float(scy), float(scz),
            float(res_atom_count), float(chain_rank), float(const)
        ])
        if res_atom_count <= AA_MIN_ATOMS.get(aa_code, 0):
            problems.append((resi, chain_rank))

    # Write "<ProteinName>.cor"
    cor_path = f"{pdb_path}.cor"
    with open(cor_path, "w", encoding="utf-8") as f:
        for r in result_rows:
            f.write(f"{r[0]:5.0f} {r[1]:4.0f} {r[2]:8.3f} {r[3]:8.3f} {r[4]:8.3f} "
                    f"{r[5]:8.3f} {r[6]:8.3f} {r[7]:8.3f} {r[8]:4.0f} {r[9]:2.0f} {r[10]:4.1f}\n")
    if problems:
        with open("problem", "w", encoding="utf-8") as f:
            f.write("There are problems about number of atoms on following residues\n")
            for resseq, crank in problems:
                f.write(f"{int(resseq)} of {int(crank)} chain\n")
    return cor_path, len(result_rows), len(problems)

def _copy_to_free_name(src: str, dst_base: str) -> str | None:
    """Copy src -> dst_base (or dst_base.N) without overwrite. Return final path or None."""
    if not os.path.isfile(src):
        return None
    dst = dst_base
    k = 1
    while os.path.exists(dst):
        dst = f"{dst_base}.{k}"
        k += 1
    shutil.copy2(src, dst)
    return dst

# -------------------- main UI --------------------
def launch(
    defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
    show_title="MCPath-style Parameters"
):
    """Launch the MCPath parameter form (Python backend only: readpdb)."""
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=15).text)
    FN = cfg["field_names"]
    SAVE_DIR = cfg.get("save_dir", "/content")
    os.makedirs(SAVE_DIR, exist_ok=True)

    logo = _logo_widget(cfg.get("branding", {}))
    DESC = {'description_width': '260px'}
    wide = W.Layout(width="460px", min_width="420px")

    # --- PDB: code OR upload ---
    pdb_code   = W.Text(value=str(cfg.get("pdb_code", "")),
                        description="PDB code:",
                        layout=wide, style=DESC)
    or_lbl     = W.HTML("<b>&nbsp;&nbsp;or&nbsp;&nbsp;</b>")
    pdb_upload = W.FileUpload(accept=".pdb", multiple=False, description="Choose File")
    file_lbl   = W.Label("No file chosen")
    def _on_upload_change(change):
        if pdb_upload.value:
            fname = next(iter(pdb_upload.value.keys()))
            file_lbl.value = fname
        else:
            file_lbl.value = "No file chosen"
    pdb_upload.observe(_on_upload_change, names="value")

    # --- Always-present fields ---
    chain_id   = W.Text(value=str(cfg.get("chain_id", "")),
                        description="Chain ID:",
                        layout=wide, style=DESC)
    email      = W.Text(value=str(cfg.get("email", "")),
                        description="Email (opt):",
                        layout=wide, style=DESC)

    # --- (Keep radio but only 'functional' used here) ---
    pred_type = W.RadioButtons(
        options=[
            ("Functional Residues (readpdb only)", "functional"),
        ],
        value="functional",
        description="Prediction:",
        layout=W.Layout(width="auto"),
        style=DESC
    )

    # Big path length widget (kept for compatibility; not used by Python reader)
    big_opts = cfg.get("path_length_options", [100000, 200000, 300000, 400000, 500000])
    big_default = int(cfg.get("path_length", big_opts[0] if big_opts else 100000))
    (pl_mode_big, pl_container_big, pl_dd_big, pl_txt_big, get_big_len) = _list_or_custom(
        label="Path length", options=big_opts, default_value=big_default,
        minv=1, maxv=10_000_000, step=1000, desc_style=DESC
    )
    functional_box = W.VBox([pl_mode_big, pl_container_big])

    # Actions & output
    btn_submit = W.Button(description="Run (Python readpdb)", button_style="success", icon="play")
    btn_clear  = W.Button(description="Clear",  button_style="warning", icon="trash")
    out        = W.Output()

    pdb_row = W.HBox([pdb_code, W.HTML("&nbsp;"), or_lbl, pdb_upload, file_lbl],
                     layout=W.Layout(align_items="center", justify_content="flex-start",
                                     flex_flow="row wrap", gap="10px"))

    children = []
    if logo: children.append(logo)
    children += [
        W.HTML(f"<h3>{show_title}</h3>"),
        pdb_row,
        W.HTML("<b>Prediction of:</b>"),
        pred_type,
        functional_box,
        chain_id, email,
        W.HBox([btn_submit, btn_clear]),
        W.HTML("<hr>"), out
    ]
    display(W.VBox(children, layout=W.Layout(width="auto")))

    def on_clear(_):
        pdb_code.value = ""
        try:
            pdb_upload.value.clear()
        except Exception:
            pdb_upload.value = {}
        file_lbl.value = "No file chosen"
        chain_id.value = ""
        email.value = ""
        pl_mode_big.value = "list"
        try:
            pl_dd_big.value = pl_dd_big.options[0]; pl_txt_big.value = int(pl_dd_big.options[0])
        except Exception:
            pass
        with out: clear_output()

    def _collect_pdb_bytes():
        if pdb_upload.value:
            (fname, meta) = next(iter(pdb_upload.value.items()))
            return meta["content"], fname
        code = pdb_code.value.strip()
        if not _is_valid_pdb_code(code):
            raise ValueError("PDB code must be exactly 4 alphanumeric characters (or upload a file).")
        return _fetch_rcsb(code), f"{code.upper()}.pdb"

    def on_submit(_):
        with out:
            clear_output()
            try:
                print("▶ Validating inputs…")
                chain_global = chain_id.value.strip()
                if not _is_valid_chain(chain_global):
                    raise ValueError("Chain ID must be a single character (e.g., A).")
                if not _is_valid_email(email.value.strip()):
                    raise ValueError("Invalid email format.")

                print("▶ Getting PDB (upload or RCSB)…")
                pdb_bytes, pdb_name = _collect_pdb_bytes()

                # Save PDB
                os.makedirs(SAVE_DIR, exist_ok=True)
                save_path = os.path.join(SAVE_DIR, pdb_name)
                with open(save_path, "wb") as f:
                    f.write(pdb_bytes)
                print(f"Saved local copy: {save_path}")

                # Build mcpath_input.txt (for compatibility/logging)
                input_path = os.path.join(os.path.dirname(save_path), "mcpath_input.txt")
                rows = [
                    "1",              # mode (functional)
                    pdb_name,         # pdb file name
                    chain_global,     # global chain
                    str(get_big_len()),
                    email.value.strip() or "-"
                ]
                with open(input_path, "w") as f:
                    for r in rows:
                        f.write(str(r).strip() + "\n")
                print(f"Input file saved: {input_path}")

                # ---- Run Python readpdb replacement ----
                print("▶ Running Python readpdb…")
                cor_path, nrows, nprob = readpdb_py(save_path, chain_global)
                # Normalize to coor_file (without overwrite)
                coor_norm = _copy_to_free_name(cor_path, os.path.join(os.path.dirname(save_path), "coor_file"))

                print("✅ readpdb completed.")
                print(f"Rows written: {nrows}")
                if nprob:
                    print(f"⚠️ problem file written with {nprob} entries.")

                # Preview first 8 lines of coor_file (or the .cor)
                preview_path = coor_norm if coor_norm and os.path.exists(coor_norm) else cor_path
                try:
                    print("\n--- Preview (first 8 lines) ---")
                    with open(preview_path, "r") as f:
                        for i, line in enumerate(f):
                            if i >= 8: break
                            print(line.rstrip("\n"))
                except Exception as e:
                    print("Preview unavailable:", e)

                print(f"\nOutputs:")
                print(f"  • .cor file : {cor_path}")
                if coor_norm:
                    print(f"  • coor_file : {coor_norm}")
                if os.path.exists(os.path.join(os.path.dirname(save_path), "problem")):
                    print(f"  • problem   : {os.path.join(os.path.dirname(save_path), 'problem')}")

            except Exception as e:
                print("❌", e)

    btn_clear.on_click(on_clear)
    btn_submit.on_click(on_submit)
