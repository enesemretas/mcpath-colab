# mcpath/ui.py
import os, re, requests, yaml, importlib, shutil, sys, glob, subprocess
from IPython.display import display, HTML, IFrame
import ipywidgets as W

# ---------- singletons (so New Job/Submit always uses the same log panel) ----------
_FORM_ROOT = None
_LOG_OUT   = None

# -------------------- validators & helpers --------------------
def _is_valid_pdb_code(c):
    return bool(re.fullmatch(r"[0-9A-Za-z]{4}", c.strip()))

def _is_valid_chain(ch):
    """Single-chain validator (for init/final chains etc.)."""
    return bool(re.fullmatch(r"[A-Za-z0-9]", ch.strip()))

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
    return (not s.strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s.strip()))

def _fetch_rcsb(code: str) -> bytes:
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.content

def _logo_widget(branding: dict):
    url = (branding.get("logo_url") or "").strip()
    if not url:
        return None
    height = int(branding.get("logo_height", 96))
    align  = (branding.get("logo_align") or "center").lower()
    try:
        img_bytes = requests.get(url, timeout=20).content
        jc = {"center":"center", "left":"flex-start", "right":"flex-end"}.get(align, "center")
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

# -------------------- Peak -> B-factor PDB + viewer helpers --------------------
_NUM_RE = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")
# matches "16A", "1913Q", and also "16 A"

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
    Extract (chain, number) pairs from lines like:
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
    Parse lines containing e.g.:
      16A -> closeness = 0.184519
    Returns list of (value, chain, number)
    """
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = _RESCHAIN_RE.search(line)
            if not m:
                continue
            rn = int(m.group(1))
            ch = m.group(2)

            nums = _NUM_RE.findall(line)
            if not nums:
                continue
            try:
                val = float(nums[-1])
            except Exception:
                continue
            rows.append((val, ch, rn))
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
            chain_order.append(ch)  # first-seen order in the PDB/selection
        chain_to_ordered[ch].append((r["resi_raw"], r["resi_num"]))
        if r["resi_num"] is not None:
            chain_to_resnums.setdefault(ch, set()).add(r["resi_num"])

    return chain_to_ordered, chain_to_resnums, chain_order

def _normalize_chain_token(ch, chain_to_ordered, chain_order):
    """
    Convert numeric chain tokens ("1","2",...) into actual PDB chain IDs ("A","B",...)
    based on chain order in residues_in_order.

    If the PDB really has a chain named "1", we keep it as-is.
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

    Also supports peaks files where chain is given as numeric index (1,2,3,...).
    Returns peak_keys = {(chain, resi_raw), ...}
    """
    if not pairs:
        return set()

    chain_to_ordered, chain_to_resnums, chain_order = _build_chain_maps(residues_in_order)

    # group numbers per chain (after normalizing chain token)
    per_chain = {}
    for ch, n in pairs:
        ch_norm = _normalize_chain_token(ch, chain_to_ordered, chain_order)
        per_chain.setdefault(ch_norm, []).append(int(n))

    peak_keys = set()

    for ch, nums in per_chain.items():
        if ch not in chain_to_ordered:
            # still unknown chain after normalization
            continue

        ordered = chain_to_ordered[ch]  # list of (resi_raw, resi_num)
        resnum_set = chain_to_resnums.get(ch, set())
        L = len(ordered)

        # count how many look like resnums vs indices
        as_resnum = sum(1 for x in nums if x in resnum_set)
        as_index  = sum(1 for x in nums if 1 <= x <= L)

        # IMPORTANT: tie-break goes to "resnum" (safer when both look plausible)
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
    # if most ints are <= N, likely indices
    if leN >= max(1, int(0.7 * len(ints_list))):
        keys = set()
        for idx in ints_list:
            if 1 <= idx <= N:
                r = residues_in_order[idx - 1]
                keys.add((r["chain"], r["resi_raw"]))
        return keys

    # else: residue numbers without chain -> mark any chain matching that resi_num
    want_nums = set(ints_list)
    keys = set()
    for r in residues_in_order:
        if r["resi_num"] in want_nums:
            keys.add((r["chain"], r["resi_raw"]))
    return keys

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
                b_str = f"{b:6.2f}"

                if len(line) < 66:
                    line = line.rstrip("\n").ljust(66) + "\n"
                line = line[:60] + b_str + line[66:]
            fout.write(line)

    return pdb_out

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

def _iframe_src_for(path: str) -> str:
    ap = os.path.abspath(path)
    # Colab serves /content via /files/
    if ap.startswith("/content/"):
        rel = ap[len("/content/"):].lstrip("/")
        return f"/files/{rel}"
    # fallback
    return ap


def _view_bfactor_pdb_py3dmol(pdb_path: str, title: str, sphere_radius: float = 0.9, show_sticks: bool = False):
    """
    Robust viewer:
      - tries display(view)
      - if blank (common in Output widgets), writes html and shows via IFrame
    No download links.
    """
    if not _ensure_py3dmol():
        print("Warning: py3Dmol not available; skipping viewer.")
        return

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

    # Attempt direct render
    try:
        display(view)
        return
    except Exception:
        pass

    # Fallback: iframe (most reliable inside ipywidgets.Output)
    html_path = os.path.splitext(pdb_path)[0] + "_py3dmol.html"
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(view._make_html())

    display(IFrame(src=_iframe_src_for(html_path), width=870, height=560))

def _write_pymol_pml(pml_path: str, pdb_close: str, pdb_betw: str):
    pml = f"""
reinitialize
bg_color white
set ray_opaque_background, off

load {pdb_close}, closeness
load {pdb_betw}, betweenness

hide everything, all

# --- closeness view ---
enable closeness
disable betweenness
show cartoon, closeness
color gray80, closeness
select c_peaks, closeness and b > 50
color red, c_peaks
show spheres, c_peaks and name CA
set sphere_scale, 0.6, c_peaks
orient closeness
scene closeness, store

# --- betweenness view ---
disable closeness
enable betweenness
show cartoon, betweenness
color gray80, betweenness
select b_peaks, betweenness and b > 50
color red, b_peaks
show spheres, b_peaks and name CA
set sphere_scale, 0.6, b_peaks
orient betweenness
scene betweenness, store
"""
    os.makedirs(os.path.dirname(pml_path) or ".", exist_ok=True)
    with open(pml_path, "w", encoding="utf-8") as f:
        f.write(pml.strip() + "\n")
    return pml_path

# --- Single row: "<Label> : [List | Custom]  <value-widget>"
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
    pdb_upload = W.FileUpload(accept=".pdb", multiple=False, description="Choose file")
    file_lbl   = W.Label("No file chosen")

    def _on_upload_change(_):
        if pdb_upload.value:
            file_lbl.value = next(iter(pdb_upload.value.keys()))
        else:
            file_lbl.value = "No file chosen"
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
    init_chain = W.Text(value="", description="Chain of initial residue:", placeholder="A", layout=wide, style=DESC)

    short_len_opts = [5, 8, 10, 13, 15, 20, 25, 30]
    row_short, short_ctrl = _list_or_custom_row("Length of Paths", short_len_opts, 5, 1, 10_000, 1)
    get_short_len = short_ctrl['get']

    num_paths_opts_mode2 = [1000, 2000, 3000, 5000, 10000]
    row_np2, np2_ctrl = _list_or_custom_row("Number of Paths", num_paths_opts_mode2, 1000, 1, 10_000_000, 100)
    get_num_paths_2 = np2_ctrl['get']

    final_idx = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1,
                                 description="Index of final residue:", layout=wide, style=DESC)
    final_chain = W.Text(value="", description="Chain of final residue:", placeholder="B", layout=wide, style=DESC)

    num_paths_opts_mode3 = [1000, 2000, 3000, 5000, 10000]
    row_np3, np3_ctrl = _list_or_custom_row("Number of Paths", num_paths_opts_mode3, 1000, 1, 10_000_000, 100)
    get_num_paths_3 = np3_ctrl['get']

    btn_submit = W.Button(description="Submit", button_style="success", icon="paper-plane")
    btn_new_job = W.Button(description="New Job",  button_style="info", icon="plus")

    if _LOG_OUT is None:
        _LOG_OUT = W.Output()
    out = _LOG_OUT

    pdb_row = W.HBox([pdb_code, W.HTML("&nbsp;"), or_lbl, pdb_upload, file_lbl],
                     layout=W.Layout(align_items="center", justify_content="flex-start", flex_flow="row wrap", gap="10px"))
    functional_box = W.VBox([row_big])
    mode2_box = W.VBox([init_idx, init_chain, row_short, row_np2])
    mode3_box = W.VBox([init_idx, init_chain, final_idx, final_chain, row_np3])

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
        W.HTML("<b>Prediction of:</b>"),
        pred_type,
        functional_box,
        mode2_box,
        mode3_box,
        chain_id, email,
        W.HBox([btn_submit, btn_new_job]),
        W.HTML("<hr>"),
        out
    ]
    root = W.VBox(body_children, layout=W.Layout(width="auto"))
    _FORM_ROOT = root
    display(root)

    # -------------------- handlers --------------------
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
        if pdb_upload.value:
            (_, meta) = next(iter(pdb_upload.value.items()))
            return meta["content"], meta.get("metadata", {}).get("name", "upload.pdb")
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
        _LOG_OUT.clear_output(wait=True)
        with _LOG_OUT:
            try:
                total_steps = 5

                chain_global_raw = chain_id.value.strip()
                if chain_global_raw and (not _is_valid_chain_list(chain_global_raw)):
                    raise ValueError("Chain ID must be empty (all chains) or a comma-separated list of single characters (e.g., A or A,B).")
                chain_global = chain_global_raw  # may be ""

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
                    if not _is_valid_chain(init_chain.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    rows = ["2", pdb_name, chain_global, str(get_short_len()), str(get_num_paths_2()),
                            str(int(init_idx.value)), (init_chain.value or "").strip(), (email.value.strip() or "-")]
                else:
                    if not _is_valid_chain(init_chain.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    if not _is_valid_chain(final_chain.value or ""):
                        raise ValueError("Chain of final residue must be a single character.")
                    rows = ["3", pdb_name, chain_global,
                            str(int(init_idx.value)), (init_chain.value or "").strip(),
                            str(int(final_idx.value)), (final_chain.value or "").strip(),
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
                        cor_path = run_readpdb(pdb_path=save_path, chain=(chain_global or "A"))
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
                if mode == "functional" and cor_path and os.path.isfile(cor_path):
                    try:
                        inf_mod = importlib.import_module("mcpath.infinite")
                        importlib.reload(inf_mod)

                        work_dir = os.path.dirname(save_path)
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

                            # ---- Step 4b: write B-factor PDBs + viewer ----
                            try:
                                base = os.path.splitext(os.path.basename(save_path))[0]

                                # Build CA mapping consistent with user's chain selection:
                                # if chain_global is empty -> all chains; else only those chains
                                residues_in_order = _ca_residues_in_order(save_path, chain_str=chain_global)

                                # -------- CLOSENESS peaks: parse chain+number and auto-map (resi vs index) --------
                                close_peak_keys = set()
                                close_src = _find_first_existing(work_dir, [
                                    "closeness_peaks", "closeness_peaks.txt",
                                    "closeness_chain_labels.txt"
                                ])
                                if close_src:
                                    rows = _parse_chain_resnum_value_lines(close_src)
                                    if rows:
                                        rows.sort(key=lambda x: x[0], reverse=True)
                                        pairs = [(ch, n) for (val, ch, n) in rows[:30]]
                                    else:
                                        pairs = _parse_reschain_pairs(close_src)
                                    close_peak_keys = _pairs_to_peak_keys_auto(pairs, residues_in_order)
                                    print(f"[closeness] parsed {len(pairs)} peaks from {os.path.basename(close_src)} -> mapped {len(close_peak_keys)} residues")
                                else:
                                    print("Warning: no closeness peaks file found.")

                                close_pdb_out = os.path.join(work_dir, f"{base}_CLOSENESS_peaks_bfac.pdb")
                                if close_peak_keys:
                                    _write_bfactor_peak_pdb(save_path, close_pdb_out, close_peak_keys, peak_b=100.0, other_b=0.0)
                                    print(f"[PDB] Wrote closeness peaks B-factor PDB: {close_pdb_out}")
                                    _view_bfactor_pdb_py3dmol(close_pdb_out, "Closeness peaks (B=100 → RED spheres)", sphere_radius=0.9, show_sticks=False)
                                else:
                                    print("Warning: closeness peak set is empty; skipping closeness PDB write.")

                                # -------- BETWEENNESS peaks: parse chain+number; fallback to pure ints; auto-map --------
                                betw_peak_keys = set()
                                betw_src = _find_first_existing(work_dir, [
                                    "betweenness_peaks", "betweenness_peaks.txt",
                                    "betw_peaks", "betw_peaks.txt",
                                    "betweenness_chain_labels.txt"
                                ])
                                if betw_src:
                                    rows = _parse_chain_resnum_value_lines(betw_src)
                                    if rows:
                                        rows.sort(key=lambda x: x[0], reverse=True)
                                        pairs = [(ch, n) for (val, ch, n) in rows[:30]]
                                        betw_peak_keys = _pairs_to_peak_keys_auto(pairs, residues_in_order)
                                        print(f"[betweenness] parsed {len(pairs)} peaks from {os.path.basename(betw_src)} -> mapped {len(betw_peak_keys)} residues")
                                    else:
                                        pairs = _parse_reschain_pairs(betw_src)
                                        if pairs:
                                            betw_peak_keys = _pairs_to_peak_keys_auto(pairs, residues_in_order)
                                            print(f"[betweenness] parsed {len(pairs)} peaks from {os.path.basename(betw_src)} -> mapped {len(betw_peak_keys)} residues")
                                        else:
                                            ints_list = _read_pure_ints(betw_src)
                                            betw_peak_keys = _ints_to_peak_keys_auto(ints_list, residues_in_order)
                                            print(f"[betweenness] parsed {len(ints_list)} ints from {os.path.basename(betw_src)} -> mapped {len(betw_peak_keys)} residues")
                                else:
                                    print("Warning: no betweenness peaks file found.")

                                betw_pdb_out = os.path.join(work_dir, f"{base}_BETWEENNESS_peaks_bfac.pdb")
                                if betw_peak_keys:
                                    _write_bfactor_peak_pdb(save_path, betw_pdb_out, betw_peak_keys, peak_b=100.0, other_b=0.0)
                                    print(f"[PDB] Wrote betweenness peaks B-factor PDB: {betw_pdb_out}")
                                    _view_bfactor_pdb_py3dmol(betw_pdb_out, "Betweenness peaks (B=100 → RED spheres)", sphere_radius=0.9, show_sticks=False)
                                else:
                                    print("Warning: betweenness peak set is empty; skipping betweenness PDB write.")

                                # -------- PyMOL script if both exist --------
                                if close_peak_keys and betw_peak_keys:
                                    pml_path = os.path.join(work_dir, f"{base}_peaks_bfac_views.pml")
                                    _write_pymol_pml(pml_path, close_pdb_out, betw_pdb_out)
                                    print(f"[PyMOL] Wrote script: {pml_path}")
                                    _show_download(pml_path, label=os.path.basename(pml_path))

                            except Exception as e_bfac:
                                print(f"Warning: B-factor peak PDB generation/view failed: {e_bfac}")

                        finally:
                            os.chdir(old_cwd2)

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
                        "chain_initial": (init_chain.value or "").strip(),
                        "index_final": int(final_idx.value),
                        "chain_final": (final_chain.value or "").strip(),
                        "number_paths": int(get_num_paths_3())
                    })

                files = {FN["pdb_file"]: (pdb_name, pdb_bytes, "chemical/x-pdb")}
                if not target_url:
                    _progress(5, total_steps, "Local processing completed (no remote submission configured).")
                    return

                try:
                    r = requests.post(target_url, data=data, files=files, timeout=180)
                    _progress(5, total_steps, f"Remote submission completed (HTTP {r.status_code}).")
                except Exception as e_sub:
                    print(f"Warning: remote submission failed: {e_sub}")

            except Exception as e:
                print("Error:", e)

    btn_new_job.on_click(on_new_job)
    btn_submit.on_click(on_submit)
