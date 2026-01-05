# mcpath/ui.py
import os, re, requests, yaml, importlib, shutil, sys, glob
from IPython.display import display, clear_output
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
        return True  # empty -> all chains
    tokens = [t.strip() for t in re.split(r"[,\s]+", ch_str) if t.strip()]
    if not tokens:
        return True
    return all(_is_valid_chain(tok) for tok in tokens)

def _is_valid_email(s):    
    return (not s.strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s.strip()))

def _fetch_rcsb(code: str) -> bytes:
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    r = requests.get(url, timeout=60); r.raise_for_status()
    return r.content

def _logo_widget(branding: dict):
    url = (branding.get("logo_url") or "").strip()
    if not url: return None
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
    """
    Minimal progress message helper.
    Example: [1/5] Input validated and PDB saved.
    """
    print(f"[{step}/{total}] {message}")

# -------------------- Peak -> B-factor PDB + viewer helpers --------------------
_NUM_RE = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")

def _ensure_py3dmol():
    try:
        import py3Dmol  # noqa: F401
        return True
    except Exception:
        try:
            import subprocess, sys
            subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "py3Dmol"])
            import py3Dmol  # noqa: F401
            return True
        except Exception as e:
            print(f"Warning: py3Dmol not available and install failed: {e}")
            return False

def _parse_chain_tokens(chain_str: str):
    chain_str = (chain_str or "").strip()
    if not chain_str:
        return None  # all chains
    toks = [t.strip() for t in re.split(r"[,\s]+", chain_str) if t.strip()]
    return set(toks) if toks else None

def _ca_residues_in_order_keys(pdb_path: str, chain_str: str = ""):
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

def _read_all_ints(path: str):
    ints = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            for n in _NUM_RE.findall(line):
                try:
                    ints.append(int(float(n)))
                except Exception:
                    pass
    return ints

def _parse_closeness_chain_labels_top_resnums(labels_path: str, top_n: int = 30):
    """
    From closeness_chain_labels.txt choose top_n residues by (last float) value.
    Returns: list of (chain, resi_num)
    """
    rows = []
    with open(labels_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue

            toks = re.split(r"\s+", s.replace(":", " "))
            # find chain token (single alnum char)
            chain = None
            chain_pos = None
            for i, t in enumerate(toks):
                if len(t) == 1 and t.isalnum():
                    chain = t
                    chain_pos = i
                    break
            if chain is None:
                continue

            # numeric tokens
            nums = _NUM_RE.findall(s)
            if not nums:
                continue

            # value = last float
            try:
                val = float(nums[-1])
            except Exception:
                continue

            # try to pick residue number = first int token AFTER the chain token, else last int in line
            resi_num = None
            for t in toks[(chain_pos + 1):]:
                if re.fullmatch(r"\d+", t):
                    resi_num = int(t)
                    break
            if resi_num is None:
                # fallback: last integer-like among nums
                int_nums = []
                for n in nums:
                    try:
                        int_nums.append(int(float(n)))
                    except Exception:
                        pass
                if int_nums:
                    # often [index, resi] -> pick second if exists
                    resi_num = int_nums[1] if len(int_nums) >= 2 else int_nums[-1]

            if resi_num is None:
                continue

            rows.append((val, chain, int(resi_num)))

    rows.sort(key=lambda x: x[0], reverse=True)
    rows = rows[:max(1, int(top_n))]
    return [(ch, rn) for _, ch, rn in rows]

def _peaks_to_keys_from_indices_or_resnums(peaks_ints, residues_in_order):
    """
    peaks_ints could be:
      - 1-based indices into residues_in_order
      - or residue numbers
    Returns set of keys: {(chain, resi_raw), ...}
    """
    keys = set()
    if not peaks_ints:
        return keys

    N = len(residues_in_order)
    if N <= 0:
        return keys

    # Heuristic: if most values are <= N, treat as indices
    leN = sum(1 for x in peaks_ints if 1 <= x <= N)
    if leN >= max(1, int(0.7 * len(peaks_ints))):
        # indices
        for idx in peaks_ints:
            if 1 <= idx <= N:
                r = residues_in_order[idx - 1]
                keys.add((r["chain"], r["resi_raw"]))
        return keys

    # else treat as residue numbers: highlight all residues that match those resi_num (across chains)
    want_nums = set(int(x) for x in peaks_ints)
    for r in residues_in_order:
        if r["resi_num"] in want_nums:
            keys.add((r["chain"], r["resi_raw"]))
    return keys

def _keys_from_chain_resnums(chain_resnums, residues_in_order):
    """
    chain_resnums: list[(chain, resi_num)]
    Returns keys matching those (chain, resi_num) in the PDB (handles insertion codes by matching resi_num).
    """
    want = {}
    for ch, rn in chain_resnums:
        want.setdefault(ch, set()).add(int(rn))

    keys = set()
    for r in residues_in_order:
        if r["chain"] in want and r["resi_num"] in want[r["chain"]]:
            keys.add((r["chain"], r["resi_raw"]))
    return keys

def _write_bfactor_peak_pdb(pdb_in: str, pdb_out: str, peak_keys: set,
                            peak_b: float = 100.0, other_b: float = 0.0):
    """
    Writes a PDB where B-factor is set to peak_b for residues in peak_keys, else other_b.
    peak_keys contains (chain, resi_raw) pairs.
    """
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

                # ensure line is long enough
                if len(line) < 66:
                    line = line.rstrip("\n").ljust(66) + "\n"

                # B-factor columns 61-66 => [60:66]
                line = line[:60] + b_str + line[66:]
            fout.write(line)

    return pdb_out

def _extract_peak_resis_from_bfactor_pdb(pdb_path: str, b_threshold: float = 50.0):
    """
    Reads CA atoms and collects residues whose B >= threshold.
    Returns dict: {chain: [resi_num, ...]}
    """
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

def _view_bfactor_pdb_py3dmol(pdb_path: str, title: str):
    """
    Viewer that highlights residues RED based on B-factor >= 50 in the given PDB.
    """
    if not _ensure_py3dmol():
        return

    import py3Dmol
    from IPython.display import display, HTML, FileLink

    pdb_text = open(pdb_path, "r", encoding="utf-8", errors="ignore").read()
    peaks = _extract_peak_resis_from_bfactor_pdb(pdb_path, b_threshold=50.0)

    view = py3Dmol.view(width=850, height=520)
    view.addModel(pdb_text, "pdb")
    view.setStyle({}, {"cartoon": {"color": "lightgray"}})

    # Highlight based on B-factor-coded peaks
    for ch, resis in peaks.items():
        if not resis:
            continue
        view.addStyle({"chain": ch, "resi": resis}, {"stick": {"color": "red"}})
        view.addStyle({"chain": ch, "resi": resis, "atom": "CA"}, {"sphere": {"color": "red", "radius": 0.75}})

    view.zoomTo()

    # Save HTML always (reliable inside widget Output)
    html_path = os.path.splitext(pdb_path)[0] + ".html"
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(view._make_html())

    display(HTML(f"<b>{title}</b>"))
    display(FileLink(html_path))
    view.show()

def _write_pymol_pml(pml_path: str, pdb_close: str, pdb_betw: str):
    """
    PyMOL script to color peaks red using B-factors (b > 50) for each object.
    """
    pml = f"""
reinitialize
bg_color white
set ray_opaque_background, off

load {pdb_close}, closeness
load {pdb_betw}, betweenness

disable betweenness
hide everything, all

# --- closeness view ---
enable closeness
show cartoon, closeness
color gray80, closeness
select c_peaks, closeness and b > 50
color red, c_peaks
show sticks, c_peaks
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
show sticks, b_peaks
show spheres, b_peaks and name CA
set sphere_scale, 0.6, b_peaks
orient betweenness
scene betweenness, store

# Tip: Scene menu -> recall closeness / betweenness
"""
    os.makedir






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

    # --- Configuration loaded here and used for defaults/resetting ---
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=30).text)
    target_url = cfg.get("target_url", "").strip()
    FN = cfg["field_names"]

    SAVE_DIR = cfg.get("save_dir", "/content")
    os.makedirs(SAVE_DIR, exist_ok=True)

    logo = _logo_widget(cfg.get("branding", {}))

    DESC = {'description_width': '180px'}
    wide = W.Layout(width="420px", min_width="360px")

    # ---------- PDB: code OR upload ----------
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

    # Chain ID: now supports comma-separated list OR empty for "all chains"
    chain_id = W.Text(
        value=str(cfg.get("chain_id", "")),
        description="Chain ID:",
        placeholder="A or A,B (empty = all chains)",
        layout=wide, style=DESC
    )
    email      = W.Text(value=str(cfg.get("email", "")),   description="Email (opt):", layout=wide, style=DESC)

    # Prediction type
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

    # ---------- Functional mode: big path_length ----------
    big_opts = cfg.get("path_length_options", [100000, 200000, 300000, 400000, 500000, 750000, 1000000, 2000000, 3000000, 4000000])
    big_default = int(cfg.get("path_length", big_opts[0] if big_opts else 100000))
    row_big, big_ctrl = _list_or_custom_row("Path length", big_opts, big_default, 1, 10_000_000, 1000)
    get_big_len = big_ctrl['get']

    # ---------- Mode 2 ----------
    init_idx_default = 1
    init_idx   = W.BoundedIntText(value=init_idx_default, min=1, max=1_000_000, step=1,
                                  description="Index of initial residue:", layout=wide, style=DESC)
    init_chain_default = ""
    init_chain = W.Text(value=init_chain_default, description="Chain of initial residue:", placeholder="A", layout=wide, style=DESC)

    short_len_opts = [5, 8, 10, 13, 15, 20, 25, 30]
    short_len_default = 5
    row_short, short_ctrl = _list_or_custom_row("Length of Paths", short_len_opts, short_len_default, 1, 10_000, 1)
    get_short_len = short_ctrl['get']

    num_paths_opts_mode2 = [1000, 2000, 3000, 5000, 10000, 20000, 30000, 40000, 50000]
    num_paths_mode2_default = 1000
    row_np2, np2_ctrl = _list_or_custom_row("Number of Paths", num_paths_opts_mode2, num_paths_mode2_default, 1, 10_000_000, 100)
    get_num_paths_2 = np2_ctrl['get']

    # ---------- Mode 3 ----------
    final_idx_default = 1
    final_idx   = W.BoundedIntText(value=final_idx_default, min=1, max=1_000_000, step=1,
                                   description="Index of final residue:", layout=wide, style=DESC)
    final_chain_default = ""
    final_chain = W.Text(value=final_chain_default, description="Chain of final residue:", placeholder="B", layout=wide, style=DESC)

    num_paths_opts_mode3 = [1000, 2000, 3000, 5000, 10000, 30000, 50000]
    num_paths_mode3_default = 1000
    row_np3, np3_ctrl = _list_or_custom_row("Number of Paths", num_paths_opts_mode3, num_paths_mode3_default, 1, 10_000_000, 100)
    get_num_paths_3 = np3_ctrl['get']

    # Actions & shared output (singleton)
    btn_submit = W.Button(description="Submit", button_style="success", icon="paper-plane")
    btn_new_job = W.Button(description="New Job",  button_style="info", icon="plus")

    if _LOG_OUT is None:
        _LOG_OUT = W.Output()
    out = _LOG_OUT  # reuse the same Output every time

    # Grouped layouts
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

    # Root container (singleton root so UI re-launch doesn’t multiply outputs)
    body_children = []
    if logo: body_children.append(logo)
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
    _FORM_ROOT = root  # track latest

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
            return meta["content"], meta["metadata"]["name"] if "metadata" in meta and "name" in meta["metadata"] else "upload.pdb"
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

    # helper to copy with auto-suffix if destination exists
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

                # -------- Step 1: validate inputs + save PDB --------
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

                # -------- Step 2: generate .cor via readpdb_strict --------
                run_readpdb = _try_import_readpdb()
                cor_path = None
                if run_readpdb is None:
                    print("Warning: readpdb_strict not found; skipping .cor generation.")
                else:
                    try:
                        cor_path = run_readpdb(input_path=input_path)
                    except TypeError:
                        # Fallback to direct mode (rare)
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

                        coor_fixed = _copy_unique(cor_src,  "coor_file", work_dir=os.path.dirname(cor_src))
                        atom_fixed = _copy_unique(atom_src, "atom_file", work_dir=os.path.dirname(atom_src))
                        path_fixed = _copy_unique(path_src, "path_file", work_dir=os.path.dirname(path_src))

                        close_mod = importlib.import_module("mcpath.closeness")
                        importlib.reload(close_mod)

                        betw_mod = importlib.import_module("mcpath.betweenness")
                        importlib.reload(betw_mod)

                        old_cwd2 = os.getcwd()
                        try:
                            os.chdir(work_dir)
                            close_mod.main()
                            betw_mod.main()

                            # ---- Step 4b: write B-factor PDBs (peaks=100) + show via viewer ----
                            try:
                                work_dir = os.path.dirname(save_path)
                                base = os.path.splitext(os.path.basename(save_path))[0]

                                # CA order mapping (used for indices->residues)
                                residues_in_order = _ca_residues_in_order_keys(save_path, chain_global)

                                # -------- Closeness peaks: from closeness_chain_labels.txt (top-N by value) --------
                                close_labels = os.path.join(work_dir, "closeness_chain_labels.txt")
                                close_peak_keys = set()
                                if os.path.isfile(close_labels):
                                    top_chain_resnums = _parse_closeness_chain_labels_top_resnums(close_labels, top_n=30)
                                    close_peak_keys = _keys_from_chain_resnums(top_chain_resnums, residues_in_order)
                                else:
                                    print("Warning: closeness_chain_labels.txt not found; cannot build closeness peak set.")

                                close_pdb_out = os.path.join(work_dir, f"{base}_CLOSENESS_peaks_bfac.pdb")
                                if close_peak_keys:
                                    _write_bfactor_peak_pdb(save_path, close_pdb_out, close_peak_keys, peak_b=100.0, other_b=0.0)
                                    print(f"[PDB] Wrote closeness peaks B-factor PDB: {close_pdb_out}")
                                else:
                                    print("Warning: closeness peak set is empty; skipping closeness PDB write.")

                                # -------- Betweenness peaks: from betweenness_peaks (indices or residue numbers) --------
                                betw_peaks_path = None
                                for cand in ["betweenness_peaks", "betweenness_peaks.txt", "betw_peaks", "betw_peaks.txt"]:
                                    p = os.path.join(work_dir, cand)
                                    if os.path.isfile(p):
                                        betw_peaks_path = p
                                        break

                                betw_peak_keys = set()
                                if betw_peaks_path:
                                    betw_ints = _read_all_ints(betw_peaks_path)
                                    betw_peak_keys = _peaks_to_keys_from_indices_or_resnums(betw_ints, residues_in_order)
                                else:
                                    print("Warning: betweenness_peaks file not found; cannot build betweenness peak set.")

                                betw_pdb_out = os.path.join(work_dir, f"{base}_BETWEENNESS_peaks_bfac.pdb")
                                if betw_peak_keys:
                                    _write_bfactor_peak_pdb(save_path, betw_pdb_out, betw_peak_keys, peak_b=100.0, other_b=0.0)
                                    print(f"[PDB] Wrote betweenness peaks B-factor PDB: {betw_pdb_out}")
                                else:
                                    print("Warning: betweenness peak set is empty; skipping betweenness PDB write.")

                                # -------- PyMOL script (real PyMOL) --------
                                if close_peak_keys and betw_peak_keys:
                                    pml_path = os.path.join(work_dir, f"{base}_peaks_bfac_views.pml")
                                    _write_pymol_pml(pml_path, close_pdb_out, betw_pdb_out)
                                    print(f"[PyMOL] Wrote script: {pml_path}")
                                    print(f"[PyMOL] Run: pymol {pml_path}")

                                # -------- Show in notebook (viewer reads B-factors and paints RED) --------
                                if os.path.isfile(close_pdb_out):
                                    _view_bfactor_pdb_py3dmol(close_pdb_out, "Closeness peaks (B=100 -> RED)")

                                if os.path.isfile(betw_pdb_out):
                                    _view_bfactor_pdb_py3dmol(betw_pdb_out, "Betweenness peaks (B=100 -> RED)")

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
