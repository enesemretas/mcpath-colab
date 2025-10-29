# mcpath/ui.py
import os, re, io, requests, yaml
from IPython.display import display, clear_output
import ipywidgets as W

# Optional: if you have readpdb_strict.py, we can attempt to import it
_rps = None
try:
    from . import readpdb_strict as _rps  # package mode
except Exception:
    try:
        import readpdb_strict as _rps     # flat mode
    except Exception:
        _rps = None

# -------------------- validators & helpers --------------------
def _is_valid_pdb_code(c): return bool(re.fullmatch(r"[0-9A-Za-z]{4}", (c or "").strip()))
def _is_valid_chain(ch):   return bool(re.fullmatch(r"[A-Za-z0-9]", (ch or "").strip()))
def _is_valid_email(s):    return (not (s or "").strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", (s or "").strip()))

def _fetch_rcsb(code: str) -> bytes:
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    r = requests.get(url, timeout=60); r.raise_for_status()
    return r.content

def _list_or_custom(label: str, options, default_value, minv, maxv=None, step=1, **kwargs):
    """Factory for a (List|Custom) toggle + Dropdown + BoundedIntText trio."""
    # Back-compat: allow callers to pass max=...
    if maxv is None and "max" in kwargs:
        maxv = kwargs["max"]

    options = sorted({int(x) for x in options})
    default_value = int(default_value)
    if default_value not in options:
        options = sorted(options + [default_value])

    mode = W.ToggleButtons(options=[("List","list"), ("Custom","custom")], value="list", description="Input:")
    dd   = W.Dropdown(options=options, value=default_value, description=label)
    txt  = W.BoundedIntText(value=default_value, min=int(minv), max=int(maxv or minv), step=int(step), description=label)

    def _sync(*_):
        if mode.value == "list":
            dd.layout.display = ""
            txt.layout.display = "none"
        else:
            dd.layout.display = "none"
            txt.layout.display = ""
    _sync()
    mode.observe(_sync, names="value")

    def get_value():
        return int(dd.value if mode.value == "list" else txt.value)

    return mode, dd, txt, get_value

def _logo_widget(branding: dict):
    url = (branding.get("logo_url") or "").strip()
    if not url:
        return None
    height = int(branding.get("logo_height", 96))
    align  = (branding.get("logo_align") or "center").lower()
    try:
        img_bytes = requests.get(url, timeout=20).content
        img = W.Image(value=img_bytes, format="png", layout=W.Layout(height=f"{height}px"))
        jc = {"center":"center", "left":"flex-start", "right":"flex-end"}.get(align, "center")
        box = W.HBox([img], layout=W.Layout(justify_content=jc))
        return box
    except Exception:
        return None

# -------------------- FileUpload (v7 & v8 compatible) --------------------
def _fu_has_file(upload_widget: W.FileUpload) -> bool:
    v = getattr(upload_widget, "value", None)
    if v is None:
        return False
    if isinstance(v, dict):
        return len(v) > 0                       # ipywidgets v7
    if isinstance(v, (tuple, list)):
        return len(v) > 0                       # ipywidgets v8
    return False

def _fu_first(upload_widget: W.FileUpload):
    """
    Returns (filename:str|None, content:bytes|None) for the first uploaded file.
    Supports ipywidgets v7 (dict) and v8 (tuple/list of UploadedFile).
    """
    v = getattr(upload_widget, "value", None)

    # v7: dict -> {"name": {"content": b"...", "metadata": ...}, ...}
    if isinstance(v, dict) and v:
        fname = next(iter(v.keys()))
        meta = v[fname]
        content = meta.get("content") if isinstance(meta, dict) else getattr(meta, "content", None)
        return fname, content

    # v8: tuple/list of UploadedFile objects (or dict-like)
    if isinstance(v, (tuple, list)) and v:
        uf = v[0]
        fname = getattr(uf, "name", None) or (uf.get("name") if isinstance(uf, dict) else None)
        content = getattr(uf, "content", None) or (uf.get("content") if isinstance(uf, dict) else None)
        return fname, content

    return None, None

def _make_fileupload(on_change_cb):
    """Create a fresh FileUpload widget + filename label, wired to on_change_cb."""
    fu = W.FileUpload(
        accept=".pdb,.ent,.PDB,.ENT,.txt,.gz",  # a bit looser for Colab quirks
        multiple=False,
        description="Choose file"
    )
    lbl = W.Label("No file chosen")

    def _on_change(change):
        if _fu_has_file(fu):
            fname, _ = _fu_first(fu)
            lbl.value = fname or "1 file"
        else:
            lbl.value = "No file chosen"
        if callable(on_change_cb):
            on_change_cb(change)

    fu.observe(_on_change, names="value")
    return fu, lbl

# -------------------- Optional local Python readpdb runner --------------------
def _maybe_run_local_readpdb(pdb_path: str, chain: str, enable: bool):
    if not enable:
        return
    if _rps is None:
        print("ℹ️ readpdb_strict.py not found in environment; skipping local parse.")
        return
    try:
        runner = None
        for name in ("run", "run_readpdb", "readpdb"):
            if hasattr(_rps, name):
                runner = getattr(_rps, name)
                break
        if runner is None:
            print("ℹ️ readpdb_strict has no callable (run/run_readpdb/readpdb); skipping.")
            return

        print("▶ Running Python readpdb_strict…")
        try:
            runner(pdb_path, chain, out_basename=os.path.basename(pdb_path), strict_matlab_mode=True)
        except TypeError:
            runner(pdb_path, chain)

        cor_file = pdb_path + ".cor"
        if os.path.exists(cor_file):
            with open(cor_file, "r") as f:
                head = "".join([next(f) for _ in range(8) if True])
            print(f"✅ Wrote: {cor_file}\n\n— .cor preview (first lines) —\n{head}")
        else:
            print("⚠️ No .cor file produced by readpdb_strict.")
    except Exception as e:
        print("⚠️ readpdb_strict failed:", e)

# -------------------- main UI --------------------
def launch(
    defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
    show_title="MCPath-style Parameters"
):
    """Launches the MCPath parameter form (UI only)."""
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=30).text)
    target_url = (cfg.get("target_url") or "").strip()
    FN = cfg["field_names"]
    run_python_readpdb = bool(cfg.get("run_python_readpdb_after_submit", False))

    # ---------- Branding / logo ----------
    logo = _logo_widget(cfg.get("branding", {}))

    # ---------- PDB: code OR upload ----------
    pdb_code = W.Text(value=str(cfg.get("pdb_code", "")), description="PDB code:")

    # build upload row with a replaceable upload widget
    upload_container = W.HBox()  # we will set its children below
    # create the first upload
    pdb_upload, file_lbl = _make_fileupload(on_change_cb=None)
    or_lbl = W.HTML("<b>&nbsp;&nbsp;or&nbsp;&nbsp;</b>")
    upload_container.children = (pdb_code, W.HTML("&nbsp;"), or_lbl, pdb_upload, file_lbl)

    # ---------- Always-present fields ----------
    chain_id   = W.Text(value=str(cfg.get("chain_id", "")), description="Chain ID:")
    email      = W.Text(value=str(cfg.get("email", "")), description="Email (opt):")

    # ---------- Prediction type ----------
    pred_type = W.RadioButtons(
        options=[
            ("Functional Residues", "functional"),
            ("Allosteric Paths (initial residue + path length)", "paths_init_len"),
            ("Allosteric Paths (initial & final residues)", "paths_init_final"),
        ],
        value="functional",
        description="Prediction:",
        layout=W.Layout(width="auto")
    )

    # ---------- Functional mode: big path_length (YAML-driven List/Custom) ----------
    big_opts = cfg.get("path_length_options", [
        100000, 200000, 300000, 400000, 500000, 750000,
        1000000, 2000000, 3000000, 4000000
    ])
    big_default = int(cfg.get("path_length", big_opts[0] if big_opts else 100000))
    (pl_mode_big, pl_dd_big, pl_txt_big, get_big_len) = _list_or_custom(
        label="Path length:", options=big_opts, default_value=big_default,
        minv=1, maxv=10_000_000, step=1000
    )

    # ---------- Mode 2 ----------
    init_idx   = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1, description="Index of initial residue")
    init_chain = W.Text(value="", description="Chain of initial residue", placeholder="A", layout=W.Layout(width="260px"))
    short_len_opts = [5, 8, 10, 13, 15, 20, 25, 30]
    (pl_mode_short, pl_dd_short, pl_txt_short, get_short_len) = _list_or_custom(
        label="Length of Paths:", options=short_len_opts, default_value=5,
        minv=1, maxv=10_000, step=1
    )
    num_paths_opts_mode2 = [1000, 2000, 3000, 5000, 10000, 20000, 30000, 40000, 50000]
    (np_mode_2, np_dd_2, np_txt_2, get_num_paths_2) = _list_or_custom(
        label="Number of Paths:", options=num_paths_opts_mode2, default_value=1000,
        minv=1, maxv=10_000_000, step=100
    )

    # ---------- Mode 3 ----------
    final_idx   = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1, description="Index of final residue")
    final_chain = W.Text(value="", description="Chain of final residue", placeholder="B", layout=W.Layout(width="260px"))
    num_paths_opts_mode3 = [1000, 2000, 3000, 5000, 10000, 30000, 50000]
    (np_mode_3, np_dd_3, np_txt_3, get_num_paths_3) = _list_or_custom(
        label="Number of Paths:", options=num_paths_opts_mode3, default_value=1000,
        minv=1, maxv=10_000_000, step=100
    )

    # ---------- Actions & output ----------
    btn_submit = W.Button(description="Submit", button_style="success", icon="paper-plane")
    btn_clear  = W.Button(description="Clear",  button_style="warning", icon="trash")
    out        = W.Output()

    functional_box = W.VBox([pl_mode_big, pl_dd_big, pl_txt_big])
    mode2_box = W.VBox([
        init_idx, init_chain,
        pl_mode_short, pl_dd_short, pl_txt_short,
        np_mode_2, np_dd_2, np_txt_2
    ])
    mode3_box = W.VBox([
        init_idx, init_chain,
        final_idx, final_chain,
        np_mode_3, np_dd_3, np_txt_3
    ])

    def _sync_mode(*_):
        functional_box.layout.display = ""
        mode2_box.layout.display = "none"
        mode3_box.layout.display = "none"
        if pred_type.value == "paths_init_len":
            functional_box.layout.display = "none"
            mode2_box.layout.display = ""
        elif pred_type.value == "paths_init_final":
            functional_box.layout.display = "none"
            mode2_box.layout.display = "none"
            mode3_box.layout.display = ""
    _sync_mode()
    pred_type.observe(_sync_mode, names="value")

    # ---------- Render the full UI ----------
    children = []
    if logo: children.append(logo)
    children += [
        W.HTML(f"<h3>{show_title}</h3>"),
        upload_container,
        W.HTML("<b>Prediction of:</b>"),
        pred_type,
        functional_box,
        mode2_box,
        mode3_box,
        chain_id, email,
        W.HBox([btn_submit, btn_clear]),
        W.HTML("<hr>"), out
    ]
    display(W.VBox(children))

    # -------------------- handlers --------------------
    def on_clear(_):
        # Reset text fields
        pdb_code.value = ""
        chain_id.value = ""
        email.value = ""
        pred_type.value = "functional"
        # reset mode widgets
        pl_mode_big.value = "list"; pl_dd_big.value = pl_dd_big.options[0]; pl_txt_big.value = int(pl_dd_big.options[0])
        init_idx.value = 1; init_chain.value = ""
        pl_mode_short.value = "list"; pl_dd_short.value = pl_dd_short.options[0]; pl_txt_short.value = int(pl_dd_short.options[0])
        np_mode_2.value = "list"; np_dd_2.value = np_dd_2.options[0]; np_txt_2.value = int(np_dd_2.options[0])
        final_idx.value = 1; final_chain.value = ""
        np_mode_3.value = "list"; np_dd_3.value = np_dd_3.options[0]; np_txt_3.value = int(np_dd_3.options[0])

        # Recreate FileUpload (v8 value is read-only)
        nonlocal pdb_upload, file_lbl
        new_upload, new_label = _make_fileupload(on_change_cb=None)
        pdb_upload = new_upload
        file_lbl = new_label
        upload_container.children = (pdb_code, W.HTML("&nbsp;"), or_lbl, pdb_upload, file_lbl)

        with out: clear_output()

    def _collect_pdb_bytes():
        # Prefer uploaded file if present
        if _fu_has_file(pdb_upload):
            fname, content = _fu_first(pdb_upload)
            if not content:
                raise ValueError("Uploaded file is empty or unreadable.")
            return content, fname
        # fallback to PDB code
        code = (pdb_code.value or "").strip()
        if not _is_valid_pdb_code(code):
            raise ValueError("PDB code must be exactly 4 alphanumeric characters (or upload a file).")
        return _fetch_rcsb(code), f"{code.upper()}.pdb"

    def on_submit(_):
        with out:
            clear_output()
            try:
                # Validate chain + email format early (but not code yet)
                chain_global = (chain_id.value or "").strip()
                if not _is_valid_chain(chain_global):
                    raise ValueError("Chain ID must be a single character (e.g., A).")
                if not _is_valid_email(email.value or ""):
                    raise ValueError("Invalid email format.")

                pdb_bytes, pdb_name = _collect_pdb_bytes()  # this checks code only if upload empty

                # build payload
                data = {
                    "prediction_mode": pred_type.value,
                    FN["chain_id"]: chain_global,
                }
                if (pdb_code.value or "").strip():
                    data[FN["pdb_code"]] = pdb_code.value.strip().upper()
                if (email.value or "").strip():
                    data[FN["email"]] = email.value.strip()

                files = {FN["pdb_file"]: (pdb_name, pdb_bytes, "chemical/x-pdb")}

                # mode-specific fields
                if pred_type.value == "functional":
                    data[FN["path_length"]] = str(get_big_len())
                elif pred_type.value == "paths_init_len":
                    if not _is_valid_chain(init_chain.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    data.update({
                        "index_initial": int(init_idx.value),
                        "chain_initial": (init_chain.value or "").strip(),
                        "length_paths":  int(get_short_len()),
                        "number_paths":  int(get_num_paths_2()),
                    })
                elif pred_type.value == "paths_init_final":
                    if not _is_valid_chain(init_chain.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    if not _is_valid_chain(final_chain.value or ""):
                        raise ValueError("Chain of final residue must be a single character.")
                    data.update({
                        "index_initial": int(init_idx.value),
                        "chain_initial": (init_chain.value or "").strip(),
                        "index_final":   int(final_idx.value),
                        "chain_final":   (final_chain.value or "").strip(),
                        "number_paths":  int(get_num_paths_3()),
                    })

                # Save a local copy of the PDB for local parsing if needed
                with open(pdb_name, "wb") as f:
                    f.write(pdb_bytes)
                print(f"Saved local copy: {os.getcwd()}/{pdb_name}")

                # (Optional) run Python readpdb_strict locally and preview .cor
                _maybe_run_local_readpdb(os.path.abspath(pdb_name), chain_global, enable=run_python_readpdb)

                if not target_url:
                    # preview only
                    print("\n(No target_url set) — preview only payload below:\n")
                    preview = dict(data)
                    preview["attached_file"] = pdb_name
                    print(preview)
                    return

                # POST to server
                print(f"Submitting to {target_url} …")
                r = requests.post(target_url, data=data, files=files, timeout=180)
                print("HTTP", r.status_code)
                try:
                    print("JSON:", r.json())
                except Exception:
                    print("Response (≤800 chars):\n", r.text[:800])

            except Exception as e:
                print("❌", e)

    btn_clear.on_click(on_clear)
    btn_submit.on_click(on_submit)
