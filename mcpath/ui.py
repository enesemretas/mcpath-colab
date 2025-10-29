# mcpath/ui.py
import os, re, io, requests, yaml
from IPython.display import display, clear_output
import ipywidgets as W

# Try to import the local Python parser; optional
_rps = None
try:
    # when installed as a package
    from . import readpdb_strict as _rps  # type: ignore
except Exception:
    try:
        # when running from a flat directory
        import readpdb_strict as _rps  # type: ignore
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
    """Create a centered/left/right logo if branding['logo_url'] is set."""
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
        # If the URL fails, just skip the logo silently
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

# -------------------- main UI --------------------
def launch(
    defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
    show_title="MCPath-style Parameters"
):
    """Launches the MCPath parameter form (UI only)."""
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=30).text)
    target_url = (cfg.get("target_url") or "").strip()
    FN = cfg["field_names"]

    # If present in your YAML and True, we will run the local Python readpdb_strict after submit
    run_python_readpdb = bool(cfg.get("run_python_readpdb_after_submit", False))

    # ---------- Branding / logo ----------
    logo = _logo_widget(cfg.get("branding", {}))

    # ---------- PDB: code OR upload (with "or" & filename label) ----------
    pdb_code   = W.Text(value=str(cfg.get("pdb_code", "")), description="PDB code:")
    or_lbl     = W.HTML("<b>&nbsp;&nbsp;or&nbsp;&nbsp;</b>")
    pdb_upload = W.FileUpload(accept=".pdb", multiple=False, description="Choose .pdb")
    file_lbl   = W.Label("No file chosen")

    def _on_upload_change(change):
        if _fu_has_file(pdb_upload):
            fname, _ = _fu_first(pdb_upload)
            file_lbl.value = fname or "1 file"
        else:
            file_lbl.value = "No file chosen"
    pdb_upload.observe(_on_upload_change, names="value")

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

    # ---------- Mode 2: initial residue + short path length + number of paths ----------
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

    # ---------- Mode 3: initial & final residues + number of paths ----------
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

    # ---------- Grouped layouts per section ----------
    pdb_row = W.HBox([pdb_code, W.HTML("&nbsp;"), or_lbl, pdb_upload, file_lbl])

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
        pdb_row,
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
        pdb_code.value = ""
        # reset upload widget value for both ipywidgets v8 and v7
        try:
            pdb_upload.value = ()     # v8
        except Exception:
            pdb_upload.value = {}     # v7
        file_lbl.value = "No file chosen"

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
        with out: clear_output()

    def _collect_pdb_bytes():
        if _fu_has_file(pdb_upload):
            fname, content = _fu_first(pdb_upload)
            if not content:
                raise ValueError("Uploaded file is empty or unreadable.")
            return content, fname
        code = pdb_code.value.strip()
        if not _is_valid_pdb_code(code):
            raise ValueError("PDB code must be exactly 4 alphanumeric characters (or upload a file).")
        return _fetch_rcsb(code), f"{code.upper()}.pdb"

    def _maybe_run_local_readpdb(pdb_path: str, chain: str):
        """Optionally run Python readpdb_strict (if available) and preview .cor head."""
        if not run_python_readpdb:
            return
        if _rps is None:
            print("ℹ️ readpdb_strict.py not found in environment; skipping local parse.")
            return

        try:
            # discover callable
            runner = None
            for name in ("run", "run_readpdb", "readpdb"):
                if hasattr(_rps, name):
                    runner = getattr(_rps, name)
                    break

            if runner is None:
                print("ℹ️ readpdb_strict has no callable (run/run_readpdb/readpdb); skipping.")
                return

            print("▶ Running Python readpdb_strict…")
            # Most implementations accept (pdb_path, chain_id, out_basename=None, strict_matlab_mode=True/False)
            try:
                runner(pdb_path, chain, out_basename=os.path.basename(pdb_path), strict_matlab_mode=True)
            except TypeError:
                # Fallback: maybe only (pdb_path, chain)
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

    def on_submit(_):
        with out:
            clear_output()
            try:
                # basic validation
                chain_global = chain_id.value.strip()
                if not _is_valid_chain(chain_global):
                    raise ValueError("Chain ID must be a single character (e.g., A).")
                if not _is_valid_email(email.value.strip()):
                    raise ValueError("Invalid email format.")

                pdb_bytes, pdb_name = _collect_pdb_bytes()

                # build payload
                data = {
                    "prediction_mode": pred_type.value,   # for your backend/router
                    FN["chain_id"]: chain_global,
                }
                if pdb_code.value.strip():
                    data[FN["pdb_code"]] = pdb_code.value.strip().upper()
                if email.value.strip():
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

                # Save a local copy of the PDB for the user (and for local parser)
                with open(pdb_name, "wb") as f:
                    f.write(pdb_bytes)
                print(f"Saved local copy: {os.getcwd()}/{pdb_name}")

                # (Optional) run Python readpdb_strict locally and preview .cor
                _maybe_run_local_readpdb(os.path.abspath(pdb_name), chain_global)

                if not target_url:
                    # preview only
                    print("\n(No target_url set) — preview only payload below:\n")
                    preview = dict(data)
                    preview["attached_file"] = pdb_name
                    print(preview)
                    return

                # real POST
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
