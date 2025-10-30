# mcpath/ui.py
import os, re, requests, yaml, importlib, shutil, sys
from IPython.display import display
import ipywidgets as W

# -------------------- validators & helpers --------------------
def _is_valid_pdb_code(c): return bool(re.fullmatch(r"[0-9A-Za-z]{4}", c.strip()))
def _is_valid_chain(ch):   return bool(re.fullmatch(r"[A-Za-z0-9]", ch.strip()))
def _is_valid_email(s):    return (not s.strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s.strip()))

def _fetch_rcsb(code: str) -> bytes:
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    r = requests.get(url, timeout=60); r.raise_for_status()
    return r.content

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
        return W.HBox([img], layout=W.Layout(justify_content=jc))
    except Exception:
        return None

# --- Inline 'List | Custom | Value' row ---
def _list_or_custom(label: str, options, default_value, minv, maxv, step=1, desc_style=None):
    options = sorted({int(x) for x in options})
    default_value = int(default_value)
    if default_value not in options:
        options = sorted(options + [default_value])

    btns = W.ToggleButtons(
        options=[("List", "list"), ("Custom", "custom")],
        value="list",
        description="",
        layout=W.Layout(width="180px", margin="0 12px 0 0"),
        style=desc_style
    )
    dd = W.Dropdown(
        options=options, value=default_value,
        description=(label if label.endswith(":") else f"{label}:"),
        layout=W.Layout(width="420px", min_width="360px"),
        style=desc_style
    )
    itxt = W.BoundedIntText(
        value=default_value, min=minv, max=maxv, step=step,
        description=(label if label.endswith(":") else f"{label}:"),
        layout=W.Layout(width="420px", min_width="360px"),
        style=desc_style
    )
    value_box = W.Box([dd], layout=W.Layout(align_items="center"))

    def _on_mode(ch):
        value_box.children = [dd] if ch["new"] == "list" else [itxt]
    btns.observe(_on_mode, names="value")

    def get_value():
        return int(dd.value if btns.value == "list" else itxt.value)

    row = W.HBox([btns, value_box], layout=W.Layout(align_items="center", gap="12px"))
    return row, btns, dd, itxt, get_value

# -------------------- main UI --------------------
def launch(
    defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
    show_title="MCPath-style Parameters"
):
    """Launches the MCPath parameter form and runs readpdb_strict + atomistic."""
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=30).text)
    target_url = cfg.get("target_url", "").strip()
    FN = cfg["field_names"]

    # Where to save files (PDB + mcpath_input.txt)
    SAVE_DIR = cfg.get("save_dir", "/content")
    os.makedirs(SAVE_DIR, exist_ok=True)

    # Branding / logo
    logo = _logo_widget(cfg.get("branding", {}))

    DESC = {'description_width': '240px'}
    wide = W.Layout(width="420px", min_width="360px")

    # ---------- PDB: code OR upload ----------
    pdb_code   = W.Text(value=str(cfg.get("pdb_code", "")),
                        description="PDB code:", layout=wide, style=DESC)
    or_lbl     = W.HTML("<b>&nbsp;&nbsp;or&nbsp;&nbsp;</b>")
    pdb_upload = W.FileUpload(accept=".pdb", multiple=False, description="Choose file")
    file_lbl   = W.Label("No file chosen")

    def _on_upload_change(change):
        if pdb_upload.value:
            fname = next(iter(pdb_upload.value.keys()))
            file_lbl.value = fname
        else:
            file_lbl.value = "No file chosen"
    pdb_upload.observe(_on_upload_change, names="value")

    # Always-present fields
    chain_id   = W.Text(value=str(cfg.get("chain_id", "")),
                        description="Chain ID:", layout=wide, style=DESC)
    email      = W.Text(value=str(cfg.get("email", "")),
                        description="Email (opt):", layout=wide, style=DESC)

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
    big_opts = cfg.get("path_length_options", [
        100000, 200000, 300000, 400000, 500000, 750000,
        1000000, 2000000, 3000000, 4000000
    ])
    big_default = int(cfg.get("path_length", big_opts[0] if big_opts else 100000))
    row_big, pl_mode_big, pl_dd_big, pl_txt_big, get_big_len = _list_or_custom(
        label="Path length", options=big_opts, default_value=big_default,
        minv=1, maxv=10_000_000, step=1000, desc_style=DESC
    )

    # ---------- Mode 2 ----------
    init_idx   = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1,
                                  description="Index of initial residue:",
                                  layout=wide, style=DESC)
    init_chain = W.Text(value="", description="Chain of initial residue:",
                        placeholder="A", layout=wide, style=DESC)

    short_len_opts = [5, 8, 10, 13, 15, 20, 25, 30]
    row_short, pl_mode_short, pl_dd_short, pl_txt_short, get_short_len = _list_or_custom(
        label="Length of Paths", options=short_len_opts, default_value=5,
        minv=1, maxv=10_000, step=1, desc_style=DESC
    )
    num_paths_opts_mode2 = [1000, 2000, 3000, 5000, 10000, 20000, 30000, 40000, 50000]
    row_np2, np_mode_2, np_dd_2, np_txt_2, get_num_paths_2 = _list_or_custom(
        label="Number of Paths", options=num_paths_opts_mode2, default_value=1000,
        minv=1, maxv=10_000_000, step=100, desc_style=DESC
    )

    # ---------- Mode 3 ----------
    final_idx   = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1,
                                   description="Index of final residue:",
                                   layout=wide, style=DESC)
    final_chain = W.Text(value="", description="Chain of final residue:",
                         placeholder="B", layout=wide, style=DESC)
    num_paths_opts_mode3 = [1000, 2000, 3000, 5000, 10000, 30000, 50000]
    row_np3, np_mode_3, np_dd_3, np_txt_3, get_num_paths_3 = _list_or_custom(
        label="Number of Paths", options=num_paths_opts_mode3, default_value=1000,
        minv=1, maxv=10_000_000, step=100, desc_style=DESC
    )

    # Actions & output
    btn_submit = W.Button(description="Submit", button_style="success", icon="paper-plane")
    btn_clear  = W.Button(description="Clear",  button_style="warning", icon="trash")
    out        = W.Output()

    # Grouped layouts
    pdb_row = W.HBox([pdb_code, W.HTML("&nbsp;"), or_lbl, pdb_upload, file_lbl],
                     layout=W.Layout(align_items="center", justify_content="flex-start",
                                     flex_flow="row wrap", gap="10px"))

    functional_box = W.VBox([row_big])
    mode2_box = W.VBox([init_idx, init_chain, row_short, row_np2])
    mode3_box = W.VBox([init_idx, init_chain, final_idx, final_chain, row_np3])

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

    # Render
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
    display(W.VBox(children, layout=W.Layout(width="auto")))

    # -------------------- handlers --------------------
    def on_clear(_):
        # Reset text fields
        pdb_code.value = ""
        chain_id.value = ""
        email.value = ""

        # Reset file upload label (keep widget)
        file_lbl.value = "No file chosen"
        try:
            pdb_upload.value = ()
        except Exception:
            try:
                pdb_upload.value = {}
            except Exception:
                pass  # leave as-is if trait is read-only

        # Reset prediction mode & visibility
        pred_type.value = "functional"

        # Reset toggle rows to 'List'
        pl_mode_big.value = "list"
        pl_mode_short.value = "list"
        np_mode_2.value = "list"
        np_mode_3.value = "list"

        # Reset numeric choices to first option in dropdowns
        pl_dd_big.value    = pl_dd_big.options[0];    pl_txt_big.value    = int(pl_dd_big.options[0])
        pl_dd_short.value  = pl_dd_short.options[0];  pl_txt_short.value  = int(pl_dd_short.options[0])
        np_dd_2.value      = np_dd_2.options[0];      np_txt_2.value      = int(np_dd_2.options[0])
        np_dd_3.value      = np_dd_3.options[0];      np_txt_3.value      = int(np_dd_3.options[0])

        # Clear ONLY the log/output box (not the whole cell)
        out.clear_output(wait=True)

    def _collect_pdb_bytes():
        if pdb_upload.value:
            (fname, meta) = next(iter(pdb_upload.value.items()))
            return meta["content"], fname
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

    def on_submit(_):
        # Clear only the log/output area; keep the form
        out.clear_output(wait=True)
        with out:
            try:
                chain_global = chain_id.value.strip()
                if not _is_valid_chain(chain_global):
                    raise ValueError("Chain ID must be a single character (e.g., A).")
                if not _is_valid_email(email.value.strip()):
                    raise ValueError("Invalid email format.")

                # obtain PDB bytes & name
                pdb_bytes, pdb_name = _collect_pdb_bytes()

                # save PDB
                save_path = os.path.join(SAVE_DIR, pdb_name)
                with open(save_path, "wb") as f:
                    f.write(pdb_bytes)
                print(f"Saved local copy: {save_path}")

                # build mcpath_input.txt
                input_path = os.path.join(os.path.dirname(save_path), "mcpath_input.txt")
                mode = pred_type.value
                if mode == "functional":
                    rows = ["1", pdb_name, chain_global, str(get_big_len()), (email.value.strip() or "-")]
                elif mode == "paths_init_len":
                    if not _is_valid_chain(init_chain.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    rows = [
                        "2", pdb_name, chain_global,
                        str(get_short_len()), str(get_num_paths_2()),
                        str(int(init_idx.value)), (init_chain.value or "").strip(),
                        (email.value.strip() or "-")
                    ]
                else:
                    if not _is_valid_chain(init_chain.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    if not _is_valid_chain(final_chain.value or ""):
                        raise ValueError("Chain of final residue must be a single character.")
                    rows = [
                        "3", pdb_name, chain_global,
                        str(int(init_idx.value)), (init_chain.value or "").strip(),
                        str(int(final_idx.value)), (final_chain.value or "").strip(),
                        (email.value.strip() or "-")
                    ]
                with open(input_path, "w") as f:
                    for r in rows: f.write(str(r).strip() + "\n")
                print(f"Input file saved: {input_path}")
                print(f"▶ Using input file: {input_path}")
                print(f"   Mode={mode}, PDB={save_path}, Chain={chain_global}")

                # run readpdb_strict
                run_readpdb = _try_import_readpdb()
                cor_path = None
                if run_readpdb is None:
                    print("ℹ️ readpdb_strict not found; skipping .cor generation.")
                else:
                    try:
                        cor_path = run_readpdb(input_path=input_path)
                    except TypeError:
                        cor_path = run_readpdb(pdb_path=save_path, chain=chain_global)
                    print(f"✔ Wrote: {cor_path}")
                    print(f"✔ COR written: {cor_path}")
                    print(f"   COR exists? {os.path.isfile(cor_path)}\n")
                    try:
                        with open(cor_path, "r") as f:
                            head = "".join([next(f) for _ in range(5)])
                        print("First lines of COR:\n" + head)
                    except Exception:
                        pass

                # Copy *.pdb.cor -> *.cor (atomistic expects <base>.cor)
                try:
                    if cor_path and os.path.isfile(cor_path):
                        base_no_ext = os.path.splitext(save_path)[0]
                        want_cor = f"{base_no_ext}.cor"
                        if os.path.abspath(cor_path) != os.path.abspath(want_cor):
                            shutil.copyfile(cor_path, want_cor)
                        print(f"↪ Copied COR to: {want_cor} (for atomistic)")

                        here = os.path.dirname(os.path.abspath(__file__))
                        param_path = os.path.join(here, "vdw_cns.param")
                        top_path   = os.path.join(here, "pdb_cns.top")
                        print(f"   param_path: {param_path} exists? {os.path.isfile(param_path)}")
                        print(f"   top_path:   {top_path} exists? {os.path.isfile(top_path)}")

                        # ensure mcpath importable
                        pkg_dir = os.path.dirname(os.path.abspath(__file__))
                        root_dir = os.path.dirname(pkg_dir)
                        if root_dir not in sys.path:
                            sys.path.append(root_dir)

                        atom_mod = importlib.import_module("mcpath.atomistic")
                        importlib.reload(atom_mod)

                        base_for_atom = os.path.splitext(save_path)[0]
                        print(f"▶ Running atomistic LJ on base: {base_for_atom}")
                        norm = atom_mod.atomistic(
                            base_for_atom,
                            param_file=param_path,
                            top_file=top_path,
                            rcut=5.0,
                            kT=1.0,
                            save_txt=True
                        )
                        print(f"✔ Saved probability matrix: {base_for_atom}_atomistic.out")
                    else:
                        print("ℹ️ No .cor available; skipping atomistic LJ step.")
                except Exception as e:
                    print("❌ atomistic run failed:", e)

                # Optional POST (preview if no target_url)
                data = {"prediction_mode": mode, FN["chain_id"]: chain_global}
                if pdb_code.value.strip():
                    data[FN["pdb_code"]] = pdb_code.value.strip().upper()
                if email.value.strip():
                    data[FN["email"]] = email.value.strip()
                if mode == "functional":
                    data[FN["path_length"]] = str(get_big_len())
                elif mode == "paths_init_len":
                    data.update({
                        "length_paths":  int(get_short_len()),
                        "number_paths":  int(get_num_paths_2()),
                        "index_initial": int(init_idx.value),
                        "chain_initial": (init_chain.value or "").strip(),
                    })
                else:
                    data.update({
                        "index_initial": int(init_idx.value),
                        "chain_initial": (init_chain.value or "").strip(),
                        "index_final":   int(final_idx.value),
                        "chain_final":   (final_chain.value or "").strip(),
                        "number_paths":  int(get_num_paths_3()),
                    })

                files = {FN["pdb_file"]: (pdb_name, pdb_bytes, "chemical/x-pdb")}

                if not target_url:
                    print("\n(No target_url set) — preview only payload below:\n")
                    preview = dict(data); preview["attached_file"] = pdb_name
                    print(preview)
                    return

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
