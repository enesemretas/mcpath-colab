# mcpath/ui.py
import os, re, sys, yaml, importlib, requests, traceback
from IPython.display import display, clear_output
import ipywidgets as W

# -------------------- validators & helpers --------------------
def _is_valid_pdb_code(c): return bool(re.fullmatch(r"[0-9A-Za-z]{4}", c.strip()))
def _is_valid_chain(ch):   return bool(re.fullmatch(r"[A-Za-z0-9]", ch.strip()))
def _is_valid_email(s):    return (not s.strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s.strip()))

def _fetch_rcsb(code: str) -> bytes:
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    r = requests.get(url, timeout=60); r.raise_for_status()
    return r.content

# ---------- “List | Custom | Value” control ----------
def list_custom_row(label: str, options, default_value, minv, maxv, step=1, desc_style=None):
    """
    Row layout:
      [ <label> ] [ List | Custom ] [ values widget (dropdown or BoundedIntText) ]

    Returns:
      row_box, get_value_fn, (mode_toggle, container, dropdown, inttext)
    """
    opts = [] if options is None else list(options)
    has_list = len(opts) > 0
    if has_list:
        # normalize to unique ints, keep default in list
        try:
            opts = sorted({int(x) for x in opts})
        except Exception:
            opts = sorted({int(float(x)) for x in opts})
    try:
        default_value = int(default_value)
    except Exception:
        default_value = minv

    # Widgets
    lab = W.Label(f"{label}:", layout=W.Layout(width="180px", min_width="140px"))
    mode = W.ToggleButtons(
        options=[("List", "list"), ("Custom", "custom")] if has_list else [("Custom", "custom")],
        value="list" if has_list else "custom",
        layout=W.Layout(width="200px"),
        style=desc_style
    )
    dd = W.Dropdown(
        options=(opts + ([default_value] if (has_list and default_value not in opts) else [])) if has_list else [],
        value=(default_value if has_list else None),
        layout=W.Layout(width="260px", min_width="220px"),
        style=desc_style
    )
    txt = W.BoundedIntText(
        value=default_value, min=minv, max=maxv, step=step,
        layout=W.Layout(width="260px", min_width="220px"),
        style=desc_style
    )

    container = W.Box([dd if has_list else txt])
    container.layout = W.Layout(align_items="center")

    def _on_mode(change):
        if change["new"] == "list" and has_list:
            container.children = [dd]
        else:
            container.children = [txt]
    mode.observe(_on_mode, names="value")

    def get_value():
        if mode.value == "list" and has_list:
            return int(dd.value)
        return int(txt.value)

    row = W.HBox(
        [lab, mode, container],
        layout=W.Layout(align_items="center", gap="12px", flex_flow="row wrap")
    )
    return row, get_value, (mode, container, dd, txt)

def _logo_widget(branding: dict):
    url = (branding.get("logo_url") or "").strip()
    if not url: return None
    height = int(branding.get("logo_height", 96))
    align  = (branding.get("logo_align") or "center").lower()
    try:
        img_bytes = requests.get(url, timeout=20).content
        img = W.Image(value=img_bytes, format="png", layout=W.Layout(height=f"{height}px"))
        jc = {"center":"center","left":"flex-start","right":"flex-end"}.get(align,"center")
        return W.HBox([img], layout=W.Layout(justify_content=jc))
    except Exception:
        return None

# -------------------- main UI --------------------
def launch(
    defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
    show_title="MCPath-style Parameters"
):
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=30).text)
    target_url = cfg.get("target_url", "").strip()
    FN = cfg["field_names"]

    SAVE_DIR = cfg.get("save_dir", "/content")
    os.makedirs(SAVE_DIR, exist_ok=True)

    logo = _logo_widget(cfg.get("branding", {}))
    DESC = {'description_width': '240px'}
    wide = W.Layout(width="420px", min_width="360px")

    # ---------- PDB: code OR upload ----------
    pdb_code   = W.Text(value=str(cfg.get("pdb_code", "")),
                        description="PDB code:", layout=wide, style=DESC)
    or_lbl     = W.HTML("<b>&nbsp;&nbsp;or&nbsp;&nbsp;</b>")
    pdb_upload = W.FileUpload(accept=".pdb", multiple=False, description="Choose file")
    file_lbl   = W.Label("No file chosen")
    def _on_upload_change(_):
        try:
            file_lbl.value = next(iter(pdb_upload.value.keys())) if pdb_upload.value else "No file chosen"
        except Exception:
            file_lbl.value = "No file chosen"
    pdb_upload.observe(_on_upload_change, names="value")

    # ---------- Always-present fields ----------
    chain_id = W.Text(value=str(cfg.get("chain_id", "")),
                      description="Chain ID:", layout=wide, style=DESC)
    email    = W.Text(value=str(cfg.get("email", "")),
                      description="Email (opt):", layout=wide, style=DESC)

    # ---------- Prediction type ----------
    pred_type = W.RadioButtons(
        options=[
            ("Functional Residues", "functional"),
            ("Allosteric Paths (initial residue + path length)", "paths_init_len"),
            ("Allosteric Paths (initial & final residues)", "paths_init_final"),
        ],
        value="functional", description="Prediction:", layout=W.Layout(width="auto"), style=DESC
    )

    # ---------- Parameter rows (List | Custom | Value) ----------
    big_opts = cfg.get("path_length_options",
        [100000, 200000, 300000, 400000, 500000, 750000, 1000000, 2000000])
    big_default = int(cfg.get("path_length", big_opts[0] if big_opts else 100000))
    row_big, get_big_len, _ = list_custom_row(
        "Path length", big_opts, big_default, minv=1, maxv=10_000_000, step=1000, desc_style=DESC
    )

    short_len_opts = [5, 8, 10, 13, 15, 20, 25, 30]
    row_short, get_short_len, _ = list_custom_row(
        "Length of Paths", short_len_opts, 5, minv=1, maxv=10_000, step=1, desc_style=DESC
    )

    num_paths_opts_mode2 = [1000, 2000, 3000, 5000, 10000, 20000, 30000, 40000, 50000]
    row_np2, get_num_paths_2, _ = list_custom_row(
        "Number of Paths", num_paths_opts_mode2, 1000, minv=1, maxv=10_000_000, step=100, desc_style=DESC
    )

    num_paths_opts_mode3 = [1000, 2000, 3000, 5000, 10000, 30000, 50000]
    row_np3, get_num_paths_3, _ = list_custom_row(
        "Number of Paths", num_paths_opts_mode3, 1000, minv=1, maxv=10_000_000, step=100, desc_style=DESC
    )

    # ---------- Other inputs for modes 2/3 ----------
    init_idx   = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1,
                                  description="Index of initial residue:", layout=wide, style=DESC)
    init_chain = W.Text(value="", description="Chain of initial residue:", placeholder="A",
                        layout=wide, style=DESC)
    final_idx   = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1,
                                   description="Index of final residue:", layout=wide, style=DESC)
    final_chain = W.Text(value="", description="Chain of final residue:", placeholder="B",
                         layout=wide, style=DESC)

    # ---------- Actions & output ----------
    btn_submit = W.Button(description="Submit", button_style="success", icon="paper-plane")
    btn_clear  = W.Button(description="Clear",  button_style="warning", icon="trash")
    out        = W.Output()

    # ---------- Grouped layouts per section ----------
    pdb_row = W.HBox([pdb_code, W.HTML("&nbsp;"), or_lbl, pdb_upload, file_lbl],
                     layout=W.Layout(align_items="center", justify_content="flex-start",
                                     flex_flow="row wrap", gap="10px"))

    functional_box = W.VBox([W.HTML("<b>Input:</b>"), row_big])
    mode2_box = W.VBox([init_idx, init_chain, W.HTML("<b>Input:</b>"), row_short, row_np2])
    mode3_box = W.VBox([init_idx, init_chain, final_idx, final_chain, W.HTML("<b>Input:</b>"), row_np3])

    def _sync_mode(*_):
        functional_box.layout.display = ""
        mode2_box.layout.display = "none"
        mode3_box.layout.display = "none"
        if pred_type.value == "paths_init_len":
            functional_box.layout.display = "none"; mode2_box.layout.display = ""
        elif pred_type.value == "paths_init_final":
            functional_box.layout.display = "none"; mode3_box.layout.display = ""
    _sync_mode(); pred_type.observe(_sync_mode, names="value")

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
    display(W.VBox(children, layout=W.Layout(width="auto")))

    # -------------------- helpers --------------------
    def _collect_pdb_bytes():
        if pdb_upload.value:
            try:
                (fname, meta) = next(iter(pdb_upload.value.items()))
                return meta["content"], fname
            except Exception:
                pass
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

    # -------------------- submit handler --------------------
    def on_submit(_):
        with out:
            clear_output()
            try:
                chain_global = chain_id.value.strip()
                if not _is_valid_chain(chain_global):
                    raise ValueError("Chain ID must be a single character (e.g., A).")
                if not _is_valid_email(email.value.strip()):
                    raise ValueError("Invalid email format.")

                # obtain PDB bytes & name
                pdb_bytes, pdb_name = _collect_pdb_bytes()

                # ---- Save PDB to SAVE_DIR ----
                save_path = os.path.join(SAVE_DIR, pdb_name)
                with open(save_path, "wb") as f:
                    f.write(pdb_bytes)
                print(f"Saved local copy: {save_path}")

                # ---- Build input file rows based on mode ----
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
                else:  # paths_init_final
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
                    f.write("\n".join([str(r).strip() for r in rows]))
                print(f"Input file saved: {input_path}")

                # ---- Run readpdb_strict immediately (if available) ----
                run_readpdb = _try_import_readpdb()
                cor_path = None
                if run_readpdb is None:
                    print("ℹ️ readpdb_strict not found; skipping .cor generation.")
                else:
                    try:
                        cor_path = run_readpdb(input_path=input_path)
                    except TypeError:
                        cor_path = run_readpdb(pdb_path=save_path, chain=chain_global)
                    print(f"✔ COR written: {cor_path}")
                    print(f"   COR exists? {os.path.isfile(cor_path)}")

                    # quick preview
                    try:
                        with open(cor_path, "r") as f:
                            head = "".join([next(f) for _ in range(5)])
                        print("\nFirst lines of COR:\n" + head)
                    except Exception:
                        pass

                # ---- Normalize COR name for atomistic (atomistic expects <base>.cor) ----
                desired_cor = None
                if cor_path and os.path.isfile(cor_path):
                    desired_cor = os.path.splitext(save_path)[0] + ".cor"   # /content/XXXX.cor
                    if cor_path != desired_cor:
                        import shutil
                        try:
                            shutil.copyfile(cor_path, desired_cor)
                            print(f"↪ Copied COR to: {desired_cor} (for atomistic)")
                        except Exception as ce:
                            print("   (Could not copy COR):", ce)
                            desired_cor = cor_path

                # ---- Run atomistic (creates <pdb>_atomistic.out) ----
                try:
                    if not desired_cor or not os.path.isfile(desired_cor):
                        print("ℹ️ No .cor available; skipping atomistic LJ step.")
                    else:
                        # ensure mcpath is importable (useful in Colab)
                        pkg_dir  = os.path.dirname(os.path.abspath(__file__))
                        root_dir = os.path.dirname(pkg_dir)
                        if root_dir not in sys.path:
                            sys.path.append(root_dir)
                            print("   Added to sys.path:", root_dir)
                        importlib.invalidate_caches()

                        try:
                            atom_mod = importlib.import_module("mcpath.atomistic")
                        except ModuleNotFoundError as me:
                            print("❌ Could not import mcpath.atomistic:", me)
                            print("   Files in pkg_dir:", os.listdir(pkg_dir))
                            raise
                        atom_mod = importlib.reload(atom_mod)

                        # Absolute paths for param/top inside mcpath
                        param_path = os.path.join(pkg_dir, "vdw_cns.param")
                        top_path   = os.path.join(pkg_dir, "pdb_cns.top")
                        print("   param_path:", param_path, "exists?", os.path.isfile(param_path))
                        print("   top_path:  ", top_path,   "exists?", os.path.isfile(top_path))
                        if not os.path.isfile(param_path): raise FileNotFoundError(param_path)
                        if not os.path.isfile(top_path):   raise FileNotFoundError(top_path)

                        # Pass the base WITHOUT extension (same base as the .cor)
                        protein_base = os.path.splitext(save_path)[0]  # /content/XXXX
                        print("▶ Running atomistic LJ on base:", protein_base)
                        norm = atom_mod.atomistic(
                            protein_base,
                            param_file=param_path,
                            top_file=top_path,
                            rcut=5.0,
                            kT=1.0,
                            save_txt=True
                        )
                        out_path = f"{protein_base}_atomistic.out"
                        print(f"✔ Saved probability matrix: {out_path} (exists? {os.path.isfile(out_path)})")
                        if os.path.isfile(out_path):
                            import numpy as np
                            preview = np.loadtxt(out_path, max_rows=5)
                            print("Preview (first 5 rows):")
                            print(preview)
                except Exception as e:
                    print("❌ atomistic run failed:", e)
                    traceback.print_exc()

                # ---- Optional POST ----
                data = {"prediction_mode": mode, FN["chain_id"]: chain_global}
                if pdb_code.value.strip(): data[FN["pdb_code"]] = pdb_code.value.strip().upper()
                if email.value.strip():    data[FN["email"]] = email.value.strip()
                if mode == "functional":
                    data[FN["path_length"]] = str(get_big_len())
                elif mode == "paths_init_len":
                    data.update({"length_paths": int(get_short_len()), "number_paths": int(get_num_paths_2())})
                else:
                    data.update({"index_initial": int(init_idx.value), "index_final": int(final_idx.value),
                                 "number_paths": int(get_num_paths_3())})

                if not target_url:
                    print("\n(No target_url set) — preview only payload below:\n")
                    print(dict(data, attached_file=pdb_name))
                    return

                print(f"Submitting to {target_url} …")
                r = requests.post(target_url, data=data,
                                  files={FN["pdb_file"]: (pdb_name, pdb_bytes, "chemical/x-pdb")},
                                  timeout=180)
                print("HTTP", r.status_code)
                try: print("JSON:", r.json())
                except Exception: print("Response:", r.text[:600])

            except Exception as e:
                print("❌", e)
                traceback.print_exc()

    btn_submit.on_click(on_submit)
    btn_clear.on_click(lambda _: clear_output())
