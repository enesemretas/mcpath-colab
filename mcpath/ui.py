# mcpath/ui.py
import os, re, requests, yaml, importlib, shutil, sys, glob
from IPython.display import display, clear_output
import ipywidgets as W

# ---------- singletons (so New Job/Submit always uses the same log panel) ----------
_FORM_ROOT = None
_LOG_OUT   = None

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

    # Always-present fields
    chain_id   = W.Text(value=str(cfg.get("chain_id", "")), description="Chain ID:",  layout=wide, style=DESC)
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

    # --- HARD RESET & RERUN: rebuilds the entire UI fresh (including upload) ---
    def on_new_job(_):
        global _FORM_ROOT, _LOG_OUT
        try:
            if _LOG_OUT:
                _LOG_OUT.clear_output(wait=True)
        except Exception:
            pass
        try:
            if _FORM_ROOT:
                _FORM_ROOT.close()  # remove old widget tree from front-end
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
        # add _1, _2, ...
        k = 1
        while True:
            cand = f"{dst}_{k}"
            if not os.path.exists(cand):
                shutil.copyfile(src_path, cand)
                return cand
            k += 1

    def on_submit(_):
        _LOG_OUT.clear_output(wait=True)  # always clear the shared area first
        with _LOG_OUT:
            try:
                chain_global = chain_id.value.strip()
                if not _is_valid_chain(chain_global):
                    raise ValueError("Chain ID must be a single character (e.g., A).")
                if not _is_valid_email(email.value.strip()):
                    raise ValueError("Invalid email format.")

                pdb_bytes, pdb_name = _collect_pdb_bytes()

                save_path = os.path.join(SAVE_DIR, pdb_name)
                with open(save_path, "wb") as f: f.write(pdb_bytes)
                print(f"Saved local copy: {save_path}")

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
                    for r in rows: f.write(str(r).strip() + "\n")
                print(f"Input file saved: {input_path}")
                print(f"▶ Using input file: {input_path}")
                print(f"    Mode={mode}, PDB={save_path}, Chain={chain_global}")

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
                    print(f"    COR exists? {os.path.isfile(cor_path)}\n")
                    try:
                        with open(cor_path, "r") as f:
                            head = "".join([next(f) for _ in range(5)])
                        print("First lines of COR:\n" + head)
                    except Exception:
                        pass

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
                        print(f"    param_path: {param_path} exists? {os.path.isfile(param_path)}")
                        print(f"    top_path:   {top_path} exists? {os.path.isfile(top_path)}")

                        pkg_dir = os.path.dirname(os.path.abspath(__file__))
                        root_dir = os.path.dirname(pkg_dir)
                        if root_dir not in sys.path: sys.path.append(root_dir)

                        atom_mod = importlib.import_module("mcpath.atomistic")
                        importlib.reload(atom_mod)

                        base_for_atom = os.path.splitext(save_path)[0]
                        print(f"▶ Running atomistic LJ on base: {base_for_atom}")
                        norm = atom_mod.atomistic(base_for_atom, param_file=param_path, top_file=top_path,
                                                  rcut=5.0, kT=1.0, save_txt=True)
                        print(f"✔ Saved probability matrix: {base_for_atom}_atomistic.out")

                        # If mode == functional, run infinite.py and copy to fixed names, then run closeness.py
                        if mode == "functional":
                            try:
                                # 1) Run infinite (long path generator)
                                inf_mod = importlib.import_module("mcpath.infinite")
                                importlib.reload(inf_mod)

                                work_dir = os.path.dirname(save_path)
                                old_cwd = os.getcwd()
                                try:
                                    os.chdir(work_dir)
                                    steps = int(get_big_len())
                                    print(f"▶ Running infinite path generator: ProteinName={os.path.basename(base_for_atom)}, steps={steps}")
                                    path_arr = inf_mod.infinite(os.path.basename(base_for_atom), path_length=steps, pottype='1')
                                    pl = path_arr.shape[1]
                                    path_src = f"{os.path.basename(base_for_atom)}_atomistic_{pl}steps_infinite.path"
                                    path_src = os.path.join(work_dir, path_src)
                                finally:
                                    os.chdir(old_cwd)

                                # 2) Copy sources to fixed names (with auto-suffix if exist)
                                cor_src  = f"{base_for_atom}.cor"
                                atom_src = f"{base_for_atom}_atomistic.out"

                                coor_fixed = _copy_unique(cor_src,  "coor_file", work_dir=os.path.dirname(cor_src))
                                atom_fixed = _copy_unique(atom_src, "atom_file", work_dir=os.path.dirname(atom_src))
                                path_fixed = _copy_unique(path_src, "path_file", work_dir=os.path.dirname(path_src))

                                print(f"✔ Fixed copies:")
                                print(f"   {cor_src}  →  {coor_fixed}")
                                print(f"   {atom_src} →  {atom_fixed}")
                                print(f"   {path_src} →  {path_fixed}")

                                # 3) Run closeness (shortest paths & centrality)
                                close_mod = importlib.import_module("mcpath.closeness")
                                importlib.reload(close_mod)

                                try:
                                    os.chdir(work_dir)
                                    print("▶ Running closeness (shortest paths + closeness centrality)…")
                                    close_mod.main()
                                    print("✔ Wrote: shortest_paths_all_pairs_chain1")
                                    print("✔ Wrote: closeness_float")
                                finally:
                                    os.chdir(old_cwd)

                            except Exception as e_inf:
                                print("❌ infinite/closeness pipeline failed:", e_inf)

                    else:
                        print("ℹ️ No .cor available; skipping atomistic LJ step.")
                except Exception as e:
                    print("❌ atomistic run failed:", e)

                # (Submission preview / POST unchanged)
                data = {"prediction_mode": mode, FN["chain_id"]: chain_global}
                if pdb_code.value.strip():
                    data[FN["pdb_code"]] = pdb_code.value.strip().upper()
                if email.value.strip():
                    data[FN["email"]] = email.value.strip()
                if mode == "functional":
                    data[FN["path_length"]] = str(get_big_len())
                elif mode == "paths_init_len":
                    data.update({"length_paths": int(get_short_len()),
                                 "number_paths": int(get_num_paths_2()),
                                 "index_initial": int(init_idx.value),
                                 "chain_initial": (init_chain.value or "").strip()})
                else:
                    data.update({"index_initial": int(init_idx.value),
                                 "chain_initial": (init_chain.value or "").strip(),
                                 "index_final": int(final_idx.value),
                                 "chain_final": (final_chain.value or "").strip(),
                                 "number_paths": int(get_num_paths_3())})

                files = {FN["pdb_file"]: (pdb_name, pdb_bytes, "chemical/x-pdb")}
                if not target_url:
                    print("\n(No target_url set) — preview only payload below:\n")
                    preview = dict(data); preview["attached_file"] = pdb_name
                    print(preview)
                    return

                print(f"Submitting to {target_url} …")
                r = requests.post(target_url, data=data, files=files, timeout=180)
                print("HTTP", r.status_code)
                try: print("JSON:", r.json())
                except Exception: print("Response (≤800 chars):\n", r.text[:800])

            except Exception as e:
                print("❌", e)

    btn_new_job.on_click(on_new_job)
    btn_submit.on_click(on_submit)
