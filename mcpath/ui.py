# mcpath/ui.py
import os, re, requests, yaml, importlib, shutil, sys
from IPython.display import display
from IPython.display import display, clear_output # <-- IMPORTED clear_output
import ipywidgets as W

# --- Try to import Colab-specific output clearing ---
try:
    from google.colab import output as colab_output
except ImportError:
    colab_output = None
# ---------------------------------------------------

# ---------- singletons (so New Job/Submit always uses the same log panel) ----------
_FORM_ROOT = None
_LOG_OUT    = None

def _logo_widget(branding: dict):
    try:
        url = branding.get("logo_url")
        if not url: return None
        height = branding.get("logo_height_px", 60)
        align = branding.get("logo_align", "center")
        img_bytes = requests.get(url, timeout=20).content
        jc = {"center":"center", "left":"flex-start", "right":"flex-end"}.get(align, "center")
        return W.HBox([W.Image(value=img_bytes, format="png", layout=W.Layout(height=f"{height}px"))],
                        layout=W.Layout(justify_content=jc))
    except Exception:
        return None

def _make_list_custom_row(description: str,
                        default_val: int,
                        list_vals: list,
                        config: dict,
                        desc_width="auto"):
    DESC = {'description_width': desc_width}
    lbl = W.Label(description, layout=W.Layout(width="auto"))

    # --- Option 1: Dropdown with predefined list ---
    dropdown = W.Dropdown(options=list_vals, value=default_val, layout=W.Layout(width="auto", flex="1 1 auto"))

    # --- Option 2: IntText for custom value ---
    intbox = W.IntText(value=default_val, layout=W.Layout(width="auto", flex="1 1 auto"))

    # --- Toggle button ---
    toggle = W.ToggleButtons(
        options=['List', 'Custom'],
        value='List' if default_val in list_vals else 'Custom',
        layout=W.Layout(width="auto"),
        tooltips=['Select from a predefined list', 'Enter a custom integer value'],
        style={'button_width': '60px'}
    )

    # Box to hold either dropdown or intbox
    value_box = W.Box([dropdown], layout=W.Layout(width="auto", flex="1 1 auto"))

    def _sync_toggle(change):
        if change['new'] == 'List':
            value_box.children = [dropdown]
            if intbox.value in list_vals:
                dropdown.value = intbox.value
            else:
                dropdown.value = list_vals[0] # Default to first in list
        else: # 'Custom'
            value_box.children = [intbox]
            intbox.value = dropdown.value # Carry over value
    toggle.observe(_sync_toggle, 'value')

    # Initialize correct view
    _sync_toggle({'new': toggle.value})

    def get_value():
        return dropdown.value if toggle.value == 'List' else intbox.value

    def set_list(val):
        """Sets the widget state to a specific value (int) from the original list."""
        if val not in list_vals: val = list_vals[0]
        toggle.value = 'List'
        dropdown.value = val
        intbox.value = val
        value_box.children = [dropdown]

    def set_custom(val):
        """Sets the widget state to a specific custom value (int)."""
        toggle.value = 'Custom'
        dropdown.value = list_vals[0] # Reset dropdown to default
        intbox.value = val
        value_box.children = [intbox]

    row = W.HBox([lbl, toggle, value_box], layout=W.Layout(align_items="center", gap="12px", width="auto"))
    return row, {'toggle': toggle, 'dropdown': dropdown, 'intbox': intbox,
                 'get': get_value, 'set_list': set_list, 'set_custom': set_custom}

# -------------------- main UI --------------------
def launch(
    config_path: str,
    target_url: str = None,
    config_overrides: dict = {},
    field_names: dict = {},
    branding: dict = {}
):
    """
    Renders the McPath submit UI.
    config_path: Path to YAML config file.
    target_url: If provided, shows a "Submit" button that POSTs data.
    """
    global _FORM_ROOT, _LOG_OUT
    cfg = {}
    if config_path and os.path.isfile(config_path):
        with open(config_path, "r") as f:
            cfg = yaml.safe_load(f)
    cfg.update(config_overrides)

    FN = { # Default field names for POST request
        "pdb_file": "pdb_file",
        "pdb_code": "pdb_code",
        "chain_id": "chain_id",
        "email": "email",
    }
    FN.update(field_names)

    DESC = {'description_width': 'initial'} #
    wide = W.Layout(width='auto', flex='1 1 auto')

    # -------------------- PDB input --------------------
    pdb_code = W.Text(value=str(cfg.get("pdb_code", "")), description="PDB Code:", placeholder="1ABC", layout=wide, style=DESC)
    chain_id = W.Text(value=str(cfg.get("chain_id", "")), description="Chain ID:", placeholder="A", layout=wide, style=DESC)
    file_lbl = W.Label("or Upload PDB:")
    pdb_upload = W.FileUpload(accept='.pdb', multiple=False, layout=W.Layout(width='auto'))

    def _on_upload_change(_):
        if pdb_upload.value:
            file_info = next(iter(pdb_upload.value.values()))
            file_lbl.value = f"Selected: {file_info['metadata']['name']}"
            pdb_code.value = "" # Clear PDB code
        else:
            file_lbl.value = "or Upload PDB:"
    pdb_upload.observe(_on_upload_change, 'value')

    # -------------------- Mode selection --------------------
    pred_type = W.Dropdown(
        options=[
            ('Functional (all paths)', 'functional'),
            ('Paths from initial residue', 'paths_init'),
            ('Paths from initial to final', 'paths_init_final')
        ],
        value='functional',
        description='Prediction Mode:',
        layout=wide,
        style=DESC
    )

    # -------------------- Mode-specific controls --------------------
    # (These are shown/hidden by _sync_mode)

    # ---------- Mode 2 & 3 ----------
    init_idx_default = 1 # Added to use for resetting
    init_idx = W.BoundedIntText(value=init_idx_default, min=1, max=1_000_000, step=1,
                                description="Index of initial residue:", layout=wide, style=DESC)
    init_chain_default = "" # Added to use for resetting
    init_chain = W.Text(value=init_chain_default, description="Chain of initial residue:", placeholder="A", layout=wide, style=DESC)

    # ---------- Mode 3 ----------
    final_idx_default = 1 # Added to use for resetting
    final_idx    = W.BoundedIntText(value=final_idx_default, min=1, max=1_000_000, step=1,
                                     description="Index of final residue:", layout=wide, style=DESC)
    final_chain_default = "" # Added to use for resetting
    final_chain = W.Text(value=final_chain_default, description="Chain of final residue:", placeholder="B", layout=wide, style=DESC)


    # --- List/Custom rows ---
    desc_width = "180px"
    big_default = cfg.get("big_cutoff_default", 50)
    big_list = cfg.get("big_cutoff_list", [10, 20, 30, 40, 50, 60, 70])
    (big_row, big_ctrl) = _make_list_custom_row("Big cutoff (Å):", big_default, big_list, cfg, desc_width)

    short_len_default = cfg.get("short_path_len_default", 5)
    short_len_list = cfg.get("short_path_len_list", [3, 4, 5, 6, 7])
    (short_row, short_ctrl) = _make_list_custom_row("Short path length:", short_len_default, short_len_list, cfg, desc_width)

    num_paths_mode2_default = cfg.get("num_paths_mode2_default", 500)
    num_paths_mode2_list = cfg.get("num_paths_mode2_list", [100, 200, 500, 1000])
    (np2_row, np2_ctrl) = _make_list_custom_row("Number of paths (Mode 2):", num_paths_mode2_default, num_paths_mode2_list, cfg, desc_width)

    num_paths_mode3_default = cfg.get("num_paths_mode3_default", 500)
    num_paths_mode3_list = cfg.get("num_paths_mode3_list", [100, 200, 500, 1000])
    (np3_row, np3_ctrl) = _make_list_custom_row("Number of paths (Mode 3):", num_paths_mode3_default, num_paths_mode3_list, cfg, desc_width)

    # Group mode-specific controls for easy show/hide
    mode2_controls = W.VBox([init_idx, init_chain, short_row, np2_row])
    mode3_controls = W.VBox([init_idx, init_chain, final_idx, final_chain, np3_row])
    mode1_controls = W.VBox([big_row])


    # -------------------- Other fields --------------------
    email = W.Text(value=str(cfg.get("email", "")), description="Email (optional):", placeholder="user@example.com", layout=wide, style=DESC)

    # -------------------- Buttons --------------------
    btn_new_job = W.Button(description="New Job", button_style='warning', icon='eraser', layout=W.Layout(width='auto'))
    btn_submit = W.Button(description="Submit", button_style='primary', icon='paper-plane', layout=W.Layout(width='auto'))
    if not target_url:
        btn_submit.description = "Preview (no URL)"
        btn_submit.button_style = ''
        btn_submit.icon = 'search'

    # -------------------- Output log --------------------
    if _LOG_OUT is None: # Create singleton on first launch
        _LOG_OUT = W.Output(layout=W.Layout(width='auto', height='auto', min_height='200px', border='1px solid #ccc', padding='8px'))

    # -------------------- Layout --------------------
    logo_widget = _logo_widget(branding)
    layout = W.VBox([
        *([logo_widget] if logo_widget else []),
        W.HBox([pdb_code, chain_id, file_lbl, pdb_upload], layout=W.Layout(gap="12px", align_items="center")),
        pred_type,
        mode1_controls,
        mode2_controls,
        mode3_controls,
        email,
        W.HBox([btn_new_job, btn_submit], layout=W.Layout(gap="12px", margin="10px 0 0 0")),
        _LOG_OUT
    ], layout=W.Layout(width='100%', max_width='800px', display='flex', flex_flow='column nowrap'))
    _FORM_ROOT = layout # Store reference for 'New Job'

    # -------------------- sync UI --------------------
    def _sync_mode(*_):
        mode = pred_type.value
        mode1_controls.layout.display = 'flex' if mode == 'functional' else 'none'
        mode2_controls.layout.display = 'flex' if mode == 'paths_init' else 'none'
        mode3_controls.layout.display = 'flex' if mode == 'paths_init_final' else 'none'

        # Show/hide based on mode
        with _LOG_OUT:
            print(f"Mode set to: {mode}")
    pred_type.observe(_sync_mode, 'value')
    _sync_mode() # Initial sync

    # -------------------- handlers --------------------

    # --- MODIFIED HANDLER: Clears entire cell output and re-displays UI ---
    def on_new_job(_):
        # 1. Clear the entire Colab cell output area
        if colab_output:
            colab_output.clear_output(wait=True)
        else:
            # Fallback for standard Jupyter/IPython
            clear_output(wait=True) 

        # 2. Reset simple text fields
        pdb_code.value = str(cfg.get("pdb_code", ""))
        chain_id.value = str(cfg.get("chain_id", ""))
        email.value    = str(cfg.get("email", ""))

        # 3. Reset file upload
        file_lbl.value = "No file chosen"
        # Clear the upload widget value
        for attempt in ((), {}):  
            try: pdb_upload.value = attempt; break
            except Exception: pass # Colab/ipywidgets quirk handling

        # 4. Reset prediction type (will trigger _sync_mode)
        pred_type.value = "functional"  

        # 5. Reset mode-specific fields
        init_idx.value = init_idx_default
        init_chain.value = init_chain_default
        final_idx.value = final_idx_default
        final_chain.value = final_chain_default

        # 6. Reset list/custom rows to their *original* default values
        big_ctrl['set_list'](big_default)
        short_ctrl['set_list'](short_len_default)  
        np2_ctrl['set_list'](num_paths_mode2_default)
        np3_ctrl['set_list'](num_paths_mode3_default)

        # 7. Clear the internal log widget's content
        if _LOG_OUT:
            _LOG_OUT.clear_output(wait=True)
            
        # 8. Re-display the entire form
        # This is necessary after clearing the cell output
        if _FORM_ROOT:
            display(_FORM_ROOT)
        
        # 9. Print reset message *inside* the log widget
        if _LOG_OUT:
            with _LOG_OUT:
                print("UI successfully reset for a new job.")

    def _collect_pdb_bytes():
        if pdb_upload.value:
            # Note: FileUpload value is a dict in ipywidgets 8
            # In older versions it was a list of dicts. This handles the common case.
            file_info = next(iter(pdb_upload.value.values()))
            return file_info["content"], file_info["metadata"]["name"]
        
        code = pdb_code.value.strip()
        if not _is_valid_pdb_code(code):
            raise ValueError("PDB code must be exactly 4 alphanumeric characters (or upload a file).")
        
        url = f"https://files.rcsb.org/download/{code}.pdb"
        print(f"Downloading PDB: {url}")
        r = requests.get(url, timeout=20)
        r.raise_for_status()
        print("Download complete.")
        return r.content, f"{code}.pdb"

    def _is_valid_pdb_code(code: str):
        return code and len(code) == 4 and code.isalnum()

    def _is_valid_chain(chain: str):
        return len(chain) == 1 and chain.strip() # Must be one char, not space

    def _try_import_readpdb():
        try:
            # Try to find it relative to the module 'mcpath'
            if 'mcpath' in sys.modules:
                pkg_dir = os.path.dirname(sys.modules['mcpath'].__file__)
                if pkg_dir not in sys.path:
                    sys.path.append(pkg_dir)
            
            # Try to find it relative to *this* file's location
            # This is fragile in notebooks, but a good fallback
            try:
                here = os.path.dirname(os.path.abspath(__file__))
                if here not in sys.path: sys.path.append(here)
                root = os.path.dirname(here)
                if root not in sys.path: sys.path.append(root)
            except NameError:
                pass # __file__ not defined

            mod = importlib.import_module("mcpath.readpdb")
            importlib.reload(mod) # ensure latest is used
            return mod.run_readpdb
        except Exception as e:
            print(f"Warning: could not import readpdb: {e}")
            return None

    def on_submit(_):
        with _LOG_OUT:
            try:
                _LOG_OUT.clear_output(wait=True)
                print("Starting job...")

                pdb_bytes, pdb_name = _collect_pdb_bytes()
                chain_global = chain_id.value.strip()
                if not _is_valid_chain(chain_global):
                    raise ValueError("Chain ID must be exactly one character.")
                
                # --- Create input files for readpdb/atomistic ---
                # (This logic is preserved from the original script)
                job_dir = cfg.get("job_dir", ".")
                os.makedirs(job_dir, exist_ok=True)
                
                # Save PDB
                save_path = os.path.join(job_dir, pdb_name)
                with open(save_path, "wb") as f: f.write(pdb_bytes)
                print(f"PDB file saved: {save_path}")

                # Save INPUT file
                input_path = os.path.join(job_dir, "INPUT")
                mode = pred_type.value

                rows = []
                if mode == "functional":
                    rows = ["1", pdb_name, chain_global, str(big_ctrl['get']()), (email.value.strip() or "-")]
                elif mode == "paths_init":
                    if not _is_valid_chain(init_chain.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    rows = ["2", pdb_name, chain_global, str(short_ctrl['get']()), str(np2_ctrl['get']()),
                            str(int(init_idx.value)), (init_chain.value or "").strip(), (email.value.strip() or "-")]
                else: # paths_init_final
                    if not _is_valid_chain(init_chain.value or ""):
                        raise ValueError("Chain of initial residue must be a single character.")
                    if not _is_valid_chain(final_chain.value or ""):
                        raise ValueError("Chain of final residue must be a single character.")
                    rows = ["3", pdb_name, chain_global,
                            str(int(init_idx.value)), (init_chain.value or "").strip(),
                            str(int(final_idx.value)), (final_chain.value or "").strip(),
                            str(np3_ctrl['get']()), # <-- ADDED num_paths_3 FOR MODE 3
                            (email.value.strip() or "-")]
                    
                    # --- CHECK: Original code for mode 3 seemed to be missing number_paths ---
                    # Original:
                    # rows = ["3", pdb_name, chain_global,
                    #         str(int(init_idx.value)), (init_chain.value or "").strip(),
                    #         str(int(final_idx.value)), (final_chain.value or "").strip(),
                    #         (email.value.strip() or "-")]
                    #
                    # This looked like an error, as mode 3 has a "Number of Paths" widget. 
                    # I've corrected it to include `get_num_paths_3()` as seen above.
                    # If this was intentional, you can revert to the original lines.
                    # --- END CHECK ---

                with open(input_path, "w") as f:
                    for r in rows: f.write(str(r).strip() + "\n")
                print(f"▶ Using input file: {input_path}")
                print(f"   Mode={mode}, PDB={save_path}, Chain={chain_global}")

                run_readpdb = _try_import_readpdb()
                cor_path = None
                if run_readpdb:
                    print(f"▶ Running readpdb...")
                    try:
                        cor_path = run_readpdb(input_path=input_path)
                    except TypeError:
                        # Fallback for older interface
                        cor_path = run_readpdb(pdb_path=save_path, chain=chain_global)
                    print(f"✔ COR written: {cor_path}")
                    if cor_path and os.path.isfile(cor_path):
                        print(f"   COR exists? True\n")
                        try:
                            with open(cor_path, "r") as f:
                                head = "".join([next(f) for _ in range(5)])
                            print("First lines of COR:\n" + head)
                        except Exception:
                            pass
                    else:
                        print(f"   COR exists? False\n")

                try:
                    if cor_path and os.path.isfile(cor_path):
                        base_no_ext = os.path.splitext(save_path)[0]
                        want_cor = f"{base_no_ext}.cor"
                        if os.path.abspath(cor_path) != os.path.abspath(want_cor):
                            shutil.copyfile(cor_path, want_cor)
                            print(f"↪ Copied COR to: {want_cor} (for atomistic)")

                        # __file__ might not be defined in Colab, assume relative paths
                        pkg_dir = os.path.dirname(sys.modules['mcpath'].__file__ if 'mcpath' in sys.modules else ".")
                        
                        param_path = os.path.join(pkg_dir, "vdw_cns.param")
                        top_path   = os.path.join(pkg_dir, "pdb_cns.top")
                        print(f"   param_path: {param_path} exists? {os.path.isfile(param_path)}")
                        print(f"   top_path:   {top_path} exists? {os.path.isfile(top_path)}")

                        atom_mod = importlib.import_module("mcpath.atomistic")
                        importlib.reload(atom_mod) # ensure latest is used

                        base_for_atom = os.path.splitext(save_path)[0]
                        print(f"▶ Running atomistic LJ on base: {base_for_atom}")
                        atom_mod.run_atomistic(base_for_atom, param_path, top_path)
                    else:
                        print("ℹ️ No .cor available; skipping atomistic LJ step.")
                except Exception as e:
                    print(f"❌ atomistic run failed: {e}")
                    import traceback
                    traceback.print_exc()

                # --- Prepare web request ---
                data = {"prediction_mode": mode, FN["chain_id"]: chain_global}
                if pdb_code.value.strip():
                    data[FN["pdb_code"]] = pdb_code.value.strip()
                if email.value.strip():
                    data[FN["email"]] = email.value.strip()
                
                if mode == "functional":
                    data.update({"big_cutoff": int(big_ctrl['get']())})
                elif mode == "paths_init":
                    data.update({"short_path_len": int(short_ctrl['get']()),
                                 "number_paths": int(np2_ctrl['get']()),
                                 "index_initial": int(init_idx.value),
                                 "chain_initial": (init_chain.value or "").strip()})
                else: # paths_init_final
                    data.update({"index_initial": int(init_idx.value),
                                 "chain_initial": (init_chain.value or "").strip(),
                                 "index_final": int(final_idx.value),
                                 "chain_final": (final_chain.value or "").strip(),
                                 "number_paths": int(np3_ctrl['get']())}) # <-- Corrected this payload

                files = {FN["pdb_file"]: (pdb_name, pdb_bytes, "chemical/x-pdb")}
                if not target_url:
                    print("\n(No target_url set) — preview only payload below:\n")
                    preview = dict(data); preview["attached_file"] = pdb_name
                    print(yaml.dump(preview, default_flow_style=False))
                    return

                print(f"Submitting to {target_url} …")
                r = requests.post(target_url, data=data, files=files, timeout=30)
                print(f"Response [{r.status_code}]:")
                print(r.text)

            except Exception as e:
                print("❌", e)
                import traceback
                traceback.print_exc()

    btn_new_job.on_click(on_new_job)
    btn_submit.on_click(on_submit)

    display(layout)
