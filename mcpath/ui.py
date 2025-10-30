# mcpath/ui.py
# UI-only form for MCPath with a single ipywidgets FileUpload (no Colab uploader).
# Run once in Colab before calling launch():
#   !pip -q install "ipywidgets==8.*" pyyaml requests
#   from google.colab import output; output.enable_custom_widget_manager()

import os, re, requests, yaml
from IPython.display import display, clear_output
import ipywidgets as W

# -------------------- validators & helpers --------------------
def _is_valid_pdb_code(c): return bool(re.fullmatch(r"[0-9A-Za-z]{4}", (c or "").strip()))
def _is_valid_chain(ch):   return bool(re.fullmatch(r"[A-Za-z0-9]", (ch or "").strip()))
def _is_valid_email(s):    return (not (s or "").strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s.strip()))

def _fetch_rcsb(code: str) -> bytes:
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    r = requests.get(url, timeout=60); r.raise_for_status()
    return r.content

def _list_or_custom(label: str, options, default_value, minv, maxv, step=1):
    options = sorted({int(x) for x in options})
    default_value = int(default_value)
    if default_value not in options:
        options = sorted(options + [default_value])

    mode = W.ToggleButtons(options=[("List","list"), ("Custom","custom")], value="list", description="Input:")
    dd   = W.Dropdown(options=options, value=default_value, description=label)
    txt  = W.BoundedIntText(value=default_value, min=minv, max=maxv, step=step, description=label)

    def _sync(*_):
        dd.layout.display  = ""   if mode.value == "list"   else "none"
        txt.layout.display = ""   if mode.value == "custom" else "none"
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
        return W.HBox([img], layout=W.Layout(justify_content=jc))
    except Exception:
        return None

# -------------------- main UI --------------------
def launch(
    defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
    show_title="MCPath (Python readpdb)"
):
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=30).text)
    target_url = cfg.get("target_url", "").strip()
    FN = cfg["field_names"]

    # Branding
    logo = _logo_widget(cfg.get("branding", {}))

    # --- PDB: code OR single-file upload (ipywidgets only, no Colab uploader) ---
    pdb_code = W.Text(value=str(cfg.get("pdb_code", "")), description="PDB code:", layout=W.Layout(width="350px"))
    or_lbl   = W.HTML("<b>&nbsp;&nbsp;or&nbsp;&nbsp;</b>")

    uploader = W.FileUpload(
        accept=".pdb",
        multiple=False,
        description="Choose file (1)",
        icon="upload",
        layout=W.Layout(width="240px", height="38px"),
        disabled=False,
    )
    file_lbl = W.HTML("No file chosen")
    picked = {"name": None, "bytes": None, "saved_path": None}

    def _on_upload_change(change):
        if not uploader.value:
            file_lbl.value = "No file chosen"
            picked.update(name=None, bytes=None, saved_path=None)
            return
        (fname, meta) = next(iter(uploader.value.items()))
        # force .pdb extension if missing
        name = fname if fname.lower().endswith(".pdb") else (os.path.splitext(fname)[0] + ".pdb")
        data = meta["content"]
        with open(name, "wb") as f:
            f.write(data)
        picked.update(name=name, bytes=data, saved_path=os.path.join(os.getcwd(), name))
        file_lbl.value = f"{name} (uploaded)"
        # IMPORTANT: do NOT assign to uploader.value (read-only); leave it as-is

    uploader.observe(_on_upload_change, names="value")

    # --- Always-present fields ---
    chain_id = W.Text(value=str(cfg.get("chain_id", "A")), description="Chain ID:", layout=W.Layout(width="350px"))
    email    = W.Text(value=str(cfg.get("email", "")), description="Email (opt):", layout=W.Layout(width="350px"))

    # --- Prediction type ---
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

    # --- Functional (long path length) ---
    big_opts = cfg.get("path_length_options", [
        100000, 200000, 300000, 400000, 500000, 750000,
        1000000, 2000000, 3000000, 4000000
    ])
    big_default = int(cfg.get("path_length", big_opts[0] if big_opts else 100000))
    (pl_mode_big, pl_dd_big, pl_txt_big, get_big_len) = _list_or_custom(
        label="Path length:", options=big_opts, default_value=big_default,
        minv=1, maxv=10_000_000, step=1000
    )

    # --- Mode 2: initial residue + short path length + number of paths ---
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

    # --- Mode 3: initial & final residues + number of paths ---
    final_idx   = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1, description="Index of final residue")
    final_chain = W.Text(value="", description="Chain of final residue", placeholder="B", layout=W.Layout(width="260px"))
    num_paths_opts_mode3 = [1000, 2000, 3000, 5000, 10000, 30000, 50000]
    (np_mode_3, np_dd_3, np_txt_3, get_num_paths_3) = _list_or_custom(
        label="Number of Paths:", options=num_paths_opts_mode3, default_value=1000,
        minv=1, maxv=10_000_000, step=100
    )

    # --- Actions & output ---
    btn_submit = W.Button(description="Submit", button_style="success", icon="paper-plane",
                          layout=W.Layout(width="220px", height="40px"))
    btn_clear  = W.Button(description="Clear",  button_style="warning", icon="trash",
                          layout=W.Layout(width="220px", height="40px"))
    out        = W.Output()
    status_bar = W.HTML("")

    # --- Grouped layouts ---
    pdb_row = W.HBox([pdb_code, W.HTML("&nbsp;"), or_lbl, uploader, W.HTML("&nbsp;"), file_lbl])

    functional_box = W.VBox([pl_mode_big, pl_dd_big, pl_txt_big])
    mode2_box = W.VBox([init_idx, init_chain, pl_mode_short, pl_dd_short, pl_txt_short, np_mode_2, np_dd_2, np_txt_2])
    mode3_box = W.VBox([init_idx, init_chain, final_idx, final_chain, np_mode_3, np_dd_3, np_txt_3])

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
        W.HBox([btn_submit, W.HTML("&nbsp;"), btn_clear]),
        status_bar,
        W.HTML("<hr>"), out
    ]
    display(W.VBox(children))

    # -------------------- handlers --------------------
    def on_clear(_):
        pdb_code.value = ""
        chain_id.value = str(cfg.get("chain_id", "A"))
        email.value = ""
        pred_type.value = "functional"
        pl_mode_big.value = "list"; pl_dd_big.value = pl_dd_big.options[0]; pl_txt_big.value = int(pl_dd_big.options[0])
        init_idx.value = 1; init_chain.value = ""
        pl_mode_short.value = "list"; pl_dd_short.value = pl_dd_short.options[0]; pl_txt_short.value = int(pl_dd_short.options[0])
        np_mode_2.value = "list"; np_dd_2.value = np_dd_2.options[0]; np_txt_2.value = int(np_dd_2.options[0])
        final_idx.value = 1; final_chain.value = ""
        np_mode_3.value = "list"; np_dd_3.value = np_dd_3.options[0]; np_txt_3.value = int(np_dd_3.options[0])
        picked.update(name=None, bytes=None, saved_path=None)
        file_lbl.value = "No file chosen"
        with out: clear_output()
        status_bar.value = ""

    def _collect_pdb_bytes_and_name():
        # Priority: uploaded file if present, otherwise RCSB using code
        if picked["bytes"] and picked["name"]:
            return picked["bytes"], picked["name"]
        code = pdb_code.value.strip()
        if not _is_valid_pdb_code(code):
            raise ValueError("PDB code must be exactly 4 alphanumeric characters (or upload a file).")
        return _fetch_rcsb(code), f"{code.upper()}.pdb"

    def on_submit(_):
        with out:
            clear_output()
            try:
                chain_global = chain_id.value.strip()
                if not _is_valid_chain(chain_global):
                    raise ValueError("Chain ID must be a single character (e.g., A).")
                if not _is_valid_email(email.value.strip()):
                    raise ValueError("Invalid email format.")

                pdb_bytes, pdb_name = _collect_pdb_bytes_and_name()

                # Save a local copy
                with open(pdb_name, "wb") as f:
                    f.write(pdb_bytes)
                print(f"Saved local copy: {os.getcwd()}/{pdb_name}")

                # OPTIONAL: call mcpath.readpdb_strict.run(...) if available
                ran_local = False
                try:
                    import importlib
                    rmod = importlib.import_module("mcpath.readpdb_strict")
                    runner = None
                    for cand in ("run", "run_readpdb", "readpdb"):
                        if hasattr(rmod, cand):
                            runner = getattr(rmod, cand); break
                    if callable(runner):
                        print("▶ Running readpdb_strict …")
                        res = runner(pdb_name, chain_global)  # expected signature
                        if res is not None:
                            print("readpdb_strict result keys:", list(res.keys()))
                        ran_local = True
                    else:
                        print("ℹ️ readpdb_strict has no callable (run/run_readpdb/readpdb); skipping.")
                except Exception as e:
                    print("ℹ️ readpdb_strict not executed:", e)

                # Build payload (if remote)
                data = {"prediction_mode": pred_type.value, FN["chain_id"]: chain_global}
                if pdb_code.value.strip(): data[FN["pdb_code"]] = pdb_code.value.strip().upper()
                if email.value.strip():    data[FN["email"]] = email.value.strip()

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

                if not target_url:
                    print("\n(No target_url set) — preview only payload below:\n")
                    preview = dict(data); preview["attached_file"] = pdb_name
                    print(preview)
                    if not ran_local:
                        print("\nTip: add mcpath/readpdb_strict.py with run(pdb_path, chain_id) to parse locally.")
                    return

                files = {FN["pdb_file"]: (pdb_name, pdb_bytes, "chemical/x-pdb")}
                print(f"Submitting to {target_url} …")
                r = requests.post(target_url, data=data, files=files, timeout=180)
                print("HTTP", r.status_code)
                try: print("JSON:", r.json())
                except Exception: print("Response (≤800 chars):\n", r.text[:800])

            except Exception as e:
                print("❌", e)

    btn_clear.on_click(on_clear)
    btn_submit.on_click(on_submit)
