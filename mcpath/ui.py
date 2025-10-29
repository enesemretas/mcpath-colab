# mcpath/ui.py
import os, re, requests, yaml, importlib
from IPython.display import display, clear_output
import ipywidgets as W

# Prefer Colab's picker; fall back gracefully
try:
    from google.colab import files as colab_files
    _HAS_COLAB = True
except Exception:
    _HAS_COLAB = False

# -------------------- validators & helpers --------------------
def _is_valid_pdb_code(c): return bool(re.fullmatch(r"[0-9A-Za-z]{4}", c.strip()))
def _is_valid_chain(ch):   return bool(re.fullmatch(r"[A-Za-z0-9]", ch.strip()))
def _is_valid_email(s):    return (not s.strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s.strip()))

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
    run_python_readpdb = bool(cfg.get("run_python_readpdb_after_submit", True))

    # ---------- Branding ----------
    logo = _logo_widget(cfg.get("branding", {}))

    # ---------- PDB: code OR single-button upload ----------
    pdb_code   = W.Text(value=str(cfg.get("pdb_code", "")), description="PDB code:")
    or_lbl     = W.HTML("<b>&nbsp;&nbsp;or&nbsp;&nbsp;</b>")
    btn_pick   = W.Button(description="Choose file", icon="upload")
    file_lbl   = W.Label("No file chosen")

    # Store chosen file (name, bytes)
    picked = {"name": None, "bytes": None}

    def _on_pick(_):
        if not _HAS_COLAB:
            file_lbl.value = "Colab picker not available in this environment."
            return
        try:
            result = colab_files.upload()
            if not result:
                return
            name, data = next(iter(result.items()))
            picked["name"]  = name
            picked["bytes"] = data
            file_lbl.value  = f"{name}"
        except Exception as e:
            file_lbl.value = f"Upload failed: {e}"

    btn_pick.on_click(_on_pick)

    # ---------- Always fields ----------
    chain_id   = W.Text(value=str(cfg.get("chain_id", "A")), description="Chain ID:")
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

    # ---------- Functional mode ----------
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

    # ---------- Actions ----------
    btn_submit = W.Button(description="Submit", button_style="success", icon="paper-plane")
    btn_clear  = W.Button(description="Clear",  button_style="warning", icon="trash")
    out        = W.Output()

    # ---------- Layout ----------
    pdb_row = W.HBox([pdb_code, W.HTML("&nbsp;"), or_lbl, btn_pick, file_lbl])
    functional_box = W.VBox([pl_mode_big, pl_dd_big, pl_txt_big])
    mode2_box = W.VBox([init_idx, init_chain, pl_mode_short, pl_dd_short, pl_txt_short, np_mode_2, np_dd_2, np_txt_2])
    mode3_box = W.VBox([init_idx, init_chain, final_idx, final_chain, np_mode_3, np_dd_3, np_txt_3])

    def _sync_mode(*_):
        functional_box.layout.display = ""
        mode2_box.layout.display = "none"
        mode3_box.layout.display = "none"
        if pred_type.value == "paths_init_len":
            functional_box.layout.display = "none"; mode2_box.layout.display = ""
        elif pred_type.value == "paths_init_final":
            functional_box.layout.display = "none"; mode3_box.layout.display = ""
    _sync_mode()
    pred_type.observe(_sync_mode, names="value")

    children = []
    if logo: children.append(logo)
    children += [
        W.HTML(f"<h3>{show_title}</h3>"),
        pdb_row,
        W.HTML("<b>Prediction of:</b>"),
        pred_type,
        functional_box, mode2_box, mode3_box,
        chain_id, email,
        W.HBox([btn_submit, btn_clear]),
        W.HTML("<hr>"), out
    ]
    display(W.VBox(children))

    # -------------------- handlers --------------------
    def on_clear(_):
        pdb_code.value = ""
        picked["name"] = None; picked["bytes"] = None
        file_lbl.value = "No file chosen"
        chain_id.value = str(cfg.get("chain_id","A"))
        email.value = ""
        pred_type.value = "functional"
        pl_mode_big.value = "list"; pl_dd_big.value = pl_dd_big.options[0]; pl_txt_big.value = int(pl_dd_big.options[0])
        init_idx.value = 1; init_chain.value = ""
        pl_mode_short.value = "list"; pl_dd_short.value = pl_dd_short.options[0]; pl_txt_short.value = int(pl_dd_short.options[0])
        np_mode_2.value = "list"; np_dd_2.value = np_dd_2.options[0]; np_txt_2.value = int(np_dd_2.options[0])
        final_idx.value = 1; final_chain.value = ""
        np_mode_3.value = "list"; np_dd_3.value = np_dd_3.options[0]; np_txt_3.value = int(np_dd_3.options[0])
        with out: clear_output()

    def _collect_pdb_bytes():
        if picked["bytes"] is not None:
            return picked["bytes"], picked["name"]
        code = pdb_code.value.strip()
        if not _is_valid_pdb_code(code):
            raise ValueError("PDB code must be exactly 4 alphanumeric characters (or upload a file).")
        return _fetch_rcsb(code), f"{code.upper()}.pdb"

    def _run_python_readpdb(saved_path: str, chain: str):
        try:
            from mcpath import readpdb_strict as rps
            importlib.reload(rps)
        except Exception as e:
            print("ℹ️ readpdb_strict could not be imported:", e)
            return None
        if not hasattr(rps, "run") or not callable(rps.run):
            print("ℹ️ readpdb_strict has no callable (run); skipping.")
            return None
        try:
            out_path = rps.run(saved_path, chain, strict_matlab_mode=True)
            return out_path
        except Exception as e:
            print("❌ readpdb_strict.run failed:", e)
            return None

    def on_submit(_):
        with out:
            clear_output()
            try:
                chain_global = chain_id.value.strip()
                if not _is_valid_chain(chain_global):
                    raise ValueError("Chain ID must be a single character (e.g., A).")
                if not _is_valid_email(email.value.strip()):
                    raise ValueError("Invalid email format.")

                pdb_bytes, pdb_name = _collect_pdb_bytes()
                save_path = os.path.join(os.getcwd(), pdb_name)
                with open(save_path, "wb") as f:
                    f.write(pdb_bytes)
                print(f"Saved local copy: {save_path}")

                cor_path = None
                if run_python_readpdb:
                    cor_path = _run_python_readpdb(save_path, chain_global)
                    if cor_path and os.path.exists(cor_path):
                        print(f"✅ Created: {cor_path}")
                        try:
                            with open(cor_path, "r") as cf:
                                head = "".join([next(cf) for _ in range(5)])
                            print("\n--- .cor (first 5 lines) ---\n" + head)
                        except Exception:
                            pass

                data = {"prediction_mode": pred_type.value, FN["chain_id"]: chain_global}
                if pdb_code.value.strip():
                    data[FN["pdb_code"]] = pdb_code.value.strip().upper()
                if email.value.strip():
                    data[FN["email"]] = email.value.strip()
                files = {FN["pdb_file"]: (pdb_name, pdb_bytes, "chemical/x-pdb")}

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
