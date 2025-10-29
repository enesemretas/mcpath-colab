# ui.py
import os, re, requests, shutil, importlib
from IPython.display import display, clear_output
import ipywidgets as W
import yaml

# Try to import local readpdb_strict.py
try:
    import readpdb_strict
except Exception as e:
    raise ImportError("Place readpdb_strict.py in the same folder as ui.py") from e

def _is_valid_pdb_code(c): return bool(re.fullmatch(r"[0-9A-Za-z]{4}", c.strip()))
def _is_valid_chain(ch):   return bool(re.fullmatch(r"[A-Za-z0-9]", ch.strip()))
def _is_valid_email(s):    return (not s.strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s.strip()))

def _fetch_rcsb(code: str) -> bytes:
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    print(f"Fetching PDB from RCSB: {url} (≤15s timeout)")
    r = requests.get(url, timeout=15); r.raise_for_status()
    return r.content

def _list_or_custom(label: str, options, default_value, minv, maxv, step=1, desc_style=None):
    options = sorted({int(x) for x in options})
    default_value = int(default_value)
    if default_value not in options:
        options = sorted(options + [default_value])
    common_layout = W.Layout(width="460px", min_width="420px")
    mode = W.ToggleButtons(
        options=[("List", "list"), ("Custom", "custom")],
        value="list", description="Input:", layout=W.Layout(width="260px"), style=desc_style
    )
    dd = W.Dropdown(
        options=options, value=default_value,
        description=(label if label.endswith(":") else f"{label}:"),
        layout=common_layout, style=desc_style
    )
    txt = W.BoundedIntText(
        value=default_value, min=minv, max=maxv, step=step,
        description=(label if label.endswith(":") else f"{label}:"),
        layout=common_layout, style=desc_style
    )
    container = W.VBox([dd])
    def _on_mode(change): container.children = [dd] if change["new"] == "list" else [txt]
    mode.observe(_on_mode, names="value")
    def get_value(): return int(dd.value if mode.value == "list" else txt.value)
    return mode, container, dd, txt, get_value

def _logo_widget(branding: dict):
    url = (branding.get("logo_url") or "").strip()
    if not url: return None
    height = int(branding.get("logo_height", 96))
    align  = (branding.get("logo_align") or "center").lower()
    try:
        img_bytes = requests.get(url, timeout=10).content
        img = W.Image(value=img_bytes, format="png", layout=W.Layout(height=f"{height}px"))
        jc = {"center":"center", "left":"flex-start", "right":"flex-end"}.get(align, "center")
        return W.HBox([img], layout=W.Layout(justify_content=jc))
    except Exception:
        return None

def _copy_to_free_name(src: str, dst_base: str) -> str | None:
    if not os.path.isfile(src):
        return None
    dst = dst_base; k = 1
    while os.path.exists(dst):
        dst = f"{dst_base}.{k}"; k += 1
    shutil.copy2(src, dst)
    return dst

def launch(
    defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
    show_title="MCPath Parameters (Python readpdb)"
):
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=15).text)
    FN = cfg["field_names"]
    SAVE_DIR = cfg.get("save_dir", "/content")
    os.makedirs(SAVE_DIR, exist_ok=True)

    logo = _logo_widget(cfg.get("branding", {}))
    DESC = {'description_width': '260px'}
    wide = W.Layout(width="460px", min_width="420px")

    # PDB code or upload
    pdb_code   = W.Text(value=str(cfg.get("pdb_code", "")),
                        description="PDB code:", layout=wide, style=DESC)
    or_lbl     = W.HTML("<b>&nbsp;&nbsp;or&nbsp;&nbsp;</b>")
    pdb_upload = W.FileUpload(accept=".pdb", multiple=False, description="Choose File")
    file_lbl   = W.Label("No file chosen")
    def _on_upload_change(change):
        if pdb_upload.value:
            fname = next(iter(pdb_upload.value.keys()))
            file_lbl.value = fname
        else:
            file_lbl.value = "No file chosen"
    pdb_upload.observe(_on_upload_change, names="value")

    chain_id   = W.Text(value=str(cfg.get("chain_id", "")),
                        description="Chain ID:", layout=wide, style=DESC)
    email      = W.Text(value=str(cfg.get("email", "")),
                        description="Email (opt):", layout=wide, style=DESC)

    # Keep path_length UI (compat/logging)
    big_opts = cfg.get("path_length_options", [100000, 200000, 300000])
    big_default = int(cfg.get("path_length", big_opts[0] if big_opts else 100000))
    (pl_mode_big, pl_container_big, pl_dd_big, pl_txt_big, get_big_len) = _list_or_custom(
        label="Path length", options=big_opts, default_value=big_default,
        minv=1, maxv=10_000_000, step=1000, desc_style=DESC
    )
    pred_type = W.RadioButtons(
        options=[("Functional Residues (readpdb only)", "functional")],
        value="functional", description="Prediction:", layout=W.Layout(width="auto"), style=DESC
    )

    btn_submit = W.Button(description="Run (Python readpdb)", button_style="success", icon="play")
    btn_clear  = W.Button(description="Clear",  button_style="warning", icon="trash")
    out        = W.Output()

    pdb_row = W.HBox([pdb_code, W.HTML("&nbsp;"), or_lbl, pdb_upload, file_lbl],
                     layout=W.Layout(align_items="center", justify_content="flex-start",
                                     flex_flow="row wrap", gap="10px"))

    display(W.VBox([
        *( [logo] if logo else [] ),
        W.HTML(f"<h3>{show_title}</h3>"),
        pdb_row,
        W.HTML("<b>Prediction of:</b>"),
        pred_type,
        W.VBox([pl_mode_big, pl_container_big]),
        chain_id, email,
        W.HBox([btn_submit, btn_clear]),
        W.HTML("<hr>"), out
    ], layout=W.Layout(width="auto")))

    def on_clear(_):
        pdb_code.value = ""
        try: pdb_upload.value.clear()
        except Exception: pdb_upload.value = {}
        file_lbl.value = "No file chosen"
        chain_id.value = ""
        email.value = ""
        pl_mode_big.value = "list"
        with out: clear_output()

    def _collect_pdb_bytes():
        if pdb_upload.value:
            (fname, meta) = next(iter(pdb_upload.value.items()))
            return meta["content"], fname
        code = pdb_code.value.strip()
        if not _is_valid_pdb_code(code):
            raise ValueError("PDB code must be exactly 4 alphanumeric characters (or upload a file).")
        return _fetch_rcsb(code), f"{code.upper()}.pdb"

    def on_submit(_):
        with out:
            clear_output()
            try:
                # reload module to pick up edits without kernel restart (handy in Colab)
                importlib.reload(readpdb_strict)

                print("▶ Validating inputs…")
                chain_global = chain_id.value.strip()
                if not _is_valid_chain(chain_global):
                    raise ValueError("Chain ID must be a single character (e.g., A).")
                if not _is_valid_email(email.value.strip()):
                    raise ValueError("Invalid email format.")

                print("▶ Getting PDB (upload or RCSB)…")
                pdb_bytes, pdb_name = _collect_pdb_bytes()

                os.makedirs(SAVE_DIR, exist_ok=True)
                save_path = os.path.join(SAVE_DIR, pdb_name)
                with open(save_path, "wb") as f:
                    f.write(pdb_bytes)
                print(f"Saved local copy: {save_path}")

                # mcpath_input.txt for provenance/logging
                input_path = os.path.join(os.path.dirname(save_path), "mcpath_input.txt")
                rows = ["1", pdb_name, chain_global, str(get_big_len()), (email.value.strip() or "-")]
                with open(input_path, "w") as f:
                    for r in rows:
                        f.write(str(r).strip() + "\n")
                print(f"Input file saved: {input_path}")

                # ---- Call strict reader
                print("▶ Running Python readpdb (strict MATLAB mode)…")
                cor_path, nrows, nprob = readpdb_strict.readpdb_py(save_path, chain_global)

                # Normalize to coor_file
                coor_norm = _copy_to_free_name(cor_path, os.path.join(os.path.dirname(save_path), "coor_file"))

                print("✅ readpdb completed.")
                print(f"Rows written: {nrows}")
                if nprob:
                    print(f"⚠️ problem file written with {nprob} entries.")

                # Preview
                preview_path = coor_norm if coor_norm and os.path.exists(coor_norm) else cor_path
                try:
                    print("\n--- Preview (first 8 lines) ---")
                    with open(preview_path, "r") as f:
                        for i, line in enumerate(f):
                            if i >= 8: break
                            print(line.rstrip("\n"))
                except Exception as e:
                    print("Preview unavailable:", e)

                print("\nOutputs:")
                print(f"  • .cor file : {cor_path}")
                if coor_norm: print(f"  • coor_file : {coor_norm}")
                prob_path = os.path.join(os.path.dirname(save_path), "problem")
                if os.path.exists(prob_path):
                    print(f"  • problem   : {prob_path}")

            except Exception as e:
                print("❌", e)

    btn_clear.on_click(on_clear)
    btn_submit.on_click(on_submit)
