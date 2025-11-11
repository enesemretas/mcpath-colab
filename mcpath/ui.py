# mcpath/ui.py — fixed NameError and improved structure

import os, re, requests, yaml, importlib, shutil, sys, traceback
from IPython.display import display
import ipywidgets as W

_FORM_ROOT = None
_LOG_OUT   = None

def _is_valid_pdb_code(c):
    return bool(re.fullmatch(r"[0-9A-Za-z]{4}", (c or "").strip()))

def _is_valid_chain(ch):
    return bool(re.fullmatch(r"[A-Za-z0-9]", (ch or "").strip()))

def _is_valid_email(s):
    s = (s or "").strip()
    return (not s) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s))

def _fetch_rcsb(code: str) -> bytes:
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.content

def _logo_widget(branding: dict):
    url = (branding.get("logo_url") or "").strip()
    if not url:
        return None
    height = int(branding.get("logo_height", 96))
    align  = (branding.get("logo_align") or "center").lower()
    try:
        img_bytes = requests.get(url, timeout=20).content
        jc = {"center": "center", "left": "flex-start", "right": "flex-end"}.get(align, "center")
        return W.HBox(
            [W.Image(value=img_bytes, format="png", layout=W.Layout(height=f"{height}px"))],
            layout=W.Layout(justify_content=jc)
        )
    except Exception:
        return None

def _list_or_custom_row(label: str, options, default_value, minv, maxv, step=1):
    options = sorted({int(x) for x in options})
    default_value = int(default_value)
    if default_value not in options:
        options = sorted(options + [default_value])
    lbl = W.Label(f"{label if label.endswith(':') else label + ':'}", layout=W.Layout(width="180px"))
    toggle = W.ToggleButtons(
        options=[("List", "list"), ("Custom", "custom")], value="list",
        style={'button_width': '120px'},
        layout=W.Layout(display="flex", flex_flow="row", width="260px", margin="0 12px 0 0")
    )
    dropdown = W.Dropdown(options=options, value=default_value,
                          layout=W.Layout(width="280px", min_width="220px"))
    intbox = W.BoundedIntText(value=default_value, min=minv, max=maxv, step=step,
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
    return row, {'get': get_value, 'set_list': set_list, 'set_custom': set_custom}

def launch(defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml", show_title="MCPath-style Parameters"):
    global _FORM_ROOT, _LOG_OUT
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=30).text)
    target_url = (cfg.get("target_url") or "").strip()
    FN = cfg["field_names"]
    SAVE_DIR = cfg.get("save_dir", "/content")
    os.makedirs(SAVE_DIR, exist_ok=True)
    logo = _logo_widget(cfg.get("branding", {}))
    DESC = {'description_width': '180px'}
    wide = W.Layout(width="420px", min_width="360px")
    pdb_code = W.Text(value=str(cfg.get("pdb_code", "")), description="PDB code:", layout=wide, style=DESC)
    or_lbl = W.HTML("<b>&nbsp;&nbsp;or&nbsp;&nbsp;</b>")
    pdb_upload = W.FileUpload(accept=".pdb", multiple=False, description="Choose file")
    file_lbl = W.Label("No file chosen")
    def _on_upload_change(_):
        if pdb_upload.value:
            file_lbl.value = next(iter(pdb_upload.value.keys()))
        else:
            file_lbl.value = "No file chosen"
    pdb_upload.observe(_on_upload_change, names="value")
    chain_id = W.Text(value=str(cfg.get("chain_id", "")), description="Chain ID:", layout=wide, style=DESC)
    email = W.Text(value=str(cfg.get("email", "")), description="Email (opt):", layout=wide, style=DESC)
    pred_type = W.RadioButtons(
        options=[
            ("Functional Residues", "functional"),
            ("Allosteric Paths (initial residue + path length)", "paths_init_len"),
            ("Allosteric Paths (initial & final residues)", "paths_init_final"),
        ],
        value="functional", description="Prediction:", layout=W.Layout(width="auto"), style=DESC)
    big_opts = cfg.get("path_length_options", [100000,200000,300000])
    big_default = int(cfg.get("path_length", big_opts[0]))
    row_big, big_ctrl = _list_or_custom_row("Path length", big_opts, big_default, 1, 1_000_000, 1000)
    get_big_len = big_ctrl['get']
    init_idx = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1, description="Index of initial residue:", layout=wide, style=DESC)
    init_chain = W.Text(value="", description="Chain of initial residue:", placeholder="A", layout=wide, style=DESC)
    row_short, short_ctrl = _list_or_custom_row("Length of Paths", [5,10,20], 5, 1, 10_000, 1)
    get_short_len = short_ctrl['get']
    row_np2, np2_ctrl = _list_or_custom_row("Number of Paths", [1000,2000,3000], 1000, 1, 10_000_000, 100)
    get_num_paths_2 = np2_ctrl['get']
    final_idx = W.BoundedIntText(value=1, min=1, max=1_000_000, step=1, description="Index of final residue:", layout=wide, style=DESC)
    final_chain = W.Text(value="", description="Chain of final residue:", placeholder="B", layout=wide, style=DESC)
    row_np3, np3_ctrl = _list_or_custom_row("Number of Paths", [1000,2000,3000], 1000, 1, 10_000_000, 100)
    get_num_paths_3 = np3_ctrl['get']
    btn_submit = W.Button(description="Submit", button_style="success", icon="paper-plane")
    btn_new_job = W.Button(description="New Job", button_style="info", icon="plus")
    if _LOG_OUT is None:
        _LOG_OUT = W.Output()
    out = _LOG_OUT
    pdb_row = W.HBox([pdb_code, or_lbl, pdb_upload, file_lbl], layout=W.Layout(align_items="center", gap="10px"))
    functional_box = W.VBox([row_big])
    mode2_box = W.VBox([init_idx, init_chain, row_short, row_np2])
    mode3_box = W.VBox([init_idx, init_chain, final_idx, final_chain, row_np3])
    def _sync_mode(*_):
        functional_box.layout.display = ""; mode2_box.layout.display = "none"; mode3_box.layout.display = "none"
        if pred_type.value == "paths_init_len": mode2_box.layout.display = ""
        elif pred_type.value == "paths_init_final": mode3_box.layout.display = ""
    pred_type.observe(_sync_mode, names="value")
    body_children = [logo] if logo else []
    body_children += [W.HTML(f"<h3>{show_title}</h3>"), pdb_row, pred_type, functional_box, mode2_box, mode3_box, chain_id, email, W.HBox([btn_submit, btn_new_job]), out]
    _FORM_ROOT = W.VBox(body_children)
    display(_FORM_ROOT)
    def on_new_job(_):
        pdb_code.value = str(cfg.get("pdb_code", "")); chain_id.value = str(cfg.get("chain_id", "")); email.value = str(cfg.get("email", ""))
        pdb_upload.value = {}; file_lbl.value = "No file chosen"; pred_type.value = "functional"; _LOG_OUT.clear_output(wait=True)
        print("UI reset for new job.")
    def on_submit(_):
        _LOG_OUT.clear_output(wait=True)
        with _LOG_OUT:
            print("Submit clicked — implement request logic here.")
    btn_new_job.on_click(on_new_job)
    btn_submit.on_click(on_submit)
