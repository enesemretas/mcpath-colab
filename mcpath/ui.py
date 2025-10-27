# mcpath/ui.py
import os, re, requests, yaml
from IPython.display import display, clear_output
import ipywidgets as W

def _is_valid_pdb_code(c): return bool(re.fullmatch(r"[0-9A-Za-z]{4}", c.strip()))
def _is_valid_chain(ch):   return bool(re.fullmatch(r"[A-Za-z0-9]", ch.strip()))
def _is_valid_email(s):    return (not s.strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s.strip()))

def _fetch_rcsb(code: str) -> bytes:
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    r = requests.get(url, timeout=60); r.raise_for_status()
    return r.content

def launch(
    defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
    show_title="MCPath-style Parameters"
):
    """Launches the MCPath parameter form (UI only)."""
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=30).text)
    target_url = cfg.get("target_url", "").strip()
    FN = cfg["field_names"]

    # --- Widgets: PDB / Chain ---
    pdb_code   = W.Text(value=str(cfg.get("pdb_code", "")), description="PDB code:")
    pdb_upload = W.FileUpload(accept=".pdb", multiple=False, description="Upload PDB")
    chain_id   = W.Text(value=str(cfg.get("chain_id", "")), description="Chain ID:")

    # --- Path length: list OR custom ---
    # Options from YAML (fallback list). Ensure unique, int, sorted.
    opts = cfg.get("path_length_options", [
        100000, 200000, 300000, 400000, 500000, 750000,
        1000000, 2000000, 3000000, 4000000
    ])
    opts = sorted({int(x) for x in opts})
    default_len = int(cfg.get("path_length", opts[0] if opts else 100000))
    if default_len not in opts:
        opts = sorted(opts + [default_len])

    pl_mode = W.ToggleButtons(
        options=[("List", "list"), ("Custom", "custom")],
        value="list" if default_len in opts else "custom",
        description="Input:"
    )
    pl_dropdown = W.Dropdown(
        options=opts, value=default_len if default_len in opts else opts[0],
        description="Path length:"
    )
    pl_custom = W.BoundedIntText(
        value=default_len, min=1, max=10_000_000, step=1000,
        description="Path length:"
    )

    def _sync_visibility(*_):
        if pl_mode.value == "list":
            pl_dropdown.layout.display = ""
            pl_custom.layout.display = "none"
        else:
            pl_dropdown.layout.display = "none"
            pl_custom.layout.display = ""
    _sync_visibility()
    pl_mode.observe(_sync_visibility, names="value")

    def get_path_len() -> int:
        return int(pl_dropdown.value if pl_mode.value == "list" else pl_custom.value)

    # --- Email & actions ---
    email      = W.Text(value=str(cfg.get("email", "")), description="Email (opt):")
    btn_submit = W.Button(description="Submit", button_style="success", icon="paper-plane")
    btn_clear  = W.Button(description="Clear",  button_style="warning", icon="trash")
    out        = W.Output()

    # --- Layout ---
    display(W.VBox([
        W.HTML(f"<h3>{show_title}</h3>"),
        pdb_code, pdb_upload, chain_id,
        pl_mode, pl_dropdown, pl_custom,
        email,
        W.HBox([btn_submit, btn_clear]),
        W.HTML("<hr>"), out
    ]))

    # --- Handlers ---
    def on_clear(_):
        pdb_code.value = ""
        chain_id.value = ""
        email.value = ""
        pdb_upload.value.clear()
        # Reset path length widgets
        pl_mode.value = "list"
        pl_dropdown.value = opts[0]
        pl_custom.value = default_len
        _sync_visibility()
        with out: clear_output()

    def on_submit(_):
        with out:
            clear_output()
            try:
                code  = pdb_code.value.strip()
                chain = chain_id.value.strip()
                mail  = email.value.strip()

                if not _is_valid_chain(chain):
                    raise ValueError("Chain ID must be a single character (e.g., A).")
                if not _is_valid_email(mail):
                    raise ValueError("Invalid email format.")

                have_file = bool(pdb_upload.value)
                have_code = bool(code)
                if not (have_file or have_code):
                    raise ValueError("Provide a PDB code or upload a .pdb file.")

                # Obtain PDB bytes
                if have_file:
                    (fname, meta) = next(iter(pdb_upload.value.items()))
                    pdb_bytes, pdb_name = meta["content"], fname
                else:
                    if not _is_valid_pdb_code(code):
                        raise ValueError("PDB code must be exactly 4 alphanumeric characters.")
                    print("Fetching PDB from RCSB…")
                    pdb_bytes, pdb_name = _fetch_rcsb(code), f"{code.upper()}.pdb"

                chosen_len = get_path_len()
                data = {
                    FN["chain_id"]:   chain,
                    FN["path_length"]: str(chosen_len),
                }
                if code: data[FN["pdb_code"]] = code.upper()
                if mail: data[FN["email"]]    = mail
                files = {FN["pdb_file"]: (pdb_name, pdb_bytes, "chemical/x-pdb")}

                # Save local copy (convenient for users)
                with open(pdb_name, "wb") as f:
                    f.write(pdb_bytes)
                print(f"Saved local copy: {os.getcwd()}/{pdb_name}")

                if not target_url:
                    print("\n(No target_url set) — preview only:")
                    print({**data, "attached_file": pdb_name})
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
