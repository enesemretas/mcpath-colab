# mcpath/ui.py
import os, re, io, requests, yaml
from IPython.display import display, HTML, clear_output
import ipywidgets as W

def _is_valid_pdb_code(c): return bool(re.fullmatch(r"[0-9A-Za-z]{4}", c.strip()))
def _is_valid_chain(ch):   return bool(re.fullmatch(r"[A-Za-z0-9]", ch.strip()))
def _is_valid_email(s):    return (not s.strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s.strip()))

def _fetch_rcsb(code):
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.content

def launch(defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
           show_title="MCPath-style Parameters"):
    """Launches the MCPath parameter form."""
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=30).text)
    target_url = cfg.get("target_url","").strip()
    FN = cfg["field_names"]

    pdb_code   = W.Text(value=str(cfg.get("pdb_code","")), description="PDB code:")
    pdb_upload = W.FileUpload(accept=".pdb", multiple=False, description="Upload PDB")
    chain_id   = W.Text(value=str(cfg.get("chain_id","")), description="Chain ID:")
    # Read options from YAML (falls back to a sane list if absent)
    opts = cfg.get("path_length_options", [
    100000, 200000, 300000, 400000, 500000, 750000,
    1000000, 2000000, 3000000, 4000000])
    # make sure all are ints and unique, sorted
    opts = sorted({int(x) for x in opts})

    default_len = int(cfg.get("path_length", opts[0] if opts else 100000))
    # if the default isn't in the list, include it and re-sort so it can be selected
    if default_len not in opts:
    opts = sorted(opts + [default_len])

    path_len = W.Dropdown(options=opts, value=default_len, description="Path length:")


    email      = W.Text(value=str(cfg.get("email","")), description="Email (opt):")
    btn_submit = W.Button(description="Submit", button_style="success", icon="paper-plane")
    btn_clear  = W.Button(description="Clear",  button_style="warning", icon="trash")
    out        = W.Output()

    display(W.VBox([
        W.HTML(f"<h3>{show_title}</h3>"),
        pdb_code, pdb_upload, chain_id, path_len, email,
        W.HBox([btn_submit, btn_clear]),
        W.HTML("<hr>"), out
    ]))

    def on_clear(_):
        pdb_code.value=""; chain_id.value=""; path_len.value=100_000; email.value=""; pdb_upload.value.clear()
        with out: clear_output()

    def on_submit(_):
        with out:
            clear_output()
            try:
                code=pdb_code.value.strip(); chain=chain_id.value.strip(); mail=email.value.strip()
                if not _is_valid_chain(chain): raise ValueError("Chain ID must be a single character.")
                if not _is_valid_email(mail):  raise ValueError("Invalid email format.")
                have_file,have_code=bool(pdb_upload.value),bool(code)
                if not (have_file or have_code): raise ValueError("Provide a PDB code or upload a file.")

                if have_file:
                    (fname,meta)=next(iter(pdb_upload.value.items()))
                    pdb_bytes, pdb_name=meta["content"], fname
                else:
                    if not _is_valid_pdb_code(code): raise ValueError("PDB code must be 4 chars.")
                    print("Fetching PDB from RCSB…")
                    pdb_bytes, pdb_name=_fetch_rcsb(code), f"{code.upper()}.pdb"

                data={FN["chain_id"]:chain, FN["path_length"]:str(path_len.value)}
                if code: data[FN["pdb_code"]]=code.upper()
                if mail: data[FN["email"]]=mail
                files={FN["pdb_file"]:(pdb_name,pdb_bytes,"chemical/x-pdb")}

                with open(pdb_name,"wb") as f: f.write(pdb_bytes)
                print(f"Saved local copy: {os.getcwd()}/{pdb_name}")

                if not target_url:
                    print("\n(No target_url set) — preview only:\n", {**data,"attached_file":pdb_name})
                    return

                print(f"Submitting to {target_url} …")
                r=requests.post(target_url,data=data,files=files,timeout=180)
                print("HTTP",r.status_code)
                try: print("JSON:",r.json())
                except Exception: print("Response (≤800 chars):\n",r.text[:800])
            except Exception as e:
                print("❌",e)

    btn_clear.on_click(on_clear)
    btn_submit.on_click(on_submit)
