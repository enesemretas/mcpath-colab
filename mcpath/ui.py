import os, re, requests, shutil, importlib, yaml
from IPython.display import display, clear_output
import ipywidgets as W

# import the new reader
from mcpath import readpdb_strict

def _is_valid_pdb_code(c): return bool(re.fullmatch(r"[0-9A-Za-z]{4}", c.strip()))
def _is_valid_chain(ch):   return bool(re.fullmatch(r"[A-Za-z0-9]", ch.strip()))
def _is_valid_email(s):    return (not s.strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s.strip()))

def _fetch_rcsb(code):
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    print(f"Fetching PDB: {url}")
    r = requests.get(url, timeout=15); r.raise_for_status()
    return r.content

def launch(defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
           show_title="MCPath (Python readpdb)"):
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=15).text)
    SAVE_DIR = cfg.get("save_dir", "/content")
    os.makedirs(SAVE_DIR, exist_ok=True)

    DESC={'description_width':'260px'}
    wide=W.Layout(width="460px")

    pdb_code=W.Text(value=str(cfg.get("pdb_code","")),description="PDB code:",layout=wide,style=DESC)
    pdb_up=W.FileUpload(accept=".pdb",multiple=False,description="Upload PDB")
    file_lbl=W.Label("No file chosen")
    def _on_up(change):
        file_lbl.value=next(iter(pdb_up.value.keys())) if pdb_up.value else "No file chosen"
    pdb_up.observe(_on_up,names="value")

    chain_id=W.Text(value=str(cfg.get("chain_id","")),description="Chain ID:",layout=wide,style=DESC)
    email=W.Text(value=str(cfg.get("email","")),description="Email (opt):",layout=wide,style=DESC)

    btn=W.Button(description="Run readpdb (Python)",button_style="success",icon="play")
    out=W.Output()

    display(W.VBox([W.HTML(f"<h3>{show_title}</h3>"),
                    W.HBox([pdb_code,pdb_up,file_lbl]),
                    chain_id,email,btn,W.HTML("<hr>"),out]))

    def _collect():
        if pdb_up.value:
            (fname,meta)=next(iter(pdb_up.value.items()))
            return meta["content"],fname
        c=pdb_code.value.strip()
        if not _is_valid_pdb_code(c): raise ValueError("Invalid PDB code.")
        return _fetch_rcsb(c),f"{c.upper()}.pdb"

    def on_run(_):
        with out:
            clear_output()
            try:
                pdb_bytes,pdb_name=_collect()
                save_path=os.path.join(SAVE_DIR,pdb_name)
                with open(save_path,"wb") as f: f.write(pdb_bytes)
                print("Saved:",save_path)
                print("▶ Running Python readpdb...")
                cor,nrow,nprob=readpdb_strict.readpdb_py(save_path,chain_id.value)
                print("✅ Done —",nrow,"rows,",nprob,"problems")
                with open(cor) as f:
                    for i,l in enumerate(f):
                        if i>=8: break
                        print(l.rstrip())
            except Exception as e:
                print("❌",e)
    btn.on_click(on_run)
