# @title MCPath (Python readpdb only) — install & launch
!pip -q install pyyaml ipywidgets requests

import sys, shutil, pathlib, subprocess, textwrap, os

REPO_URL = "https://github.com/enesemretas/mcpath-colab.git"
REPO_DIR = pathlib.Path("/content/mcpath-colab")
PKG_DIR  = REPO_DIR / "mcpath"

# Fresh clone
if REPO_DIR.exists():
    shutil.rmtree(REPO_DIR)
subprocess.run(["git", "clone", "--depth", "1", REPO_URL, str(REPO_DIR)], check=True)
PKG_DIR.mkdir(parents=True, exist_ok=True)

# --- write readpdb_strict.py ---
readpdb_strict_py = r"""
import numpy as np

AA_CODE = {
    "GLY":1,"ALA":2,"VAL":3,"ILE":4,"LEU":5,"SER":6,"THR":7,"ASP":8,"ASN":9,"GLU":10,
    "GLN":11,"LYS":12,"ARG":13,"CYS":14,"MET":15,"PHE":16,"TYR":17,"TRP":18,"HIS":19,"PRO":20
}
AA_CONST = {1:75.0,2:89.1,3:117.2,4:131.2,5:131.2,6:105.1,7:119.1,8:133.1,9:132.1,10:147.1,
            11:146.0,12:146.2,13:174.2,14:121.2,15:149.2,16:165.2,17:181.2,18:204.2,19:155.2,20:115.1}
AA_MIN_ATOMS = {1:1,2:4,3:6,4:7,5:7,6:5,7:5,8:7,9:7,10:8,11:8,12:8,13:10,14:5,15:6,16:10,17:11,18:13,19:9,20:6}
AA_SC_RULES = {
 1:("index",[2]), 2:("index",[5]), 3:("mean",[6,7]), 4:("index",[8]), 5:("mean",[7,8]),
 6:("index",[6]), 7:("index",[6]), 8:("mean",[7,8]), 9:("mean",[7,8]), 10:("mean",[8,9]),
 11:("mean",[8,9]), 12:("index",[9]), 13:("mean",[8,10,11]), 14:("index",[6]), 15:("index",[7]),
 16:("mean",[6,7,8,9,10,11]), 17:("mean",[6,7,8,9,10,11,12]), 18:("mean",[7,8,9,10,11,12,13,14]),
 19:("mean",[6,7,8,9,10]), 20:("mean",[5,6,7]),
}
ALTNAME_NORMALIZE = {"HIP":"HIS","HIE":"HIS","HID":"HIS","HSD":"HIS","HSE":"HIS","HSP":"HIS",
                     "LYN":"LYS","GLH":"GLU","CYM":"CYS","CYX":"CYS","ASH":"ASN","TYM":"TYR"}

def _safe_sub(s,a,b):
    a = max(1,a); b = min(len(s),b)
    return s[a-1:b] if a<=b else ""

def _trim_to_first_model(lines):
    has = any(L.startswith("MODEL") for L in lines)
    if not has: return lines
    out = []; in_m = False
    for L in lines:
        if L.startswith("MODEL"):
            if in_m: break
            in_m = True; continue
        if L.startswith("ENDMDL") and in_m: break
        if not in_m:
            if L.startswith(("ATOM","HETATM","TER")): out.append(L)
        else:
            out.append(L)
    return out

def _parse_atoms_from_path(pdb_path):
    with open(pdb_path,"r",encoding="utf-8",errors="ignore") as f:
        lines=f.read().splitlines()
    lines=_trim_to_first_model(lines)
    atoms=[]
    for L in lines:
        if not L.startswith(("ATOM","HETATM")): continue
        if len(L)<54: continue
        try:
            atoms.append({
                "AtomSerNo":int(_safe_sub(L,7,11).strip() or "0"),
                "AtomName":_safe_sub(L,13,16).strip(),
                "altLoc":_safe_sub(L,17,17).strip(),
                "resName":_safe_sub(L,18,20).strip(),
                "chainID":_safe_sub(L,22,22).strip(),
                "resSeq":int(_safe_sub(L,23,26).strip() or "0"),
                "iCode":_safe_sub(L,27,27).strip(),
                "X":float(_safe_sub(L,31,38)), "Y":float(_safe_sub(L,39,46)), "Z":float(_safe_sub(L,47,54)),
                "element":_safe_sub(L,77,78).strip(),
            })
        except: pass
    return atoms

def _normalize_resname(nm): 
    u=(nm or "").upper(); 
    return ALTNAME_NORMALIZE.get(u,u)

def _matlab_h_skip(atom_name, element):
    nm=(atom_name or "").upper(); el=(element or "").upper()
    if el=="H" or nm in {"H","D","E","G","HN"}: return True
    if ("H" in nm) and ("N" not in nm): return True
    return False

def _chain_rank_map_stream(atoms):
    ranks={}; r=0; last=None
    for a in atoms:
        ch=(a["chainID"] or "").upper()
        if ch!=last:
            if ch not in ranks: r+=1; ranks[ch]=r
            last=ch
    return ranks

def _build_residues_stream(atoms):
    for a in atoms:
        a["resName"]=_normalize_resname(a["resName"])
        a["chainID"]=(a["chainID"] or "").upper()
    ordered=[]; cur=None; buf=[]
    for a in atoms:
        key=(a["chainID"],a["resSeq"],a["iCode"])
        if key!=cur:
            if buf: ordered.append((cur,buf))
            cur=key; buf=[a]
        else:
            buf.append(a)
    if buf: ordered.append((cur,buf))
    res_list=[]
    for key, group in ordered:
        keep=[True]*len(group)
        for i in range(len(group)-1):
            a=group[i]; b=group[i+1]
            if (a.get("altLoc") or "") and (b.get("altLoc") or "") and (a["altLoc"]!=b["altLoc"]):
                keep[i]=False
        heavy=[]
        for i,a in enumerate(group):
            if not keep[i]: continue
            if _matlab_h_skip(a.get("AtomName",""), a.get("element","")): continue
            heavy.append(a)
        res_list.append({
            "key":key,"chain":key[0],"resSeq":key[1],"iCode":key[2],
            "resName":group[0]["resName"] if group else "", "atoms":group, "heavy":heavy
        })
    return res_list

def _pick_sc_xyz(aa_code, heavy):
    min_need=AA_MIN_ATOMS.get(aa_code,10**9)
    if len(heavy)<=min_need: return None
    mode, idxs = AA_SC_RULES[aa_code]
    coords=[]
    for idx in idxs:
        j=idx-1
        if 0<=j<len(heavy):
            a=heavy[j]; coords.append((a["X"],a["Y"],a["Z"]))
        else:
            return None
    if mode=="index": return coords[0]
    arr=np.array(coords,float); c=arr.mean(0)
    return float(c[0]), float(c[1]), float(c[2])

def readpdb_py(pdb_path, chainID):
    chainID=(chainID or "").strip()[:1].upper()
    atoms_all=_parse_atoms_from_path(pdb_path)
    ranks_full=_chain_rank_map_stream(atoms_all)
    atoms_use=[a for a in atoms_all if not chainID or (a["chainID"] or "").upper()==chainID] or atoms_all
    residues=_build_residues_stream(atoms_use)

    rows=[]; problems=[]
    for res in residues:
        rn=res["resName"]
        if rn not in AA_CODE: continue
        heavy=res["heavy"]
        ca=next((a for a in heavy if (a.get("AtomName","").upper()=="CA")), None) \
           or next((a for a in res["atoms"] if (a.get("AtomName","").upper()=="CA")), None)
        if ca is None: continue
        aa=AA_CODE[rn]; chain_rank=float(ranks_full.get(res["chain"],0))
        cnt=float(len(heavy))
        sc=_pick_sc_xyz(aa, heavy)
        if sc is None: scx,scy,scz = float(ca["X"]),float(ca["Y"]),float(ca["Z"])
        else: scx,scy,scz = sc
        const=AA_CONST[aa]
        rows.append([
            float(res["resSeq"]), float(aa),
            float(ca["X"]),float(ca["Y"]),float(ca["Z"]),
            float(scx),float(scy),float(scz),
            cnt, chain_rank, float(const)
        ])
        if cnt <= AA_MIN_ATOMS.get(aa,0):
            problems.append((res["resSeq"], int(chain_rank)))

    cor_path=f"{pdb_path}.cor"
    with open(cor_path,"w",encoding="utf-8") as f:
        for r in rows:
            f.write(f"{r[0]:5.0f} {r[1]:4.0f} {r[2]:8.3f} {r[3]:8.3f} {r[4]:8.3f} "
                    f"{r[5]:8.3f} {r[6]:8.3f} {r[7]:8.3f} {r[8]:4.0f} {r[9]:2.0f} {r[10]:4.1f}\n")
    if problems:
        with open(os.path.join(os.path.dirname(pdb_path),"problem"),"w",encoding="utf-8") as f:
            f.write("There are problems about number of atoms on following residues\n")
            for resseq,crank in problems:
                f.write(f"{int(resseq)} of {int(crank)} chain\n")
    return cor_path, len(rows), len(problems)
"""
with open(PKG_DIR / "readpdb_strict.py", "w", encoding="utf-8") as f:
    f.write(readpdb_strict_py)

# --- write a minimal ui.py that calls readpdb_strict ---
ui_py = r"""
import os, re, requests, shutil, importlib, yaml
from IPython.display import display, clear_output
import ipywidgets as W

try:
    import readpdb_strict
except Exception as e:
    raise ImportError("readpdb_strict.py missing next to ui.py") from e

def _is_valid_pdb_code(c): return bool(re.fullmatch(r"[0-9A-Za-z]{4}", c.strip()))
def _is_valid_chain(ch):   return bool(re.fullmatch(r"[A-Za-z0-9]", ch.strip()))
def _is_valid_email(s):    return (not s.strip()) or bool(re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", s.strip()))

def _fetch_rcsb(code: str) -> bytes:
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    print(f"Fetching PDB from RCSB: {url} (≤15s timeout)")
    r = requests.get(url, timeout=15); r.raise_for_status()
    return r.content

def _logo_widget(branding: dict):
    url = (branding.get("logo_url") or "").strip()
    if not url: return None
    height = int(branding.get("logo_height", 96))
    align  = (branding.get("logo_align") or "center").lower()
    try:
        img_bytes = requests.get(url, timeout=10).content
        img = W.Image(value=img_bytes, format="png", layout=W.Layout(height=f"{height}px"))
        jc = {"center":"center","left":"flex-start","right":"flex-end"}.get(align,"center")
        return W.HBox([img], layout=W.Layout(justify_content=jc))
    except Exception:
        return None

def _copy_to_free_name(src: str, dst_base: str) -> str | None:
    if not os.path.isfile(src): return None
    dst = dst_base; k = 1
    while os.path.exists(dst):
        dst = f"{dst_base}.{k}"; k += 1
    shutil.copy2(src, dst); return dst

def launch(defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
           show_title="MCPath (Python readpdb)"):
    cfg = yaml.safe_load(requests.get(defaults_url, timeout=15).text)
    SAVE_DIR = cfg.get("save_dir", "/content")
    os.makedirs(SAVE_DIR, exist_ok=True)
    logo = _logo_widget(cfg.get("branding", {}))
    DESC = {'description_width':'260px'}
    wide = W.Layout(width="460px", min_width="420px")

    pdb_code = W.Text(value=str(cfg.get("pdb_code","")), description="PDB code:", layout=wide, style=DESC)
    or_lbl   = W.HTML("<b>&nbsp;&nbsp;or&nbsp;&nbsp;</b>")
    pdb_up   = W.FileUpload(accept=".pdb", multiple=False, description="Choose File")
    file_lbl = W.Label("No file chosen")
    def _on_up(change):
        file_lbl.value = next(iter(pdb_up.value.keys())) if pdb_up.value else "No file chosen"
    pdb_up.observe(_on_up, names="value")

    chain_id = W.Text(value=str(cfg.get("chain_id","")), description="Chain ID:", layout=wide, style=DESC)
    email    = W.Text(value=str(cfg.get("email","")), description="Email (opt):", layout=wide, style=DESC)
    path_len = W.IntText(value=int(cfg.get("path_length",100000)), description="Path length (for log):",
                         layout=wide, style=DESC)

    btn_run  = W.Button(description="Run (Python readpdb)", button_style="success", icon="play")
    btn_clr  = W.Button(description="Clear", button_style="warning", icon="trash")
    out      = W.Output()

    row = W.HBox([pdb_code, W.HTML("&nbsp;"), or_lbl, pdb_up, file_lbl],
                 layout=W.Layout(align_items="center", justify_content="flex-start",
                                 flex_flow="row wrap", gap="10px"))

    display(W.VBox([
        *( [logo] if logo else []),
        W.HTML(f"<h3>{show_title}</h3>"),
        row, chain_id, email, path_len,
        W.HBox([btn_run, btn_clr]), W.HTML("<hr>"), out
    ]))

    def _collect():
        if pdb_up.value:
            (fname, meta) = next(iter(pdb_up.value.items()))
            return meta["content"], fname
        code = pdb_code.value.strip()
        if not _is_valid_pdb_code(code):
            raise ValueError("PDB code must be exactly 4 alphanumeric characters (or upload a file).")
        return _fetch_rcsb(code), f"{code.upper()}.pdb"

    def on_clear(_):
        pdb_code.value=""; chain_id.value=""; email.value=""; path_len.value=100000
        try: pdb_up.value.clear()
        except Exception: pdb_up.value={}
        file_lbl.value="No file chosen"
        with out: clear_output()

    def on_run(_):
        with out:
            clear_output()
            try:
                importlib.reload(readpdb_strict)
                ch = chain_id.value.strip()
                if not _is_valid_chain(ch): raise ValueError("Chain ID must be a single character (e.g., A).")
                if not _is_valid_email(email.value.strip()): raise ValueError("Invalid email format.")
                pdb_bytes, pdb_name = _collect()

                save_path = os.path.join(SAVE_DIR, pdb_name)
                with open(save_path,"wb") as f: f.write(pdb_bytes)
                print("Saved:", save_path)

                # provenance file (like original)
                with open(os.path.join(SAVE_DIR,"mcpath_input.txt"),"w") as f:
                    for r in ["1", pdb_name, ch, str(int(path_len.value)), (email.value.strip() or "-")]:
                        f.write(str(r).strip()+"\n")

                print("▶ Running Python readpdb (strict mode)…")
                cor_path, nrows, nprob = readpdb_strict.readpdb_py(save_path, ch)

                coor_norm = _copy_to_free_name(cor_path, os.path.join(SAVE_DIR,"coor_file"))
                print("✅ Done. Rows:", nrows)
                if nprob: print(f"⚠️ problem file with {nprob} entries.")

                prev = coor_norm if coor_norm and os.path.exists(coor_norm) else cor_path
                print("\n--- Preview (first 8 lines) ---")
                try:
                    with open(prev,"r") as f:
                        for i,line in enumerate(f):
                            if i>=8: break
                            print(line.rstrip())
                except Exception as e:
                    print("Preview unavailable:", e)

                print("\nOutputs:")
                print("  •", cor_path)
                if coor_norm: print("  •", coor_norm)
                prob = os.path.join(SAVE_DIR,"problem")
                if os.path.exists(prob): print("  •", prob)
            except Exception as e:
                print("❌", e)

    btn_run.on_click(on_run)
    btn_clr.on_click(on_clear)
"""
with open(PKG_DIR / "ui.py", "w", encoding="utf-8") as f:
    f.write(ui_py)

# Import & launch
sys.path.append(str(REPO_DIR))
from mcpath.ui import launch
from google.colab import output as colab_output
colab_output.enable_custom_widget_manager()
launch(
    defaults_url="https://raw.githubusercontent.com/enesemretas/mcpath-colab/main/config/defaults.yaml",
    show_title="MCPath (Python readpdb)"
)
