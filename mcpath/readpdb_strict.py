import os, numpy as np

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

def _safe_sub(s,a,b): a=max(1,a); b=min(len(s),b); return s[a-1:b] if a<=b else ""

def _trim_to_first_model(lines):
    has=any(L.startswith("MODEL") for L in lines)
    if not has: return lines
    out=[]; in_m=False
    for L in lines:
        if L.startswith("MODEL"):
            if in_m: break
            in_m=True; continue
        if L.startswith("ENDMDL") and in_m: break
        if not in_m:
            if L.startswith(("ATOM","HETATM","TER")): out.append(L)
        else: out.append(L)
    return out

def _parse_atoms_from_path(pdb_path):
    with open(pdb_path,"r",encoding="utf-8",errors="ignore") as f: lines=f.read().splitlines()
    lines=_trim_to_first_model(lines)
    atoms=[]
    for L in lines:
        if not L.startswith(("ATOM","HETATM")) or len(L)<54: continue
        try:
            atoms.append({
                "AtomSerNo":int(_safe_sub(L,7,11).strip() or "0"),
                "AtomName":_safe_sub(L,13,16).strip(),
                "altLoc":_safe_sub(L,17,17).strip(),
                "resName":_safe_sub(L,18,20).strip(),
                "chainID":_safe_sub(L,22,22).strip(),
                "resSeq":int(_safe_sub(L,23,26).strip() or "0"),
                "iCode":_safe_sub(L,27,27).strip(),
                "X":float(_safe_sub(L,31,38)),
                "Y":float(_safe_sub(L,39,46)),
                "Z":float(_safe_sub(L,47,54)),
                "element":_safe_sub(L,77,78).strip(),
            })
        except: pass
    return atoms

def _normalize_resname(nm): u=(nm or "").upper(); return ALTNAME_NORMALIZE.get(u,u)

def _matlab_h_skip(atom_name,element):
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
        else: buf.append(a)
    if buf: ordered.append((cur,buf))
    res=[]
    for key,group in ordered:
        keep=[True]*len(group)
        for i in range(len(group)-1):
            a=group[i]; b=group[i+1]
            if (a.get("altLoc") or "") and (b.get("altLoc") or "") and (a["altLoc"]!=b["altLoc"]):
                keep[i]=False
        heavy=[]
        for i,a in enumerate(group):
            if not keep[i]: continue
            if _matlab_h_skip(a.get("AtomName",""),a.get("element","")): continue
            heavy.append(a)
        res.append({"key":key,"chain":key[0],"resSeq":key[1],"iCode":key[2],
                    "resName":group[0]["resName"],"atoms":group,"heavy":heavy})
    return res

def _pick_sc_xyz(aa,heavy):
    if len(heavy)<=AA_MIN_ATOMS.get(aa,1e9): return None
    mode,idxs=AA_SC_RULES[aa]; coords=[]
    for idx in idxs:
        j=idx-1
        if 0<=j<len(heavy): a=heavy[j]; coords.append((a["X"],a["Y"],a["Z"]))
        else: return None
    if mode=="index": return coords[0]
    arr=np.array(coords,float); c=arr.mean(0)
    return float(c[0]),float(c[1]),float(c[2])

def readpdb_py(pdb_path,chainID):
    chainID=(chainID or "").strip()[:1].upper()
    atoms_all=_parse_atoms_from_path(pdb_path)
    ranks=_chain_rank_map_stream(atoms_all)
    atoms_use=[a for a in atoms_all if not chainID or (a["chainID"] or "").upper()==chainID] or atoms_all
    residues=_build_residues_stream(atoms_use)
    rows=[]; problems=[]
    for res in residues:
        rn=res["resName"]
        if rn not in AA_CODE: continue
        heavy=res["heavy"]
        ca=next((a for a in heavy if a["AtomName"].upper()=="CA"),None)\
            or next((a for a in res["atoms"] if a["AtomName"].upper()=="CA"),None)
        if not ca: continue
        aa=AA_CODE[rn]; chrank=float(ranks.get(res["chain"],0)); cnt=float(len(heavy))
        sc=_pick_sc_xyz(aa,heavy)
        scx,scy,scz=(float(ca["X"]),float(ca["Y"]),float(ca["Z"])) if sc is None else sc
        rows.append([float(res["resSeq"]),float(aa),
                     float(ca["X"]),float(ca["Y"]),float(ca["Z"]),
                     scx,scy,scz,cnt,chrank,float(AA_CONST[aa])])
        if cnt<=AA_MIN_ATOMS.get(aa,0): problems.append((res["resSeq"],int(chrank)))
    cor_path=f"{pdb_path}.cor"
    with open(cor_path,"w",encoding="utf-8") as f:
        for r in rows:
            f.write(f"{r[0]:5.0f} {r[1]:4.0f} {r[2]:8.3f} {r[3]:8.3f} {r[4]:8.3f} "
                    f"{r[5]:8.3f} {r[6]:8.3f} {r[7]:8.3f} {r[8]:4.0f} {r[9]:2.0f} {r[10]:4.1f}\n")
    if problems:
        with open(os.path.join(os.path.dirname(pdb_path),"problem"),"w",encoding="utf-8") as f:
            f.write("There are problems about number of atoms on following residues\n")
            for resseq,cr in problems: f.write(f"{int(resseq)} of {int(cr)} chain\n")
    return cor_path,len(rows),len(problems)
