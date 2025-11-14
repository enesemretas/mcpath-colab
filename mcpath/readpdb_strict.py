# readpdb_strict.py
"""
Strict-MATLAB-mode PDB → .cor converter (pure Python) with mcpath_input.txt support.

Primary entry:
    run(pdb_path=None, chain=None, out_basename=None, strict_matlab_mode=True, input_path=None)

Behavior:
- If mcpath_input.txt exists (path given or auto-located), it is parsed and used:
    mode=1 : functional → lines: [1, pdb_name, chain, path_len, email]
    mode=2 : paths_init_len → [2, pdb_name, chain, L, N, idx_init, chain_init, email]
    mode=3 : paths_init_final → [3, pdb_name, chain, idx_i, chain_i, idx_f, chain_f, email]
  Only (pdb_name, chain) are needed for .cor creation; others are ignored here.
- If no input file is found, falls back to direct arguments (pdb_path, chain).

Chain semantics (updated):
- chain may be:
    "" (empty)         → use ALL chains.
    "A"                → only chain A (by rank).
    "A,B" or "A,B,C"   → only those chains (by rank).
- If the PDB has no chain IDs at all, ALL residues are used regardless of chain.
- Writes "<pdb_name>.cor" next to the PDB and returns the .cor path.
"""

from __future__ import annotations
import os, re
from collections import defaultdict
from typing import Optional, Tuple, Dict, List

# ---------------- Constants & remaps (MATLAB-compatible) ----------------
AA_SET = {
    "GLY":1, "ALA":2, "VAL":3, "ILE":4, "LEU":5, "SER":6, "THR":7, "ASP":8,
    "ASN":9, "GLU":10, "GLN":11, "LYS":12, "ARG":13, "CYS":14, "MET":15,
    "PHE":16, "TYR":17, "TRP":18, "HIS":19, "PRO":20
}
RESNAME_REMAP = {
    "HIP":"HIS","HIE":"HIS","HID":"HIS","HSD":"HIS","HSE":"HIS","HSP":"HIS",
    "LYN":"LYS","GLH":"GLU","CYM":"CYS","CYX":"CYS","ASH":"ASN","TYM":"TYR"
}
RES_MASS = {  # same numbers you use in MATLAB output column 11
    1:75,   2:89.1, 3:117.2,4:131.2,5:131.2,6:105.1,7:119.1,8:133.1,9:132.1,10:147.1,
    11:146, 12:146.2,13:174.2,14:121.2,15:149.2,16:165.2,17:181.2,18:204.2,19:155.2,20:115.1
}

# ---------------- Small helpers ----------------
def _safesub(s: str, a: int, b: int) -> str:
    a = max(1, a); b = min(len(s), b)
    return "" if a > b else s[a-1:b]

def _trim_to_first_model(pdb_text: List[str]) -> List[str]:
    has_model = any(l.startswith("MODEL") for l in pdb_text)
    if not has_model:
        return pdb_text
    out, in_model, seen_first = [], False, False
    for L in pdb_text:
        if L.startswith("MODEL"):
            if seen_first: break
            in_model = True; seen_first = True
            continue
        if L.startswith("ENDMDL") and in_model:
            break
        if not in_model:
            if L.startswith(("ATOM  ","HETATM","TER")):
                out.append(L)
        else:
            out.append(L)
    return out

def _parse_pdb_atoms(lines: List[str]) -> List[Dict]:
    atoms = []
    for L in lines:
        if not (L.startswith("ATOM") or L.startswith("HETATM")):
            continue
        AtomSerNo = _safesub(L, 7, 11).strip()
        AtomName  = _safesub(L,13, 16).strip()
        altLoc    = _safesub(L,17, 17).strip()
        resName   = _safesub(L,18, 20).strip().upper()
        chainID   = _safesub(L,22, 22).strip().upper()
        resSeq    = _safesub(L,23, 26).strip()
        iCode     = _safesub(L,27, 27).strip()
        X         = _safesub(L,31, 38).strip()
        Y         = _safesub(L,39, 46).strip()
        Z         = _safesub(L,47, 54).strip()
        elem      = _safesub(L,77, 78).strip()
        if not elem:
            nm = re.sub(r"[^A-Za-z]","", AtomName)
            elem = (nm[:1] or "").upper()
        try:
            atoms.append(dict(
                AtomSerNo=int(AtomSerNo or "0"),
                AtomName=AtomName, altLoc=altLoc, resName=resName, chainID=chainID,
                resSeq=int(resSeq or "0"), iCode=iCode,
                X=float(X or "nan"), Y=float(Y or "nan"), Z=float(Z or "nan"),
                element=elem
            ))
        except Exception:
            pass
    return atoms

def _is_hydrogen(atomname: str, element: str) -> bool:
    an = atomname.upper()
    if element.upper()=="H": return True
    if an in {"H","D","E","G","HN"}: return True
    if ("H" in an) and ("N" not in an): return True
    return False

def _stable_unique(seq):
    seen=set(); out=[]
    for x in seq:
        if x not in seen:
            out.append(x); seen.add(x)
    return out

def _chain_rank(all_atoms, req_chain):
    chains = _stable_unique([a["chainID"] for a in all_atoms if a["chainID"]])
    try:
        return chains.index(req_chain.upper())+1
    except ValueError:
        return None

def _avg(coords):
    n = len(coords)
    if n==0: return (0.0,0.0,0.0)
    sx = sum(c[0] for c in coords); sy = sum(c[1] for c in coords); sz = sum(c[2] for c in coords)
    return (sx/n, sy/n, sz/n)

def _rename_res(resname: str) -> str:
    return RESNAME_REMAP.get(resname.upper(), resname.upper())

# ---------------- Side-chain centroid rules (index-based; MATLAB-like) ----------------
def _sidechain_centroid_by_rule(rescode: int, per_res_atoms: List[Dict], ca_xyz):
    xyz = [(a["X"],a["Y"],a["Z"]) for a in per_res_atoms]
    n = len(xyz)
    def ok(k): return 1 <= k <= n

    if rescode==1:        # GLY -> 2nd atom if exists
        if ok(2): return xyz[1]
    elif rescode==2:      # ALA -> #5
        if ok(5): return xyz[4]
    elif rescode==3:      # VAL -> mean #6,#7
        if ok(7): return ((xyz[5][0]+xyz[6][0])/2.0, (xyz[5][1]+xyz[6][1])/2.0, (xyz[5][2]+xyz[6][2])/2.0)
    elif rescode==4:      # ILE -> #8
        if ok(8): return xyz[7]
    elif rescode==5:      # LEU -> mean #7,#8
        if ok(8): return ((xyz[6][0]+xyz[7][0])/2.0, (xyz[6][1]+xyz[7][1])/2.0, (xyz[6][2]+xyz[7][2])/2.0)
    elif rescode==6:      # SER -> #6
        if ok(6): return xyz[5]
    elif rescode==7:      # THR -> #6
        if ok(6): return xyz[5]
    elif rescode==8:      # ASP -> mean #7,#8
        if ok(8): return ((xyz[6][0]+xyz[7][0])/2.0, (xyz[6][1]+xyz[7][1])/2.0, (xyz[6][2]+xyz[7][2])/2.0)
    elif rescode==9:      # ASN -> mean #7,#8
        if ok(8): return ((xyz[6][0]+xyz[7][0])/2.0, (xyz[6][1]+xyz[7][1])/2.0, (xyz[6][2]+xyz[7][2])/2.0)
    elif rescode==10:     # GLU -> mean #8,#9
        if ok(9): return ((xyz[7][0]+xyz[8][0])/2.0, (xyz[7][1]+xyz[8][1])/2.0, (xyz[7][2]+xyz[8][2])/2.0)
    elif rescode==11:     # GLN -> mean #8,#9
        if ok(9): return ((xyz[7][0]+xyz[8][0])/2.0, (xyz[7][1]+xyz[8][1])/2.0, (xyz[7][2]+xyz[8][2])/2.0)
    elif rescode==12:     # LYS -> #9
        if ok(9): return xyz[8]
    elif rescode==13:     # ARG -> mean #8,#10,#11
        if ok(11):
            xs = xyz[7][0]+xyz[9][0]+xyz[10][0]
            ys = xyz[7][1]+xyz[9][1]+xyz[10][1]
            zs = xyz[7][2]+xyz[9][2]+xyz[10][2]
            return (xs/3.0, ys/3.0, zs/3.0)
    elif rescode==14:     # CYS -> #6
        if ok(6): return xyz[5]
    elif rescode==15:     # MET -> #7
        if ok(7): return xyz[6]
    elif rescode==16:     # PHE -> mean #6..#11
        if ok(11): return _avg(xyz[5:11])
    elif rescode==17:     # TYR -> mean #6..#12
        if ok(12): return _avg(xyz[5:12])
    elif rescode==18:     # TRP -> mean #7..#14
        if ok(14): return _avg(xyz[6:14])
    elif rescode==19:     # HIS -> mean #6..#10
        if ok(10): return _avg(xyz[5:10])
    elif rescode==20:     # PRO -> mean #5..#7
        if ok(7): return _avg(xyz[4:7])
    return ca_xyz

# ---------------- Core converter ----------------
def readpdb_to_cor(pdb_path: str, chain: str, out_basename: Optional[str]=None, strict_matlab_mode: bool=True) -> str:
    """
    Convert PDB to .cor.

    chain string semantics:
      ""           -> all chains
      "A"          -> only chain A (by rank)
      "A,B"        -> chains A and B (by rank)
      "A,B,C" ...  -> multiple
    If PDB has NO chain IDs, all residues are kept regardless of chain.
    """
    chain_raw = (chain or "").strip().upper()

    with open(pdb_path, "r") as f:
        raw = f.read().splitlines()

    lines = _trim_to_first_model(raw)
    atoms_all = _parse_pdb_atoms(lines)
    if not atoms_all:
        raise RuntimeError("No atoms parsed from PDB.")

    # Detect if PDB has any chain IDs at all
    has_any_chain = any(a["chainID"] for a in atoms_all)

    # Determine allowed chain ranks (or None for "all")
    allowed_ranks = None
    if has_any_chain and chain_raw:
        tokens = [t.strip() for t in re.split(r"[,\s]+", chain_raw) if t.strip()]
        ranks: List[int] = []
        for tk in tokens:
            rk = _chain_rank(atoms_all, tk)
            if rk is not None and rk not in ranks:
                ranks.append(rk)
        if ranks:
            allowed_ranks = set(ranks)
        else:
            # Requested chain(s) not found → fall back to ALL chains
            allowed_ranks = None
    else:
        # Either no chain IDs in PDB or empty chain string: use ALL residues
        allowed_ranks = None

    # Clean (drop hydrogens; simple altLoc filter; rename residues)
    cleaned = []
    i = 0
    while i < len(atoms_all):
        a = atoms_all[i]
        if _is_hydrogen(a["AtomName"], a["element"]):
            i += 1; continue
        if i+1 < len(atoms_all):
            b = atoms_all[i+1]
            if a["altLoc"] and b["altLoc"] and a["altLoc"] != b["altLoc"]:
                i += 1; continue
        cleaned.append(dict(
            AtomSerNo = a["AtomSerNo"],
            AtomName  = a["AtomName"].upper(),
            resName   = _rename_res(a["resName"]),
            chainID   = a["chainID"].upper(),
            resSeq    = a["resSeq"],
            X=a["X"], Y=a["Y"], Z=a["Z"],
            element=a["element"].upper()
        ))
        i += 1

    # Group by (chain,resSeq)
    res_blocks = defaultdict(list)
    for at in cleaned:
        res_blocks[(at["chainID"], at["resSeq"])] .append(at)

    out_rows = []
    for (chainID, resSeq) in sorted(res_blocks.keys(), key=lambda t: (t[0], t[1])):
        perres = res_blocks[(chainID, resSeq)]
        resname = _rename_res(perres[0]["resName"])
        rescode = AA_SET.get(resname, 0)
        atom_count = len(perres)
        CA = next((a for a in perres if a["AtomName"]=="CA"), None)
        if rescode==0 or CA is None:
            continue
        CAxyz = (CA["X"], CA["Y"], CA["Z"])
        SCxyz = _sidechain_centroid_by_rule(rescode, perres, CAxyz)
        mass  = RES_MASS.get(rescode, 0.0)
        chain_rank_val = _chain_rank(atoms_all, chainID) or 1
        out_rows.append((
            int(resSeq),
            int(rescode),
            float(CAxyz[0]), float(CAxyz[1]), float(CAxyz[2]),
            float(SCxyz[0]), float(SCxyz[1]), float(SCxyz[2]),
            int(atom_count),
            int(chain_rank_val),
            float(mass)
        ))

    # Filter by allowed_ranks if any
    if allowed_ranks is not None:
        out_rows = [r for r in out_rows if r[9] in allowed_ranks]

    # copy CA to SC if SC is zero (MATLAB fallback)
    fixed = []
    for r in out_rows:
        r = list(r)
        if r[5]==0.0 and r[6]==0.0 and r[7]==0.0:
            r[5], r[6], r[7] = r[2], r[3], r[4]
        fixed.append(tuple(r))
    out_rows = fixed

    out_base = out_basename or os.path.basename(pdb_path)
    out_path = os.path.join(os.path.dirname(pdb_path), f"{out_base}.cor")
    with open(out_path, "w") as f:
        for r in out_rows:
            f.write(f"{r[0]:5.0f} {r[1]:4.0f} {r[2]:8.3f} {r[3]:8.3f} {r[4]:8.3f} "
                    f"{r[5]:8.3f} {r[6]:8.3f} {r[7]:8.3f} {r[8]:4.0f} {r[9]:2.0f} {r[10]:4.1f}\n")
    print(f"✔ Wrote: {out_path}")
    return out_path

# ---------------- mcpath_input.txt parsing & runner ----------------
def _find_default_input_candidates(start_dir: Optional[str]) -> List[str]:
    cands = []
    if start_dir:
        cands.append(os.path.join(start_dir, "mcpath_input.txt"))
    cands += [
        os.path.join(os.getcwd(), "mcpath_input.txt"),
        "/content/mcpath_input.txt",
    ]
    seen=set(); out=[]
    for p in cands:
        if p not in seen:
            out.append(p); seen.add(p)
    return out

def parse_mcpath_input(input_path: str) -> Dict[str,str]:
    with open(input_path, "r") as f:
        lines = [ln.strip() for ln in f.readlines()]
    if not lines or not lines[0].isdigit():
        raise ValueError("mcpath_input.txt malformed: first line must be mode 1/2/3")
    mode = int(lines[0])
    if mode == 1 and len(lines) >= 3:
        return {"mode": "functional", "pdb_name": lines[1], "chain": lines[2]}
    if mode == 2 and len(lines) >= 8:
        return {"mode": "paths_init_len", "pdb_name": lines[1], "chain": lines[2]}
    if mode == 3 and len(lines) >= 8:
        return {"mode": "paths_init_final", "pdb_name": lines[1], "chain": lines[2]}
    raise ValueError("mcpath_input.txt malformed: insufficient lines for selected mode.")

def run_from_input(input_path: str, base_dir: Optional[str]=None) -> str:
    info = parse_mcpath_input(input_path)
    pdb_name = info["pdb_name"]
    chain    = info["chain"].strip().upper()
    root = base_dir or os.path.dirname(os.path.abspath(input_path))
    pdb_path = pdb_name if os.path.isabs(pdb_name) else os.path.join(root, pdb_name)
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB not found: {pdb_path}")
    print(f"▶ Using input file: {input_path}")
    print(f"   Mode={info['mode']}, PDB={pdb_path}, Chain='{chain or 'ALL'}'")
    return readpdb_to_cor(pdb_path, chain, out_basename=os.path.basename(pdb_name))

# ---------------- Unified entrypoint (used by UI) ----------------
def run(pdb_path: Optional[str]=None,
        chain: Optional[str]=None,
        out_basename: Optional[str]=None,
        strict_matlab_mode: bool=True,
        input_path: Optional[str]=None) -> str:
    """
    Preferred behavior:
      - If input_path is given (or auto-found), use mcpath_input.txt.
      - Else, require pdb_path & chain and convert directly.

    chain may be "" or comma-separated as in readpdb_to_cor.
    """
    if input_path:
        return run_from_input(input_path)

    start_dir = os.path.dirname(os.path.abspath(pdb_path)) if pdb_path else None
    for cand in _find_default_input_candidates(start_dir):
        if os.path.isfile(cand):
            return run_from_input(cand, base_dir=start_dir)

    if not pdb_path:
        raise ValueError("No mcpath_input.txt found and pdb_path not provided.")
    pdb_path = os.path.abspath(pdb_path)
    chain_str = (chain or "").strip().upper()
    print(f"▶ No input file found; converting directly: {pdb_path}, chain='{chain_str or 'ALL'}'")
    return readpdb_to_cor(pdb_path, chain_str, out_basename=out_basename or os.path.basename(pdb_path))
