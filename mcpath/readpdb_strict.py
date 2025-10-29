"""
Strict-MATLAB-mode PDB → .cor converter (pure Python).

- Parses ATOM/HETATM lines from a PDB file.
- Replicates your MATLAB readpdb.m filtering for chain, MODEL trimming,
  hydrogen skipping, residue renames (HIP/HIE/HID→HIS, etc.), and
  side-chain centroid logic.
- Writes "<input>.pdb.cor" with the **same 11-column layout** your MATLAB code writes.

Entry point expected by the UI:
    run(pdb_path: str, chain: str, out_basename: str|None=None, strict_matlab_mode: bool=True)
"""

from __future__ import annotations
import os, re
from collections import defaultdict

AA_SET = {
    "GLY":1, "ALA":2, "VAL":3, "ILE":4, "LEU":5, "SER":6, "THR":7, "ASP":8,
    "ASN":9, "GLU":10, "GLN":11, "LYS":12, "ARG":13, "CYS":14, "MET":15,
    "PHE":16, "TYR":17, "TRP":18, "HIS":19, "PRO":20
}
RESNAME_REMAP = {
    "HIP":"HIS","HIE":"HIS","HID":"HIS","HSD":"HIS","HSE":"HIS","HSP":"HIS",
    "LYN":"LYS","GLH":"GLU","CYM":"CYS","CYX":"CYS","ASH":"ASN","TYM":"TYR"
}
RES_MASS = { # column 11 in your output (same constants you use)
    1:75,   2:89.1, 3:117.2,4:131.2,5:131.2,6:105.1,7:119.1,8:133.1,9:132.1,10:147.1,
    11:146, 12:146.2,13:174.2,14:121.2,15:149.2,16:165.2,17:181.2,18:204.2,19:155.2,20:115.1
}

def _safesub(s: str, a: int, b: int) -> str:
    a = max(1, a); b = min(len(s), b)
    return "" if a > b else s[a-1:b]

def _trim_to_first_model(pdb_text: list[str]) -> list[str]:
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
            # keep ATOM/HETATM/TER that appear before MODEL (rare)
            if L.startswith(("ATOM  ","HETATM","TER")):
                out.append(L)
        else:
            out.append(L)
    return out

def _parse_pdb_atoms(lines: list[str]) -> list[dict]:
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
            # guess element from AtomName like MATLAB fallback
            nm = re.sub(r"[^A-Za-z]","", AtomName)
            elem = (nm[:1] or "").upper()
        try:
            d = dict(
                AtomSerNo=int(AtomSerNo or "0"),
                AtomName=AtomName, altLoc=altLoc, resName=resName, chainID=chainID,
                resSeq=int(resSeq or "0"), iCode=iCode,
                X=float(X or "nan"), Y=float(Y or "nan"), Z=float(Z or "nan"),
                element=elem
            )
            atoms.append(d)
        except Exception:
            # ignore malformed rows
            pass
    return atoms

def _is_hydrogen(atomname: str, element: str) -> bool:
    # MATLAB criteria: element=='H' or atomname=='H' or name in {'D','E','G','HN'} or contains 'H' but not 'N'
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

def _sidechain_centroid_by_rule(rescode: int, per_res_atoms: list[dict], ca_xyz: tuple[float,float,float]):
    """
    Re-creates your residue-specific rules:
    - Use defined atom indices to average; if insufficient atom count in residue,
      fall back to CA coordinates.
    """
    # Build list in residue order
    xyz = [(a["X"],a["Y"],a["Z"]) for a in per_res_atoms]
    n = len(xyz)

    def ok(k): return k<=n and k>=1

    if rescode==1:        # GLY -> take the 2nd atom if exists, else CA
        if ok(2): return xyz[1]
    elif rescode==2:      # ALA -> atom #5 else CA
        if ok(5): return xyz[4]
    elif rescode==3:      # VAL -> mean of atoms #6 and #7
        if ok(7): return ((xyz[5][0]+xyz[6][0])/2.0, (xyz[5][1]+xyz[6][1])/2.0, (xyz[5][2]+xyz[6][2])/2.0)
    elif rescode==4:      # ILE -> atom #8
        if ok(8): return xyz[7]
    elif rescode==5:      # LEU -> mean of #7 and #8
        if ok(8): return ((xyz[6][0]+xyz[7][0])/2.0, (xyz[6][1]+xyz[7][1])/2.0, (xyz[6][2]+xyz[7][2])/2.0)
    elif rescode==6:      # SER -> atom #6
        if ok(6): return xyz[5]
    elif rescode==7:      # THR -> atom #6
        if ok(6): return xyz[5]
    elif rescode==8:      # ASP -> mean of #7 and #8
        if ok(8): return ((xyz[6][0]+xyz[7][0])/2.0, (xyz[6][1]+xyz[7][1])/2.0, (xyz[6][2]+xyz[7][2])/2.0)
    elif rescode==9:      # ASN -> mean of #7 and #8
        if ok(8): return ((xyz[6][0]+xyz[7][0])/2.0, (xyz[6][1]+xyz[7][1])/2.0, (xyz[6][2]+xyz[7][2])/2.0)
    elif rescode==10:     # GLU -> mean of #8 and #9
        if ok(9): return ((xyz[7][0]+xyz[8][0])/2.0, (xyz[7][1]+xyz[8][1])/2.0, (xyz[7][2]+xyz[8][2])/2.0)
    elif rescode==11:     # GLN -> mean of #8 and #9
        if ok(9): return ((xyz[7][0]+xyz[8][0])/2.0, (xyz[7][1]+xyz[8][1])/2.0, (xyz[7][2]+xyz[8][2])/2.0)
    elif rescode==12:     # LYS -> atom #9
        if ok(9): return xyz[8]
    elif rescode==13:     # ARG -> mean of #8,#10,#11
        if ok(11): 
            xs = xyz[7][0]+xyz[9][0]+xyz[10][0]
            ys = xyz[7][1]+xyz[9][1]+xyz[10][1]
            zs = xyz[7][2]+xyz[9][2]+xyz[10][2]
            return (xs/3.0, ys/3.0, zs/3.0)
    elif rescode==14:     # CYS -> atom #6
        if ok(6): return xyz[5]
    elif rescode==15:     # MET -> atom #7
        if ok(7): return xyz[6]
    elif rescode==16:     # PHE -> mean of #6..#11
        if ok(11): return _avg(xyz[5:11])
    elif rescode==17:     # TYR -> mean of #6..#12
        if ok(12): return _avg(xyz[5:12])
    elif rescode==18:     # TRP -> mean of #7..#14
        if ok(14): return _avg(xyz[6:14])
    elif rescode==19:     # HIS -> mean of #6..#10
        if ok(10): return _avg(xyz[5:10])
    elif rescode==20:     # PRO -> mean of #5..#7
        if ok(7): return _avg(xyz[4:7])

    return ca_xyz  # fallback: CA

def _rename_res(resname: str) -> str:
    return RESNAME_REMAP.get(resname.upper(), resname.upper())

def readpdb_to_cor(pdb_path: str, chain: str, out_basename: str|None=None, strict_matlab_mode: bool=True):
    """
    Main converter. Writes "<pdb_path>.cor".
    Returns output file path.
    """
    chain = (chain or "").strip().upper()
    if len(chain)!=1:
        raise ValueError("chain must be a single letter/char, e.g., 'A'.")

    with open(pdb_path, "r") as f:
        raw = f.read().splitlines()

    # Trim to first MODEL (MATLAB behavior)
    lines = _trim_to_first_model(raw)

    # Parse atoms
    atoms = _parse_pdb_atoms(lines)
    if not atoms:
        raise RuntimeError("No atoms parsed from PDB.")

    # Filter by requested chain rank (like MATLAB using rank in result col10)
    rank = _chain_rank(atoms, chain)  # 1-based within unique chains order
    if rank is None:
        # warn but continue (MATLAB emits a warning and keeps all)
        pass

    # Build per-residue lists in order (skipping hydrogens & conflicting altLoc pairs)
    cleaned = []
    i = 0; n = len(atoms); mult_count = 0
    while i < n:
        a = atoms[i]
        # skip H
        if _is_hydrogen(a["AtomName"], a["element"]):
            i += 1
            continue
        # skip altLoc pair mismatch (if consecutive same atom w/ different altLoc)
        if i+1 < n:
            a2 = atoms[i+1]
            if a["altLoc"] and a2["altLoc"] and (a["altLoc"] != a2["altLoc"]):
                i += 1
                continue

        # normalize residue names
        resname = _rename_res(a["resName"])

        cleaned.append(dict(
            AtomSerNo=a["AtomSerNo"],
            AtomName=a["AtomName"].upper(),
            resName=resname,
            chainID=a["chainID"].upper(),
            resSeq=a["resSeq"],
            X=a["X"], Y=a["Y"], Z=a["Z"],
            element=a["element"].upper()
        ))
        i += 1

    # walk and build output rows per residue (CA rows for standard AAs)
    out_rows = []  # rows for writing later
    current_chain_numeric = 1
    chain_end_marks = []
    # compute chain numeric indices by boundaries (MATLAB newchainid)
    for j in range(1, len(cleaned)):
        if cleaned[j-1]["chainID"] != cleaned[j]["chainID"]:
            chain_end_marks.append(j)
    chain_end_marks.append(len(cleaned))  # last end

    # map residue → atoms aggregated
    res_blocks = defaultdict(list)
    for a in cleaned:
        res_blocks[(a["chainID"], a["resSeq"])].append(a)

    # produce per-residue outputs
    # result columns (11): resID, resCode(1..20), CAx, CAy, CAz, SCx, SCy, SCz, atom_count_in_res, chain_rank, mass
    # We also keep track to compute "atom_count_in_res" like MATLAB
    for (chainID, resSeq) in sorted(res_blocks.keys(), key=lambda t: (t[0], t[1])):
        perres = res_blocks[(chainID, resSeq)]
        # residue type code
        resname = _rename_res(perres[0]["resName"])
        rescode = AA_SET.get(resname, 0)

        # find CA in perres and atom count
        atom_count = len(perres)
        CA = next((a for a in perres if a["AtomName"]=="CA"), None)
        if rescode==0 or CA is None:
            continue  # keep only standard AA with CA

        CAxyz = (CA["X"], CA["Y"], CA["Z"])
        SCxyz = _sidechain_centroid_by_rule(rescode, perres, CAxyz)
        mass  = RES_MASS.get(rescode, 0.0)

        # Determine chain numeric like MATLAB "newchainid": assign incremental numbers per boundary
        # Use stable order of unique chains; MATLAB uses boundary-based numeric but we can emulate by rank when available
        # For compatibility with your later filters, we store both:
        #   col10 = rank (1-based) among unique chains (if available)
        #   col9  = per-res atom count
        chain_rank_val = _chain_rank(atoms, chainID) or 1

        out_rows.append((
            int(resSeq),                 # 1 resID
            int(rescode),                # 2 resCode 1..20
            float(CAxyz[0]),             # 3
            float(CAxyz[1]),             # 4
            float(CAxyz[2]),             # 5
            float(SCxyz[0]),             # 6
            float(SCxyz[1]),             # 7
            float(SCxyz[2]),             # 8
            int(atom_count),             # 9 atom count in residue
            int(chain_rank_val),         # 10 chain rank (for filter)
            float(mass)                  # 11 residue mass constant
        ))

    # keep only requested chain (same as MATLAB branch using 'rank')
    if rank is not None:
        out_rows = [r for r in out_rows if r[9] == rank]

    # side-chain missing fallback (zeros → copy CA) — already avoided above,
    # but we keep exact MATLAB behavior if any zeros slipped in:
    fixed_rows = []
    for r in out_rows:
        r = list(r)
        if (abs(r[5])==0 and abs(r[6])==0 and abs(r[7])==0):
            r[5],r[6],r[7] = r[2],r[3],r[4]
        fixed_rows.append(tuple(r))
    out_rows = fixed_rows

    # write output
    out_base = out_basename or os.path.basename(pdb_path)
    out_path = os.path.join(os.path.dirname(pdb_path), f"{out_base}.cor")
    with open(out_path, "w") as f:
        for r in out_rows:
            f.write(f"{r[0]:5.0f} {r[1]:4.0f} {r[2]:8.3f} {r[3]:8.3f} {r[4]:8.3f} "
                    f"{r[5]:8.3f} {r[6]:8.3f} {r[7]:8.3f} {r[8]:4.0f} {r[9]:2.0f} {r[10]:4.1f}\n")
    return out_path

# ---------- UI entrypoint wrapper (DON'T CHANGE NAME) ----------
def _call_first_existing_fn(fn_names, *args, **kwargs):
    g = globals()
    for name in fn_names:
        fn = g.get(name)
        if callable(fn):
            try:
                return fn(*args, **kwargs)
            except TypeError:
                return fn(*args)
    raise RuntimeError("readpdb_strict: no suitable callable found.")

def run(pdb_path: str, chain: str, out_basename: str | None = None, strict_matlab_mode: bool = True):
    """
    Standard callable expected by the UI. Keeps behavior identical to MATLAB.
    """
    pdb_path = os.path.abspath(pdb_path)
    if out_basename is None:
        out_basename = os.path.basename(pdb_path)
    return _call_first_existing_fn(
        ("readpdb_to_cor",), pdb_path, chain, out_basename=out_basename, strict_matlab_mode=strict_matlab_mode
    )
