# readpdb_strict.py
# STRICT MATLAB MODE reader that mirrors your readpdb.m behavior.
import os
import numpy as np

AA_CODE = {
    "GLY": 1, "ALA": 2, "VAL": 3, "ILE": 4, "LEU": 5,
    "SER": 6, "THR": 7, "ASP": 8, "ASN": 9, "GLU":10,
    "GLN":11, "LYS":12, "ARG":13, "CYS":14, "MET":15,
    "PHE":16, "TYR":17, "TRP":18, "HIS":19, "PRO":20
}
AA_CONST = {
    1:  75.0,  2:  89.1,  3: 117.2,  4: 131.2,  5: 131.2,
    6: 105.1,  7: 119.1,  8: 133.1,  9: 132.1, 10: 147.1,
   11: 146.0, 12: 146.2, 13: 174.2, 14: 121.2, 15: 149.2,
   16: 165.2, 17: 181.2, 18: 204.2, 19: 155.2, 20: 115.1,
}
# thresholds (must be strictly greater than) copied from MATLAB branches
AA_MIN_ATOMS = {
    1:  1,  2:  4,  3:  6,  4:  7,  5:  7,
    6:  5,  7:  5,  8:  7,  9:  7, 10:  8,
   11:  8, 12:  8, 13: 10, 14:  5, 15:  6,
   16: 10, 17: 11, 18: 13, 19:  9, 20:  6,
}
# side-chain index/average recipes (1-based indices on "heavy" list)
AA_SC_RULES = {
     1: ("index", [2]),                         # GLY
     2: ("index", [5]),                         # ALA
     3: ("mean",  [6,7]),                       # VAL
     4: ("index", [8]),                         # ILE
     5: ("mean",  [7,8]),                       # LEU
     6: ("index", [6]),                         # SER
     7: ("index", [6]),                         # THR
     8: ("mean",  [7,8]),                       # ASP
     9: ("mean",  [7,8]),                       # ASN
    10: ("mean",  [8,9]),                       # GLU
    11: ("mean",  [8,9]),                       # GLN
    12: ("index", [9]),                         # LYS
    13: ("mean",  [8,10,11]),                   # ARG
    14: ("index", [6]),                         # CYS
    15: ("index", [7]),                         # MET
    16: ("mean",  [6,7,8,9,10,11]),             # PHE
    17: ("mean",  [6,7,8,9,10,11,12]),          # TYR
    18: ("mean",  [7,8,9,10,11,12,13,14]),      # TRP
    19: ("mean",  [6,7,8,9,10]),                # HIS
    20: ("mean",  [5,6,7]),                     # PRO
}
ALTNAME_NORMALIZE = {
    "HIP":"HIS", "HIE":"HIS", "HID":"HIS", "HSD":"HIS", "HSE":"HIS", "HSP":"HIS",
    "LYN":"LYS", "GLH":"GLU", "CYM":"CYS", "CYX":"CYS", "ASH":"ASN", "TYM":"TYR"
}

def _safe_sub(s, a, b):
    a = max(1, a); b = min(len(s), b)
    return s[a-1:b] if a <= b else ""

def _trim_to_first_model(pdb_lines):
    has_model = any(l.startswith("MODEL") for l in pdb_lines)
    if not has_model: return pdb_lines
    new_lines = []; in_model = False
    for L in pdb_lines:
        if L.startswith("MODEL"):
            if in_model: break
            in_model = True; continue
        if L.startswith("ENDMDL") and in_model: break
        if not in_model:
            if L.startswith("ATOM") or L.startswith("HETATM") or L.startswith("TER"):
                new_lines.append(L)
        else:
            new_lines.append(L)
    return new_lines

def _parse_atoms_from_path(pdb_path):
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.read().splitlines()
    lines = _trim_to_first_model(lines)
    atoms = []
    for L in lines:
        if not (L.startswith("ATOM") or L.startswith("HETATM")): continue
        if len(L) < 54: continue
        try:
            atoms.append({
                "AtomSerNo": int(_safe_sub(L,7,11).strip() or "0"),
                "AtomName" : _safe_sub(L,13,16).strip(),
                "altLoc"   : _safe_sub(L,17,17).strip(),
                "resName"  : _safe_sub(L,18,20).strip(),
                "chainID"  : _safe_sub(L,22,22).strip(),
                "resSeq"   : int(_safe_sub(L,23,26).strip() or "0"),
                "iCode"    : _safe_sub(L,27,27).strip(),
                "X"        : float(_safe_sub(L,31,38)),
                "Y"        : float(_safe_sub(L,39,46)),
                "Z"        : float(_safe_sub(L,47,54)),
                "element"  : _safe_sub(L,77,78).strip(),
                "raw"      : L,
            })
        except:
            continue
    return atoms

def _normalize_resname(nm):
    u = (nm or "").upper()
    return ALTNAME_NORMALIZE.get(u, u)

def _matlab_h_skip(atom_name: str, element: str) -> bool:
    nm = (atom_name or "").upper()
    el = (element or "").upper()
    if el == "H" or nm in {"H","D","E","G","HN"}:
        return True
    if ("H" in nm) and ("N" not in nm):
        return True
    return False

def _chain_rank_map_stream(atoms):
    ranks = {}; r = 0; last = None
    for a in atoms:
        ch = (a["chainID"] or "").upper()
        if ch != last:
            if ch not in ranks:
                r += 1; ranks[ch] = r
            last = ch
    return ranks

def _build_residues_stream(atoms):
    # normalize names + chain case
    for a in atoms:
        a["resName"] = _normalize_resname(a["resName"])
        a["chainID"] = (a["chainID"] or "").upper()

    # group by (chain,resSeq,iCode) in file order
    res_list = []
    cur_key = None
    cur_atoms = []
    ordered = []
    for a in atoms:
        key = (a["chainID"], a["resSeq"], a["iCode"])
        if key != cur_key:
            if cur_atoms:
                ordered.append((cur_key, cur_atoms))
            cur_key = key
            cur_atoms = [a]
        else:
            cur_atoms.append(a)
    if cur_atoms:
        ordered.append((cur_key, cur_atoms))

    # apply neighbor altLoc rule and hydrogen skip
    for key, atoms_in_res in ordered:
        keep = [True] * len(atoms_in_res)
        for i in range(len(atoms_in_res)-1):
            a = atoms_in_res[i]; b = atoms_in_res[i+1]
            if (a.get("altLoc") or "") and (b.get("altLoc") or "") and (a["altLoc"] != b["altLoc"]):
                keep[i] = False
        heavy = []
        for i, a in enumerate(atoms_in_res):
            if not keep[i]:
                continue
            if _matlab_h_skip(a.get("AtomName",""), a.get("element","")):
                continue
            heavy.append(a)
        ch, resi, ic = key
        res_list.append({
            "key": key, "chain": ch, "resSeq": resi, "iCode": ic,
            "resName": atoms_in_res[0]["resName"] if atoms_in_res else "",
            "atoms": atoms_in_res, "heavy": heavy
        })
    return res_list

def _pick_sc_xyz_matlab_like(aa_code, heavy):
    min_need = AA_MIN_ATOMS.get(aa_code, 10**9)
    if len(heavy) <= min_need:
        return None
    mode, idxs = AA_SC_RULES[aa_code]
    coords = []
    for idx1 in idxs:   # 1-based
        j = idx1 - 1
        if 0 <= j < len(heavy):
            a = heavy[j]
            coords.append((a["X"], a["Y"], a["Z"]))
        else:
            return None  # out-of-range → fallback to CA
    if mode == "index":
        return coords[0]
    arr = np.array(coords, dtype=float)
    c = arr.mean(axis=0)
    return float(c[0]), float(c[1]), float(c[2])

def readpdb_py(pdb_path: str, chainID: str | None):
    """
    STRICT MATLAB MODE.
    Returns: (cor_path, n_rows, n_problems)
    Also writes optional 'problem' file like MATLAB.
    """
    chainID = (chainID or "").strip()[:1].upper()
    atoms_all = _parse_atoms_from_path(pdb_path)
    ranks_full = _chain_rank_map_stream(atoms_all)

    # chain filter after rank calculation
    if chainID:
        atoms_use = [a for a in atoms_all if (a["chainID"] or "").upper() == chainID]
        if not atoms_use:
            print(f'WARNING: requested chain "{chainID}" not found. No chain filter applied.')
            atoms_use = atoms_all
    else:
        atoms_use = atoms_all

    residues = _build_residues_stream(atoms_use)

    result_rows = []
    problems = []
    for res in residues:
        rn = res["resName"]
        if rn not in AA_CODE:
            continue
        heavy = res["heavy"]
        # CA from heavy first, then fallback to all atoms
        ca = next((a for a in heavy if (a.get("AtomName","").upper() == "CA")), None)
        if ca is None:
            ca = next((a for a in res["atoms"] if (a.get("AtomName","").upper() == "CA")), None)
        if ca is None:
            continue

        aa_code = AA_CODE[rn]
        chain_rank = float(ranks_full.get(res["chain"], 0))
        res_atom_count = float(len(heavy))

        sc = _pick_sc_xyz_matlab_like(aa_code, heavy)
        if sc is None:
            scx, scy, scz = float(ca["X"]), float(ca["Y"]), float(ca["Z"])
        else:
            scx, scy, scz = sc

        const = AA_CONST[aa_code]
        result_rows.append([
            float(res["resSeq"]), float(aa_code),
            float(ca["X"]), float(ca["Y"]), float(ca["Z"]),
            float(scx), float(scy), float(scz),
            res_atom_count, chain_rank, float(const)
        ])

        if res_atom_count <= AA_MIN_ATOMS.get(aa_code, 0):
            problems.append((res["resSeq"], int(chain_rank)))

    cor_path = f"{pdb_path}.cor"
    with open(cor_path, "w", encoding="utf-8") as f:
        for r in result_rows:
            f.write(f"{r[0]:5.0f} {r[1]:4.0f} {r[2]:8.3f} {r[3]:8.3f} {r[4]:8.3f} "
                    f"{r[5]:8.3f} {r[6]:8.3f} {r[7]:8.3f} {r[8]:4.0f} {r[9]:2.0f} {r[10]:4.1f}\n")

    if problems:
        with open("problem", "w", encoding="utf-8") as f:
            f.write("There are problems about number of atoms on following residues\n")
            for resseq, crank in problems:
                f.write(f"{int(resseq)} of {int(crank)} chain\n")

    return cor_path, len(result_rows), len(problems)
