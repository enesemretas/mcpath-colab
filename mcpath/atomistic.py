# -*- coding: utf-8 -*-
"""
PROGRAM FOR RESIDUAL INTERACTION PROBABILITY CALCULATION
USING ATOMISTIC VAN DER WAALS POTENTIALS (Python / Colab version)

Inputs (same as MATLAB version):
- 'vdw_cns.param' : parameter file (for atom-type -> epsilon,sigma)
- 'pdb_cns.top'   : topology file  (for (res,atom) -> parameter-type mapping)
- '<ProteinName>.pdb' : the structure file
- '<ProteinName>.cor' : coordinates/atom-block info (from your readpdb workflow)
Output:
- '<ProteinName>_atomistic.out' : normalized probability matrix (tab-separated)

Notes:
- kT is set to 1.0 to reproduce your MATLAB behavior exactly
- Hydrogens & altLoc handling are matched to your MATLAB rules
- For r < rmin, energy uses the constant value with rmin in the denominator (same as your code)
- For r > rcut, energy = 0
- Diagonal energies forced to 0; diagonal of normalized matrix forced to 0
"""

import os
import re
import math
import numpy as np

# ---- Install Biopython if missing (for clean PDB parsing) ----
try:
    import Bio  # noqa
except ImportError:
    import sys, subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython"])
from Bio.PDB import PDBParser


# -------------------------- Helpers --------------------------

def _trim(s):
    """Trim and return a plain Python str even if the input is bytes/NumPy scalar."""
    if isinstance(s, bytes):
        s = s.decode("utf-8", "ignore")
    return str(s).strip()


def parse_vdw_params(param_path):
    """
    Parse 'vdw_cns.param' lines that start with ' NONBonded' (1-based in MATLAB).
    We expect fixed columns similar to the MATLAB slicing:
      partype  = cols 13-16
      pareps   = cols 21-26
      parsigma = cols 30-35
    Returns: dict[partype] = (eps, sigma)
    """
    partype_to_epssigma = {}
    with open(param_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            # keep exact behavior: line[1:10] in MATLAB == Python line[1:10]
            if len(line) >= 10 and line[1:10] == "NONBonded":
                # Use safe slicing (pad if short)
                seg = line.rstrip("\n")
                def field(seg, i, j):
                    return _trim(seg[i:j]) if len(seg) >= j else _trim(seg[i:])

                partype  = field(seg, 12, 16)  # MATLAB 13:16 -> Python 12:16
                pareps_s = field(seg, 20, 26)  # MATLAB 21:26 -> Python 20:26
                parsig_s = field(seg, 29, 35)  # MATLAB 30:35 -> Python 29:35

                if partype:
                    try:
                        eps   = float(pareps_s)
                        sigma = float(parsig_s)
                        partype_to_epssigma[partype] = (eps, sigma)
                    except ValueError:
                        # skip malformed numeric rows silently,
                        # mirroring tolerance of the MATLAB script
                        pass
    return partype_to_epssigma


def parse_topology(top_path):
    """
    Parse 'pdb_cns.top' with MATLAB columns:
      parrestype = cols 1:3
      paratomtype= cols 11:13
      partypecor = cols 21:24  # (map to vdw param 'partype')
      charge     = cols 34:38  # (unused here, but parsed anyway)
    Returns: dict[(res, atom)] -> mapped_type
    """
    mapping = {}
    with open(top_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            seg = line.rstrip("\n")
            def field(seg, i, j):
                return _trim(seg[i:j]) if len(seg) >= j else _trim(seg[i:])

            res  = field(seg, 0, 3)    # 1:3 -> 0:3
            atm  = field(seg, 10, 13)  # 11:13 -> 10:13
            corr = field(seg, 20, 24)  # 21:24 -> 20:24
            # charge = field(seg, 33, 38)  # 34:38 -> 33:38 (unused)

            if res and atm and corr:
                mapping[(res.upper(), atm.upper())] = corr
    return mapping


def load_pdb_atoms(pdb_path):
    """
    Load atoms from PDB (like your pdbreader).
    Filters:
      - Only atoms with iCode='' and altLoc in {'', ' ', 'A'} (choose first conformer)
      - Skip hydrogens via (element=='H') or special MATLAB-style name rules
    Returns dict with arrays for matched heavy atoms:
        {
          'serial': np.array(int),
          'atom':   list[str],
          'res':    list[str],
          'chain':  list[str],
          'resseq': np.array(int),
          'xyz':    (N,3) float,
          'element': list[str],
          'altloc': list[str],
          'icode':  list[str],
        }
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", pdb_path)

    serials, atom_names, resnames, chains, resseqs = [], [], [], [], []
    xs, ys, zs = [], [], []
    elements, altlocs, icodes = [], [], []

    # Simple first-occurrence for altLocs: prefer '' or 'A'
    seen_altloc_key = set()

    for model in structure:
        for chain in model:
            for residue in chain:
                # residue id: (' ', resseq, icode)
                het, resseq, icode = residue.id
                resname = _trim(residue.get_resname()).upper()

                for atom in residue:
                    atom_name = _trim(atom.get_name()).upper()
                    element = _trim(atom.element).upper() if atom.element else ""
                    altloc = _trim(getattr(atom, "altloc", "")).upper()
                    serial = int(atom.get_serial_number())

                    # Keep only main conformer (similar to MATLAB's strict altLoc logic)
                    # Build a key: same position (chain, resseq, atom_name)
                    key = (chain.id, resseq, icode, atom_name)
                    if key in seen_altloc_key:
                        continue
                    # accept blank or 'A' first
                    if altloc not in ("", " ", "A"):
                        # if we haven't seen any altloc for this key yet, skip these non-primary ones
                        # (mimics dropping altLoc mismatches)
                        continue
                    seen_altloc_key.add(key)

                    # Hydrogen / name-based hydrogen rules (match MATLAB)
                    # - If element == 'H' -> skip
                    # - If atom name starts with 'H' and does NOT contain 'N' anywhere -> skip
                    if element == "H":
                        continue
                    if atom_name.startswith("H") and ("N" not in atom_name):
                        continue
                    # Also skip specific atom codes 'D','E','G','HN' (mirrors your checks)
                    if atom_name in ("D", "E", "G", "HN"):
                        continue

                    # Normalize non-standard residue names to canonical
                    resname_norm = resname
                    if resname_norm in ("HIP","HIE","HID","HSD","HSE","HSP"):
                        resname_norm = "HIS"
                    elif resname_norm == "LYN":
                        resname_norm = "LYS"
                    elif resname_norm == "GLH":
                        resname_norm = "GLU"
                    elif resname_norm in ("CYM","CYX"):
                        resname_norm = "CYS"
                    elif resname_norm == "ASH":
                        resname_norm = "ASN"
                    elif resname_norm == "TYM":
                        resname_norm = "TYR"

                    # Accept atom
                    x, y, z = atom.coord
                    serials.append(serial)
                    atom_names.append(atom_name)
                    resnames.append(resname_norm)
                    chains.append(_trim(chain.id).upper())
                    resseqs.append(int(resseq))
                    xs.append(float(x)); ys.append(float(y)); zs.append(float(z))
                    elements.append(element)
                    altlocs.append(altloc)
                    icodes.append(_trim(icode))

    xyz = np.column_stack([xs, ys, zs]).astype(float)
    return {
        "serial": np.array(serials, dtype=int),
        "atom":   atom_names,
        "res":    resnames,
        "chain":  chains,
        "resseq": np.array(resseqs, dtype=int),
        "xyz":    xyz,
        "element": elements,
        "altloc": altlocs,
        "icode":  icodes,
    }


def map_eps_sigma(atoms, top_map, vdw_map):
    """
    Map each atom to (eps, sigma) via:
       (res, atom) --(pdb_cns.top)--> partype --> (vdw_cns.param) --> (eps,sigma)
    Also perform the MATLAB 'O' normalizations for atom names:
       if atom in {'OXT','O2','O1'} or element=='O' => atom = 'O'
    Return filtered atoms (only those with mapped params) and arrays eps, sigma.
    """
    res = atoms["res"]
    atm = atoms["atom"]
    elem = atoms["element"]

    # Normalize atom name for O variants
    norm_atom = []
    for a, el in zip(atm, elem):
        if a in ("OXT", "O2", "O1") or el == "O":
            norm_atom.append("O")
        else:
            norm_atom.append(a)

    eps_list, sig_list, keep_idx = [], [], []

    for i, (rname, aname) in enumerate(zip(res, norm_atom)):
        key = (rname.upper(), aname.upper())
        ptype = top_map.get(key, None)
        if ptype is None:
            continue
        es = vdw_map.get(ptype, None)
        if es is None:
            continue
        eps, sig = es
        # ignore zero params just like MATLAB collects "errors"
        if (eps == 0.0) or (sig == 0.0):
            continue
        keep_idx.append(i)
        eps_list.append(float(eps))
        sig_list.append(float(sig))

    if not keep_idx:
        raise RuntimeError("No atoms found with valid (eps, sigma) mapping. "
                           "Check that pdb_cns.top and vdw_cns.param match your PDB.")

    keep_idx = np.array(keep_idx, dtype=int)
    # Slice atoms down to only those successfully mapped
    filtered_atoms = {
        k: (val[keep_idx] if isinstance(val, np.ndarray) else [val[j] for j in keep_idx])
        if k in ("serial", "resseq", "xyz") else
        ([atoms[k][j] for j in keep_idx] if isinstance(atoms[k], list) else atoms[k])
        for k in atoms.keys()
    }
    # Fix xyz selection (np.ndarray)
    filtered_atoms["xyz"] = atoms["xyz"][keep_idx, :]
    return filtered_atoms, np.array(eps_list, dtype=float), np.array(sig_list, dtype=float)


def compute_lj_energy(coords, eps, sigma, rcut=5.0, verbose_every=1000):
    """
    Compute full ds x ds LJ 12-6 energy with:
      - Lorentz-Berthelot mixing: sigma_mix = (sigma_i + sigma_j)/2
                                   eps_mix   = sqrt(eps_i * eps_j)
      - rmin = 2^(1/6) * sigma_mix
      - if i==j -> 0
      - if r < rmin -> 4*eps_mix*((sigma_mix/rmin)^12 - (sigma_mix/rmin)^6)  (NOTE: constant cap, to match MATLAB)
      - if r > rcut -> 0
      - else       -> 4*eps_mix*((sigma_mix/r)^12   - (sigma_mix/r)^6)
    Returns energy matrix (float64).
    """
    xyz = coords.astype(float)
    ds  = xyz.shape[0]

    # Pairwise distances via broadcasting
    # memory: O(N^2) – same as MATLAB; for huge N consider block-chunking
    diff = xyz[:, None, :] - xyz[None, :, :]
    r = np.linalg.norm(diff, axis=2)

    sigma_i = sigma.reshape(-1, 1)
    sigma_j = sigma.reshape(1, -1)
    eps_i   = eps.reshape(-1, 1)
    eps_j   = eps.reshape(1, -1)

    sigma_mix = (sigma_i + sigma_j) / 2.0
    eps_mix   = np.sqrt(eps_i * eps_j)
    rmin      = (2.0 ** (1.0 / 6.0)) * sigma_mix

    energy = np.zeros((ds, ds), dtype=float)

    # Masks
    eye_mask = np.eye(ds, dtype=bool)
    with np.errstate(divide="ignore", invalid="ignore"):
        # Region r <= rcut and not on diagonal
        within = (r <= rcut) & (~eye_mask)

        # Submasks
        lt_rmin = within & (r < rmin)   # use constant capped value
        ge_rmin = within & (r >= rmin)  # use normal LJ

        # Constant-capped energy for r < rmin
        if np.any(lt_rmin):
            # (sigma_mix/rmin)^n -> equals (1 / 2^(n/6)) since rmin = 2^(1/6) * sigma_mix
            # But we compute it literally to mirror MATLAB:
            ratio_const = (sigma_mix[lt_rmin] / rmin[lt_rmin])
            e_const = 4.0 * eps_mix[lt_rmin] * (ratio_const**12 - ratio_const**6)
            energy[lt_rmin] = e_const

        # Normal LJ for r >= rmin
        if np.any(ge_rmin):
            ratio = sigma_mix[ge_rmin] / r[ge_rmin]
            e = 4.0 * eps_mix[ge_rmin] * (ratio**12 - ratio**6)
            energy[ge_rmin] = e

    # i==j -> 0
    np.fill_diagonal(energy, 0.0)
    return energy


def block_sum_residual_energy(A, rcor):
    """
    Reproduce the residue-level energy block sum from MATLAB.

    A: (ds x ds) energy matrix (atom-atom)
    rcor: np.ndarray from '<ProteinName>.cor' (expected: column 9 has atom count per residue)

    For residues l,k:
      Sum A[ii, jj] over
        ii in [kalani, kalani+rcor(l,9)-1], capped to A.shape[0]
        jj in [kalanj, kalanj+rcor(k,9)-1], capped to A.shape[1]
      Then divide by 2, and zero diagonal.

    Returns resenergy: (R x R)
    """
    Asz1, Asz2 = A.shape
    # MATLAB 1-based column 9 → Python index 8
    if rcor.shape[1] < 9:
        raise ValueError(".cor file must have at least 9 columns; column 9 is atom count per residue.")
    counts = rcor[:, 8].astype(int)
    totresno = counts.shape[0]

    resenergy = np.zeros((totresno, totresno), dtype=float)

    kalani = 0  # Python 0-based
    for l in range(totresno):
        kalanj = 0
        for k in range(totresno):
            i_start = kalani
            i_end   = kalani + counts[l]  # exclusive
            j_start = kalanj
            j_end   = kalanj + counts[k]  # exclusive

            # Clamp to A bounds
            i_start_eff = min(max(i_start, 0), Asz1)
            i_end_eff   = min(max(i_end,   0), Asz1)
            j_start_eff = min(max(j_start, 0), Asz2)
            j_end_eff   = min(max(j_end,   0), Asz2)

            if (i_start_eff < i_end_eff) and (j_start_eff < j_end_eff):
                block = A[i_start_eff:i_end_eff, j_start_eff:j_end_eff]
                s = float(block.sum())
            else:
                s = 0.0

            val = 0.5 * s
            if l == k:
                val = 0.0
            resenergy[l, k] = val

            if k != totresno - 1:
                kalanj += counts[k]
        kalani += counts[l]

    return resenergy


def row_normalize_weights(resenergy, kT=1.0):
    """
    weight = exp(-E/kT); row-normalize; zero diagonal.
    Returns normalized (same shape).
    """
    W = np.exp(-resenergy / float(kT))
    # row-wise normalize
    row_sums = W.sum(axis=1, keepdims=True)
    # avoid divide-by-zero
    row_sums[row_sums == 0.0] = 1.0
    N = W / row_sums
    # zero diagonal
    np.fill_diagonal(N, 0.0)
    return N


def atomistic(ProteinName, param_file="vdw_cns.param", top_file="pdb_cns.top", rcut=5.0, kT=1.0, save_txt=True):
    """
    Main driver: replicates the MATLAB function y = atomistic(ProteinName)

    Creates '<ProteinName>_atomistic.out' with the normalized matrix (tab-separated).
    Returns the normalized matrix as a NumPy array.
    """
    # 1) Read parameter & topology maps
    vdw_map = parse_vdw_params(param_file)
    top_map = parse_topology(top_file)
    if not vdw_map:
        raise RuntimeError(f"No NONBonded entries parsed from {param_file} — check file format.")
    if not top_map:
        raise RuntimeError(f"No (res,atom)->type mappings parsed from {top_file} — check file format.")

    # 2) Load PDB atoms, filter
    pdb_path = ProteinName if ProteinName.lower().endswith(".pdb") else f"{ProteinName}.pdb"
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB not found: {pdb_path}")
    atoms = load_pdb_atoms(pdb_path)

    # 3) Map to (eps, sigma); keep only successfully mapped atoms
    atoms_f, eps, sigma = map_eps_sigma(atoms, top_map, vdw_map)
    xyz = atoms_f["xyz"]
    ds = xyz.shape[0]
    if ds == 0:
        raise RuntimeError("After mapping, zero atoms remain. Check your mapping and filters.")

    # 4) LJ energy
    energy = compute_lj_energy(xyz, eps, sigma, rcut=rcut)

    # 5) Residue-block summation using '<ProteinName>.cor'
    cor_path = f"{os.path.splitext(pdb_path)[0]}.cor"
    if not os.path.isfile(cor_path):
        raise FileNotFoundError(f".cor file not found: {cor_path} — it must be created by your 'readpdb' workflow.")
    rcor = np.loadtxt(cor_path)  # mirrors MATLAB dlmread
    if rcor.ndim == 1:
        rcor = rcor.reshape(1, -1)

    resenergy = block_sum_residual_energy(energy, rcor)

    # 6) Normalize weights row-wise (e^{-E/kT})
    normalized = row_normalize_weights(resenergy, kT=kT)

    # 7) Save
    if save_txt:
        out_path = f"{os.path.splitext(pdb_path)[0]}_atomistic.out"
        np.savetxt(out_path, normalized, fmt="%.8f", delimiter="\t")
        print(f"Saved probability matrix to: {out_path}")

    return normalized


# -------------------------- Example usage --------------------------
# Uncomment and adapt the line below in Colab once your input files are present:
# norm = atomistic("YOUR_PDB_BASENAME_OR_PATH_WITHOUT_EXT", param_file="vdw_cns.param", top_file="pdb_cns.top", rcut=5.0, kT=1.0)
# For example, if your files are:
#   ./vdw_cns.param
#   ./pdb_cns.top
#   ./4CFR.pdb
#   ./4CFR.cor
# then run:
# norm = atomistic("4CFR")
