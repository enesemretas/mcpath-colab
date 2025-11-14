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
    import Bio  # noqa: F401
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
            # MATLAB used tline(2:10) == 'NONBonded' (1-based),
            # which equals Python line[1:10]
            if len(line) >= 10 and line[1:10] == "NONBonded":
                seg = line.rstrip("\n")

                def field(seg, i, j):
                    # Python slices are 0-based, end-exclusive
                    return _trim(seg[i:j]) if len(seg) >= j else _trim(seg[i:])

                partype  = field(seg, 12, 16)  # 13:16 in MATLAB
                pareps_s = field(seg, 20, 26)  # 21:26
                parsig_s = field(seg, 29, 35)  # 30:35

                if partype:
                    try:
                        eps   = float(pareps_s)
                        sigma = float(parsig_s)
                        partype_to_epssigma[partype] = (eps, sigma)
                    except ValueError:
                        # tolerate malformed numerics like MATLAB script
                        pass
    return partype_to_epssigma


def parse_topology(top_path):
    """
    Parse 'pdb_cns.top' with MATLAB columns:
      parrestype = cols 1:3
      paratomtype= cols 11:13
      partypecor = cols 21:24  # (map to vdw param 'partype')
      charge     = cols 34:38  # (unused here)
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

            if res and atm and corr:
                mapping[(res.upper(), atm.upper())] = corr
    return mapping


def load_pdb_atoms(pdb_path):
    """
    Load atoms from PDB (like your pdbreader).
    Filters:
      - Only atoms with iCode='' and altLoc in {'', ' ', 'A'} (choose first conformer)
      - Skip hydrogens via (element=='H') or name rules
    Returns dict with arrays for matched heavy atoms.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", pdb_path)

    serials, atom_names, resnames, chains, resseqs = [], [], [], [], []
    xs, ys, zs = [], [], []
    elements, altlocs, icodes = [], [], []

    # choose first altLoc seen for each (chain, resseq, icode, atom_name)
    seen_altloc_key = set()

    for model in structure:
        for chain in model:
            for residue in chain:
                het, resseq, icode = residue.id
                resname = _trim(residue.get_resname()).upper()
                for atom in residue:
                    atom_name = _trim(atom.get_name()).upper()
                    element = _trim(atom.element).upper() if atom.element else ""
                    altloc = _trim(getattr(atom, "altloc", "")).upper()
                    serial = int(atom.get_serial_number())
                    key = (chain.id, resseq, icode, atom_name)

                    if key in seen_altloc_key:
                        continue
                    if altloc not in ("", " ", "A"):
                        continue
                    seen_altloc_key.add(key)

                    # hydrogen rules (match MATLAB)
                    if element == "H":
                        continue
                    if atom_name.startswith("H") and ("N" not in atom_name):
                        continue
                    if atom_name in ("D", "E", "G", "HN"):
                        continue

                    # residue normalization
                    resname_norm = resname
                    if resname_norm in ("HIP", "HIE", "HID", "HSD", "HSE", "HSP"):
                        resname_norm = "HIS"
                    elif resname_norm == "LYN":
                        resname_norm = "LYS"
                    elif resname_norm == "GLH":
                        resname_norm = "GLU"
                    elif resname_norm in ("CYM", "CYX"):
                        resname_norm = "CYS"
                    elif resname_norm == "ASH":
                        resname_norm = "ASN"
                    elif resname_norm == "TYM":
                        resname_norm = "TYR"

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
    res  = atoms["res"]
    atm  = atoms["atom"]
    elem = atoms["element"]

    # Normalize atom name for O variants
    norm_atom = []
    for a, el in zip(atm, elem):
        norm_atom.append("O" if (a in ("OXT", "O2", "O1") or el == "O") else a)

    eps_list, sig_list, keep_idx = [], [], []
    for i, (rname, aname) in enumerate(zip(res, norm_atom)):
        key   = (rname.upper(), aname.upper())
        ptype = top_map.get(key)
        if ptype is None:
            continue
        es = vdw_map.get(ptype)
        if es is None:
            continue
        eps, sig = es
        if (eps == 0.0) or (sig == 0.0):
            continue
        keep_idx.append(i)
        eps_list.append(float(eps))
        sig_list.append(float(sig))

    if not keep_idx:
        raise RuntimeError("No atoms found with valid (eps, sigma) mapping. "
                           "Check that pdb_cns.top and vdw_cns.param match your PDB.")

    keep_idx = np.array(keep_idx, dtype=int)

    # --- FIXED: build filtered atom dict without the undefined `val` name ---
    filtered_atoms = {}
    for k, v in atoms.items():
        if k == "xyz":
            # handled separately for 2D slice
            continue
        if isinstance(v, np.ndarray):
            filtered_atoms[k] = v[keep_idx]
        elif isinstance(v, list):
            filtered_atoms[k] = [v[j] for j in keep_idx]
        else:
            # metadata scalars/unused fields are copied as-is
            filtered_atoms[k] = v
    # xyz is (N,3)
    filtered_atoms["xyz"] = atoms["xyz"][keep_idx, :]

    return filtered_atoms, np.array(eps_list, dtype=float), np.array(sig_list, dtype=float)


def compute_lj_energy(coords, eps, sigma, rcut=5.0, verbose_every=1000):
    """
    Compute full ds x ds LJ 12-6 energy with:
      - Lorentz-Berthelot mixing: sigma_mix = (sigma_i + sigma_j)/2
                                   eps_mix   = sqrt(eps_i * eps_j)
      - rmin = 2^(1/6) * sigma_mix
      - if i==j -> 0
      - if r < rmin -> 4*eps_mix*((sigma_mix/rmin)^12 - (sigma_mix/rmin)^6)  (constant cap)
      - if r > rcut -> 0
      - else       -> 4*eps_mix*((sigma_mix/r)^12   - (sigma_mix/r)^6)
    """
    xyz = coords.astype(float)
    ds  = xyz.shape[0]

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
    eye_mask = np.eye(ds, dtype=bool)

    with np.errstate(divide="ignore", invalid="ignore"):
        within = (r <= rcut) & (~eye_mask)
        lt_rmin = within & (r < rmin)
        ge_rmin = within & (r >= rmin)

        if np.any(lt_rmin):
            ratio_const = (sigma_mix[lt_rmin] / rmin[lt_rmin])
            e_const = 4.0 * eps_mix[lt_rmin] * (ratio_const**12 - ratio_const**6)
            energy[lt_rmin] = e_const

        if np.any(ge_rmin):
            ratio = sigma_mix[ge_rmin] / r[ge_rmin]
            e = 4.0 * eps_mix[ge_rmin] * (ratio**12 - ratio**6)
            energy[ge_rmin] = e

    np.fill_diagonal(energy, 0.0)
    return energy


def block_sum_residual_energy(A, rcor):
    """
    Residue-level energy block sum (MATLAB replica).
    rcor: column 9 (index 8) is atom count per residue.
    """
    Asz1, Asz2 = A.shape
    if rcor.shape[1] < 9:
        raise ValueError(".cor file must have at least 9 columns; column 9 is atom count per residue.")
    counts = rcor[:, 8].astype(int)
    totresno = counts.shape[0]

    resenergy = np.zeros((totresno, totresno), dtype=float)
    kalani = 0
    for l in range(totresno):
        kalanj = 0
        for k in range(totresno):
            i_start = kalani; i_end = kalani + counts[l]
            j_start = kalanj; j_end = kalanj + counts[k]

            i_start = min(max(i_start, 0), Asz1)
            i_end   = min(max(i_end,   0), Asz1)
            j_start = min(max(j_start, 0), Asz2)
            j_end   = min(max(j_end,   0), Asz2)

            s = float(A[i_start:i_end, j_start:j_end].sum()) if (i_end > i_start and j_end > j_start) else 0.0
            resenergy[l, k] = 0.0 if l == k else 0.5 * s

            if k != totresno - 1:
                kalanj += counts[k]
        kalani += counts[l]
    return resenergy


def row_normalize_weights(resenergy, kT=1.0):
    """weight = exp(-E/kT); row-normalize; zero diagonal."""
    W = np.exp(-resenergy / float(kT))
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0.0] = 1.0
    N = W / row_sums
    np.fill_diagonal(N, 0.0)
    return N


def atomistic(ProteinName, param_file="vdw_cns.param", top_file="pdb_cns.top",
              rcut=5.0, kT=1.0, save_txt=True):
    """
    Main driver: replicates MATLAB y = atomistic(ProteinName).
    Creates '<ProteinName>_atomistic.out' and returns the normalized matrix (NumPy array).
    """
    vdw_map = parse_vdw_params(param_file)
    top_map = parse_topology(top_file)
    if not vdw_map:
        raise RuntimeError(f"No NONBonded entries parsed from {param_file}.")
    if not top_map:
        raise RuntimeError(f"No (res,atom)->type mappings parsed from {top_file}.")

    pdb_path = ProteinName if ProteinName.lower().endswith(".pdb") else f"{ProteinName}.pdb"
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB not found: {pdb_path}")
    atoms = load_pdb_atoms(pdb_path)

    atoms_f, eps, sigma = map_eps_sigma(atoms, top_map, vdw_map)
    xyz = atoms_f["xyz"]
    if xyz.shape[0] == 0:
        raise RuntimeError("After mapping, zero atoms remain. Check mapping and filters.")

    energy = compute_lj_energy(xyz, eps, sigma, rcut=rcut)

    cor_path = f"{os.path.splitext(pdb_path)[0]}.cor"
    if not os.path.isfile(cor_path):
        raise FileNotFoundError(f".cor file not found: {cor_path}")
    rcor = np.loadtxt(cor_path)
    if rcor.ndim == 1:
        rcor = rcor.reshape(1, -1)

    resenergy  = block_sum_residual_energy(energy, rcor)
    normalized = row_normalize_weights(resenergy, kT=kT)

    if save_txt:
        out_path = f"{os.path.splitext(pdb_path)[0]}_atomistic.out"
        np.savetxt(out_path, normalized, fmt="%.8f", delimiter="\t")

    return normalized

# Example (in Colab):
# norm = atomistic("4CFR", param_file="mcpath/vdw_cns.param", top_file="mcpath/pdb_cns.top", rcut=5.0, kT=1.0)
