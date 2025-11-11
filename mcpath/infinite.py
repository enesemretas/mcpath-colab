# Infinite path generator — Python/NumPy version (Colab-ready)
# Save this cell and run. Then call:
#   path = infinite("4CFR", path_length=500_000, pottype='1')
#
# Files expected in the working dir:
#   <ProteinName>.cor
#   <ProteinName>_<tip>.out  where tip from pottype:
#       '1' -> _atomistic, '2' -> _BJ, '3' -> _TD, '4' -> _dist, '5' -> _markov
#
# Output written:
#   <ProteinName>_<tip>_<steps>steps_infinite.path   (tab-delimited)

import os
import numpy as np

def infinite(ProteinName: str, path_length: int, pottype: str = '1'):
    """
    Generate a single long residue path using the normalized transition
    matrix in <ProteinName>_<tip>.out and residue info from <ProteinName>.cor.

    Parameters
    ----------
    ProteinName : str
        Base name of files (no extension).
    path_length : int
        Number of steps to generate (e.g., 500_000).
    pottype : {'1','2','3','4','5'}
        1='atomistic', 2='BJ', 3='TD', 4='dist', 5='markov'

    Returns
    -------
    path : np.ndarray (shape: 1 x path_length)
        Encoded path values: resID + chainID/10
    """

    # --- map pottype to suffixes (kept identical to your MATLAB)
    tip_map = {
        '1': ('_atomistic', 'atomistic'),
        '2': ('_BJ',        'BJ'),
        '3': ('_TD',        'TD'),
        '4': ('_dist',      'dist'),
        '5': ('_markov',    'markov'),
    }
    if pottype not in tip_map:
        raise ValueError("pottype must be one of '1','2','3','4','5'")

    Type, tip = tip_map[pottype]

    # --- load inputs
    cor_file = f"{ProteinName}.cor"
    out_file = f"{ProteinName}{Type}.out"

    if not os.path.isfile(cor_file):
        raise FileNotFoundError(f"Missing {cor_file}")
    if not os.path.isfile(out_file):
        raise FileNotFoundError(f"Missing {out_file}")

    # coor: numeric matrix; uses columns 1 and 10 in MATLAB (-> 0 and 9 in Python)
    coor = np.loadtxt(cor_file)
    # normalized: probability/weight matrix (resno x resno)
    normalized = np.loadtxt(out_file)

    # --- defaults replicated
    baslangic  = 1   # residue ID (as in coor column 1 in MATLAB)
    baslangic1 = 1   # chain ID number (as in coor column 10 in MATLAB)
    pathlength = int(path_length)

    resno = normalized.shape[0]
    assert normalized.shape[1] == resno, "normalized must be square"

    # --- find the initial row index (MATLAB set 'initial' but then used startres=1;
    #     here we prefer to actually start from the found index if it exists)
    initial = None
    # In MATLAB: coor(i,10)==baslangic1 && coor(i,1)==baslangic
    # Python (0-based): coor[i, 9] and coor[i, 0].
    for i in range(resno):
        if int(coor[i, 9]) == int(baslangic1) and int(coor[i, 0]) == int(baslangic):
            initial = i
            break
    if initial is None:
        # fallback to first row if not found (MATLAB effectively started at 1 anyway)
        initial = 0

    # --- build newmtrx like MATLAB (size: (resno+1) x (resno+1), first row/col are headers 0..resno-1)
    newmtrx = np.zeros((resno + 1, resno + 1), dtype=float)
    header = np.arange(0, resno, dtype=float)
    newmtrx[0, 1:] = header
    newmtrx[1:, 0] = header
    newmtrx[1:, 1:] = np.array(normalized, dtype=float)

    upper = 1.0  # used as an upper threshold in selection
    # Note: cap matrix is written in MATLAB but never used later; we keep structure.
    cap = np.ones((resno, resno), dtype=float)

    # --- NEW MATRIX GENERATION (row normalization / thresholding) ---
    # MATLAB loops by row (i=2..resno+1); 'sum' is re-used per row
    if pottype == '1':
        # atomistic branch
        for i in range(1, resno + 1):
            row_sum = np.sum(newmtrx[i, 1:])
            if row_sum <= 0:
                # avoid divide-by-zero (if an all-zero row)
                continue
            newmtrx[i, 1:] = newmtrx[i, 1:] / row_sum
            # set diagonal to 1 temporarily
            newmtrx[i, i] = 1.0

            # zero out minimum entry in the row and diagonal
            row_min = np.min(newmtrx[i, :])
            # minimum among all entries, MATLAB uses entire row; keep same
            newmtrx[i, newmtrx[i, :] == row_min] = 0.0
            newmtrx[i, i] = 0.0

            # renormalize after zeroing
            row_sum = np.sum(newmtrx[i, 1:])
            if row_sum > 0:
                newmtrx[i, 1:] = newmtrx[i, 1:] / row_sum
                newmtrx[i, i] = 0.0
        lower = 0.0
    else:
        # BJ / TD / dist / markov branch
        for i in range(1, resno + 1):
            row_sum = np.sum(newmtrx[i, 1:])
            if row_sum <= 0:
                continue
            average = row_sum / resno

            newmtrx[i, 1:] = newmtrx[i, 1:] / row_sum
            newmtrx[i, i] = 1.0

            # threshold: values below average -> 0; diagonal -> 0
            below = newmtrx[i, :] < average
            newmtrx[i, below] = 0.0
            newmtrx[i, i] = 0.0

            # renormalize to sum=1 (excluding col 0 header)
            row_sum = np.sum(newmtrx[i, 1:])
            if row_sum > 0:
                newmtrx[i, 1:] = newmtrx[i, 1:] / row_sum
                newmtrx[i, i] = 0.0

            # (cap filling was used to store non-zeros; not used later)
            nz = newmtrx[i, 1:] > 0
            nz_vals = newmtrx[i, 1:][nz]
            count = nz_vals.size
            if count > 0:
                cap[i - 1, :count] = nz_vals
        lower = 0.0

    # --- PATH GENERATION (single path, like the MATLAB script effectively does) ---
    rng = np.random.default_rng()  # modern RNG
    taban = max(1, pathlength // 10)

    # MATLAB started with 'startres=1' (1-based). We use 'initial' found above (0-based).
    startres = int(initial)  # index in [0..resno-1]

    sequence = np.zeros((pathlength,), dtype=int)
    path = np.zeros((1, pathlength), dtype=float)

    # step 1
    sequence[0] = startres
    # encode resID + chainID/10  (MATLAB: coor(row,1)+coor(row,10)/10)
    path[0, 0] = float(coor[startres, 0]) + float(coor[startres, 9]) / 10.0

    for step in range(1, pathlength):
        # Collect neighbors from row (startres+1) and cols 2..resno+1
        probs = newmtrx[startres + 1, 1:].copy()
        # zero diagonal and any values outside (lower, upper)
        probs[startres] = 0.0
        mask = (probs > lower) & (probs < upper)
        candidates = np.nonzero(mask)[0]  # indices 0..resno-1
        if candidates.size == 0:
            # If we hit a dead-end row, fallback to uniform among non-diagonal indices
            candidates = np.delete(np.arange(resno), startres)
            probs = np.ones_like(candidates, dtype=float) / max(1, candidates.size)
            next_idx = int(rng.choice(candidates, p=probs))
        else:
            p_raw = probs[candidates]
            s = p_raw.sum()
            if s <= 0:
                # safeguard: uniform among candidates
                p = np.ones_like(candidates, dtype=float) / candidates.size
            else:
                p = p_raw / s
            next_idx = int(rng.choice(candidates, p=p))

        startres = next_idx
        sequence[step] = startres
        path[0, step] = float(coor[startres, 0]) + float(coor[startres, 9]) / 10.0

        if (step + 1) % taban == 0:
            print(f"Number of steps generated so far: {step + 1}")

    # --- Trim trailing all-zero columns (mirror of MATLAB's end-trim pass)
    # (Here we always generate exactly 'pathlength' steps; kept for parity.)
    # If any trailing zeros exist (shouldn't unless pathlength==0), drop them.
    last_nonzero = np.max(np.nonzero(np.abs(path[0]) > 0)) if np.any(path) else -1
    if last_nonzero >= 0:
        path = path[:, :last_nonzero + 1]

    # --- Write .path (tab-delimited)
    pl = path.shape[1]
    out_name = f"{ProteinName}{Type}_{pl}steps_infinite.path"
    np.savetxt(out_name, path, fmt="%.10g", delimiter="\t")
    print(f"Saved path to: {out_name}")

    # --- stats (same idea as MATLAB block)
    coor_enc = coor[:, 0] + coor[:, 9] / 10.0
    stats = np.column_stack([coor_enc, np.zeros_like(coor_enc)])
    # count visits for first (and only) path row
    for val in path[0, :]:
        hits = np.where(np.isclose(stats[:, 0], val, atol=1e-12))[0]
        if hits.size:
            stats[hits[0], 1] += 1

    return path
