# infinite.py
# Python translation of MATLAB "infinite" function
# (single long path; same logic as the original code)

import os
import numpy as np


def infinite(ProteinName: str, path_length: int, pottype: str = "1"):
    """
    Generate a long residue path using the normalized matrix in
    <ProteinName><Type>.out and coordinates in <ProteinName>.cor

    Parameters
    ----------
    ProteinName : str
        Base name of the files (no extension).
    path_length : int
        Number of steps to generate (MATLAB: pathlength).
    pottype : {'1','2','3','4','5'}
        1='atomistic', 2='BJ', 3='TD', 4='dist', 5='markov'.

    Returns
    -------
    path : np.ndarray, shape (1, L)
        Encoded path values: resID + chainID/10
    """

    # ----- map pottype -> Type suffix (as in MATLAB) -----
    type_map = {
        "1": "_atomistic",
        "2": "_BJ",
        "3": "_TD",
        "4": "_dist",
        "5": "_markov",
    }
    if pottype not in type_map:
        raise ValueError("pottype must be one of '1','2','3','4','5'")
    Type = type_map[pottype]

    # ----- load inputs -----
    cor_file = f"{ProteinName}.cor"
    out_file = f"{ProteinName}{Type}.out"

    if not os.path.isfile(cor_file):
        raise FileNotFoundError(cor_file)
    if not os.path.isfile(out_file):
        raise FileNotFoundError(out_file)

    # coor: same usage as MATLAB: col1=resID, col10=chainID
    coor = np.loadtxt(cor_file)
    normalized = np.loadtxt(out_file)

    # MATLAB defaults:
    baslangic = 1    # residue number (unused for start in MATLAB infinite)
    baslangic1 = 1   # chain number
    pathlength = int(path_length)

    # size of probability matrix
    resno = normalized.shape[0]
    assert normalized.shape == (resno, resno)

    # ----- build newmtrx (resno+1 x resno+1) with header row/col -----
    # MATLAB:
    # for i=1:resno+1
    #   newmtrx(1,i)=i-1;
    #   newmtrx(i,1)=i-1;
    # end
    newmtrx = np.zeros((resno + 1, resno + 1), dtype=float)
    for i in range(resno + 1):
        newmtrx[0, i] = i  # 0..resno
        newmtrx[i, 0] = i
    # newmtrx(i+1,j+1)=normalized(i,j)
    newmtrx[1:, 1:] = normalized

    upper = 1.0
    lower = 0.0

    # cap matrix (not used later, but we keep it for completeness)
    cap = np.ones((resno, resno), dtype=float)

    # ----- NEW MATRIX GENERATION -----
    # We follow the MATLAB loops explicitly.

    if int(pottype) == 1:
        # atomistic
        for i in range(1, resno + 1):  # MATLAB: i=2:resno+1
            # sum over columns 2..resno+1
            row_sum = 0.0
            for j in range(1, resno + 1):
                row_sum += newmtrx[i, j]

            if row_sum == 0:
                continue

            # first normalization + diagonal = 1
            for j in range(1, resno + 1):
                newmtrx[i, j] = newmtrx[i, j] / row_sum
                if i == j:
                    newmtrx[i, j] = 1.0

            # find minimum over whole row (including col 0)
            row_min = np.min(newmtrx[i, :])

            # zero minimum and diagonal
            for j in range(1, resno + 1):
                if newmtrx[i, j] == row_min:
                    newmtrx[i, j] = 0.0
                if i == j:
                    newmtrx[i, j] = 0.0

            # renormalize (excluding col 0)
            row_sum = 0.0
            for j in range(1, resno + 1):
                row_sum += newmtrx[i, j]

            if row_sum == 0:
                continue

            for j in range(1, resno + 1):
                newmtrx[i, j] = newmtrx[i, j] / row_sum
                if i == j:
                    newmtrx[i, j] = 0.0

    else:
        # BJ / TD / dist / markov
        for i in range(1, resno + 1):
            # original sum (unnormalized) over j=2..resno+1
            row_sum = 0.0
            for j in range(1, resno + 1):
                row_sum += newmtrx[i, j]

            if row_sum == 0:
                continue

            average = row_sum / resno

            # normalize, set diagonal=1
            for j in range(1, resno + 1):
                newmtrx[i, j] = newmtrx[i, j] / row_sum
                if i == j:
                    newmtrx[i, j] = 1.0

            # threshold by average (note: average from *pre*-normalization)
            for j in range(1, resno + 1):
                if newmtrx[i, j] < average:
                    newmtrx[i, j] = 0.0
                if i == j:
                    newmtrx[i, j] = 0.0

            # renormalize and fill cap
            row_sum = 0.0
            for j in range(1, resno + 1):
                row_sum += newmtrx[i, j]

            if row_sum == 0:
                continue

            counter = 0
            for j in range(1, resno + 1):
                newmtrx[i, j] = newmtrx[i, j] / row_sum
                if i == j:
                    newmtrx[i, j] = 0.0
                if newmtrx[i, j] != 0:
                    counter += 1
                    cap[i - 1, counter - 1] = newmtrx[i, j]

    # ----- PATH GENERATION -----
    # MATLAB:
    # startres=1; runs=runs+1; newmtrx=permtrx; counter=1; sequence(counter,runs)=startres; ...

    startres = 1  # 1-based residue index
    taban = max(1, pathlength // 10)
    rng = np.random.default_rng()

    sequence = np.zeros(pathlength, dtype=int)   # 1-based residue indices
    path = np.zeros((1, pathlength), dtype=float)

    # first step
    counter = 1
    sequence[counter - 1] = startres
    row0 = startres - 1  # convert to 0-based index into coor
    path[0, counter - 1] = float(coor[row0, 0]) + float(coor[row0, 9]) / 10.0

    while counter != pathlength:
        counter += 1

        # --- build P (indices & probabilities) ---
        # MATLAB:
        # j=0; P=0; for i=2:resno+1, if newmtrx(startres+1,i)>lower && <upper ...
        j = 0
        P_idx = []
        P_prob = []

        row_idx = startres  # since MATLAB row = startres+1 (1-based) => index 'startres'
        for col in range(1, resno + 1):  # 1..resno corresponds to col 2..resno+1
            val = newmtrx[row_idx, col]
            if (val > lower) and (val < upper):
                j += 1
                P_idx.append(int(newmtrx[0, col]))   # header = residue index
                P_prob.append(float(val))

        if j == 0:
            # MATLAB behavior: no update to startres; path stays at same residue
            sequence[counter - 1] = startres
        else:
            P_prob = np.array(P_prob, dtype=float)
            # Q(2,i) = P(2,i) / sum; cumulative in gen(2,i)
            s = P_prob.sum()
            if s <= 0:
                # should not happen; fall back to uniform among P_idx
                q = np.ones(j, dtype=float) / j
            else:
                q = P_prob / s

            cdf = np.cumsum(q)
            r = rng.random()
            k = int(np.searchsorted(cdf, r, side="right"))
            if k >= j:
                k = j - 1
            startres = int(P_idx[k])     # 1-based residue index

            sequence[counter - 1] = startres

        # append encoded residue to path
        row0 = startres - 1
        path[0, counter - 1] = float(coor[row0, 0]) + float(coor[row0, 9]) / 10.0

        if counter % taban == 0:
            print(f"Number of steps generated so far: {counter}")

    # ----- trim trailing zero columns (mirrors MATLAB, mostly no-op here) -----
    # MATLAB does a check for max(abs(path(:,pl+1-i)))~=0; here path has exactly
    # 'pathlength' entries, all filled, so we just keep them.

    pl = path.shape[1]

    # ----- write .path file -----
    out_name = f"{ProteinName}{Type}_{pl}steps_infinite.path"
    np.savetxt(out_name, path, fmt="%.10g", delimiter="\t")
    print(f"Saved path to: {out_name}")

    # ----- stats block (same idea as MATLAB; optional) -----
    coor_enc = coor[:, 0] + coor[:, 9] / 10.0
    stats = np.column_stack([coor_enc, np.zeros_like(coor_enc)])
    for val in path[0, :]:
        hits = np.where(np.isclose(stats[:, 0], val, atol=1e-12))[0]
        if hits.size > 0:
            stats[hits[0], 1] += 1

    return path
