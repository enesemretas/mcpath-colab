# closeness.py
# Direct translation of MATLAB "closeness" function
# - Reads latest path_file[_N], coor_file[_N], atom_file[_N] in current dir
# - Computes shortest_paths_all_pairs_chain1 and closeness_float
# - Uses a MATLAB-like PEAKDET algorithm to find peaks in closeness
#   with delta = std(closeness_vals) on the selected chain.
# - Creates a plot of residue labels ('1A', '2A', ...) vs closeness
#   and a label file with lines like: '1A' , 0.014814

import os
import re
import numpy as np
import matplotlib.pyplot as plt

# --------------------------- configuration ---------------------------
IN_PATH_BASE = "path_file"      # tab-separated; first row = long path
IN_COOR_BASE = "coor_file"      # numeric; col1=resID, cols6-8=CA xyz, col10=chainID
IN_ATOM_BASE = "atom_file"      # optional, just reported

PATHS_FILE      = "shortest_paths_all_pairs_chain1"
CLOSENESS_FILE  = "closeness_float"
PEAKS_FILE      = "closeness_peaks"            # output for peak residues
PLOT_FILE       = "closeness_chain_plot.png"   # plot for selected chain
LABELS_FILE     = "closeness_chain_labels.txt" # '1A' , 0.014814 style

CHAIN_ID_SELECTED = 1
TOL = 1e-9
# --------------------------------------------------------------------


def _pick_latest(basename: str, directory: str = ".") -> str:
    """
    Pick the latest file among basename, basename_2, basename_3, ... in `directory`.
    'Latest' = highest numeric suffix (no suffix = 0).
    """
    pattern = re.compile(rf"^{re.escape(basename)}(?:_(\d+))?$")
    best = None
    best_k = -1
    for fname in os.listdir(directory):
        m = pattern.match(fname)
        if m:
            k = int(m.group(1)) if m.group(1) is not None else 0
            if k > best_k:
                best_k = k
                best = fname
    if best is None:
        raise FileNotFoundError(f"No file found for base '{basename}' (including numbered variants).")
    return os.path.join(directory, best)


def _load_path_first_row(pathfile: str) -> np.ndarray:
    # MATLAB: P = readmatrix(...); P(~isfinite(P))=0; seq = P(1,:); seq = seq(seq~=0);
    P = np.loadtxt(pathfile, delimiter="\t", ndmin=2)
    if P.ndim == 1:
        P = P.reshape(1, -1)
    P[~np.isfinite(P)] = 0.0
    seq = P[0, :]
    seq = seq[seq != 0]
    if seq.size < 2:
        raise ValueError("Path too short (need at least 2 nodes).")
    return seq


def _load_coor_matrix(coorfile: str) -> np.ndarray:
    coor = np.loadtxt(coorfile)
    if coor.ndim != 2 or coor.shape[1] < 10:
        raise ValueError("coor_file must have at least 10 columns (got shape %s)" % (coor.shape,))
    return coor


def encode_node(resid: float, chainid: float) -> float:
    # MATLAB: encode = @(rid,cid) rid + cid/10;
    return resid + chainid / 10.0


# ------------------- MATLAB-like PEAKDET implementation -------------------
def peakdet(v, delta):
    """
    Python version of the MATLAB 'peakdet' function by Eli Billauer.

    v: 1D array-like of values
    delta: positive scalar threshold
    Returns:
      maxtab: array([[idx_max, value_max], ...])
      mintab: array([[idx_min, value_min], ...])
    where idx_* are indices (0-based, like Python).
    """
    v = np.asarray(v, dtype=float).flatten()
    if v.size == 0:
        return np.empty((0, 2)), np.empty((0, 2))

    if np.ndim(delta) != 0:
        raise ValueError("Input argument DELTA must be a scalar")
    if delta <= 0:
        raise ValueError("Input argument DELTA must be positive")

    maxtab = []
    mintab = []

    mn = np.inf
    mx = -np.inf
    mnpos = np.nan
    mxpos = np.nan

    lookformax = True

    for i in range(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = i
        if this < mn:
            mn = this
            mnpos = i

        if lookformax:
            if this < mx - delta:
                # maximum peak found
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = i
                lookformax = False
        else:
            if this > mn + delta:
                # minimum peak found
                mintab.append((mnpos, mn))
                mx = this
                mxpos = i
                lookformax = True

    if len(maxtab) == 0:
        maxtab_arr = np.empty((0, 2))
    else:
        maxtab_arr = np.asarray(maxtab, dtype=float)

    if len(mintab) == 0:
        mintab_arr = np.empty((0, 2))
    else:
        mintab_arr = np.asarray(mintab, dtype=float)

    return maxtab_arr, mintab_arr
# -------------------------------------------------------------------------


def main():
    # --- pick latest files (MATLAB had fixed names; we extend to numbered versions) ---
    in_pathfile = _pick_latest(IN_PATH_BASE)
    in_coorfile = _pick_latest(IN_COOR_BASE)
    try:
        in_atomfile = _pick_latest(IN_ATOM_BASE)
    except FileNotFoundError:
        in_atomfile = None

    print("Using files:")
    print(f"  path_file  -> {in_pathfile}")
    print(f"  coor_file  -> {in_coorfile}")
    if in_atomfile:
        print(f"  atom_file  -> {in_atomfile}")

    # ------------------------ Load data (MATLAB block) ------------------------
    seq = _load_path_first_row(in_pathfile)
    n = seq.size
    assert n >= 2, "Path too short."

    coor = _load_coor_matrix(in_coorfile)
    resID_list   = coor[:, 0]        # col1
    chainID_list = coor[:, 9]        # col10
    CAxyz        = coor[:, 5:8]      # col6-8

    # Decode nodes v = resID + chainID/10
    resID_seq   = np.floor(seq + 1e-6)
    chainID_seq = np.round((seq - resID_seq) * 10)

    # --- Residues on selected chain (SORTED) ---
    mask_chain = (chainID_seq == CHAIN_ID_SELECTED)
    res_on_chain = np.unique(resID_seq[mask_chain])
    Nchain = res_on_chain.size
    assert Nchain > 0, f"No residues from chain {CHAIN_ID_SELECTED} found in the selected path."

    # --- Precompute step geometry & validity (CA distances) ---
    keyCoor = (resID_list.astype(np.int64) * 100) + chainID_list.astype(np.int64)
    coorIdxMap = {int(k): int(i) for i, k in enumerate(keyCoor)}

    dist = np.zeros(n - 1, dtype=float)
    invalid = np.zeros(n - 1, dtype=bool)

    for k in range(n - 1):  # MATLAB k = 1:n-1
        ra = resID_seq[k];   ca = chainID_seq[k]
        rb = resID_seq[k+1]; cb = chainID_seq[k+1]
        k1 = int(ra) * 100 + int(ca)
        k2 = int(rb) * 100 + int(cb)
        if k1 not in coorIdxMap or k2 not in coorIdxMap:
            invalid[k] = True
            dist[k] = 0.0
        else:
            va = CAxyz[coorIdxMap[k1], :]
            vb = CAxyz[coorIdxMap[k2], :]
            dist[k] = float(np.linalg.norm(vb - va))

    # S = [0; cumsum(dist)], C = [0; cumsum(invalid)]
    S = np.zeros(n, dtype=float)
    C = np.zeros(n, dtype=int)
    if n > 1:
        S[1:] = np.cumsum(dist)
        C[1:] = np.cumsum(invalid.astype(int))
    maxC = int(np.max(C))

    # --- For each residue, indices where it occurs in the path ---
    occIdx = {}
    for rid in res_on_chain:
        rid_val = float(rid)
        occ = np.where((resID_seq == rid_val) & (chainID_seq == CHAIN_ID_SELECTED))[0]
        occIdx[float(rid_val)] = occ

    # --- Fast shortest search for all ordered pairs (i != j) ---
    res_on_chain = res_on_chain.reshape(1, -1).ravel()
    N = res_on_chain.size

    bestI   = np.full((N, N), np.nan, dtype=float)
    bestJ   = np.full((N, N), np.nan, dtype=float)
    bestD   = np.full((N, N), np.inf, dtype=float)
    bestLen = np.zeros((N, N), dtype=int)

    dBestForJ = np.full(n, np.inf, dtype=float)
    iBestForJ = np.full(n, np.nan, dtype=float)
    lenForJ   = np.zeros(n, dtype=int)

    bestS_forC = np.full(maxC + 1, -np.inf, dtype=float)  # index by C(t) (0..maxC)
    bestI_forC = np.full(maxC + 1, np.nan, dtype=float)

    for ai in range(N):  # MATLAB ai = 1:N
        sRID = float(res_on_chain[ai])
        sIdx = occIdx[sRID]  # sorted indices where sRID occurs

        # Reset per-j results
        dBestForJ[:] = np.inf
        iBestForJ[:] = np.nan
        lenForJ[:]   = 0
        bestS_forC[:] = -np.inf
        bestI_forC[:] = np.nan

        # Single sweep across the whole path
        p = 0  # pointer into sIdx (0-based)
        sIdx_len = sIdx.size

        for t in range(n):  # MATLAB t = 1:n
            # If t is an eligible start index for sRID, update tables at its C(t)
            while p < sIdx_len and sIdx[p] == t:
                c = int(C[t])   # NOTE: we use C(t) (0..maxC), not C(t)+1
                s_val = S[t]
                ib = bestI_forC[c]
                sb = bestS_forC[c]
                if (s_val > sb + TOL) or (abs(s_val - sb) <= TOL and (np.isnan(ib) or t > ib)):
                    bestS_forC[c] = s_val
                    bestI_forC[c] = float(t)
                p += 1

            # If t can be an end (on selected chain), try best compatible start
            if chainID_seq[t] == CHAIN_ID_SELECTED:
                c = int(C[t])
                ib = bestI_forC[c]
                if not np.isnan(ib):
                    ib_int = int(ib)
                    d   = S[t] - bestS_forC[c]
                    length = t - ib_int + 1
                    if (d + TOL < dBestForJ[t]) or (abs(d - dBestForJ[t]) <= TOL and length < lenForJ[t]):
                        dBestForJ[t] = d
                        iBestForJ[t] = float(ib_int)
                        lenForJ[t]   = length

        # For each end residue b, pick best among its j occurrences
        for bi in range(N):
            if ai == bi:
                continue
            eRID = float(res_on_chain[bi])
            js = occIdx[eRID]
            if js.size == 0:
                continue

            dvals = dBestForJ[js]
            pos = int(np.argmin(dvals))
            dmin = dvals[pos]
            if np.isfinite(dmin):
                jStar = js[pos]
                # tie by shortest len
                ties = np.where(np.abs(dBestForJ[js] - dmin) <= TOL)[0]
                if ties.size > 1:
                    tie_js   = js[ties]
                    tie_lens = lenForJ[tie_js]
                    kmin = int(np.argmin(tie_lens))
                    jStar = tie_js[kmin]
                    dmin  = dBestForJ[jStar]
                bestD[ai, bi]   = float(dmin)
                bestI[ai, bi]   = float(iBestForJ[jStar])
                bestJ[ai, bi]   = float(jStar)
                bestLen[ai, bi] = int(lenForJ[jStar])

    # --- Build output matrix for paths: [start end distance nodes...] ---
    maxNodes = 0
    for a in range(N):
        for b in range(N):
            if a == b:
                continue
            if np.isnan(bestI[a, b]) or np.isnan(bestJ[a, b]):
                continue
            ib = int(bestI[a, b])
            jb = int(bestJ[a, b])
            maxNodes = max(maxNodes, jb - ib + 1)

    totalRows = N * (N - 1)
    M = np.zeros((totalRows, 3 + maxNodes), dtype=float)
    row = 0

    for a in range(N):
        for b in range(N):
            if a == b:
                continue
            row += 1  # MATLAB 1-based row; Python row-1 index
            startVal = encode_node(float(res_on_chain[a]), float(CHAIN_ID_SELECTED))
            endVal   = encode_node(float(res_on_chain[b]), float(CHAIN_ID_SELECTED))
            M[row-1, 0] = startVal
            M[row-1, 1] = endVal

            if np.isnan(bestI[a, b]) or np.isnan(bestJ[a, b]):
                M[row-1, 2] = np.inf
            else:
                ib = int(bestI[a, b])
                jb = int(bestJ[a, b])
                seg = seq[ib:jb+1]
                M[row-1, 2] = bestD[a, b]
                M[row-1, 3:3+seg.size] = seg

    # Sort rows by start then end
    order = np.lexsort((M[:, 1], M[:, 0]))
    M = M[order, :]

    # --- Write paths file with 6-decimal distance in column 3 ---
    with open(PATHS_FILE, "w") as fid:
        for r in range(M.shape[0]):
            start_node = M[r, 0]
            end_node   = M[r, 1]
            dist_val   = M[r, 2]
            fid.write(f" {start_node:4.1f}\t{end_node:4.1f}\t")
            if np.isfinite(dist_val):
                fid.write(f"{dist_val:.6f}")
            else:
                fid.write("Inf")
            tail = M[r, 3:]
            nz = np.nonzero(tail != 0)[0]
            if nz.size > 0:
                tail = tail[:nz[-1] + 1]
                for v in tail:
                    fid.write(f"\t{v:4.1f}")
            fid.write("\n")

    # --- Build symmetric "shortest of the pair" matrix for closeness ---
    start_nodes = M[:, 0]
    end_nodes   = M[:, 1]
    distances   = M[:, 2]

    residues = np.unique(np.concatenate([start_nodes, end_nodes]))
    Nres = residues.size
    idxMap = {str(resid): i for i, resid in enumerate(residues)}

    distMat = np.full((Nres, Nres), np.inf, dtype=float)
    for k in range(distances.size):
        d = distances[k]
        if not np.isfinite(d):
            continue
        i = idxMap[str(start_nodes[k])]
        j = idxMap[str(end_nodes[k])]
        if d < distMat[i, j]:
            distMat[i, j] = d
        if d < distMat[j, i]:
            distMat[j, i] = d
    np.fill_diagonal(distMat, 0.0)
    distMat = np.minimum(distMat, distMat.T)

    # --- Closeness centrality ---
    closenessVals = np.zeros(Nres, dtype=float)
    for i in range(Nres):
        rowVals = distMat[i, :]
        finiteVals = rowVals[np.isfinite(rowVals)]
        s = float(np.sum(finiteVals))
        if s > 0:
            closenessVals[i] = (Nres - 1) / s
        else:
            closenessVals[i] = 0.0

    # --- Save EXACT format: "  xx.x\t0.xxxxxx" ---
    with open(CLOSENESS_FILE, "w") as fid:
        for i in range(Nres):
            fid.write(f" {residues[i]:4.1f}\t{closenessVals[i]:.6f}\n")

    # ----------------- Find peak residues via peakdet -----------------
    # Decode residues -> (resID, chainID)
    res_int  = np.floor(residues + 1e-6).astype(int)
    chain_id = np.round((residues - res_int) * 10).astype(int)

    # Focus only on selected chain
    mask_sel_chain = (chain_id == CHAIN_ID_SELECTED)
    res_int_chain        = res_int[mask_sel_chain]
    closeness_chain      = closenessVals[mask_sel_chain]
    residues_chain_enc   = residues[mask_sel_chain]

    # Sort along sequence (residue number)
    order_chain = np.argsort(res_int_chain)
    res_int_chain      = res_int_chain[order_chain]
    closeness_chain    = closeness_chain[order_chain]
    residues_chain_enc = residues_chain_enc[order_chain]

    peaks_idx = []

    if res_int_chain.size > 0:
        v = closeness_chain
        std_v = float(np.std(v))
        if std_v > 0:
            delta = std_v  # PEAK DELTA = STANDARD DEVIATION
            try:
                maxtab, mintab = peakdet(v, delta)
            except ValueError as e:
                print("Peak detection warning:", e)
                maxtab = np.empty((0, 2))

            if maxtab.size > 0:
                # maxtab[:,0] are indices (0-based)
                peaks_idx = maxtab[:, 0].astype(int).tolist()
        else:
            # zero variance -> flat signal, no peaks
            maxtab = np.empty((0, 2))

    # Write peaks to file
    with open(PEAKS_FILE, "w") as fid:
        fid.write("# resID\tchainID\tencoded_node\tcloseness_peak\n")
        for idx in peaks_idx:
            fid.write(
                f"{res_int_chain[idx]:4d}\t"
                f"{CHAIN_ID_SELECTED:d}\t"
                f"{residues_chain_enc[idx]:4.1f}\t"
                f"{closeness_chain[idx]:.6f}\n"
            )

    # ----------------- Plot & label file for closeness on chain -----------------
    if res_int_chain.size > 0:
        # Map numeric chain ID to letter if possible
        if 1 <= CHAIN_ID_SELECTED <= 26:
            chain_letter = chr(ord("A") + CHAIN_ID_SELECTED - 1)
        else:
            chain_letter = str(CHAIN_ID_SELECTED)

        labels = [f"{rid}{chain_letter}" for rid in res_int_chain]

        # Text file: '1A' , 0.014814
        with open(LABELS_FILE, "w") as lf:
            for lab, val in zip(labels, closeness_chain):
                lf.write(f"'{lab}' , {val:.6f}\n")

        # PLOT: all residues + peaks
        x = np.arange(len(res_int_chain))

        plt.figure(figsize=(12, 4))
        # full profile
        plt.plot(x, closeness_chain, marker="o", linestyle="-", label="Closeness")

        # overlay peaks if any
        if peaks_idx:
            peak_x = np.array(peaks_idx, dtype=int)
            peak_y = closeness_chain[peak_x]
            plt.scatter(peak_x, peak_y, s=50, marker="s", label="Peaks")  # larger markers

        plt.xticks(x, labels, rotation=90, fontsize=6)
        plt.xlabel("Residue (number + chain)")
        plt.ylabel("Closeness centrality")
        plt.title(f"Closeness centrality along chain {chain_letter}")
        if peaks_idx:
            plt.legend()
        plt.tight_layout()
        plt.savefig(PLOT_FILE, dpi=300)
        plt.close()

        print(f"Saved closeness plot to '{PLOT_FILE}' and labels to '{LABELS_FILE}'.")
    else:
        print("No residues on selected chain for plotting / label export.")

    print(f'Wrote {PATHS_FILE} (distances with 6 decimals), '
          f'{CLOSENESS_FILE} (closeness with 6 decimals), '
          f'{PEAKS_FILE} (peak residues via peakdet, delta = std(closeness)).')

    if peaks_idx:
        print("Peak residues (resID on chain", CHAIN_ID_SELECTED, "):")
        for idx in peaks_idx:
            print(f"  {res_int_chain[idx]}  -> closeness = {closeness_chain[idx]:.6f}")
    else:
        print("No peaks found on the selected chain using peakdet with delta = std(closeness).")


if __name__ == "__main__":
    main()
