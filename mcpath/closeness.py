# closeness.py
# Closeness based on a single "infinite" path (step-count distances)
# ---------------------------------------------------------------
# - Reads latest path_file[_N] and coor_file[_N] in current dir
# - path_file: first row is a long path of encoded nodes v = resID + chainID/10
# - Distance l_ij = minimal number of steps along this path from i to j
#                    (index difference q - p, with q > p)
# - Closeness:
#       O_i = (N-1) / sum_{j != i} l_ij
#   where N is the number of residues on the selected chain.
# - Outputs:
#   * shortest_paths_all_pairs_chain1 : rows [start end distance nodes...]
#   * closeness_float                  : "  xx.x\t0.xxxxxx"
#   * closeness_peaks                  : peak residues from closeness profile
#   * closeness_chain_labels.txt       : "'1A' , 0.096112" style
#   * closeness_chain_plot.png         : profile + peaks

import os
import re
import numpy as np
import matplotlib.pyplot as plt

# --------------------------- configuration ---------------------------
IN_PATH_BASE = "path_file"      # tab-separated; first row = long path
IN_COOR_BASE = "coor_file"      # numeric; col1=resID, col10=chainID (only for info)

PATHS_FILE      = "shortest_paths_all_pairs_chain1"
CLOSENESS_FILE  = "closeness_float"
PEAKS_FILE      = "closeness_peaks"            # peak residues
PLOT_FILE       = "closeness_chain_plot.png"   # plot for selected chain
LABELS_FILE     = "closeness_chain_labels.txt" # '1A' , 0.096112 style

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
        raise FileNotFoundError(
            f"No file found for base '{basename}' (including numbered variants)."
        )
    return os.path.join(directory, best)


def _load_path_first_row(pathfile: str) -> np.ndarray:
    """
    Load first row of path_file as a vector seq, strip zeros.
    """
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
    """
    Load coor_file (only for sanity/info; not used in distance).
    """
    coor = np.loadtxt(coorfile)
    if coor.ndim != 2 or coor.shape[1] < 10:
        raise ValueError(
            "coor_file must have at least 10 columns (got shape %s)" % (coor.shape,)
        )
    return coor


def encode_node(resid: float, chainid: float) -> float:
    # MATLAB-style encoding: v = resid + chainid/10
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
    where idx_* are indices (0-based).
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
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = i
                lookformax = False
        else:
            if this > mn + delta:
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


def _build_pos_lists(resID_seq, chainID_seq, chain_id_selected):
    """
    Build a dict: resid(int) -> sorted np.array of indices where that residue
    appears on the selected chain. Single pass over the path (O(n)).
    """
    pos_dict = {}
    for idx in range(resID_seq.size):
        if chainID_seq[idx] != chain_id_selected:
            continue
        rid = int(resID_seq[idx])
        if rid not in pos_dict:
            pos_dict[rid] = []
        pos_dict[rid].append(idx)

    # convert lists to np.array and sort (they are already increasing but sort for safety)
    for rid in pos_dict:
        pos_arr = np.asarray(pos_dict[rid], dtype=int)
        pos_arr.sort()
        pos_dict[rid] = pos_arr

    return pos_dict


def _min_forward_steps(posA, posB):
    """
    Compute minimal q - p with q>p between two sorted index arrays posA and posB,
    using a two-pointer linear scan (O(lenA + lenB)).
    Returns (best_len, best_p, best_q), where best_len=np.inf if no q>p.
    """
    i = 0
    j = 0
    lenA = posA.size
    lenB = posB.size
    best_len = np.inf
    best_p = -1
    best_q = -1

    # We want, for each p in posA, the smallest q in posB such that q>p.
    # Two-pointer strategy.
    while i < lenA and j < lenB:
        p = posA[i]
        q = posB[j]
        if q <= p:
            j += 1
        else:
            d = q - p
            if d < best_len - TOL:
                best_len = d
                best_p = p
                best_q = q
            # move to next p; maybe there's an even closer q for larger p
            i += 1

    return best_len, best_p, best_q


def _chain_to_int(ch):
    """
    Convert a chain ID from the UI to numeric:
    'A' -> 1, 'B' -> 2, ..., 'Z' -> 26
    numeric strings ('1','2',...) are passed through.
    """
    if ch is None:
        return None
    if isinstance(ch, str):
        ch = ch.strip()
        if not ch:
            return None
        # single letter like A,B,C,...
        if ch.isalpha() and len(ch) == 1:
            return ord(ch.upper()) - ord('A') + 1
        # digit(s): "1", "2", ...
        if ch.isdigit():
            return int(ch)
    # fallback: try int
    return int(ch)


def main(chain_id=None):
    global CHAIN_ID_SELECTED

    # if UI passes something, convert it
    cid = _chain_to_int(chain_id)
    if cid is not None:
        CHAIN_ID_SELECTED = cid
        
    # --- pick latest files ---
    in_pathfile = _pick_latest(IN_PATH_BASE)
    in_coorfile = _pick_latest(IN_COOR_BASE)

    print("Using files:")
    print(f"  path_file  -> {in_pathfile}")
    print(f"  coor_file  -> {in_coorfile}")

    # ------------------------ Load data ------------------------
    seq = _load_path_first_row(in_pathfile)  # long path
    n = seq.size
    assert n >= 2, "Path too short."

    # coor is not used for distances, only for info/sanity if needed
    coor = _load_coor_matrix(in_coorfile)

    # Decode nodes v = resID + chainID/10
    resID_seq   = np.floor(seq + 1e-6)
    chainID_seq = np.round((seq - resID_seq) * 10)

    # --- Build position lists in ONE pass ---
    pos_dict = _build_pos_lists(resID_seq, chainID_seq, CHAIN_ID_SELECTED)

    # --- Residues on selected chain (SORTED by residue number) ---
    res_on_chain = np.array(sorted(pos_dict.keys()), dtype=int)
    N = res_on_chain.size
    assert N > 0, f"No residues from chain {CHAIN_ID_SELECTED} found in the selected path."

    # map residue index in res_on_chain -> its pos array
    pos_list = [pos_dict[rid] for rid in res_on_chain]

    # --- Distance matrix L: minimal index difference along the path (q - p, q>p) ---
    # Also store best start/end indices in seq for reconstructing the actual subpath.
    L = np.full((N, N), np.inf, dtype=float)
    start_idx = np.full((N, N), -1, dtype=int)
    end_idx   = np.full((N, N), -1, dtype=int)

    # main double loop over residues (N typically small, e.g. <= 300)
    for a in range(N):
        posA = pos_list[a]
        if posA.size == 0:
            continue
        for b in range(N):
            if a == b:
                continue
            posB = pos_list[b]
            if posB.size == 0:
                continue

            best_len, best_p, best_q = _min_forward_steps(posA, posB)
            if best_len < np.inf:
                L[a, b] = float(best_len)
                start_idx[a, b] = best_p
                end_idx[a, b] = best_q

    # set diagonal to zero
    np.fill_diagonal(L, 0.0)

    # --- Build PATHS_FILE: [start end distance nodes...] ---
    maxNodes = 0
    for a in range(N):
        for b in range(N):
            if a == b:
                continue
            ib = start_idx[a, b]
            jb = end_idx[a, b]
            if ib >= 0 and jb >= ib:
                maxNodes = max(maxNodes, jb - ib + 1)

    totalRows = N * (N - 1)
    M = np.zeros((totalRows, 3 + maxNodes), dtype=float)
    row = 0

    for a in range(N):
        ridA = float(res_on_chain[a])
        for b in range(N):
            if a == b:
                continue
            row += 1
            ridB = float(res_on_chain[b])
            startVal = encode_node(ridA, float(CHAIN_ID_SELECTED))
            endVal   = encode_node(ridB, float(CHAIN_ID_SELECTED))
            M[row-1, 0] = startVal
            M[row-1, 1] = endVal

            ib = start_idx[a, b]
            jb = end_idx[a, b]
            if ib < 0 or jb < 0:
                M[row-1, 2] = np.inf
            else:
                seg = seq[ib:jb+1]
                M[row-1, 2] = L[a, b]
                M[row-1, 3:3+seg.size] = seg

    # Sort rows by start then end
    order = np.lexsort((M[:, 1], M[:, 0]))
    M = M[order, :]

    # --- Write PATHS_FILE with 6-decimal distance in column 3 ---
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

    # --- Closeness centrality (vectorized) ---
    # distances are directed (i -> j). For each i, sum finite l_ij, j != i.
    # Build a mask of finite entries, then zero out infinities.
    finite_mask = np.isfinite(L)
    # exclude self-distances from the sum
    np.fill_diagonal(finite_mask, False)
    L_for_sum = np.where(finite_mask, L, 0.0)
    sum_dists = np.sum(L_for_sum, axis=1)  # shape (N,)

    # (N-1)/sum_dists; where sum_dists==0 => closeness=0
    closenessVals = np.zeros(N, dtype=float)
    nonzero = sum_dists > 0
    closenessVals[nonzero] = (N - 1) / sum_dists[nonzero]

    # Encoded nodes (resID + chainID/10)
    residues_encoded = np.array(
        [encode_node(float(rid), float(CHAIN_ID_SELECTED)) for rid in res_on_chain],
        dtype=float
    )

    # --- Save CLOSENESS_FILE: "  xx.x\t0.xxxxxx" ---
    with open(CLOSENESS_FILE, "w") as fid:
        for i in range(N):
            fid.write(f" {residues_encoded[i]:4.1f}\t{closenessVals[i]:.6f}\n")

    # ----------------- Peak residues via peakdet (on selected chain) -----------------
    res_int_chain   = res_on_chain.copy()
    closeness_chain = closenessVals.copy()

    peaks_idx = []
    if closeness_chain.size > 0:
        v = closeness_chain
        std_v = float(np.std(v))
        if std_v > 0:
            try:
                maxtab, mintab = peakdet(v, std_v)
            except ValueError as e:
                print("Peak detection warning:", e)
                maxtab = np.empty((0, 2))
            if maxtab.size > 0:
                peaks_idx = maxtab[:, 0].astype(int).tolist()

    # Write peaks file
    with open(PEAKS_FILE, "w") as fid:
        fid.write("# resID\tchainID\tencoded_node\tcloseness_peak\n")
        for idx in peaks_idx:
            fid.write(
                f"{res_int_chain[idx]:4d}\t"
                f"{CHAIN_ID_SELECTED:d}\t"
                f"{residues_encoded[idx]:4.1f}\t"
                f"{closeness_chain[idx]:.6f}\n"
            )

    # ----------------- Plot & label file for closeness on chain -----------------
    if res_int_chain.size > 0:
        # chain letter (for "1A", "2A", ...)
        if 1 <= CHAIN_ID_SELECTED <= 26:
            chain_letter = chr(ord("A") + CHAIN_ID_SELECTED - 1)
        else:
            chain_letter = str(CHAIN_ID_SELECTED)

        labels = [f"{rid}{chain_letter}" for rid in res_int_chain]

        # Text file: '1A' , 0.096112
        with open(LABELS_FILE, "w") as lf:
            for lab, val in zip(labels, closeness_chain):
                lf.write(f"'{lab}' , {val:.6f}\n")

        # PLOT: all residues + peaks
        x = np.arange(len(res_int_chain))

        plt.figure(figsize=(12, 4))
        plt.plot(x, closeness_chain, marker="o", linestyle="-", label="Closeness")

        if peaks_idx:
            peak_x = np.array(peaks_idx, dtype=int)
            peak_y = closeness_chain[peak_x]
            plt.scatter(peak_x, peak_y, s=50, marker="s", label="Peaks")

        plt.xticks(x, labels, rotation=90, fontsize=6)
        plt.xlabel("Residue (number + chain)")
        plt.ylabel("Closeness centrality")
        plt.title(f"Closeness centrality along chain {chain_letter} (step-count metric)")
        if peaks_idx:
            plt.legend()
        plt.tight_layout()
        plt.savefig(PLOT_FILE, dpi=300)
        plt.close()

        print(f"Saved closeness plot to '{PLOT_FILE}' and labels to '{LABELS_FILE}'.")
    else:
        print("No residues on selected chain for plotting / label export.")

    print(f'Wrote {PATHS_FILE}, {CLOSENESS_FILE}, {PEAKS_FILE} (peaks), and {LABELS_FILE}.')

    if peaks_idx:
        print("Peak residues (resID on chain", CHAIN_ID_SELECTED, "):")
        for idx in peaks_idx:
            print(f"  {res_int_chain[idx]}  -> closeness = {closeness_chain[idx]:.6f}")
    else:
        print("No peaks found on the selected chain with delta = std(closeness).")


if __name__ == "__main__":
    main()
