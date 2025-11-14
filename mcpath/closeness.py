# closeness.py
# Closeness based on a single "infinite" path (step-count distances)
# ---------------------------------------------------------------
# - Reads latest path_file[_N] and coor_file[_N] in current dir
# - path_file: first row is a long path of encoded nodes v = resID + chainID/10
# - coor_file: columns 1=resID, 10=chainID (numeric rank, as assigned by readpdb_strict)
# - Distance l_ij = minimal number of steps along this path from i to j
#                    (index difference q - p, with q > p)
# - Closeness (for all residues present in both .cor and path):
#       O_i = (N-1) / sum_{j != i} l_ij
#   where N is the number of residues included.
# - Outputs:
#   * shortest_paths_all_pairs_chain1 : rows [start end distance nodes...]
#   * closeness_float                  : "  xx.x\t0.xxxxxx"
#   * closeness_peaks                  : peak residues from closeness profile
#   * closeness_chain_labels.txt       : "'12A' , 0.096112" style (using PDB chain IDs)
#   * closeness_chain_plot.png         : profile + peaks

import os
import re
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

# --------------------------- configuration ---------------------------
IN_PATH_BASE = "path_file"      # tab-separated; first row = long path
IN_COOR_BASE = "coor_file"      # numeric; col1=resID, col10=chainID (numeric rank)

PATHS_FILE      = "shortest_paths_all_pairs_chain1"
CLOSENESS_FILE  = "closeness_float"
PEAKS_FILE      = "closeness_peaks"            # peak residues
PLOT_FILE       = "closeness_chain_plot.png"   # plot for all residues
LABELS_FILE     = "closeness_chain_labels.txt" # '12A' , 0.096112 style

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


def _find_latest_pdb(directory: str = ".") -> str | None:
    """
    Find the most recently modified .pdb file in `directory`.
    Returns its full path, or None if none found.
    """
    pdb_files = [f for f in os.listdir(directory) if f.lower().endswith(".pdb")]
    if not pdb_files:
        return None
    pdb_files_full = [os.path.join(directory, f) for f in pdb_files]
    pdb_files_full.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return pdb_files_full[0]


def _build_chain_mapping_from_pdb(coor: np.ndarray, directory: str = ".") -> dict[int, str]:
    """
    Build mapping from numeric chain ID (as stored in .cor/.coor) to actual
    chain letter in the PDB, using the order of first appearance in the PDB.

    Steps:
    - Find latest .pdb in `directory`.
    - Read lines starting with ATOM/HETATM.
    - Extract chainID at column 22 (PDB 1-based indexing).
    - Take unique chain IDs in order of first appearance → list of letters.
    - Numeric chain IDs in coor (1..N) are mapped to this list:
          1 -> first letter, 2 -> second, ...
    If anything fails, fall back to mapping k -> str(k).
    """
    # numeric chain IDs present in coor
    chain_nums = sorted(set(coor[:, 9].astype(int)))

    pdb_path = _find_latest_pdb(directory)
    if pdb_path is None:
        # Safe fallback: no PDB → use numeric IDs as strings
        return {k: str(k) for k in chain_nums}

    chain_letters_ordered: list[str] = []
    seen = set()

    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                if len(line) < 22:
                    continue
                ch = line[21].strip()  # column 22 in 1-based indexing
                if not ch:
                    continue
                if ch not in seen:
                    seen.add(ch)
                    chain_letters_ordered.append(ch)
    except Exception:
        # If reading/parsing PDB fails, fall back to numeric-as-string
        return {k: str(k) for k in chain_nums}

    if not chain_letters_ordered:
        # No chain letters parsed → fallback
        return {k: str(k) for k in chain_nums}

    # Now build the mapping: numeric rank (as used in .cor) → actual chain letter
    mapping: dict[int, str] = {}
    for k in chain_nums:
        # ranks are 1-based; if out of range, fallback to numeric string
        idx = k - 1
        if 0 <= idx < len(chain_letters_ordered):
            mapping[k] = chain_letters_ordered[idx]
        else:
            mapping[k] = str(k)

    return mapping


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
    Load coor_file (only for residue + chain info; not used in distance).
    """
    coor = np.loadtxt(coorfile)
    if coor.ndim != 2 or coor.shape[1] < 10:
        raise ValueError(
            "coor_file must have at least 10 columns (got shape %s)" % (coor.shape,)
        )
    return coor


def encode_node_float(resid: float, chain_num: float) -> float:
    """
    MATLAB-style encoding: v = resid + chain/10  (float).
    """
    return resid + chain_num / 10.0


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


def main():
    # --- pick latest files ---
    in_pathfile = _pick_latest(IN_PATH_BASE)
    in_coorfile = _pick_latest(IN_COOR_BASE)
    work_dir = os.path.dirname(os.path.abspath(in_coorfile)) or "."

    # print("Using files:")
    # print(f"  path_file  -> {in_pathfile}")
    # print(f"  coor_file  -> {in_coorfile}")

    # ------------------------ Load data ------------------------
    seq = _load_path_first_row(in_pathfile)  # long path (float encoded as resid + chain/10)
    n = seq.size
    assert n >= 2, "Path too short."

    coor = _load_coor_matrix(in_coorfile)

    # Build mapping numeric_chainID -> actual chain letter from PDB
    chain_map = _build_chain_mapping_from_pdb(coor, directory=work_dir)

    # From coor_file: residue & chain (numeric)
    resID_cor  = coor[:, 0].astype(int)
    chain_cor  = coor[:, 9].astype(int)  # column 10 (0-based index 9)

    # Build integer "node key" for residues in .cor:
    # key = 10*resID + chainNum  (since encoding is resid + chain/10)
    node_keys_cor = set(int(r * 10 + c) for r, c in zip(resID_cor, chain_cor))

    # Convert path sequence to integer keys the same way:
    # seq (float) = resid + chain/10 → key = round(10*seq)
    seq_keys = np.round(seq * 10.0).astype(int)

    # Build position list for all nodes appearing in the path
    pos_dict_all = defaultdict(list)
    for idx, k in enumerate(seq_keys):
        pos_dict_all[k].append(idx)
    for k in pos_dict_all:
        arr = np.asarray(pos_dict_all[k], dtype=int)
        arr.sort()
        pos_dict_all[k] = arr

    # Restrict to residues that are both in .cor and in the path
    common_keys = [k for k in node_keys_cor if k in pos_dict_all]
    if not common_keys:
        raise RuntimeError("No residues found that are common to .cor and path_file.")

    # Sort residues grouped by chain ID first, then by residue ID
    # key = 10*resID + chainNum → resID = key // 10, chainNum = key % 10
    common_keys_sorted = sorted(common_keys, key=lambda k: (k % 10, k // 10))
    node_keys = np.array(common_keys_sorted, dtype=int)

    # Decode back to resID and chainNum from key = 10*resID + chainNum
    resIDs = (node_keys // 10).astype(int)
    chainNums = (node_keys % 10).astype(int)

    # position lists for each node (in the new chain-grouped order)
    pos_list = [pos_dict_all[k] for k in node_keys]

    N = len(node_keys)

    # --- Distance matrix L: minimal index difference along the path (q - p, q>p) ---
    # Also store best start/end indices in seq for reconstructing the actual subpath.
    L = np.full((N, N), np.inf, dtype=float)
    start_idx = np.full((N, N), -1, dtype=int)
    end_idx   = np.full((N, N), -1, dtype=int)

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
    # start/end are encoded float nodes (resID + chain/10), segments from original seq
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
        start_enc = node_keys[a] / 10.0  # float encoded
        for b in range(N):
            if a == b:
                continue
            row += 1
            end_enc = node_keys[b] / 10.0
            M[row-1, 0] = start_enc
            M[row-1, 1] = end_enc

            ib = start_idx[a, b]
            jb = end_idx[a, b]
            if ib < 0 or jb < 0:
                M[row-1, 2] = np.inf
            else:
                seg = seq[ib:jb+1]
                M[row-1, 2] = L[a, b]
                M[row-1, 3:3+seg.size] = seg

    # Sort rows by start then end (encodings)
    order_paths = np.lexsort((M[:, 1], M[:, 0]))
    M = M[order_paths, :]

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
    finite_mask = np.isfinite(L)
    np.fill_diagonal(finite_mask, False)   # exclude self-distances
    L_for_sum = np.where(finite_mask, L, 0.0)
    sum_dists = np.sum(L_for_sum, axis=1)  # shape (N,)

    closenessVals = np.zeros(N, dtype=float)
    nonzero = sum_dists > 0
    closenessVals[nonzero] = (N - 1) / sum_dists[nonzero]

    # Encoded nodes v = resID + chain/10
    residues_encoded = node_keys.astype(float) / 10.0

    # --- Save CLOSENESS_FILE: "  xx.x\t0.xxxxxx" ---
    with open(CLOSENESS_FILE, "w") as fid:
        for i in range(N):
            fid.write(f" {residues_encoded[i]:4.1f}\t{closenessVals[i]:.6f}\n")

    # ----------------- Peak residues via peakdet (for all residues) -----------------
    res_int_all   = resIDs.copy()
    chain_int_all = chainNums.copy()
    closeness_all = closenessVals.copy()

    peaks_idx = []
    if closeness_all.size > 0:
        v = closeness_all
        std_v = float(np.std(v))
        if std_v > 0:
            try:
                maxtab, _ = peakdet(v, std_v)
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
                f"{res_int_all[idx]:4d}\t"
                f"{chain_int_all[idx]:d}\t"
                f"{residues_encoded[idx]:4.1f}\t"
                f"{closeness_all[idx]:.6f}\n"
            )

    # ----------------- Plot & label file for closeness -----------------
    if res_int_all.size > 0:
        # Human-readable labels (residue + real chain letter)
        labels = []
        for r, c in zip(res_int_all, chain_int_all):
            ch = chain_map.get(c, str(c))
            labels.append(f"{r}{ch}")

        # Text file: '12A' , 0.096112
        with open(LABELS_FILE, "w") as lf:
            for lab, val in zip(labels, closeness_all):
                lf.write(f"'{lab}' , {val:.6f}\n")

        # PLOT: all residues + peaks (x-axis grouped by chain)
        x = np.arange(len(res_int_all))

        plt.figure(figsize=(14, 4))
        plt.plot(x, closeness_all, marker="o", linestyle="-", linewidth=1.0, markersize=2.5, label="Closeness")

        if peaks_idx:
            peak_x = np.array(peaks_idx, dtype=int)
            peak_y = closeness_all[peak_x]
            plt.scatter(peak_x, peak_y, s=40, marker="s", label="Peaks")

        # --- Make x-axis readable: show at most ~80 labels ---
        Nres = len(labels)
        if Nres <= 80:
            tick_positions = x
        else:
            step = max(1, Nres // 80)
            tick_positions = x[::step]
        tick_labels = [labels[i] for i in tick_positions]

        plt.xticks(tick_positions, tick_labels, rotation=90, fontsize=6)

        # --- Visual chain grouping: vertical lines and chain letters on top ---
        unique_chains, first_idx, counts = np.unique(chain_int_all, return_index=True, return_counts=True)
        ymax = float(np.max(closeness_all)) if closeness_all.size > 0 else 1.0
        for i in range(1, len(first_idx)):
            # vertical separator between chains
            plt.axvline(first_idx[i] - 0.5, linestyle="--", linewidth=0.5, alpha=0.5)

        for ch_num, start, count in zip(unique_chains, first_idx, counts):
            center = start + (count - 1) / 2.0
            ch_letter = chain_map.get(int(ch_num), str(int(ch_num)))
            plt.text(center, ymax * 1.02, ch_letter,
                     ha="center", va="bottom", fontsize=8)

        plt.xlabel("Residue (grouped by chain: number + chain)")
        plt.ylabel("Closeness centrality")
        plt.title("Closeness centrality along all residues (step-count metric)")

        if peaks_idx:
            plt.legend()
        plt.tight_layout()
        plt.savefig(PLOT_FILE, dpi=300, bbox_inches="tight")
        plt.close()

        # print(f"Saved closeness plot to '{PLOT_FILE}' and labels to '{LABELS_FILE}'.")
    else:
        print("No residues for plotting / label export.")

    print(f"Wrote {PATHS_FILE}, {CLOSENESS_FILE}, {PEAKS_FILE} (peaks), and {LABELS_FILE}.")

    if peaks_idx:
        print("Peak residues (resID + chain):")
        for idx in peaks_idx:
            ch = chain_map.get(chain_int_all[idx], str(chain_int_all[idx]))
            print(
                f"  {res_int_all[idx]}{ch}  -> closeness = {closeness_all[idx]:.6f}"
            )
    else:
        print("No peaks found with delta = std(closeness).")


if __name__ == "__main__":
    main()
