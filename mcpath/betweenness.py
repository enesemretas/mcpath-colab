# betweenness.py  (C-code compatible)
# ----------------------------------
# Matches the provided C code behavior as closely as possible using
# shortest_paths_all_pairs_chain1 + coor_file.
#
# C behavior:
# - Residue universe/order = coor_file order (N = res_num)
# - For each residue i:
#     bcount[i] = number of times coor[i] appears in my_path[k][l].shorty over ALL ordered pairs (k,l), k!=l
#     NOTE: counts endpoints too; counts multiple occurrences within same path.
# - betweenness[i] = bcount[i] / ((N-1)*N)
#
# Outputs:
#   * betweenness_float
#   * betweenness_peaks
#   * betweenness_chain_labels.txt
#   * betweenness_chain_plot.png

import os
import re
import numpy as np
import matplotlib.pyplot as plt

# ------------- configuration -------------
PATHS_BASE       = "shortest_paths_all_pairs_chain1"
IN_COOR_BASE     = "coor_file"

BETWEENNESS_FILE = "betweenness_float"
PEAKS_FILE       = "betweenness_peaks"
LABELS_FILE      = "betweenness_chain_labels.txt"
PLOT_FILE        = "betweenness_chain_plot.png"

# C uses E = 1e-6 for float comparisons
E_C = 1e-6
# ----------------------------------------


def _pick_latest(basename: str, directory: str = ".") -> str:
    """Pick basename, basename_2, basename_3, ... with highest suffix."""
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


def _find_latest_pdb(directory: str = ".") -> str | None:
    pdb_files = [f for f in os.listdir(directory) if f.lower().endswith(".pdb")]
    if not pdb_files:
        return None
    full = [os.path.join(directory, f) for f in pdb_files]
    full.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return full[0]


def _build_chain_mapping_from_pdb(coor: np.ndarray, directory: str = ".") -> dict[int, str]:
    """Map numeric chain ranks (coor col10) -> PDB chain letters by first appearance."""
    chain_nums = sorted(set(coor[:, 9].astype(int)))
    pdb_path = _find_latest_pdb(directory)
    if pdb_path is None:
        return {k: str(k) for k in chain_nums}

    chain_letters_ordered = []
    seen = set()
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                if len(line) < 22:
                    continue
                ch = line[21].strip()
                if not ch:
                    continue
                if ch not in seen:
                    seen.add(ch)
                    chain_letters_ordered.append(ch)
    except Exception:
        return {k: str(k) for k in chain_nums}

    if not chain_letters_ordered:
        return {k: str(k) for k in chain_nums}

    mapping = {}
    for k in chain_nums:
        idx = k - 1
        mapping[k] = chain_letters_ordered[idx] if 0 <= idx < len(chain_letters_ordered) else str(k)
    return mapping


def _load_coor_matrix(coorfile: str) -> np.ndarray:
    coor = np.loadtxt(coorfile, dtype=float)
    if coor.ndim != 2 or coor.shape[1] < 10:
        raise ValueError(f"coor_file must have at least 10 columns (got shape {coor.shape})")
    return coor


def _node_key_from_encoded_float(v: float) -> int:
    """
    v = resid + chain/10  -> key = round(10*v)
    Using rounding makes membership robust and consistent.
    """
    return int(np.round(v * 10.0))


# ------------------- MATLAB-like PEAKDET implementation -------------------
def peakdet(v, delta):
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

    maxtab_arr = np.asarray(maxtab, dtype=float) if maxtab else np.empty((0, 2))
    mintab_arr = np.asarray(mintab, dtype=float) if mintab else np.empty((0, 2))
    return maxtab_arr, mintab_arr
# -------------------------------------------------------------------------


def main():
    # --- pick latest files ---
    paths_file = _pick_latest(PATHS_BASE)
    coor_file  = _pick_latest(IN_COOR_BASE, directory=os.path.dirname(paths_file) or ".")
    work_dir = os.path.dirname(os.path.abspath(paths_file)) or "."

    # --- load coor: defines residue universe/order (C behavior) ---
    coor = _load_coor_matrix(coor_file)
    resID = coor[:, 0].astype(int)
    chainNum = coor[:, 9].astype(int)
    encoded = resID.astype(float) + chainNum.astype(float) / 10.0  # coor[i] in C
    N = int(encoded.size)

    # map encoded -> index in coor order
    keys = np.array([int(r * 10 + c) for r, c in zip(resID, chainNum)], dtype=int)
    key_to_i = {k: i for i, k in enumerate(keys)}

    # counts as in C: bcount[i]
    bcount = np.zeros(N, dtype=np.int64)

    # --- read shortest paths and count occurrences (C behavior) ---
    # C loops: for i in residues:
    #            for k in residues:
    #              for l in residues (l!=k):
    #                for m in path[k][l]:
    #                   if node==coor[i]: bcount[i]++
    #
    # Equivalent: for each path(k,l), for each node in nodes: bcount[node] += occurrences
    # (endpoints included, repeated occurrences count multiple times)

    with open(paths_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue

            dist_str = parts[2].strip()
            if dist_str == "Inf":
                # no path; C would still have some default length, but in practice
                # your paths_file should omit nodes for such cases. We skip.
                continue

            # nodes are columns 4..end in your file
            # (and may be empty for weird lines)
            if len(parts) < 4:
                continue

            # Parse node list; count each occurrence
            for token in parts[3:]:
                token = token.strip()
                if not token:
                    continue
                try:
                    v = float(token)
                except ValueError:
                    continue

                k = _node_key_from_encoded_float(v)
                idx = key_to_i.get(k, None)
                if idx is not None:
                    bcount[idx] += 1

    # --- normalize exactly like C ---
    denom = (N - 1) * N
    if denom <= 0:
        raise RuntimeError("Invalid N from coor_file; cannot normalize.")

    betw_vals = bcount.astype(float) / float(denom)

    # --- write betweenness_float (same style) ---
    with open(BETWEENNESS_FILE, "w") as fid:
        for v, val in zip(encoded, betw_vals):
            fid.write(f" {v:5.1f}\t{val:.6f}\n")

    # --- chain letter mapping for labels/plot ---
    chain_map = _build_chain_mapping_from_pdb(coor, directory=work_dir)
    labels = [f"{r}{chain_map.get(int(c), str(int(c)))}" for r, c in zip(resID, chainNum)]

    with open(LABELS_FILE, "w") as lf:
        for lab, val in zip(labels, betw_vals):
            lf.write(f"'{lab}' , {val:.6f}\n")

    # --- peaks ---
    peaks_idx = []
    if betw_vals.size > 0:
        std_v = float(np.std(betw_vals))
        if std_v > 0:
            maxtab, _ = peakdet(betw_vals, std_v)
            if maxtab.size > 0:
                peaks_idx = maxtab[:, 0].astype(int).tolist()

    with open(PEAKS_FILE, "w") as pf:
        pf.write("# resID\tchainID_numeric\tencoded_node\tbetweenness_peak\n")
        for idx in peaks_idx:
            pf.write(f"{resID[idx]:4d}\t{chainNum[idx]:d}\t{encoded[idx]:5.1f}\t{betw_vals[idx]:.6f}\n")

    # --- plot ---
    x = np.arange(N)
    plt.figure(figsize=(14, 4))
    plt.plot(x, betw_vals, marker="o", linestyle="-", linewidth=1.0, markersize=2.5, label="Betweenness (C-style)")

    if peaks_idx:
        px = np.array(peaks_idx, dtype=int)
        py = betw_vals[px]
        plt.scatter(px, py, s=40, marker="s", label="Peaks")

    # show at most ~80 x labels
    if N <= 80:
        tick_positions = x
    else:
        step = max(1, N // 80)
        tick_positions = x[::step]
    tick_labels = [labels[i] for i in tick_positions]
    plt.xticks(tick_positions, tick_labels, rotation=90, fontsize=6)

    ymax = float(np.max(betw_vals)) if betw_vals.size else 1.0
    plt.ylim(0, ymax * 1.10 if ymax > 0 else 1.0)

    plt.xlabel("Residue (coor_file order)")
    plt.ylabel("Betweenness (C normalization)")
    plt.title("Betweenness centrality (C-style counting: endpoints included, repeats counted)")

    if peaks_idx:
        plt.legend()

    plt.tight_layout()
    plt.savefig(PLOT_FILE, dpi=300, bbox_inches="tight")
    plt.close()

    if peaks_idx:
        print("Peak residues (resID + chain):")
        for idx in peaks_idx:
            ch = chain_map.get(int(chainNum[idx]), str(int(chainNum[idx])))
            print(f"  {resID[idx]}{ch}  -> betweenness = {betw_vals[idx]:.6f}")
    else:
        print("No peaks found with delta = std(betweenness).")


if __name__ == "__main__":
    main()
