# betweenness.py
# --------------------------
# Computes betweenness centrality on a chain using the
# shortest paths stored in shortest_paths_all_pairs_chain1.
#
# For each unordered pair (j,k) with a finite shortest path:
#   - Denominator (global): total_paths = number of such shortest paths.
#   - Numerator for residue i: count how many of those paths have i as an
#     INTERNAL node (not start or end). If residue i appears multiple times
#     in the same path, it still counts only once for that path.
#
# betweenness(i) = ( #paths where i is internal ) / (total_paths)
#
# Paths and nodes are encoded as v = resID + chainID/10.
#
# Outputs (analogous to closeness.py):
#   * betweenness_float              : "  xx.x\t0.xxxxxx"
#   * betweenness_peaks              : peak residues
#   * betweenness_chain_labels.txt   : "'1A' , 0.096112" style
#   * betweenness_chain_plot.png     : profile + peaks

import os
import re
import numpy as np
import matplotlib.pyplot as plt

# ------------- configuration -------------
PATHS_BASE       = "shortest_paths_all_pairs_chain1"
BETWEENNESS_FILE = "betweenness_float"
PEAKS_FILE       = "betweenness_peaks"
LABELS_FILE      = "betweenness_chain_labels.txt"
PLOT_FILE        = "betweenness_chain_plot.png"

TOL = 1e-9

# default – can be overridden via main(chain_id=...)
CHAIN_ID_SELECTED = 1
# ----------------------------------------


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


def _chain_to_int(ch):
    """
    Convert chain ID from UI to a numeric index.
    'A' -> 1, 'B' -> 2, ..., 'Z' -> 26
    '1','2',... stay numeric.
    """
    if ch is None:
        return None
    if isinstance(ch, str):
        ch = ch.strip()
        if not ch:
            return None
        if ch.isalpha() and len(ch) == 1:
            return ord(ch.upper()) - ord("A") + 1
        if ch.isdigit():
            return int(ch)
    return int(ch)


def _decode_node(v: float):
    """
    Decode MATLAB-style node encoding v = resid + chainID/10
    Returns (resid:int, chainID:int)
    """
    resid = int(np.floor(v + 1e-6))
    chain = int(round((v - resid) * 10))
    return resid, chain


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


def main(chain_id=None):
    """
    Compute betweenness on the chain specified by chain_id (e.g. 'A' or '1'),
    using shortest_paths_all_pairs_chain1[_N] in the current directory.
    """
    global CHAIN_ID_SELECTED

    # --- handle chain ID coming from Colab UI ---
    cid = _chain_to_int(chain_id) if chain_id is not None else None
    if cid is not None:
        CHAIN_ID_SELECTED = cid

    # --- open latest shortest_paths file ---
    paths_file = _pick_latest(PATHS_BASE)
    print(f"Using shortest-paths file: {paths_file}")
    print(f"Betweenness will be computed for chain ID = {CHAIN_ID_SELECTED}")

    # betweenness counts: encoded node (res + chain/10) -> count of paths where it's internal
    betw_counts = {}
    # all residues on this chain (start, end or internal) so they are in output even if 0
    residues_on_chain = set()

    # To avoid double counting j→k and k→j, use unordered pair keys
    seen_pairs = set()

    # Global denominator: total number of shortest paths (unordered j,k) with finite distance
    total_paths = 0

    # --- read file line by line ---
    with open(paths_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue

            # first two columns: start and end encoded nodes
            try:
                s = float(parts[0])
                t = float(parts[1])
            except ValueError:
                continue

            dist_str = parts[2].strip()
            if dist_str == "Inf":
                # no valid shortest path; skip
                continue

            # columns 4..end: path nodes (encoded) – may be empty
            nodes = []
            for p in parts[3:]:
                p = p.strip()
                if not p:
                    continue
                try:
                    nodes.append(float(p))
                except ValueError:
                    pass
            if not nodes:
                # something odd; skip
                continue

            # record start/end residues on this chain (betweenness can be 0 but they should appear)
            for v in (s, t):
                resid, chain = _decode_node(v)
                if chain == CHAIN_ID_SELECTED:
                    residues_on_chain.add(v)
                    betw_counts.setdefault(v, 0.0)

            # unordered pair key
            u = min(s, t)
            v2 = max(s, t)
            pair_key = (u, v2)
            if pair_key in seen_pairs:
                # we've already processed this j–k pair
                continue
            seen_pairs.add(pair_key)

            # a valid shortest path for this unordered pair exists → update denominator
            total_paths += 1

            # internal nodes: exclude first and last
            if len(nodes) <= 2:
                continue  # no internal residues for this path

            internal_nodes = nodes[1:-1]

            # each residue should count at most once per path, even if it appears multiple times
            internal_unique = set(internal_nodes)

            for enc in internal_unique:
                resid, chain = _decode_node(enc)
                if chain != CHAIN_ID_SELECTED:
                    continue  # only count internal residues on the selected chain
                residues_on_chain.add(enc)
                betw_counts.setdefault(enc, 0.0)
                betw_counts[enc] += 1.0

    if not residues_on_chain:
        print("No residues found on selected chain in shortest-paths file; nothing to do.")
        return

    if total_paths == 0:
        print("No finite shortest paths found; betweenness will be zero for all residues.")
        total_paths = 1  # avoid division by zero; all numerators are zero anyway

    # Sort residues by encoded value (resID + chain/10) – OK for single chain
    residues_sorted = sorted(residues_on_chain)
    numer_vals = np.array([betw_counts.get(r, 0.0) for r in residues_sorted], dtype=float)
    residues_encoded = np.array(residues_sorted, dtype=float)

    # Normalize by total_paths
    betw_vals = numer_vals / float(total_paths)

    # --- write betweenness file in same style as closeness_float ---
    with open(BETWEENNESS_FILE, "w") as fid:
        for enc, val in zip(residues_encoded, betw_vals):
            fid.write(f" {enc:4.1f}\t{val:.6f}\n")

    print(f"Wrote betweenness centrality to '{BETWEENNESS_FILE}'.")
    print(f"Total number of shortest paths (denominator) = {total_paths}")

    # ----------------- Peak detection & reporting -----------------
    peaks_idx = []
    if betw_vals.size > 0:
        std_v = float(np.std(betw_vals))
        if std_v > 0:
            try:
                maxtab, _ = peakdet(betw_vals, std_v)
            except ValueError as e:
                print("Peak detection warning:", e)
                maxtab = np.empty((0, 2))
            if maxtab.size > 0:
                peaks_idx = maxtab[:, 0].astype(int).tolist()
        else:
            maxtab = np.empty((0, 2))

    # Determine chain letter (for labels like "1A")
    chain_num = CHAIN_ID_SELECTED
    if 1 <= chain_num <= 26:
        chain_letter = chr(ord("A") + chain_num - 1)
    else:
        chain_letter = str(chain_num)

    # Decode residue numbers (integers)
    res_ids = np.array([_decode_node(enc)[0] for enc in residues_encoded], dtype=int)
    labels = [f"{rid}{chain_letter}" for rid in res_ids]

    # --- labels file: '1A' , 0.096112
    with open(LABELS_FILE, "w") as lf:
        for lab, val in zip(labels, betw_vals):
            lf.write(f"'{lab}' , {val:.6f}\n")

    # --- peaks file (detailed) ---
    with open(PEAKS_FILE, "w") as pf:
        pf.write("# resID\tchainID\tencoded_node\tbetweenness_peak\n")
        for idx in peaks_idx:
            pf.write(
                f"{res_ids[idx]:4d}\t"
                f"{chain_num:d}\t"
                f"{residues_encoded[idx]:4.1f}\t"
                f"{betw_vals[idx]:.6f}\n"
            )

    # ----------------- Plot -----------------
    x = np.arange(len(res_ids))

    plt.figure(figsize=(12, 4))
    plt.plot(x, betw_vals, marker="o", linestyle="-", label="Betweenness")

    if peaks_idx:
        peak_x = np.array(peaks_idx, dtype=int)
        peak_y = betw_vals[peak_x]
        plt.scatter(peak_x, peak_y, s=50, marker="s", label="Peaks")

    plt.xticks(x, labels, rotation=90, fontsize=6)
    plt.xlabel("Residue (number + chain)")
    plt.ylabel("Betweenness centrality")
    plt.title(f"Betweenness centrality along chain {chain_letter}")
    if peaks_idx:
        plt.legend()
    plt.tight_layout()
    plt.savefig(PLOT_FILE, dpi=300)
    plt.close()

    print(f"Saved betweenness plot to '{PLOT_FILE}' and labels to '{LABELS_FILE}'.")
    print(f"Peaks written to '{PEAKS_FILE}'.")

    # Small text summary
    print("Summary (first few residues):")
    for enc, num, val in list(zip(residues_encoded, numer_vals, betw_vals))[:10]:
        resid, ch = _decode_node(enc)
        print(
            f"  residue {resid} chain {ch} (encoded {enc:4.1f}) "
            f"-> paths_through_i = {num:.0f}, betweenness = {val:.6f}"
        )

    if peaks_idx:
        print("Peak residues (resID on chain", chain_letter, "):")
        for idx in peaks_idx:
            print(f"  {res_ids[idx]}  -> betweenness = {betw_vals[idx]:.6f}")
    else:
        print("No peaks found on the selected chain with delta = std(betweenness).")


if __name__ == "__main__":
    main()
