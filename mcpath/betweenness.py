# betweenness.py
# --------------------------
# Computes betweenness centrality for all residues using the
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
# Paths and nodes are encoded as v = resID + chainID/10,
# where chainID is the numeric rank assigned by readpdb_strict.
#
# Outputs (analogous to closeness.py):
#   * betweenness_float              : "  xx.x\t0.xxxxxx"
#   * betweenness_peaks              : peak residues
#   * betweenness_chain_labels.txt   : "'12A' , 0.096112" style
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


def _find_latest_pdb(directory: str = ".") -> str | None:
    """
    Find the most recently modified .pdb file in `directory`.
    Returns full path or None if not found.
    """
    pdb_files = [f for f in os.listdir(directory) if f.lower().endswith(".pdb")]
    if not pdb_files:
        return None
    pdb_files_full = [os.path.join(directory, f) for f in pdb_files]
    pdb_files_full.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return pdb_files_full[0]


def _build_chain_mapping_from_pdb(chain_nums: list[int], directory: str = ".") -> dict[int, str]:
    """
    Build mapping from numeric chain ID (rank used in encoding) to actual PDB chain letter.

    - Find latest .pdb in `directory`.
    - Read ATOM/HETATM records and extract chainID at column 22.
    - Collect unique chain IDs in order of first appearance → list of letters.
    - Map:
          1 -> first letter, 2 -> second, ...
      for any chain_nums in the provided list.

    If anything fails or PDB is missing, fall back to k -> str(k).
    """
    chain_nums_sorted = sorted(set(chain_nums))
    pdb_path = _find_latest_pdb(directory)
    if pdb_path is None:
        return {k: str(k) for k in chain_nums_sorted}

    chain_letters_ordered: list[str] = []
    seen = set()

    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                if len(line) < 22:
                    continue
                ch = line[21].strip()  # column 22 (1-based)
                if not ch:
                    continue
                if ch not in seen:
                    seen.add(ch)
                    chain_letters_ordered.append(ch)
    except Exception:
        # On any error, fall back
        return {k: str(k) for k in chain_nums_sorted}

    if not chain_letters_ordered:
        return {k: str(k) for k in chain_nums_sorted}

    mapping: dict[int, str] = {}
    for k in chain_nums_sorted:
        idx = k - 1  # ranks are 1-based
        if 0 <= idx < len(chain_letters_ordered):
            mapping[k] = chain_letters_ordered[idx]
        else:
            mapping[k] = str(k)

    return mapping


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


def main():
    """
    Compute betweenness for all residues present in the shortest-paths file,
    without requiring a chain ID. Chain letters in labels are taken from
    the PDB via first-appearance order (consistent with readpdb_strict).
    """
    # --- open latest shortest_paths file ---
    paths_file = _pick_latest(PATHS_BASE)
    work_dir = os.path.dirname(os.path.abspath(paths_file)) or "."
    # print(f"Using shortest-paths file: {paths_file}")

    # betweenness counts: encoded node (res + chain/10) -> count of paths where it's internal
    betw_counts: dict[float, float] = {}
    # all residues (start, end or internal) so they are in output even if 0
    residues_encoded_set = set()
    # track which chain numeric IDs appear at all
    chain_nums_set = set()

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

            # decode start/end to track chain IDs and residues
            for v in (s, t):
                resid, chain = _decode_node(v)
                residues_encoded_set.add(v)
                chain_nums_set.add(chain)
                betw_counts.setdefault(v, 0.0)

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
                residues_encoded_set.add(enc)
                chain_nums_set.add(chain)
                betw_counts.setdefault(enc, 0.0)
                betw_counts[enc] += 1.0

    if not residues_encoded_set:
        print("No residues found in shortest-paths file; nothing to do.")
        return

    if total_paths == 0:
        print("No finite shortest paths found; betweenness will be zero for all residues.")
        total_paths = 1  # avoid division by zero; all numerators are zero anyway

    # --- sort residues grouped by chainID (numeric rank) and then residue ID ---
    residues_sorted = sorted(
        residues_encoded_set,
        key=lambda v: (_decode_node(v)[1], _decode_node(v)[0])  # (chainID, resid)
    )
    residues_encoded = np.array(residues_sorted, dtype=float)
    numer_vals = np.array([betw_counts.get(r, 0.0) for r in residues_encoded], dtype=float)

    # Decode residue numbers and chain numeric ranks (in the new grouped order)
    res_ids = np.array([_decode_node(enc)[0] for enc in residues_encoded], dtype=int)
    chain_nums = np.array([_decode_node(enc)[1] for enc in residues_encoded], dtype=int)

    # Build mapping numeric chainID -> PDB chain letter from PDB
    chain_map = _build_chain_mapping_from_pdb(list(chain_nums_set), directory=work_dir)

    # Normalize by total_paths
    betw_vals = numer_vals / float(total_paths)

    # --- write betweenness file in same style as closeness_float ---
    with open(BETWEENNESS_FILE, "w") as fid:
        for enc, val in zip(residues_encoded, betw_vals):
            fid.write(f" {enc:4.1f}\t{val:.6f}\n")

    # print(f"Wrote betweenness centrality to '{BETWEENNESS_FILE}'.")
    # print(f"Total number of shortest paths (denominator) = {total_paths}")

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

    # Labels with real chain letters from PDB (same grouped order)
    labels = []
    for rid, cnum in zip(res_ids, chain_nums):
        ch = chain_map.get(cnum, str(cnum))
        labels.append(f"{rid}{ch}")

    # --- labels file: '12A' , 0.096112
    with open(LABELS_FILE, "w") as lf:
        for lab, val in zip(labels, betw_vals):
            lf.write(f"'{lab}' , {val:.6f}\n")

    # --- peaks file (detailed) ---
    with open(PEAKS_FILE, "w") as pf:
        pf.write("# resID\tchainID_numeric\tencoded_node\tbetweenness_peak\n")
        for idx in peaks_idx:
            pf.write(
                f"{res_ids[idx]:4d}\t"
                f"{chain_nums[idx]:d}\t"
                f"{residues_encoded[idx]:4.1f}\t"
                f"{betw_vals[idx]:.6f}\n"
            )

    # ----------------- Plot -----------------
    x = np.arange(len(res_ids))

    plt.figure(figsize=(14, 4))
    plt.plot(x, betw_vals, marker="o", linestyle="-",
             linewidth=1.0, markersize=2.5, label="Betweenness")

    if peaks_idx:
        peak_x = np.array(peaks_idx, dtype=int)
        peak_y = betw_vals[peak_x]
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
    if betw_vals.size > 0:
        ymax = float(np.max(betw_vals))
    else:
        ymax = 1.0

    unique_chains, first_idx, counts = np.unique(chain_nums, return_index=True, return_counts=True)

    # vertical separators between chains
    for i in range(1, len(first_idx)):
        plt.axvline(first_idx[i] - 0.5, linestyle="--", linewidth=0.5, alpha=0.5)

    # chain letters above each chain block
    for ch_num, start, count in zip(unique_chains, first_idx, counts):
        center = start + (count - 1) / 2.0
        ch_letter = chain_map.get(int(ch_num), str(int(ch_num)))
        plt.text(center, ymax * 1.02, ch_letter,
                 ha="center", va="bottom", fontsize=8)

    plt.ylim(0, ymax * 1.10 if ymax > 0 else 1.0)

    plt.xlabel("Residue (grouped by chain: number + chain)")
    plt.ylabel("Betweenness centrality")
    plt.title("Betweenness centrality along all residues")

    if peaks_idx:
        plt.legend()
    plt.tight_layout()
    plt.savefig(PLOT_FILE, dpi=300, bbox_inches="tight")
    plt.close()

    # print(f"Saved betweenness plot to '{PLOT_FILE}' and labels to '{LABELS_FILE}'.")
    # print(f"Peaks written to '{PEAKS_FILE}'.")

    # Small text summary
    # print("Summary (first few residues):")
    for enc, num, val, rid, cnum in list(zip(residues_encoded, numer_vals, betw_vals, res_ids, chain_nums))[:10]:
        ch = chain_map.get(cnum, str(cnum))
    #    print(
    #        f"  residue {rid}{ch} (encoded {enc:4.1f}) "
    #        f"-> paths_through_i = {num:.0f}, betweenness = {val:.6f}"
    #    )

    if peaks_idx:
        print("Peak residues (resID + chain):")
        for idx in peaks_idx:
            ch = chain_map.get(chain_nums[idx], str(chain_nums[idx]))
            print(f"  {res_ids[idx]}{ch}  -> betweenness = {betw_vals[idx]:.6f}")
    else:
        print("No peaks found with delta = std(betweenness).")


if __name__ == "__main__":
    main()
