# betweenness.py
# --------------------------
# Computes betweenness centrality on a chain using the
# shortest paths stored in shortest_paths_all_pairs_chain1.
#
# Formula (paper):
#   b_i = sum_{jk} g_jik / g_jk
# where:
#   g_jk   = number of shortest paths between j and k
#   g_jik  = number of those paths that pass through i
#
# In our current pipeline, for each pair (j,k) we keep
# exactly ONE shortest path, so g_jk = 1.
# Thus, for each unordered pair (j,k), if residue i is an
# internal node on that path, it contributes +1 to b_i.

import os
import re
import numpy as np

# ------------- configuration -------------
PATHS_BASE       = "shortest_paths_all_pairs_chain1"
BETWEENNESS_FILE = "betweenness_float"
TOL = 1e-9

# default – will be overridden if main(chain_id=...) is used
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


def main(chain_id=None):
    global CHAIN_ID_SELECTED

    # --- handle chain ID coming from Colab UI ---
    cid = _chain_to_int(chain_id) if chain_id is not None else None
    if cid is not None:
        CHAIN_ID_SELECTED = cid

    # --- open latest shortest_paths file ---
    paths_file = _pick_latest(PATHS_BASE)
    print(f"Using shortest-paths file: {paths_file}")
    print(f"Betweenness will be computed for chain ID = {CHAIN_ID_SELECTED}")

    # betweenness counts: encoded node (res + chain/10) -> count
    betw_counts = {}
    # also track all residues on this chain (even if betw=0) so they appear in output
    residues_on_chain = set()

    # To avoid double counting j→k and k→j, use unordered pair keys
    seen_pairs = set()

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
                # no valid path; skip
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

            # ensure start/end also recorded as residues on this chain (even if no betweenness)
            for v in (s, t):
                resid, chain = _decode_node(v)
                if chain == CHAIN_ID_SELECTED:
                    residues_on_chain.add(v)
                    betw_counts.setdefault(v, 0.0)

            # unordered pair key
            u = min(s, t)
            v = max(s, t)
            pair_key = (u, v)
            if pair_key in seen_pairs:
                # we've already processed this unordered pair; skip the reverse row
                continue
            seen_pairs.add(pair_key)

            # internal nodes: exclude first and last node in the path
            if len(nodes) <= 2:
                continue  # no internal residues

            internal_nodes = nodes[1:-1]

            # Update betweenness: for each internal residue on selected chain, +1
            for enc in internal_nodes:
                resid, chain = _decode_node(enc)
                if chain != CHAIN_ID_SELECTED:
                    continue
                residues_on_chain.add(enc)
                betw_counts.setdefault(enc, 0.0)
                betw_counts[enc] += 1.0

    if not residues_on_chain:
        print("No residues found on selected chain in shortest-paths file; nothing to do.")
        return

    # Sort residues by encoded value (resID + chain/10) – OK for single chain
    residues_sorted = sorted(residues_on_chain)
    values = [betw_counts.get(r, 0.0) for r in residues_sorted]

    # Optional: normalization (commented; you can enable if needed)
    # N = len(residues_sorted)
    # max_pairs = (N - 1) * (N - 2) / 2.0 if N > 2 else 1.0
    # values = [v / max_pairs for v in values]

    # --- write betweenness file in same style as closeness_float ---
    with open(BETWEENNESS_FILE, "w") as fid:
        for enc, val in zip(residues_sorted, values):
            fid.write(f" {enc:4.1f}\t{val:.6f}\n")

    print(f"Wrote betweenness centrality to '{BETWEENNESS_FILE}'.")
    print("Summary (first few residues):")
    for enc, val in list(zip(residues_sorted, values))[:10]:
        resid, ch = _decode_node(enc)
        print(f"  residue {resid} chain {ch} (encoded {enc:4.1f}) -> betweenness {val:.6f}")


if __name__ == "__main__":
    main()
