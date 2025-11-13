# closeness.py
# ---------------------------------------------------------------
# 1) Reproduces the MATLAB "closeness" function for
#    shortest_paths_all_pairs_chain1:
#       - Uses Euclidean CA distances along the long path.
#       - Same sweep / tie-breaking as your MATLAB code.
#       - Writes paths file with Euclidean distance in column 3.
# 2) Then builds a *step-count* distance matrix:
#       d_ij = (#residues on Euclidean-shortest path between i and j) - 1
#    using the minimum of i->j and j->i (symmetric).
# 3) Computes closeness centrality from these step-count distances:
#       O_i = (N-1) / sum_{j!=i} d_ij
#    and writes "closeness_float" as in MATLAB.

import os
import re
import numpy as np

# --------------------------- configuration ---------------------------
IN_PATH_BASE = "path_file"      # tab-separated; first row = long path
IN_COOR_BASE = "coor_file"      # numeric; col1=resID, cols6-8=CA xyz, col10=chainID
IN_ATOM_BASE = "atom_file"      # optional, just reported

PATHS_FILE      = "shortest_paths_all_pairs_chain1"
CLOSENESS_FILE  = "closeness_float"

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
    """Load first row of path_file as a vector seq, strip zeros (MATLAB style)."""
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
        raise ValueError(
            "coor_file must have at least 10 columns (got shape %s)" % (coor.shape,)
        )
    return coor


def encode_node(resid: float, chainid: float) -> float:
    # MATLAB: encode = @(rid,cid) rid + cid/10;
    return resid + chainid / 10.0


def main():
    # --- pick latest files ---
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

    # --- Load data (as in MATLAB) ---
    seq = _load_path_first_row(in_pathfile)
    n = seq.size
    assert n >= 2, "Path too short."

    coor = _load_coor_matrix(in_coorfile)
    resID_list   = coor[:, 0]
    chainID_list = coor[:, 9]
    CAxyz        = coor[:, 5:8]   # CA coordinates (x,y,z)

    # Decode nodes v = resID + chainID/10
    resID_seq   = np.floor(seq + 1e-6)
    chainID_seq = np.round((seq - resID_seq) * 10)

    # --- Residues on selected chain (SORTED) ---
    res_on_chain = np.unique(resID_seq[chainID_seq == CHAIN_ID_SELECTED])
    Nchain = res_on_chain.size
    assert Nchain > 0, (
        f"No residues from chain {CHAIN_ID_SELECTED} found in the selected path."
    )

    # --- Precompute step geometry & validity (CA distances) ---
    keyCoor = (resID_list.astype(np.int64) * 100) + chainID_list.astype(np.int64)
    coorIdxMap = {int(k): int(i) for i, k in enumerate(keyCoor)}

    dist = np.zeros(n - 1, dtype=float)
    invalid = np.zeros(n - 1, dtype=bool)

    for k in range(n - 1):
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

    S = np.zeros(n, dtype=float)
    C = np.zeros(n, dtype=int)
    if n > 1:
        S[1:] = np.cumsum(dist)
        C[1:] = np.cumsum(invalid.astype(int))
    maxC = int(np.max(C))

    # --- For each residue on selected chain, indices where it occurs in the path ---
    occIdx = {}
    for rid in res_on_chain:
        rid_val = float(rid)
        occ = np.where((resID_seq == rid_val) & (chainID_seq == CHAIN_ID_SELECTED))[0]
        occIdx[rid_val] = occ

    # --- Fast shortest search for all ordered pairs (i != j) ---
    res_on_chain = res_on_chain.reshape(1, -1).ravel()
    N = res_on_chain.size

    bestI   = np.full((N, N), np.nan, dtype=float)
    bestJ   = np.full((N, N), np.nan, dtype=float)
    bestD   = np.full((N, N), np.inf, dtype=float)   # Euclidean shortest distance
    bestLen = np.zeros((N, N), dtype=int)            # #nodes of that path

    dBestForJ = np.full(n, np.inf, dtype=float)
    iBestForJ = np.full(n, np.nan, dtype=float)
    lenForJ   = np.zeros(n, dtype=int)

    bestS_forC = np.full(maxC + 1, -np.inf, dtype=float)
    bestI_forC = np.full(maxC + 1, np.nan, dtype=float)

    for ai in range(N):
        sRID = float(res_on_chain[ai])
        sIdx = occIdx[sRID]          # sorted indices where sRID occurs

        # Reset per-j results
        dBestForJ[:] = np.inf
        iBestForJ[:] = np.nan
        lenForJ[:]   = 0
        bestS_forC[:] = -np.inf
        bestI_forC[:] = np.nan

        # Single sweep across the whole path
        p = 0
        sIdx_len = sIdx.size

        for t in range(n):
            # If t is an eligible start index for sRID, update tables at its C(t)
            while p < sIdx_len and sIdx[p] == t:
                c = C[t]           # 0..maxC
                s_val = S[t]
                ib = bestI_forC[c]
                sb = bestS_forC[c]
                if (s_val > sb + TOL) or (
                    abs(s_val - sb) <= TOL and (np.isnan(ib) or t > ib)
                ):
                    bestS_forC[c] = s_val
                    bestI_forC[c] = float(t)
                p += 1

            # If t can be an end (on selected chain), try best compatible start
            if chainID_seq[t] == CHAIN_ID_SELECTED:
                c = C[t]
                ib = bestI_forC[c]
                if not np.isnan(ib):
                    ib_int = int(ib)
                    d   = S[t] - bestS_forC[c]
                    length = t - ib_int + 1
                    if (d + TOL < dBestForJ[t]) or (
                        abs(d - dBestForJ[t]) <= TOL and length < lenForJ[t]
                    ):
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

    # ---------------- Build output matrix M exactly like MATLAB ----------------
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
            row += 1
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
                M[row-1, 2] = bestD[a, b]                 # Euclidean distance
                M[row-1, 3:3+seg.size] = seg             # full path nodes

    # Sort rows by start then end (ascending)
    order = np.lexsort((M[:, 1], M[:, 0]))
    M = M[order, :]

    # --- Write PATHS_FILE with 6-decimal Euclidean distance in column 3 ---
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

    # ================== Step-count distance matrix for closeness ==================
    # For each residue index a,b (0..N-1), define
    #   d_ij = (#residues on Euclidean-shortest path) - 1
    #         = bestLen[a,b] - 1
    # Use min(d_ij, d_ji) to make it symmetric.
    distMat = np.full((N, N), np.inf, dtype=float)
    for a in range(N):
        distMat[a, a] = 0.0
        for b in range(a + 1, N):
            len_ab = bestLen[a, b]
            len_ba = bestLen[b, a]
            d = np.inf
            if len_ab > 0:
                d = min(d, float(len_ab - 1))
            if len_ba > 0:
                d = min(d, float(len_ba - 1))
            distMat[a, b] = distMat[b, a] = d
    # =====================================================================

    # --- Closeness centrality from step-count distances ---
    closenessVals = np.zeros(N, dtype=float)
    for i in range(N):
        rowVals = distMat[i, :]
        mask = np.isfinite(rowVals)
        mask[i] = False                     # exclude self (distance 0)
        finiteVals = rowVals[mask]
        if finiteVals.size > 0:
            s = float(np.sum(finiteVals))
            closenessVals[i] = (N - 1) / s
        else:
            closenessVals[i] = 0.0

    # Encoded residues for selected chain
    residues_encoded = np.array(
        [encode_node(float(rid), float(CHAIN_ID_SELECTED)) for rid in res_on_chain],
        dtype=float,
    )

    # --- Save CLOSENESS_FILE in MATLAB's "  xx.x\t0.xxxxxx" format ---
    with open(CLOSENESS_FILE, "w") as fid:
        for i in range(N):
            fid.write(f" {residues_encoded[i]:4.1f}\t{closenessVals[i]:.6f}\n")

    print(
        f"Wrote {PATHS_FILE} (Euclidean distances) and "
        f"{CLOSENESS_FILE} (step-count closeness with 6 decimals)."
    )
    print("Closeness min / max:", float(np.min(closenessVals)), float(np.max(closenessVals)))


if __name__ == "__main__":
    main()
