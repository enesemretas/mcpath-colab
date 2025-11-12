# closeness.py
# Efficient all-pairs shortest contiguous paths + closeness centrality
# Uses most recent coor_file[_N], atom_file[_N], path_file[_N] in current directory.

import os
import re
import math
import numpy as np

# --------------------------- configuration ---------------------------
IN_PATH_BASE = "path_file"                  # tab-separated; first row = long path
IN_COOR_BASE = "coor_file"                  # numeric; col1=resID, cols6-8=CA xyz, col10=chainID
IN_ATOM_BASE = "atom_file"                  # optional, only reported

PATHS_OUT     = "shortest_paths_all_pairs_chain1"
CLOSENESS_OUT = "closeness_float"

CHAIN_ID_SELECTED = 1                       # which chain to analyze
TOL = 1e-9                                   # tie tolerance
# --------------------------------------------------------------------


def _pick_latest(basename: str, directory: str = ".") -> str:
    """
    Pick the latest file among basename, basename_2, basename_3, ... in `directory`.
    'Latest' = the one with the highest numeric suffix (no suffix = 0).
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
    # Load full matrix (variable width), keep only first row, coerce non-finite to 0
    P = np.loadtxt(pathfile, delimiter="\t", ndmin=2)
    if P.ndim == 1:
        P = P.reshape(1, -1)
    P[~np.isfinite(P)] = 0.0
    seq = P[0, :]
    seq = seq[seq != 0]  # drop zeros
    if seq.size < 2:
        raise ValueError("Path too short (need at least 2 nodes).")
    return seq


def _load_coor_matrix(coorfile: str) -> np.ndarray:
    coor = np.loadtxt(coorfile)
    if coor.ndim != 2 or coor.shape[1] < 10:
        raise ValueError("coor_file must have at least 10 columns.")
    return coor


def _optional_load_atom(atomfile: str):
    try:
        _ = np.loadtxt(atomfile)
    except Exception:
        print(f"⚠️  Could not read atom_file '{atomfile}'. Continuing.")


def encode_node(resid: int, chainid: int) -> float:
    # MATLAB formatting: rid + cid/10 with one decimal
    return float(resid) + float(chainid) / 10.0


def main():
    # --- pick latest files
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

    # --- load data
    seq = _load_path_first_row(in_pathfile)
    n = seq.size

    coor = _load_coor_matrix(in_coorfile)
    # coor: col1=resID, col6-8=CAxyz, col10=chainID  (MATLAB-style)
    resID_list   = coor[:, 0].astype(int)
    chainID_list = coor[:, 9].astype(int)
    CAxyz        = coor[:, 5:8]

    if in_atomfile:
        _optional_load_atom(in_atomfile)

    # Decode nodes v = resID + chainID/10
    resID_seq   = np.floor(seq + 1e-6).astype(int)
    chainID_seq = np.round((seq - resID_seq) * 10).astype(int)

    # --- residues on selected chain (sorted unique)
    mask_chain = (chainID_seq == CHAIN_ID_SELECTED)
    res_on_chain = np.unique(resID_seq[mask_chain])
    Nchain = res_on_chain.size
    if Nchain == 0:
        raise ValueError(f"No residues from chain {CHAIN_ID_SELECTED} found in the selected path.")

    # --- mapping (resID, chainID) -> row index in coor
    keyCoor = resID_list.astype(np.int64) * 100 + chainID_list.astype(np.int64)
    coorIdxMap = {int(k): int(i) for i, k in enumerate(keyCoor)}

    # --- step geometry & validity across the path
    dist = np.zeros(n - 1, dtype=float)
    invalid = np.zeros(n - 1, dtype=bool)
    for k in range(n - 1):
        ra, ca = int(resID_seq[k]), int(chainID_seq[k])
        rb, cb = int(resID_seq[k + 1]), int(chainID_seq[k + 1])
        k1 = ra * 100 + ca
        k2 = rb * 100 + cb
        if k1 not in coorIdxMap or k2 not in coorIdxMap:
            invalid[k] = True
            dist[k] = 0.0
        else:
            va = CAxyz[coorIdxMap[k1], :]
            vb = CAxyz[coorIdxMap[k2], :]
            dist[k] = float(np.linalg.norm(vb - va))

    # prefix sums S (distance) and C (invalid count)
    S = np.empty(n, dtype=float)
    S[0] = 0.0
    if n > 1:
        S[1:] = np.cumsum(dist)
    C = np.empty(n, dtype=int)
    C[0] = 0
    if n > 1:
        C[1:] = np.cumsum(invalid.astype(int))
    maxC = int(np.max(C))

    # --- occurrence indices for each residue on the selected chain
    occIdx = {}
    for rid in res_on_chain:
        occ = np.where((resID_seq == rid) & (chainID_seq == CHAIN_ID_SELECTED))[0]
        # Defensive clamp (shouldn't be needed but safe)
        occ = occ[(occ >= 0) & (occ < n)]
        occIdx[rid] = occ

    # --- residues → 0..N-1 index
    residues_chain = res_on_chain.copy()
    N = residues_chain.size
    resid2idx = {int(rid): i for i, rid in enumerate(residues_chain)}

    # --- containers for best segments (all ordered pairs)
    bestD   = np.full((N, N), np.inf, dtype=float)
    bestLen = np.zeros((N, N), dtype=int)
    bestI   = np.full((N, N), -1, dtype=int)   # path start index
    bestJ   = np.full((N, N), -1, dtype=int)   # path end index

    # Pre-allocated arrays reused for each start residue
    dBestForJ = np.full(n, np.inf, dtype=float)
    iBestForJ = np.full(n, -1, dtype=int)
    lenForJ   = np.zeros(n, dtype=int)

    bestS_forC = np.full(maxC + 1, -np.inf, dtype=float)
    bestI_forC = np.full(maxC + 1, -1, dtype=int)

    # ================= core O(N * n) algorithm =================
    for ai, sRID in enumerate(residues_chain):
        sIdx = occIdx[sRID]
        if sIdx.size == 0:
            continue

        # reset per-j bests and C-level tables
        dBestForJ[:] = np.inf
        iBestForJ[:] = -1
        lenForJ[:]   = 0
        bestS_forC[:] = -np.inf
        bestI_forC[:] = -1

        p = 0  # pointer into sIdx
        for t in range(n):
            # 1) If t is a start occurrence of sRID, update "best start so far" for C-level
            while p < sIdx.size and sIdx[p] == t:
                c = int(C[t])  # 0..maxC
                if c < 0 or c > maxC:
                    p += 1
                    continue
                s_val = S[t]
                sb    = bestS_forC[c]
                ib    = bestI_forC[c]
                if (s_val > sb + TOL) or (abs(s_val - sb) <= TOL and (ib < 0 or t > ib)):
                    bestS_forC[c] = s_val
                    bestI_forC[c] = t
                p += 1

            # 2) If t can be an end (on selected chain), try best compatible start at same C(t)
            if chainID_seq[t] != CHAIN_ID_SELECTED:
                continue
            c = int(C[t])
            if c < 0 or c > maxC:
                continue
            ib = bestI_forC[c]
            if ib < 0 or ib >= n:
                continue

            # candidate segment [ib .. t]
            d   = S[t] - bestS_forC[c]
            ln  = t - ib + 1
            # compare against per-position best for t (tie by shorter len)
            if (d + TOL < dBestForJ[t]) or (abs(d - dBestForJ[t]) <= TOL and ln < lenForJ[t]):
                dBestForJ[t] = d
                iBestForJ[t] = ib
                lenForJ[t]   = ln

        # 3) Finalize best per end residue eRID for this start residue sRID
        for bi, eRID in enumerate(residues_chain):
            if ai == bi:
                continue
            js = occIdx[eRID]
            if js.size == 0:
                continue

            # clamp occurrence indices just in case
            js = js[(js >= 0) & (js < n)]
            if js.size == 0:
                continue

            dvals = dBestForJ[js]
            pos = int(np.argmin(dvals))
            dmin = float(dvals[pos])
            if not np.isfinite(dmin):
                continue

            # tie-break on length
            ties = np.where(np.abs(dvals - dmin) <= TOL)[0]
            if ties.size > 1:
                tie_js   = js[ties]
                tie_lens = lenForJ[tie_js]
                pos2 = int(np.argmin(tie_lens))
                jStar = int(tie_js[pos2])
            else:
                jStar = int(js[pos])

            if jStar < 0 or jStar >= n:
                continue
            ib = iBestForJ[jStar]
            if ib < 0 or ib >= n:
                continue

            bestD[ai, bi]   = dmin
            bestLen[ai, bi] = int(lenForJ[jStar])
            bestI[ai, bi]   = ib
            bestJ[ai, bi]   = jStar
    # ===========================================================

    # --- build output matrix M: rows = N*(N-1), columns = [start end dist nodes...]
    maxNodes = 0
    for a in range(N):
        for b in range(N):
            if a == b:
                continue
            ib = bestI[a, b]
            jb = bestJ[a, b]
            if ib < 0 or jb < 0 or ib >= n or jb >= n or ib > jb:
                continue
            maxNodes = max(maxNodes, jb - ib + 1)

    if maxNodes == 0:
        maxNodes = 1  # avoid zero-width arrays

    totalRows = N * (N - 1)
    M = np.zeros((totalRows, 3 + maxNodes), dtype=float)
    row = 0
    for a in range(N):
        for b in range(N):
            if a == b:
                row += 1
                continue

            startVal = encode_node(int(residues_chain[a]), CHAIN_ID_SELECTED)
            endVal   = encode_node(int(residues_chain[b]), CHAIN_ID_SELECTED)
            M[row, 0] = startVal
            M[row, 1] = endVal

            ib = bestI[a, b]
            jb = bestJ[a, b]
            if ib < 0 or jb < 0 or ib >= n or jb >= n or ib > jb or not np.isfinite(bestD[a, b]):
                M[row, 2] = np.inf
            else:
                seg = seq[ib:jb+1]
                M[row, 2] = bestD[a, b]
                M[row, 3:3+seg.size] = seg
            row += 1

    # sort rows by start then end ascending
    order = np.lexsort((M[:, 1], M[:, 0]))
    M = M[order, :]

    # --- write paths_file with 6-decimal distance (col 3), nodes 1 decimal
    with open(PATHS_OUT, "w") as f:
        for r in range(M.shape[0]):
            start_node = M[r, 0]
            end_node   = M[r, 1]
            dist_val   = M[r, 2]

            f.write(f" {start_node:4.1f}\t{end_node:4.1f}\t")
            if np.isfinite(dist_val):
                f.write(f"{dist_val:.6f}")
            else:
                f.write("Inf")

            tail = M[r, 3:]
            nz = np.nonzero(tail != 0)[0]
            if nz.size > 0:
                tail = tail[:nz[-1] + 1]
                for v in tail:
                    f.write(f"\t{v:4.1f}")
            f.write("\n")

    # --- build symmetric “shortest of the pair” distance matrix for closeness
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

    # closeness centrality
    closenessVals = np.zeros(Nres, dtype=float)
    for i in range(Nres):
        rowVals = distMat[i, :]
        finiteVals = rowVals[np.isfinite(rowVals)]  # includes 0 on diagonal
        s = float(np.sum(finiteVals))
        closenessVals[i] = (Nres - 1) / s if s > 0 else 0.0

    # --- write closeness file
    with open(CLOSENESS_OUT, "w") as f:
        for i in range(Nres):
            f.write(f" {residues[i]:4.1f}\t{closenessVals[i]:.6f}\n")

    print(f"Wrote {PATHS_OUT} (distances with 6 decimals) and {CLOSENESS_OUT} (closeness with 6 decimals).")


if __name__ == "__main__":
    main()
