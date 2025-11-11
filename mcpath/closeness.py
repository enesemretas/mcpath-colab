# closeness.py
# Reads the most recent coor_file / atom_file / path_file (with optional _N suffix)
# Computes all-pairs shortest segments on a single long path for a selected chain,
# writes:
#   - "shortest_paths_all_pairs_chain1"  (distance with 6 decimals, nodes as 1 decimal)
#   - "closeness_float"                  (node, closeness with 6 decimals)

import os
import re
import math
import numpy as np

# --------------------------- configuration ---------------------------
IN_PATH_BASE = "path_file"                  # tab-separated; first row = long path
IN_COOR_BASE = "coor_file"                  # numeric; col1=resID, cols6-8=CA xyz, col10=chainID
IN_ATOM_BASE = "atom_file"                  # not used in math here, validated if present

PATHS_OUT     = "shortest_paths_all_pairs_chain1"
CLOSENESS_OUT = "closeness_float"

CHAIN_ID_SELECTED = 1                       # set which chain to analyze (e.g., 1)
TOL = 1e-9                                   # tie tolerance
# --------------------------------------------------------------------

def _pick_latest(basename: str, directory: str = ".") -> str:
    """
    Pick the latest file among basename, basename_2, basename_3, ... in `directory`.
    'Latest' = the one with the highest numeric suffix (no suffix = 0).
    If none exist, raises FileNotFoundError.
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
    # load full matrix (variable width), keep only first row, coerce non-finite to 0
    P = np.loadtxt(pathfile, delimiter="\t", ndmin=2)
    if P.ndim == 1:  # single row can become 1D
        P = P.reshape(1, -1)
    P[~np.isfinite(P)] = 0.0
    seq = P[0, :]
    # drop trailing zeros
    seq = seq[seq != 0]
    if seq.size < 2:
        raise ValueError("Path too short (need at least 2 nodes).")
    return seq

def _load_coor_matrix(coorfile: str) -> np.ndarray:
    coor = np.loadtxt(coorfile)
    if coor.ndim != 2 or coor.shape[1] < 10:
        raise ValueError("coor_file must have at least 10 columns.")
    return coor

def _optional_load_atom(atomfile: str):
    # Not used in the computation; we just confirm it exists/readable if present
    try:
        _ = np.loadtxt(atomfile)
    except Exception:
        # If it can't be read, just warn to console; not critical for outputs
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

    print(f"Using files:")
    print(f"  path_file  -> {in_pathfile}")
    print(f"  coor_file  -> {in_coorfile}")
    if in_atomfile:
        print(f"  atom_file  -> {in_atomfile}")

    # --- load data
    seq = _load_path_first_row(in_pathfile)
    n = seq.size

    coor = _load_coor_matrix(in_coorfile)
    # Columns: (MATLAB 1-based) col1=resID, col6-8=CAxyz, col10=chainID
    resID_list   = coor[:, 0].astype(int)
    chainID_list = coor[:, 9].astype(int)
    CAxyz        = coor[:, 5:8]  # (x,y,z)

    if in_atomfile:
        _optional_load_atom(in_atomfile)

    # Decode nodes v = resID + chainID/10
    resID_seq   = np.floor(seq + 1e-6).astype(int)
    chainID_seq = np.round((seq - resID_seq) * 10).astype(int)

    # --- residues on selected chain (sorted unique)
    res_on_chain = np.unique(resID_seq[chainID_seq == CHAIN_ID_SELECTED])
    Nchain = res_on_chain.size
    if Nchain == 0:
        raise ValueError(f"No residues from chain {CHAIN_ID_SELECTED} found in the selected path.")

    # --- build mapping from (resID, chainID) -> row index in coor
    # key = resID*100 + chainID  (matches MATLAB int64/100 trick)
    keyCoor = (resID_list.astype(np.int64) * 100 + chainID_list.astype(np.int64))
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
    S = np.concatenate([[0.0], np.cumsum(dist)])              # prefix distances (len n)
    C = np.concatenate([[0], np.cumsum(invalid).astype(int)]) # prefix invalid-counts (len n)
    maxC = int(np.max(C))

    # Occurrence indices for each residue on the selected chain
    occIdx = {}
    for rid in res_on_chain:
        occIdx[rid] = np.where((resID_seq == rid) & (chainID_seq == CHAIN_ID_SELECTED))[0]

    # --- fast shortest search for all ordered pairs (i != j)
    N = Nchain
    bestI   = np.full((N, N), np.nan)
    bestJ   = np.full((N, N), np.nan)
    bestD   = np.full((N, N), np.inf)
    bestLen = np.zeros((N, N), dtype=int)

    dBestForJ = np.full(n, np.inf)
    iBestForJ = np.full(n, np.nan)
    lenForJ   = np.zeros(n, dtype=int)

    bestS_forC = np.full(maxC + 1, -np.inf)  # best start S for each C level
    bestI_forC = np.full(maxC + 1, np.nan)   # index of best start

    for ai, sRID in enumerate(res_on_chain):
        sIdx = occIdx[sRID]   # sorted indices where sRID occurs

        # reset per-j results
        dBestForJ[:] = np.inf
        iBestForJ[:] = np.nan
        lenForJ[:]   = 0
        bestS_forC[:] = -np.inf
        bestI_forC[:] = np.nan

        p = 0
        for t in range(n):
            # update candidate starts (if t is an occurrence of sRID)
            while p < sIdx.size and sIdx[p] == t:
                c = int(C[t])  # 0..maxC
                s = S[t]
                sb = bestS_forC[c]
                ib = bestI_forC[c]
                if (s > sb + TOL) or (abs(s - sb) <= TOL and (np.isnan(ib) or t > ib)):
                    bestS_forC[c] = s
                    bestI_forC[c] = float(t)
                p += 1

            # if t is an end on selected chain, try best compatible start at same C(t)
            if chainID_seq[t] == CHAIN_ID_SELECTED:
                c = int(C[t])
                ib = bestI_forC[c]
                if not np.isnan(ib):
                    ib_i = int(ib)
                    d   = S[t] - bestS_forC[c]
                    ln  = t - ib_i + 1
                    if (d + TOL < dBestForJ[t]) or (abs(d - dBestForJ[t]) <= TOL and ln < lenForJ[t]):
                        dBestForJ[t] = d
                        iBestForJ[t] = float(ib_i)
                        lenForJ[t]   = ln

        # finalize best for each end residue
        for bi, eRID in enumerate(res_on_chain):
            if ai == bi:
                continue
            js = occIdx[eRID]
            if js.size == 0:
                continue

            dvals = dBestForJ[js]
            pos = int(np.argmin(dvals))
            dmin = dvals[pos]
            if not np.isfinite(dmin):
                continue

            # tie by shortest len
            ties = np.where(np.abs(dvals - dmin) <= TOL)[0]
            if ties.size > 1:
                tie_js = js[ties]
                tie_lens = lenForJ[tie_js]
                pos2 = int(np.argmin(tie_lens))
                jStar = int(tie_js[pos2])
                dmin  = float(dBestForJ[jStar])
            else:
                jStar = int(js[pos])

            bestD[ai, bi]   = float(dmin)
            bestI[ai, bi]   = float(int(iBestForJ[jStar]))
            bestJ[ai, bi]   = float(jStar)
            bestLen[ai, bi] = int(lenForJ[jStar])

    # --- build output matrix M: rows = N*(N-1), columns = [start end dist nodes...]
    maxNodes = 0
    for a in range(N):
        for b in range(N):
            if a == b or math.isnan(bestI[a, b]) or math.isnan(bestJ[a, b]):
                continue
            maxNodes = max(maxNodes, int(bestJ[a, b] - bestI[a, b] + 1))

    totalRows = N * (N - 1)
    M = np.zeros((totalRows, 3 + maxNodes), dtype=float)
    row = 0
    for a in range(N):
        for b in range(N):
            if a == b:
                row += 1
                continue
            startVal = encode_node(int(res_on_chain[a]), CHAIN_ID_SELECTED)
            endVal   = encode_node(int(res_on_chain[b]), CHAIN_ID_SELECTED)
            M[row, 0] = startVal
            M[row, 1] = endVal

            if math.isnan(bestI[a, b]) or math.isnan(bestJ[a, b]):
                M[row, 2] = np.inf
            else:
                ib = int(bestI[a, b])
                jb = int(bestJ[a, b])
                seg = seq[ib:jb+1]
                M[row, 2] = bestD[a, b]
                M[row, 3:3+seg.size] = seg
            row += 1

    # sort rows by start then end ascending
    # argsort by column 0 then column 1
    order = np.lexsort((M[:,1], M[:,0]))
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
    distMat = np.minimum(distMat, distMat.T)  # enforce symmetry

    # closeness centrality: (Nres-1)/sum_j d(i,j) over finite distances
    closenessVals = np.zeros(Nres, dtype=float)
    for i in range(Nres):
        rowVals = distMat[i, :]
        finiteVals = rowVals[np.isfinite(rowVals)]  # includes 0 on diagonal
        s = float(np.sum(finiteVals))
        closenessVals[i] = (Nres - 1) / s if s > 0 else 0.0

    # --- write closeness file: "  xx.x\t0.xxxxxx"
    with open(CLOSENESS_OUT, "w") as f:
        for i in range(Nres):
            f.write(f" {residues[i]:4.1f}\t{closenessVals[i]:.6f}\n")

    print(f"Wrote {PATHS_OUT} (distances with 6 decimals) and {CLOSENESS_OUT} (closeness with 6 decimals).")

if __name__ == "__main__":
    main()
