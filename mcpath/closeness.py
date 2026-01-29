# closeness.py  (C kodundaki mantıkla uyumlu sürüm)
# -------------------------------------------------
# C ile uyum:
# - N: coor_file'daki tüm residüler (sıra bozulmaz)
# - l_ij: path üzerinde i->j için min(q-p) (q>p); yoksa l_ij = pl
# - closeness_i = (N-1) / sum_{j!=i} l_ij
# - path_file: C'deki gibi "flat" okunur (tüm sayılar), sadece ilk satır değil.
# - shortest_paths_all_pairs_chain1: nodes segmenti C gibi finali dahil etmez (uzunluk = q-p)

import os
import re
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

# --------------------------- configuration ---------------------------
IN_PATH_BASE = "path_file"
IN_COOR_BASE = "coor_file"

PATHS_FILE      = "shortest_paths_all_pairs_chain1"
CLOSENESS_FILE  = "closeness_float"
PEAKS_FILE      = "closeness_peaks"
PLOT_FILE       = "closeness_chain_plot.png"
LABELS_FILE     = "closeness_chain_labels.txt"

# C kodundaki E = 1e-6
E_C = 1e-6
# --------------------------------------------------------------------


def _pick_latest(basename: str, directory: str = ".") -> str:
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
    pdb_files_full = [os.path.join(directory, f) for f in pdb_files]
    pdb_files_full.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return pdb_files_full[0]


def _build_chain_mapping_from_pdb(coor: np.ndarray, directory: str = ".") -> dict[int, str]:
    chain_nums = sorted(set(coor[:, 9].astype(int)))
    pdb_path = _find_latest_pdb(directory)
    if pdb_path is None:
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

    mapping: dict[int, str] = {}
    for k in chain_nums:
        idx = k - 1
        mapping[k] = chain_letters_ordered[idx] if 0 <= idx < len(chain_letters_ordered) else str(k)
    return mapping


def _load_path_flat(pathfile: str) -> np.ndarray:
    """
    C kodu gibi: dosyadaki tüm float'ları (satır/kolon fark etmeden) sırayla oku.
    """
    P = np.loadtxt(pathfile, dtype=float)
    P = np.asarray(P, dtype=float).ravel()
    P[~np.isfinite(P)] = 0.0
    if P.size < 2:
        raise ValueError("Path too short (need at least 2 nodes).")
    return P


def _load_coor_matrix(coorfile: str) -> np.ndarray:
    coor = np.loadtxt(coorfile, dtype=float)
    if coor.ndim != 2 or coor.shape[1] < 10:
        raise ValueError(f"coor_file must have at least 10 columns (got shape {coor.shape})")
    return coor


def _node_key_from_encoded_float(v: float) -> int:
    """
    C: v = resid + chain/10
    Güvenli eşleşme için: key = round(10*v)
    (E=1e-6 << 0.1 olduğu için bu, C'deki fabs< E ile aynı amaca hizmet eder)
    """
    return int(np.round(v * 10.0))


def _min_forward_steps_C(posA: np.ndarray, posB: np.ndarray, pl: int):
    """
    C mantığı:
    - min(q-p) (q>p)
    - bulunamazsa minimum = pl
    Döndürür: (best_len, best_p, best_q, found)
    """
    if posA.size == 0 or posB.size == 0:
        return pl, -1, -1, False

    i = 0
    j = 0
    best_len = pl
    best_p = -1
    best_q = -1
    found = False

    # posA/posB zaten sorted
    while i < posA.size and j < posB.size:
        p = posA[i]
        q = posB[j]
        if q <= p:
            j += 1
        else:
            d = q - p
            if d < best_len:
                best_len = d
                best_p = p
                best_q = q
                found = True
            i += 1

    return best_len, best_p, best_q, found


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
    in_pathfile = _pick_latest(IN_PATH_BASE)
    in_coorfile = _pick_latest(IN_COOR_BASE)
    work_dir = os.path.dirname(os.path.abspath(in_coorfile)) or "."

    # ---- Load files ----
    pathi = _load_path_flat(in_pathfile)      # C gibi flat oku
    pl = int(pathi.size)

    coor = _load_coor_matrix(in_coorfile)     # tüm residüler
    chain_map = _build_chain_mapping_from_pdb(coor, directory=work_dir)

    # C'deki coor[i] = cor_1 + cor_10/10
    resID_cor = coor[:, 0].astype(int)
    chain_cor = coor[:, 9].astype(int)        # col10 (0-based 9)

    coor_enc = resID_cor.astype(float) + chain_cor.astype(float) / 10.0  # length N
    N = int(coor_enc.size)

    # Path sequence -> integer keys
    seq_keys = np.round(pathi * 10.0).astype(int)

    # Positions of each node key in the path
    pos_dict = defaultdict(list)
    for idx, k in enumerate(seq_keys):
        pos_dict[k].append(idx)
    for k in list(pos_dict.keys()):
        arr = np.asarray(pos_dict[k], dtype=int)
        arr.sort()
        pos_dict[k] = arr

    # Node keys for residues in coor (same encoding)
    node_keys = np.array([int(r * 10 + c) for r, c in zip(resID_cor, chain_cor)], dtype=int)
    pos_list = [pos_dict.get(k, np.asarray([], dtype=int)) for k in node_keys]

    # ---- Distances + indices (for PATHS_FILE) ----
    L = np.full((N, N), float(pl), dtype=float)   # C gibi default = pl
    start_idx = np.full((N, N), -1, dtype=int)
    end_idx   = np.full((N, N), -1, dtype=int)

    # diagonal = 0
    np.fill_diagonal(L, 0.0)

    for a in range(N):
        posA = pos_list[a]
        for b in range(N):
            if a == b:
                continue
            posB = pos_list[b]
            best_len, best_p, best_q, found = _min_forward_steps_C(posA, posB, pl=pl)
            L[a, b] = float(best_len)
            if found:
                start_idx[a, b] = best_p
                end_idx[a, b] = best_q
            # found değilse C'de best_initial/best_final tanımsız olurdu;
            # biz güvenli şekilde -1 bırakıyoruz ama mesafe pl kalıyor.

    # ---- Write PATHS_FILE (C segment mantığıyla: final dahil değil) ----
    # maxNodes = max distance (q-p) among found pairs
    maxNodes = 0
    for a in range(N):
        for b in range(N):
            if a == b:
                continue
            if start_idx[a, b] >= 0 and end_idx[a, b] > start_idx[a, b]:
                d = int(L[a, b])
                if d > maxNodes:
                    maxNodes = d

    totalRows = N * (N - 1)
    M = np.zeros((totalRows, 3 + maxNodes), dtype=float)
    row = 0
    for a in range(N):
        for b in range(N):
            if a == b:
                continue
            start_node = coor_enc[a]
            end_node   = coor_enc[b]
            dist_val   = L[a, b]

            M[row, 0] = start_node
            M[row, 1] = end_node
            M[row, 2] = dist_val

            ib = start_idx[a, b]
            if ib >= 0 and np.isfinite(dist_val) and dist_val > 0:
                d = int(dist_val)
                seg = pathi[ib:ib + d]  # C gibi: length = d (final düğüm yok)
                M[row, 3:3 + seg.size] = seg
            row += 1

    # (isteğe bağlı) sort: start then end
    order_paths = np.lexsort((M[:, 1], M[:, 0]))
    M = M[order_paths, :]

    with open(PATHS_FILE, "w") as fid:
        for r in range(M.shape[0]):
            start_node = M[r, 0]
            end_node   = M[r, 1]
            dist_val   = M[r, 2]
            fid.write(f" {start_node:4.1f}\t{end_node:4.1f}\t{dist_val:.6f}")
            tail = M[r, 3:]
            nz = np.nonzero(tail != 0)[0]
            if nz.size > 0:
                tail = tail[:nz[-1] + 1]
                for v in tail:
                    fid.write(f"\t{v:4.1f}")
            fid.write("\n")

    # ---- C-style closeness ----
    # sum_{j!=i} L[i,j]  (artık hepsi sonlu: default pl)
    sum_dists = np.sum(L, axis=1)  # includes diagonal (0)
    closenessVals = np.zeros(N, dtype=float)
    nonzero = sum_dists > 0
    closenessVals[nonzero] = (N - 1) / sum_dists[nonzero]

    # ---- Write closeness_float (C formatına yakın) ----
    with open(CLOSENESS_FILE, "w") as fid:
        for i in range(N):
            fid.write(f" {coor_enc[i]:5.1f}\t{closenessVals[i]:.6f}\n")

    # ---- Peaks / labels / plot (mevcut özellikler korunuyor) ----
    peaks_idx = []
    if closenessVals.size > 0:
        std_v = float(np.std(closenessVals))
        if std_v > 0:
            maxtab, _ = peakdet(closenessVals, std_v)
            if maxtab.size > 0:
                peaks_idx = maxtab[:, 0].astype(int).tolist()

    with open(PEAKS_FILE, "w") as fid:
        fid.write("# resID\tchainID\tencoded_node\tcloseness_peak\n")
        for idx in peaks_idx:
            fid.write(f"{resID_cor[idx]:4d}\t{chain_cor[idx]:d}\t{coor_enc[idx]:5.1f}\t{closenessVals[idx]:.6f}\n")

    # labels file: '12A' , 0.096112
    labels = []
    for r, c in zip(resID_cor, chain_cor):
        ch = chain_map.get(int(c), str(int(c)))
        labels.append(f"{int(r)}{ch}")

    with open(LABELS_FILE, "w") as lf:
        for lab, val in zip(labels, closenessVals):
            lf.write(f"'{lab}' , {val:.6f}\n")

    # Plot
    x = np.arange(N)
    plt.figure(figsize=(14, 4))
    plt.plot(x, closenessVals, marker="o", linestyle="-", linewidth=1.0, markersize=2.5, label="Closeness")
    if peaks_idx:
        peak_x = np.array(peaks_idx, dtype=int)
        peak_y = closenessVals[peak_x]
        plt.scatter(peak_x, peak_y, s=40, marker="s", label="Peaks")

    # x labels (sadeleştirilmiş)
    if N <= 80:
        tick_positions = x
    else:
        step = max(1, N // 80)
        tick_positions = x[::step]
    tick_labels = [labels[i] for i in tick_positions]
    plt.xticks(tick_positions, tick_labels, rotation=90, fontsize=6)

    plt.xlabel("Residue (coor_file order)")
    plt.ylabel("Closeness centrality (C-style)")
    plt.title("Closeness centrality (min forward steps; default=pl if not found)")
    if peaks_idx:
        plt.legend()
    plt.tight_layout()
    plt.savefig(PLOT_FILE, dpi=300, bbox_inches="tight")
    plt.close()

    if peaks_idx:
        print("Peak residues (resID + chain):")
        for idx in peaks_idx:
            ch = chain_map.get(int(chain_cor[idx]), str(int(chain_cor[idx])))
            print(f"  {resID_cor[idx]}{ch}  -> closeness = {closenessVals[idx]:.6f}")
    else:
        print("No peaks found with delta = std(closeness).")


if __name__ == "__main__":
    main()
