"""Sum of per-plane low x high core ranks of an ANGLE_k delta table (lowbits 4).

core rank of a plane = GF(2) rank of the H x 15 matrix of non-constant low-ANF
coefficients (Moebius transform over the 4 low address bits, one row per
16-address high block). Padded cells are set to 0 (upper bound).
Time O(planes * H * 15), memory O(rows * 75).
"""
import sys
import numpy as np

def mobius16(m):  # m: (H, 16) uint8 bits -> ANF coefficients over 4 low bits
    m = m.copy()
    for b in range(4):
        s = 1 << b
        for x in range(16):
            if x & s:
                m[:, x] ^= m[:, x ^ s]
    return m

def gf2_rank(rows):  # rows: iterable of ints
    basis = []
    for r in rows:
        for b in basis:
            r = min(r, r ^ b)
        if r:
            basis.append(r)
    return len(basis)

def core_ranks(words, J=5, brot=8):
    n, ncol = words.shape
    K = -(-ncol // J)
    padcols = K * J - ncol
    w = np.concatenate([words, np.zeros((n, padcols), dtype=words.dtype)], axis=1)
    H = -(-n // 16)
    out = []
    prev = np.zeros((n, J), dtype=np.int64)
    for k in range(K):
        cur = w[:, k * J:(k + 1) * J].astype(np.int64)
        delta = cur ^ prev
        prev = cur
        ranks = []
        for j in range(J):
            for bit in range(brot):
                p = ((delta[:, j] >> bit) & 1).astype(np.uint8)
                P = np.zeros(H * 16, dtype=np.uint8); P[:n] = p
                C = mobius16(P.reshape(H, 16))[:, 1:]
                ints = (C.astype(np.int64) << np.arange(15)).sum(1)
                ranks.append(gf2_rank(int(v) for v in ints if v))
        out.append(ranks)
    return np.array(out)  # (K, J*brot)

if __name__ == "__main__":
    d = np.load(sys.argv[1])
    A = d["angle_words"].astype(np.int64) & 255
    R, B = d["W"].shape[:2]
    rng = np.random.default_rng(0)
    def report(tag, tab):
        r = core_ranks(tab)
        tot = r.sum()
        print(f"{tag:42s} rows={tab.shape[0]:5d} H={-(-tab.shape[0]//16):3d} "
              f"sum_core_rank={tot:5d} per_loader={tot/r.shape[0]:6.1f} "
              f"model_T/pass={15*0 + 0.42*tot + r.shape[0]*(-(-tab.shape[0]//16)):7.0f} "
              f"full15={int((r==15).sum())} zero={int((r==0).sum())}")
        return r
    base = report("shipped table", A)
    print("per-bit mean rank (bit0..7):", np.round(base.reshape(15, 5, 8).mean(axis=(0, 1)), 2))
    report("random control (same shape)", rng.integers(0, 256, A.shape))
    one, two = A[:76], A[76:].reshape(R, B, -1)
    report("drop last group (7 groups, re-packed)", np.concatenate([one, two[:-1].reshape(-1, 75)]))
    report("drop 2 groups (6 groups)", np.concatenate([one, two[:-2].reshape(-1, 75)]))
    keep = np.sort(np.concatenate([rng.choice(B, 129, replace=False) for _ in range(R)]).reshape(R, 129), axis=1)
    report("B 143->129 (14 random rows/group dropped)", np.concatenate([one, np.stack([two[g, keep[g]] for g in range(R)]).reshape(-1, 75)]))
    report("first 512 rows only", A[:512])
