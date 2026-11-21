"""r53_chem_gauge_group (design seat R53-CHEM-GAUGE): THE SIGN GAUGE OF A ROW AS IT ACTS ON WORDS, counted on a table.

The verifier decodes a row from 75 words along the tree (theta = pi * word / 2^(brot-1); edge t = (p, c), parent before
child): (amplitude left at p, amplitude sent to c) = r * (cos theta, sin theta), r the amplitude at p before the edge.
The table uses u u^T only. So the sign of every INTERMEDIATE amplitude is free: one bit per edge (74 edges whose
incoming amplitude is not the root's 1, plus the global sign) = 75 bits.

THE GROUP [P; checked numerically below on every row]: write h_t = 1 when word t takes a HALF TURN (word XOR 2^(brot-1))
and n_t = 1 when it is NEGATED (word -> -word mod 2^brot). Every gauge move is word_t -> (-1)^n_t word_t + 2^(brot-1) h_t
with h ANY pattern of 75 bits and
        n_t = sum of h over the LATER edges of the same parent p  +  sum of h over ALL edges whose parent is c   (mod 2)
(n = N h; the map h -> gauge move is a bijection). Negation keeps the lowest set bit of a word and complements every
bit above it: an odd word keeps bit 0 and flips bits 1..brot-1; a half turn flips the top bit only.
THE SIGNED RANGE (top two bits equal, theta in [-pi/2, pi/2)) IS h = 0: THE ONE GAUGE WITHOUT A HALF TURN. So no gauge
move at all is free for the loaders as they stand: each h_t = 1 breaks the copy 'bit 6 = bit 5' of word t.

    python r53_chem_gauge_group.py TABLE.npz [--brot 7] [--loaders 4,6,8,10,12,14] [--rand 20]
"""
import argparse
import sys

import numpy as np


def decode(words, edges, brot):
    th = np.pi * np.asarray(words, np.float64) / (1 << (brot - 1))
    v = np.zeros((th.shape[0], edges.shape[0] + 1))
    v[:, 0] = 1.0
    for t, (p, c) in enumerate(edges):
        a, b = v[:, p].copy(), v[:, c].copy()
        v[:, p] = np.cos(th[:, t]) * a - np.sin(th[:, t]) * b
        v[:, c] = np.sin(th[:, t]) * a + np.cos(th[:, t]) * b
    return v


def gauge_matrix(edges):
    """N (T x T over GF(2)): n = N h."""
    T = len(edges)
    N = np.zeros((T, T), np.uint8)
    for t, (p, c) in enumerate(edges):
        for u in range(T):
            if (u > t and edges[u][0] == p) or edges[u][0] == c:
                N[t, u] = 1
    return N


def rank2(M):
    M = (np.array(M, np.uint8) & 1).copy()
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]:
                piv = i
                break
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
        if r == rows:
            break
    return r


def flips(w, b, brot):
    """1 where negation flips bit b of the word (b above the lowest set bit)."""
    w = np.asarray(w, np.int64) % (1 << brot)
    low = w & -w                                          # lowest set bit (0 for w = 0)
    return ((low > 0) & ((1 << b) > low)).astype(np.uint8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--brot", type=int, default=7)
    ap.add_argument("--loaders", default="4,6,8,10,12,14")
    ap.add_argument("--rand", type=int, default=20)
    ap.add_argument("--order", default="", help="positions -> columns (default: the tree's own)")
    a = ap.parse_args()
    z = np.load(a.table)
    words = np.asarray(z["angle_words"], np.int64)
    edges = np.asarray(z["tree_edges"], np.int64)
    brot = a.brot
    full, half, q = 1 << brot, 1 << (brot - 1), 1 << (brot - 2)
    words = words % full
    R_, T = words.shape
    N = gauge_matrix(edges)
    print(f"table {a.table}: {R_} rows x {T} words, brot {brot}")
    # ---- anatomy of the tree
    par = edges[:, 0]
    chi = edges[:, 1]
    nedge_of = np.bincount(par, minlength=T + 1)
    last_of_parent = np.array([not np.any(par[t + 1:] == par[t]) for t in range(T)])
    leaf_child = nedge_of[chi] == 0
    both_final = last_of_parent & leaf_child
    print(f"tree: {int((nedge_of > 0).sum())} inner nodes, {int((nedge_of == 0).sum())} leaves; edges that are the LAST of their parent {int(last_of_parent.sum())}, "
          f"into a leaf {int(leaf_child.sum())}, both (their word cannot be negated: both outputs are final amplitudes) {int(both_final.sum())}")
    print(f"rank of N over GF(2) (independent negation patterns) {rank2(N)} of {T}; columns of N per word (words negated by one half turn): "
          f"min {int(N.sum(0).min())} mean {N.sum(0).mean():.2f} max {int(N.sum(0).max())}")
    for k in range(15):
        cols = np.arange(5 * k, 5 * k + 5)
        print(f"  batch {k:2d}: negatable words {int((~both_final[cols]).sum())} of 5; half turns that negate a word of this batch: {int((N[cols].sum(0) > 0).sum())}")
    # ---- exactness of the group on every row, random h
    rng = np.random.default_rng(0)
    v0 = decode(words, edges, brot)
    worst = 0.0
    for _ in range(a.rand):
        h = rng.integers(0, 2, (R_, T)).astype(np.int64)
        n = (h @ N.T.astype(np.int64)) & 1
        w1 = (np.where(n == 1, -words, words) + half * h) % full
        v1 = decode(w1, edges, brot)
        sg = np.sign((v0 * v1).sum(1))
        worst = max(worst, float(np.abs(v1 * sg[:, None] - v0).max()))
    print(f"EXACTNESS [M]: {a.rand} random half-turn patterns per row, n = N h: max |decoded new - (+-) decoded old| = {worst:.2e}")
    # a negation WITHOUT its half turns is not a gauge move: one example
    w1 = words.copy()
    odd_cols = np.where(~both_final)[0]
    w1[:, odd_cols[0]] = (-w1[:, odd_cols[0]]) % full
    v1 = decode(w1, edges, brot)
    sg = np.sign((v0 * v1).sum(1))
    print(f"control: word {int(odd_cols[0])} negated alone: max |change of the vector| = {float(np.abs(v1 * sg[:, None] - v0).max()):.3f} (not a gauge move)")
    # ---- signed range and special words
    top, nxt = (words >> (brot - 1)) & 1, (words >> (brot - 2)) & 1
    print(f"signed range: words with top two bits different {int((top != nxt).sum())} of {words.size}; words = 0: {int((words == 0).sum())}; "
          f"words = +-{q}: {int(((words == q) | (words == full - q)).sum())}")
    # amplitude entering each word (dead words: no amplitude: free)
    th = np.pi * words / half
    amp = np.zeros((R_, T))
    v = np.zeros((R_, T + 1))
    v[:, 0] = 1.0
    for t, (p, c) in enumerate(edges):
        amp[:, t] = np.abs(v[:, p])
        aa = v[:, p].copy()
        v[:, p] = np.cos(th[:, t]) * aa
        v[:, c] = np.sin(th[:, t]) * aa
    dead = amp < 1e-12
    print(f"dead words (no amplitude enters: any value decodes the same): {int(dead.sum())} of {dead.size} = {dead.mean() * 100:.2f}%; rows with any {int(dead.any(1).sum())}; "
          f"by batch {[int(dead[:, 5 * k:5 * k + 5].sum()) for k in range(15)]}")
    print(f"word parity: odd {((words & 1) == 1).mean():.4f}; = 2 mod 4 {((words & 3) == 2).mean():.4f}; = 0 mod 4 {((words & 3) == 0).mean():.4f}")
    # ---- how many bit-b relations of the pairs of the wired loaders a row can set by gauge
    if a.order:
        order = np.load(a.order).astype(np.int64)
    else:
        order = np.arange(T)
    loaders = [int(x) for x in a.loaders.split(",")]
    pairs = [(order[5 * k - 5 + j], order[5 * k + j]) for k in loaders for j in range(5)]
    for bits in ((1,), (1, 2), (1, 2, 3)):
        ranks_free, ranks_cols = [], []
        for r in range(R_):
            w = words[r]
            rows_ = []
            for (e, f) in pairs:
                for b in bits:
                    row = (flips(w[e], b, brot) * N[e]) ^ (flips(w[f], b, brot) * N[f])
                    rows_.append(row)
            M = np.array(rows_, np.uint8)
            ranks_free.append(rank2(M))
        ranks_free = np.array(ranks_free)
        nrel = len(pairs) * len(bits)
        hist = np.bincount(ranks_free, minlength=nrel + 1)
        nz = np.nonzero(hist)[0]
        print(f"bits {bits} on loaders {loaders}: {nrel} relations a row; RANK over GF(2) of the gauge on them, every half turn allowed "
              f"(each one breaks a copy bit 6 = bit 5): mean {ranks_free.mean():.2f} (two-body rows {ranks_free[76:].mean():.2f}, one-body {ranks_free[:76].mean():.2f}), "
              f"min {ranks_free.min()}, max {ranks_free.max()}; share of the relations {ranks_free.mean() / nrel:.3f}; histogram {dict((int(i), int(hist[i])) for i in nz)}")
    print("WITH THE SIGNED RANGE KEPT (h = 0): rank 0 on every row [P]: the gauge is spent on the 75 copies.")
    # the trade: a half-turn pattern h costs one generic output per (loader, slot) whose pair has h_e XOR h_f not affine in the address
    # count for the gauge that sets the most bit-1 relations: minimal number of half turns for the rank found (greedy basis)
    cost = []
    for r in rng.choice(np.arange(76, R_), size=min(200, R_ - 76), replace=False):
        w = words[r]
        M = np.array([(flips(w[e], 1, brot) * N[e]) ^ (flips(w[f], 1, brot) * N[f]) for (e, f) in pairs], np.uint8)
        # greedy: columns (half turns) by weight until the rank is reached
        full_rank = rank2(M)
        chosen = []
        cur = 0
        for c in np.argsort(-M.sum(0), kind="stable"):
            if cur == full_rank:
                break
            if rank2(M[:, chosen + [int(c)]]) > cur:
                chosen.append(int(c))
                cur += 1
        cost.append((full_rank, len(chosen), len(set(cc // 5 for cc in chosen))))
    cost = np.array(cost)
    print(f"THE TRADE [M 200 two-body rows]: to set rank {cost[:, 0].mean():.1f} bit-1 relations a row needs {cost[:, 1].mean():.1f} half turns "
          f"(= that many words leaving the signed range, in {cost[:, 2].mean():.1f} batches): one copy lost per relation gained, never fewer [P: rank].")


if __name__ == "__main__":
    main()
