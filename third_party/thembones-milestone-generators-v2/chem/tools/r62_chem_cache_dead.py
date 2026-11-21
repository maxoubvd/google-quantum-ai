"""r62_chem_cache_dead: EXACT DON'T-CARE WORDS of a table (seat R62-CHEM-CACHE).

The verifier decodes a row by Givens rotations along the tree, in edge order, starting from e_0
(scoreinfra/dfthc_math.py reconstruct_u_vectors): parent <- cos(theta) parent, child <- sin(theta) parent.
A word whose PARENT amplitude is exactly zero when its edge is applied has no effect on the decoded vector:
its seven bits are free at zero error. That happens below an edge whose word is 0 or 64 (sin = 0: the child
subtree is dead) and after an edge whose word is 32 or 96 (cos = 0: the parent keeps nothing).
This prints, per column and per batch, how many words of the table are such exact don't-cares, and how many
MORE would be if every word with |signed value| <= t were rounded to zero (t = 1, 2): the price side of that
is the rounding builder's (not measured here).
Usage: python -I tools/r62_chem_cache_dead.py <submission.npz> [--brot 7]
"""
from __future__ import annotations

import argparse

import numpy as np


def dead_mask(w, edges, brot, thresh=0):
    """Boolean (rows, cols): True where the word is an exact don't-care (after optional zeroing |v| <= thresh)."""
    A, ncol = w.shape
    half = 1 << (brot - 1)          # pi
    quarter = half >> 1             # pi / 2
    N = ncol + 1
    alive = np.zeros((A, N), bool)
    alive[:, 0] = True
    dead = np.zeros((A, ncol), bool)
    full = 1 << brot
    for c, (p, ch) in enumerate(edges):
        v = w[:, c] % full
        par = alive[:, p].copy()
        dead[:, c] = ~par
        sv = ((v + quarter) % half) - quarter          # value mod pi in [-pi/2, pi/2)
        sin0 = (np.abs(sv) <= thresh) & (sv != -quarter)
        cos0 = (v % half) == quarter
        alive[:, ch] = par & ~sin0
        alive[:, p] = par & ~cos0
    return dead


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--brot", type=int, default=7)
    args = ap.parse_args()
    with np.load(args.table, allow_pickle=False) as z:
        w = np.array(z["angle_words"]).astype(np.int64)
        edges = np.array(z["tree_edges"]).astype(np.int64)
    A, ncol = w.shape
    print(f"table {args.table}: rows {A}, columns {ncol}")
    # tree shape: subtree sizes
    kids = {}
    for c, (p, ch) in enumerate(edges):
        kids.setdefault(int(p), []).append((c, int(ch)))
    below = np.zeros(ncol, int)      # number of later edges that die if column c's child is dead or parent dies

    def sub(v):
        n = 0
        for c, ch in kids.get(v, []):
            n += 1 + sub(ch)
        return n
    for c, (p, ch) in enumerate(edges):
        below[c] = sub(int(ch))
    print("edges under each column's child (subtree size):", below.tolist())
    for t in (0, 1, 2):
        d = dead_mask(w, edges, args.brot, t)
        tot = int(d.sum())
        per_batch = [int(d[:, 5 * k:5 * k + 5].sum()) for k in range(ncol // 5)]
        rows_any = int(d.any(axis=1).sum())
        tw = d[76:]
        print(f"threshold {t}: exact don't-care words {tot} of {d.size} ({100.0 * tot / d.size:.2f}%), rows with one "
              f"or more {rows_any}; one-body rows {int(d[:76].sum())}, two-body rows {int(tw.sum())}")
        print(f"   per batch (5 x {A} words each): {per_batch}")
        if t == 0:
            vals, cnt = np.unique(w % (1 << args.brot), return_counts=True)
            z0 = {int(v): int(c) for v, c in zip(vals, cnt) if int(v) in (0, 64, 32, 96)}
            print(f"   words equal to 0 / 64 / 32 / 96 in the table: {z0}")
            sv = ((w + 32) % 64) - 32
            print(f"   words with |value mod pi| <= 1: {int((np.abs(sv) <= 1).sum())}, <= 2: {int((np.abs(sv) <= 2).sum())}, "
                  f"<= 3: {int((np.abs(sv) <= 3).sum())}")
            print(f"   per column, words equal to zero mod pi: {(sv == 0).sum(axis=0).tolist()}")
        # the point stage's view: dead (row, word) pairs at the three low points D of a kept-product loader k.
        # A pair of loader k counts when the word is dead in batch k (the loader may then write anything there).
        import itertools
        low = np.arange(A) & 15
        cnt = np.bincount(d.sum(axis=1), minlength=76)
        top = np.sort(d.sum(axis=1))[::-1]
        print(f"   dead words a row: rows with 0: {int(cnt[0])}, 1-4: {int(cnt[1:5].sum())}, 5-9: {int(cnt[5:10].sum())}, "
              f"10-19: {int(cnt[10:20].sum())}, 20+: {int(cnt[20:].sum())}; the 228 fullest rows hold {int(top[:228].sum())}")
        for k in range(ncol // 5):
            dk = d[:, 5 * k:5 * k + 5]
            if k:
                dk = dk | d[:, 5 * k - 5:5 * k]          # either batch of the loader
            per_pt = np.array([dk[low == p].sum() for p in range(16)])
            rows_pt = np.array([(low == p).sum() for p in range(16)])
            best = max(itertools.combinations(range(16), 3), key=lambda D: per_pt[list(D)].sum())
            n = int(per_pt[list(best)].sum())
            den = int(rows_pt[list(best)].sum()) * 5
            print(f"   loader k{k}: dead pairs (either batch) {int(dk.sum())}; at D = 13,14,15: "
                  f"{int(per_pt[13:16].sum())} of {int(rows_pt[13:16].sum()) * 5}; best D {best}: {n} of {den} "
                  f"({100.0 * n / den:.1f}%)")


if __name__ == "__main__":
    main()
