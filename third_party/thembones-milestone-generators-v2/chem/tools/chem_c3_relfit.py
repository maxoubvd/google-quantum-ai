"""chem_c3_relfit.py: builder C3 (fifth holder, Nov 18 05:45Z). KILL TEST of a POST-HOC bit-1 relation on a FINISHED,
GAUGED table: for loader k, the relation 'XOR over slots j of output (j, 1) = T(address)' with T = f(address & 15) XOR
g(address >> 4) (the family F, free for the loader [M R52]). The rows' present parity p(r) = XOR_j bit1(w[r,5k+j]) XOR
bit1(w[r,5k-5+j]) is fitted by the best (f, g) (exhaustive over f's 2^16 tables, g by majority per block); every row left
mismatched needs ONE bit-1 flip of one of its ten words: free when one of them is DEAD (r62_chem_cache_dead.dead_mask:
never read by the decoder; bit 1 leaves the carrier bits 5-6 and the wire bit 0 alone), else a LIVE flip (priced later).
Prints per loader: mismatched rows, rows repaired by a dead word, rows that need a live flip.
    python -I chem_c3_relfit.py SUBMISSION.npz --ks 4,6,8,10,12,14 [--exclude-words 5k:5k+5 ...]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r62_chem_cache_dead import dead_mask  # noqa: E402  (read-only import)


def fit(p, addr):
    """best f (16 bits, exhaustive) and g (per block, majority) for parity p over rows with addresses addr."""
    lo = addr & 15
    blk = addr >> 4
    nb = int(blk.max()) + 1
    best = None
    fs = ((np.arange(1 << 16)[:, None] >> np.arange(16)[None, :]) & 1).astype(np.int8)    # 65536 x 16
    # q[f, r] = p[r] XOR f[lo[r]]; per block the count of ones; g(b) = majority -> mismatches = min(ones, n - ones)
    q = (p[None, :] ^ fs[:, lo]).astype(np.int16)
    ones = np.zeros((fs.shape[0], nb), np.int32)
    for b in range(nb):
        m = blk == b
        ones[:, b] = q[:, m].sum(axis=1)
    nrow = np.bincount(blk, minlength=nb)
    mis = np.minimum(ones, nrow[None, :] - ones).sum(axis=1)
    i = int(np.argmin(mis))
    g = (ones[i] * 2 > nrow).astype(np.int8)
    t = fs[i][lo] ^ g[blk]
    return int(mis[i]), fs[i], g, t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sub")
    ap.add_argument("--ks", default="4,6,8,10,12,14")
    a = ap.parse_args()
    z = np.load(a.sub)
    w = np.asarray(z["angle_words"], np.int64)
    n = 1219
    w = w[:n]
    dead = dead_mask(np.asarray(z["angle_words"], np.int64), np.asarray(z["tree_edges"], np.int64), 7)[:n]
    addr = np.arange(n)
    tot_live = 0
    for k in [int(x) for x in a.ks.split(",")]:
        cols = [5 * k + j for j in range(5)] + [5 * (k - 1) + j for j in range(5)]
        p = np.zeros(n, np.int8)
        for c in cols:
            p ^= ((w[:, c] >> 1) & 1).astype(np.int8)
        mis, f, g, t = fit(p, addr)
        bad = (p ^ t).astype(bool)
        has_dead = dead[:, cols].any(axis=1)
        free = int((bad & has_dead).sum())
        live = int((bad & ~has_dead).sum())
        tot_live += live
        print(f"k{k}: rows {n}, mismatched after the best F target {mis} ({100 * mis / n:.1f}%), repaired by a dead word {free}, "
              f"need a LIVE flip {live}; rows with a dead word among the ten {int(has_dead.sum())}")
    print(f"total live flips {tot_live}")


if __name__ == "__main__":
    main()
