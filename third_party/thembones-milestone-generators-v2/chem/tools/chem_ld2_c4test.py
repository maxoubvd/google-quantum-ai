"""chem_ld2_c4test: does a table admit a 68-AND high cache on loader k? (builder C2, fifth holder; item C-4)

Necessary condition [P]: a 68-AND cache spans a 76-dim space H of functions of the 77 high blocks that holds the
constant and the 7 address bits, so H = {f : sum_{p in V} f(p) = 0} for a nonempty set V of blocks with |V| even and
XOR of its block indices 0. Every product factor's high part and every pure-high correction lies in H, so the V-sum
of every product (l1 + h1)(l2 + h2) is the constant sigma(h1 h2), and the V-sum of every pure term is 0. Hence for
every payload bit b of loader k's delta (words of batch k-1 XOR batch k; batch -1 = 0):
    S_b(l) = sum_{p in V} F_b(p, l)  is the same for all 16 low addresses l,
where F_b(76, l) for l >= 4 (addresses 1220..1231, never reached) is free. This tool solves that linear system over
GF(2) for V, in two cases (block 76 in V: the padding rows are free; not in V: all rows count), and prints the
dimension of the solution space (0 = no 68-AND cache exists for any plan on this table).

Usage: python -I chem_ld2_c4test.py <table.npz> [--bits 7] [--k 0-14]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chem_ld2_hb import HB, LAST  # noqa: E402  (tenth holder: 77 / 76 unless CHEM_HB is set)


def nullspace_gf2(M):
    """Basis of {v : M v = 0} over GF(2) (rows of the returned uint8 array)."""
    M = (np.array(M, dtype=np.uint8) & 1).copy()
    r, c = M.shape
    piv = []
    row = 0
    for col in range(c):
        if row >= r:
            break
        nz = np.nonzero(M[row:, col])[0]
        if len(nz) == 0:
            continue
        p = row + nz[0]
        if p != row:
            M[[row, p]] = M[[p, row]]
        others = np.nonzero(M[:, col])[0]
        others = others[others != row]
        M[others] ^= M[row]
        piv.append(col)
        row += 1
    free = [j for j in range(c) if j not in piv]
    basis = []
    for f in free:
        v = np.zeros(c, np.uint8)
        v[f] = 1
        for i, pc in enumerate(piv):
            v[pc] = M[i, f]
        basis.append(v)
    return np.array(basis, dtype=np.uint8).reshape(len(basis), c)


def delta_bits(words, k, nbits):
    """F[b, p, l] for b < 5*nbits over the 77 x 16 grid (padding rows = 0, marked by `valid`)."""
    A = words.shape[0]
    cur = words[:, 5 * k:5 * k + 5].astype(np.int64)
    prev = words[:, 5 * k - 5:5 * k].astype(np.int64) if k >= 1 else np.zeros_like(cur)
    d = prev ^ cur
    F = np.zeros((5 * nbits, HB * 16), np.uint8)
    for j in range(5):
        for b in range(nbits):
            F[j * nbits + b, :A] = (d[:, j] >> b) & 1
    valid = np.zeros(HB * 16, bool)
    valid[:A] = True
    return F.reshape(-1, HB, 16), valid.reshape(HB, 16)


def solve(F, valid, with76):
    """Rows of the system in v (HB unknowns): even size, XOR 0, V-sum constant in l on the counted rows."""
    rows = [np.ones(HB, np.uint8)]
    pidx = np.arange(HB)
    for j in range(7):
        rows.append(((pidx >> j) & 1).astype(np.uint8))
    nb = F.shape[0]
    for b in range(nb):
        for l in range(1, 16):
            if with76 and not valid[LAST, l]:
                continue                          # padding row of the last block is free: equation dropped
            rows.append(F[b, :, l] ^ F[b, :, 0])
    M = np.array(rows, np.uint8)
    if with76:
        # force v_last = 1 by solving the affine system: move the last column to the right-hand side
        N = nullspace_gf2(M)
        sols = [v for v in N]
        # the solution space of M v = 0 restricted to v_last = 1 is nonempty iff some basis vector has v_last = 1
        return len(N), any(v[LAST] for v in N)
    M0 = M.copy()
    M0[:, LAST] = 0
    M0 = np.vstack([M0, np.eye(HB, dtype=np.uint8)[LAST]])   # v_last = 0
    N = nullspace_gf2(M0)
    return len(N), len(N) > 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--bits", type=int, default=7)
    ap.add_argument("--k", default="0-14")
    a = ap.parse_args()
    with np.load(a.table, allow_pickle=False) as z:
        words = np.array(z["angle_words"]).astype(np.int64)
    lo, hi = (int(x) for x in a.k.split("-")) if "-" in a.k else (int(a.k), int(a.k))
    for k in range(lo, hi + 1):
        F, valid = delta_bits(words, k, a.bits)
        n1, ok1 = solve(F, valid, True)
        n0, ok0 = solve(F, valid, False)
        print(f"k{k}: V with block 76: nullity {n1}, solution with v76=1: {ok1}; V without 76: nullity {n0}",
              flush=True)


if __name__ == "__main__":
    sys.exit(main())
