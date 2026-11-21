#!/usr/bin/env python3
"""r33_inv-wire2_rowlaw.py -- exhaustive check of the row law used by r33_inv-wire2_report.py (own code).

One non-restoring row with an odd divisor d = 2Q + 1 (Q in [2^m, 2^(m+1)): m variable bits between the constant
lowest bit and the constant leading one):  res = reg - s d,  s = +1 if reg >= 0 else -1,  with the invariant
|reg| <= 2 d before the row (so |res| <= d after it).

Claim [P]: the row needs m + 1 ANDs: the carry out of the lowest position is the register's own bit (the
divisor's lowest bit is the constant 1), one majority per variable divisor bit gives the carry into the
leading-one position, ONE more gives the carry out of it, and every result bit above is an XOR of wires on
all reachable states. With m ANDs only (no carry out of the leading-one position) the top bits are NOT affine
(printed as "m ANDs: False"). This is the line's feed row: 12 + 1 feed carry at m = 12 [M F1]; its no-feed
row costs m because there the lowest carry is a root wire. The check enumerates every (d, reg) for the given m
and solves, over GF(2), for each top result bit as an affine function of: the register's bits above the last
computed carry, that carry, the sign.

  r33_inv-wire2_rowlaw.py [m_lo m_hi]        default 3..9
Also for tapered rows: divisor T = rhu(d / 2^k) (lowest bit variable, the rounding bit is the carry-in):
m + 1 - k ANDs, the same statement for the top bits.
"""
import sys

import numpy as np


def gf2_affine(X, y):
    """is y an affine function of the columns of X over GF(2)? (X: n x k 0/1, y: n 0/1)"""
    A = np.concatenate([X, np.ones((X.shape[0], 1), dtype=np.uint8), y[:, None]], axis=1).astype(np.uint8)
    A = np.unique(A, axis=0)
    k = A.shape[1] - 1
    M = A.copy()
    r = 0
    for c in range(k):
        piv = None
        for i in range(r, M.shape[0]):
            if M[i, c]:
                piv = i
                break
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(M.shape[0]):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
    # inconsistent iff a row is 0..0 | 1
    bad = (M[:, :k].sum(axis=1) == 0) & (M[:, k] == 1)
    return not bool(bad.any())


def check_exact(m):
    W = m + 5                                     # register bits used here (two's complement)
    rows = []
    for Q in range(1 << m, 1 << (m + 1)):
        d = 2 * Q + 1
        reg = np.arange(-2 * d, 2 * d + 1, dtype=np.int64)
        s = np.where(reg >= 0, 1, -1)
        res = reg - s * d
        assert np.abs(res).max() <= d
        rows.append(np.stack([np.full_like(reg, d), reg, s, res], axis=1))
    A = np.concatenate(rows)
    d, reg, s, res = A[:, 0], A[:, 1], A[:, 2], A[:, 3]
    mask = (1 << W) - 1
    ru = reg & mask
    add = (-s * d) & mask                         # the addend as a W-bit pattern
    # carries: c[i] = carry into bit i of ru + add
    c = np.zeros((W + 1, len(ru)), dtype=np.int64)
    for i in range(W):
        c[i + 1] = (((ru >> i) & 1) + ((add >> i) & 1) + c[i]) >> 1
    # carry into bit 1 is the register's bit 0 (no AND)
    free_c1 = bool(np.array_equal(c[1], ru & 1))
    # ANDs: carries into bits 2 .. m+1 (one per variable divisor bit) -> m ANDs; + the carry into bit m+2
    out = {}
    for extra in (0, 1):
        top = m + 1 + extra                       # first position whose carry-in is the last one computed
        feats = [(ru >> i) & 1 for i in range(top, W)] + [c[top], (s < 0).astype(np.int64)]
        X = np.stack(feats, axis=1).astype(np.uint8)
        ok = True
        for i in range(top, W):
            y = ((res & mask) >> i) & 1
            ok &= gf2_affine(X, y.astype(np.uint8))
        out[extra] = ok
    return free_c1, out, len(ru)


def check_tapered(m, k):
    W = m + 5
    rows = []
    for Q in range(1 << m, 1 << (m + 1)):
        d = 2 * Q + 1
        T = max(1, (d + (1 << (k - 1))) >> k)
        Tp = max(1, (d + (1 << (k - 2))) >> (k - 1)) if k >= 2 else d     # the previous row's divisor
        reg = np.arange(-Tp, Tp + 1, dtype=np.int64)                        # |reg| <= previous divisor
        s = np.where(reg >= 0, 1, -1)
        res = reg - s * T
        rows.append(np.stack([np.full_like(reg, d), reg, s, res, np.full_like(reg, T)], axis=1))
    A = np.concatenate(rows)
    d, reg, s, res, T = A[:, 0], A[:, 1], A[:, 2], A[:, 3], A[:, 4]
    mask = (1 << W) - 1
    ru = reg & mask
    add = (-s * T) & mask
    c = np.zeros((W + 1, len(ru)), dtype=np.int64)
    for i in range(W):
        c[i + 1] = (((ru >> i) & 1) + ((add >> i) & 1) + c[i]) >> 1
    nvar = max(0, m + 1 - k)                      # variable bits of d >> k below its leading one
    best = None
    for extra in (0, 1, 2):
        top = nvar + extra                        # ANDs: carries into bits 1 .. top
        feats = [(ru >> i) & 1 for i in range(top, W)] + [c[top], (s < 0).astype(np.int64)]
        X = np.stack(feats, axis=1).astype(np.uint8)
        ok = True
        for i in range(top, W):
            y = ((res & mask) >> i) & 1
            ok &= gf2_affine(X, y.astype(np.uint8))
        if ok:
            best = top
            break
    return best, nvar


def main():
    lo, hi = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (3, 9)
    for m in range(lo, hi + 1):
        f, ok, n = check_exact(m)
        line = (f"m {m}: exact row: {n} states, carry out of the lowest bit = register bit: {f}; top bits affine "
                f"with {m} ANDs: {ok[0]}, with {m + 1} ANDs: {ok[1]}")
        tap = []
        for k in range(1, m + 3):
            best, nv = check_tapered(m, k)
            tap.append(f"k{k}:{best}")
        print(line + " | tapered rows (k: ANDs needed, None = more than variable bits + 2): " + " ".join(tap), flush=True)


if __name__ == "__main__":
    main()
