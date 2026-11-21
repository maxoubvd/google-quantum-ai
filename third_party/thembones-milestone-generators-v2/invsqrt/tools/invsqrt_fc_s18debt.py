#!/usr/bin/env python3
"""invsqrt_fc_s18debt.py -- FC eighteenth holder: THE FORWARD'S CARRY DEBTS of MAIN row 24 in a FUSED UNDO
(R65-INV-CHAINALL-reprice's item: ROW9's chain block + row 24 undone as one map).  Row 24's narrow forward X-measures
its carries (FC's NW.fwd: debts[i] = hmr(carry)), so any undo that does not recompute them must pay each debt
(-1)^(m c_p) as a phase of the wires live then.  This tool computes, on ALL inputs (R54's integer model through the
designer's tools/r65_inv_chainall_carry.chain, read-only import), the carry c_p (p = 1..9) of row 24's framed add
(X + d, X = framed window before the row) and the EXACT GF(2) degree (<= 3, else '>= 4') of c_p as a function of:
   L1: v0..v6, sg, w0..w9         (R24's window, before the fused undo's step 1)
   L2: v0..v6, sg, u0..u3, w4..w9 (after step 1: U's four bits back in place)
   L3: v0..v6, sg, X0..X(p-1)     (the add's own inputs below p: what today's undo sees)
Degree 1-2 = free (Clifford repair), 3 = at least one non-Clifford gate under the debt bit (>= 0.5 a shot), >= 4 = at
least two (>= 1.0) [P: R65-INV-CHAINALL-reprice's degree rule].
  invsqrt_fc_s18debt.py SRC HTAB.npz [row=24]
Safe to run: reads only, prints only."""
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


class Basis:
    def __init__(self):
        self.b = {}

    def reduce(self, v):
        while v:
            w = self.b.get(v.bit_length())
            if w is None:
                return v
            v ^= w
        return 0

    def add(self, v):
        v = self.reduce(v)
        if v:
            self.b[v.bit_length()] = v


def fn_rows(cols, target):
    """distinct input patterns -> value; None if not a function."""
    key = np.zeros(len(target), dtype=np.int64)
    for j, c in enumerate(cols):
        key |= (c.astype(np.int64) & 1) << j
    order = np.argsort(key, kind="stable")
    ks, ts = key[order], target[order].astype(np.int64)
    st = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
    mx = np.maximum.reduceat(ts, st)
    mn = np.minimum.reduceat(ts, st)
    if not np.array_equal(mx, mn):
        return None
    return list(zip(ks[st].tolist(), mn.tolist()))


def degree(rows, n, dmax=3):
    cols = [0] * n
    t = 0
    for r, (p, y) in enumerate(rows):
        for j in range(n):
            if (p >> j) & 1:
                cols[j] |= 1 << r
        if y:
            t |= 1 << r
    ones = (1 << len(rows)) - 1
    B = Basis()
    B.add(ones)
    for j in range(n):
        B.add(cols[j])
    if B.reduce(t) == 0:
        return 1
    for deg in range(2, dmax + 1):
        for comb in itertools.combinations(range(n), deg):
            v = ones
            for j in comb:
                v &= cols[j]
            B.add(v)
        if B.reduce(t) == 0:
            return deg
    return dmax + 1


def main():
    src = os.path.abspath(os.path.expanduser(sys.argv[1]))
    hpath = os.path.expanduser(sys.argv[2])
    opt = dict(a.split("=", 1) for a in sys.argv[3:])
    row = int(opt.get("row", "24"))
    import r65_inv_chainall_carry as CT
    d, U = CT.chain(src, hpath)
    v = (d - 257) >> 1
    vb = [(v >> j) & 1 for j in range(7)]
    N = (np.int64(1) << 35) + (U << 8)
    R = N.copy()
    for b in range(26, row - 1, -1):
        neg = (R < 0)
        Rn = R - np.where(neg, -1, 1) * (d << b)
        if b == row:
            sg = neg.astype(np.int64)
            Wb = (R >> b) & 1023
            X = np.where(neg, Wb, (~Wb) & 1023)
            Xp = (X + d) & 1023
            Wa = (Rn >> b) & 1023
            assert np.array_equal(np.where(neg, Xp, (~Xp) & 1023), Wa)
            L = 28 - b
            ub = [((N >> (b + j)) & 1) for j in range(L)]
            wb = [(Wa >> p) & 1 for p in range(10)]
            xb = [(X >> p) & 1 for p in range(10)]
            car = [np.zeros(len(d), dtype=np.int64)]
            for p in range(9):
                dp = (d >> p) & 1
                car.append((dp & xb[p]) | (dp & car[p]) | (xb[p] & car[p]))
            print("row %d: %d inputs; carries of the framed add X + d (c_p = carry into position p)" % (b, len(d)))
            sets = {
                "L1 (v, sg, w0..w9)": vb + [sg] + wb,
                "L2 (v, sg, u0..u%d, w%d..w9)" % (L - 1, L): vb + [sg] + ub + wb[L:],
            }
            for p in range(1, 10):
                out = ["c%d" % p]
                for nm, cols in sets.items():
                    rows = fn_rows(cols, car[p])
                    out.append("%s: %s" % (nm.split()[0], "not a function" if rows is None else
                                           "deg %s (%d patterns)" % (
                                               {4: ">=4"}.get(degree(rows, len(cols)), degree(rows, len(cols))),
                                               len(rows))))
                cols3 = vb + [sg] + xb[:p]
                rows3 = fn_rows(cols3, car[p])
                out.append("L3: %s" % ("not a function" if rows3 is None else "deg %s" % (
                    {4: ">=4"}.get(degree(rows3, len(cols3)), degree(rows3, len(cols3))))))
                print("  " + "  |  ".join(out), flush=True)
        R = Rn


if __name__ == "__main__":
    main()
