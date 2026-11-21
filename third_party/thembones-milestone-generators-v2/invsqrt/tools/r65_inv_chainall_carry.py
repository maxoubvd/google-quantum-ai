#!/usr/bin/env python3
"""r65_inv_chainall_carry.py -- R65-INV-CHAINALL (designer; model only): THE CARRIES OF MAIN's TOP ROWS IN THE UNDO.

A row b of MAIN adds or subtracts d on the window (positions b..b+9) under the sign wire; its undo makes the same
eight carries c_2..c_9 again (c_(p+1) = MAJ(d_p, X_p, c_p) on the framed window X, c_1 = X_0).  In the undo game a
carry is erased for free only while its predecessor is live, which is what costs clean wires (rows 24, 23 pay 12.125
at two clean wires, 22: 10.25, 21: 9.5 against 8).  In the top rows the window holds 28 - b bits of U in ten wires,
so a carry may be ONE AND of two XORs of wires that are live when it has to go, WITHOUT its predecessor: then it is
X-measured and repaired by one CZ, at no Toffoli and no clean wire.
For every row b and carry i this tool writes the minimum-AND problem (F3's spec format) of c_i over the wires live
in the ascending two-wire schedule of the undo, at the moment c_i is no longer needed:
     v0..v6, sg (the row's control), n0..ni (window positions 0..i already restored to the pre-row value, framed),
     o(i+1)..o9 (positions not yet restored: the post-row value, framed), cn = c_(i+1) (just made; not for i = 9)
on ALL 3,133,452 inputs (R54's integer model of the chain, read-only import).
  r65_inv_chainall_carry.py SRC HTAB.npz OUTDIR [rows=24,23,22,21,20,19]
Safe to run: reads only; writes the specs under OUTDIR (my own folder)."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def chain(src, hpath):
    sys.argv[1] = src
    sys.path.insert(0, HERE)
    import r54_inv_tail_dc as DC
    RM, F8 = DC.RM, DC.F8
    x = np.arange(1, RM.MAXX + 1, dtype=np.int64)
    n, k, Q, d, delta = F8.ref_chain(x, 7)
    sg1, q1, _ = DC.div_trace(delta.astype(np.int64) << 5, d, DC.A_TOP, DC.S1)
    sg2, P, _ = DC.div2_true(q1, d)
    h, seen = DC.htab(hpath)
    y = -((P + 1) >> 1)
    T = y >> 9
    U = y + h[T + 1024]
    return d.astype(np.int64), U.astype(np.int64)


def main():
    src = os.path.abspath(os.path.expanduser(sys.argv[1]))
    hpath = os.path.expanduser(sys.argv[2])
    outdir = os.path.expanduser(sys.argv[3])
    opt = dict(a.split("=", 1) for a in sys.argv[4:])
    rows = [int(t) for t in opt.get("rows", "24,23,22,21,20,19").split(",")]
    os.makedirs(outdir, exist_ok=True)
    d, U = chain(src, hpath)
    v = (d - 257) >> 1
    N = (np.int64(1) << 35) + (U << 8)
    R = N.copy()
    summary = []
    for b in range(26, min(rows) - 1, -1):
        neg = (R < 0)
        Rn = R - np.where(neg, -1, 1) * (d << b)
        if b in rows:
            sg = neg.astype(np.int64)
            Wb = (R >> b) & 1023
            X = np.where(neg, Wb, (~Wb) & 1023)                 # framed window before the row's add
            Xp = (X + d) & 1023                                 # after the add (framed)
            Wa = np.where(neg, Xp, (~Xp) & 1023)
            assert np.array_equal(Wa, (Rn >> b) & 1023), "row %d: the framed add is not the row" % b
            c = [np.zeros(len(d), dtype=np.int64)]
            for p in range(9):
                dp = (d >> p) & 1
                xp = (X >> p) & 1
                c.append((dp & xp) | (dp & c[p]) | (xp & c[p]))
            # ---- mode "upper" (upper=1): the window AFTER the row (the undo's input), framed: its wire p as a
            # function of v, sg and wires below it (A: all of w0..w(p-1); B, for the top wire only: just the low
            # wires at register positions <= 27, i.e. w0..w(27-b)): what the erase of that wire costs when the
            # undo runs on a truncated window
            if opt.get("upper"):
                nlow = 28 - b
                for p in range(9, nlow - 1, -1):
                    for mode, nin in (("A", p), ("B", nlow)):
                        if mode == "B" and (p != 9 or nin == p):
                            continue
                        names = ["v%d" % j for j in range(7)] + ["sg"] + ["w%d" % q for q in range(nin)]
                        key = v | (sg << 7)
                        for q in range(nin):
                            key = key | (((Xp >> q) & 1) << (8 + q))
                        tv = (Xp >> p) & 1
                        order = np.argsort(key, kind="stable")
                        ks, ts = key[order], tv[order]
                        st = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
                        mx = np.maximum.reduceat(ts, st)
                        mn = np.minimum.reduceat(ts, st)
                        tag = "up%s_b%d_w%d" % (mode, b, p)
                        if not np.array_equal(mx, mn):
                            print("row %d wire w%d (%s, %d wires below): NOT a function (%d of %d patterns mixed)"
                                  % (b, p, mode, nin, int((mx != mn).sum()), len(st)), flush=True)
                            continue
                        rws = [[int(a_), int(b_)] for a_, b_ in zip(ks[st], mn)]
                        json.dump(dict(inputs=names, rows=rws), open(os.path.join(outdir, tag + ".json"), "w"))
                        print("row %d wire w%d (%s, %d wires below): %6d patterns, ones %6d -> %s.json" % (
                            b, p, mode, nin, len(rws), int(sum(r_[1] for r_ in rws)), tag), flush=True)
                R = Rn
                continue
            for i in range(2, 10):
                names = ["v%d" % j for j in range(7)] + ["sg"] + ["n%d" % p for p in range(i + 1)] + \
                        ["o%d" % p for p in range(i + 1, 10)] + (["cn"] if i < 9 else [])
                key = v.copy()
                pos = 7
                key |= sg << pos
                pos += 1
                for p in range(i + 1):
                    key |= ((X >> p) & 1) << pos
                    pos += 1
                for p in range(i + 1, 10):
                    key |= ((Xp >> p) & 1) << pos
                    pos += 1
                if i < 9:
                    key |= c[i + 1] << pos
                    pos += 1
                assert pos == len(names)
                order = np.argsort(key, kind="stable")
                ks, ts = key[order], c[i][order]
                st = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
                en = np.r_[st[1:], len(ks)]
                mx = np.maximum.reduceat(ts, st)
                mn = np.minimum.reduceat(ts, st)
                assert np.array_equal(mx, mn), "row %d carry %d is not a function of its context" % (b, i)
                rws = [[int(a_), int(b_)] for a_, b_ in zip(ks[st], mn)]
                # quick exact pre-checks on the care rows: constant / affine
                A = np.array([[(p_ >> j) & 1 for j in range(len(names))] + [1] for p_, _ in rws], dtype=np.uint8)
                t = np.array([f_ for _, f_ in rws], dtype=np.uint8)
                aff = gf2_solvable(A, t)
                tag = "car_b%d_c%d" % (b, i)
                json.dump(dict(inputs=names, rows=rws), open(os.path.join(outdir, tag + ".json"), "w"))
                ones = int(t.sum())
                summary.append((b, i, len(rws), ones, aff))
                print("row %d carry c%d: %6d patterns, ones %6d, affine in the context: %s -> %s.json" % (
                    b, i, len(rws), ones, aff, tag), flush=True)
        R = Rn
    json.dump(summary, open(os.path.join(outdir, "summary.json"), "w"))


def gf2_solvable(A, t):
    """is A z = t solvable over GF(2)?  (rows packed into Python ints)"""
    n = A.shape[1]
    rows = []
    for r in range(A.shape[0]):
        val = 0
        for j in range(n):
            if A[r, j]:
                val |= 1 << j
        rows.append(val | (int(t[r]) << n))
    basis = {}
    for val in set(rows):
        while val & ((1 << n) - 1):
            h = (val & ((1 << n) - 1)).bit_length()
            w = basis.get(h)
            if w is None:
                basis[h] = val
                break
            val ^= w
        else:
            if val >> n:
                return False
    return True


if __name__ == "__main__":
    main()
