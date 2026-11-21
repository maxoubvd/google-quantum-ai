"""chem_ld2_e17cap.py: builder C2 (seventeenth holder, Nov 17). A CAPPED LINE TABLE FOR MY ENGINE (route (a)).

Why: on a wired loader (25 generic outputs) the block cut is 4 / 4 / 4 / 3 low coordinates; a pair of points of one
4-coordinate block shares a fan-out space of dimension d up to 13 at codimension 48-51, and the enumeration keeps
d + d(d-1)/2 + (d-2) lines for it: 27,644 lines, and the exact scans stall for hours [M sixteenth holder].
What: after chem_own_engine._lines has enumerated the lines, every run of lines of one point pair (i, j) longer than
K = env C2_LINECAP keeps its d basis fan-outs and K - d of its other fan-outs (seeded by the pair and C2_CAPSEED);
the arrays are compacted in place. Every kept line is a real line, so the price stays an UPPER bound on the products
(the Lovasz rank with fewer lines can only be lower) and every extracted plan is valid; what is lost is choice.
env C2_MAXLINES raises the enumeration limit (60,000) so that late pairs are not cut off before the cap.
Off unless C2_LINECAP > 0 (chem_ld2_engine.install() calls patch()).
"""
from __future__ import annotations

import os

import numpy as np

_INV = {d * (d + 1) // 2 + d - 2: d for d in range(5, 64)}   # run length of a pair with d > 4 -> d
_SMALL = {1: 1, 3: 2, 7: 3, 15: 4}                            # d <= 4: all 2^d - 1 fan-outs, m = 1 .. 2^d - 1


def patch(eng):
    K = int(os.environ.get("C2_LINECAP", "0"))
    if K <= 0 or getattr(eng, "_e17_patched", False):
        return eng
    ml = int(os.environ.get("C2_MAXLINES", "0"))
    if ml > eng.MAXLINES:
        eng.MAXLINES = ml
    seed = int(os.environ.get("C2_CAPSEED", "0"))
    orig = eng._lines
    stats = {"calls": 0, "before": 0, "after": 0}

    def capped(pt_blk, pt_ymask, pt_outs, pt_phase, ncoord, Rm, ecount, nb, li, lj, ls, deg):
        n = orig(pt_blk, pt_ymask, pt_outs, pt_phase, ncoord, Rm, ecount, nb, li, lj, ls, deg)
        stats["calls"] += 1
        stats["before"] += n
        if n == 0:
            return n
        a, b = li[:n], lj[:n]
        brk = np.flatnonzero((a[1:] != a[:-1]) | (b[1:] != b[:-1])) + 1
        starts = np.r_[0, brk]
        lens = np.diff(np.r_[starts, n])
        big = np.flatnonzero(lens > K)
        if len(big) == 0:
            stats["after"] += n
            return n
        keep = np.ones(n, np.bool_)
        for q in big:
            st, ln = int(starts[q]), int(lens[q])
            if ln in _SMALL:
                basis = [0, 1, 3, 7][:_SMALL[ln]]
            elif ln in _INV:
                basis = list(range(_INV[ln]))
            else:
                basis = []                      # a run cut by the enumeration limit: sample only
            rest = np.array([x for x in range(ln) if x not in set(basis)], np.int64)
            m = max(0, K - len(basis))
            chosen = set(basis)
            if m and len(rest):
                r = np.random.default_rng((int(a[st]) * 1000003 + int(b[st]) * 7919 + seed) & 0x7fffffff)
                chosen |= set(r.choice(rest, size=min(m, len(rest)), replace=False).tolist())
            drop = [x for x in range(ln) if x not in chosen]
            keep[st + np.array(drop, np.int64)] = False
        idx = np.flatnonzero(keep)
        m = len(idx)
        li[:m] = li[idx]
        lj[:m] = lj[idx]
        ls[:m] = ls[idx]
        stats["after"] += m
        return m

    eng._lines = capped
    eng._e17_patched = True
    eng._e17_stats = stats
    return eng
