#!/usr/bin/env python3
"""r89_chem_cascade_dead.py (seat R89-CHEM-CASCADE): a FORCED relation space in the 23-output loader's cells-bound
group 4, paid by the pool kernel and the batch's dead cells, so that the cubic block's REMAINDER can be tiled.

FORM. C4's dead23 imposes m RANDOM column relations on group 4 (4 outputs, a = 14 live coordinates, 56 cells), then one
rank-additive cubic block (C4's tools/chem_c4_cubic.py) and tiles the rest. A random relation space of dim c folds onto
every fan-out pattern f with rank c, so each pattern keeps 14 - c dims; after the block (6 more generic dims) a remainder
of rank 44 (codim 12) keeps 2 dims a pattern = 15 planes for 22 needed (C4 measured 28 / 30 there). Here the relations
are a FULL SECOND COLUMN KILL on a live coordinate t2 (one relation e_o (x) t2 per output: the four fold to rank 1 on
EVERY pattern) plus g random relations. With g = 2 and the block the remainder's 12-dim relation space folds to rank
<= 1 + 2 + 6 = 9: every pattern keeps 5 dims (2 planes, 30 for 22) [P]: group 4 = 2 + 22 = 24.

USE (wrapper, C4's build unchanged): r89_chem_cascade_build.py replaces chem_c4_dead.solve_dead by solve_struct before
chem_c4_build.main() runs. Env: R89_G (generic relations, default 2), R89_T2 (comma list of candidate t2, default all
live), R89_TRIALS (draws per t2, default 4). Every solution is re-checked: the residual rank is recomputed (<= cells - m).
"""
import os

import numpy as np

import chem_loaders_core as core
import chem_c4_dead as cd

bits = core.bits
STATS = {"calls": 0, "solved": 0, "fail": 0, "t2": {}, "consec": 0}


class R89Abort(Exception):
    """Raised after R89_MAXFAIL consecutive misses: NOT caught by the engine's closure retry (RuntimeError /
    AssertionError), so a pool that cannot pay the form ends the build at once (the caller tries the next plan / pair)."""


def solve_struct(dec, outs, live, killed, R0, N, K, P, usable, sols, kern, nr, b, rows_of, ncell):
    """Drop-in for chem_c4_dead.solve_dead: same unknowns (kernel + pool-compensated dead cells), STRUCTURED relations
    (a full kill of coordinate t2 on every output + g random ones). Returns (chosen combos, picked cells) or None."""
    import chem_c4_drop as dr
    STATS["calls"] += 1
    cells_all = cd.HOOK["cells"]
    an = dec.an
    t = killed[0]
    liveu = np.uint32(live)
    gi_of = {o: gi for gi, o in enumerate(outs)}
    cells = [c for c in cells_all if c[0] in gi_of]
    if not cells:
        return None
    # pool t-columns (the first killed coordinate): C4's compensation, copied
    piv = {}
    for pi in range(len(usable)):
        v = 0
        for r in range(nr):
            if (int(P[pi][r]) >> t) & 1:
                v |= 1 << r
        tag = 1 << pi
        while v:
            hh = v.bit_length() - 1
            if hh not in piv:
                piv[hh] = (v, tag)
                break
            v ^= piv[hh][0]
            tag ^= piv[hh][1]

    def express(c):
        tag = 0
        while c:
            hh = c.bit_length() - 1
            if hh not in piv:
                return None
            c ^= piv[hh][0]
            tag ^= piv[hh][1]
        return tag
    unk = {gi: [] for gi in range(b)}
    pend = {}
    for c in cells:
        gi = gi_of[c[0]]
        E = cd._effect(an, nr, c[1], c[2], c[4])
        cvec = 0
        for r in range(nr):
            if (int(E[r]) >> t) & 1:
                cvec |= 1 << r
        tag = express(cvec)
        group = [c]
        if tag is None:
            if gi not in pend:
                pend[gi] = (E, cvec, c)
                continue
            E0, c0, cc0 = pend.pop(gi)
            E = E ^ E0
            tag = express(cvec ^ c0)
            group = [cc0, c]
            if tag is None:
                continue
        Ec = E.copy()
        for pi in bits(tag):
            Ec ^= P[pi]
        assert not ((Ec >> np.uint32(t)) & 1).any()
        unk[gi].append((E, tag, group, Ec & liveu))
    rng = dec.rng
    g = int(os.environ.get("R89_G", "2"))
    lc = [x for x in bits(live)]
    t2s = [int(x) for x in os.environ["R89_T2"].split(",")] if os.environ.get("R89_T2") else list(lc)
    t2s = [x for x in t2s if x in lc]
    rng.shuffle(t2s)
    t2s.sort(key=lambda x: -STATS["t2"].get(x, 0))     # the kill columns that paid before go first
    offs, n = [], 0
    for gi in range(b):
        offs.append(n)
        n += K + len(unk[gi])

    def rel_rows(vs):
        rows = []
        for r in range(nr):
            mask, rhs = 0, 0
            for gi in range(b):
                v = vs[gi]
                if not v:
                    continue
                rhs ^= cd.parity(int(R0[gi][r]) & v)
                for j in range(K):
                    if cd.parity(int(N[j][r]) & v):
                        mask |= 1 << (offs[gi] + j)
                for i, u in enumerate(unk[gi]):
                    if cd.parity(int(u[3][r]) & v):
                        mask |= 1 << (offs[gi] + K + i)
            rows.append((mask, rhs))
        return rows
    m_tot = b + g
    for t2 in t2s:
        kill = []
        for gi in range(b):
            vs = [0] * b
            vs[gi] = 1 << t2
            kill += rel_rows(vs)
        for _trial in range(int(os.environ.get("R89_TRIALS", "4"))):
            cons = list(kill)
            for _ in range(g):
                cons += rel_rows([rng.getrandbits(32) & live for _ in range(b)])
            x = cd._solve(cons, n, rng)
            if x is None:
                continue
            cols, chosen, picked = [], [], []
            for gi in range(b):
                col = R0[gi].copy()
                comb = sols[gi]
                xs = (x >> offs[gi]) & ((1 << (K + len(unk[gi]))) - 1)
                for j in range(K):
                    if xs >> j & 1:
                        col ^= N[j]
                        comb ^= kern[j]
                for i, u in enumerate(unk[gi]):
                    if xs >> (K + i) & 1:
                        col ^= u[3]
                        comb ^= u[1]
                        picked.append((outs[gi], u[0], u[2]))
                cols.append(col)
                chosen.append(comb)
            # exact re-checks: column t2 zero on every output; rank at most cells - m
            assert all(not ((cc >> np.uint32(t2)) & 1).any() for cc in cols), "kill of t2 broke"
            if dr._rank(rows_of(cols)) <= ncell - m_tot:
                STATS["solved"] += 1
                STATS["consec"] = 0
                STATS["t2"][t2] = STATS["t2"].get(t2, 0) + 1
                STATS["last"] = dict(t2=t2, g=g, rank=dr._rank(rows_of(cols)), toggles=len(picked))
                return chosen, picked
    STATS["fail"] += 1
    STATS["consec"] += 1
    # R89_STRICT=1 (default): a miss fails this attempt (the engine retries the group) instead of the silent fallback to
    # two random relations (group 4 then reads 26, worse than C4's dead23); R89_MAXFAIL misses in a row end the build
    if os.environ.get("R89_STRICT", "1") == "1":
        if STATS["consec"] >= int(os.environ.get("R89_MAXFAIL", "40")):
            raise R89Abort(f"r89: {STATS['consec']} consecutive misses of the structured relation space")
        raise RuntimeError("r89: no structured solve on this attempt")
    return None
