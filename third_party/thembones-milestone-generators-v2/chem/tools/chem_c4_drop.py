#!/usr/bin/env python3
"""chem_c4_drop.py: builder C4 (holder 3). THE KILL SYSTEM'S FREE BITS BUY ONE MORE ROW DEPENDENCY.

The engine's kill step (chem_loaders_core.Decomp.apply_kills, read-only) solves, per output of a group, the column
kills and d random row kills with the pool's fan-out bits and keeps ONE particular solution. When the pool has more
products than equations (a 24-output loader's group 4: 114 pool products against 76 + 2 x 14 equations) every output
keeps K free directions (the kernel of the pool's equation matrix, the same for every output). This module keeps the
kill step exact and spends those free bits on one more linear dependency among the group's residual rows: dim W drops
by one (73 -> 72 on a 24-output group 4: 37 -> 36 products, as on a 25-output loader whose group 4 gets a third row
kill from group 3's 42 products).

Search (exact, no estimate): fix random kernel choices x_o for the first b-1 outputs; the left kernel L of their
residual columns (rows x cells) holds the dependencies those outputs allow; for each lambda in L, the last output
needs lambda.R_last + sum_j x_j lambda.N_j = 0: a 14-bit linear system in its K free bits. A solution whose full
residual rank is below the base rank is applied (fan-out bits of the pool products, as apply_kills does).
Every result is re-checked: the kill equations hold and rank(W) is recomputed from the applied residual.

Use: `import chem_c4_drop; chem_c4_drop.install()` before the engine runs (chem_c4_build.py does it when C4_DROP=1).
C4_DROP_PAR=1 (default): try only when the base rank is odd (a drop then saves a product); C4_DROP_TRIALS (default 60).
"""
import os
import random

import numpy as np

import chem_loaders_core as core

STATS = {"calls": 0, "tried": 0, "dropped": 0, "nokernel": 0}
DEAD_ROW_HOOK = None      # holder 4: chem_c4_dead.solve_dead_row (row-bound groups)
DEAD_HOOK = None          # holder 4: chem_c4_dead.solve_dead (dead toggles as extra unknowns)


def _elim_with_comb(vecs):
    """Gaussian elimination over GF(2) of int vectors. Returns (pivots {bit: (vec, comb)}, kernel combos)."""
    piv, kern = {}, []
    for i, v in enumerate(vecs):
        c = 1 << i
        while v:
            h = v.bit_length() - 1
            e = piv.get(h)
            if e is None:
                break
            v ^= e[0]
            c ^= e[1]
        if v:
            piv[v.bit_length() - 1] = (v, c)
        else:
            kern.append(c)
    return piv, kern


def _reduce(piv, v):
    c = 0
    while v:
        h = v.bit_length() - 1
        e = piv.get(h)
        if e is None:
            return v, c
        v ^= e[0]
        c ^= e[1]
    return 0, c


def _rank(rows):
    piv = {}
    r = 0
    for v in rows:
        while v:
            h = v.bit_length() - 1
            e = piv.get(h)
            if e is None:
                piv[h] = v
                r += 1
                break
            v ^= e
    return r


def _xor_rows(P, comb):
    """XOR of the pool rows P[i] (uint32 arrays over block rows) for the set bits of comb."""
    idx = core.bits(comb)
    if not idx:
        return np.zeros(P.shape[1], np.uint32)
    return np.bitwise_xor.reduce(P[idx], axis=0)


def apply_kills_drop(self, outs, coords, dl, d, shift=0):
    """core.Decomp.apply_kills with the same equations, plus the rank-drop search on the kernel."""
    STATS["calls"] += 1
    an, nr = self.an, self.nr
    order = core.kill_order(coords, an)
    killed = order[shift:shift + dl] if shift + dl <= len(order) else order[:dl]
    live = coords
    for t in killed:
        live &= ~(1 << t)
    usable = [i for i, p in enumerate(self.prods) if not ((p[0] | p[2]) & ~coords)]
    if not usable:
        return None
    P = self.PM[usable]
    cols = [0] * len(usable)
    tgts = [0] * len(outs)
    Xg = self.X[:, outs].T
    width = 0
    for t in killed:
        for i, v in enumerate(core.rows_to_ints((P >> np.uint32(t)) & 1)):
            cols[i] |= v << width
        for i, v in enumerate(core.rows_to_ints((Xg >> np.uint32(t)) & 1)):
            tgts[i] |= v << width
        width += nr
    for _ in range(d):
        z = np.array([self.rng.random() < 0.5 for _ in range(nr)])
        if not z.any():
            z[self.rng.randrange(nr)] = True
        red = np.bitwise_xor.reduce(P[:, z], axis=1) & np.uint32(live)
        for i, v in enumerate(red.tolist()):
            cols[i] |= v << width
        redx = np.bitwise_xor.reduce(Xg[:, z], axis=1) & np.uint32(live)
        for i, v in enumerate(redx.tolist()):
            tgts[i] |= v << width
        width += an.nl
    piv, kern = _elim_with_comb(cols)
    sols = []
    for t in tgts:
        r, c = _reduce(piv, t)
        if r:
            return None
        sols.append(c)
    b = len(outs)
    liveu = np.uint32(live)
    # residual columns of the group with the particular solutions (uint32 per block row, live bits only)
    R0 = [((self.X[:, o] ^ _xor_rows(P, c)) & liveu) for o, c in zip(outs, sols)]
    K = len(kern)
    N = [(_xor_rows(P, k) & liveu) for k in kern]

    def rows_of(cols_):
        return [sum(int(cols_[g][M]) << (16 * g) for g in range(len(cols_))) for M in range(nr)]

    base = _rank(rows_of(R0))
    chosen = list(sols)
    ncell = b * bin(live).count("1")
    # holder 4: one line per kill system (group size, pool, free kernel K, base rank, live cells)
    STATS.setdefault("log", []).append(dict(b=b, pool=len(usable), K=K, base=base, ncell=ncell,
                                           killed=list(killed), d=d))
    if os.environ.get("C4_DROPDUMP") and base == ncell and not os.path.exists(os.environ["C4_DROPDUMP"]):
        # holder 4: dump a cells-bound group's kill system (particular residual R0, shared kernel N) for offline tests
        np.savez(os.environ["C4_DROPDUMP"], outs=np.array(outs), live=live, R0=np.stack(R0), N=np.stack(N) if N else
                 np.zeros((0, nr), np.uint32), killed=np.array(killed), nr=nr, P=P, usable=np.array(usable),
                 sols=np.array([str(c) for c in sols]), kern=np.array([str(c) for c in kern]))
    if K > 0 and base == ncell and b >= 2:
        # CELLS-BOUND group (rank = number of live cells < rows): the free bits impose COLUMN relations instead (m random
        # relations v: X'(x) v = 0 is linear in x, nr equations each); W = cells - m, products ceil(W / 2).
        STATS["tried"] += 1
        mmax = (K * b) // (nr - 1)
        m = 0
        for mm in range(1, mmax + 1):
            if (ncell - mm + 1) // 2 < (ncell - m + 1) // 2:
                m = mm
        if DEAD_HOOK is not None:
            self.dead_toggles = [c for c in getattr(self, "dead_toggles", []) if c[0] not in outs]
            res = DEAD_HOOK(self, outs, live, killed, R0, N, K, P, usable, sols, kern, nr, b, rows_of, ncell)
            if res is not None:
                chosen, picked = res
                for (o, E, group) in picked:            # the toggles change the target itself
                    self.X[:, o] ^= E
                    self.dead_toggles.extend(group)
                STATS["dropped"] += 1
                m = 0                                   # the plain column drop is replaced
        if m:
            got = _column_drop(self.rng, R0, N, K, b, live, nr, m, rows_of, ncell)
            if got is not None:
                STATS["dropped"] += 1
                STATS["coldrop"] = STATS.get("coldrop", 0) + m
                chosen = []
                for c, x in zip(sols, got):
                    cc = c
                    for j in core.bits(x):
                        cc ^= kern[j]
                    chosen.append(cc)
        want = False
    else:
        want = K > 0 and (base % 2 == 1 or os.environ.get("C4_DROP_PAR", "1") == "0")
        if DEAD_ROW_HOOK is not None and b >= 2:
            self.dead_toggles = [c for c in getattr(self, "dead_toggles", []) if c[0] not in outs]
            res = DEAD_ROW_HOOK(self, outs, live, killed, R0, N, K, P, usable, sols, kern, nr, b, rows_of, base)
            if res is not None:
                chosen, picked = res
                for (o, E, group) in picked:
                    self.X[:, o] ^= E
                    self.dead_toggles.extend(group)
                STATS["dropped"] += 1
                want = False
    if K == 0:
        STATS["nokernel"] += 1
    if want and b >= 2:
        STATS["tried"] += 1
        rng = self.rng
        trials = int(os.environ.get("C4_DROP_TRIALS", "60"))
        found = None
        for _ in range(trials):
            xs = [rng.getrandbits(K) for _ in range(b - 1)]
            cols_ = []
            for g in range(b - 1):
                v = R0[g].copy()
                for j in core.bits(xs[g]):
                    v ^= N[j]
                cols_.append(v)
            A = rows_of(cols_)
            _, lk = _elim_with_comb(A)            # left kernel: row combinations whose first b-1 outputs vanish
            if len(lk) > 14:
                lk = lk[:14]
            last = R0[b - 1]
            for m in range(1, 1 << len(lk)):
                lam = 0
                for q in core.bits(m):
                    lam ^= lk[q]
                rows_l = core.bits(lam)
                v = int(np.bitwise_xor.reduce(last[rows_l])) if rows_l else 0
                us = [int(np.bitwise_xor.reduce(Nj[rows_l])) for Nj in N]
                pv, _ = _elim_with_comb(us)
                r, c = _reduce(pv, v)
                if r:
                    continue
                # c combines the reduced u's; map back to the original u indices: re-solve directly
                x_last = _solve_comb(us, v)
                if x_last is None:
                    continue
                v_last = last.copy()
                for j in core.bits(x_last):
                    v_last ^= N[j]
                rk = _rank(rows_of(cols_ + [v_last]))
                if rk < base:
                    found = (xs + [x_last], rk)
                    break
            if found:
                break
        if found:
            STATS["dropped"] += 1
            xs, rk = found
            chosen = []
            for c, x in zip(sols, xs):
                cc = c
                for j in core.bits(x):
                    cc ^= kern[j]
                chosen.append(cc)
    for o, c in zip(outs, chosen):
        idx = core.bits(c)
        if idx:
            self.X[:, o] ^= np.bitwise_xor.reduce(P[idx], axis=0)
            for i in idx:
                self.prods[usable[i]][4] ^= 1 << o
    # exact re-check of the kill equations on the applied residual
    for o in outs:
        for t in killed:
            assert not ((self.X[:, o] >> np.uint32(t)) & 1).any(), "drop: a column kill broke"
    return live


def _column_drop(rng, R0, N, K, b, live, nr, m, rows_of, ncell, trials=8):
    """m random column relations v (over the b x a live cells) imposed exactly: unknowns x[o][j] (K a output), equations
    sum_o R0_o v_o + sum_{o,j} x[o][j] N_j v_o = 0 on every block row. Returns the list x[o] (int masks) or None."""
    lc = core.bits(live)
    for _ in range(trials):
        vs = [[rng.getrandbits(32) & live for _ in range(b)] for _ in range(m)]
        eqs = []                                    # (coefficient mask over b*K unknowns, rhs bit)
        for v in vs:
            # per block row r: parity(R0_o[r] & v_o) summed over o, and the same for each N_j
            r0 = np.zeros(nr, np.uint8)
            for o in range(b):
                r0 ^= _parity_rows(R0[o] & np.uint32(v[o]))
            colsj = []
            for o in range(b):
                for j in range(K):
                    colsj.append(_parity_rows(N[j] & np.uint32(v[o])))
            M = np.array(colsj, np.uint8).T          # (nr, b*K)
            for r in range(nr):
                eqs.append((int("".join("1" if t else "0" for t in M[r][::-1]), 2) if M[r].any() else 0, int(r0[r])))
        # Gaussian elimination with the rhs as the top bit
        nb = b * K
        piv = {}
        ok = True
        for coef, rhs in eqs:
            v = coef | (rhs << nb)
            while v & ((1 << nb) - 1):
                h = (v & ((1 << nb) - 1)).bit_length() - 1
                e = piv.get(h)
                if e is None:
                    piv[h] = v
                    break
                v ^= e
            else:
                if v >> nb & 1:
                    ok = False
                    break
        if not ok:
            continue
        # back-substitution: a particular solution (free unknowns 0)
        x = 0
        for h in sorted(piv):
            v = piv[h]
            rhs = v >> nb & 1
            rest = v & ((1 << nb) - 1) & ~(1 << h)
            val = rhs ^ (bin(rest & x).count("1") & 1)
            if val:
                x |= 1 << h
        xs = [(x >> (o * K)) & ((1 << K) - 1) for o in range(b)]
        cols_ = []
        for o in range(b):
            col = R0[o].copy()
            for j in core.bits(xs[o]):
                col ^= N[j]
            cols_.append(col)
        if _rank(rows_of(cols_)) <= ncell - m:
            return xs
    return None


def _parity_rows(arr):
    """parity of the set bits of each uint32 entry -> uint8 array."""
    a = arr.astype(np.uint32)
    a ^= a >> 16
    a ^= a >> 8
    a ^= a >> 4
    a ^= a >> 2
    a ^= a >> 1
    return (a & 1).astype(np.uint8)


def _solve_comb(us, v):
    """x with XOR_{j in x} us[j] == v, or None."""
    piv, _ = _elim_with_comb(us)
    r, c = _reduce(piv, v)
    if r:
        return None
    return c


def install():
    core.Decomp.apply_kills = apply_kills_drop
