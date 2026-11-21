"""chem_ld2_basis: a change of basis of the low space for the block engine of chem_own_engine (ours).

Idea (public, maoxin123 notebook "QCC v22 | 4,089-Toffoli chemistry loaders": "split the 15 low monomials into
five blocks of three and search the subspaces of each block"; read in their circuits by
tools/chem_ld2_dissect.py + chem_ld2_yfan.py: their five blocks are generic 3-dim subspaces of the low space, and
in coded circuits they keep the degree filtration). Nothing here reads a public file: the basis is drawn at random
(seeded) and chosen by our own search.

Conventions. The analysis `an` of chem_loaders_core expresses the mixed tensor of circuit k in low coordinates t
(basis functions f_t: address bits, carriers, computed monomials); an.X0[M, o] is a mask over t. A new basis is
g_u = sum_t P[u, t] f_t, given by the rows P[u] (masks over t). Then
    X0_new = X0_old P^-1   (row vectors),      y_old = y_new P   (a product's low part, for emission).
New coordinates keep the phase structure: g_u is live in phase p iff every f_t it uses is. Levels: level(u) =
max level of the f_t it uses. A P that is block lower-triangular in the levels ("filt") keeps every coded output
(mixed part of degree <= d) inside the new coordinates of level <= d, so codes keep their small supports.
"""
from __future__ import annotations

import copy
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
from chem_loaders_core import bits  # noqa: E402


_orig_make_blocks = core.make_blocks


def _make_blocks(an, supp, nb):
    """core.make_blocks, unless the analysis carries a fixed partition (blocks of a searched basis)."""
    fixed = getattr(an, "blocks_override", None)
    if fixed is not None:
        return [list(b) for b in fixed]
    return _orig_make_blocks(an, supp, nb)


core.make_blocks = _make_blocks


def mat_inv(rows, n):
    """Inverse of the GF(2) matrix with rows `rows` (masks over n columns); raises if singular."""
    A = list(rows)
    Inv = [1 << i for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if (A[r] >> col) & 1), None)
        if piv is None:
            raise ValueError("singular")
        A[col], A[piv] = A[piv], A[col]
        Inv[col], Inv[piv] = Inv[piv], Inv[col]
        for r in range(n):
            if r != col and (A[r] >> col) & 1:
                A[r] ^= A[col]
                Inv[r] ^= Inv[col]
    return Inv


def vec_mul(y, P):
    """Row vector y (mask) times matrix with rows P."""
    out = 0
    for u in bits(y):
        out ^= P[u]
    return out


def phase_sig(an, t):
    return tuple(bool(ym >> t & 1) for ym in an.ymasks)


def make_P(an, seed: int, mode: str = "filt", lower: float = 0.5):
    """Rows P[u] of a random basis change compatible with the phases of `an`.
    mode 'id': identity. 'filt': mixing inside (level, phase class) groups, plus random lower-level terms from
    coordinates live wherever u is live. 'full': mixing inside phase classes (levels ignored), plus terms from
    coordinates live wherever u is live."""
    n = an.nl
    if mode == "id":
        return [1 << u for u in range(n)]
    rng = random.Random(seed)
    groups = {}
    for t in range(n):
        key = (phase_sig(an, t),) + ((an.level[t],) if mode == "filt" else ())
        groups.setdefault(key, []).append(t)
    P = [0] * n
    for key, ts in groups.items():
        m = len(ts)
        while True:                       # random invertible m x m
            M = [rng.getrandbits(m) for _ in range(m)]
            try:
                mat_inv(M, m)
                break
            except ValueError:
                continue
        for i, u in enumerate(ts):
            v = 0
            for j in bits(M[i]):
                v |= 1 << ts[j]
            P[u] = v
    group_of = {t: key for key, ts in groups.items() for t in ts}
    if lower > 0:
        # elementary row additions P[u] += P[t] (always invertible); 'filt' only from lower levels (keeps the
        # filtration); both modes only rows whose support is live wherever u is live (keeps the phases)
        for u in range(n):
            su = phase_sig(an, u)
            for t in range(n):
                if t == u or group_of[t] == group_of[u]:
                    continue
                if mode == "filt" and an.level[t] >= an.level[u]:
                    continue
                live_ok = all(all(b or not a for a, b in zip(su, phase_sig(an, x))) for x in bits(P[t]))
                if not live_ok:
                    continue
                if rng.random() < lower:
                    P[u] ^= P[t]
    mat_inv(P, n)                          # check invertible
    return P


def dont_care_low(an):
    """Old-coordinate masks spanning the free part of the top row: the table ends inside the last block of
    2**lb addresses (nv valid addresses, nv a power of two), so the ANF coefficients of the top high monomial
    on the low monomials m >= nv are free (the circuit's values on the invalid addresses are don't-cares)."""
    A, L, H = an.A, an.L, an.H
    nv = A - L * (H - 1)
    if nv <= 0 or nv >= L:
        return []
    if nv & (nv - 1) == 0:
        return [an.mono_coords[m] for m in range(nv, L)]
    # tenth holder: any nv (even B leaves 12 valid rows in the top block): the coefficient of low monomial m is free
    # iff no valid address l < nv contains m (f(l) = sum over m <= l of c_m); for a power of two this is m >= nv
    return [an.mono_coords[m] for m in range(1, L) if all((m & l) != m for l in range(nv))]


def transform(an, P, free_top=False, blocks=None):
    """Copy of `an` in the new low coordinates (for the engine). The copy keeps the old analysis in .an_old.
    free_top: the top row of the mixed tensor is removed from what the engine sees (its relation pools are
    computed without it); the true row and the don't-care space are kept (.top_true, .dlow) so that the plan
    extraction (chem_ld2_engine.Engine2) can choose the don't-cares that put the row into the features' span."""
    n = an.nl
    Q = mat_inv(P, n)
    new = copy.copy(an)
    X = an.X0
    Xn = np.zeros_like(X)
    for t in range(n):
        bit = ((X >> np.uint32(t)) & np.uint32(1)).astype(bool)
        Xn[bit] ^= np.uint32(Q[t])
    new.free_top = False
    new.blocks_override = [list(b) for b in blocks] if blocks else None
    new.X0 = Xn
    new.level = [max(an.level[t] for t in bits(P[u])) for u in range(n)]
    new.ymasks = [sum(1 << u for u in range(n) if (P[u] & ~ym) == 0) for ym in an.ymasks]
    for u in range(n):
        assert any(ym >> u & 1 for ym in new.ymasks), "a new coordinate is live in no phase"
    new.P, new.Q, new.an_old = list(P), Q, an
    new.top_true = Xn[an.nr - 1].copy()
    new.top_info = None
    if free_top:
        new.top_info = fix_top_row(new, [vec_mul(v, Q) for v in dont_care_low(an)])
    return new


def fix_top_row(an2, dlow):
    """Choose the don't-cares of the top row (in the new coordinates) so that, in EVERY block of the engine,
    the top row lies in the span of the other rows: then each block's relation pool is that of the other rows
    alone (one relation more than with the zero-filled row), and any plan of the engine reproduces the table
    on every valid address. One linear solve: unknowns are the don't-care vectors z (x) e_o (z in the
    don't-care low space, only existing cells), equations r(top + d) = 0 for every pool relation r of the other
    rows of every block. Writes the new top row into an2.X0 and returns (equations, unknowns, solved)."""
    nr = an2.nr
    X = an2.X0
    nb = X.shape[1]
    bd = core.BlockDecomp(an2, 0, 1)
    bd.prepare()
    supp = bd.supp
    unknowns = []                                   # (o, z) with every cell (u, o) existing
    for o in range(nb):
        sp = core.Span()
        for z in dlow:
            zz = z
            if zz & ~supp[o]:
                continue
            if sp.add(zz):
                unknowns.append((o, zz))
    rows_eq = []                                    # (mask over unknowns, rhs)
    for bi in bd.binfo:
        cell = bi["cell"]
        ts = set(bi["ts"])
        W = core.Span()
        for v in bi["rows"][:nr - 1]:
            W.add(v)
        rel = core.nullspace([p[0] for p in W.piv.values()], bi["amb"])
        top = 0
        for o in range(nb):
            for u in bits(int(X[nr - 1, o])):
                if u in ts:
                    top |= 1 << cell[(u, o)]
        uvec = []
        for o, z in unknowns:
            v = 0
            for u in bits(z):
                if u in ts:
                    v |= 1 << cell[(u, o)]
            uvec.append(v)
        for r in rel:
            m = 0
            for i, v in enumerate(uvec):
                if bin(r & v).count("1") & 1:
                    m |= 1 << i
            rows_eq.append((m, bin(r & top).count("1") & 1))
    # solve m . c = rhs over GF(2)
    piv = {}
    for m, rhs in rows_eq:
        while m:
            h = m.bit_length() - 1
            if h in piv:
                m ^= piv[h][0]
                rhs ^= piv[h][1]
            else:
                piv[h] = (m, rhs)
                break
        if not m and rhs:
            return (len(rows_eq), len(unknowns), False)
    sol = 0
    for h in sorted(piv):                           # back substitution, lowest pivot first
        m, rhs = piv[h]
        rest = m & ~(1 << h)
        val = rhs ^ (bin(rest & sol).count("1") & 1)
        if val:
            sol |= 1 << h
    for i, (o, z) in enumerate(unknowns):
        if (sol >> i) & 1:
            X[nr - 1, o] ^= np.uint32(z)
    # check
    bd2 = core.BlockDecomp(an2, 0, 1)
    bd2.prepare()
    for bi in bd2.binfo:
        W = core.Span()
        for v in bi["rows"][:nr - 1]:
            W.add(v)
        assert W.reduce(bi["rows"][nr - 1])[0] == 0, "top row not in the span after the solve"
    return (len(rows_eq), len(unknowns), True)


def map_back(prods, P):
    """Products [y1, A1, y2, A2, s, phase] from new to old low coordinates."""
    out = []
    for y1, a1, y2, a2, s, ph in prods:
        out.append([vec_mul(y1, P), a1, vec_mul(y2, P) if y2 else 0, a2, s, ph])
    return out


def check_equivalent(an, an_new):
    """X0 of an_new times P equals X0 of an (round trip)."""
    P = an_new.P
    X = an_new.X0.copy()
    if getattr(an_new, "top_info", None) is not None:
        X[an.nr - 1] = an_new.top_true            # the top row's don't-cares were chosen; compare the rest
    back = np.zeros_like(X)
    for u in range(an.nl):
        bit = ((X >> np.uint32(u)) & np.uint32(1)).astype(bool)
        back[bit] ^= np.uint32(P[u])
    return bool(np.array_equal(back, an.X0))
