#!/usr/bin/env python3
"""chem_c4_cubic.py: builder C4 (holder 4). A RANK-ADDITIVE CUBIC BLOCK inside one generic group of a wired loader.

THE BLOCK (two Toffolis): P0 = (A1 + y2)(A2 + y1) on a spare wire (P0 feeds no output), P1 = (A + y)(B + z + P0) fanned
into a pattern s of the group's outputs, then P0 is erased by measurement (HMR + CZ, no Toffoli). Expanding P1 (products
of functions pointwise; u' = u without its constant) gives the mixed part
    C = A (x) z' + B' (x) y + (A A1) (x) y1 + (A A2) (x) y2 + A1 (x) (y y1) + A2 (x) (y y2),   z' = z + y1 y2, B' = B + A1 A2,
a rank-6 matrix from two Toffolis (a plain product gives rank 2 a Toffoli).

RANK-ADDITIVE [P]: with the group's residual R = X Y^T (X a basis of its column space, Y of its row space W), C = U V^T is
rank-additive (rank(R - C) = rank(R) - 6) iff U = X alpha, V = Y beta and beta^T alpha = I_6. Construction, exact and
linear step by step: the six flats (z', y, y1, y2, y y1, y y2 on s) in W; then A, A1, A2, B' each by one linear system
(column space, the high cache's H = {V-sum 0}, the six core conditions, A A1 and A A2 in the column space); P1's sigma (the
closure's V-sum of its pure-high part A B') is set to 0, so the plan's closure is untouched. The remainder R - C is then
tiled by the engine's own group solve. Every step is re-checked: P1's mixed part recomputed from its truth table over all
77 x 16 addresses, and rank(R - C) recomputed.

Use: C4_CUBIC=<generic group index> in tools/chem_c4_build.py (installs the hook after chem_c4_drop); emission by
emit_c4x() below (a copy of builder C2's chem_ld2_c4.emit_c4, read-only, with the block added).
"""
import os
import random

import numpy as np

import chem_loaders_core as core
import chem_ld2_c4 as c4

bits = core.bits
HB = c4.HB
_HT = c4._HT                       # truth of x^M over the HB blocks
STATS = {"calls": 0, "built": 0, "fail": {}}


# ------------------------------------------------------------------------------------------------- GF(2) helpers
def _eliminate(vecs):
    """Echelon form with combination tags: (pivots {bit: (vec, tag)}, list of kernel tags)."""
    piv, kern = {}, []
    for i, v in enumerate(vecs):
        tag = 1 << i
        while v:
            h = v.bit_length() - 1
            e = piv.get(h)
            if e is None:
                break
            v ^= e[0]
            tag ^= e[1]
        if v:
            piv[v.bit_length() - 1] = (v, tag)
        else:
            kern.append(tag)
    return piv, kern


def _reduce(piv, v):
    tag = 0
    while v:
        h = v.bit_length() - 1
        e = piv.get(h)
        if e is None:
            return v, tag
        v ^= e[0]
        tag ^= e[1]
    return 0, tag


def _rank(vecs):
    return len(_eliminate(vecs)[0])


def _solve_affine(cons, n, rng):
    """cons: list of (mask over n unknowns, rhs bit). A uniformly random solution (int), or None if inconsistent."""
    piv = {}
    for mask, rhs in cons:
        v = mask | (rhs << n)
        while v & ((1 << n) - 1):
            h = (v & ((1 << n) - 1)).bit_length() - 1
            e = piv.get(h)
            if e is None:
                piv[h] = v
                break
            v ^= e
        else:
            if v >> n & 1:
                return None
    free = [i for i in range(n) if i not in piv]
    x = 0
    for i in free:
        if rng.random() < 0.5:
            x |= 1 << i
    for h in sorted(piv):                       # pivots in increasing order: each row's other bits are lower or free
        v = piv[h]
        rest = v & ((1 << n) - 1) & ~(1 << h)
        val = (v >> n & 1) ^ (bin(rest & x).count("1") & 1)
        if val:
            x |= 1 << h
    for mask, rhs in cons:                      # exact re-check
        assert (bin(mask & x).count("1") & 1) == rhs
    return x


def _dot(a, b):
    return bin(a & b).count("1") & 1


# ------------------------------------------------------------------------------------------ functions of h and l
def _mobius_int(t, nbits):
    """ANF of a truth table given as an int (bit p = value at point p, p < 2^nbits); returns the coefficient int."""
    arr = np.array([(t >> p) & 1 for p in range(1 << nbits)], np.uint8)
    core.mobius(arr, nbits)
    return core.to_int(arr)


def h_truth(coef):
    """Truth table (HB-bit int) of the high function with coefficient vector coef (bit r = monomial r + 1)."""
    return c4.hmask_truth(coef)


def h_coef(t):
    """Coefficient vector (bit r = monomial r + 1) and constant of an HB-bit truth table (down-set {0..HB-1})."""
    a = _mobius_int(t, 7)
    const = a & 1
    coef = 0
    for M in range(1, HB):
        if a >> M & 1:
            coef |= 1 << (M - 1)
    return coef, const


class LowSide:
    """The loader's 15 low basis functions: truth tables over l = 0..15 and coordinates of any function of l."""

    def __init__(self, an):
        self.an = an
        self.nl = an.nl
        self.T = []
        for bs in an.basis:
            t = 0
            for m in bits(bs["anf"]):
                for l in range(16):
                    if (l & m) == m:
                        t ^= 1 << l
            self.T.append(t)

    def truth(self, y):
        t = 0
        for i in bits(y):
            t ^= self.T[i]
        return t

    def coords(self, t):
        """Coordinates (mask over the basis) of the low function with truth table t, constant dropped."""
        a = _mobius_int(t, 4)
        y = 0
        for m in range(1, 16):
            if a >> m & 1:
                y ^= self.an.mono_coords[m]
        return y

    def mult(self, y, v):
        return self.coords(self.truth(y) & self.truth(v))


def mixed_rows(an, low, F):
    """The engine's representation (nr rows of coordinate masks) of the mixed part of F, F[h][l] over HB x 16."""
    nr = an.nr
    G = np.zeros((128, 16), np.uint8)
    G[:HB, :] = F
    G = G.T.copy()                              # (16, 128): Moebius along h for every l
    core.mobius(G, 7)
    rows = []
    for M in range(1, nr + 1):
        t = 0
        for l in range(16):
            if G[l, M]:
                t |= 1 << l
        rows.append(low.coords(t))
    return rows


# ------------------------------------------------------------------------------------------------- the block
def make_block(dec, outs, live, rng, tries=200):
    """Construct a rank-additive, closure-neutral cubic block on the group (outs, live) of Decomp dec.
    Returns dict(block=..., rows=..., s=..., rank0, rank1) or None."""
    STATS["calls"] += 1
    an = dec.an
    nr = an.nr
    lc = bits(live)
    a = len(lc)
    b = len(outs)
    loc = {t: i for i, t in enumerate(lc)}
    low = LowSide(an)
    Xg = dec.X[:, outs]
    rows = []
    for r in range(nr):
        v = 0
        for gi in range(b):
            x = int(Xg[r, gi])
            for t in bits(x):
                v |= 1 << (gi * a + loc[t])
        rows.append(v)
    ncell = a * b
    piv_r, _ = _eliminate(rows)
    rank0 = len(piv_r)
    if rank0 < 8:
        STATS["fail"]["rank"] = STATS["fail"].get("rank", 0) + 1
        return None
    # row space W: basis Y (independent rows) and its relations (W-perp)
    Ybasis = [v for v, _ in piv_r.values()]
    pivY, _ = _eliminate(Ybasis)
    rel = core.nullspace(Ybasis, ncell)
    # ONE rank factorization R = X Y^T: Y = the basis Ybasis of W, X[r] = coordinates of row r in Ybasis; the columns of X
    # are then a basis of R's column space (nr-bit ints)
    Xrow = []
    for r in range(nr):
        rem, tag = _reduce(pivY, rows[r])
        assert rem == 0
        Xrow.append(tag)
    Xbasis = []
    for i in range(rank0):
        v = 0
        for r in range(nr):
            if Xrow[r] >> i & 1:
                v |= 1 << r
        Xbasis.append(v)
    pivX, kx = _eliminate(Xbasis)
    assert not kx, "X columns dependent"
    lker = core.nullspace(Xbasis, nr)          # z with z . x = 0 for every column: u in colspace iff z . u = 0
    # coordinates in X: pick rank0 pivot rows where X is invertible; coords(u) = Xp^-1 u_P (valid on the column space)
    Pset = sorted(pivX.keys())                  # pivot bit positions of an echelon basis of colspace
    E = [pivX[p][0] for p in Pset]              # echelon basis vectors (combinations of Xbasis given by tags)
    Etag = [pivX[p][1] for p in Pset]
    # in the echelon basis E (each E_i has leading bit Pset[i]), coordinates of u are found by reduction; as linear forms:
    # form_i(u) for u in colspace = the coefficient of E_i, computed as a mask over u's bits by back substitution
    m = len(E)
    order = sorted(range(m), key=lambda i: -Pset[i])   # leading bits from high to low
    form = [0] * m
    for i in order:                             # coefficient of E_i = bit Pset[i] of u minus contributions of higher E_j
        f = 1 << Pset[i]
        for j in order:
            if j == i:
                break
            if E[j] >> Pset[i] & 1:
                f ^= form[j]
        form[i] = f
    # convert forms from the echelon basis E to the basis Xbasis: u = sum_i c_i E_i, E_i = sum_k tag_ik Xbasis_k
    xform = [0] * m                             # xform[k]: mask over u's bits giving the coefficient of Xbasis_k
    for i in range(m):
        for k in bits(Etag[i]):
            xform[k] ^= form[i]

    def coordsX_form(beta):
        """Linear form (mask over nr bits) u -> beta . coordsX(u), valid for u in the column space."""
        f = 0
        for k in bits(beta):
            f ^= xform[k]
        return f

    def coordsY(flat):
        r, tag = _reduce(pivY, flat)
        if r:
            return None
        return tag                               # tag: combination of Ybasis indices

    # high side helpers
    Vset = dec.V if getattr(dec, "V", None) is not None else 0
    hv = 0                                       # V-sum functional on coefficient vectors
    for r in range(nr):
        if bin(_HT[r + 1] & Vset).count("1") & 1:
            hv |= 1 << r
    HTr = [_HT[r + 1] for r in range(nr)]

    def mult_matrix(A):
        """Columns of u -> coef(T_A * T_u) (constant dropped), as nr ints."""
        TA = h_truth(A)
        colsM = []
        for r in range(nr):
            co, _ = h_coef(TA & HTr[r])
            colsM.append(co)
        return colsM

    def pull(mask, colsM):
        """Linear form f(M u) as a mask on u."""
        out = 0
        for r in range(nr):
            if _dot(mask, colsM[r]):
                out |= 1 << r
        return out

    for attempt in range(tries):
        # ---- low side
        s = rng.randrange(1, 1 << b)
        gis = bits(s)
        smask = []                               # allowed v (over live local coords): v . sum_{gi in s} r_gi = 0
        for rr in rel:
            mm = 0
            for gi in gis:
                mm ^= (rr >> (gi * a)) & ((1 << a) - 1)
            smask.append(mm)
        nonlive = sum(1 << t for t in range(an.nl) if t not in loc)

        def local(v):                            # global coordinate mask (live only) -> local a-bit vector
            return sum(1 << loc[t] for t in bits(v))

        def allowed_cons_global(extra_map=None):
            """Constraints on a global coordinate vector v (15 bits): zero off the live set; flat in W."""
            cons = [(1 << t, 0) for t in bits(nonlive)]
            for mm in smask:
                g = sum(1 << lc[ci] for ci in bits(mm))
                cons.append((g, 0))
            return cons

        base = allowed_cons_global()
        y = _solve_affine(base, an.nl, rng)
        if not y:
            STATS["fail"]["y"] = STATS["fail"].get("y", 0) + 1
            continue

        def cons_mult(yy):                       # constraints on v such that (yy * v) is allowed too
            Mcol = [low.mult(yy, 1 << t) for t in range(an.nl)]
            out = []
            for cm, rhs in base:
                f = 0
                for t in range(an.nl):
                    if _dot(cm, Mcol[t]):
                        f |= 1 << t
                out.append((f, rhs))
            return out

        cm_y = cons_mult(y)
        y1 = _solve_affine(base + cm_y, an.nl, rng)
        y2 = _solve_affine(base + cm_y, an.nl, rng)
        zp = _solve_affine(base, an.nl, rng)
        if not (y1 and y2 and zp):
            STATS["fail"]["low"] = STATS["fail"].get("low", 0) + 1
            continue
        yy1, yy2 = low.mult(y, y1), low.mult(y, y2)
        vlist = [zp, y, y1, y2, yy1, yy2]
        flats = [sum(local(v) << (gi * a) for gi in gis) for v in vlist]
        if _rank(flats) < 6:
            STATS["fail"]["flats"] = STATS["fail"].get("flats", 0) + 1
            continue
        betas = [coordsY(f) for f in flats]
        if any(bt is None for bt in betas):
            STATS["fail"]["notinW"] = STATS["fail"].get("notinW", 0) + 1
            continue
        forms = [[coordsX_form(betas[j]) for j in range(6)]]   # forms[0][j]: u -> beta_j . alpha(u)
        F = forms[0]
        colcons = [(z, 0) for z in lker]

        def core_cons(k):                       # beta_j . alpha(u) = delta_jk
            return [(F[j], 1 if j == k else 0) for j in range(6)]

        # ---- high side: A (term 0: A (x) z'), then A1 (term 4: A1 (x) y y1; term 2: A A1 (x) y1),
        #      A2 (term 5; term 3: A A2 (x) y2), B' (term 1: B' (x) y)
        A = _solve_affine(colcons + [(hv, 0)] + core_cons(0), nr, rng)
        if not A:
            STATS["fail"]["A"] = STATS["fail"].get("A", 0) + 1
            continue
        MA = mult_matrix(A)
        lk_pull = [(pull(z, MA), 0) for z in lker]
        A1 = _solve_affine(colcons + [(hv, 0)] + core_cons(4) + lk_pull +
                           [(pull(F[j], MA), 1 if j == 2 else 0) for j in range(6)], nr, rng)
        A2 = _solve_affine(colcons + [(hv, 0)] + core_cons(5) + lk_pull +
                           [(pull(F[j], MA), 1 if j == 3 else 0) for j in range(6)], nr, rng)
        if not (A1 and A2):
            STATS["fail"]["A12"] = STATS["fail"].get("A12", 0) + 1
            continue
        TA, TA1, TA2 = h_truth(A), h_truth(A1), h_truth(A2)
        c12, k12 = h_coef(TA1 & TA2)
        assert k12 == 0
        # B' : in colspace, core term 1, B = B' + A1 A2 in H, sigma = V-sum(T_A * T_B') = 0
        sig = 0
        for r in range(nr):
            if bin(TA & HTr[r] & Vset).count("1") & 1:
                sig |= 1 << r
        vs12 = bin(TA1 & TA2 & Vset).count("1") & 1
        Bp = _solve_affine(colcons + core_cons(1) + [(hv, vs12), (sig, 0)], nr, rng)
        if Bp is None:
            STATS["fail"]["B"] = STATS["fail"].get("B", 0) + 1
            continue
        Bc = Bp ^ c12                             # B = B' + A1 A2 (as functions; no constants anywhere)
        z = zp ^ low.mult(y1, y2)                 # z = z' + y1 y2
        # ---- exact check 1: P1's mixed part from its truth table over HB x 16 addresses
        Th = {k: h_truth(v) for k, v in (("A", A), ("B", Bc), ("A1", A1), ("A2", A2))}
        Tl = {k: low.truth(v) for k, v in (("y", y), ("z", z), ("y1", y1), ("y2", y2))}
        Fh = np.zeros((HB, 16), np.uint8)
        for h in range(HB):
            for l in range(16):
                def hv_(k):
                    return (Th[k] >> h) & 1

                def lv_(k):
                    return (Tl[k] >> l) & 1
                p0 = (hv_("A1") ^ lv_("y2")) & (hv_("A2") ^ lv_("y1"))
                Fh[h, l] = (hv_("A") ^ lv_("y")) & (hv_("B") ^ lv_("z") ^ p0)
        mrows = mixed_rows(an, low, Fh)
        # the same from the formula: sum_k u_k (x) v_k
        TAA1, _ = h_coef(TA & TA1)
        TAA2, _ = h_coef(TA & TA2)
        us = [A, Bp, TAA1, TAA2, A1, A2]
        frows = [0] * nr
        for uk, vk in zip(us, vlist):
            for r in bits(uk):
                frows[r] ^= vk
        if mrows != frows:
            STATS["fail"]["formula"] = STATS["fail"].get("formula", 0) + 1
            continue
        # ---- exact check 2: rank(R - C) = rank(R) - 6
        newrows = []
        for r in range(nr):
            v = rows[r]
            for gi in gis:
                v ^= local(mrows[r]) << (gi * a)
            newrows.append(v)
        rank1 = _rank(newrows)
        if rank1 != rank0 - 6:
            STATS["fail"]["rank1"] = STATS["fail"].get("rank1", 0) + 1
            continue
        if any(mrows[r] & nonlive for r in range(nr)):
            STATS["fail"]["nonlive"] = STATS["fail"].get("nonlive", 0) + 1
            continue
        STATS["built"] += 1
        sglob = sum(1 << outs[gi] for gi in gis)
        return dict(A=A, B=Bc, A1=A1, A2=A2, y=y, z=z, y1=y1, y2=y2, s=sglob, rows=mrows,
                    rank0=rank0, rank1=rank1, attempt=attempt)
    return None


def apply_block(dec, blk):
    """Subtract the block's mixed part from the group's outputs in the pattern."""
    rowsarr = np.array(blk["rows"], np.uint32)
    for o in bits(blk["s"]):
        dec.X[:, o] ^= rowsarr


def install(target_gi):
    """Wrap core.Decomp.apply_kills: after the kills of generic group target_gi, build and subtract one block."""
    inner = core.Decomp.apply_kills

    def wrapped(self, outs, coords, dl, d, shift=0):
        live = inner(self, outs, coords, dl, d, shift)
        if live is None or getattr(self, "_gi", -1) != target_gi:
            return live
        self.cubic_blocks = []
        for _ in range(int(os.environ.get("C4_CUBIC_N", "1"))):     # blocks one after the other, each rank-additive
            blk = make_block(self, outs, live, self.rng, int(os.environ.get("C4_CUBIC_TRIES", "300")))
            if not blk:
                break
            apply_block(self, blk)
            self.cubic_blocks.append(blk)
        return live
    core.Decomp.apply_kills = wrapped


# --------------------------------------------------------------------------------------- emission (C2's emit_c4)
def emit_c4x(an, prods, V, blocks):
    """kmx text and stats: chem_ld2_c4.emit_c4 (builder C2, read-only; copied here) with the cubic blocks added."""
    bits_, BROT, J, QUBIT_CAP = core.bits, core.BROT, core.J, core.QUBIT_CAP
    A, lb, L, nl, n_sel = an.A, an.lb, an.L, an.nl, an.n_sel
    assert an.H == c4.HB and an.hbits == 7
    aa = np.arange(A)
    lidx, hidx = aa & (L - 1), aa >> lb
    ALL = (1 << A) - 1
    LBf = []
    for bs in an.basis:
        tt = np.zeros(L, np.uint8)
        for m in bits_(bs["anf"]):
            tt[m] = 1
        tt = core.mobius(tt, lb)
        LBf.append(core.to_int(tt[lidx]))

    def addr_truth(g):
        bv = np.array([(g >> p) & 1 for p in range(c4.HB)], np.uint8)
        return core.to_int(bv[hidx])
    HMf = [None] + [core.to_int((hidx & M) == M) for M in range(1, an.H)]

    ops_c, truths = c4.build_hcache(V)
    nwc = len(truths)
    CT = c4.TSpan()
    CT.add((1 << c4.HB) - 1, 1 << nwc)
    for w, f in enumerate(truths):
        CT.add(f, 1 << w)

    def hwires(mask):
        r, tag = CT.reduce(c4.hmask_truth(mask))
        assert r == 0, "a product high part is outside H"
        return tag

    def wire_q(w):
        return n_sel + BROT * w[0] + w[1]

    ops = []
    free = sorted(wire_q(c["wire"]) for c in an.carriers if c["clean"]) + list(range(n_sel + J * BROT, QUBIT_CAP))
    used_max = n_sel + J * BROT - 1

    def take():
        nonlocal used_max
        if not free:
            raise ValueError("out of qubits")
        q = free.pop(0)
        used_max = max(used_max, q)
        return q

    def give(q):
        free.append(q)
        free.sort()

    lowq = [None] * nl
    for t, bs in enumerate(an.basis):
        if bs["kind"] == "lin":
            lowq[t] = bs["bit"]
        elif bs["kind"] == "car":
            lowq[t] = wire_q(bs["wire"])
    prep = []
    for c in an.carriers:
        q = wire_q(c["wire"])
        if c["conj"]:
            prep.append(f"CX q{q - 1} q{q}")
        if c["const"]:
            prep.append(f"X q{q}")
        if c["clean"]:
            for t in c["combo"]:
                prep.append(f"CX q{lowq[t]} q{q}")
    ops += prep

    def controls(s1, s2, const1=False, const2=False):
        s1, s2 = set(s1), set(s2)
        assert s1 and s2 and s1 != s2, "degenerate controls"
        p1 = min(s1 - s2) if s1 - s2 else min(s1)
        pre = [f"CX q{q} q{p1}" for q in sorted(s1) if q != p1]
        if p1 in s2:
            s2 = (s1 ^ s2) | {p1}
            const2 = const2 ^ const1
        if const1:
            pre.append(f"X q{p1}")
        p2 = min(q for q in s2 if q != p1)
        pre += [f"CX q{q} q{p2}" for q in sorted(s2) if q != p2]
        if const2:
            pre.append(f"X q{p2}")
        return pre, p1, p2

    cq = [lb + i for i in range(7)]
    made_high = []
    for (t1, t2) in ops_c:
        pre, p1, p2 = controls([cq[w] for w in t1], [cq[w] for w in t2])
        t = take()
        ops.extend(pre)
        ops.append(f"CCX q{p1} q{p2} q{t}")
        ops.extend(reversed(pre))
        cq.append(t)
        made_high.append((pre, p1, p2, t))

    def low_and(t):
        par, i = an.recipe[t]
        qs = [lowq[x] for x in par]
        piv = next(q for q in qs if q != i)
        pre = [f"CX q{q} q{piv}" for q in qs if q != piv]
        tq = take()
        ops.extend(pre)
        ops.append(f"CCX q{piv} q{i} q{tq}")
        ops.extend(reversed(pre))
        lowq[t] = tq
        return (pre, piv, i, tq, t)

    def low_del(rec):
        pre, piv, i, tq, t = rec
        ops.extend(pre)
        ops.append(f"HMR q{tq} b0")
        ops.append(f"CZ q{piv} q{i} if b0")
        ops.extend(reversed(pre))
        lowq[t] = None
        give(tq)

    nw = len(an.out_wires)
    wq = [wire_q(w) for w in an.out_wires]

    def func(low_mask, high_mask):
        f = 0
        for t in bits_(low_mask):
            f ^= LBf[t]
        for M in bits_(high_mask):
            f ^= HMf[M + 1]
        return f

    acc = [0] * nw
    targets_of = []
    for a1, b1, a2, b2, s, _ in prods:
        P = func(a1, b2) & func(a2, b1)
        tg = [q for q in range(nw) if (an.G[q] & s).bit_count() & 1]
        assert tg, "product without a target"
        targets_of.append(tg)
        for q in tg:
            acc[q] ^= P
    btg = []
    for blk in blocks:                           # holder 4: the cubic blocks' P1 into their targets
        T0 = func(blk["y2"], blk["A1"]) & func(blk["y1"], blk["A2"])
        T1 = func(blk["y"], blk["A"]) & (func(blk["z"], blk["B"]) ^ T0)
        tg = [q for q in range(nw) if (an.G[q] & blk["s"]).bit_count() & 1]
        assert tg, "block without a target"
        btg.append(tg)
        for q in tg:
            acc[q] ^= T1
    pure = core.Span()
    pure.add(ALL, 1)
    for t in range(nl):
        pure.add(LBf[t], 1 << (1 + t))
    for w, f in enumerate(truths):
        pure.add(addr_truth(f), 1 << (1 + nl + w))
    corr = []
    for q in range(nw):
        r, c = pure.reduce(an.delta_int[q] ^ acc[q])
        assert r == 0, f"wire {an.out_wires[q]}: pure-high part outside H (tau != kappa) on circuit {an.k}"
        corr.append(c)

    core_made = [low_and(t) for t in an.core]
    by_phase = [[] for _ in an.excl]
    for i, p in enumerate(prods):
        by_phase[p[5]].append(i)
    n_emitted = 0
    mask_w = (1 << nwc) - 1

    def do_product(i):
        nonlocal n_emitted
        a1, b1, a2, b2, s, _ = prods[i]
        h2, h1 = hwires(b2), hwires(b1)
        s1 = [lowq[t] for t in bits_(a1)] + [cq[w] for w in bits_(h2 & mask_w)]
        s2 = [lowq[t] for t in bits_(a2)] + [cq[w] for w in bits_(h1 & mask_w)]
        assert None not in s1 and None not in s2, "product uses a low function that is not cached in its phase"
        pre, p1, p2 = controls(s1, s2, bool(h2 >> nwc & 1), bool(h1 >> nwc & 1))
        tg = targets_of[i]
        first = wq[tg[0]]
        fan = [f"CX q{first} q{wq[q]}" for q in tg[1:]]
        ops.extend(pre)
        ops.extend(fan)
        ops.append(f"CCX q{p1} q{p2} q{first}")
        ops.extend(reversed(fan))
        ops.extend(reversed(pre))
        n_emitted += 1

    def do_block(j):
        """P0 on a fresh wire; P1 = (A + y)(B + z + P0) into the targets; P0 erased by HMR + CZ."""
        nonlocal n_emitted
        blk = blocks[j]
        hA1, hA2, hA, hB = hwires(blk["A1"]), hwires(blk["A2"]), hwires(blk["A"]), hwires(blk["B"])
        s1 = [lowq[t] for t in bits_(blk["y2"])] + [cq[w] for w in bits_(hA1 & mask_w)]
        s2 = [lowq[t] for t in bits_(blk["y1"])] + [cq[w] for w in bits_(hA2 & mask_w)]
        assert None not in s1 and None not in s2
        pre0, p1, p2 = controls(s1, s2, bool(hA1 >> nwc & 1), bool(hA2 >> nwc & 1))
        t0 = take()
        ops.extend(pre0)
        ops.append(f"CCX q{p1} q{p2} q{t0}")
        ops.extend(reversed(pre0))
        u1 = [lowq[t] for t in bits_(blk["y"])] + [cq[w] for w in bits_(hA & mask_w)]
        u2 = [lowq[t] for t in bits_(blk["z"])] + [cq[w] for w in bits_(hB & mask_w)] + [t0]
        assert None not in u1 and None not in u2
        pre1, q1, q2 = controls(u1, u2, bool(hA >> nwc & 1), bool(hB >> nwc & 1))
        tg = btg[j]
        first = wq[tg[0]]
        fan = [f"CX q{first} q{wq[q]}" for q in tg[1:]]
        ops.extend(pre1)
        ops.extend(fan)
        ops.append(f"CCX q{q1} q{q2} q{first}")
        ops.extend(reversed(fan))
        ops.extend(reversed(pre1))
        ops.extend(pre0)
        ops.append(f"HMR q{t0} b0")
        ops.append(f"CZ q{p1} q{p2} if b0")
        ops.extend(reversed(pre0))
        give(t0)
        n_emitted += 2

    def low_corrections(ts):
        for q in range(nw):
            for t in ts:
                if corr[q] >> (1 + t) & 1:
                    ops.append(f"CX q{lowq[t]} q{wq[q]}")

    excl_all = {t for ex in an.excl for t in ex}
    for ph, ex in enumerate(an.excl):
        made = [low_and(t) for t in ex]
        for i in by_phase[ph]:
            do_product(i)
        if ph == len(an.excl) - 1:
            for j in range(len(blocks)):
                do_block(j)
        low_corrections(ex)
        if ph == len(an.excl) - 1:
            low_corrections([t for t in range(nl) if t not in excl_all])
            for q in range(nw):
                c = corr[q]
                if c & 1:
                    ops.append(f"X q{wq[q]}")
                for w in range(nwc):
                    if c >> (1 + nl + w) & 1:
                        ops.append(f"CX q{cq[w]} q{wq[q]}")
        for rec in reversed(made):
            low_del(rec)
    for rec in reversed(core_made):
        low_del(rec)
    for pre, p1, p2, t in reversed(made_high):
        ops.extend(pre)
        ops.append(f"HMR q{t} b0")
        ops.append(f"CZ q{p1} q{p2} if b0")
        ops.extend(reversed(pre))
    ops.extend(reversed(prep))

    head = []
    for i in range(n_sel):
        head.append(f"APPEND_TO_REGISTER q{i} r0")
    head.append('REGISTER r0 "selection"')
    for j in range(J):
        for b in range(BROT):
            head.append(f"APPEND_TO_REGISTER q{n_sel + BROT * j + b} r{j + 1}")
        head.append(f'REGISTER r{j + 1} "word_{j}"')
    infra = len(ops_c) + len(an.comp)
    stats = dict(k=an.k, infra=infra, products=n_emitted, toffoli=infra + n_emitted, qubits=used_max + 1,
                 phases=len(an.excl), outputs=len(an.out_level), ops=len(ops), cache="c4-68", blocks=len(blocks))
    return "\n".join(head + ops) + "\n", stats
