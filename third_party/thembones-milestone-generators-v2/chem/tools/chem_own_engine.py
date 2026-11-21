"""chem_own_engine: product plans for one ANGLE_k lookup circuit, priced exactly by linear matroid parity.

The circuit scheme (address split, caches, carriers, emission) is the one of tools/chem_loaders_core.py (ours,
imported unchanged).  This file replaces its pairing step.

Setting (see research/ws_chem_own.md).  The low coordinates are cut into blocks; block b has a cell for every
(low coordinate, output) pair, its table rows span W_b, and a plan keeps a relation space R_b inside the
annihilator of W_b.  A feature is (point, s): a nonzero vector y of one block and a fan-out s with y (x) s in
F_b = R_b-perp; two features with the same fan-out share one Toffoli.  With every feature paired the product
count is (cells - sum_b dim R_b) / 2, so the whole game is to keep as many relations as possible while a
perfect pairing still exists.

What is implemented here, and where each idea comes from:
  * exact pairing count (ours, after the public description "plans priced by matroid parity" in marculera's
    notebook "DF-THC 11.2x: Codes on the XOR of Five Words"; the notebook gives no code): the maximum number
    of pairs is rank(Y) / 2 for Y = sum over candidate lines of x_l (u v^T + v u^T) with random x_l in
    GF(2^12) (Lovasz 1979; one-sided error, a reported rank is always achievable);
  * the search "drop spare classes, then swap at a fixed class count until the pairing is perfect again"
    (same notebook, finding 2), written here on the relation spans (Engine.level_up): add one relation
    (the best of all single relations by exact count), then repair with two kinds of single exchanges,
    Engine.trial (open a direction, close with the best relation of any block) and Engine.swap2 (add a
    relation, open the best direction of any block), until every feature is paired again;
  * the target "five blocks of rank 105 - 19" (maoxin123's notebook "LB 25.303 | Combined invsqrt + refit
    chem": 217 -> 215 products) is only used as a yardstick: the search is free to reach any codimension
    (it reaches 94 = 216 products on the public words; 95 was not reached);
  * ours, not in the notebooks: all single-relation additions of a block priced from ONE elimination
    (bordered-matrix identity rank(Y - K) = rho + rank(C_K) - 2|K|, _transform / _scan_block / _scan_exact);
    candidate closings ranked by the features that survive them, for all 2^m candidates with one
    Walsh-Hadamard transform (_block_caps); all 2^e openings of a block ranked by the point pairs that gain
    a common fan-out, also one transform (_open_hist); the don't-care addresses of the last partial address
    block used to shrink every block's row space (relax_last_row);
  * extraction of the pairing itself (ours, after Cheung, Kwok, Lau 2014): keep Y^-1 and delete every line
    whose removal leaves Y nonsingular (rank-2 Woodbury updates); what is left is a parity basis, which is
    then checked exactly over GF(2) before emission (unpaired features go through dummy coordinates).
Tried and switched off by default (measured, research/ws_chem_own.md): the pair-obstruction filter
(pair_filter) and the equal-price plateau walk (walk in level_up).
No plan, span or fan-out is read from any public circuit: everything is computed from the table.
"""
from __future__ import annotations

import os
import random
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
from chem_loaders_core import Span, bits  # noqa: E402

EMAX = 48          # largest relation space of a block that is stored
CMAX = 4           # largest block (low coordinates)
GF_BITS = 12
GF_POLY = 0x1053   # x^12 + x^6 + x^4 + x + 1, primitive (checked at import)
MAXLINES = 60000


def _gf_tables(nbits: int, poly: int):
    """exp (length 4(q-1): three periods of the powers of the generator, then q-1 zeros) and log.
    LOGZ = 3(q-1) is used as the "logarithm of zero": exp[LOGZ + t] = 0 for 0 <= t < q-1."""
    q = 1 << nbits
    exp = np.zeros(4 * (q - 1), np.uint16)
    log = np.zeros(q, np.int64)
    x = 1
    seen = set()
    for i in range(q - 1):
        exp[i] = x
        log[x] = i
        seen.add(x)
        x <<= 1
        if x & q:
            x ^= poly
    assert x == 1 and len(seen) == q - 1, "GF polynomial is not primitive"
    exp[q - 1:2 * (q - 1)] = exp[:q - 1]
    exp[2 * (q - 1):3 * (q - 1)] = exp[:q - 1]
    return exp, log


GF_EXP, GF_LOG = _gf_tables(GF_BITS, GF_POLY)


# ------------------------------------------------------------------------------------------ numba kernels
@njit(cache=True)
def _nullspace(cons, ncons, support, nb, out):
    """Basis of {s inside support : parity(s & c) = 0 for all constraints c}; returns its dimension."""
    prow = np.zeros(64, np.int64)
    pbit = np.zeros(64, np.int64)
    ispiv = np.zeros(64, np.bool_)
    npiv = 0
    for j in range(ncons):
        v = cons[j] & support
        for q in range(npiv):
            if (v >> pbit[q]) & 1:
                v ^= prow[q]
        if v != 0:
            pb = 0
            t = v
            while t > 1:
                t >>= 1
                pb += 1
            for q in range(npiv):
                if (prow[q] >> pb) & 1:
                    prow[q] ^= v
            prow[npiv] = v
            pbit[npiv] = pb
            ispiv[pb] = True
            npiv += 1
    cnt = 0
    for f in range(nb):
        if (support >> f) & 1 and not ispiv[f]:
            x = np.int64(1) << f
            for q in range(npiv):
                if (prow[q] >> f) & 1:
                    x |= np.int64(1) << pbit[q]
            out[cnt] = x
            cnt += 1
    return cnt


@njit(cache=True)
def _point_cons(pt_blk, pt_ymask, pt_outs, ncoord, Rm, ecount):
    npnt = len(pt_blk)
    cons = np.zeros((npnt, EMAX), np.int64)
    for i in range(npnt):
        b = pt_blk[i]
        for j in range(ecount[b]):
            c = np.int64(0)
            for ci in range(ncoord[b]):
                if (pt_ymask[i] >> ci) & 1:
                    c ^= Rm[b, j, ci]
            cons[i, j] = c & pt_outs[i]
    return cons


@njit(cache=True)
def _lines(pt_blk, pt_ymask, pt_outs, pt_phase, ncoord, Rm, ecount, nb, li, lj, ls, deg):
    """All candidate lines (i, j, s): s a nonzero fan-out valid at points i and j. Returns the count."""
    npnt = len(pt_blk)
    cons = _point_cons(pt_blk, pt_ymask, pt_outs, ncoord, Rm, ecount)
    both = np.zeros(2 * EMAX, np.int64)
    basis = np.zeros(64, np.int64)
    maxl = len(li)
    n = 0
    for i in range(npnt):
        deg[i] = 0
    for i in range(npnt):
        ei = ecount[pt_blk[i]]
        for j in range(i + 1, npnt):
            if (pt_phase[i] & pt_phase[j]) == 0:
                continue
            sup = pt_outs[i] & pt_outs[j]
            if sup == 0:
                continue
            ej = ecount[pt_blk[j]]
            for q in range(ei):
                both[q] = cons[i, q]
            for q in range(ej):
                both[ei + q] = cons[j, q]
            d = _nullspace(both, ei + ej, sup, nb, basis)
            if d == 0:
                continue
            deg[i] += d
            deg[j] += d
            if d <= 4:
                for m in range(1, 1 << d):
                    s = np.int64(0)
                    for q in range(d):
                        if (m >> q) & 1:
                            s ^= basis[q]
                    if n < maxl:
                        li[n] = i
                        lj[n] = j
                        ls[n] = s
                        n += 1
            else:
                # the basis, all sums of two, and sums of three consecutive basis vectors
                for q in range(d):
                    if n < maxl:
                        li[n] = i
                        lj[n] = j
                        ls[n] = basis[q]
                        n += 1
                for q in range(d):
                    for q2 in range(q + 1, d):
                        if n < maxl:
                            li[n] = i
                            lj[n] = j
                            ls[n] = basis[q] ^ basis[q2]
                            n += 1
                for q in range(d - 2):
                    if n < maxl:
                        li[n] = i
                        lj[n] = j
                        ls[n] = basis[q] ^ basis[q + 1] ^ basis[q + 2]
                        n += 1
    return n


@njit(cache=True)
def _singles(pt_blk, pt_ymask, pt_outs, ncoord, Rm, ecount, nb, ndummy, start, li, lj, ls):
    """Lines (feature, dummy coordinate) for a basis of the fan-outs of every point; appended at `start`."""
    npnt = len(pt_blk)
    cons = _point_cons(pt_blk, pt_ymask, pt_outs, ncoord, Rm, ecount)
    basis = np.zeros(64, np.int64)
    n = start
    for i in range(npnt):
        d = _nullspace(cons[i], ecount[pt_blk[i]], pt_outs[i], nb, basis)
        for q in range(d):
            for k in range(ndummy):
                if n < len(li):
                    li[n] = i
                    lj[n] = -1 - k
                    ls[n] = basis[q]
                    n += 1
    return n


@njit(cache=True)
def _feat(i, s, pt_blk, pt_ymask, ncoord, cellid, freeidx, nb, out):
    b = pt_blk[i]
    ym = pt_ymask[i]
    n = 0
    for ci in range(ncoord[b]):
        if (ym >> ci) & 1:
            for o in range(nb):
                if (s >> o) & 1:
                    f = freeidx[cellid[b, ci, o]]
                    if f >= 0:
                        out[n] = f
                        n += 1
    return n


@njit(cache=True)
def _build_y(n, li, lj, ls, lx, pt_blk, pt_ymask, ncoord, cellid, freeidx, nb, r, size):
    Y = np.zeros((size, size), np.uint16)
    i1 = np.zeros(CMAX * 64, np.int64)
    i2 = np.zeros(CMAX * 64, np.int64)
    for ell in range(n):
        n1 = _feat(li[ell], ls[ell], pt_blk, pt_ymask, ncoord, cellid, freeidx, nb, i1)
        if lj[ell] >= 0:
            n2 = _feat(lj[ell], ls[ell], pt_blk, pt_ymask, ncoord, cellid, freeidx, nb, i2)
        else:
            i2[0] = r + (-1 - lj[ell])
            n2 = 1
        x = lx[ell]
        for a in range(n1):
            ra = i1[a]
            for b in range(n2):
                cb = i2[b]
                Y[ra, cb] ^= x
                Y[cb, ra] ^= x
    return Y


@njit(cache=True)
def gf_rank(Y, exp, log):
    """Rank over GF(2^k) (Y is destroyed)."""
    n = Y.shape[0]
    qm1 = len(exp) // 4
    logz = 3 * qm1
    lg = np.empty(n, np.int64)
    row = 0
    for col in range(n):
        piv = -1
        for i in range(row, n):
            if Y[i, col] != 0:
                piv = i
                break
        if piv < 0:
            continue
        if piv != row:
            for j in range(col, n):
                t = Y[row, j]
                Y[row, j] = Y[piv, j]
                Y[piv, j] = t
        lp = log[Y[row, col]]
        for j in range(col + 1, n):
            v = Y[row, j]
            if v != 0:
                t = log[v] - lp
                if t < 0:
                    t += qm1
                lg[j] = t
            else:
                lg[j] = logz
        for i in range(row + 1, n):
            c = Y[i, col]
            if c != 0:
                lc = log[c]
                for j in range(col + 1, n):
                    Y[i, j] ^= exp[lc + lg[j]]
        row += 1
        if row == n:
            break
    return row


@njit(cache=True)
def gf_inverse(Y, exp, log):
    """(ok, inverse) by Gauss-Jordan over GF(2^k)."""
    n = Y.shape[0]
    qm1 = len(exp) // 4
    A = np.zeros((n, 2 * n), np.uint16)
    for i in range(n):
        for j in range(n):
            A[i, j] = Y[i, j]
        A[i, n + i] = 1
    cols = np.empty(2 * n, np.int64)
    lgs = np.empty(2 * n, np.int64)
    for col in range(n):
        piv = -1
        for i in range(col, n):
            if A[i, col] != 0:
                piv = i
                break
        if piv < 0:
            return False, A[:, n:].copy()
        if piv != col:
            for j in range(2 * n):
                t = A[col, j]
                A[col, j] = A[piv, j]
                A[piv, j] = t
        lp = log[A[col, col]]
        cnt = 0
        for j in range(col, 2 * n):
            v = A[col, j]
            if v != 0:
                v2 = exp[log[v] + qm1 - lp]
                A[col, j] = v2
                if j != col:
                    cols[cnt] = j
                    lgs[cnt] = log[v2]
                    cnt += 1
        for i in range(n):
            if i != col:
                c = A[i, col]
                if c != 0:
                    lc = log[c]
                    for q in range(cnt):
                        A[i, cols[q]] ^= exp[lc + lgs[q]]
                    A[i, col] = 0
    return True, A[:, n:].copy()


@njit(cache=True)
def _extract(N, n, li, lj, ls, lx, order, pt_blk, pt_ymask, ncoord, cellid, freeidx, nb, r, exp, log):
    """Delete every line whose removal keeps Y nonsingular (N = Y^-1 is updated in place).
    Returns keep flags; the kept lines are a parity basis with high probability."""
    size = N.shape[0]
    qm1 = len(exp) // 4
    keep = np.ones(n, np.bool_)
    i1 = np.zeros(CMAX * 64, np.int64)
    i2 = np.zeros(CMAX * 64, np.int64)
    w1 = np.zeros(size, np.uint16)
    w2 = np.zeros(size, np.uint16)
    l1 = np.zeros(size, np.int64)
    l2 = np.zeros(size, np.int64)
    for t in range(n):
        ell = order[t]
        n1 = _feat(li[ell], ls[ell], pt_blk, pt_ymask, ncoord, cellid, freeidx, nb, i1)
        if lj[ell] >= 0:
            n2 = _feat(lj[ell], ls[ell], pt_blk, pt_ymask, ncoord, cellid, freeidx, nb, i2)
        else:
            i2[0] = r + (-1 - lj[ell])
            n2 = 1
        for u in range(size):
            w2[u] = 0
        for b in range(n2):
            row = i2[b]
            for u in range(size):
                w2[u] ^= N[row, u]
        a = 0
        for q in range(n1):
            a ^= w2[i1[q]]
        x = lx[ell]
        xa = 0
        if a != 0:
            xa = exp[log[x] + log[a]]
        if xa == 1:
            continue                       # removing this line would make Y singular: it stays
        for u in range(size):
            w1[u] = 0
        for b in range(n1):
            row = i1[b]
            for u in range(size):
                w1[u] ^= N[row, u]
        den = 1 ^ xa
        lc = (log[x] + qm1 - log[den]) % qm1       # log of c = x / (1 + x a)
        for u in range(size):
            l1[u] = log[w1[u]] if w1[u] != 0 else -1
            l2[u] = log[w2[u]] if w2[u] != 0 else -1
        for u in range(size):
            a1 = l1[u]
            a2 = l2[u]
            if a1 < 0 and a2 < 0:
                continue
            for v in range(size):
                acc = 0
                if a1 >= 0 and l2[v] >= 0:
                    acc ^= exp[lc + a1 + l2[v]]
                if a2 >= 0 and l1[v] >= 0:
                    acc ^= exp[lc + a2 + l1[v]]
                N[u, v] ^= acc
        keep[ell] = False
    return keep


@njit(cache=True)
def _parity(x):
    x ^= x >> 32
    x ^= x >> 16
    x ^= x >> 8
    x ^= x >> 4
    x ^= x >> 2
    x ^= x >> 1
    return x & 1


@njit(cache=True)
def _transform(Y, exp, log):
    """Reduced row echelon form of [Y | I] over GF(2^k).
    Returns rho = rank Y, ET (r x r, ET[i] = column i of the row transform E, E Y = [[R], [0]]),
    Rnp (rho x d, the non-pivot columns of R), pivrow (coordinate -> pivot row or -1), npidx
    (coordinate -> index among the non-pivot coordinates or -1)."""
    r = Y.shape[0]
    w = 2 * r
    qm1 = len(exp) // 4
    A = np.zeros((r, w), np.uint16)
    for i in range(r):
        for j in range(r):
            A[i, j] = Y[i, j]
        A[i, r + i] = 1
    pivrow = np.full(r, -1, np.int64)
    cols = np.empty(w, np.int64)
    lgs = np.empty(w, np.int64)
    row = 0
    for col in range(r):
        piv = -1
        for i in range(row, r):
            if A[i, col] != 0:
                piv = i
                break
        if piv < 0:
            continue
        if piv != row:
            for j in range(col, w):
                t = A[row, j]
                A[row, j] = A[piv, j]
                A[piv, j] = t
        lp = log[A[row, col]]
        cnt = 0
        for j in range(col, w):
            v = A[row, j]
            if v != 0:
                v2 = exp[log[v] + qm1 - lp]
                A[row, j] = v2
                if j != col:
                    cols[cnt] = j
                    lgs[cnt] = log[v2]
                    cnt += 1
        for i in range(r):
            if i != row:
                c = A[i, col]
                if c != 0:
                    lc = log[c]
                    for q in range(cnt):
                        A[i, cols[q]] ^= exp[lc + lgs[q]]
                    A[i, col] = 0
        pivrow[col] = row
        row += 1
    rho = row
    d = r - rho
    npidx = np.full(r, -1, np.int64)
    Rnp = np.zeros((rho, d), np.uint16)
    t = 0
    for col in range(r):
        if pivrow[col] < 0:
            npidx[col] = t
            for i in range(rho):
                Rnp[i, t] = A[i, col]
            t += 1
    ET = np.zeros((r, r), np.uint16)
    for i in range(r):
        for j in range(r):
            ET[j, i] = A[i, r + j]
    return rho, ET, Rnp, pivrow, npidx


@njit(cache=True)
def _scan_block(rho, ET, Rnp, pivrow, npidx, sel, li, lj, ls, lx, pt_blk, pt_ymask, ncoord, cellid, freeidx, nb,
                blk, Qm, nq, exp, log):
    """Pieces of the Schur complement for the lines `sel` (those touching block blk); see Engine.scan.

    With U the feature vectors of those lines (two columns per line) and D = diag blocks (1/x)[[0,1],[1,0]],
    rank(Y - lines K) = rank [[Y, U_K], [U_K^T, D_K]] - 2|K| = rho + rank(C_K) - 2|K| with
    C_K = [[0, Zr_K], [Zc_K, D_K + T_KK]] (T: 2n x 2n, Zr: d x 2n, Zc: 2n x d, d = r - rho).
    h1/h2: for each line, bit t = <Q_t, its feature in block blk> (0 when that feature is in another block)."""
    r = ET.shape[0]
    n = len(sel)
    d = r - rho
    qm1 = len(exp) // 4
    i1 = np.zeros(CMAX * 64, np.int64)
    h1 = np.zeros(n, np.int64)
    h2 = np.zeros(n, np.int64)
    dinv = np.zeros(n, np.uint16)
    PT = np.zeros((2 * n, r), np.uint16)          # PT[c] = E u_c
    T = np.zeros((2 * n, 2 * n), np.uint16)
    Zc = np.zeros((2 * n, d), np.uint16)
    Zr = np.zeros((d, 2 * n), np.uint16)
    sup_start = np.zeros(2 * n + 1, np.int64)
    sup = np.zeros(2 * n * CMAX * 64, np.int64)
    for x in range(n):
        ell = sel[x]
        for side in range(2):
            pt = li[ell] if side == 0 else lj[ell]
            cnt = _feat(pt, ls[ell], pt_blk, pt_ymask, ncoord, cellid, freeidx, nb, i1)
            c = 2 * x + side
            sup_start[c + 1] = sup_start[c] + cnt
            for q in range(cnt):
                coord = i1[q]
                sup[sup_start[c] + q] = coord
                for u in range(r):
                    PT[c, u] ^= ET[coord, u]
                pr = pivrow[coord]
                if pr >= 0:
                    for t in range(d):
                        Zc[c, t] ^= Rnp[pr, t]
                else:
                    Zc[c, npidx[coord]] ^= 1
            if pt_blk[pt] == blk:
                hm = np.int64(0)
                for t in range(nq):
                    cv = np.int64(0)
                    for ci in range(ncoord[blk]):
                        if (pt_ymask[pt] >> ci) & 1:
                            cv ^= Qm[t, ci]
                    if _parity(cv & ls[ell]):
                        hm |= np.int64(1) << t
                if side == 0:
                    h1[x] = hm
                else:
                    h2[x] = hm
        dinv[x] = exp[qm1 - log[lx[ell]]]
    for c in range(2 * n):
        for i in range(d):
            Zr[i, c] = PT[c, rho + i]
    for c in range(2 * n):
        for q in range(sup_start[c], sup_start[c + 1]):
            pr = pivrow[sup[q]]
            if pr >= 0:
                for c2 in range(2 * n):
                    T[c, c2] ^= PT[c2, pr]
    return T, Zr, Zc, dinv, h1, h2


@njit(cache=True)
def _essential(T, Zr, Zc, dinv):
    """Flags of the lines that every maximum pairing contains (their own Schur block C_l vanishes)."""
    n = len(dinv)
    d = Zr.shape[0]
    out = np.zeros(n, np.bool_)
    for x in range(n):
        ok = True
        for side in range(2):
            c = 2 * x + side
            for i in range(d):
                if Zr[i, c] != 0 or Zc[c, i] != 0:
                    ok = False
        if ok:
            if T[2 * x, 2 * x] != 0 or T[2 * x + 1, 2 * x + 1] != 0:
                ok = False
            if T[2 * x, 2 * x + 1] != dinv[x] or T[2 * x + 1, 2 * x] != dinv[x]:
                ok = False
        out[x] = ok
    return out


@njit(cache=True)
def _pair_obstructions(T, Zr, Zc, dinv, ess, exp, log):
    """obs[x] = set (62-bit words) of the lines x' that cannot be removed together with x without losing a
    pair: rank of the Schur block C_{x,x'} below 4. Essential lines are skipped (they can never be removed)."""
    n = len(dinv)
    d = Zr.shape[0]
    nw = (n + 61) // 62
    obs = np.zeros((n, nw), np.int64)
    size = d + 4
    C = np.zeros((size, size), np.uint16)
    idx = np.zeros(2, np.int64)
    for x in range(n):
        if ess[x]:
            continue
        for x2 in range(x + 1, n):
            if ess[x2]:
                continue
            idx[0] = x
            idx[1] = x2
            for i in range(size):
                for j in range(size):
                    C[i, j] = 0
            for u in range(2):
                xu = idx[u]
                for su in range(2):
                    cu = 2 * xu + su
                    ru = d + 2 * u + su
                    for i in range(d):
                        C[i, ru] = Zr[i, cu]
                        C[ru, i] = Zc[cu, i]
                    for v in range(2):
                        xv = idx[v]
                        C[ru, d + 2 * v] = T[cu, 2 * xv]
                        C[ru, d + 2 * v + 1] = T[cu, 2 * xv + 1]
                C[d + 2 * u, d + 2 * u + 1] ^= dinv[xu]
                C[d + 2 * u + 1, d + 2 * u] ^= dinv[xu]
            if gf_rank(C, exp, log) < 4:
                obs[x, x2 // 62] |= np.int64(1) << (x2 % 62)
                obs[x2, x // 62] |= np.int64(1) << (x % 62)
    return obs


@njit(cache=True)
def _popcount(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@njit(cache=True)
def _filter_cands(h1, h2, ess, obs, ncand, vmax, out):
    """Candidates a ranked by violated necessary conditions for zero loss: v = 16 x (essential lines killed)
    + (obstructing pairs killed). A test over GF(2^12) reports a false obstruction about once in a thousand
    pairs, so candidates with up to `vmax` violations are kept (the exact count decides afterwards).
    Writes (a, v, killed) triples to out[0], out[1], out[2] and returns their number."""
    n = len(h1)
    nw = obs.shape[1]
    kmask = np.zeros(nw, np.int64)
    klist = np.zeros(n, np.int64)
    cap = out.shape[1]
    cnt = 0
    for a in range(1, ncand):
        k = 0
        v = 0
        for q in range(nw):
            kmask[q] = 0
        for x in range(n):
            if _parity(a & h1[x]) | _parity(a & h2[x]):
                if ess[x]:
                    v += 16
                    if v > vmax:
                        break
                klist[k] = x
                k += 1
                kmask[x // 62] |= np.int64(1) << (x % 62)
        if v > vmax:
            continue
        twice = 0
        for t in range(k):
            x = klist[t]
            for q in range(nw):
                twice += _popcount(obs[x, q] & kmask[q])
            if v + twice // 2 > vmax:
                break
        v += twice // 2
        if v <= vmax and cnt < cap:
            out[0, cnt] = a
            out[1, cnt] = v
            out[2, cnt] = k
            cnt += 1
    return cnt


@njit(cache=True)
def _open_hist(blk, pt_blk, pt_ymask, pt_outs, pt_phase, ncoord, Rm, ecount, nb):
    """gain[lam] * 256 for every direction lam (bit j = coefficient on relation j of block blk; opening lam
    keeps the relations sum a_j r_j with even lam.a): the number of point pairs (p in blk, p' elsewhere) whose
    fan-out spaces gain a common dimension. A pair gains iff lam lies in L = C_p(sigma_p'), the image of the
    partner's fan-outs under the constraint map of p; the indicator of L is the average of the characters of
    its annihilator, so all 2^e directions are counted with one transform."""
    npnt = len(pt_blk)
    e = ecount[blk]
    size = 1 << e
    hist = np.zeros(size, np.int64)
    cons = _point_cons(pt_blk, pt_ymask, pt_outs, ncoord, Rm, ecount)
    sig = np.zeros(64, np.int64)
    img = np.zeros(64, np.int64)
    ann = np.zeros(64, np.int64)
    full = (np.int64(1) << e) - 1
    for p in range(npnt):
        if pt_blk[p] != blk:
            continue
        for q in range(npnt):
            if pt_blk[q] == blk:
                continue
            if (pt_phase[p] & pt_phase[q]) == 0:
                continue
            sup = pt_outs[p] & pt_outs[q]
            if sup == 0:
                continue
            d = _nullspace(cons[q], ecount[pt_blk[q]], sup, nb, sig)
            if d == 0:
                continue
            for t in range(d):
                v = np.int64(0)
                for j in range(e):
                    if _parity(cons[p, j] & sig[t]):
                        v |= np.int64(1) << j
                img[t] = v
            da = _nullspace(img, d, full, e, ann)
            if da > 8:
                continue
            w = np.int64(256) >> da
            for m in range(1 << da):
                x = np.int64(0)
                for t in range(da):
                    if (m >> t) & 1:
                        x ^= ann[t]
                hist[x] += w
    h = 1
    while h < size:
        for s in range(0, size, 2 * h):
            for q in range(s, s + h):
                x = hist[q]
                y = hist[q + h]
                hist[q] = x + y
                hist[q + h] = x - y
        h *= 2
    return hist


@njit(cache=True)
def _block_caps(n, li, lj, ls, pt_blk, pt_ymask, ncoord, blk, Qm, nq):
    """cap2[a] = 2 x (number of block-blk features carried by the lines that survive candidate a).
    A line dies when one of its features in the block has odd a.h (h = class of the feature in the basis Q).
    Counted for all a at once with one Walsh-Hadamard transform."""
    size = 1 << nq
    hist = np.zeros(size, np.int64)
    for ell in range(n):
        i = li[ell]
        j = lj[ell]
        bi = pt_blk[i] == blk
        bj = pt_blk[j] == blk
        if not (bi or bj):
            continue
        ha = np.int64(0)
        hb = np.int64(0)
        for side in range(2):
            pt = i if side == 0 else j
            if pt_blk[pt] != blk:
                continue
            hm = np.int64(0)
            for t in range(nq):
                cv = np.int64(0)
                for ci in range(ncoord[blk]):
                    if (pt_ymask[pt] >> ci) & 1:
                        cv ^= Qm[t, ci]
                if _parity(cv & ls[ell]):
                    hm |= np.int64(1) << t
            if side == 0:
                ha = hm
            else:
                hb = hm
        if bi and bj:                      # both features in the block: survives iff a.ha = a.hb = 0, weight 2
            hist[0] += 1
            hist[ha] += 1
            hist[hb] += 1
            hist[ha ^ hb] += 1
        elif bi:
            hist[0] += 1
            hist[ha] += 1
        else:
            hist[0] += 1
            hist[hb] += 1
    h = 1
    while h < size:
        for s in range(0, size, 2 * h):
            for q in range(s, s + h):
                x = hist[q]
                y = hist[q + h]
                hist[q] = x + y
                hist[q + h] = x - y
        h *= 2
    return hist


@njit(cache=True)
def _scan_counts(h1, h2, ncand):
    """Number of lines killed by every candidate a (a line dies when a feature of it has odd a.h)."""
    out = np.zeros(ncand, np.int64)
    n = len(h1)
    for a in range(1, ncand):
        c = 0
        for x in range(n):
            if _parity(a & h1[x]) | _parity(a & h2[x]):
                c += 1
        out[a] = c
    return out


@njit(cache=True)
def _scan_exact(T, Zr, Zc, dinv, h1, h2, cands, rho, stop_at, exp, log):
    """Exact pairing number after each candidate (in order; stops early once `stop_at` is reached)."""
    n = len(h1)
    d = Zr.shape[0]
    out = np.full(len(cands), -1, np.int64)
    idx = np.zeros(n, np.int64)
    for q in range(len(cands)):
        a = cands[q]
        k = 0
        for x in range(n):
            if _parity(a & h1[x]) | _parity(a & h2[x]):
                idx[k] = x
                k += 1
        size = d + 2 * k
        C = np.zeros((size, size), np.uint16)
        for u in range(k):
            xu = idx[u]
            for su in range(2):
                cu = 2 * xu + su
                ru = d + 2 * u + su
                for i in range(d):
                    C[i, ru] = Zr[i, cu]
                    C[ru, i] = Zc[cu, i]
                for v in range(k):
                    xv = idx[v]
                    C[ru, d + 2 * v] = T[cu, 2 * xv]
                    C[ru, d + 2 * v + 1] = T[cu, 2 * xv + 1]
            C[d + 2 * u, d + 2 * u + 1] ^= dinv[xu]
            C[d + 2 * u + 1, d + 2 * u] ^= dinv[xu]
        rk = gf_rank(C, exp, log)
        out[q] = (rho + rk - 2 * k) // 2
        if out[q] >= stop_at:
            break
    return out


# ------------------------------------------------------------------------------ don't-care addresses
def relax_last_row(an: core.Analysis) -> int:
    """Use the don't-care addresses of the last, partial address block (ours; the public notes only say that
    padding rows are don't-cares).

    The table ends inside the last block of 2**lb addresses (4 valid addresses of 16 for 1,220 rows). The
    values a circuit produces on the other addresses of that block are free, and they enter exactly one row
    of the mixed tensor, the row of the top high monomial, in the low monomials that contain an address bit
    beyond the valid ones. chem_loaders_core sets those coefficients to zero. Here they are chosen, block by
    block, so that this row becomes a combination of the other rows inside every block that has a free
    coordinate: the row space of the block loses a dimension, i.e. the block gains one relation to choose
    from (which matters most where a block has no spare relation left).
    Only existing cells are written (no output gains a low coordinate). Returns the number of blocks whose
    row space shrank. Modifies an.X0 in place; emission needs no change (it checks valid addresses only)."""
    A, L, H, nr = an.A, an.L, an.H, an.nr
    nv = A - L * (H - 1)
    if nv <= 0 or nv >= L or nv & (nv - 1) or an.X0.shape[1] == 0:
        return 0
    mono = []
    for bs in an.basis:
        anf = bs["anf"]
        if anf & (anf - 1):
            return 0                        # a low coordinate that is not a single monomial: not handled
        mono.append(anf.bit_length() - 1)
    free = sum(1 << t for t, m in enumerate(mono) if m >= nv)
    X = an.X0
    nb = X.shape[1]
    last = nr - 1
    assert not (int(np.bitwise_or.reduce(X[last])) & free), "free coefficients of the last row are not zero"
    bd = core.BlockDecomp(an, 0, 1)
    bd.prepare()
    changed = 0
    for bi in bd.binfo:
        ts, cell = bi["ts"], bi["cell"]
        fixed_cells = sum(1 << idx for (t, o), idx in cell.items() if not free >> t & 1)
        if fixed_cells == (1 << bi["amb"]) - 1:
            continue                        # no free cell in this block
        rows = bi["rows"]
        sp = Span()
        for i, v in enumerate(rows[:last]):
            r, c = sp.reduce(v & fixed_cells, 1 << i)
            if r:
                sp.piv[r.bit_length() - 1] = (r, c)
        r, combo = sp.reduce(rows[last] & fixed_cells)
        if r:
            continue                        # the fixed part is not reachable: keep the row as it is
        new = 0
        for i in bits(combo):
            new ^= rows[i]
        assert (new & fixed_cells) == (rows[last] & fixed_cells)
        full = Span()
        before = sum(1 for v in rows if full.add(v))
        if new == rows[last]:
            continue
        for (t, o), idx in cell.items():
            if free >> t & 1:
                if new >> idx & 1:
                    X[last, o] |= np.uint32(1 << t)
                else:
                    X[last, o] &= np.uint32(~(1 << t) & 0xFFFFFFFF)
        full = Span()
        after = sum(1 for i, v in enumerate(rows) if full.add(new if i == last else v))
        changed += before - after
    return changed


# ---------------------------------------------------------------------------------------------- engine
class Engine:
    """Relation spans of the blocks of one circuit, their exact price, the swap search and the final plan."""

    def __init__(self, an: core.Analysis, seed: int = 0):
        self.an = an
        bd = core.BlockDecomp(an, seed, 1)
        bd.prepare()
        self.bd, self.binfo, self.nb = bd, bd.binfo, bd.nb
        assert self.nb <= 62
        self.rng = random.Random(1000003 * seed + 7919 * an.k + 4242)
        nblk = len(self.binfo)
        self.nblk = nblk
        self.ncoord = np.array([len(bi["ts"]) for bi in self.binfo], np.int64)
        assert self.ncoord.max() <= CMAX
        self.ncell = sum(bi["amb"] for bi in self.binfo)
        self.cellid = np.full((nblk, CMAX, 64), self.ncell, np.int64)   # index ncell = "no cell"
        self.loc = []
        for b, bi in enumerate(self.binfo):
            ci_of = {t: ci for ci, t in enumerate(bi["ts"])}
            loc = [None] * bi["amb"]
            for (t, o), idx in bi["cell"].items():
                self.cellid[b, ci_of[t], o] = bi["offset"] + idx
                loc[idx] = (ci_of[t], o)
            self.loc.append(loc)
        pts = []
        for b, bi in enumerate(self.binfo):
            ts = bi["ts"]
            ph = sum(1 << p for p in bi["phases"])
            for ymask in range(1, 1 << len(ts)):
                yt = [ts[i] for i in bits(ymask)]
                outs = sum(1 << o for o in range(self.nb) if all(bd.supp[o] >> t & 1 for t in yt))
                if outs:
                    pts.append((b, ymask, outs, ph, sum(1 << t for t in yt), yt))
        self.points = pts
        self.pt_blk = np.array([p[0] for p in pts], np.int64)
        self.pt_ymask = np.array([p[1] for p in pts], np.int64)
        self.pt_outs = np.array([p[2] for p in pts], np.int64)
        self.pt_phase = np.array([p[3] for p in pts], np.int64)
        self.R = [[] for _ in range(nblk)]
        self.Rm = np.zeros((nblk, EMAX, CMAX), np.int64)
        self.ecount = np.zeros(nblk, np.int64)
        self.pivmask = [None] * nblk
        self.dirty = [True] * nblk
        self.li = np.zeros(MAXLINES, np.int64)
        self.lj = np.zeros(MAXLINES, np.int64)
        self.ls = np.zeros(MAXLINES, np.int64)
        self.deg = np.zeros(len(pts), np.int64)
        xr = np.random.default_rng(977 * seed + an.k + 31)
        self.lx = xr.integers(1, 1 << GF_BITS, MAXLINES).astype(np.uint16)
        self.evals = 0
        self.scans = 0
        self.trials = 0
        self.swaps2 = 0
        self.filtered = 0
        self.mix = 0.25
        self.vmax = 4
        self.pair_filter = False
        self.exacts = 0
        self.laterals = 0
        self.last_free = 0
        self.mismatch = 0
        self.log = []

    # ---- spans
    def max_e(self, b):
        return min(EMAX, len(self.binfo[b]["rel"]))

    def random_relation(self, b, sp: Span):
        rel = self.binfo[b]["rel"]
        rng = self.rng
        for _ in range(200):
            v = 0
            for r in rel:
                if rng.random() < 0.5:
                    v ^= r
            if v and sp.reduce(v)[0]:
                return v
        return None

    def init_spans(self, evec):
        for b in range(self.nblk):
            sp = Span()
            R = []
            while len(R) < min(evec[b], self.max_e(b)):
                v = self.random_relation(b, sp)
                if v is None:
                    break
                sp.add(v)
                R.append(v)
            self.R[b] = R
            self.dirty[b] = True

    def set_block(self, b, R):
        self.R[b] = R
        self.dirty[b] = True

    def _sync(self):
        for b in range(self.nblk):
            if not self.dirty[b]:
                continue
            bi = self.binfo[b]
            R = self.R[b]
            self.Rm[b] = 0
            sp = Span()
            loc = self.loc[b]
            for j, v in enumerate(R):
                ok = sp.add(v)
                assert ok, "dependent relations"
                for idx in bits(v):
                    ci, o = loc[idx]
                    self.Rm[b, j, ci] |= 1 << o
            self.ecount[b] = len(R)
            m = np.ones(bi["amb"], np.bool_)
            for p in sp.piv:
                m[p] = False
            self.pivmask[b] = m
            self.dirty[b] = False
        free = np.concatenate(self.pivmask + [np.zeros(1, np.bool_)])
        idx = np.cumsum(free) - 1
        idx[~free] = -1
        self.freeidx = idx.astype(np.int64)
        self.r = int(free.sum())

    @property
    def codim(self):
        return [len(R) for R in self.R]

    # ---- exact price
    def price(self):
        """(products, lines): products = r - nu with nu = rank(Y)/2 the exact maximum number of pairs
        (a lower bound with probability about 1/4096 of being one short)."""
        self._sync()
        n = _lines(self.pt_blk, self.pt_ymask, self.pt_outs, self.pt_phase, self.ncoord, self.Rm, self.ecount,
                   self.nb, self.li, self.lj, self.ls, self.deg)
        Y = _build_y(n, self.li, self.lj, self.ls, self.lx, self.pt_blk, self.pt_ymask, self.ncoord, self.cellid,
                     self.freeidx, self.nb, self.r, self.r)
        rank = gf_rank(Y, GF_EXP, GF_LOG)
        self.evals += 1
        self.nlines = n
        return self.r - rank // 2, n

    # ---- moves
    def swap(self, b):
        """A random neighbour of R_b in the Grassmann graph: drop one direction, add a new relation."""
        R = self.R[b]
        e = len(R)
        rng = self.rng
        j = rng.randrange(e)
        new = [R[i] ^ (R[j] if rng.random() < 0.5 else 0) for i in range(e) if i != j]
        sp = Span()
        for v in new:
            sp.add(v)
        v = self.random_relation(b, sp)
        if v is None:
            return None
        new.append(v)
        return new

    def grow(self, b):
        R = self.R[b]
        if len(R) >= self.max_e(b):
            return None
        sp = Span()
        for v in R:
            sp.add(v)
        v = self.random_relation(b, sp)
        if v is None:
            return None
        return R + [v]

    def climb(self, iters, cur=None, stop_at=None):
        """Swap search at fixed codimensions; accepts moves that do not raise the price."""
        if cur is None:
            cur, _ = self.price()
        rng = self.rng
        for it in range(iters):
            if stop_at is not None and cur <= stop_at:
                break
            b = rng.randrange(self.nblk)
            if not self.R[b]:
                continue
            new = self.swap(b)
            if new is None:
                continue
            old = self.R[b]
            self.set_block(b, new)
            p, _ = self.price()
            if p <= cur:
                cur = p
            else:
                self.set_block(b, old)
        return cur

    def snapshot(self):
        return [list(R) for R in self.R]

    def restore(self, snap):
        for b, R in enumerate(snap):
            self.set_block(b, list(R))

    # ---- exact scans of all single-relation additions
    def _complement(self, b, R, cap=13):
        """Vectors spanning the relations of block b modulo span(R); a random `cap`-dimensional part if more."""
        sp = Span()
        for v in R:
            sp.add(v)
        out = []
        for v in self.binfo[b]["rel"]:
            if sp.reduce(v)[0]:
                sp.add(v)
                out.append(v)
        if len(out) > cap:
            rng = self.rng
            mixed, sp2 = [], Span()
            while len(mixed) < cap:
                v = 0
                for q in out:
                    if rng.random() < 0.5:
                        v ^= q
                if v and sp2.add(v):
                    mixed.append(v)
            out = mixed
        return out

    def prepare(self):
        """Lines, Lovasz matrix and its echelon transform for the current spans. Returns nu."""
        self._sync()
        n = _lines(self.pt_blk, self.pt_ymask, self.pt_outs, self.pt_phase, self.ncoord, self.Rm, self.ecount,
                   self.nb, self.li, self.lj, self.ls, self.deg)
        Y = _build_y(n, self.li, self.lj, self.ls, self.lx, self.pt_blk, self.pt_ymask, self.ncoord, self.cellid,
                     self.freeidx, self.nb, self.r, self.r)
        rho, ET, Rnp, pivrow, npidx = _transform(Y, GF_EXP, GF_LOG)
        self.tr = (n, rho, ET, Rnp, pivrow, npidx)
        self.evals += 1
        self.scans += 1
        self.nlines = n
        return rho // 2

    def scan(self, b, Q, need=None, ntop=64, cap=256, cands=None):
        """Exact pairing numbers after adding the relation sum_{t in a} Q[t] to block b (call prepare() first).
        cands given: exactly those candidates. Otherwise, need = nu of the current spans: zero-loss mode,
        every candidate that keeps all essential lines of the block is evaluated (at most `cap`, fewest killed
        lines first; exhaustive when fewer exist); else the `ntop` candidates that kill the fewest lines.
        Returns [(nu', a, killed)] in evaluation order (stops at the first candidate reaching `need`)."""
        n, rho, ET, Rnp, pivrow, npidx = self.tr
        tb = self.pt_blk
        sel = np.flatnonzero((tb[self.li[:n]] == b) | (tb[self.lj[:n]] == b)).astype(np.int64)
        nq = len(Q)
        Qm = self._qm(b, Q)
        T, Zr, Zc, dinv, h1, h2 = _scan_block(rho, ET, Rnp, pivrow, npidx, sel, self.li, self.lj, self.ls, self.lx,
                                              self.pt_blk, self.pt_ymask, self.ncoord, self.cellid, self.freeidx,
                                              self.nb, b, Qm, nq, GF_EXP, GF_LOG)
        zero_loss = need is not None and need >= rho // 2
        if cands is not None:
            cands = [int(z) for z in cands]
            killed = None
        elif zero_loss:
            killed = _scan_counts(h1, h2, 1 << nq)
            ess = _essential(T, Zr, Zc, dinv)
            if self.pair_filter:
                obs = _pair_obstructions(T, Zr, Zc, dinv, ess, GF_EXP, GF_LOG)
            else:
                obs = np.zeros((len(dinv), (len(dinv) + 61) // 62 + 1), np.int64)
            buf = np.zeros((3, 8192), np.int64)
            cnt = _filter_cands(h1, h2, ess, obs, 1 << nq, self.vmax, buf)
            self.last_free = cnt
            self.filtered += cnt
            order = sorted(range(cnt), key=lambda q: (int(buf[1, q]), int(buf[2, q])))
            cands = [int(buf[0, q]) for q in order[:cap]]
        else:
            killed = _scan_counts(h1, h2, 1 << nq)
            order = np.argsort(killed[1:], kind="stable") + 1
            cands = [int(z) for z in order[:ntop]]
        if not cands:
            return []
        stop_at = need if need is not None else 10 ** 9
        nus = _scan_exact(T, Zr, Zc, dinv, h1, h2, np.array(cands, np.int64), rho, stop_at, GF_EXP, GF_LOG)
        self.exacts += int((nus >= 0).sum())
        return [(int(nu), z, int(killed[z]) if killed is not None else 0) for nu, z in zip(nus, cands) if nu >= 0]

    def _vec(self, Q, a):
        v = 0
        for t in bits(a):
            v ^= Q[t]
        return v

    def grow_scan(self, ntop=96, pick=0):
        """The single relation (over all blocks) whose addition keeps the most pairs, by exact count.
        Applies it and returns the new price, or None when no block can grow. pick > 0 takes a random one
        of the best `pick + 1` (used on retries)."""
        nu0 = self.prepare()
        need = min(nu0, (self.r - 1) // 2)
        found = []
        for b in range(self.nblk):
            if len(self.R[b]) >= self.max_e(b):
                continue
            Q = self._complement(b, self.R[b])
            if not Q:
                continue
            res = self.scan(b, Q, need=need, ntop=ntop)
            if not res and need >= nu0:
                res = self.scan(b, Q, need=None, ntop=ntop)
            for nu, a, kl in res:
                found.append((-nu, kl, self.rng.random(), b, a, Q))
        if not found:
            return None
        found.sort(key=lambda z: z[:3])
        _, _, _, b, a, Q = found[self.rng.randrange(min(len(found), pick + 1)) if pick else 0]
        self.set_block(b, self.R[b] + [self._vec(Q, a)])
        p, _ = self.price()
        return p

    def _qm(self, b, Q):
        Qm = np.zeros((max(1, len(Q)), CMAX), np.int64)
        loc = self.loc[b]
        for t, v in enumerate(Q):
            for idx in bits(v):
                ci, o = loc[idx]
                Qm[t, ci] |= 1 << o
        return Qm

    def caps(self, b, Q):
        """For every relation sum_{t in a} Q[t] that could be added to block b: the number of block-b features
        left on the candidate lines (an upper bound on what a pairing can use there). Needs fresh lines."""
        return _block_caps(self.nlines, self.li, self.lj, self.ls, self.pt_blk, self.pt_ymask, self.ncoord, b,
                           self._qm(b, Q), len(Q)) // 2

    def fresh_lines(self):
        self._sync()
        self.nlines = _lines(self.pt_blk, self.pt_ymask, self.pt_outs, self.pt_phase, self.ncoord, self.Rm,
                             self.ecount, self.nb, self.li, self.lj, self.ls, self.deg)
        return self.nlines

    def block_dim(self, b):
        return self.binfo[b]["amb"] - len(self.R[b])

    def trial(self, b, cur, ntop=16, lateral=0.0):
        """Open one random direction of R_b, then look at every relation of every block that could close the
        spans again. If the opening raised the pairing number by exactly what is needed, ALL closings that
        keep it are found: candidates that kill an essential line or an obstructing pair of lines are
        discarded (exact necessary tests), the rest are priced exactly from one elimination. Otherwise the
        `ntop` candidates of each block that keep the most features are priced. Never raises the price; with
        probability `lateral` an equal-price closing is taken when no better one exists. Returns the new price."""
        R = self.R[b]
        e = len(R)
        if e == 0:
            return cur
        rng = self.rng
        j = rng.randrange(e)
        H = [R[i] ^ (R[j] if rng.random() < 0.5 else 0) for i in range(e) if i != j]
        self._sync()
        r0 = self.r
        nu0 = r0 - cur
        need = min(nu0 + 1, r0 // 2)
        slack = r0 - 2 * need              # deficiency allowed after closing
        snap = self.snapshot()
        self.set_block(b, H)
        nu1 = self.prepare()
        self.trials += 1
        go_lateral = bool(lateral) and rng.random() < lateral
        if nu1 < need and not go_lateral:
            self.restore(snap)
            return cur
        order = list(range(self.nblk))
        rng.shuffle(order)
        best = None
        side = []
        top = min(nu1, r0 // 2)
        for b2 in order:
            if len(self.R[b2]) >= self.max_e(b2):
                continue
            if b2 == b:
                Q = [R[j]] + self._complement(b2, H + [R[j]], cap=12)
            else:
                Q = self._complement(b2, self.R[b2], cap=12)
            if not Q:
                continue
            res = []
            if nu1 == need:
                res = self.scan(b2, Q, need=nu1, cap=ntop)
                if b2 == b:
                    res = [x for x in res if x[1] != 1]
            elif nu1 > need:
                cap = self.caps(b2, Q)
                cap[0] = -1
                if b2 == b:
                    cap[1] = -1
                thresh = self.block_dim(b2) - 1 - slack
                good = [int(z) for z in np.argsort(-cap, kind="stable")[:ntop] if cap[z] >= max(thresh, 0)]
                if good:
                    res = self.scan(b2, Q, need=top, cands=good)
            for nu, z, _ in res:
                if nu >= need and (best is None or nu > best[0]):
                    best = (nu, b2, self._vec(Q, z))
            if best is not None and best[0] >= top:
                break
            if best is None and go_lateral and len(side) < 4:
                cap = self.caps(b2, Q)
                cap[0] = -1
                if b2 == b:
                    cap[1] = -1
                good = [int(z) for z in np.argsort(-cap, kind="stable")[:6] if cap[z] >= 0]
                for nu, z, _ in (self.scan(b2, Q, need=None, cands=good) if good else []):
                    if nu == nu0:
                        side.append((b2, self._vec(Q, z)))
        if best is None and go_lateral and side:
            b2, v = side[rng.randrange(len(side))]
            self.set_block(b2, self.R[b2] + [v])
            p, _ = self.price()
            if p == cur:
                self.laterals += 1
                return p
            self.restore(snap)
            return cur
        if best is None:
            self.restore(snap)
            return cur
        _, b2, v = best
        self.set_block(b2, self.R[b2] + [v])
        p, _ = self.price()
        if p > cur:                        # the two exact counts disagree (field too small): keep the old spans
            self.mismatch += 1
            self.restore(snap)
            return cur
        return p

    def grow_zero_loss(self):
        """Exhaustive: a relation of any block whose addition keeps EVERY pair (essential-line and
        obstructing-pair filters, then exact counts). Applies it and returns the new price, else None."""
        nu0 = self.prepare()
        r0 = self.r
        order = list(range(self.nblk))
        self.rng.shuffle(order)
        for b in order:
            if len(self.R[b]) >= self.max_e(b):
                continue
            Q = self._complement(b, self.R[b], cap=12)
            if not Q:
                continue
            for nu, z, _ in self.scan(b, Q, need=nu0, cap=48):
                if nu >= nu0:
                    old = self.R[b]
                    self.set_block(b, old + [self._vec(Q, z)])
                    p, _ = self.price()
                    if p <= r0 - 1 - nu0:
                        return p
                    self.mismatch += 1
                    self.set_block(b, old)
        return None

    def lateral_step(self, cur, ntop=10):
        """One random step on the plateau of equal price: open a direction of a random block, close with a
        different relation (of any block) that gives the same pairing number. Returns True if the spans moved."""
        rng = self.rng
        b = rng.randrange(self.nblk)
        R = self.R[b]
        e = len(R)
        if e == 0:
            return False
        j = rng.randrange(e)
        H = [R[i] ^ (R[j] if rng.random() < 0.5 else 0) for i in range(e) if i != j]
        self._sync()
        r0 = self.r
        nu0 = r0 - cur
        snap = self.snapshot()
        self.set_block(b, H)
        nu1 = self.prepare()
        self.trials += 1
        order = list(range(self.nblk))
        rng.shuffle(order)
        for b2 in order:
            if len(self.R[b2]) >= self.max_e(b2):
                continue
            if b2 == b:
                Q = [R[j]] + self._complement(b2, H + [R[j]], cap=12)
            else:
                Q = self._complement(b2, self.R[b2], cap=12)
            if not Q:
                continue
            if nu1 == nu0:
                res = self.scan(b2, Q, need=nu0, cap=24)
                if b2 == b:
                    res = [x for x in res if x[1] != 1]
            else:
                cap = self.caps(b2, Q)
                cap[0] = -1
                if b2 == b:
                    cap[1] = -1
                good = [int(z) for z in np.argsort(-cap, kind="stable")[:ntop] if cap[z] >= 0]
                res = self.scan(b2, Q, need=None, cands=good) if good else []
            ok = [z for nu, z, _ in res if nu >= nu0]
            if ok:
                z = ok[rng.randrange(len(ok))]
                self.set_block(b2, self.R[b2] + [self._vec(Q, z)])
                p, _ = self.price()
                if p <= cur:
                    self.laterals += 1
                    return True
                self.mismatch += 1
                break
        self.restore(snap)
        return False

    def total_dim(self):
        """Sum over point pairs of the dimension of their common fan-outs (needs fresh lines)."""
        return int(self.deg.sum()) // 2

    def swap2(self, cur, tot=None, nopen=3, nexact=6):
        """Close-then-open exchange. A relation that keeps the most features is added to a random block; then
        EVERY direction of every block's relation space is ranked by the number of point pairs that gain a
        common fan-out when it is opened (one transform per block), and the best few are priced exactly.
        Accepted if the price falls, or stays while the total common-fan-out dimension grows.
        Returns (price, total dimension)."""
        rng = self.rng
        if tot is None:
            self.fresh_lines()
            tot = self.total_dim()
        snap = self.snapshot()
        b = rng.randrange(self.nblk)
        if len(self.R[b]) >= self.max_e(b):
            return cur, tot
        self.fresh_lines()
        Q = self._complement(b, self.R[b], cap=12)
        if not Q:
            return cur, tot
        cap = self.caps(b, Q)
        cap[0] = -1
        top = [int(z) for z in np.argsort(-cap, kind="stable")[:3] if cap[z] >= 0]
        if not top:
            return cur, tot
        a = top[rng.randrange(len(top))]
        added = self._vec(Q, a)
        self.set_block(b, self.R[b] + [added])
        self._sync()
        self.swaps2 += 1
        cands = []
        for b2 in range(self.nblk):
            e2 = len(self.R[b2])
            if e2 < 2 or e2 > 22:
                continue
            hist = _open_hist(b2, self.pt_blk, self.pt_ymask, self.pt_outs, self.pt_phase, self.ncoord,
                              self.Rm, self.ecount, self.nb)
            hist[0] = -1
            if b2 == b:
                hist[1 << (e2 - 1)] = -1        # that direction would only undo the relation just added
            for lam in np.argsort(-hist, kind="stable")[:nopen]:
                cands.append((int(hist[lam]), rng.random(), b2, int(lam)))
        cands.sort(reverse=True)
        best = None
        for gain, _, b2, lam in cands[:nexact]:
            R2 = self.R[b2]
            j0 = (lam & -lam).bit_length() - 1
            opened = [R2[i] ^ (R2[j0] if (lam >> i) & 1 else 0) for i in range(len(R2)) if i != j0]
            self.set_block(b2, opened)
            p, _ = self.price()
            self.exacts += 1
            t2 = self.total_dim()
            if best is None or (p, -t2) < (best[0], -best[1]):
                best = (p, t2, b2, opened)
            self.set_block(b2, R2)
            if p < cur:
                break
        if best is not None and (best[0] < cur or (best[0] == cur and best[1] > tot)):
            self.set_block(best[2], best[3])
            self._sync()
            if best[0] == cur:
                self.laterals += 1
            return best[0], best[1]
        self.restore(snap)
        return cur, tot

    def level_up(self, cur, quick=60, trials=150, pick=0, lateral=0.15, walk=60, verbose=False):
        """One more relation with a perfect pairing. Returns (price, perfect?). Three routes, cheapest first:
        1. a few random relations and random swaps;
        2. the best relation by exact scan, then open/close trials until the pairing is perfect again;
        3. (from a state with one unpaired feature, where the new relation must lose nothing) walk the
           plateau of equal price and test every state exhaustively for such a relation."""
        rng = self.rng
        base = self.snapshot()
        self._sync()
        odd = self.r % 2 == 1 and cur == (self.r + 1) // 2
        cands = []
        for b in range(self.nblk):
            for _ in range(2):
                new = self.grow(b)
                if new is None:
                    continue
                old = self.R[b]
                self.set_block(b, new)
                p, _ = self.price()
                cands.append((p, rng.random(), b, new))
                self.set_block(b, old)
        if not cands:
            return cur, False
        cands.sort()
        p, _, b, new = cands[0]
        self.set_block(b, new)
        self._sync()
        goal = (self.r + 1) // 2
        if p > goal and quick and p <= goal + 5:
            p = self.climb(quick, p, stop_at=goal)
        if p == goal:
            return p, True
        self.restore(base)
        p = self.grow_scan(ntop=96, pick=pick)
        if p is None:
            return cur, False
        if verbose:
            print(f"    grow scan: price {p} goal {goal} codim {self.codim} lines {self.nlines}", flush=True)
        t = 0
        last = 0
        tot = None
        while p > goal and t - last < trials:
            if rng.random() < self.mix:
                q, tot = self.swap2(p, tot)
                kind = "swap2"
            else:
                b = rng.randrange(self.nblk)
                q = self.trial(b, p, lateral=lateral)
                tot = None
                kind = "trial"
            if q < p:
                last = t
                if verbose:
                    print(f"    {kind} {t}: price {p} -> {q} codim {self.codim} lines {self.nlines} "
                          f"exact {self.exacts} laterals {self.laterals}", flush=True)
            p = q
            t += 1
        if p == goal or not odd or not walk:
            return p, p == goal
        lossy = (p, self.snapshot())
        self.restore(base)
        for step in range(walk):
            q = self.grow_zero_loss()
            if q is not None:
                if verbose:
                    print(f"    zero-loss relation after {step} plateau steps: price {q} codim {self.codim} "
                          f"lines {self.nlines} exact {self.exacts} filtered {self.filtered}", flush=True)
                return q, q == goal
            for _ in range(8):
                if self.lateral_step(cur):
                    break
        if verbose:
            print(f"    no zero-loss relation on {walk} plateau states (laterals {self.laterals})", flush=True)
        self.restore(lossy[1])
        return lossy[0], False

    # ---- final plan
    def plan(self, tries=4):
        """Products [y1, A1, y2, A2, fan-out, phase] for the current spans (exact GF(2) check included)."""
        an, binfo = self.an, self.binfo
        self._sync()
        r = self.r
        last = None
        for attempt in range(tries):
            xr = np.random.default_rng(1009 * attempt + 17 * an.k + 5)
            n = _lines(self.pt_blk, self.pt_ymask, self.pt_outs, self.pt_phase, self.ncoord, self.Rm, self.ecount,
                       self.nb, self.li, self.lj, self.ls, self.deg)
            lx = xr.integers(1, 1 << GF_BITS, MAXLINES).astype(np.uint16)
            Y = _build_y(n, self.li, self.lj, self.ls, lx, self.pt_blk, self.pt_ymask, self.ncoord, self.cellid,
                         self.freeidx, self.nb, r, r)
            nu = gf_rank(Y.copy(), GF_EXP, GF_LOG) // 2
            d = r - 2 * nu
            ntot = n
            if d:
                ntot = _singles(self.pt_blk, self.pt_ymask, self.pt_outs, self.ncoord, self.Rm, self.ecount,
                                self.nb, d, n, self.li, self.lj, self.ls)
            Y = _build_y(ntot, self.li, self.lj, self.ls, lx, self.pt_blk, self.pt_ymask, self.ncoord,
                         self.cellid, self.freeidx, self.nb, r, r + d)
            ok, N = gf_inverse(Y, GF_EXP, GF_LOG)
            if not ok:
                last = f"attempt {attempt}: Y singular (nu {nu}, deficiency {d})"
                continue
            order = xr.permutation(ntot).astype(np.int64)
            keep = _extract(N, ntot, self.li, self.lj, self.ls, lx, order, self.pt_blk, self.pt_ymask,
                            self.ncoord, self.cellid, self.freeidx, self.nb, r, GF_EXP, GF_LOG)
            kept = np.flatnonzero(keep)
            pairs = [(int(self.li[q]), int(self.lj[q]), int(self.ls[q])) for q in kept if self.lj[q] >= 0]
            singles = [(int(self.li[q]), int(self.ls[q])) for q in kept if self.lj[q] < 0]
            if 2 * len(pairs) + len(singles) != r:
                last = f"attempt {attempt}: kept {len(pairs)} pairs + {len(singles)} singles for rank {r}"
                continue
            prods = self._products(pairs, singles)
            if prods is None:
                last = f"attempt {attempt}: extracted features are dependent over GF(2)"
                continue
            self.log.append(f"plan: codim {self.codim} cells {self.ncell} rank {r} lines {n} pairs {len(pairs)} "
                            f"singles {len(singles)} -> {len(prods)} products")
            return prods
        raise RuntimeError(f"no plan extracted: {last}")

    def _flat(self, i, s):
        b, _, _, _, _, yt = self.points[i]
        bi = self.binfo[b]
        cell = bi["cell"]
        v = 0
        for t in yt:
            for o in bits(s):
                v |= 1 << cell[(t, o)]
        return v << bi["offset"]

    def _products(self, pairs, singles):
        binfo, points = self.binfo, self.points
        feats = []
        for n, (i, j, s) in enumerate(pairs):
            feats.append((self._flat(i, s), ("p", n, 0)))
            feats.append((self._flat(j, s), ("p", n, 1)))
        for n, (i, s) in enumerate(singles):
            feats.append((self._flat(i, s), ("s", n)))
        fs = Span()
        for idx, (v, _) in enumerate(feats):
            if not fs.add(v, 1 << idx):
                return None
        coef = [0] * len(feats)
        for bi in binfo:
            for M, v in enumerate(bi["rows"]):
                rr, c = fs.reduce(v << bi["offset"])
                if rr:
                    return None
                for idx in bits(c):
                    coef[idx] |= 1 << M
        A = {key: coef[idx] for idx, (_, key) in enumerate(feats)}

        def phase_of(*pts):
            common = None
            for p in pts:
                ph = set(binfo[points[p][0]]["phases"])
                common = ph if common is None else common & ph
            return min(common)

        prods = []
        for n, (i, j, s) in enumerate(pairs):
            a1, a2 = A[("p", n, 0)], A[("p", n, 1)]
            if a1 and a2:
                prods.append([points[i][4], a1, points[j][4], a2, s, phase_of(i, j)])
            elif a1:
                prods.append([points[i][4], a1, 0, 0, s, phase_of(i)])
            elif a2:
                prods.append([points[j][4], a2, 0, 0, s, phase_of(j)])
        for n, (i, s) in enumerate(singles):
            a1 = A[("s", n)]
            if a1:
                prods.append([points[i][4], a1, 0, 0, s, phase_of(i)])
        return prods


def search(an: core.Analysis, seed: int = 0, effort: int = 1, start=None, verbose=False, engine=None,
           max_levels: int = 10 ** 6, mix=None):
    """Raise the span codimensions one relation at a time; after each step repair the pairing (random swaps,
    then exact scans). Deterministic for given (table, k, seed, effort).
    Returns (engine at its best spans, best price, history)."""
    t0 = time.time()
    quick, trials, walk, max_fail = 60 * effort, 150 * effort, 0, 1 + effort
    if engine is None:
        eng = Engine(an, seed)
        evec = list(start) if start is not None else list(eng.bd.default_e())
        eng.init_spans(evec)
        cur, _ = eng.price()
        cur = eng.climb(100 * effort, cur, stop_at=(eng.r + 1) // 2)
    else:
        eng = engine
        cur, _ = eng.price()
    if mix is not None:
        eng.mix = mix
    best = (cur, eng.snapshot())
    hist = [(sum(eng.codim), eng.r, cur, eng.evals, eng.scans, round(time.time() - t0, 1))]
    if verbose:
        print("start", hist[-1], eng.codim, flush=True)
    fails = 0
    levels = 0
    while fails < max_fail and levels < max_levels:
        base = eng.snapshot()
        base_price = cur
        p, perfect = eng.level_up(cur, quick=quick, trials=trials, pick=3 * fails, walk=walk, verbose=verbose)
        hist.append((sum(eng.codim), eng.r, p, eng.evals, eng.scans, round(time.time() - t0, 1)))
        if verbose:
            print("level", hist[-1], eng.codim, "perfect" if perfect else "deficient", flush=True)
        if sum(eng.codim) == sum(len(R) for R in base):
            break                           # nothing could be added
        if p < best[0] or (p == best[0] and perfect):
            best = (p, eng.snapshot())
        if perfect:
            fails = 0
            cur = p
            levels += 1
        else:                               # back to the last perfect spans, try another relation
            fails += 1
            eng.restore(base)
            cur = base_price
    eng.restore(best[1])
    return eng, best[0], hist
