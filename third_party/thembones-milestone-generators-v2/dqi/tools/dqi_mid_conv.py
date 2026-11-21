#!/usr/bin/env python3
"""dqi_mid_conv.py -- word-level classical model of the kit's converter (workstream dqi_mid, items D1 and D3).

Transcribed, function by function, from the kit the 12-bit line builds with (read-only copy ~/qcc/ws/dqi_x12c/kit):
  final_decoder/conversion.py          convert_safe(monic_root=True) -> _bez_monic -> look_root_cached
  exp/convert_short_hankel/kernel.py   the convert_safe skeleton (square roots, top reduction, swaps, short Hankel)
  exp/safe_scalar/converter.py         look_solve, pre, post_look, q_post, base, canonical_permute, masked_product
  exp/convert/repair.py                rotate
  exp/convert/kernels.py               tri, hankel, convolution_xor, vandermonde, cross_terms
Every word is a numpy array over shots: the control flow does not depend on the data, and a swap under a quantum
flag is an elementwise select, so one call runs a whole batch.

Modes:
  'monic'    the kit as built: B monic (B_35 = 1 implicit); the root quotient q1 = 1/a lives in its own cached
             12-qubit word through the whole recursion.
  'uncached' same values; the root quotient is not cached.  Its only use inside the recursion is the root's
             child's q_post (the borrowed quotient of a degree gap), where it is rebuilt as 1/a_child (in a root gap
             the root's leader has moved into the child's g[-1]); the model checks on every shot that the rebuilt
             word equals the cached one wherever the gap flag uses it.  The root's post rebuilds 1/a itself.
  'gauge'    D1 as designed (SIG + converter scale in the gauge): the pair is lam times a monic pair and B's leading
             coefficient lam is explicit in the gauge word G.  After the square root G holds O's leader; the root
             runs the kit's generic level with f[-1] = G (q1 = G/a written into G: no cached word).  The kit uses
             the implicit 1 in three more places, which get their factor here: the top and E_17 reductions (one
             inverse of G before the Euclid) and the Hankel diagonal (34 leaf scalings, one inverse of G after).
Counters (static, per call): 'prod' (42 Toffolis: mul_xor, scale, convolution and Hankel rank updates; 'scale' is
the subset that are in-place scalings), inverses ('inv' at the levels, 'qinv' in q_post, 'ginv' of the gauge,
'rinv' rebuilt root quotients), 'zflag' zero tests, 'cswap' word swaps under a flag, 'mprod' masked products.
Usage: python dqi_mid_conv.py selftest [shots=512] [seed=1]
"""
import sys
from collections import Counter
from contextlib import contextmanager
from functools import lru_cache

import numpy as np

from dqi_x12c_gf import I64, EXP, LOG, ORD, mul, inv0

PRICE = dict(prod=42, inv=65, qinv=65, ginv=65, rinv=65, zflag=11, cswap=12, mprod=66)


# ---- scalar field helpers (the kit's conventions) -------------------------------------------------------------
def spow(a, n):
    """The kit's power(): power(0, 0) = 1."""
    if a == 0:
        return int(n == 0)
    return int(EXP[(LOG[a] * n) % ORD])


def smul(a, b):
    return int(EXP[LOG[a] + LOG[b]]) if a and b else 0


def sinv(a):
    assert a
    return int(EXP[(ORD - LOG[a]) % ORD])


def vpow(x, n):
    x = np.asarray(x, dtype=I64)
    return np.where(x == 0, 0, EXP[(LOG[x] * n) % ORD])


def mat_inv(a):
    """The kit's mi(): Gauss-Jordan inverse over GF(2^12)."""
    n = len(a)
    a = [list(r) + [int(i == j) for j in range(n)] for i, r in enumerate(a)]
    for j in range(n):
        p = next(i for i in range(j, n) if a[i][j])
        a[p], a[j] = a[j], a[p]
        t = sinv(a[j][j])
        a[j] = [smul(t, v) for v in a[j]]
        for i in range(n):
            if i != j and a[i][j]:
                t = a[i][j]
                a[i] = [x ^ smul(t, y) for x, y in zip(a[i], a[j])]
    return [r[n:] for r in a]


@lru_cache(None)
def vandermonde(n):
    V = [[spow(t, j) for j in range(n)] for t in range(n)]
    return V, mat_inv(V)


@lru_cache(None)
def cross_terms(k, m):
    n = k + m - 1
    V, Wm = vandermonde(n)
    return [(V[t][:k][::-1], [Wm[j][t] for j in range(n)], V[t][:m]) for t in range(n)]


class Wd:
    """One 12-bit word (a numpy array over shots)."""
    __slots__ = ('v',)

    def __init__(self, v):
        self.v = np.array(v, dtype=I64)


def words(arr):
    """(shots, n) array -> list of n words."""
    return [Wd(arr[:, i]) for i in range(arr.shape[1])]


def values(ws):
    return np.stack([w.v for w in ws], axis=1)


# ---- the machine -------------------------------------------------------------------------------------------
class Conv:
    def __init__(self, n, mode='monic'):
        assert mode in ('monic', 'uncached', 'gauge')
        self.n = n
        self.mode = mode
        self.c = Counter()
        self.checks = Counter()
        self.bad_rows = np.zeros(n, dtype=bool)        # 'uncached': shots where a gap used a wrong rebuilt word

    # primitive operations
    def lin(self, coeffs, ws):
        acc = np.zeros(self.n, dtype=I64)
        for c, w in zip(coeffs, ws):
            if c:
                acc ^= mul(c, w.v)
        return acc

    def mul_xor(self, x, y, out):
        out.v = out.v ^ mul(x.v, y.v)
        self.c['prod'] += 1

    def scale(self, x, a, ai, divide):
        x.v = mul(x.v, ai.v if divide else a.v)
        self.c['prod'] += 1
        self.c['scale'] += 1

    @staticmethod
    def swap(x, y):
        x.v, y.v = y.v, x.v

    def cswap(self, F, x, y):
        xv = np.where(F, y.v, x.v)
        y.v = np.where(F, x.v, y.v)
        x.v = xv
        self.c['cswap'] += 1

    def zflag(self, w):
        self.c['zflag'] += 1
        return w.v == 0

    def rank_update(self, aw, ac, bw, bc, tw, tc):
        p = mul(self.lin(ac, aw), self.lin(bc, bw))
        self.c['prod'] += 1
        for c, t in zip(tc, tw):
            if c:
                t.v = t.v ^ mul(c, p)

    def convolution_xor(self, a, c, out, reverse=False):
        N = len(a) + len(c) - 1
        V, Wm = vandermonde(N)
        for t in (range(N - 1, -1, -1) if reverse else range(N)):
            tc = [Wm[j][t] for j in range(len(out))]
            if any(tc):
                self.rank_update(a, V[t][:len(a)], c, V[t][:len(c)], out, tc)

    def masked_product(self, F, x, c, target):
        target.v = target.v ^ np.where(F, mul(x.v, c.v), 0)
        self.c['mprod'] += 1

    @contextmanager
    def inverse(self, a, tag='inv'):
        """Zero-safe inverse context (inverse_apply_safe / q_inverse): inside it a reads 1 where it is 0, ai = 1/a."""
        self.c[tag] += 1
        Z = a.v == 0
        a.v = np.where(Z, 1, a.v)
        ai = Wd(inv0(a.v))
        yield ai
        assert np.array_equal(ai.v, inv0(a.v)), 'inverse context: a changed inside its own context'
        a.v = np.where(Z, 0, a.v)

    def runops(self, ops, a, ai, undo=False):
        for op in (reversed(ops) if undo else ops):
            if op[0] == 's':
                self.scale(op[1], a, ai, op[2] ^ undo)
            elif op[0] == 'm':
                self.mul_xor(op[1], op[2], op[3])
            elif op[0] == 'c':
                self.convolution_xor(op[1], op[2], op[3], undo)
            else:
                raise ValueError(op[0])

    def frob(self, ws, power):
        for w in ws:
            w.v = vpow(w.v, power)

    # exp/convert/repair.py rotate, exp/safe_scalar canonical_permute
    def rotate(self, v, F, reverse=False):
        for j in (range(1, len(v)) if reverse else range(len(v) - 1, 0, -1)):
            self.cswap(F, v[j], v[j - 1])

    def canonical_permute(self, f, g, F, reverse=False):
        cyc = []
        for i in range(len(g) - 1):
            cyc += [g[i], f[i]]
        cyc += [g[-1]]
        for i in (range(len(cyc) - 2, -1, -1) if reverse else range(len(cyc) - 1)):
            self.cswap(F, cyc[i], cyc[i + 1])

    # exp/safe_scalar/converter.py
    def pre(self, f, g, payloads, undo=False):
        m = len(f) - 1
        n = m - 1
        a = g[-1]
        q0, q1 = f[-2], f[-1]
        with self.inverse(a) as ai:
            ops = [('s', v[-1], True) for v in payloads]
            ops += [('m', g[i], v[-1], v[i]) for v in payloads for i in range(n)]
            ops += [('s', q1, True), ('m', q1, g[-2], q0), ('s', q0, True), ('c', g[:-1], [q0, q1], f[:-2])]
            self.runops(ops, a, ai, undo)

    def q_post(self, q, payloads, parent, reverse=False, a=None):
        if parent is None:
            with self.inverse(q, 'qinv') as qi:
                for v in payloads:
                    self.scale(v[-1], q, qi, not reverse)
            return
        tops, p, c = parent
        rebuilt = None
        if isinstance(p, tuple):                       # 'uncached': the root quotient rebuilt from this level's a
            shadow = p[1]
            p = Wd(inv0(np.where(a.v == 0, 1, a.v)))
            self.c['rinv'] += 2                        # rebuilt before the first and before the last swap
            rebuilt = (p.v.copy(), shadow)
        F = self.zflag(q)
        if rebuilt is not None:
            used = F
            self.checks['root quotient used by a gap'] += int(used.sum())
            wrong = used & (rebuilt[0] != rebuilt[1].v)
            self.checks['rebuilt root quotient wrong'] += int(wrong.sum())
            self.bad_rows |= wrong
        self.cswap(F, q, p)
        if not reverse:
            for v, t in zip(payloads, tops):
                self.cswap(F, v[-1], t)
        with self.inverse(q, 'qinv') as qi:
            if not reverse:
                for v in payloads:
                    self.scale(v[-1], q, qi, True)
                for v, t in zip(payloads, tops):
                    self.masked_product(F, v[-1], c, t)
            else:
                for v, t in zip(payloads, tops):
                    self.masked_product(F, v[-1], c, t)
                for v in payloads:
                    self.scale(v[-1], q, qi, False)
        if reverse:
            for v, t in zip(payloads, tops):
                self.cswap(F, v[-1], t)
        self.cswap(F, q, p)
        if rebuilt is not None:                        # the borrowed word must be back before its erasure
            assert np.array_equal(p.v, rebuilt[0])
        assert np.array_equal(q.v == 0, F), 'q_post: flag not restored'

    def post_look(self, f, g, payloads, parent, reverse=False):
        m = len(f) - 1
        n = m - 1
        a = g[-1]
        q0, q1 = f[-2], f[-1]

        def dots():
            for v in payloads:
                for i in range(n):
                    self.mul_xor(g[i], v[i], v[-1])

        def last(undo):
            with self.inverse(a) as ai:
                ops = [('s', v[-1], True) for v in payloads]
                ops += [('c', g[:-1], [q0, q1], f[:-2]), ('s', q0, False), ('m', q1, g[-2], q0), ('s', q1, False)]
                self.runops(ops, a, ai, undo)
        if not reverse:
            self.q_post(q1, payloads, parent, a=a)
            dots()
            last(False)
        else:
            last(True)
            dots()
            self.q_post(q1, payloads, parent, True, a=a)

    def base(self, f, g, payloads, parent, reverse=False):
        a = g[0]
        q = f[-1]

        def first(undo):
            with self.inverse(a) as ai:
                self.runops([('s', v[0], True) for v in payloads] + [('s', q, True)], a, ai, undo)

        def last(undo):
            with self.inverse(a) as ai:
                self.runops([('s', v[0], True) for v in payloads] + [('s', q, False)], a, ai, undo)
        if not reverse:
            first(False)
            self.q_post(q, payloads, parent, a=a)
            last(False)
        else:
            last(True)
            self.q_post(q, payloads, parent, True, a=a)
            first(True)

    def look_solve(self, f, g, payloads, parent=None, reverse=False):
        m = len(f) - 1
        if m == 1:
            return self.base(f, g, payloads, parent, reverse)
        F = self.zflag(g[-1])
        self.rotate(g, F)
        if not reverse:
            self.pre(f, g, payloads)
        else:
            self.post_look(f, g, payloads, parent, True)
        self.canonical_permute(f, g, F)
        assert np.array_equal(g[-1].v == 0, F), 'look_solve: flag after the permutation'
        child = ([v[-1] for v in payloads], f[-1], f[-2])
        self.look_solve(g, f[:-2], [v[:-1] for v in payloads], child, reverse)
        F2 = self.zflag(g[-1])
        assert np.array_equal(F2, F)
        self.canonical_permute(f, g, F2, True)
        if not reverse:
            self.post_look(f, g, payloads, parent)
        else:
            self.pre(f, g, payloads, True)
        self.rotate(g, F2, True)
        assert np.array_equal(g[-1].v == 0, F2), 'look_solve: flag at the end'

    # final_decoder/conversion.py: the cached root (and its 'uncached' twin)
    def cached_root_action(self, f, g, payloads, ai, Z, post, undo):
        a = g[-1]
        q0, q1 = f[-2], f[-1]
        n = len(g) - 1
        if post:
            ops = [('s', v[-1], False) for v in payloads]
            ops += [('m', g[i], v[i], v[-1]) for v in payloads for i in range(n)]
            ops += [('s', v[-1], True) for v in payloads]
            ops += [('c', g[:-1], [q0, q1], f[:-2]), ('s', q0, False), ('m', q1, g[-2], q0)]
        else:
            ops = [('s', v[-1], True) for v in payloads]
            ops += [('m', g[i], v[-1], v[i]) for v in payloads for i in range(n)]
            ops += [('m', q1, g[-2], q0), ('s', q0, True), ('c', g[:-1], [q0, q1], f[:-2])]
        a.v = np.where(Z, 1, a.v)                      # cx(Z, a[0]): a reads 1 where it is 0
        self.runops(ops, a, ai, undo)
        a.v = np.where(Z, 0, a.v)

    def look_root_cached(self, O, g, payloads, reverse=False):
        cached = self.mode == 'monic'
        F = self.zflag(g[-1])
        self.rotate(g, F)
        a = g[-1]
        self.c['inv'] += 1                             # inverse_zero(a): ai and Z = [a = 0]
        Z = a.v == 0
        ai = Wd(np.where(Z, 1, inv0(a.v)))             # cx(Z, ai[0])
        f = list(O) + [ai]
        self.cached_root_action(f, g, payloads, ai, Z, post=reverse, undo=reverse)
        self.canonical_permute(f, g, F)
        p = f[-1] if cached else ('rebuild', Wd(ai.v.copy()))
        self.look_solve(g, f[:-2], [v[:-1] for v in payloads], ([v[-1] for v in payloads], p, f[-2]), reverse)
        F2 = self.zflag(g[-1])
        assert np.array_equal(F2, F)
        self.canonical_permute(f, g, F2, True)
        Z2 = self.zflag(a)
        assert np.array_equal(Z2, Z)
        if not cached:
            self.c['inv'] += 1                         # the root's post rebuilds 1/a (its mask came from the pre)
        assert np.array_equal(ai.v, np.where(Z, 1, inv0(a.v)))
        self.cached_root_action(f, g, payloads, ai, Z2, post=not reverse, undo=reverse)
        self.rotate(g, F2, True)

    # exp/convert/kernels.py tri / hankel, with an optional diagonal (the gauge)
    def tri(self, q, x, reverse=False, g0=None, g0i=None):
        n = len(x)
        if n == 1:
            if g0 is not None:
                self.scale(x[0], g0, g0i, reverse)     # forward: x0 *= g, reverse: x0 /= g
            return
        k = n // 2
        m = n - k
        if not reverse:
            self.tri(q[:m - 1], x[k:], False, g0, g0i)
        else:
            self.tri(q[:k - 1], x[:k], True, g0, g0i)
        ts = cross_terms(k, m)
        for ac, bc, tc in (reversed(ts) if reverse else ts):
            self.rank_update(x[:k], ac, q[:n - 1], bc, x[k:], tc)
        if not reverse:
            self.tri(q[:k - 1], x[:k], False, g0, g0i)
        else:
            self.tri(q[:m - 1], x[k:], True, g0, g0i)

    def hankel(self, Bs, As, reverse=False, g0=None, g0i=None):
        assert len(Bs) == len(As)
        q = list(reversed(Bs[1:]))
        if reverse:
            for i in range(len(As) // 2):
                self.swap(As[i], As[-1 - i])
            self.tri(q, As, True, g0, g0i)
        else:
            self.tri(q, As, False, g0, g0i)
            for i in range(len(As) // 2):
                self.swap(As[i], As[-1 - i])

    def bez(self, O, E, payloads, reverse, G=None):
        def reduce():
            for i in range(len(O)):
                self.mul_xor(E[-1], O[i], E[i])
        reduce()
        if self.mode == 'gauge':
            self.look_solve(O + [G], E[:-1], payloads, None, reverse)
        else:
            self.look_root_cached(O, E[:-1], payloads, reverse)
        reduce()

    # exp/convert_short_hankel/kernel.py convert_safe (monic_root=True), and the gauge variant
    def convert(self, A, B, reverse=False, G=None):
        assert len(A) == len(B) == 35
        gauge = self.mode == 'gauge'
        assert gauge == (G is not None)
        allw = A + B + ([G] if gauge else [])
        self.frob(allw, 2048)
        E, O, u, w, top = B[::2], B[1::2], A[1::2], A[:34:2], A[-1]
        if not reverse:
            if gauge:
                with self.inverse(G, 'ginv') as gi:
                    self.scale(top, G, gi, True)
                    self.scale(E[-1], G, gi, True)
            for i in range(17):
                self.mul_xor(top, O[i], w[i])
            self.bez(O, E, [u, w], False, G)
            if gauge:
                with self.inverse(G, 'ginv') as gi:
                    self.scale(E[-1], G, gi, False)
                    for i in range(17):
                        self.mul_xor(E[i + 1], w[i], top)
                    for i in range(17):
                        self.swap(A[2 * i], A[2 * i + 1])
                    self.hankel(B[1:], A[:34], False, G, gi)
            else:
                for i in range(17):
                    self.mul_xor(E[i + 1], w[i], top)
                for i in range(17):
                    self.swap(A[2 * i], A[2 * i + 1])
                self.hankel(B[1:], A[:34])
            for i in reversed(range(34)):
                self.swap(A[i], A[i + 1])
        else:
            for i in range(34):
                self.swap(A[i], A[i + 1])
            if gauge:
                with self.inverse(G, 'ginv') as gi:
                    self.hankel(B[1:], A[:34], True, G, gi)
                    for i in range(17):
                        self.swap(A[2 * i], A[2 * i + 1])
                    for i in range(17):
                        self.mul_xor(E[i + 1], w[i], top)
                    self.scale(E[-1], G, gi, True)
            else:
                self.hankel(B[1:], A[:34], True)
                for i in range(17):
                    self.swap(A[2 * i], A[2 * i + 1])
                for i in range(17):
                    self.mul_xor(E[i + 1], w[i], top)
            self.bez(O, E, [u, w], True, G)
            for i in range(17):
                self.mul_xor(top, O[i], w[i])
            if gauge:
                with self.inverse(G, 'ginv') as gi:
                    self.scale(E[-1], G, gi, False)
                    self.scale(top, G, gi, False)
        self.frob(allw, 2)

    def toffolis(self):
        return sum(PRICE[k] * v for k, v in self.c.items() if k in PRICE)


# ---- polynomial helpers for the checks (arrays (shots, coefficients), lowest degree first) -------------------
def pmul_v(a, b):
    out = np.zeros((a.shape[0], a.shape[1] + b.shape[1] - 1), dtype=I64)
    for i in range(a.shape[1]):
        out[:, i:i + b.shape[1]] ^= mul(a[:, i:i + 1], b)
    return out


def pmod_monic_v(r, B):
    """r mod B for B monic of degree d = B.shape[1] - 1 (B[:, d] == 1)."""
    r = r.copy()
    d = B.shape[1] - 1
    for k in range(r.shape[1] - 1, d - 1, -1):
        t = r[:, k].copy()
        r[:, k - d:k + 1] ^= mul(t[:, None], B)
    return r[:, :d]


def deriv_v(P):
    D = np.zeros_like(P)
    D[:, 0:-1:2] = P[:, 1::2]
    return D


def residue_ok(A, Bfull, C):
    """C * B' == A mod B (B monic, 36 coefficients; A, C 35)."""
    lhs = pmod_monic_v(pmul_v(C, deriv_v(Bfull)), Bfull)
    return (lhs == A[:, :35]).all(axis=1)


def coprime_halves(Bfull):
    """gcd(E, O) == 1 per shot (scalar Euclid; slow, used only in the self-test)."""
    from dqi_x12c_gf import smul as sm, sinv as si
    out = np.zeros(Bfull.shape[0], dtype=bool)
    for s in range(Bfull.shape[0]):
        b = [int(x) for x in vpow(Bfull[s], 2048)]
        E = b[0::2]
        O = b[1::2]

        def trim(p):
            while p and p[-1] == 0:
                p = p[:-1]
            return p
        f, g = trim(O), trim(E)
        while g:
            while len(f) >= len(g):
                c = sm(f[-1], si(g[-1]))
                sh = len(f) - len(g)
                for i, y in enumerate(g):
                    f[i + sh] ^= sm(c, y)
                f = trim(f)
                if not f:
                    break
            f, g = g, f
        out[s] = len(f) == 1
    return out


# ---- self-test ------------------------------------------------------------------------------------------------
def random_inputs(n, rng, root_gap=False):
    A = rng.integers(0, 4096, size=(n, 35), dtype=I64)
    if not root_gap:
        Bl = rng.integers(0, 4096, size=(n, 35), dtype=I64)
        return A, Bl
    # B = E^2 + z O^2 with O monic of degree 17 and E mod O of degree <= 15 (a zero leader at the root)
    O = rng.integers(0, 4096, size=(n, 18), dtype=I64)
    O[:, 17] = 1
    R = rng.integers(0, 4096, size=(n, 17), dtype=I64)
    R[:, 16] = 0
    e17 = rng.integers(0, 4096, size=n, dtype=I64)
    E = mul(e17[:, None], O)
    E[:, :17] ^= R
    Bfull = np.zeros((n, 36), dtype=I64)
    Bfull[:, 0::2] = vpow(E, 2)
    Bfull[:, 1::2] = vpow(O, 2)
    return A, Bfull[:, :35]


def run_mode(mode, A, Bl, lam=None):
    n = A.shape[0]
    cv = Conv(n, mode)
    if mode == 'gauge':
        Aw, Bw, G = words(mul(lam[:, None], A)), words(mul(lam[:, None], Bl)), Wd(lam)
    else:
        Aw, Bw, G = words(A), words(Bl), None
    cv.convert(Aw, Bw, False, G)
    fwd = dict(cv.c)
    C = values(Aw)
    Bout = values(Bw)
    Gout = None if G is None else G.v.copy()
    cv.c = Counter()
    cv.convert(Aw, Bw, True, G)
    rev = dict(cv.c)
    back = (values(Aw), values(Bw), None if G is None else G.v.copy())
    return C, Bout, Gout, back, fwd, rev, cv.checks


def selftest(n, seed):
    rng = np.random.default_rng(seed)
    for label, gap in (('random monic B', False), ('forced root gap (E mod O of degree <= 15)', True)):
        A, Bl = random_inputs(n, rng, gap)
        Bfull = np.concatenate([Bl, np.ones((n, 1), dtype=I64)], axis=1)
        cop = coprime_halves(Bfull)
        lam = rng.integers(1, 4096, size=n, dtype=I64)
        res = {}
        for mode in ('monic', 'uncached', 'gauge'):
            res[mode] = run_mode(mode, A, Bl, lam if mode == 'gauge' else None)
        C0 = res['monic'][0]
        ok_true = residue_ok(A, Bfull, C0)
        print('== %s: %d shots, gcd(E, O) = 1 in %d' % (label, n, cop.sum()))
        print('   monic: C * B\' = A mod B on %d of %d coprime shots (all shots %d)'
              % ((ok_true & cop).sum(), cop.sum(), ok_true.sum()))
        for mode in ('monic', 'uncached', 'gauge'):
            C, Bout, Gout, back, fwd, rev, checks = res[mode]
            if mode == 'gauge':
                inA, inB = mul(lam[:, None], A), mul(lam[:, None], Bl)
                same_in = (back[0] == inA).all(axis=1) & (back[1] == inB).all(axis=1) & (back[2] == lam)
                bkeep = (Bout == inB).all(axis=1) & (Gout == lam)
            else:
                same_in = (back[0] == A).all(axis=1) & (back[1] == Bl).all(axis=1)
                bkeep = (Bout == Bl).all(axis=1)
            sameC = (C == C0).all(axis=1)
            tof_f = sum(PRICE[k] * v for k, v in fwd.items() if k in PRICE)
            tof_r = sum(PRICE[k] * v for k, v in rev.items() if k in PRICE)
            print('   %-8s C equal to monic %d/%d; B (and G) kept %d; round trip exact %d; checks %s'
                  % (mode, sameC.sum(), n, bkeep.sum(), same_in.sum(), dict(checks)))
            print('            forward %s -> %d Toffolis; reverse %s -> %d'
                  % ({k: fwd[k] for k in sorted(fwd)}, tof_f, {k: rev[k] for k in sorted(rev)}, tof_r))


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'selftest'
    shots = int(sys.argv[2]) if len(sys.argv) > 2 else 512
    sd = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    if cmd == 'selftest':
        selftest(shots, sd)
