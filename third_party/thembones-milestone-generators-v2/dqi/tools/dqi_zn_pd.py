#!/usr/bin/env python3
"""dqi_zn_pd.py -- builder Z1 (Nov 14): word-level model of the converter's half-Euclid (bez step) with PAIRED
payload digits (R33-DQI-CONV), against an exact word-level model of the kit's own bez step.

K-model = the kit's exp/safe_scalar/converter.py bez_lookahead, forward direction, word by word: reduce, borrowed
lead word = 1, look_solve recursion (zero flag of the pivot, rotate, pre with the zero-pivot clip a := 1, canonical
permute, child, reverse permute, post_look = q_post with the parent gap route + dots + last), base level, reduce.
Every product / scaling / conditional word operation is counted by role (pay = payload, loc = locator).
Unconvert is the gate-by-gate inverse of convert, so the forward map decides values and the product counts.

Registers are numpy int64 arrays over shots (one GF(2^12) word per array); flags are boolean arrays.
"""
import collections
import os
import sys

import numpy as np

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
from dqi_x12c_gf import mul, inv0  # noqa: E402  (modulus 0x1053, the kit's field)

I64 = np.int64


class Mach:
    """Word registers + counters. Handles are ints; lists of handles stand for the kit's word lists."""

    def __init__(self, n):
        self.n = n
        self.R = []
        self.cnt = collections.Counter()

    def new(self, val=None):
        self.R.append(np.zeros(self.n, dtype=I64) if val is None else np.array(val, dtype=I64))
        return len(self.R) - 1

    def zero(self, h):
        return self.R[h] == 0

    # ---- counted operations (role, kind)
    def mx(self, x, y, t, role):
        """t ^= x * y (one product, 42)."""
        self.R[t] = self.R[t] ^ mul(self.R[x], self.R[y])
        self.cnt[(role, 'prod')] += 1

    def scl(self, x, s, role):
        """x := x * s for a value array s (one in-place scaling)."""
        self.R[x] = mul(self.R[x], s)
        self.cnt[(role, 'scale')] += 1

    def cswap(self, F, x, y, role):
        a, b = self.R[x], self.R[y]
        self.R[x], self.R[y] = np.where(F, b, a), np.where(F, a, b)
        self.cnt[(role, 'cswap')] += 1

    def mmx(self, F, x, y, t, role):
        """t ^= F x y (masked product: 12 ANDs + product + 12 CCZ)."""
        self.R[t] = self.R[t] ^ np.where(F, mul(self.R[x], self.R[y]), 0)
        self.cnt[(role, 'mprod')] += 1

    def conv(self, a, c, out, role):
        """out[i] ^= coefficient i of (a * c), i < len(out) (the kit's convolution_xor: len(a)+len(c)-1 points)."""
        R = self.R
        for i in range(len(out)):
            acc = R[out[i]]
            for j in range(len(c)):
                k = i - j
                if 0 <= k < len(a):
                    acc = acc ^ mul(R[a[k]], R[c[j]])
            R[out[i]] = acc
        self.cnt[(role, 'conv_pts')] += len(a) + len(c) - 1

    def clip(self, h):
        """Zero-pivot clip of the kit's safe inverses: during the action the word holds 1 if it was 0.
        Returns (F, a_eff, ai_eff); call unclip(h, F) after the action."""
        F = self.R[h] == 0
        self.R[h] = np.where(F, 1, self.R[h])
        self.cnt[('loc', 'inverse')] += 1
        return F, self.R[h].copy(), inv0(self.R[h])

    def unclip(self, h, F):
        self.R[h] = np.where(F, 0, self.R[h])


def rotate(M, v, F, reverse=False):
    """repair.rotate: conditioned on F, v := w * v cyclically (the zero top word goes to v[0])."""
    rng = range(1, len(v)) if reverse else range(len(v) - 1, 0, -1)
    for j in rng:
        M.cswap(F, v[j], v[j - 1], 'route')


def canonical_permute(M, f, g, F, reverse=False):
    cyc = []
    for i in range(len(g) - 1):
        cyc += [g[i], f[i]]
    cyc += [g[-1]]
    seq = range(len(cyc) - 2, -1, -1) if reverse else range(len(cyc) - 1)
    for i in seq:
        M.cswap(F, cyc[i], cyc[i + 1], 'route')


# ------------------------------------------------------------------------------------------- K-model (the kit)
class Kit:
    """The kit's bez step, forward. Hooks (pre_pay, post_pay) let a variant replace the payload operations."""

    def __init__(self, M):
        self.M = M
        self.trace = []          # (level, event, payload snapshot) when tracing

    # payload part of pre (inside the pivot's inverse context): digit and rank-one updates
    def pre_pay(self, m, f, g, P, a, ai, parent=None):
        M = self.M
        for v in P:
            M.scl(v[-1], ai, 'pay')
        for v in P:
            for i in range(m - 1):
                M.mx(g[i], v[-1], v[i], 'pay')

    def dots(self, m, f, g, P):
        M = self.M
        for v in P:
            for i in range(m - 1):
                M.mx(g[i], v[i], v[-1], 'pay')

    def pre(self, f, g, P, parent=None):
        M = self.M
        m = len(f) - 1
        F0, a, ai = M.clip(g[-1])
        self.pre_pay(m, f, g, P, a, ai, parent)
        q0, q1 = f[-2], f[-1]
        M.scl(q1, ai, 'loc')
        M.mx(q1, g[-2], q0, 'loc')
        M.scl(q0, ai, 'loc')
        M.conv(g[:-1], [q0, q1], f[:-2], 'loc')
        M.unclip(g[-1], F0)

    def last(self, f, g, P):
        M = self.M
        F0, a, ai = M.clip(g[-1])
        for v in P:
            M.scl(v[-1], ai, 'pay')
        q0, q1 = f[-2], f[-1]
        M.conv(g[:-1], [q0, q1], f[:-2], 'loc')
        M.scl(q0, a, 'loc')
        M.mx(q1, g[-2], q0, 'loc')
        M.scl(q1, a, 'loc')
        M.unclip(g[-1], F0)

    def q_post(self, q, P, m, parent):
        M = self.M
        if parent is None:
            F0, a, ai = M.clip(q)
            for v in P:
                M.scl(v[-1], ai, 'pay')
            M.unclip(q, F0)
            return
        tops, p, c = parent
        F = M.zero(q)
        M.cswap(F, q, p, 'route')
        for v, t in zip(P, tops):
            M.cswap(F, v[-1], t, 'route')
        F0, a, ai = M.clip(q)
        for v in P:
            M.scl(v[-1], ai, 'pay')
        for v, t in zip(P, tops):
            M.mmx(F, v[-1], c, t, 'route')
        M.unclip(q, F0)
        M.cswap(F, q, p, 'route')
        assert (M.zero(q) == F).all()

    def base(self, f, g, P, parent):
        M = self.M
        q = f[-1]
        F0, a, ai = M.clip(g[0])
        for v in P:
            M.scl(v[0], ai, 'pay')
        M.scl(q, ai, 'loc')
        M.unclip(g[0], F0)
        self.at_bottom(f, g, P)
        self.snap_bottom()
        self.q_post(q, P, 1, parent)
        F0, a, ai = M.clip(g[0])
        for v in P:
            M.scl(v[0], ai, 'pay')
        M.scl(q, a, 'loc')
        M.unclip(g[0], F0)

    def at_bottom(self, f, g, P):
        pass

    def snap_bottom(self):
        """Payload words (17 per payload) and every level's flag, right after the conversion at the bottom."""
        self.bottom = np.stack([self.M.R[h] for v in self.P0 for h in v], 1)

    def post_look(self, f, g, P, parent):
        m = len(f) - 1
        self.q_post(f[-1], P, m, parent)
        self.dots(m, f, g, P)
        self.last(f, g, P)

    def look(self, f, g, P, parent=None):
        M = self.M
        m = len(f) - 1
        if m == 1:
            return self.base(f, g, P, parent)
        F = M.zero(g[-1])
        rotate(M, g, F)
        self.lv[m] = dict(f=f, g=g, P=P, parent=parent, F=F)
        self.pre(f, g, P, parent)
        canonical_permute(M, f, g, F)
        assert (M.zero(g[-1]) == F).all()
        self.flags[m] = F
        child = ([v[-1] for v in P], f[-1], f[-2])
        self.look(g, f[:-2], [v[:-1] for v in P], child)
        F2 = M.zero(g[-1])
        assert (F2 == F).all()
        canonical_permute(M, f, g, F, True)
        self.post_look(f, g, P, parent)
        rotate(M, g, F, True)
        assert (M.zero(g[-1]) == F).all()

    def bez(self, O, E, P):
        """O: 17 handles (monic, O_17 = 1 implicit), E: 18 handles, P: payload handle lists (17 each)."""
        M = self.M
        self.flags = {}
        self.lv = {}
        self.P0 = P
        for i in range(len(O)):
            M.mx(E[-1], O[i], E[i], 'loc')
        lead = M.new(np.ones(M.n, dtype=I64))
        self.look(O + [lead], E[:-1], P)
        assert (M.R[lead] == 1).all()
        for i in range(len(O)):
            M.mx(E[-1], O[i], E[i], 'loc')


# ------------------------------------------------------------------------------------------- P-model (paired)
UPPER = (17, 15, 13, 11, 9, 7, 5, 3)


class Paired(Kit):
    """Paired payload digits on the down sweep, the kit's return sweep unchanged.
    Block = (upper level U, lower level U - 1). At U: no payload operation. At U - 1, inside its pivot's inverse
    context (a = its divisor's top, clipped), with t = the upper level's top word (the parent's top):
      D1 = t / a;  D0 = (v[-1] + D1 g[-2]) / a;  rest v[i] ^= D0 g[i] + D1 g[i-1] for i < m - 1
    (point form: one product per rest word, plus the one for D0). Conversion at the bottom, while every quotient
    is alive: X = D1 / clip(q1'), Y = D0 + X q0' (q1', q0' = the lower level's quotient words);
    kit digits d_U = X + X(block above), d_L = Y; the base digit ^= X(last block)."""

    def __init__(self, M, upper=UPPER, route='generic'):
        super().__init__(M)
        self.upper = set(upper)
        self.lower = {u - 1 for u in upper}
        self.route = route

    def pre_pay(self, m, f, g, P, a, ai, parent=None):
        M = self.M
        if m in self.upper:
            return
        if m not in self.lower:
            return Kit.pre_pay(self, m, f, g, P, a, ai, parent)
        R = M.R
        for v, t in zip(P, parent[0]):
            M.scl(t, ai, 'pay')                       # D1
            M.mx(t, g[-2], v[-1], 'pay')
            M.scl(v[-1], ai, 'pay')                   # D0
            for i in range(m - 1):
                acc = mul(R[v[-1]], R[g[i]])
                if i >= 1:
                    acc = acc ^ mul(R[t], R[g[i - 1]])
                R[v[i]] = R[v[i]] ^ acc
            M.cnt[('pay', 'prod')] += m - 1           # point form: m - 1 rest words, one product each

    def at_bottom(self, f, g, P):
        if self.route == 'full':
            return self.at_bottom_full(f, g, P)
        M = self.M
        blocks = sorted(self.upper, reverse=True)
        for U in blocks:
            lu, ll = self.lv[U], self.lv[U - 1]
            q1, q0 = ll['f'][-1], ll['f'][-2]
            F0, qa, qi = M.clip(q1)
            for vu in lu['P']:
                M.scl(vu[-1], qi, 'pay')              # X = D1 / q1'
            M.unclip(q1, F0)
            for vu, vl in zip(lu['P'], ll['P']):
                M.mx(vu[-1], q0, vl[-1], 'pay')       # Y = D0 + X q0'
        # corrections, bottom-up (free XORs in the generic chart): X of block (U, U-1) goes into the digit of
        # level U - 2 (the next block's X word, a kit level's digit, or the base word)
        R = M.R
        for U in sorted(self.upper):
            T = U - 2
            tg = [v[0] for v in P] if T == 1 else [v[-1] for v in self.lv[T]['P']]
            for t, vu in zip(tg, self.lv[U]['P']):
                R[t] = R[t] ^ R[vu[-1]]


    def tops(self, m, P):
        """Digit words of level m (its payload top words), or the base words for m = 1."""
        return [v[0] for v in P] if m == 1 else [v[-1] for v in self.lv[m]['P']]

    def qflag(self, m):
        """[q1 of level m == 0] at the bottom (a zero test, 11 + 11 Toffolis). Equals F_(m+1), the rotation flag of
        level m + 1 (its divisor's top was zero), since level m's dividend top is that word."""
        self.M.cnt[('route', 'flag')] += 1
        return self.M.R[self.lv[m]['f'][-1]] == 0

    def at_bottom_full(self, f, g, P):
        """Conversion with exact gap routes (one gap per block at most; double gaps are outside the kit's domain).
        Block (U, L): B = F_U = [q1_L == 0] (gap at U: the block's digits ARE the kit's, no carry);
        C = F_L = [q1_(L-1) == 0] (gap at L: the carry enters the next pair of levels through f_(U-2)).
        Carry of block j: x = [not B] X. Into levels (U', L') = (U-2, U-3):
          after a block without C: x g_U' (entry) -> d_U' += x if not F_U', else d_L' += x;
          after a block with C:    x f_U'         -> d_U' += x q0_U' and d_L' += x."""
        M = self.M
        R = M.R
        blocks = sorted(self.upper, reverse=True)
        info = {}
        for U in blocks:
            L = U - 1
            lu, ll = self.lv[U], self.lv[L]
            q1, q0 = ll['f'][-1], ll['f'][-2]
            Bf = self.qflag(L)
            Cf = self.qflag(L - 1) if L - 1 >= 2 else np.zeros(M.n, dtype=bool)
            F0, qa, qi = M.clip(q1)
            for vu in lu['P']:
                M.scl(vu[-1], qi, 'pay')              # X = D1 / q1'  (D1 itself when B: q1' = 0 is clipped to 1)
            M.unclip(q1, F0)
            for vu, vl in zip(lu['P'], ll['P']):
                M.mmx(~Bf, vu[-1], q0, vl[-1], 'pay')  # Y = D0 + [not B] X q0'
            info[U] = (Bf, Cf)
        # carries, bottom-up
        for U in sorted(self.upper):
            Bf, Cf = info[U]
            Up, Lp = U - 2, U - 3
            xs = [np.where(Bf, 0, R[vu[-1]]) for vu in self.lv[U]['P']]
            assert Lp >= 2, 'full route: blocks with U >= 5'
            FU = self.qflag(Lp)
            tu, tl = self.tops(Up, P), self.tops(Lp, P)
            q0u = self.lv[Up]['f'][-2]
            for k, x in enumerate(xs):
                # no C: route by F_U'
                R[tu[k]] = R[tu[k]] ^ np.where(~Cf & ~FU, x, 0)
                R[tl[k]] = R[tl[k]] ^ np.where(~Cf & FU, x, 0)
                # C: x f_U' = x q0_U' g_U' + x g_L'
                R[tu[k]] = R[tu[k]] ^ np.where(Cf, mul(x, R[q0u]), 0)
                R[tl[k]] = R[tl[k]] ^ np.where(Cf, x, 0)
                M.cnt[('route', 'cxor')] += 3
                M.cnt[('route', 'mprod')] += 1


class KitQB(Kit):
    """The kit with every q_post moved to the bottom (levels 1..17 in the kit's order). Its return is then the
    transpose of the descent's payload operations: Ret = Down^T Q. Must equal Kit word for word."""

    def base(self, f, g, P, parent):
        M = self.M
        q = f[-1]
        F0, a, ai = M.clip(g[0])
        for v in P:
            M.scl(v[0], ai, 'pay')
        M.scl(q, ai, 'loc')
        M.unclip(g[0], F0)
        self.lv[1] = dict(f=f, g=g, P=P, parent=parent, F=None)
        self.at_bottom(f, g, P)
        self.snap_bottom()
        for m in range(1, 18):
            lv = self.lv[m]
            self.q_post(lv['f'][-1], lv['P'], m, lv['parent'])
        self.after_q(f, g, P)
        F0, a, ai = M.clip(g[0])
        for v in P:
            M.scl(v[0], ai, 'pay')
        M.scl(q, a, 'loc')
        M.unclip(g[0], F0)

    def after_q(self, f, g, P):
        pass

    def post_look(self, f, g, P, parent):
        m = len(f) - 1
        self.dots(m, f, g, P)
        self.last(f, g, P)


class PairedR(Paired, KitQB):
    """Paired descent (route 'full') and paired return: Ret = Down2^T Cv^T Q, Q = the kit's q_posts at the bottom.
    Cv^T: carries transposed (top-down), then per block t_U ^= [not B] slot(L) q0_L and t_U *= 1 / clip(q1_L).
    Down2^T at the lower level L (inside last's inverse context of a_L): slot(L) ^= <g_L low, v low> (S0),
    t_U ^= <(w g_L) low, v low> (S1): L - 1 point products for both; slot(L) *= 1 / a_L; t_U ^= slot(L) g_L[-2];
    t_U *= 1 / a_L. At the upper level U: no payload operation at all."""

    def __init__(self, M, upper=(17, 15, 13, 11), route='full'):
        Paired.__init__(self, M, upper=upper, route=route)

    base = KitQB.base

    def after_q(self, f, g, P):
        M = self.M
        R = M.R
        blocks = sorted(self.upper, reverse=True)
        info = {}
        for U in blocks:
            L = U - 1
            info[U] = (self.qflag(L), self.qflag(L - 1), self.qflag(U - 3))
        for U in blocks:                                  # carries transposed, top-down
            Bf, Cf, FU = info[U]
            Up, Lp = U - 2, U - 3
            tu, tl = self.tops(Up, P), self.tops(Lp, P)
            q0u = self.lv[Up]['f'][-2]
            for k, vu in enumerate(self.lv[U]['P']):
                src = (np.where(~Cf & ~FU, R[tu[k]], 0) ^ np.where((~Cf & FU) | Cf, R[tl[k]], 0)
                       ^ np.where(Cf, mul(R[q0u], R[tu[k]]), 0))
                R[vu[-1]] = R[vu[-1]] ^ np.where(Bf, 0, src)
                M.cnt[('route', 'cxor')] += 3
                M.cnt[('route', 'mprod')] += 1
        for U in blocks:                                  # conversion transposed
            Bf = info[U][0]
            ll, lu = self.lv[U - 1], self.lv[U]
            q1, q0 = ll['f'][-1], ll['f'][-2]
            for vu, vl in zip(lu['P'], ll['P']):
                M.mmx(~Bf, vl[-1], q0, vu[-1], 'pay')
            F0, qa, qi = M.clip(q1)
            for vu in lu['P']:
                M.scl(vu[-1], qi, 'pay')
            M.unclip(q1, F0)

    def post_look(self, f, g, P, parent):
        m = len(f) - 1
        if m not in self.upper and m not in self.lower:
            return KitQB.post_look(self, f, g, P, parent)
        M = self.M
        R = M.R
        F0, a, ai = M.clip(g[-1])
        if m in self.lower:
            for v, t in zip(P, parent[0]):
                s0 = np.zeros(M.n, dtype=I64)
                s1 = np.zeros(M.n, dtype=I64)
                for i in range(m - 1):
                    s0 = s0 ^ mul(R[g[i]], R[v[i]])
                    if i >= 1:
                        s1 = s1 ^ mul(R[g[i - 1]], R[v[i]])
                R[v[-1]] = R[v[-1]] ^ s0
                R[t] = R[t] ^ s1
                M.cnt[('pay', 'prod')] += m - 1          # transposed point form
                M.scl(v[-1], ai, 'pay')
                M.mx(v[-1], g[-2], t, 'pay')
                M.scl(t, ai, 'pay')
        q0, q1 = f[-2], f[-1]
        M.conv(g[:-1], [q0, q1], f[:-2], 'loc')
        M.scl(q0, a, 'loc')
        M.mx(q1, g[-2], q0, 'loc')
        M.scl(q1, a, 'loc')
        M.unclip(g[-1], F0)


class PairedR2(PairedR):
    """pd3: no extra inverse, no extra scaling. Descent: at U the kit's own scaling only (t := t / a_U = X, no rank
    updates); at L (locator part of pre first, so q1', q0' exist): B = [q1' == 0], c = q1' + B (clip),
    Z0 = [not B] q0'; per payload D1 = X c into a temp, slot(L) ^= D1 g[-2], slot(L) *= 1 / a_L (= D0), rest
    v_low ^= D0 g_low + D1 (w g)_low (L - 1 point products), D1 erased by measurement, Y: slot(L) ^= Z0 X.
    Carries at the bottom (bottom-up): Z = [C] q0_U' + [not C and not F], ZB = [not B] Z;
    d_U' ^= ZB X (one product), d_L' ^= [not B and (C or F)] X (one controlled XOR).
    Return = exact transpose after Q (the kit's q_posts at the bottom): carries^T top-down at the bottom; at L:
    X ^= Z0 slot(L), S0 -> slot(L), S1 -> temp, slot(L) *= 1 / a_L, temp ^= slot(L) g[-2], X ^= temp c (temp erased
    by measurement); at U: X *= 1 / a_U (the kit's own last scaling), no dots."""

    def pre(self, f, g, P, parent=None):
        m = len(f) - 1
        if m not in self.upper and m not in self.lower:
            return Kit.pre(self, f, g, P, parent)
        M = self.M
        R = M.R
        F0, a, ai = M.clip(g[-1])
        q0, q1 = f[-2], f[-1]
        if m in self.upper:
            for v in P:
                M.scl(v[-1], ai, 'pay')                       # X = t / a_U
        M.scl(q1, ai, 'loc')
        M.mx(q1, g[-2], q0, 'loc')
        M.scl(q0, ai, 'loc')
        M.conv(g[:-1], [q0, q1], f[:-2], 'loc')
        if m in self.lower:
            B = R[q1] == 0
            M.cnt[('route', 'flag')] += 1
            c = np.where(B, 1, R[q1])
            Z0 = np.where(B, 0, R[q0])
            M.cnt[('route', 'and12')] += 1
            for v, t in zip(P, parent[0]):
                D1 = mul(R[t], c)
                R[v[-1]] = R[v[-1]] ^ mul(D1, R[g[-2]])
                M.cnt[('pay', 'prod')] += 2
                M.scl(v[-1], ai, 'pay')                       # D0
                for i in range(m - 1):
                    acc = mul(R[v[-1]], R[g[i]])
                    if i >= 1:
                        acc = acc ^ mul(D1, R[g[i - 1]])
                    R[v[i]] = R[v[i]] ^ acc
                M.cnt[('pay', 'prod')] += m - 1
                R[v[-1]] = R[v[-1]] ^ mul(Z0, R[t])           # Y
                M.cnt[('pay', 'prod')] += 1
        M.unclip(g[-1], F0)

    def _route_words(self, U):
        M = self.M
        R = M.R
        L, Up, Lp = U - 1, U - 2, U - 3
        Bf, Cf, FU = self.qflag(L), self.qflag(L - 1), self.qflag(Lp)
        q0u = self.lv[Up]['f'][-2]
        C2 = self.qflag(U - 4)                              # gap at L' too: the second part moves to U - 4
        Z = np.where(Cf, R[q0u], np.where(~FU, 1, 0))
        ZB = np.where(Bf, 0, Z)
        # zero region (padded charts): there every level reads as a gap; F(U') and F(L') both set means the
        # divisor at U' is 0 (single gaps never come two in a row), and the carry vanishes
        selL = ~Bf & ~C2 & (Cf | FU)
        selU2 = ~Bf & Cf & C2 & ~FU
        if os.environ.get('ZN_PD_SU2XOR') == '1':           # sixth holder: B5's (F, C2) law, transposed (t2 ^= X)
            selU2 = ~Bf & C2 & (Cf ^ FU)
        M.cnt[('route', 'and12')] += 2
        return ZB, selL, selU2

    def at_bottom(self, f, g, P):
        M = self.M
        R = M.R
        for U in sorted(self.upper):                          # carries, bottom-up
            ZB, selL, selU2 = self._route_words(U)
            tu, tl, t2 = self.tops(U - 2, P), self.tops(U - 3, P), self.tops(U - 4, P)
            for k, vu in enumerate(self.lv[U]['P']):
                X = R[vu[-1]]
                R[tu[k]] = R[tu[k]] ^ mul(ZB, X)
                R[tl[k]] = R[tl[k]] ^ np.where(selL, X, 0)
                R[t2[k]] = R[t2[k]] ^ np.where(selU2, X, 0)
                M.cnt[('route', 'prod')] += 1
                M.cnt[('route', 'cxor')] += 2

    def after_q(self, f, g, P):
        M = self.M
        R = M.R
        for U in sorted(self.upper, reverse=True):            # carries transposed, top-down
            ZB, selL, selU2 = self._route_words(U)
            tu, tl, t2 = self.tops(U - 2, P), self.tops(U - 3, P), self.tops(U - 4, P)
            for k, vu in enumerate(self.lv[U]['P']):
                R[vu[-1]] = (R[vu[-1]] ^ mul(ZB, R[tu[k]]) ^ np.where(selL, R[tl[k]], 0)
                             ^ np.where(selU2, R[t2[k]], 0))
                M.cnt[('route', 'prod')] += 1
                M.cnt[('route', 'cxor')] += 2

    def post_look(self, f, g, P, parent):
        m = len(f) - 1
        if m not in self.upper and m not in self.lower:
            return KitQB.post_look(self, f, g, P, parent)
        M = self.M
        R = M.R
        F0, a, ai = M.clip(g[-1])
        q0, q1 = f[-2], f[-1]
        if m in self.upper:
            for v in P:
                M.scl(v[-1], ai, 'pay')                       # X^T
        else:
            B = R[q1] == 0
            M.cnt[('route', 'flag')] += 1
            c = np.where(B, 1, R[q1])
            Z0 = np.where(B, 0, R[q0])
            M.cnt[('route', 'and12')] += 1
            for v, t in zip(P, parent[0]):
                R[t] = R[t] ^ mul(Z0, R[v[-1]])               # Y^T
                s0 = np.zeros(M.n, dtype=I64)
                s1 = np.zeros(M.n, dtype=I64)
                for i in range(m - 1):
                    s0 = s0 ^ mul(R[g[i]], R[v[i]])
                    if i >= 1:
                        s1 = s1 ^ mul(R[g[i - 1]], R[v[i]])
                R[v[-1]] = R[v[-1]] ^ s0
                M.cnt[('pay', 'prod')] += 1 + (m - 1)
                M.scl(v[-1], ai, 'pay')
                s1 = s1 ^ mul(R[v[-1]], R[g[-2]])
                R[t] = R[t] ^ mul(s1, c)
                M.cnt[('pay', 'prod')] += 2
        M.conv(g[:-1], [q0, q1], f[:-2], 'loc')
        M.scl(q0, a, 'loc')
        M.mx(q1, g[-2], q0, 'loc')
        M.scl(q1, a, 'loc')
        M.unclip(g[-1], F0)


PRICE = {('pay', 'prod'): 42, ('pay', 'scale'): 42, ('pay', 'mprod'): 66, ('route', 'prod'): 42,
         ('route', 'mprod'): 66, ('route', 'cswap'): 12, ('route', 'cxor'): 12, ('route', 'flag'): 22,
         ('route', 'and12'): 12, ('loc', 'inverse'): 98, ('loc', 'prod'): 42, ('loc', 'scale'): 42,
         ('loc', 'conv_pts'): 42}


def price(cnt):
    """Toffolis per direction from the counts, at the kernel prices above [E: 42 per product / scaling as R7-1;
    66 masked product; 12 per controlled word op; 22 per zero test with its clear; 98 per 12-bit inverse]."""
    return sum(PRICE.get(k, 0) * v for k, v in cnt.items())


def compare(O, E, u, w, cls=None, **kw):
    """Run K and a variant on the same inputs; per-shot equality of every output word."""
    ko, kc, _ = run_bez(Kit, O, E, u, w)
    po, pc, _ = run_bez(cls or Paired, O, E, u, w, **kw)
    eq = np.ones(O.shape[0], dtype=bool)
    for k in ko:
        eq &= (ko[k] == po[k]).all(axis=1)
    return eq, kc, pc


# ------------------------------------------------------------------------------------------- charts and I/O
def frob_sqrt(x):
    """x -> x^2048 (the square root in GF(2^12)), elementwise."""
    y = np.array(x, dtype=I64)
    for _ in range(11):
        y = mul(y, y)
    return y


def bez_inputs(A, B):
    """(A, B) as the decode feeds convert (35 + 35 words, B35 = 1 implicit) -> (O, E, u, w) at bez entry
    (frob 2048 on both, then w[i] ^= top * O[i], as in convert_safe forward)."""
    Bs = frob_sqrt(B)
    As = frob_sqrt(A)
    E = Bs[:, 0::2]            # 18 words B0, B2, .., B34
    O = Bs[:, 1::2]            # 17 words B1, .., B33
    u = As[:, 1::2]            # 17 words
    w = As[:, 0:34:2].copy()   # 17 words
    top = As[:, 34]
    for i in range(17):
        w[:, i] ^= mul(top, O[:, i])
    return O, E, u, w


def run_bez(cls, O, E, u, w, **kw):
    n = O.shape[0]
    M = Mach(n)
    hO = [M.new(O[:, i]) for i in range(17)]
    hE = [M.new(E[:, i]) for i in range(18)]
    hu = [M.new(u[:, i]) for i in range(17)]
    hw = [M.new(w[:, i]) for i in range(17)]
    K = cls(M, **kw)
    K.bez(hO, hE, [hu, hw])
    out = dict(O=np.stack([M.R[h] for h in hO], 1), E=np.stack([M.R[h] for h in hE], 1),
               u=np.stack([M.R[h] for h in hu], 1), w=np.stack([M.R[h] for h in hw], 1))
    return out, M.cnt, K


REAL = os.path.expanduser('~/qcc/ws/dqi_zn/box_out/charts_real.npz')


def chart_sets(kinds, rng):
    """Yield (label, A, B, weights or None). Kind 'r' = the real converter charts (verifier-sampled, weights 18..35)."""
    import dqi_int_convtest as CT  # charts (B5, read-only)
    for kd in kinds:
        if kd == 'r':
            z = np.load(REAL)
            yield 'real charts', z['An'], z['Bn'], z['w']
        else:
            A, B, label = CT.charts(int(kd), rng)
            yield label, A, B, None


def main():
    kinds = (sys.argv[1] if len(sys.argv) > 1 else '0,1,2,3,4,r').split(',')
    rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    route = sys.argv[3] if len(sys.argv) > 3 else "generic"
    upper = tuple(int(x) for x in (sys.argv[4] if len(sys.argv) > 4 else "17,15,13,11,9,7,5,3").split(","))
    for label, A, B, wt in chart_sets(kinds, rng):
        O, E, u, w = bez_inputs(A, B)
        cname = sys.argv[5] if len(sys.argv) > 5 else 'Paired'
        cls = globals()[cname]
        kw = {} if cname == 'KitQB' else dict(route=route, upper=upper)
        eq, kc, pc = compare(O, E, u, w, cls=cls, **kw)
        line = '%-22s P == K on %d / %d' % (label, eq.sum(), len(eq))
        if wt is not None:
            bad = collections.Counter(int(x) for x in wt[~eq])
            line += '; failures by weight %s' % dict(sorted(bad.items()))
        print(line, flush=True)
    print('K counts', dict(kc))
    print('P counts', dict(pc))
    pk, pp = price(kc), price(pc)
    print('price per direction [E kernel prices]: kit %d, variant %d, delta %+d' % (pk, pp, pp - pk))


if __name__ == '__main__':
    main()
