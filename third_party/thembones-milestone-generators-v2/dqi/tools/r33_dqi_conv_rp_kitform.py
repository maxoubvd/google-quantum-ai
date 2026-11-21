#!/usr/bin/env python3
"""r33_dqi_conv_rp_kitform.py -- R33 re-pricing of DQI-CONV: the 2-word payload digits in the KIT's OWN form.

Reference: a classical, word-level port of the kit's half-Euclid sweep with its payloads, vectorised over shots
(playbook/kit/decoder/exp/safe_scalar/converter.py: look_solve, pre, post_look, q_post, base, rotate,
canonical_permute; convert_short_hankel/kernel.py: Frobenius, split, w += top O, the monic root word). The kit is
non-monic: the payload's top word is divided by the pivot (safe divisor: 0 -> 1), every level zero-tests its pivot,
rotates g and permutes (f, g) on a gap, and q_post takes the gap path (swaps with the parent, masked product) when
q1 = 0. Every op is executed on every shot (circuit semantics), so the counts are static.

Paired variant (the design under re-pricing), for a set of step levels m (block = kit levels m+1, m):
  descent at level m (after its locator step): free pre-correction of the low words from the two top words,
    D1 = v_m / beta, D0 = (v_(m-1) + D1 g_(m-2)) / beta, rest by m-1 rank updates at public points,
    c1 = D1 / q1, c2 = D0 + q0 c1; level m+1 does no payload op.
  correction: the next digit word (payload index m-2) gets c1_raw when the pivot of level m-1 is nonzero
    (controlled XOR), applied bottom-up after the first unpaired level below the lowest block.
  return: the blocks' Gram scalings (1/q1 of each level) are applied early, the transposed corrections top-down,
    then each block runs the exact transposes (m+1 products, 3 scalings); level m+1 does no payload op.
Compares every word (payloads and locator) after the sweep with the reference, per shot, and classifies the
mismatches by the reference's own flags.
Usage: python r33_dqi_conv_rp_kitform.py real|syn [n_syn=20000] [seed=1] [steps=16,14,12,10 ...]
"""
import os
import sys
import time
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dqi_arch_model as M  # noqa: E402

EXP, LOG, ORD = M.EXP, M.LOG, M.ORD
I64 = np.int64


def fmul(a, b):
    a = np.asarray(a, dtype=I64)
    b = np.asarray(b, dtype=I64)
    return np.where((a == 0) | (b == 0), 0, EXP[LOG[a] + LOG[b]])


def fsafe(a):
    return np.where(a == 0, 1, a)


def finv_safe(a):
    return EXP[(ORD - LOG[fsafe(a)]) % ORD]


def fsqrt(a):
    a = np.asarray(a, dtype=I64)
    return np.where(a == 0, 0, EXP[(LOG[a] * 2048) % ORD])


def cmul(c, a):
    """public constant times a word vector (free)."""
    return fmul(np.full(np.shape(a), c, dtype=I64), a)


# ---- public point tables for the m-1 point rank updates ----
PTS = [int(EXP[i]) for i in range(1, 40)]          # distinct nonzero public points


def ppow(p, k):
    return int(EXP[(LOG[p] * k) % ORD]) if k else 1


_TAB = {}


def point_tables(m):
    """For a step level m (g of degree m-1): points p_1..p_(m-1), the interpolation matrix INT (coefficients of
    the degree <= m-2 interpolant from the m-1 values) and the pre-correction vectors z^m mod Pi, z^(m-1) mod Pi."""
    if m in _TAB:
        return _TAB[m]
    n = m - 1
    p = PTS[:n]
    V = np.array([[ppow(x, k) for k in range(n)] for x in p], dtype=I64)   # values = V @ coeffs
    INT = gf_inv_matrix(V)                                                   # coeffs = INT @ values
    Pi = np.array([1], dtype=I64)
    for x in p:
        Pi = M.pmul(Pi, np.array([x, 1], dtype=I64))
    zm = np.zeros(m + 1, dtype=I64)
    zm[m] = 1
    zm1 = np.zeros(m, dtype=I64)
    zm1[m - 1] = 1
    cm = np.zeros(n, dtype=I64)
    cm1 = np.zeros(n, dtype=I64)
    r = M.pmod(zm, Pi)
    cm[:len(r)] = r
    r = M.pmod(zm1, Pi)
    cm1[:len(r)] = r
    _TAB[m] = (p, INT, cm, cm1)
    return _TAB[m]


def gf_inv_matrix(A):
    n = A.shape[0]
    A = A.copy()
    Iv = np.eye(n, dtype=I64)
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r, c])
        if piv != c:
            A[[c, piv]] = A[[piv, c]]
            Iv[[c, piv]] = Iv[[piv, c]]
        s = int(M.inv(A[c, c]))
        A[c] = M.mul(s, A[c])
        Iv[c] = M.mul(s, Iv[c])
        for r in range(n):
            if r != c and A[r, c]:
                f = int(A[r, c])
                A[r] ^= M.mul(f, A[c])
                Iv[r] ^= M.mul(f, Iv[c])
    return Iv


class Sim:
    def __init__(self, W):
        self.W = W
        self.cnt = Counter()

    def mul(self, x, y, t, tag):
        self.W[t] ^= fmul(self.W[x], self.W[y])
        self.cnt[tag + ' products'] += 1

    def scale(self, x, a, divide, tag):
        self.W[x] = fmul(self.W[x], finv_safe(self.W[a]) if divide else fsafe(self.W[a]))
        self.cnt[tag + ' scalings'] += 1

    def cswap(self, F, i, j, tag='route'):
        wi = self.W[i].copy()
        wj = self.W[j].copy()
        self.W[i] = np.where(F, wj, wi)
        self.W[j] = np.where(F, wi, wj)
        self.cnt[tag + ' cswap words'] += 1


class Conv:
    """One run of the converter's sweep (bez_lookahead with the monic root word) on all shots."""

    def __init__(self, sim, steps=()):
        self.s = sim
        self.steps = set(steps)                       # paired step levels m
        self.firsts = {m + 1 for m in self.steps}     # first level of each block: no payload op
        self.low = min(self.steps) - 1 if self.steps else None
        self.lv = {}                                  # level -> dict(f, g, F) recorded on the descent
        self.flags = {}                               # level -> zero-pivot flag (vector)
        self.term = {}                                # level -> g entirely zero (vector)

    # ---- kit pieces ----
    def rotate(self, g, F, reverse=False):
        rng = range(1, len(g)) if reverse else range(len(g) - 1, 0, -1)
        for j in rng:
            self.s.cswap(F, g[j], g[j - 1])

    def permute(self, f, g, F, reverse=False):
        cyc = []
        for i in range(len(g) - 1):
            cyc += [g[i], f[i]]
        cyc += [g[-1]]
        seq = range(len(cyc) - 2, -1, -1) if reverse else range(len(cyc) - 1)
        for i in seq:
            self.s.cswap(F, cyc[i], cyc[i + 1])

    def conv(self, f, g, q0, q1):
        W = self.s.W
        m = len(f) - 1
        for j in range(m - 1):
            t = fmul(W[q0], W[g[j]])
            if j >= 1:
                t ^= fmul(W[q1], W[g[j - 1]])
            W[f[j]] ^= t
        self.s.cnt['locator products'] += m

    def pre(self, f, g, P, parent):
        s, W = self.s, self.s.W
        m = len(f) - 1
        n = m - 1
        a, q0, q1 = g[-1], f[-2], f[-1]
        if m not in self.firsts and m not in self.steps:
            for v in P:
                s.scale(v[-1], a, True, 'payload')
            for v in P:
                for i in range(n):
                    s.mul(g[i], v[-1], v[i], 'payload')
        s.scale(q1, a, True, 'locator')
        s.mul(q1, g[-2], q0, 'locator')
        s.scale(q0, a, True, 'locator')
        self.conv(f, g, q0, q1)
        if m in self.steps:
            self.paired_down(f, g, P, parent)

    def paired_down(self, f, g, P, parent):
        s, W = self.s, self.s.W
        m = len(f) - 1
        n = m - 1
        a, q0, q1 = g[-1], f[-2], f[-1]
        p, INT, cm, cm1 = point_tables(m)
        tops = parent[0]
        gs = W[g].copy()
        gs[-1] = fsafe(gs[-1])                       # the word a reads 1 inside its safe inverse context
        gval = [np.bitwise_xor.reduce(np.stack([cmul(ppow(x, k), gs[k]) for k in range(m)]), axis=0) for x in p]
        for v, t in zip(P, tops):
            c1, c2, low = t, v[-1], v[:-1]
            for i in range(n):                       # (1) free pre-correction of the low words
                W[low[i]] ^= cmul(int(cm[i]), W[c1]) ^ cmul(int(cm1[i]), W[c2])
            s.scale(c1, a, True, 'payload')          # (2) D1
            s.mul(c1, g[-2], c2, 'payload')          # (3) D0
            s.scale(c2, a, True, 'payload')
            vals = []
            for i, x in enumerate(p):                # (4) m-1 rank updates at public points
                vals.append(fmul(cmul(x, W[c1]) ^ W[c2], gval[i]))
            s.cnt['payload products'] += n
            for j in range(n):
                acc = np.zeros_like(W[c1])
                for i in range(n):
                    if INT[j, i]:
                        acc ^= cmul(int(INT[j, i]), vals[i])
                W[low[j]] ^= acc
            s.scale(c1, q1, True, 'payload')         # (5) c1 = D1 / q1
            s.mul(q0, c1, c2, 'payload')             #     c2 = D0 + q0 c1

    def paired_up(self, f, g, P, parent):
        s, W = self.s, self.s.W
        m = len(f) - 1
        n = m - 1
        a, q0, q1 = g[-1], f[-2], f[-1]
        p, INT, cm, cm1 = point_tables(m)
        tops = parent[0]
        gs = W[g].copy()
        gs[-1] = fsafe(gs[-1])
        gval = [np.bitwise_xor.reduce(np.stack([cmul(ppow(x, k), gs[k]) for k in range(m)]), axis=0) for x in p]
        for v, t in zip(P, tops):
            c1, c2, low = t, v[-1], v[:-1]
            s.mul(q0, c2, c1, 'payload')             # (5^T) y_D1 = (y_c1 + q0 y_c2) / q1, y_D0 = y_c2
            s.scale(c1, q1, True, 'payload')
            yv = []                                  # (4^T) INT^T y_low, one product per point
            for i in range(n):
                acc = np.zeros_like(W[c1])
                for j in range(n):
                    if INT[j, i]:
                        acc ^= cmul(int(INT[j, i]), W[low[j]])
                yv.append(fmul(acc, gval[i]))
            s.cnt['payload products'] += n
            for i, x in enumerate(p):
                W[c1] ^= cmul(x, yv[i])
                W[c2] ^= yv[i]
            s.scale(c2, a, True, 'payload')          # (3^T)
            s.mul(g[-2], c2, c1, 'payload')
            s.scale(c1, a, True, 'payload')          # (2^T)
            for i in range(n):                       # (1^T) free
                W[c1] ^= cmul(int(cm[i]), W[low[i]])
                W[c2] ^= cmul(int(cm1[i]), W[low[i]])

    def q_post(self, q, P, parent, m):
        s, W = self.s, self.s.W
        skip = (m in self.firsts) or (m in self.steps)
        if parent is None:
            if not skip:
                for v in P:
                    s.scale(v[-1], q, True, 'payload')
            return
        tops, pw, c = parent
        F = W[q] == 0
        s.cswap(F, q, pw, 'route')
        if skip:
            s.cswap(F, q, pw, 'route')
            return
        for v, t in zip(P, tops):
            s.cswap(F, v[-1], t, 'route')
        for v in P:
            s.scale(v[-1], q, True, 'payload')
        for v, t in zip(P, tops):
            W[t] ^= np.where(F, fmul(W[v[-1]], W[c]), 0)
            s.cnt['payload masked products'] += 1
        s.cswap(F, q, pw, 'route')

    def post_look(self, f, g, P, parent):
        s, W = self.s, self.s.W
        m = len(f) - 1
        n = m - 1
        a, q0, q1 = g[-1], f[-2], f[-1]
        self.q_post(q1, P, parent, m)
        if m == self.low:
            self.early_gram_and_transposed_corrections(P)
        if m in self.steps:
            self.paired_up(f, g, P, parent)
        elif m not in self.firsts:
            for v in P:
                for i in range(n):
                    s.mul(g[i], v[i], v[-1], 'payload')
            for v in P:
                s.scale(v[-1], a, True, 'payload')
        self.conv(f, g, q0, q1)
        s.scale(q0, a, False, 'locator')
        s.mul(q1, g[-2], q0, 'locator')
        s.scale(q1, a, False, 'locator')

    def base(self, f, g, P, parent):
        s = self.s
        a, q = g[0], f[-1]
        for v in P:
            s.scale(v[0], a, True, 'payload')
        s.scale(q, a, True, 'locator')
        self.q_post(q, P, parent, 1)
        for v in P:
            s.scale(v[0], a, True, 'payload')
        s.scale(q, a, False, 'locator')

    # ---- the correction network of the paired variant ----
    def corrections_down(self, P):
        W = self.s.W
        P = self.Pfull
        for m in sorted(self.steps):                 # bottom-up: raw c1 words are read before they change
            sflag = ~self.flags[m - 1]
            for v in P:
                W[v[m - 2]] ^= np.where(sflag, W[v[m]], 0)
                self.s.cnt['payload controlled XOR words'] += 1

    def early_gram_and_transposed_corrections(self, P):
        s, W = self.s, self.s.W
        P = self.Pfull
        for m in sorted(self.steps):                # the blocks' Gram scalings, moved here
            for lvl, idx in ((m + 1, m), (m, m - 1)):
                q1 = self.lv[lvl]['f'][-1]
                for v in P:
                    s.scale(v[idx], q1, True, 'payload')
        for m in sorted(self.steps, reverse=True):   # transposed corrections, top-down
            sflag = ~self.flags[m - 1]
            for v in P:
                W[v[m]] ^= np.where(sflag, W[v[m - 2]], 0)
                self.s.cnt['payload controlled XOR words'] += 1

    def look_solve(self, f, g, P, parent=None):
        W = self.s.W
        m = len(f) - 1
        if m == 1:
            return self.base(f, g, P, parent)
        F = W[g[-1]] == 0
        self.flags[m] = F
        self.term[m] = np.all(W[g] == 0, axis=0)
        self.lv[m] = dict(f=f, g=g)
        self.s.cnt['zero tests'] += 1
        self.rotate(g, F)
        self.pre(f, g, P, parent)
        self.permute(f, g, F)
        bad = (W[g[-1]] == 0) != F
        self.s.cnt['clear_flag violations'] += int(bad.sum())
        child =([v[-1] for v in P], f[-1], f[-2])
        self.look_solve(g, f[:-2], [v[:-1] for v in P], child)
        F2 = W[g[-1]] == 0
        self.permute(f, g, F2, True)
        self.post_look(f, g, P, parent)
        self.rotate(g, F2, True)

    def run(self, O, E, payloads, lead):
        W = self.s.W
        for i in range(len(O)):
            W[E[i]] ^= fmul(W[E[-1]], W[O[i]])
        W[lead] = 1
        self.Pfull = payloads
        self.look_solve(O + [lead], E[:-1], payloads)
        for i in range(len(O)):
            W[E[i]] ^= fmul(W[E[-1]], W[O[i]])


class ConvPaired(Conv):
    def pre(self, f, g, P, parent):
        super().pre(f, g, P, parent)
        if len(f) - 1 == self.low:
            self.corrections_down(P)


class ConvRouted(Conv):
    """Paired variant with descent routes for interior gaps (Fa: gap at the block's first level m+1, Fb: at the step
    level m, Fc: at level m-1 right below):
      c2 = D0 + [~Fa] q0 c1                                (masked product)
      generic:  digit(m-1) ^= c1 if ~Fc, digit(m-2) ^= c1 if Fc
      Fa:       no correction (the kit's two levels already form the 2-word digit)
      Fb:       digit(m-1) ^= c1 q0_(m-1) (masked product), digit(m-2) ^= c1
    All corrections are applied bottom-up after the pre of level min(steps) - 2. Only the descent is routed here:
    the comparison is made on the payload state at the bottom of the descent."""

    def __init__(self, sim, steps=()):
        super().__init__(sim, steps)
        self.low2 = min(self.steps) - 2

    def pre(self, f, g, P, parent):
        super().pre(f, g, P, parent)
        if len(f) - 1 == self.low2:
            self.routed_corrections()

    def paired_down(self, f, g, P, parent):
        m = len(f) - 1
        self.masked_q0 = True
        super().paired_down(f, g, P, parent)

    def routed_corrections(self):
        s, W, P = self.s, self.s.W, self.Pfull
        fl = self.flags
        for m in sorted(self.steps):
            # interior gap at level L = zero pivot at L and a nonzero pivot at L-1 (an early end zeroes both)
            Fa = fl[m + 1] & ~fl[m]
            Fb = fl[m] & ~fl[m - 1]
            Fc = fl[m - 1] & ~fl[m - 2]
            q0n = self.lv[m - 1]['f'][-2]
            gen = ~Fa & ~Fb
            for v in P:
                c1 = W[v[m]]
                W[v[m - 2]] ^= np.where(gen & ~fl[m - 1], c1, 0)
                W[v[m - 3]] ^= np.where((gen & Fc) | Fb, c1, 0)
                W[v[m - 2]] ^= np.where(Fb, fmul(c1, W[q0n]), 0)
                s.cnt['payload controlled XOR words'] += 2
                s.cnt['payload masked products'] += 1


def conv_with_masked_q0(cls):
    """Wrap paired_down so that the conversion's q0 c1 term is masked by ~Fa (gap at the block's first level)."""
    orig = cls.paired_down

    def paired_down(self, f, g, P, parent):
        m = len(f) - 1
        W = self.s.W
        q0 = f[-2]
        Fa = self.flags.get(m + 1)
        if Fa is None:
            return orig(self, f, g, P, parent)
        saved = W[q0].copy()
        # the masked product is counted as a product here and as +1 masked product per payload below
        W[q0] = np.where(Fa, 0, saved)
        try:
            orig(self, f, g, P, parent)
        finally:
            W[q0] = saved
        self.s.cnt['payload masked products (conversion)'] += 2
    cls.paired_down = paired_down
    return cls


conv_with_masked_q0(ConvRouted)


class Snap:
    """Mixin: keep a copy of all words at the bottom of the descent (before the base's return part)."""

    def base(self, f, g, P, parent):
        s = self.s
        a, q = g[0], f[-1]
        for v in P:
            s.scale(v[0], a, True, 'payload')
        s.scale(q, a, True, 'locator')
        self.bottom = s.W.copy()
        self.q_post(q, P, parent, 1)
        for v in P:
            s.scale(v[0], a, True, 'payload')
        s.scale(q, a, False, 'locator')


class ConvFull(ConvRouted):
    """Routed descent + full return. The kit's return equals [all q_post Grams bottom-up] then [all transposed level
    ops bottom-up] (each q_post touches only its own top and its parent's top; each transposed level op writes its
    own top from lower words). So the return here keeps the kit's q_post at every level unchanged (applied at the
    bottom) and replaces the transposed level ops by the exact transposes of the paired and routed descent."""

    def base(self, f, g, P, parent):
        s = self.s
        a, q = g[0], f[-1]
        for v in P:
            s.scale(v[0], a, True, 'payload')
        s.scale(q, a, True, 'locator')
        self.lv[1] = dict(f=f, g=g)
        self.bottom = s.W.copy()
        for L in range(1, 18):                       # the kit's Gram stage, bottom-up, unchanged
            fL = self.lv[L]['f']
            PL = [v[:L] for v in self.Pfull]
            par = None if L == 17 else ([v[L] for v in self.Pfull], self.lv[L + 1]['f'][-1],
                                        self.lv[L + 1]['f'][-2])
            self.q_post_kit(fL[-1], PL, par)
        for v in P:
            s.scale(v[0], a, True, 'payload')
        s.scale(q, a, False, 'locator')

    def q_post_kit(self, q, P, parent):
        s, W = self.s, self.s.W
        if parent is None:
            for v in P:
                s.scale(v[-1], q, True, 'payload')
            return
        tops, pw, c = parent
        F = W[q] == 0
        s.cswap(F, q, pw, 'route')
        for v, t in zip(P, tops):
            s.cswap(F, v[-1], t, 'route')
        for v in P:
            s.scale(v[-1], q, True, 'payload')
        for v, t in zip(P, tops):
            W[t] ^= np.where(F, fmul(W[v[-1]], W[c]), 0)
            s.cnt['payload masked products'] += 1
        s.cswap(F, q, pw, 'route')

    def transposed_corrections(self):
        s, W, P, fl = self.s, self.s.W, self.Pfull, self.flags
        for m in sorted(self.steps, reverse=True):   # top-down: the reverse of the descent's bottom-up order
            Fa = fl[m + 1] & ~fl[m]
            Fb = fl[m] & ~fl[m - 1]
            Fc = fl[m - 1] & ~fl[m - 2]
            q0n = self.lv[m - 1]['f'][-2]
            gen = ~Fa & ~Fb
            for v in P:
                W[v[m]] ^= np.where(Fb, fmul(W[v[m - 2]], W[q0n]), 0)
                W[v[m]] ^= np.where((gen & Fc) | Fb, W[v[m - 3]], 0)
                W[v[m]] ^= np.where(gen & ~fl[m - 1], W[v[m - 2]], 0)
                s.cnt['payload controlled XOR words'] += 2
                s.cnt['payload masked products'] += 1

    def post_look(self, f, g, P, parent):
        s, W = self.s, self.s.W
        m = len(f) - 1
        n = m - 1
        a, q0, q1 = g[-1], f[-2], f[-1]
        if m == self.low2:
            self.transposed_corrections()
        if m in self.steps:
            Fa = self.flags.get(m + 1)
            saved = W[q0].copy()
            W[q0] = np.where(Fa, 0, saved)           # the masked conversion, transposed
            self.paired_up(f, g, P, parent)
            W[q0] = saved
            s.cnt['payload masked products (conversion)'] += 2
        elif m not in self.firsts:
            for v in P:
                for i in range(n):
                    s.mul(g[i], v[i], v[-1], 'payload')
            for v in P:
                s.scale(v[-1], a, True, 'payload')
        self.conv(f, g, q0, q1)
        s.scale(q0, a, False, 'locator')
        s.mul(q1, g[-2], q0, 'locator')
        s.scale(q1, a, False, 'locator')


class RefSnap(Snap, Conv):
    pass


class PairedSnap(Snap, ConvPaired):
    pass


class RoutedSnap(Snap, ConvRouted):
    pass


def layout():
    A = list(range(0, 35))
    B = list(range(35, 70))
    lead = 70
    return A, B, lead, 71


def load_real(path):
    d = np.load(path)
    n = int(d['shots'])
    bits = np.unpackbits(d['ab'], axis=1)[:, :n].astype(I64)
    words = np.zeros((70, n), dtype=I64)
    for w in range(70):
        for j in range(12):
            words[w] |= bits[12 * w + j] << j
    return words, d['weights'].astype(int)


def make_syn(n, rng, low=False):
    law = M.weight_law()
    law = law / law.sum()
    wts = rng.integers(14, 26, size=n) if low else rng.choice(len(law), size=n, p=law)
    words = np.zeros((70, n), dtype=I64)
    words[:35] = rng.integers(0, 4096, size=(35, n))
    for s in range(n):
        w = int(wts[s])
        t = 35 - w
        if w == 0:
            continue
        cf = rng.integers(0, 4096, size=w)
        cf[0] = rng.integers(1, 4096)
        words[35 + t:35 + t + w, s] = cf          # B = z^t B~, B~ monic of degree w (B35 = 1 implicit)
    return words, wts


def sweep(words, steps, cls):
    A, B, lead, nw = layout()
    S = words.shape[1]
    W = np.zeros((nw, S), dtype=I64)
    W[:70] = fsqrt(words)                          # Frobenius on both, as convert_safe does first
    sim = Sim(W)
    E = B[::2]
    O = B[1::2]
    u = A[1::2]
    w = A[:34:2]
    top = A[-1]
    for i in range(17):
        W[w[i]] ^= fmul(W[top], W[O[i]])
    c = cls(sim, steps)
    c.run(O, E, [u, w], lead)
    return W, sim.cnt, c


def descent_check(words, wts, steps):
    """Compare the state at the bottom of the descent: reference, plain paired, routed."""
    S = words.shape[1]
    _, cr, ref = sweep(words, (), RefSnap)
    _, cp, pp = sweep(words, steps, PairedSnap)
    _, cq, rr = sweep(words, steps, RoutedSnap)
    interior = {m: ref.flags[m] & ~ref.term[m] for m in ref.flags}
    for name, run in (('plain paired', pp), ('routed', rr)):
        mism = np.any(run.bottom != ref.bottom, axis=0)
        pats = Counter()
        for sidx in np.flatnonzero(mism):
            pats[tuple(m for m in range(17, 1, -1) if interior[m][sidx])] += 1
        gap_any = np.zeros(S, dtype=bool)
        for m in range(min(steps) - 3, 18):
            gap_any |= interior.get(m, np.zeros(S, dtype=bool))
        print('  descent %-12s steps %s: differing shots %d / %d (interior-gap shots near the blocks %d); '
              'gap levels of the differing shots: %s'
              % (name, steps, mism.sum(), S, gap_any.sum(), dict(pats.most_common(12))))
        if name == 'routed':
            ok = gap_any & ~mism
            print('    interior-gap shots repaired by the descent routes: %d of %d' % (ok.sum(), gap_any.sum()))
    return cq


def full_check(words, wts, steps):
    """Whole sweep (descent and return): kit reference against the routed full variant."""
    S = words.shape[1]
    Wr, cr, ref = sweep(words, (), Conv)
    Wf, cf, ff = sweep(words, steps, ConvFull)
    interior = {m: ref.flags[m] & ~ref.term[m] for m in ref.flags}
    mism = np.any(Wf != Wr, axis=0)
    gap_any = np.zeros(S, dtype=bool)
    for m in range(min(steps) - 3, 18):
        gap_any |= interior.get(m, np.zeros(S, dtype=bool))
    pats = Counter()
    for sidx in np.flatnonzero(mism):
        pats[tuple(m for m in range(17, 1, -1) if interior[m][sidx])] += 1
    print('  full sweep, steps %s: differing shots %d / %d; interior-gap shots near the blocks %d, of them '
          'differing %d; gap levels of the differing shots: %s; weights of the differing shots: %s'
          % (steps, mism.sum(), S, gap_any.sum(), (mism & gap_any).sum(), dict(pats.most_common(10)),
             dict(sorted(Counter(wts[mism].tolist()).items()))))
    keys = sorted(set(cr) | set(cf))
    print('    counts per convert (reference -> routed):',
          ', '.join('%s %d -> %d' % (k, cr[k], cf[k]) for k in keys if cr[k] != cf[k] or 'payload' in k))
    return mism


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
    if mode.startswith('full'):
        sub = mode.split(':')[1] if ':' in mode else 'real'
        nsyn = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
        seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
        configs = [tuple(int(x) for x in a.split(',')) for a in sys.argv[4:]] or [(16, 14, 12, 10, 8)]
        if sub == 'real':
            words, wts = load_real(os.path.expanduser('~/qcc/ws/dqi_build/r7/ab/ab_data.npz'))
        else:
            words, wts = make_syn(nsyn, np.random.default_rng(seed), low=(sub == 'low'))
        print('full check, %s, shots %d' % (sub, words.shape[1]))
        for steps in configs:
            full_check(words, wts, steps)
        return
    if mode.startswith('descent'):
        sub = mode.split(':')[1] if ':' in mode else 'real'
        nsyn = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
        seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
        configs = [tuple(int(x) for x in a.split(',')) for a in sys.argv[4:]] or [(16, 14, 12, 10, 8)]
        if sub == 'real':
            words, wts = load_real(os.path.expanduser('~/qcc/ws/dqi_build/r7/ab/ab_data.npz'))
        else:
            words, wts = make_syn(nsyn, np.random.default_rng(seed), low=(sub == 'low'))
        print('descent check, %s, shots %d' % (sub, words.shape[1]))
        for steps in configs:
            cq = descent_check(words, wts, steps)
            print('  routed counts (whole sweep, both payloads):', dict(sorted(cq.items())))
        return
    nsyn = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    configs = [tuple(int(x) for x in a.split(',')) for a in sys.argv[4:]] or [(16, 14, 12, 10), (16, 14, 12, 10, 8),
                                                                              (16, 14, 12, 10, 8, 6)]
    t0 = time.time()
    if mode == 'real':
        words, wts = load_real(os.path.expanduser('~/qcc/ws/dqi_build/r7/ab/ab_data.npz'))
    else:
        words, wts = make_syn(nsyn, np.random.default_rng(seed), low=(mode == 'low'))
    S = words.shape[1]
    print('mode %s, shots %d, weights %d..%d' % (mode, S, wts.min(), wts.max()))
    Wr, cr, ref = sweep(words, (), Conv)
    print('reference: %s  (%.1f s)' % (dict(sorted(cr.items())), time.time() - t0))
    # classes from the reference's own flags
    interior = {m: ref.flags[m] & ~ref.term[m] for m in ref.flags}
    for steps in configs:
        Wp, cp, pc = sweep(words, steps, ConvPaired)
        mism = np.any(Wp != Wr, axis=0)
        lowest = min(steps) - 1
        near = np.zeros(S, dtype=bool)
        for m in range(lowest, 18):
            near |= interior.get(m, np.zeros(S, dtype=bool))
        dp = cp['payload products'] - cr['payload products']
        dsc = cp['payload scalings'] - cr['payload scalings']
        dmk = cp['payload masked products'] - cr['payload masked products']
        dcx = cp['payload controlled XOR words']
        print('steps %s: payload products %+d, scalings %+d, masked products %+d, controlled XOR words %d '
              '(both payloads, one convert)' % (steps, dp, dsc, dmk, dcx))
        print('   mismatching shots %d / %d; of them with an interior gap at levels %d..17: %d; '
              'interior-gap shots total %d; weights of mismatches: %s'
              % (mism.sum(), S, lowest, (mism & near).sum(), near.sum(),
                 dict(sorted(Counter(wts[mism].tolist()).items()))))
        other = mism & ~near
        if other.any():
            print('   mismatches WITHOUT an interior gap near the blocks: %d, weights %s'
                  % (other.sum(), dict(sorted(Counter(wts[other].tolist()).items()))))
        lv_mis, lv_all = Counter(), Counter()
        for m in range(lowest, 18):
            gm = interior.get(m, np.zeros(S, dtype=bool))
            lv_all[m] = int(gm.sum())
            lv_mis[m] = int((gm & mism).sum())
        print('   interior gaps by kit level (mismatching / all): %s'
              % ' '.join('%d:%d/%d' % (m, lv_mis[m], lv_all[m]) for m in range(17, lowest - 1, -1) if lv_all[m]))
    print('done %.1f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
