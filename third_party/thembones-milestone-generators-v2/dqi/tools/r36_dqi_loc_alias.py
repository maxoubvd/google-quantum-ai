#!/usr/bin/env python3
"""r36_dqi_loc_alias.py -- R36 designer DQI-LOC: the "cross-word alias" form of the converter's locator update.

The kit's descent forms the quotient (q0, q1) in place on the dividend's two head words by scale / cross / scale and
then runs shifted_conv: the low m - 1 coefficients of coeff(x) g'(x) by m evaluation products (the full
m-coefficient product; its head coefficient is computed and discarded). That head coefficient is exactly the cross
product X = q1 * (g[-2] or g[-3] by F), which the descent has just XORed into the q0 slot. So the product may be
taken at m - 1 points only (i.e. modulo Pi(x) = prod (x - t)), with the alias of X x^(m-1) repaired by constants:

    out ^= corr * slot      (before the cross: slot = f_(m-1))
    slot ^= q1 * c_sel      (the kit's cross, unchanged)
    out ^= corr * slot      (after the cross: net corr * X)
    ...scale q0...          (unchanged)
    out ^= interpolation of coeff(t) g'(t) over the m - 1 points      (m - 1 products instead of m)

corr_j = sum_t W[j][t] t^(m-1) over the m - 1 points (the constants of high_full_restore's own alias repair).

Part 1 (kernel): the direct kernels' arithmetic (pivot, selected cross, shifted pair, evaluation conv) old against
new against the schoolbook product, on random and degenerate words, every level m, both charts F.
Part 2 (sweep): the kit's reference sweep (class Conv copied from tools/dqi_k1_conv.py, which copied it from
tools/r33_dqi_conv_rp_kitform.py) with the alias form in the descent AND in the in-place return, on the real
(A, B) shots of the r7 census and on synthetic shots: every word of the final state against the reference.
Part 3 (count): the line's levels at the kit's prices.

  python r36_dqi_loc_alias.py kernel [shots=4000] [seed=1]
  python r36_dqi_loc_alias.py sweep real|syn|low [n=6000] [seed=1]
  python r36_dqi_loc_alias.py count
"""
import os
import sys
from collections import Counter

import numpy as np

MOD = 4179          # x^12 + x^6 + x^4 + x + 1 (the kit's GF(2^12))
ORD = 4095
I64 = np.int64
EXP = np.zeros(2 * ORD, dtype=I64)
LOG = np.zeros(4096, dtype=I64)
_x = 1
for _i in range(ORD):
    EXP[_i] = _x
    LOG[_x] = _i
    _x <<= 1
    if _x & 4096:
        _x ^= MOD
EXP[ORD:] = EXP[:ORD]


def fmul(a, b):
    a = np.asarray(a, dtype=I64)
    b = np.asarray(b, dtype=I64)
    return np.where((a == 0) | (b == 0), 0, EXP[LOG[a] + LOG[b]])


def smul(c, a):
    """classical constant c times word array a (a free linear map in the circuit)"""
    c = int(c)
    a = np.asarray(a, dtype=I64)
    if c == 0:
        return np.zeros_like(a)
    return np.where(a == 0, 0, EXP[LOG[a] + LOG[c]])


def fsafe(a):
    return np.where(a == 0, 1, a)


def finv_safe(a):
    return EXP[(ORD - LOG[fsafe(a)]) % ORD]


def fsqrt(a):
    a = np.asarray(a, dtype=I64)
    return np.where(a == 0, 0, EXP[(LOG[a] * 2048) % ORD])


def s1(a, b):
    a, b = int(a), int(b)
    return 0 if a == 0 or b == 0 else int(EXP[LOG[a] + LOG[b]])


def spow(t, j):
    if j == 0:
        return 1
    if t == 0:
        return 0
    return int(EXP[(LOG[t] * j) % ORD])


_VW = {}


def vandermonde(n):
    """V[t][j] = t^j for the points t = 0..n-1 (as field elements), W = V^-1 (the kit's convention)."""
    if n not in _VW:
        V = [[spow(t, j) for j in range(n)] for t in range(n)]
        M = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(V)]
        for col in range(n):
            piv = next(r for r in range(col, n) if M[r][col])
            M[col], M[piv] = M[piv], M[col]
            ic = int(EXP[(ORD - LOG[M[col][col]]) % ORD])
            M[col] = [s1(ic, v) for v in M[col]]
            for r in range(n):
                if r != col and M[r][col]:
                    f = M[r][col]
                    M[r] = [v ^ s1(f, w) for v, w in zip(M[r], M[col])]
        _VW[n] = (V, [row[n:] for row in M])
    return _VW[n]


def alias_constants(n_points, power):
    """corr_j = sum_t W[j][t] t^power over the n_points points: the coefficients of x^power mod prod (x - t)."""
    V, W = vandermonde(n_points)
    return [int(np.bitwise_xor.reduce([s1(W[j][t], spow(t, power)) for t in range(n_points)]))
            for j in range(n_points)]


def conv_eval(g, coeff, N, nout, cnt=None, tag='locator products'):
    """The kit's conv: sum over the points t < N of W[j][t] * (g(t) * coeff(t)), rows j < nout. N products."""
    V, W = vandermonde(N)
    out = [np.zeros_like(g[0]) for _ in range(nout)]
    for t in range(N):
        x = np.zeros_like(g[0])
        for j, w in enumerate(g):
            x ^= smul(spow(t, j), w)
        y = np.zeros_like(g[0])
        for j, w in enumerate(coeff):
            y ^= smul(spow(t, j), w)
        p = fmul(x, y)
        if cnt is not None:
            cnt[tag] += 1
        for j in range(nout):
            out[j] ^= smul(W[j][t], p)
    return out


def school_low(g, coeff, nout):
    out = [np.zeros_like(g[0]) for _ in range(nout)]
    for i, a in enumerate(g):
        for j, c in enumerate(coeff):
            if i + j < nout:
                out[i + j] ^= fmul(a, c)
    return out


# ----------------------------------------------------------------------------------------------------- part 1
def kernel_case(m, S, rng, shifted=True):
    """One level of the direct descent (direct_pre, not undo), locator part only. Returns mismatch counts."""
    f = rng.integers(0, 4096, size=(m + 1, S)).astype(I64)
    g = rng.integers(0, 4096, size=(m, S)).astype(I64)
    r = rng.random(S)
    g[-1, r < 0.35] = 0                      # zero top word: F = 1
    if m >= 3:
        g[-2, (r < 0.10)] = 0                # top and second both zero: pivot stays zero (safe inverse)
    g[:, r < 0.03] = 0                       # early end: g entirely zero
    f[-1, (r > 0.90)] = 0                    # zero dividend head
    f[:, (r > 0.97)] = 0
    F = (g[-1] == 0) if shifted else np.zeros(S, dtype=bool)
    gs = g.copy()                            # pivot(b, g, F): swap g[-1], g[-2] under F
    gs[-1] = np.where(F, g[-2], g[-1])
    gs[-2] = np.where(F, g[-1], g[-2])
    a = gs[-1]
    ai = finv_safe(a)                        # inverse of the clipped pivot (the kit's Z flag)
    q1 = fmul(f[-1], ai)                     # scale(q1, a, ai, True)
    if m >= 3:
        csel = np.where(F, gs[-3], gs[-2])   # dk.selected(g[-2], g[-3], F)
    else:
        csel = gs[-2]
    gp = [gs[i] for i in range(m - 1)]       # g[:-1] of the swapped array
    slot0 = f[-2]

    def coeff_of(q0):
        z = np.zeros(S, dtype=I64)
        return [np.where(F, z, q0), np.where(F, q0, q1), np.where(F, q1, z)]     # shifted_pair

    # old kernel: cross, scale, conv at m points truncated to m - 1 rows
    cnt_old, cnt_new = Counter(), Counter()
    slot = slot0 ^ fmul(q1, csel)
    q0 = fmul(slot, ai)
    upd = conv_eval(gp, coeff_of(q0), m, m - 1, cnt_old)
    old = [f[j] ^ upd[j] for j in range(m - 1)]
    # exact: schoolbook low product
    ex = school_low(gp, coeff_of(q0), m - 1)
    exact = [f[j] ^ ex[j] for j in range(m - 1)]
    # new kernel: alias repair around the cross, conv at m - 1 points
    corr = alias_constants(m - 1, m - 1)
    new = [f[j] ^ smul(corr[j], slot0) for j in range(m - 1)]
    slot_n = slot0 ^ fmul(q1, csel)
    new = [new[j] ^ smul(corr[j], slot_n) for j in range(m - 1)]
    q0n = fmul(slot_n, ai)
    upd = conv_eval(gp, coeff_of(q0n), m - 1, m - 1, cnt_new)
    new = [new[j] ^ upd[j] for j in range(m - 1)]
    bad_old = np.zeros(S, dtype=bool)
    bad_new = np.zeros(S, dtype=bool)
    for j in range(m - 1):
        bad_old |= old[j] != exact[j]
        bad_new |= new[j] != exact[j]
    bad_new |= q0n != q0
    return (int(bad_old.sum()), int(bad_new.sum()), int(F.sum()), int((a == 0).sum()),
            cnt_old['locator products'], cnt_new['locator products'])


def kernel_restore_case(m, S, rng):
    """The in-place return (direct_pre undo in its ordinary form): conv, scale q0, cross, scale q1, with the alias
    repair around the cross; against the schoolbook restore of the dividend."""
    g = rng.integers(0, 4096, size=(m, S)).astype(I64)
    r = rng.random(S)
    g[-1, r < 0.35] = 0
    if m >= 3:
        g[-2, (r < 0.10)] = 0
    g[:, r < 0.03] = 0
    F = g[-1] == 0
    gs = g.copy()
    gs[-1] = np.where(F, g[-2], g[-1])
    gs[-2] = np.where(F, g[-1], g[-2])
    a = fsafe(gs[-1])
    csel = np.where(F, gs[-3], gs[-2]) if m >= 3 else gs[-2]
    gp = [gs[i] for i in range(m - 1)]
    rem = rng.integers(0, 4096, size=(m - 1, S)).astype(I64)
    q0 = rng.integers(0, 4096, size=S).astype(I64)
    q1 = rng.integers(0, 4096, size=S).astype(I64)
    q1[r > 0.9] = 0
    z = np.zeros(S, dtype=I64)
    coeff = [np.where(F, z, q0), np.where(F, q0, q1), np.where(F, q1, z)]
    ex = school_low(gp, coeff, m - 1)
    f_exact = [rem[j] ^ ex[j] for j in range(m - 1)]
    lo_exact = fmul(q0, a) ^ fmul(q1, csel)
    hi_exact = fmul(q1, a)
    cnt = Counter()
    corr = alias_constants(m - 1, m - 1)
    upd = conv_eval(gp, coeff, m - 1, m - 1, cnt)
    new = [rem[j] ^ upd[j] for j in range(m - 1)]
    slot = fmul(q0, a)                                   # scale(q0, a, ai, False)
    new = [new[j] ^ smul(corr[j], slot) for j in range(m - 1)]
    slot = slot ^ fmul(q1, csel)                         # cross
    new = [new[j] ^ smul(corr[j], slot) for j in range(m - 1)]
    hi = fmul(q1, a)                                     # scale(q1, a, ai, False)
    bad = (slot != lo_exact) | (hi != hi_exact)
    for j in range(m - 1):
        bad |= new[j] != f_exact[j]
    return int(bad.sum()), cnt['locator products']


def run_kernel(S, seed):
    rng = np.random.default_rng(seed)
    print('part 1: direct descent kernel, %d shots per level, seed %d' % (S, seed))
    print('  m   F=1  a=0   old!=exact  new!=exact  products old -> new   return(in place) bad, products')
    tot = [0, 0, 0]
    for m in range(3, 18):
        bo, bn, nf, na, po, pn = kernel_case(m, S, rng)
        br, pr = kernel_restore_case(m, S, rng)
        tot[0] += bo
        tot[1] += bn
        tot[2] += br
        print('  %2d  %4d %4d   %9d  %9d   %2d -> %2d                 %d, %d' % (m, nf, na, bo, bn, po, pn, br, pr))
    bo, bn, nf, na, po, pn = kernel_case(2, S, rng, shifted=False)
    print('   2  (no shift)  %9d  %9d   %2d -> %2d' % (bo, bn, po, pn))
    print('KERNEL %s: old mismatches %d, new mismatches %d, in-place return mismatches %d'
          % ('PASS' if tot[0] + tot[1] + tot[2] + bo + bn == 0 else 'FAIL', tot[0] + bo, tot[1] + bn, tot[2]))


# ----------------------------------------------------------------------------------------------------- part 2
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
    """The kit's reference sweep (copied from tools/dqi_k1_conv.py; probes removed)."""

    def __init__(self, sim):
        self.s = sim

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

    def head_down(self, f, g):
        s = self.s
        a, q0, q1 = g[-1], f[-2], f[-1]
        s.scale(q1, a, True, 'locator')
        s.mul(q1, g[-2], q0, 'locator')
        s.scale(q0, a, True, 'locator')
        self.conv(f, g, q0, q1)

    def head_up(self, f, g):
        s = self.s
        a, q0, q1 = g[-1], f[-2], f[-1]
        self.conv(f, g, q0, q1)
        s.scale(q0, a, False, 'locator')
        s.mul(q1, g[-2], q0, 'locator')
        s.scale(q1, a, False, 'locator')

    def pre(self, f, g, P, parent):
        s = self.s
        m = len(f) - 1
        n = m - 1
        a = g[-1]
        for v in P:
            s.scale(v[-1], a, True, 'payload')
        for v in P:
            for i in range(n):
                s.mul(g[i], v[-1], v[i], 'payload')
        self.head_down(f, g)

    def q_post(self, q, P, parent, m):
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

    def post_look(self, f, g, P, parent):
        s = self.s
        m = len(f) - 1
        n = m - 1
        a, q1 = g[-1], f[-1]
        self.q_post(q1, P, parent, m)
        for v in P:
            for i in range(n):
                s.mul(g[i], v[i], v[-1], 'payload')
        for v in P:
            s.scale(v[-1], a, True, 'payload')
        self.head_up(f, g)

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

    def look_solve(self, f, g, P, parent=None):
        W = self.s.W
        m = len(f) - 1
        if m == 1:
            return self.base(f, g, P, parent)
        F = W[g[-1]] == 0
        self.rotate(g, F)
        self.pre(f, g, P, parent)
        self.permute(f, g, F)
        child = ([v[-1] for v in P], f[-1], f[-2])
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
        self.look_solve(O + [lead], E[:-1], payloads)
        for i in range(len(O)):
            W[E[i]] ^= fmul(W[E[-1]], W[O[i]])


class ConvAlias(Conv):
    """The same sweep with the cross-word alias form of the locator update, descent and in-place return."""

    def conv_points(self, f, g, q0, q1):
        W = self.s.W
        m = len(f) - 1
        upd = conv_eval([W[x] for x in g[:-1]], [W[q0], W[q1]], m - 1, m - 1, self.s.cnt)
        for j in range(m - 1):
            W[f[j]] ^= upd[j]

    def repair(self, f, slot):
        W = self.s.W
        m = len(f) - 1
        corr = alias_constants(m - 1, m - 1)
        for j in range(m - 1):
            W[f[j]] ^= smul(corr[j], W[slot])          # free: a constant multiple of a live word

    def head_down(self, f, g):
        s = self.s
        a, q0, q1 = g[-1], f[-2], f[-1]
        s.scale(q1, a, True, 'locator')
        self.repair(f, q0)
        s.mul(q1, g[-2], q0, 'locator')
        self.repair(f, q0)
        s.scale(q0, a, True, 'locator')
        self.conv_points(f, g, q0, q1)

    def head_up(self, f, g):
        s = self.s
        a, q0, q1 = g[-1], f[-2], f[-1]
        self.conv_points(f, g, q0, q1)
        s.scale(q0, a, False, 'locator')
        self.repair(f, q0)
        s.mul(q1, g[-2], q0, 'locator')
        self.repair(f, q0)
        s.scale(q1, a, False, 'locator')


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
    """Random (A, B): B = z^t B~, B~ monic of a weight 14..35 (low: 0..25, many early ends)."""
    wts = rng.integers(0, 26, size=n) if low else rng.integers(14, 36, size=n)
    words = np.zeros((70, n), dtype=I64)
    words[:35] = rng.integers(0, 4096, size=(35, n))
    for s in range(n):
        w = int(wts[s])
        t = 35 - w
        if w == 0:
            continue
        cf = rng.integers(0, 4096, size=w)
        cf[0] = rng.integers(1, 4096)
        if rng.random() < 0.3:                           # force interior zero coefficients (gap charts)
            cf[rng.integers(0, w, size=max(1, w // 4))] = 0
            cf[0] = max(1, int(cf[0]))
        words[35 + t:35 + t + w, s] = cf
    return words, wts


def sweep(words, cls):
    A, B = list(range(0, 35)), list(range(35, 70))
    lead, nw = 70, 71
    S = words.shape[1]
    W = np.zeros((nw, S), dtype=I64)
    W[:70] = fsqrt(words)
    sim = Sim(W)
    E, O = B[::2], B[1::2]
    u, w, top = A[1::2], A[:34:2], A[-1]
    for i in range(17):
        W[w[i]] ^= fmul(W[top], W[O[i]])
    c = cls(sim)
    c.run(O, E, [u, w], lead)
    return W, sim.cnt


def run_sweep(mode, n, seed):
    rng = np.random.default_rng(seed)
    if mode == 'real':
        words, wts = load_real(os.path.expanduser('~/qcc/ws/dqi_build/r7/ab/ab_data.npz'))
    else:
        words, wts = make_syn(n, rng, low=(mode == 'low'))
    S = words.shape[1]
    W0, c0 = sweep(words, Conv)
    W1, c1 = sweep(words, ConvAlias)
    bad = np.any(W0 != W1, axis=0)
    print('part 2: sweep, mode %s, shots %d, weights %d..%d' % (mode, S, wts.min(), wts.max()))
    print('  shots whose final 71 words differ from the reference sweep: %d of %d' % (bad.sum(), S))
    print('  locator products per sweep (descent + return): reference %d, alias %d'
          % (c0['locator products'], c1['locator products']))
    print('  all other counts equal: %s' % all(c0[k] == c1[k] for k in c0 if k != 'locator products'))
    print('SWEEP %s' % ('PASS' if bad.sum() == 0 else 'FAIL'))


# ----------------------------------------------------------------------------------------------------- part 3
def run_count():
    """The line's shifted_conv calls at the kit's prices [A020: 42 m + 36 per call, m = 6..17, both directions]."""
    old = sum(42 * m + 36 for m in range(6, 18))
    new = sum(42 * (m - 1) + 36 for m in range(6, 18))
    print('part 3: shifted_conv, levels 6..17, one direction: %d -> %d (%d); two directions: %d'
          % (old, new, new - old, 2 * (new - old)))
    low_old = {3: 99, 4: 153}
    low_new = {3: 84, 4: 126}
    d = sum(low_new[m] - low_old[m] for m in low_old)
    print('low levels 3, 4 (truncated_small 99 / 153 against 2 / 3 products): %d per descent, %d for two' % (d, 2 * d))
    ret_old = sum(42 * m + 156 for m in range(6, 17))
    ret_new = sum(42 * m + 132 for m in range(6, 17))
    print('reverse return (R2), levels 6..16: high_full_restore + high_close_heads %d -> in-place alias form %d (%d)'
          % (ret_old, ret_new, ret_new - ret_old))


def main():
    a = sys.argv[1:]
    mode = a[0] if a else 'kernel'
    if mode == 'kernel':
        run_kernel(int(a[1]) if len(a) > 1 else 4000, int(a[2]) if len(a) > 2 else 1)
    elif mode == 'sweep':
        run_sweep(a[1] if len(a) > 1 else 'real', int(a[2]) if len(a) > 2 else 6000, int(a[3]) if len(a) > 3 else 1)
    else:
        run_count()


if __name__ == '__main__':
    main()
