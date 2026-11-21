#!/usr/bin/env python3
"""dqi_k1_conv.py -- K1, row C-c of the idea-10 table: the converter's ascent inverses from a carried chain.

Word-level model of the kit's half-Euclid sweep with its two payloads (classical, vectorised over shots). The
sweep (look_solve / pre / q_post / post_look / base / rotate / permute / conv) is COPIED from
tools/r33_dqi_conv_rp_kitform.py (R33 DQI-CONV re-pricer; a port of the kit's
playbook/kit/decoder/exp/safe_scalar/converter.py), reduced to the reference path, so this file imports no other
seat's live file.

Row C-c: the root's lead word is 1, so 1/a_k = q1_k * (1/a_(k+1)) with q1_k = f_top / a_k the quotient top the
descent leaves in f[-1]. Chain H_(K+1) = 1, H_k = q1_k H_(k+1) (top-down, one product per level): it gives every
1/a_k (the kit's `last` inverse) and 1/q1_k = a_k H_(k+1) (the kit's q_post inverse) without an inverse.
Bottom-up (convert's ascent): H_(k+1) = H_k * (1/q1_k), one product per level beside the q_post inverse.
This tool checks, per level and shot, whether the chain's values equal what the kit uses (finv_safe of the live
words) on every shot, generic and gap charts, and counts the operations of both forms.
  python dqi_k1_conv.py [n_syn=4000] [seed=1] [low=0]
"""
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import dqi_arch_model as M  # noqa: E402  (field tables only)

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
    """The kit's sweep (copied, reference path) with probes for row C-c."""

    def __init__(self, sim):
        self.s = sim
        self.flags = {}
        self.term = {}
        self.desc = {}          # level -> (a_k, q1_k) read at the end of pre (after rotate)
        self.asc = {}           # level -> (a used by last, q1 used by q_post, kit's 1/q1, kit's 1/a)

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
        s = self.s
        m = len(f) - 1
        n = m - 1
        a, q0, q1 = g[-1], f[-2], f[-1]
        s.cnt['inverses (pre)'] += 1
        for v in P:
            s.scale(v[-1], a, True, 'payload')
        for v in P:
            for i in range(n):
                s.mul(g[i], v[-1], v[i], 'payload')
        s.scale(q1, a, True, 'locator')
        s.mul(q1, g[-2], q0, 'locator')
        s.scale(q0, a, True, 'locator')
        self.conv(f, g, q0, q1)
        self.desc[m] = (s.W[a].copy(), s.W[q1].copy())

    def q_post(self, q, P, parent, m):
        s, W = self.s, self.s.W
        s.cnt['inverses (q_post)'] += 1
        self.asc.setdefault(m, {})['q1'] = W[q].copy()
        self.asc[m]['qi_kit'] = finv_safe(W[q])
        if parent is None:
            for v in P:
                s.scale(v[-1], q, True, 'payload')
            return
        tops, pw, c = parent
        F = W[q] == 0
        self.asc[m]['gapF'] = F.copy()
        s.cswap(F, q, pw, 'route')
        for v, t in zip(P, tops):
            s.cswap(F, v[-1], t, 'route')
        self.asc[m]['qi_kit_after_swap'] = finv_safe(W[q])
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
        for v in P:
            for i in range(n):
                s.mul(g[i], v[i], v[-1], 'payload')
        s.cnt['inverses (last)'] += 1
        self.asc[m]['a'] = W[a].copy()
        self.asc[m]['ai_kit'] = finv_safe(W[a])
        for v in P:
            s.scale(v[-1], a, True, 'payload')
        self.conv(f, g, q0, q1)
        s.scale(q0, a, False, 'locator')
        s.mul(q1, g[-2], q0, 'locator')
        s.scale(q1, a, False, 'locator')

    def base(self, f, g, P, parent):
        s, W = self.s, self.s.W
        a, q = g[0], f[-1]
        s.cnt['inverses (base)'] += 1
        for v in P:
            s.scale(v[0], a, True, 'payload')
        s.scale(q, a, True, 'locator')
        self.desc[1] = (W[a].copy(), W[q].copy())
        self.q_post(q, P, parent, 1)
        self.asc[1]['a'] = W[a].copy()
        self.asc[1]['ai_kit'] = finv_safe(W[a])
        for v in P:
            s.scale(v[0], a, True, 'payload')
        s.scale(q, a, False, 'locator')

    def look_solve(self, f, g, P, parent=None):
        W = self.s.W
        m = len(f) - 1
        if m == 1:
            return self.base(f, g, P, parent)
        F = W[g[-1]] == 0
        self.flags[m] = F
        self.term[m] = np.all(W[g] == 0, axis=0)
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


def make_syn(n, rng, low=False):
    """Random (A, B) with the verifier's weight law (copied from the kit-form model): B = z^t B~, B~ monic."""
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
        words[35 + t:35 + t + w, s] = cf
    return words, wts


def sweep(words):
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
    c = Conv(sim)
    c.run(O, E, [u, w], lead)
    return W, sim.cnt, c


def chain_check(c, S):
    """Compare the chain's values with the kit's inverses, level by level (m = 17 root ... 1 base)."""
    gap_any = np.zeros(S, dtype=bool)
    for m in c.flags:
        gap_any |= c.flags[m]
    for m in c.asc:
        if 'gapF' in c.asc[m]:
            gap_any |= c.asc[m]['gapF']
    H = np.ones(S, dtype=I64)                          # H_(18) = 1/lead = 1
    bad_last = np.zeros(S, dtype=bool)
    bad_qi = np.zeros(S, dtype=bool)
    rows = []
    for m in range(17, 0, -1):                         # top-down: unconvert's reverse ascent
        A = c.asc[m]
        qi_chain = fmul(fsafe(A['a']), H)              # 1/q1_m = a~_m H_(m+1), a~ = a xor [a = 0]
        H = fmul(fsafe(A['q1']), H)                    # H_m = q1~_m H_(m+1) = 1/a_m, q1~ = q1 xor [q1 = 0]
        ok_last = H == A['ai_kit']
        ok_qi = qi_chain == A['qi_kit']
        bad_last |= ~ok_last
        bad_qi |= ~ok_qi
        rows.append((m, int((~ok_last).sum()), int((~ok_qi).sum())))
    up_bad = np.zeros(S, dtype=bool)                   # bottom-up: convert's ascent, H_(m+1) = H_m / q1_m
    Hb = c.asc[1]['ai_kit']
    for m in range(1, 17):
        Hb = fmul(Hb, c.asc[m]['qi_kit'])
        up_bad |= Hb != c.asc[m + 1]['ai_kit']
    return rows, bad_last, bad_qi, up_bad, gap_any


def main():
    kw = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
    n = int(kw.get('n_syn', 4000))
    rng = np.random.default_rng(int(kw.get('seed', 1)))
    words, wts = make_syn(n, rng, low=bool(int(kw.get('low', 0))))
    W, cnt, c = sweep(words)
    rows, bad_last, bad_qi, up_bad, gap = chain_check(c, n)
    print('shots %d (weights %d..%d), gap shots (any zero pivot or q1) %d' % (n, wts.min(), wts.max(), gap.sum()))
    print('top-down chain (unconvert reverse ascent): shots where H_m != kit 1/a_m: %d (gap %d); where a_m H_(m+1)'
          ' != kit 1/q1_m: %d (gap %d)' % (bad_last.sum(), (bad_last & gap).sum(), bad_qi.sum(), (bad_qi & gap).sum()))
    print('bottom-up chain (convert ascent): shots where H_m * kit(1/q1_m) != kit 1/a_(m+1): %d (gap %d)'
          % (up_bad.sum(), (up_bad & gap).sum()))
    print('per level (m: last mismatches, q_post mismatches): ' + ', '.join('%d:%d/%d' % r for r in rows))
    print('kit operation counts (one sweep): ' + ', '.join('%s %d' % kv for kv in sorted(cnt.items())))


if __name__ == '__main__':
    main()
