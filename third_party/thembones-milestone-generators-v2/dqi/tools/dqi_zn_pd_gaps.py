#!/usr/bin/env python3
"""dqi_zn_pd_gaps.py -- builder Z1: the paired converter model against the kit model on every gap pattern:
a forced zero pivot at each level 2..17 alone, and every pair of gaps d = 2, 3, 4 levels apart (B5's gap_charts:
a gap at level m = the remainder after degree m has degree m - 2), plus padded synthetic charts of every weight
18..35 (B = z^(35-w) B_w, A = z^(35-w) A_w, the shape of the verifier's charts).
Usage: python dqi_zn_pd_gaps.py <class> [upper=17,15,13,11] [shots=256] [seed=1]
"""
import collections
import os
import sys

import numpy as np

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import dqi_zn_pd as PD  # noqa: E402
import dqi_int_convtest as CT  # noqa: E402  (gap_charts; B5, read-only)


def padded(rng, n, w):
    """Charts shaped like the decode's at weight w: B = z^s B_w (B_w monic of degree w with distinct nonzero roots),
    A = z^s A_w with A_w = sum Y_k prod_(j != k) (z - x_j) (the numerator of sum Y_k / (z - x_k))."""
    s = 35 - w
    A = np.zeros((n, 35), dtype=np.int64)
    B = np.zeros((n, 35), dtype=np.int64)
    for i in range(n):
        xs = rng.choice(np.arange(1, 4096), size=w, replace=False)
        ys = rng.integers(1, 4096, size=w)
        b = np.array([1], dtype=np.int64)
        for x in xs:
            b = np.concatenate([[0], b]) ^ np.concatenate([PD.mul(np.int64(x), b), [0]])
        a = np.zeros(w, dtype=np.int64)
        for k in range(w):
            p = np.array([1], dtype=np.int64)
            for j in range(w):
                if j != k:
                    p = np.concatenate([[0], p]) ^ np.concatenate([PD.mul(np.int64(xs[j]), p), [0]])
            a[:len(p)] ^= PD.mul(np.int64(ys[k]), p)
        full_b = np.zeros(36, dtype=np.int64)
        full_b[s:s + w + 1] = b
        assert full_b[35] == 1
        B[i] = full_b[:35]
        A[i, s:s + w] = a
    return A, B


def main():
    cname = sys.argv[1]
    upper = tuple(int(x) for x in (sys.argv[2] if len(sys.argv) > 2 else '17,15,13,11').split(','))
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 256
    rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    cls = getattr(PD, cname)
    kw = {} if cname in ('Kit', 'KitQB') else dict(upper=upper, route='full')
    CT.N = n
    tests = [('gap %d' % m, {m}) for m in range(17, 1, -1)]
    for d in (2, 3, 4):
        for m in range(17, 2 + d, -1):
            tests.append(('gaps %d,%d' % (m, m - d), {m, m - d}))
    bad = collections.OrderedDict()
    tot = 0
    for label, gaps in tests:
        A, B = CT.gap_charts(rng, n, gaps)
        O, E, u, w = PD.bez_inputs(A, B)
        eq, _, _ = PD.compare(O, E, u, w, cls=cls, **kw)
        tot += 1
        if not eq.all():
            bad[label] = int((~eq).sum())
    print('%s upper %s: %d gap patterns x %d shots; patterns with a mismatch: %d %s' % (
        cname, upper, tot, n, len(bad), dict(bad)), flush=True)
    badw = {}
    for w in range(18, 36):
        A, B = padded(rng, max(16, n // 8), w)
        O, E, u, w_ = PD.bez_inputs(A, B)
        eq, _, _ = PD.compare(O, E, u, w_, cls=cls, **kw)
        if not eq.all():
            badw[w] = int((~eq).sum())
    print('padded synthetic charts, weights 18..35 (%d shots each): mismatches %s' % (max(16, n // 8), badw), flush=True)


if __name__ == '__main__':
    main()
