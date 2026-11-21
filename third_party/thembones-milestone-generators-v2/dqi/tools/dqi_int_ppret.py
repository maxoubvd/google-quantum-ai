#!/usr/bin/env python3
"""dqi_int_ppret.py -- workstream dqi_int (B5): rewrite pass P2 "retained products" on an emitted .kmb.

Pairs (first i, second j) of unconditioned free CCX computing the same product (exact over many runs:
tools/dqi_int_pprep.py) get a scratch wire s that is structurally zero and untouched in the ORIGINAL circuit over the
whole gap (reset by HMR / R or never touched, next event after j); greedy, shortest gaps first, each op in one pair,
no two pairs overlapping on a wire. Then
  i: CCX(a, b, t)   -> CCX(a, b, s); CX(s, t)
  j: CCX(a', b', t'') -> CX(s, t''); HMR(s) -> m (new bit); CZ(a', b') if m
-1 Toffoli per pair, no new qubit (s is an existing idle wire), one new classical bit per pair.
Usage: python dqi_int_ppret.py <in.kmb> <pairs.npz> <out.kmb> [max pairs=all] [lo:hi of the selected list]
"""
import bisect
import dataclasses
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402

src, pairs_npz, dst = sys.argv[1], sys.argv[2], sys.argv[3]
maxp = int(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] != 'all' else None
rng = sys.argv[5] if len(sys.argv) > 5 else ':'
c = kmb.read(src)
P = np.load(pairs_npz)
first, second = P['pair_first'].astype(np.int64), P['pair_second'].astype(np.int64)
ops = np.asarray(c.ops)
n = len(ops)
base = ops & 0x7F
nq = (ops & 0x30) >> 4
W = np.full((n, 3), -1, dtype=np.int64)
for k, arrs in ((3, c.q3), (2, c.q2), (1, (c.q1,))):
    sel = np.flatnonzero(nq == k)
    for j in range(k):
        W[sel, j] = np.asarray(arrs[j], dtype=np.int64)
ev_op = np.repeat(np.arange(n, dtype=np.int64), 3)
ev_w = W.reshape(-1)
keep = ev_w >= 0
ev_op, ev_w = ev_op[keep], ev_w[keep]
o = np.lexsort((ev_op, ev_w))
ev_op, ev_w = ev_op[o], ev_w[o]
reset = np.isin(base[ev_op], (0x50, 0x12))
IS, IE = [], []                                   # per wire: zero-idle intervals (start, end), disjoint, sorted
ZI = np.load(sys.argv[6]) if len(sys.argv) > 6 else None   # tools/dqi_int_ppzero.py intervals (else structural)
for w in range(c.num_qubits):
    if ZI is not None:
        sel_ = ZI['wire'] == w
        s = ZI['start'][sel_].astype(np.int64)
        nxt = ZI['end'][sel_].astype(np.int64)
        o_ = np.argsort(s)
        IS.append(s[o_])
        IE.append(nxt[o_])
        continue
    lo_ = np.searchsorted(ev_w, w, side='left')
    hi_ = np.searchsorted(ev_w, w, side='right')
    e = ev_op[lo_:hi_]
    r = reset[lo_:hi_]
    s = np.concatenate([[-1], e[r]]).astype(np.int64)
    nxt = np.concatenate([e[:1] if len(e) else np.array([n]), np.concatenate([e[1:], [n]])[r]]).astype(np.int64)
    IS.append(s)
    IE.append(nxt)
used = {w: [] for w in range(c.num_qubits)}       # assigned (i, j) per wire, sorted by i
used_ops = set()
order = np.argsort(second - first, kind='stable')
chosen = []
for p in order:
    i, j = int(first[p]), int(second[p])
    if i in used_ops or j in used_ops:
        continue
    for w in range(c.num_qubits):
        k_ = bisect.bisect_right(IS[w], i - 1) - 1    # interval with start < i
        if k_ < 0 or IE[w][k_] <= j:
            continue
        L = used[w]
        q = bisect.bisect_left(L, (i, j))
        if (q > 0 and L[q - 1][1] >= i) or (q < len(L) and L[q][0] <= j):
            continue
        L.insert(q, (i, j))
        chosen.append((i, j, w))
        used_ops.add(i)
        used_ops.add(j)
        break
    if maxp is not None and len(chosen) >= maxp:
        break
chosen.sort()
lo, hi = rng.split(':')
chosen = chosen[int(lo) if lo else 0:int(hi) if hi else len(chosen)]
print('pairs offered %d, scheduled %d (rewriting %d)' % (len(first), len(used_ops) // 2, len(chosen)), flush=True)
nb = c.num_bits
inserts = {}
delete = set()
for idx, (i, j, s) in enumerate(chosen):
    a, b, t = (int(x) for x in W[i])
    a2, b2, t2 = (int(x) for x in W[j])
    m = nb + idx
    inserts[i] = [(0x31, (a, b, s), 0, 0), (0x21, (s, t), 0, 0)]
    inserts[j] = [(0x21, (s, t2), 0, 0), (0x50, (s,), m, 0), (0xA0, (a2, b2), 0, m)]
    delete.add(i)
    delete.add(j)
out = kw.edit(c, delete, inserts)
out = dataclasses.replace(out, num_bits=nb + len(chosen))
nqub = kw.write(dst, out)
print('qubits %d, bits %d -> %d, static Toffoli-like %d -> %d' % (
    nqub, nb, out.num_bits, int(((ops & 0x7C) == 0x30).sum()), int(((out.ops & 0x7C) == 0x30).sum())))
np.save(dst + '.pairs.npy', np.array(chosen, dtype=np.int64))
