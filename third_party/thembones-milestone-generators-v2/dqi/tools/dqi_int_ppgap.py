#!/usr/bin/env python3
"""dqi_int_ppgap.py -- workstream dqi_int (B5): for repeated-product pairs (first, second op) from
tools/dqi_int_ppfp.py, is there a wire that is structurally zero and untouched over the whole gap (reset by HMR / R,
or never touched yet, and no event until after `second`)? That is the scratch a retained product needs (no rewrite).
Usage: python dqi_int_ppgap.py <circuit.kmb> <fp.npz>
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402

c = kmb.read(sys.argv[1])
fp = np.load(sys.argv[2])
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
reset = np.isin(base[ev_op], (0x50, 0x12))          # HMR, R leave the wire at 0
# zero-idle intervals per wire: (start, end) = (reset op or -1, next event op or n)
starts, ends, wires = [], [], []
for w in range(c.num_qubits):
    lo = np.searchsorted(ev_w, w, side='left')
    hi = np.searchsorted(ev_w, w, side='right')
    e = ev_op[lo:hi]
    r = reset[lo:hi]
    s = np.concatenate([[-1], e[r]])
    nxt = np.concatenate([e[:1] if len(e) else [n], np.concatenate([e[1:], [n]])[r]])
    starts.append(s)
    ends.append(nxt)
    wires.append(np.full(len(s), w))
S = np.concatenate(starts)
E = np.concatenate(ends)
first, second = fp['pair_first'], fp['pair_second']
print('pairs', len(first), 'zero-idle intervals', len(S), flush=True)
# a pair is covered if some interval has S < first and E > second; sort intervals by start, prefix max of end
oi = np.argsort(S)
S2, E2 = S[oi], E[oi]
pm = np.maximum.accumulate(E2)
idx = np.searchsorted(S2, first, side='left') - 1
covered = (idx >= 0) & (pm[np.maximum(idx, 0)] > second)
gap = second - first
print('pairs with a zero-idle wire over the gap: %d of %d; gap median %d (covered: %d)' % (
    int(covered.sum()), len(first), int(np.median(gap)), int(np.median(gap[covered])) if covered.any() else -1))
