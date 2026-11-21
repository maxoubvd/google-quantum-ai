#!/usr/bin/env python3
"""dqi_int_ppret2.py -- workstream dqi_int (B5): rewrite pass P2b, retained products with ROVING scratch.

As tools/dqi_int_ppret.py (pairs i < j of free CCX computing the same product, exact by tools/dqi_int_pprep.py), but
the retained copy s need not stay on one wire: it hops between wires that are idle and zero in the ORIGINAL circuit,
moved by a SWAP (Clifford, free) just before the original touches its current wire. A pair is admissible when at
every moment of its gap the original leaves a free wire not taken by another retained copy; pairs are taken greedily,
shortest gaps first (segment tree over the compressed timeline). No new qubit; one new classical bit per pair.
  i: CCX(a, b, t) -> CCX(a, b, s0); CX(s0, t)        hops: SWAP(s, s') before the original touches s
  j: CCX(a', b', t'') -> CX(s, t''); HMR(s) -> m; CZ(a', b') if m
Idle-zero intervals: tools/dqi_int_ppzero.py output (structural ones only with --structural).
Usage: python dqi_int_ppret2.py <in.kmb> <pairs.npz> <zero.npz> <out.kmb> [--structural] [max pairs]
"""
import dataclasses
import heapq
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402

args = [a for a in sys.argv[1:] if not a.startswith('--')]
STRUCT = '--structural' in sys.argv
src, pairs_npz, zero_npz, dst = args[:4]
maxp = int(args[4]) if len(args) > 4 else None
c = kmb.read(src)
ops = np.asarray(c.ops)
n = len(ops)
nq = (ops & 0x30) >> 4
W = np.full((n, 3), -1, dtype=np.int64)
for k, arrs in ((3, c.q3), (2, c.q2), (1, (c.q1,))):
    sel = np.flatnonzero(nq == k)
    for j in range(k):
        W[sel, j] = np.asarray(arrs[j], dtype=np.int64)
Z = np.load(zero_npz)
zw, zs, ze = Z['wire'], Z['start'].astype(np.int64), Z['end'].astype(np.int64)
if STRUCT:
    keep = Z['structural']
    zw, zs, ze = zw[keep], zs[keep], ze[keep]
# a wire with idle interval (s, e) is free at "after-op" times s .. e-1 (touched at ops s and e)
lo_t = zs
hi_t = ze - 1
ok = hi_t >= lo_t
zw, lo_t, hi_t = zw[ok], lo_t[ok], hi_t[ok]
# compressed timeline: breakpoints where the free count changes, plus pair endpoints
P = np.load(pairs_npz)
pf, ps = P['pair_first'].astype(np.int64), P['pair_second'].astype(np.int64)
# a pair needs a free wire at after-op times i-1 .. j
plo, phi = pf - 1, ps
pts = np.unique(np.concatenate([lo_t, hi_t + 1, plo, phi + 1, [-1, n]]))
seg_of = lambda t: np.searchsorted(pts, t, side='right') - 1        # noqa: E731
M = len(pts)
diff = np.zeros(M + 1, dtype=np.int64)
np.add.at(diff, seg_of(lo_t), 1)
np.add.at(diff, seg_of(hi_t + 1), -1)
free = np.cumsum(diff)[:M]
print('timeline segments %d; free wires: min %d, median %d, max %d' % (M, free.min(), int(np.median(free)),
      free.max()), flush=True)


class SegTree:
    """range add, range min over M leaves."""
    def __init__(self, vals):
        self.n = 1
        while self.n < len(vals):
            self.n *= 2
        self.mn = np.full(2 * self.n, 1 << 40, dtype=np.int64)
        self.lz = np.zeros(2 * self.n, dtype=np.int64)
        self.mn[self.n:self.n + len(vals)] = vals
        for i in range(self.n - 1, 0, -1):
            self.mn[i] = min(self.mn[2 * i], self.mn[2 * i + 1])

    def _q(self, node, nl, nr, l, r):
        if r < nl or nr < l:
            return 1 << 40
        if l <= nl and nr <= r:
            return self.mn[node]
        mid = (nl + nr) // 2
        return self.lz[node] + min(self._q(2 * node, nl, mid, l, r), self._q(2 * node + 1, mid + 1, nr, l, r))

    def _a(self, node, nl, nr, l, r, v):
        if r < nl or nr < l:
            return
        if l <= nl and nr <= r:
            self.mn[node] += v
            self.lz[node] += v
            return
        mid = (nl + nr) // 2
        self._a(2 * node, nl, mid, l, r, v)
        self._a(2 * node + 1, mid + 1, nr, l, r, v)
        self.mn[node] = self.lz[node] + min(self.mn[2 * node], self.mn[2 * node + 1])

    def query(self, l, r):
        return self._q(1, 0, self.n - 1, l, r)

    def add(self, l, r, v):
        self._a(1, 0, self.n - 1, l, r, v)


st = SegTree(free)
order = np.argsort(ps - pf, kind='stable')
used_ops = set()
chosen = []
for p in order:
    i, j = int(pf[p]), int(ps[p])
    if i in used_ops or j in used_ops:
        continue
    a, b = int(seg_of(np.array([i - 1]))[0]), int(seg_of(np.array([j]))[0])
    if st.query(a, b) >= 1:
        st.add(a, b, -1)
        chosen.append((i, j))
        used_ops.add(i)
        used_ops.add(j)
        if maxp is not None and len(chosen) >= maxp:
            break
chosen.sort()
print('pairs offered %d, admitted %d' % (len(pf), len(chosen)), flush=True)

# sweep: assign wires, relocate a copy before the original touches its wire (one pass)
starts = {}
for w, l, h in zip(zw.tolist(), lo_t.tolist(), hi_t.tolist()):
    starts.setdefault(l, []).append((w, h))
pstart, pend = {}, {}
for idx, (i, j) in enumerate(chosen):
    pstart.setdefault(i - 1, []).append(idx)
    pend.setdefault(j, []).append(idx)
free_heap = []                     # max-heap on the interval end: (-until, wire)
holding = {}                       # wire -> (pair idx, until)
pwire, start_wire, end_wire, hops = {}, {}, {}, {}


def take_free(t):
    while free_heap:
        nu, w = heapq.heappop(free_heap)
        if -nu >= t and w not in holding:
            return w, -nu
    raise RuntimeError('no free wire at %d' % t)


queue = sorted(set(starts) | set(pstart) | set(pend))
heapq.heapify(queue)
seen = set()
while queue:
    t = heapq.heappop(queue)
    if t in seen:
        continue
    seen.add(t)
    for w, h in starts.get(t, []):
        heapq.heappush(free_heap, (-h, w))
    for idx in pend.get(t, []):
        w = pwire.pop(idx)
        end_wire[idx] = w
        until = holding.pop(w)[1]
        if until > t:
            heapq.heappush(free_heap, (-until, w))
    for idx in pstart.get(t, []):
        w, until = take_free(t + 1)
        holding[w] = (idx, until)
        pwire[idx] = w
        start_wire[idx] = w
        heapq.heappush(queue, until)
    for w in [w for w, (idx, until) in holding.items() if until == t]:
        idx, _ = holding.pop(w)
        w2, until2 = take_free(t + 1)
        holding[w2] = (idx, until2)
        pwire[idx] = w2
        hops.setdefault(t, []).append((w, w2))
        heapq.heappush(queue, until2)
print('hops (SWAPs inserted): %d' % sum(len(v) for v in hops.values()), flush=True)
nb = c.num_bits
inserts = {}
delete = set()
for t, lst in hops.items():
    inserts.setdefault(t, []).extend((0x22, (a, b), 0, 0) for a, b in lst)
for idx, (i, j) in enumerate(chosen):
    a, b, t = (int(x) for x in W[i])
    a2, b2, t2 = (int(x) for x in W[j])
    s0, s1 = start_wire[idx], end_wire[idx]
    m = nb + idx
    inserts.setdefault(i, [])
    inserts[i] = [(0x31, (a, b, s0), 0, 0), (0x21, (s0, t), 0, 0)] + inserts[i]
    inserts.setdefault(j, [])
    inserts[j] = [(0x21, (s1, t2), 0, 0), (0x50, (s1,), m, 0), (0xA0, (a2, b2), 0, m)] + inserts[j]
    delete.add(i)
    delete.add(j)
out = kw.edit(c, delete, inserts)
out = dataclasses.replace(out, num_bits=nb + len(chosen))
nqub = kw.write(dst, out)
print('qubits %d, bits %d -> %d, static Toffoli-like %d -> %d' % (
    nqub, nb, out.num_bits, int(((ops & 0x7C) == 0x30).sum()), int(((out.ops & 0x7C) == 0x30).sum())))
np.save(dst + '.pairs.npy', np.array(chosen, dtype=np.int64))
