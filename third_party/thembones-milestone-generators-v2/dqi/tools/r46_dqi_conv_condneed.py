#!/usr/bin/env python3
"""r46_dqi_conv_condneed.py -- R46-DQI-CONV: how many wires does a classically conditioned block allocate?

tools/r46_dqi_conv_zpark.py cannot evict a parked flag inside a PUSH_CONDITION block, so it evicts every park at the
block's entry when fewer than R46_ZP_MARGIN wires are free (default 40: 23 evictions in build zpsu1). This reads, for
every outermost conditioned block of convert / unconvert in the raw file, the allocation at its entry and the largest
allocation inside (allocation trace of tools/r46_dqi_conv_ztrace.py, same static sequence), and prints what each
margin would evict and whether it is safe (need <= margin for every block entered above CAP - margin).
  python r46_dqi_conv_condneed.py <raw kmb> <census npz> <trace npz> [CAP=902] [k ranges=90368-145074,189360-258566]
"""
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import kmb  # noqa: E402

c = kmb.read(sys.argv[1])
z = np.load(sys.argv[2])
T = np.load(sys.argv[3])
CAP = int(sys.argv[4]) if len(sys.argv) > 4 else 902
rng = [tuple(int(x) for x in r.split('-')) for r in (sys.argv[5] if len(sys.argv) > 5 else '90368-145074,189360-258566').split(',')]
ops = np.asarray(c.ops)
t_op = z['t_op']
nt = len(t_op)
i = [j for j in range(20) if 'total_%d' % j in T and int(T['total_%d' % j]) == nt]
K, A = T['K_%d' % i[-1]].astype(np.int64), T['A_%d' % i[-1]].astype(np.int64)
mx = np.zeros(nt + 2, np.int64)
np.maximum.at(mx, K, A)
last = np.zeros(nt + 2, np.int64)
last[K] = A
has = np.zeros(nt + 2, bool)
has[K] = True
lev = np.zeros(nt + 2, np.int64)        # level after the last event at or before k
cur = 0
for k in range(nt + 2):
    if has[k]:
        m = max(mx[k], cur)
        cur = last[k]
        mx[k] = m
    else:
        mx[k] = cur
    lev[k] = cur
push = np.flatnonzero(ops == 0x81)
pop = np.flatnonzero(ops == 0x01)
ev = np.concatenate([push, pop])
kind = np.concatenate([np.ones(len(push), np.int64), -np.ones(len(pop), np.int64)])
o = np.argsort(ev)
ev, kind = ev[o], kind[o]
depth = np.cumsum(kind)
start = ev[(kind == 1) & (depth == 1)]
end = ev[(kind == -1) & (depth == 0)]
assert len(start) == len(end)
k0 = np.searchsorted(t_op, start)
k1 = np.searchsorted(t_op, end)
rows = []
for a, b in zip(k0, k1):
    if not any(lo <= a < hi for lo, hi in rng):
        continue
    entry = int(lev[a - 1]) if a > 0 else 0
    inside = int(mx[a:b + 1].max()) if b > a else entry
    rows.append((entry, max(0, inside - entry), b - a))
rows = np.array(rows)
print('outermost conditioned blocks in the ranges: %d; with Toffolis inside: %d; largest need %d wires'
      % (len(rows), int((rows[:, 2] > 0).sum()), int(rows[:, 1].max())))
for mrg in (0, 4, 8, 12, 16, 20, 24, 32, 40):
    hit = rows[:, 0] + mrg > CAP
    unsafe = int(((rows[:, 0] + rows[:, 1] > CAP - 1) & ~hit).sum())
    print('margin %2d: blocks entered with fewer free wires than the margin %4d; blocks NOT pre-evicted whose inside '
          'peak leaves no wire for one park %d' % (mrg, int(hit.sum()), unsafe))
need = rows[:, 1]
print('need histogram (wires: blocks): ' + ', '.join('%d: %d' % (n, int((need == n).sum())) for n in sorted(set(need.tolist()))))
print('entry level of the blocks that need wires (need, entry): ' + ', '.join('%d@%d' % (n, e) for e, n, _ in rows[rows[:, 1] > 0][:60]))
