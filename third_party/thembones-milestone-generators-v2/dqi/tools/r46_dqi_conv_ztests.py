#!/usr/bin/env python3
"""r46_dqi_conv_ztests.py -- R46-DQI-CONV: every 12-bit zero test (zero_toggle, 11 ANDs) of convert / unconvert:
which word (wires read from the census file's control columns), which call site, and how often the SAME wires are
tested again (a predicate measured out and made again).
  python r46_dqi_conv_ztests.py <profdir> <one census npz> [stages=convert,unconvert] [list=0]
"""
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import r46_dqi_conv_report as R  # noqa: E402

pdir, f = sys.argv[1], sys.argv[2]
stages = sys.argv[3].split(',') if len(sys.argv) > 3 else ['convert', 'unconvert']
full = len(sys.argv) > 4 and sys.argv[4] == '1'
marks, stacks, seg, ex, shots = R.load(pdir)
z = np.load(f)
wa, wb, act, cb, n = z['wa'], z['wb'], z['act'], z['cb'], int(z['samples'])
ones = z['ones']
k0 = 0
tests = []
for i, (s, mk, c) in enumerate(seg):
    nm = R.names(stacks[s])
    if c and nm[-1] == 'zero_toggle' and marks[mk] and marks[mk][0] in stages:
        for j in range(0, c, 11):
            ks = np.arange(k0 + j, k0 + min(j + 11, c))
            wires = tuple([int(wa[ks[0]])] + [int(x) for x in wb[ks]])
            depth = sum(1 for x in nm if x in R.REC)
            # rate of the predicate = rate of the last AND's product is not stored; use cb of last (ladder) op
            tests.append((marks[mk][0], k0 + j, wires, depth, '>'.join(nm[-5:-1]), len(ks), int(ones[ks[-1]])))
    k0 += c
if len(sys.argv) > 5:                      # save (index of the last AND of each test, slot id) for r46_dqi_conv_zsame.py
    slots = {}
    rows = []
    for t in tests:
        sid = slots.setdefault(tuple(sorted(t[2])), len(slots))
        rows.append((t[1] + t[5] - 1, sid, min(t[2])))
    np.save(sys.argv[5], np.array(rows, dtype=np.int64))
    print('saved %d tests, %d slots -> %s' % (len(rows), len(slots), sys.argv[5]))
for st in stages:
    T = [t for t in tests if t[0] == st]
    byw = defaultdict(list)
    for t in T:
        byw[tuple(sorted(t[2]))].append(t)
    rep = sum(len(v) - 1 for v in byw.values())
    print('== %s: %d zero tests (%d ANDs), %d distinct wire sets, %d repeats on the same wires (%d ANDs)'
          % (st, len(T), sum(t[5] for t in T), len(byw), rep, 11 * rep))
    for w, v in sorted(byw.items(), key=lambda kv: kv[1][0][1]):
        if len(v) > 1 or full:
            print('  wires q%d..q%d x%d (hits per test: %s)' % (min(w), max(w), len(v), ','.join(str(t[6]) for t in v)))
            if full:
                for t in v:
                    print('      k%-7d d%-2d hits %-5d %s' % (t[1], t[3], t[6], t[4]))
