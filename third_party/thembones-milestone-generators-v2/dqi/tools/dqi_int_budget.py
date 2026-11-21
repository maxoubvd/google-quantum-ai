#!/usr/bin/env python3
"""dqi_int_budget.py -- workstream dqi_int (B5): THE QUBIT BUDGET trace (lead 00:19Z): tight live-qubit peaks of a DQI
.kmb per top-level stage and per innermost '@enter:' mark inside each stage.

live(t) = Q - free(t), free(t) = wires exactly zero after op t (exact zero intervals of tools/dqi_pp_zero.py, as in
tools/dqi_pp_capprof.py): allocated-but-zero wires (placeholders) count as free, so this is the TIGHT peak.
Per region: tight peak, and the number of ops with live above 890, 878 and 866 (how long the region sits there).
Usage: python dqi_int_budget.py <raw.kmb> <zero.npz> [min ops per mark row=1]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
from dqi_pp_gvn import stages  # noqa: E402


def main():
    c = kmb.read(sys.argv[1])
    n = len(c.ops)
    Z = np.load(sys.argv[2])
    minops = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    s = np.clip(Z['start'].astype(np.int64), 0, n)
    e = np.clip(Z['end'].astype(np.int64), 0, n)
    d = np.zeros(n + 2, np.int32)
    ok = e > s
    np.add.at(d, s[ok], 1)
    np.add.at(d, e[ok], -1)
    live = c.num_qubits - np.cumsum(d[:n])
    tp = c.tag_positions()
    tags = c.tags
    print('Q %d; tight peak of the whole file %d' % (c.num_qubits, int(live.max())))
    print('%-34s %6s %10s %10s %10s %10s' % ('region', 'peak', 'ops>890', 'ops>878', 'ops>866', 'ops'))
    for nm, a, b in stages(c):
        L = live[a:b]
        if len(L) == 0:
            continue
        print('%-34s %6d %10d %10d %10d %10d' % (nm, int(L.max()), int((L > 890).sum()), int((L > 878).sum()),
                                                 int((L > 866).sum()), len(L)))
        # innermost marks inside the stage (first level below the stage name)
        rows = {}
        stack = []
        last = a
        k0 = int(np.searchsorted(tp, a))
        k1 = int(np.searchsorted(tp, b))
        for k in range(k0, k1 + 1):
            p = int(tp[k]) if k < len(tp) else b
            p = min(p, b)
            key = stack[0] if stack else '(stage body)'
            seg = live[last:p]
            if len(seg):
                r = rows.setdefault(key, [0, 0, 0, 0, 0])
                r[0] = max(r[0], int(seg.max()))
                r[1] += int((seg > 890).sum())
                r[2] += int((seg > 878).sum())
                r[3] += int((seg > 866).sum())
                r[4] += len(seg)
            last = p
            if k >= len(tp) or p >= b:
                break
            t = tags[k]
            if t.startswith('@enter:') and t[7:] != nm:
                stack.append(t[7:])
            elif t.startswith('@exit:') and stack and t[6:] == stack[-1]:
                stack.pop()
        for key, r in sorted(rows.items(), key=lambda kv: -kv[1][0]):
            if r[4] >= minops and r[0] > 866:
                print('   %-31s %6d %10d %10d %10d %10d' % (key[:31], r[0], r[1], r[2], r[3], r[4]))


if __name__ == '__main__':
    main()
