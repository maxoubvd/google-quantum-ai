#!/usr/bin/env python3
"""r48_dqi_ret2_occ.py -- seat R48-DQI-RET2 (fourth instance): occupancy of the top wires, stage by stage. Read only.

For a decoder file and its exact zero intervals (tools/dqi_pp_zero.py; live(t) = Q - wires exactly zero after op t,
as tools/dqi_int_budget.py) prints per top-level stage: the tight peak, all ops and Toffoli-class ops at
live >= each threshold, and the number of maximal runs of consecutive ops at live >= the FIRST threshold (the
"runs at 902" of R46-DQI-SCAN's 901 list).

Usage:
  python r48_dqi_ret2_occ.py <file.kmb> <zero.npz> [thresholds=902,901,900,896,890]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
from dqi_pp_gvn import stages  # noqa: E402


def main():
    path, zp = sys.argv[1:3]
    th = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else '902,901,900,896,890').split(',')]
    c = kmb.read(path)
    ops = np.asarray(c.ops)
    n = len(ops)
    Z = np.load(zp)
    s = np.clip(Z['start'].astype(np.int64), 0, n)
    e = np.clip(Z['end'].astype(np.int64), 0, n)
    d = np.zeros(n + 2, np.int32)
    ok = e > s
    np.add.at(d, s[ok], 1)
    np.add.at(d, e[ok], -1)
    live = c.num_qubits - np.cumsum(d[:n])
    toff = ((ops & 0x7C) == 0x30) & (ops != 0x81)
    print('%s: Q %d, ops %d, Toffoli-class %d' % (path, c.num_qubits, n, int(toff.sum())))
    print('%-20s %5s %9s ' % ('stage', 'peak', 'runs>=%d' % th[0])
          + ' '.join('%10s %8s' % ('ops>=%d' % t, 'T>=%d' % t) for t in th))
    for nm, a, b in stages(c):
        if nm in ('mask_projection', 'initial'):
            continue
        L = live[a:b]
        T = toff[a:b]
        hi = (L >= th[0]).astype(np.int8)
        runs = int((np.diff(np.concatenate(([0], hi))) == 1).sum())
        print('%-20s %5d %9d ' % (nm, int(L.max()), runs)
              + ' '.join('%10d %8d' % (int((L >= t).sum()), int((T & (L >= t)).sum())) for t in th))


if __name__ == '__main__':
    main()
