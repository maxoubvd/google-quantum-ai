#!/usr/bin/env python3
"""dqi_pp_capprof.py -- workstream dqi_pp (P1): the free-wire profile of a DQI .kmb from its EXACT zero intervals.

free(t) = number of wires that are exactly zero after op t (intervals of tools/dqi_pp_zero.py: free at after-op times
start .. end-1). Prints, per stage, the minimum and the share of ops at each low count, and for one stage the minimum
in windows of W ops (where an extra word of workspace would fit).
Usage: python dqi_pp_capprof.py <circuit.kmb> <zero.npz> [stage=padded_unsynthesis] [W=200000]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
from dqi_pp_gvn import stages  # noqa: E402


def main():
    c = kmb.read(sys.argv[1])
    Z = np.load(sys.argv[2])
    want = sys.argv[3] if len(sys.argv) > 3 else 'padded_unsynthesis'
    W = int(sys.argv[4]) if len(sys.argv) > 4 else 200000
    n = len(c.ops)
    d = np.zeros(n + 2, np.int32)
    s = np.clip(Z['start'].astype(np.int64), 0, n)
    e = np.clip(Z['end'].astype(np.int64), 0, n)
    ok = e > s
    np.add.at(d, s[ok], 1)
    np.add.at(d, e[ok], -1)
    free = np.cumsum(d[:n])
    print('ops %d; free wires: min %d, median %d' % (n, int(free.min()), int(np.median(free))))
    for nm, a, b in stages(c):
        f = free[a:b]
        if len(f) == 0:
            continue
        print('%-20s min %3d  ops at 0: %8d  <=5: %8d  <=11: %8d  <=23: %8d of %d' % (
            nm, int(f.min()), int((f == 0).sum()), int((f <= 5).sum()), int((f <= 11).sum()), int((f <= 23).sum()),
            len(f)))
        if nm == want:
            for lo in range(a, b, W):
                ff = free[lo:min(lo + W, b)]
                print('   ops %d-%d: min free %d, p10 %d' % (lo, min(lo + W, b), int(ff.min()),
                                                           int(np.percentile(ff, 10))))


if __name__ == '__main__':
    main()
