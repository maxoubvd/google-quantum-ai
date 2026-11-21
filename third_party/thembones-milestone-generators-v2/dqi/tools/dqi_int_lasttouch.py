#!/usr/bin/env python3
"""dqi_int_lasttouch.py -- B5 (17th holder, Nov 16): read only.  The ops that touch wire w in a window [a, b) of a .kmb
(op index, op code, operands), to place a park spec's 'hmr' line (the last touch of a wire before a run).
Usage: python dqi_int_lasttouch.py <dqi.kmb> <wire> <a> <b>
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402


def main():
    path, w, a, b = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    k = kmb.read(path)
    ops = np.asarray(k.ops)
    kk = ((ops & 0x30) >> 4).astype(np.int8)
    i1 = np.cumsum(kk == 1) - 1
    i2 = np.cumsum(kk == 2) - 1
    i3 = np.cumsum(kk == 3) - 1
    M = kmb.ID_MASK
    q1 = np.asarray(k.q1)
    q2 = [np.asarray(x) for x in k.q2]
    q3 = [np.asarray(x) for x in k.q3]
    hits = []
    for t in range(a, b):
        c = int(kk[t])
        if c == 3:
            o = [int(q3[j][i3[t]]) & M for j in range(3)]
        elif c == 2:
            o = [int(q2[j][i2[t]]) & M for j in range(2)]
        elif c == 1:
            o = [int(q1[i1[t]]) & M]
        else:
            o = []
        if w in o:
            hits.append((t, hex(int(ops[t])), o))
    print('%s: wire %d touched %d times in [%d, %d)' % (path, w, len(hits), a, b))
    for h in hits[-12:]:
        print('  op %d code %s operands %s' % h)


if __name__ == '__main__':
    main()
