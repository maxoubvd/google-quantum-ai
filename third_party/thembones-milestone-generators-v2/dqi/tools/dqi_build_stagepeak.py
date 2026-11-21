#!/usr/bin/env python
"""Highest qubit index touched in each top-level stage of a .kmb file (workstream dqi_build).

The builder hands out the lowest free indices, so (highest index + 1) is the stage's peak width as the cost sees it.
Usage: python dqi_build_stagepeak.py <file.kmb> [more.kmb ...]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb

for path in sys.argv[1:]:
    c = kmb.read(path)
    ops = np.asarray(c.ops)
    k3 = np.concatenate([[0], np.cumsum((ops & 0x30) == 0x30)])
    k2 = np.concatenate([[0], np.cumsum((ops & 0x30) == 0x20)])
    k1 = np.concatenate([[0], np.cumsum((ops & 0x30) == 0x10)])
    q3 = [np.asarray(a) & kmb.ID_MASK for a in c.q3]
    q2 = [np.asarray(a) & kmb.ID_MASK for a in c.q2]
    q1 = np.asarray(c.q1) & kmb.ID_MASK
    depth, start, name, out = 0, None, None, []
    for p, t in zip(c.tag_positions().tolist(), c.tags):
        if t.startswith('@enter:'):
            if depth == 0:
                start, name = p, t[7:]
            depth += 1
        elif t.startswith('@exit:'):
            depth -= 1
            if depth == 0:
                m = 0
                for arrs, k in ((q3, k3), (q2, k2), ([q1], k1)):
                    for a in arrs:
                        seg = a[k[start]:k[p]]
                        if len(seg):
                            m = max(m, int(seg.max()) + 1)
                out.append(f'{name} {m}')
    print(os.path.basename(os.path.dirname(path)) + ': ' + '  '.join(out))
