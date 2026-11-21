#!/usr/bin/env python3
"""dqi_pp_capshift.py -- workstream dqi_pp (P1): how many of the exact pass's links survive when every moment has d
fewer free wires (a reading of the pass at Q = 902 - d on today's layout: free(t) -> max(0, free(t) - d)).

Same link set, zero intervals and admission as tools/dqi_pp_ret.py in exact mode (built by dqi_pp_windows.build).
Also per stage: links admitted at d = 0 by the stage of their second occurrence.
Usage: python dqi_pp_capshift.py <raw.kmb> <gvn.npz> <zero.npz> [d list, default 0,1,2,4,6,8,12]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
from dqi_pp_windows import build, run_admit  # noqa: E402
from dqi_pp_gvn import stages  # noqa: E402


def main():
    t0 = time.time()
    raw, gvn, zero = sys.argv[1:4]
    ds = [int(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else [0, 1, 2, 4, 6, 8, 12]
    c = kmb.read(raw)
    ops = np.asarray(c.ops)
    n = len(ops)
    B = build(c, np.load(gvn), np.load(zero), ops, n)
    st = stages(c)
    for d in ds:
        adm, _, _, _ = run_admit(B, np.maximum(B['free'] - d, 0))
        line = 'd = %2d (Q %d): links admitted %d of %d' % (d, c.num_qubits - d, int(adm.sum()), len(adm))
        if d == 0:
            by = {}
            for L in np.flatnonzero(adm).tolist():
                j = int(B['lk_hi'][L])
                name = next((nm for nm, a, b in st if a <= j <= b), '-')
                by[name] = by.get(name, 0) + 1
            line += '; by stage of the second occurrence: %s' % by
        print(line, '(%.1fs)' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
