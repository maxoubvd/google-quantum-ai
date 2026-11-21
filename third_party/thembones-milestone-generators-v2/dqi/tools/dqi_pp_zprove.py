#!/usr/bin/env python3
"""dqi_pp_zprove.py -- workstream dqi_pp (P1): prove sampled idle-zero intervals exactly, merge with the exact set.

Candidates: intervals of tools/dqi_int_ppzero.py (zero on every sample) that tools/dqi_pp_zero.py did not prove.
Each is checked by the local exhaustive prover (tools/dqi_pp_prove.py, claim class 6: the wire is zero right after
its touch at `start`). Proved intervals join the exact ones; the rest are dropped (they may be rare flags).
Usage: python dqi_pp_zprove.py <circuit.kmb> <sampled zero.npz> <exact zero.npz> <out.npz> [K=18] [L=400000]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_pp_prove as P  # noqa: E402


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    Zs = np.load(sys.argv[2])
    Ze = np.load(sys.argv[3])
    out = sys.argv[4]
    K = int(sys.argv[5]) if len(sys.argv) > 5 else 18
    L = int(sys.argv[6]) if len(sys.argv) > 6 else 400000
    ke = set(zip(Ze['wire'].tolist(), Ze['start'].tolist()))
    w, s, e = Zs['wire'].astype(np.int64), Zs['start'].astype(np.int64), Zs['end'].astype(np.int64)
    cand = np.array([k for k in range(len(w)) if s[k] >= 0 and e[k] > s[k] + 1 and (int(w[k]), int(s[k])) not in ke],
                    dtype=np.int64)
    print('sampled intervals %d, exact %d, candidates to prove %d' % (len(w), len(Ze['wire']), len(cand)), flush=True)
    ops, W0, W1, W2, BT, CB, DEP = P.decode(c)
    res, nvs = P.prove_all(ops, W0, W1, W2, BT, CB, DEP, s[cand] + 1, np.full(len(cand), 6, np.int64), w[cand], K, L,
                           c.num_qubits, c.num_bits)
    ok = res == 1
    print('proved %d, not proved %d, cone too large %d (%.1fs)' % (int(ok.sum()), int((res == 0).sum()),
          int((res == -1).sum()), time.time() - t0), flush=True)
    pw, ps_, pe = w[cand][ok], s[cand][ok], e[cand][ok]
    wire = np.concatenate([Ze['wire'], pw])
    start = np.concatenate([Ze['start'], ps_])
    end = np.concatenate([Ze['end'], pe])
    structural = np.concatenate([Ze['structural'], np.zeros(len(pw), bool)])
    np.savez(out, wire=wire, start=start, end=end, structural=structural)
    print('merged intervals %d; wire-ops added %d' % (len(wire), int((pe - ps_).sum())))


if __name__ == '__main__':
    main()
