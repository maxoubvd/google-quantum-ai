#!/usr/bin/env python3
"""dqi_pp_pprove.py -- workstream dqi_pp (P1): prove retained-product links exactly (two-time claim).

Links: consecutive occurrences (i, j) of one product fingerprint (tools/dqi_pp_pairs2.py, 192 runs) with gap < G.
Each is checked by the local exhaustive prover (tools/dqi_pp_prove.py class 7: a_i * b_i at op i equals a_j * b_j at
op j on every state of the window). Links with the same exact value-numbering form are exact without the prover.
Writes the exact links (i, j) for tools/dqi_pp_ret.py (DQI_PP_LINKS).
Usage: python dqi_pp_pprove.py <circuit.kmb> <pairs2.npz> <gvn.npz> <out.npz> [G=200000] [K=20]
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
    Pp = np.load(sys.argv[2])
    G = np.load(sys.argv[3])
    out = sys.argv[4]
    gmax = int(sys.argv[5]) if len(sys.argv) > 5 else 200000
    K = int(sys.argv[6]) if len(sys.argv) > 6 else 20
    fo, key = Pp['fo'], Pp['key']
    order = np.lexsort((fo, key))
    same = key[order][1:] == key[order][:-1]
    li, lj = fo[order[:-1][same]], fo[order[1:][same]]
    t = G['t_op']
    a, b = np.searchsorted(t, li), np.searchsorted(t, lj)
    exact = (G['t_p0'][a] == G['t_p0'][b]) & (G['t_p1'][a] == G['t_p1'][b]) & (G['t_pc'][a] == G['t_pc'][b])
    short = (lj - li) < gmax
    cand = short & ~exact
    print('links %d (gap < %d: %d), value-numbering exact %d, to prove %d' % (
        len(li), gmax, int(short.sum()), int(exact.sum()), int(cand.sum())), flush=True)
    ops, W0, W1, W2, BT, CB, DEP = P.decode(c)
    ci, cj = li[cand].astype(np.int64), lj[cand].astype(np.int64)
    res, nv = P.prove_all(ops, W0, W1, W2, BT, CB, DEP, cj, np.full(len(cj), 7, np.int64), ci, K, max(gmax + 1000,
                          400000), c.num_qubits, c.num_bits)
    print('proved %d, refuted %d, cone too large %d (%.1fs)' % (int((res == 1).sum()), int((res == 0).sum()),
          int((res == -1).sum()), time.time() - t0), flush=True)
    ok_i = np.r_[li[exact], ci[res == 1]]
    ok_j = np.r_[lj[exact], cj[res == 1]]
    np.savez(out, i=ok_i, j=ok_j, refuted_i=ci[res == 0], refuted_j=cj[res == 0])
    print('exact links written: %d' % len(ok_i))


if __name__ == '__main__':
    main()
