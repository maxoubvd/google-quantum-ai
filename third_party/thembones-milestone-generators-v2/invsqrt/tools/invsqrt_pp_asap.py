#!/usr/bin/env python3
"""invsqrt_pp_asap.py -- workstream invsqrt_pp (P1): the measure-early pass written as a .kmb, with an all-input check.

Every top-level unconditioned HMR moves up to right after the last op that touches its wire or mentions its bit
(after the POP of the block holding that op). 0 Toffolis; the moved measurement acts on the same value and its kick
commutes with every op it crosses. HMRs that could go before op 0 stay where they are.
Check: the child, simulated on all 3,133,452 inputs (zero transcript and 3 random transcripts), ends in the parent's
final state with the parent's Z-phase on every input, and every Z_POW sees the parent's wire value on every input.
Usage: python invsqrt_pp_asap.py <parent.kmb> <out.kmb>
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402
import invsqrt_pp_sim as sim  # noqa: E402
from invsqrt_pp_links import asap_order  # noqa: E402


def zpow_values(C, outcome=None):
    """Hash of the wire value at every Z_POW (in op order)."""
    out = []
    r = np.random.default_rng(1).integers(1, 2 ** 62, size=sim.W, dtype=np.int64).astype(np.uint64)

    def hook(i, S, P, bits, act):
        if C.base[i] == 0x13 and act:
            out.append(int(np.bitwise_xor.reduce(S[C.qs[i][0]] * r)))
    sim.run(C, hook=hook, outcome=outcome)
    return out


def main():
    t0 = time.time()
    src, dst = sys.argv[1], sys.argv[2]
    C = sim.load(src)
    order, nmoved = asap_order(C)
    pos = {op: k for k, op in enumerate(order)}
    # every moved HMR h follows the op just before it in the new order; express as delete + insert after anchor
    moved = [h for h in range(C.n) if C.base[h] == 0x50 and C.level[h] == 0 and C.cnd[h] < 0
             and pos[h] > 0 and order[pos[h] - 1] != h - 1]
    k = kmb.read(src)
    ops = np.asarray(k.ops)
    delete, inserts = [], {}
    for h in moved:
        anchor = order[pos[h] - 1]
        delete.append(h)
        inserts.setdefault(anchor, []).append((int(ops[h]), C.qs[h], C.bit[h], 0))
    # keep the relative order of HMRs inserted after the same anchor as in the new order
    k2 = kw.edit(k, delete, inserts)
    kw.write(dst, k2)
    D = sim.load(dst)
    print('HMRs moved', len(moved), 'ops', C.n, '->', D.n, flush=True)
    S0, P0, _, _ = sim.run(C)
    ok = True
    for seed in (0, 1, 2, 3):
        oc = None
        if seed:
            dr = np.random.default_rng(seed).integers(0, 2, size=D.n)
            oc = (lambda j, d=dr: int(d[j]))
        S1, P1, _, _ = sim.run(D, outcome=oc)
        e = bool(np.array_equal(S1, S0)) and bool(np.array_equal(P1, P0))
        print('  seed', seed, 'final state and Z-phase equal to the parent:', e, flush=True)
        ok &= e
    zp = zpow_values(C) == zpow_values(D)
    print('  every Z_POW sees the parent value:', zp)
    print('CHECK', 'PASS' if ok and zp else 'FAIL', round(time.time() - t0, 1), 's')


if __name__ == '__main__':
    main()
