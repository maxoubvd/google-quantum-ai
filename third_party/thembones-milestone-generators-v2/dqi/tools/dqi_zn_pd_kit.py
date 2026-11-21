#!/usr/bin/env python3
"""dqi_zn_pd_kit.py -- builder Z1: the kit's own bez step (exp/safe_scalar/converter.py bez_lookahead, forward)
built in kickmix on registers (O 17, E 18, u 17, w 17), run on the same charts as the word model
tools/dqi_zn_pd.py, outputs compared word by word. Validates the K-model (gap routes, clips) against the kit.
Usage: python dqi_zn_pd_kit.py [kinds=0,1,2,3,4] [seed=1]
"""
import os
import sys

import numpy as np

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
KIT = os.path.expanduser('~/qcc/ws/dqi_sw/inputs/source_cap902/playbook/kit/decoder')
sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
import kickmix as km  # noqa: E402
sys.modules['circ'] = km
sys.path.insert(0, KIT)
sys.path.insert(0, KIT + '/exp/convert')
import dqi_int_convtest as CT  # noqa: E402  (charts, bit packing; B5, read-only)
import dqi_zn_pd as PD  # noqa: E402
from exp.safe_scalar import converter as SS  # noqa: E402


def build():
    b = km.CircuitBuilder()
    regs = dict(O=[b.alloc_qubits(12) for _ in range(17)], E=[b.alloc_qubits(12) for _ in range(18)],
                u=[b.alloc_qubits(12) for _ in range(17)], w=[b.alloc_qubits(12) for _ in range(17)])
    phase = {'a': {}, 'q': {}}
    SS.bez_lookahead(b, regs['O'], regs['E'], [regs['u'], regs['w']], False, None, phase)
    return b.finish_circuit(), regs


def run(c, regs, vals, seed):
    n = vals['O'].shape[0]
    sim = km.Simulator(batch_size=n, debug_print=False, counting_operations=True, seed=seed)
    sim.use_same_registers_as(c)
    qs = {}
    for k in ('O', 'E', 'u', 'w'):
        qs[k] = km.array([q for w in regs[k] for q in w])
        sim.write_across_shots(qs[k], CT.to_bits(vals[k]))
    sim.do(c)
    out = {k: CT.from_bits(sim.read_across_shots(qs[k]), len(regs[k])) for k in qs}
    return out, sim.operation_statistics().num_toffoli_like / n


def main():
    kinds = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else '0,1,2,3,4').split(',')]
    rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    c, regs = build()
    print('kit bez forward: qubits %d, static Toffoli-like %d' % (c.num_qubits, c.operation_statistics().num_toffoli_like),
          flush=True)
    for kd in kinds:
        A, B, label = CT.charts(kd, rng)
        O, E, u, w = PD.bez_inputs(A, B)
        vals = dict(O=O, E=E, u=u, w=w)
        ko, tof = run(c, regs, vals, 7)
        mo, cnt, _ = PD.run_bez(PD.Kit, O, E, u, w)
        eq = np.ones(O.shape[0], dtype=bool)
        for k in ('O', 'E', 'u', 'w'):
            eq &= (ko[k] == mo[k]).all(axis=1)
        print('%-22s model == kit on %d / %d shots; kit executed %.1f' % (label, eq.sum(), len(eq), tof), flush=True)


if __name__ == '__main__':
    main()
