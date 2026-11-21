#!/usr/bin/env python
"""R7: executed Toffolis per profile segment (workstream dqi_build).

Runs a circuit built with the patch `r7prof` on verifier-distributed shots. At every DEBUG_PRINT marker the
simulator's executed Toffoli counter is read; the increase since the previous marker belongs to the segment that
the previous marker opened. The simulator runs the circuit once per chunk of 256 shots, so the marker sequence
repeats per chunk.
Usage: python dqi_build_r7_prof_run.py <probe.kmb> <segments.json> <out.json> [batches=8] [batch size=1024]
                                       [workers=3] [seed=1]
Output JSON: {shots, errors, executed: [sum over shots per segment], total}.
"""
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

DATA = os.path.expanduser('~/qcc/data')
sys.path[:0] = [os.path.join(DATA, 'kickmix'), os.path.join(DATA, 'verifiers')]

path, sfile, out = sys.argv[1:4]
batches = int(sys.argv[4]) if len(sys.argv) > 4 else 8
bsize = int(sys.argv[5]) if len(sys.argv) > 5 else 1024
workers = int(sys.argv[6]) if len(sys.argv) > 6 else 3
seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
NS = len(json.load(open(sfile))['segments'])
_ctx = {}


def setup():
    if _ctx:
        return
    import kickmix as km
    from dqi import verify_submission_dqi as v
    c = km.Circuit.from_file(path, format='kmb')
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=bsize)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    syn_q, mbuc_b, anc = v.validate_registers(c, spec)
    _ctx.update(km=km, v=v, c=c, spec=spec, tables=tables, syn=km.array(syn_q), mbuc=km.array(mbuc_b), anc=anc)


def run(job):
    setup()
    index, rs = job
    km, v, c, spec, tables = (_ctx[k] for k in ('km', 'v', 'c', 'spec', 'tables'))
    rng = np.random.default_rng(rs)
    syndromes, mbuc, phases = v.sample_shots(spec, tables, rng, bsize)
    acc = np.zeros(NS, dtype=np.int64)
    st = dict(k=0, prev=None, t=0)
    holder = {}

    def cb(target, vals, cond, cvals):
        if target is None:
            return
        t = int(holder['sim'].operation_statistics().num_toffoli_like)
        if st['prev'] is not None:
            acc[st['prev']] += t - st['t']
        st['prev'] = st['k'] % NS
        st['t'] = t
        st['k'] += 1

    sim = km.Simulator(batch_size=bsize, debug_print=cb, counting_operations=True, seed=rs * 7 + 3)
    holder['sim'] = sim
    sim.use_same_registers_as(c)
    bits = v.to_bits(syndromes, spec.degree)
    sim.write_across_shots(_ctx['syn'], bits)
    sim.write_across_shots(_ctx['mbuc'], v.to_bits(mbuc, spec.degree))
    sim.do(c)
    total = int(sim.operation_statistics().num_toffoli_like)
    acc[st['prev']] += total - st['t']
    assert st['k'] == NS * -(-bsize // 256), (st['k'], NS)
    flipped = sim.read_phase_flipped_across_shots()
    diff = sim.read_across_shots(_ctx['syn']) ^ bits
    anc = _ctx['anc']
    junk = sim.read_across_shots(anc) if len(anc) else np.zeros((0, bsize), dtype=bool)
    val = diff.any(axis=0) | junk.any(axis=0)
    ph = (flipped ^ phases) & ~val
    return int(val.sum()) + int(ph.sum()), acc, total


if __name__ == '__main__':
    jobs = [(i, seed * 100000 + i) for i in range(batches)]
    acc = np.zeros(NS, dtype=np.int64)
    errors = total = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for e, a, t in ex.map(run, jobs):
            errors += e
            acc += a
            total += t
    shots = batches * bsize
    with open(out, 'w') as f:
        json.dump(dict(shots=shots, errors=errors, executed=acc.tolist(), total=total), f)
    print('shots', shots, 'errors', errors, 'mean executed Toffolis', total / shots, 'charged', acc.sum() / shots)
