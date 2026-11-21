#!/usr/bin/env python3
"""dqi_pp_paired.py -- workstream dqi_pp (P1): PAIRED-SEED test of a rewritten DQI circuit against its parent.

The official verifier draws its shots from the circuit's own hash and the salt. Here the CHILD's batches (child hash,
salt) are generated exactly as the verifier does (tools/dqi_build_failshot.py), and every batch is run on the child
and on the parent with the same inputs and the same simulator seed. Per shot: value error (syndrome changed or an
ancilla left dirty) or phase error. An exact rewrite fails exactly where the parent fails: a CHILD-ONLY failure is a
failure class of the rewrite. Prints the four counts and the child-only shots.
Usage: python dqi_pp_paired.py <child.kmb> <parent.kmb> <salt|salt0> [workers=8] [batches=all]
"""
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

DATA = os.path.expanduser('~/qcc/data')
sys.path[:0] = [os.path.join(DATA, 'kickmix'), os.path.join(DATA, 'verifiers')]
import kickmix as km  # noqa: E402
from dqi import verify_submission_dqi as v  # noqa: E402

G = {}


def init(child, parent):
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    G['spec'] = spec
    G['tables'] = v.FieldTables(km.GF2Field(spec.degree))
    G['wp'] = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    G['circ'] = {}
    for name, path in (('child', child), ('parent', parent)):
        c = km.Circuit.from_file(path, format='kmb')
        syn_q, mbuc_b, anc = v.validate_registers(c, spec)
        G['circ'][name] = (c, km.array(syn_q), km.array(mbuc_b), anc)


def run(name, syndromes, mbuc, phases, sim_seed):
    c, syn_a, mbuc_a, anc = G['circ'][name]
    spec = G['spec']
    sim = km.Simulator(batch_size=len(phases), debug_print=False, seed=sim_seed)
    sim.use_same_registers_as(c)
    ins = v.to_bits(syndromes, spec.degree)
    sim.write_across_shots(syn_a, ins)
    sim.write_across_shots(mbuc_a, v.to_bits(mbuc, spec.degree))
    sim.do(c)
    flipped = np.asarray(sim.read_phase_flipped_across_shots())
    diff = np.asarray(sim.read_across_shots(syn_a)) ^ ins
    junk = np.asarray(sim.read_across_shots(anc)) if len(anc) else np.zeros((0, len(phases)), dtype=bool)
    val = diff.any(axis=0) | junk.any(axis=0)
    return val, (flipped ^ phases) & ~val


def batch(args):
    i, s0, s1 = args
    spec = G['spec']
    rng = np.random.default_rng(int.from_bytes(s0, 'little'))
    syndromes, mbuc, phases = v.sample_shots(spec, G['tables'], rng, v.DEFAULT_BATCH_SIZE, G['wp'])
    sim_seed = int.from_bytes(s1, 'little')
    vc, pc = run('child', syndromes, mbuc, phases, sim_seed)
    vp, pp = run('parent', syndromes, mbuc, phases, sim_seed)
    fc, fp = vc | pc, vp | pp
    return i, np.flatnonzero(fc & ~fp).tolist(), np.flatnonzero(fp & ~fc).tolist(), int((fc & fp).sum()), \
        int(vc.sum()), int(pc.sum())


def main():
    t0 = time.time()
    child, parent, salt = sys.argv[1], sys.argv[2], sys.argv[3]
    salt = '' if salt == 'salt0' else salt
    W = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    nb = int(sys.argv[5]) if len(sys.argv) > 5 else None
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    c = km.Circuit.from_file(os.environ.get('DQI_PP_SEEDFROM', child), format='kmb')   # fixed shots for bisection
    base_seed = v._initial_seed(v._circuit_sha256(c), salt)
    del c
    starts = range(0, spec.num_shots, v.DEFAULT_BATCH_SIZE)
    seeds = np.random.SeedSequence(base_seed, pool_size=8).generate_state(8 * len(starts),
                                                                          dtype=np.uint64).reshape(-1, 2, 4)
    jobs = [(i, seeds[i, 0].tobytes(), seeds[i, 1].tobytes()) for i in range(len(starts))][:nb]
    if os.environ.get('DQI_PP_BATCHES'):                 # only these batch indices
        keep = set(int(x) for x in os.environ['DQI_PP_BATCHES'].split(','))
        jobs = [j for j in jobs if j[0] in keep]
    co, po, both, vals, phs = [], [], 0, 0, 0
    with Pool(W, initializer=init, initargs=(child, parent)) as pool:
        for i, c_only, p_only, b, nv, nph in pool.imap_unordered(batch, jobs):
            co += [(i, j) for j in c_only]
            po += [(i, j) for j in p_only]
            both += b
            vals += nv
            phs += nph
    print('salt %r: batches %d, child failures %d (value %d, phase %d); both %d, CHILD-ONLY %d, parent-only %d '
          '(%.0fs)' % (salt, len(jobs), vals + phs, vals, phs, both, len(co), len(po), time.time() - t0))
    for i, j in sorted(co)[:20]:
        print('  child-only: batch %d shot %d' % (i, j))


if __name__ == '__main__':
    main()
