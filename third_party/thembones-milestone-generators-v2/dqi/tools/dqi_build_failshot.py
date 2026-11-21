#!/usr/bin/env python
"""Find the failing shots of a DQI circuit under the official verifier's seeding and replay them elsewhere.

Usage: python dqi_build_failshot.py find <circuit.kmb> <salt> <out.npz>
       python dqi_build_failshot.py replay <circuit.kmb> <shots.npz> [repeats=4]
'find' regenerates the verifier's batches (same seeds), records every failing shot (syndrome, masks, phase, error
weight is not stored by the verifier, so the syndrome is kept). 'replay' runs another circuit on those inputs,
padded to one batch, with several simulator seeds, and reports failures per input.
"""
import os
import sys

import numpy as np

DATA = os.path.expanduser('~/qcc/data')
sys.path[:0] = [os.path.join(DATA, 'kickmix'), os.path.join(DATA, 'verifiers')]
import kickmix as km
from dqi import verify_submission_dqi as v

mode, path = sys.argv[1], sys.argv[2]
c = km.Circuit.from_file(path, format='kmb')
spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
B = v.DEFAULT_BATCH_SIZE
syn_q, mbuc_b, anc = v.validate_registers(c, spec)
syn_a, mbuc_a = km.array(syn_q), km.array(mbuc_b)


def run(syndromes, mbuc, phases, sim_seed):
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


if mode == 'find':
    salt, out = sys.argv[3], sys.argv[4]
    weight_probs = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    tables = v.FieldTables(km.GF2Field(spec.degree))
    base_seed = v._initial_seed(v._circuit_sha256(c), salt)
    starts = range(0, spec.num_shots, B)
    seeds = np.random.SeedSequence(base_seed, pool_size=8).generate_state(8 * len(starts), dtype=np.uint64).reshape(-1, 2, 4)
    found = []
    for i, s in enumerate(starts):
        rng = np.random.default_rng(int.from_bytes(seeds[i, 0], 'little'))
        syndromes, mbuc, phases = v.sample_shots(spec, tables, rng, B, weight_probs)
        val, ph = run(syndromes, mbuc, phases, int.from_bytes(seeds[i, 1], 'little'))
        for j in np.flatnonzero(val | ph):
            found.append((syndromes[j], mbuc[j], phases[j], bool(val[j])))
            print('batch', i, 'shot', int(j), 'value error' if val[j] else 'phase error', flush=True)
    np.savez_compressed(out, syndromes=np.array([f[0] for f in found]), mbuc=np.array([f[1] for f in found]),
                        phases=np.array([f[2] for f in found]))
    print('failing shots', len(found))
else:
    d = np.load(sys.argv[3])
    reps = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    n = len(d['phases'])
    for r in range(reps):
        val, ph = run(d['syndromes'], d['mbuc'], d['phases'], 1000 + r)
        print('seed', r, 'value errors', val.astype(int).tolist(), 'phase errors', ph.astype(int).tolist())
