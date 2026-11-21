#!/usr/bin/env python3
"""dqi_int_rate.py -- workstream dqi_int (B5): FAILED SHOTS PER RUN of a DQI circuit and its parent on the SAME shots.

Copy of P1's paired-seed test (tools/dqi_pp_paired.py, credit P1) extended to many salts in one worker pool and to a
full record of every failing shot of BOTH circuits. For each salt the CHILD's batches are drawn exactly as the official
verifier draws them for the child file (seed from the child's sha256 and the salt), and every batch is run on the
child and on the parent with the same syndromes, masks, phases and simulator seed. A failure is a value error
(syndrome changed or an ancilla left dirty) or a phase error, as the verifier counts them.
Usage: python dqi_int_rate.py <child.kmb> <parent.kmb> <salts> <workers> <out.tsv>
  <salts>: comma list; 'salt0' = the empty salt; 'alt-salt-A..B' expands to alt-salt-A .. alt-salt-B.
One line per salt is appended to out.tsv (whole lines only):
  salt  child_value  child_phase  parent_value  parent_phase  both  child_only  parent_only  seconds  shots
  where shots = batch:shot:CV|CP|C-:PV|PP|P- for every shot that fails on either circuit.
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
    k, i, s0, s1 = args
    spec = G['spec']
    rng = np.random.default_rng(int.from_bytes(s0, 'little'))
    syndromes, mbuc, phases = v.sample_shots(spec, G['tables'], rng, v.DEFAULT_BATCH_SIZE, G['wp'])
    sim_seed = int.from_bytes(s1, 'little')
    vc, pc = run('child', syndromes, mbuc, phases, sim_seed)
    vp, pp = run('parent', syndromes, mbuc, phases, sim_seed)
    rows = []
    for j in np.flatnonzero(vc | pc | vp | pp):
        rows.append((int(j), 'CV' if vc[j] else ('CP' if pc[j] else 'C-'), 'PV' if vp[j] else ('PP' if pp[j] else 'P-')))
    return k, i, rows


def salts_of(spec):
    out = []
    for s in spec.split(','):
        if '..' in s:
            a, b = s.split('..')
            pre = a.rstrip('0123456789')
            out += ['%s%d' % (pre, n) for n in range(int(a[len(pre):]), int(b) + 1)]
        elif s:
            out.append(s)
    return out


def main():
    child, parent, salts, W, out = sys.argv[1], sys.argv[2], salts_of(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    c = km.Circuit.from_file(child, format='kmb')
    sha = v._circuit_sha256(c)
    del c
    starts = range(0, spec.num_shots, v.DEFAULT_BATCH_SIZE)
    jobs, meta = [], []
    for k, salt in enumerate(salts):
        s = '' if salt == 'salt0' else salt
        seeds = np.random.SeedSequence(v._initial_seed(sha, s), pool_size=8).generate_state(
            8 * len(starts), dtype=np.uint64).reshape(-1, 2, 4)
        jobs += [(k, i, seeds[i, 0].tobytes(), seeds[i, 1].tobytes()) for i in range(len(starts))]
        meta.append(salt)
    res = {k: [] for k in range(len(salts))}
    left = {k: len(starts) for k in range(len(salts))}
    t0 = {k: None for k in range(len(salts))}
    tstart = time.time()
    with Pool(W, initializer=init, initargs=(child, parent)) as pool:
        for k, i, rows in pool.imap(batch, jobs, chunksize=1):
            if t0[k] is None:
                t0[k] = time.time()
            res[k] += [(i, j, a, b) for j, a, b in rows]
            left[k] -= 1
            if left[k] == 0:
                r = res[k]
                cv = sum(a == 'CV' for _, _, a, _ in r)
                cp = sum(a == 'CP' for _, _, a, _ in r)
                pv = sum(b == 'PV' for _, _, _, b in r)
                pph = sum(b == 'PP' for _, _, _, b in r)
                both = sum(a != 'C-' and b != 'P-' for _, _, a, b in r)
                co = sum(a != 'C-' and b == 'P-' for _, _, a, b in r)
                po = sum(a == 'C-' and b != 'P-' for _, _, a, b in r)
                shots = ' '.join('%d:%d:%s:%s' % x for x in sorted(r)) or '-'
                line = '\t'.join(map(str, (meta[k], cv, cp, pv, pph, both, co, po,
                                           '%.0f' % (time.time() - (t0[k] or tstart)), shots)))
                with open(out, 'a') as f:
                    f.write(line + '\n')
                print(line, flush=True)
    print('done %d salts in %.0f s; child %s parent %s' % (len(salts), time.time() - tstart, sha[:16], parent),
          flush=True)


if __name__ == '__main__':
    main()
