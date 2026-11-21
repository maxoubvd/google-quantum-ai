#!/usr/bin/env python
"""Per-stage mean EXECUTED Toffoli-like gates of a DQI generator build (workstream dqi_build).

The builder's finish_circuit() is cumulative and repeatable, so at every top-level b.mark exit the circuit
built so far is simulated on a fixed set of verifier-style shots with operation counting; the differences of
the cumulative means are the per-stage executed means. Static counts are recorded alongside.

Usage: python dqi_build_stage_means.py <entry.py> <out.json> [--shots 512] [--patch file.py ...]
  A patch file defines install(km, module), as for dqi_extract_gen.py.
"""
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

DATA = os.path.expanduser('~/qcc/data')
sys.path[:0] = [os.path.join(DATA, 'kickmix'), os.path.join(DATA, 'verifiers')]
sys.dont_write_bytecode = True
import kickmix as km
from dqi import verify_submission_dqi as v

sys.modules['circ'] = km
args = sys.argv[1:]
entry = Path(args.pop(0)).resolve()
out = Path(args.pop(0)).resolve()
shots, patches = 512, []
while args:
    a = args.pop(0)
    if a == '--shots':
        shots = int(args.pop(0))
    elif a == '--patch':
        patches.append(Path(args.pop(0)).resolve())
    else:
        raise SystemExit('unknown argument ' + a)
BATCH = 256
assert shots % BATCH == 0

metadata = []
for i in range(70):
    metadata.append(f'REGISTER r{i} "S{i}"')
    metadata.extend(f'APPEND_TO_REGISTER q{12 * i + k} r{i}' for k in range(12))
for j in range(4095):
    metadata.append(f'REGISTER r{70 + j} "E{j}"')
    metadata.extend(f'APPEND_TO_REGISTER b{12 * j + k} r{70 + j}' for k in range(12))
header = km.Circuit('\n'.join(metadata))

spec = v.ProblemSpecDQI.benchmark_instance(num_shots=shots)
weight_probs = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
tables = v.FieldTables(km.GF2Field(spec.degree))
rows, depth, builder = [], [0], [None]
start = time.time()


def measure(tag):
    c = header + builder[0].finish_circuit()
    static = int(c.operation_statistics().num_toffoli_like)
    syn_q, mbuc_b, _ = v.validate_registers(c, spec)
    syn_a, mbuc_a = km.array(syn_q), km.array(mbuc_b)
    total = 0
    for k in range(shots // BATCH):
        rng = np.random.default_rng(20261111 + k)
        syndromes, mbuc, _ = v.sample_shots(spec, tables, rng, BATCH, weight_probs)
        sim = km.Simulator(batch_size=BATCH, debug_print=False, counting_operations=True, seed=977 + k)
        sim.use_same_registers_as(c)
        sim.write_across_shots(syn_a, v.to_bits(syndromes, spec.degree))
        sim.write_across_shots(mbuc_a, v.to_bits(mbuc, spec.degree))
        sim.do(c)
        total += int(sim.operation_statistics().num_toffoli_like)
    rows.append(dict(stage=tag, cum_static=static, cum_mean=total / shots, qubits=c.num_qubits, ops=len(c)))
    print(f'{tag:22s} cum static {static:8d}  cum mean {total / shots:11.1f}  [{time.time() - start:6.0f}s]', flush=True)
    del c


_mark = km.CircuitBuilder.mark


class Marked:
    def __init__(self, inner, tag):
        self.inner, self.tag = inner, tag

    def __enter__(self):
        depth[0] += 1
        return self.inner.__enter__()

    def __exit__(self, *a):
        r = self.inner.__exit__(*a)
        depth[0] -= 1
        if depth[0] == 0 and a[0] is None:
            measure(self.tag)
        return r


km.CircuitBuilder.mark = lambda self, tag, *a, **kw: Marked(_mark(self, tag, *a, **kw), tag)
sys.path.insert(0, str(entry.parent))
spec_m = importlib.util.spec_from_file_location('decoder_under_test', entry)
module = importlib.util.module_from_spec(spec_m)
sys.modules['decoder_under_test'] = module
spec_m.loader.exec_module(module)
for p in patches:
    ps = importlib.util.spec_from_file_location('dqi_build_patch_' + p.stem, p)
    pm = importlib.util.module_from_spec(ps)
    sys.modules[ps.name] = pm
    ps.loader.exec_module(pm)
    pm.install(km, module)
b = km.CircuitBuilder()
builder[0] = b
S = [b.alloc_qubits(12) for _ in range(70)]
E = [b.alloc_bits(12) for _ in range(4095)]
module.decode(b, S, E)
prev_s = prev_m = 0
table = []
for r in rows:
    table.append(dict(stage=r['stage'], static=r['cum_static'] - prev_s, mean=round(r['cum_mean'] - prev_m, 1)))
    prev_s, prev_m = r['cum_static'], r['cum_mean']
out.write_text(json.dumps(dict(entry=str(entry), shots=shots, qubits=rows[-1]['qubits'], total_static=prev_s,
                               total_mean=prev_m, stages=table, seconds=round(time.time() - start, 1)), indent=1) + '\n')
print('\nstage                    static       mean')
for t in table:
    if t['static']:
        print(f"{t['stage']:22s} {t['static']:8d} {t['mean']:10.1f}")
print(f"{'total':22s} {prev_s:8d} {prev_m:10.1f}   qubits {rows[-1]['qubits']}")
