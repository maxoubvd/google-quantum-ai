#!/usr/bin/env python
"""Profiled build of a DQI decoder from a generator module exposing decode(b, S, E).

Same circuit shape as the public carry-q26-892 build.py (kickmix as the kit's `circ`, registers
S0..S69 and E0..E4094). In addition every Toffoli-like builder call is charged to
(top-level b.mark stage, conditional?, Python call stack below the entry), and experiment patches
can be installed between import and construction.

Usage: python dqi_extract_gen.py <decode_module.py> <out_dir> [--patch file.py ...] [--no-circuit]
  A patch file defines install(km, module) and may return a JSON-serialisable report.
Writes <out_dir>/dqi.kmb, build.json, toffoli_profile.json.
"""
import collections
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
sys.dont_write_bytecode = True
import kickmix as km

sys.modules['circ'] = km
args = sys.argv[1:]
entry = Path(args.pop(0)).resolve()
out = Path(args.pop(0)).resolve()
patches, write_circuit = [], True
while args:
    a = args.pop(0)
    if a == '--patch':
        patches.append(Path(args.pop(0)).resolve())
    elif a == '--no-circuit':
        write_circuit = False
    else:
        raise SystemExit('unknown argument ' + a)
out.mkdir(parents=True, exist_ok=True)

PAID = ('ccx', 'ccz', 'cswap', 'init_and')
SKIP = {'<lambda>', 'paid', '<module>', '<listcomp>', '<genexpr>'}
profile = collections.Counter()
stage, depth, cond = ['<none>'], [0], [0]


def width(a, kw):
    """Toffoli-like gates emitted by one builder call: broadcast width, 0 if any operand is classical."""
    n = 1
    for x in list(a) + list(kw.values()):
        try:
            k = len(x)
        except TypeError:
            if type(x).__name__ == 'b':
                return 0
            continue
        if k and type(x[0]).__name__ == 'b':
            return 0
        n = max(n, k)
    return n


def wrap(name):
    original = getattr(km.CircuitBuilder, name)

    def paid(self, *a, **kw):
        n = width(a, kw)
        if n:
            names = []
            f = sys._getframe(1)
            while f is not None:
                co = f.f_code.co_name
                if co not in SKIP:
                    names.append(co + '@' + os.path.basename(f.f_code.co_filename)[:-3])
                f = f.f_back
            profile[(stage[0], bool(cond[0]), tuple(reversed(names)))] += n
        return original(self, *a, **kw)
    setattr(km.CircuitBuilder, name, paid)


for name in PAID:
    wrap(name)
_push, _pop, _mark = km.CircuitBuilder.push_condition, km.CircuitBuilder.pop_condition, km.CircuitBuilder.mark


def push_condition(self, *a, **kw):
    cond[0] += 1
    return _push(self, *a, **kw)


def pop_condition(self, *a, **kw):
    cond[0] -= 1
    return _pop(self, *a, **kw)


class Marked:
    def __init__(self, inner, tag):
        self.inner, self.tag = inner, tag

    def __enter__(self):
        if depth[0] == 0:
            stage[0] = self.tag
        depth[0] += 1
        return self.inner.__enter__()

    def __exit__(self, *a):
        depth[0] -= 1
        if depth[0] == 0:
            stage[0] = '<none>'
        return self.inner.__exit__(*a)


km.CircuitBuilder.push_condition = push_condition
km.CircuitBuilder.pop_condition = pop_condition
km.CircuitBuilder.mark = lambda self, tag, *a, **kw: Marked(_mark(self, tag, *a, **kw), tag)

start = time.time()
sys.path.insert(0, str(entry.parent))
spec = importlib.util.spec_from_file_location('decoder_under_test', entry)
module = importlib.util.module_from_spec(spec)
sys.modules['decoder_under_test'] = module
spec.loader.exec_module(module)
reports = {}
for p in patches:
    ps = importlib.util.spec_from_file_location('dqi_extract_patch_' + p.stem, p)
    pm = importlib.util.module_from_spec(ps)
    sys.modules[ps.name] = pm
    ps.loader.exec_module(pm)
    reports[p.name] = pm.install(km, module)
print('imported in', round(time.time() - start, 1), 's', flush=True)
b = km.CircuitBuilder()
S = [b.alloc_qubits(12) for _ in range(70)]
E = [b.alloc_bits(12) for _ in range(4095)]
module.decode(b, S, E)
print('decoded in', round(time.time() - start, 1), 's', flush=True)
metadata = []
for i in range(70):
    metadata.append(f'REGISTER r{i} "S{i}"')
    metadata.extend(f'APPEND_TO_REGISTER q{12 * i + k} r{i}' for k in range(12))
for j in range(4095):
    metadata.append(f'REGISTER r{70 + j} "E{j}"')
    metadata.extend(f'APPEND_TO_REGISTER b{12 * j + k} r{70 + j}' for k in range(12))
c = km.Circuit('\n'.join(metadata)) + b.finish_circuit()
record = dict(Q=c.num_qubits, T_static=int(c.operation_statistics().num_toffoli_like), num_ops=len(c),
              num_bits=c.num_bits, entry=str(entry), patches=[str(p) for p in patches], patch_reports=reports)
if write_circuit:
    target = out / 'dqi.kmb'
    c.to_file(str(target), format='kmb')
    with open(target, 'rb') as f:
        record['sha256'] = hashlib.file_digest(f, 'sha256').hexdigest()
record['seconds'] = round(time.time() - start, 1)
rows = [dict(stage=s, conditional=cd, stack=list(st), toffoli=n) for (s, cd, st), n in profile.items()]
rows.sort(key=lambda r: -r['toffoli'])
total = sum(r['toffoli'] for r in rows)
record['profile_total'] = total
(out / 'build.json').write_text(json.dumps(record, indent=1, default=str) + '\n')
(out / 'toffoli_profile.json').write_text(json.dumps(dict(total=total, static_T=record['T_static'], rows=rows)))
print(json.dumps({k: v for k, v in record.items() if k != 'patch_reports'}, default=str), flush=True)
assert total == record['T_static'], 'profile does not reconcile with the static count'
