#!/usr/bin/env python3
"""dqi_x12d_alloctrace.py -- qubit allocation trace of a DQI build (workstream dqi_x12d).

Builds the circuit of <entry.py> in memory (as tools/dqi_convert_build.py does, nothing written) with the builder's
alloc_qubits / free wrapped: per top-level stage mark, the maximum number of allocated qubits (the live count)
against the maximum index + 1 (what the score counts).  At the moment the highest index is first reached in each
stage, it prints the allocated registers above index 840 grouped by allocating call site (file:line of the first
frame outside kickmix), so fragmentation can be traced to the register that sits high.
Usage: python dqi_x12d_alloctrace.py <entry.py> [stage ...]   (environment as for the build)
"""
import collections
import importlib.util
import os
import sys
import traceback

sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
import kickmix as km                                                       # noqa: E402
sys.modules['circ'] = km

entry = os.path.abspath(sys.argv[1])
want = set(sys.argv[2:])
sys.path.insert(0, os.path.dirname(entry))
spec = importlib.util.spec_from_file_location('decoder_under_test', entry)
mod = importlib.util.module_from_spec(spec)
sys.modules['decoder_under_test'] = mod
spec.loader.exec_module(mod)

CB = km.CircuitBuilder
_alloc, _free, _mark = CB.alloc_qubits, CB.free, CB.mark
live = {}                      # qubit index -> call site
state = dict(stage=None, depth=0, maxlive=collections.Counter(), maxidx=collections.Counter(), snap={})


def site():
    for fr in reversed(traceback.extract_stack()[:-2]):
        if 'kickmix' not in fr.filename and 'alloctrace' not in fr.filename:
            return '%s:%d %s' % (os.path.basename(fr.filename), fr.lineno, fr.name)
    return '?'


def idx(q):
    """Index of a qubit; None for a classical bit (b.free also releases bits)."""
    if isinstance(q, int):
        return q
    s = str(q)
    return int(s[1:]) if s.startswith('q') else None


def alloc(self, n):
    regs = _alloc(self, n)
    s = site()
    for q in regs:
        live[idx(q)] = s
    st = state['stage']
    at = int(os.environ.get('DQI_X12D_TRACE_AT', '0'))         # record every distinct live set >= this count
    if st is not None and at and len(live) >= at and (not want or st in want):
        key = tuple(sorted(collections.Counter(v for k, v in live.items() if v != '?').items()))
        sets = state.setdefault('atsets', {}).setdefault(st, {})
        sets[key] = max(sets.get(key, 0), len(live))
    if st is not None:
        if len(live) > state['maxlive'][st] and (not want or st in want):
            if os.environ.get('DQI_X12D_TRACE_ALL') == '1':        # every live register but the 70 input words
                hi = collections.Counter(v for k, v in live.items() if v != '?')
            else:
                hi = collections.Counter(v for k, v in live.items() if k >= 840)
            state.setdefault('livesnap', {})[st] = (len(live), hi.most_common(18))
        state['maxlive'][st] = max(state['maxlive'][st], len(live))
        m = max(live) + 1
        if m > state['maxidx'][st]:
            state['maxidx'][st] = m
            if not want or st in want:
                hi = collections.Counter(v for k, v in live.items() if k >= 840)
                state['snap'][st] = (m, len(live), hi.most_common(12))
    return regs


def free(self, regs):
    try:
        it = list(regs)
    except TypeError:
        it = [regs]
    for q in it:
        i = idx(q)
        if i is not None:
            live.pop(i, None)
    return _free(self, regs)


class _Mark:
    def __init__(self, b, name):
        self.b, self.name, self.inner = b, name, _mark(b, name)

    def __enter__(self):
        if state['depth'] == 0:
            state['stage'] = self.name
        state['depth'] += 1
        return self.inner.__enter__()

    def __exit__(self, *a):
        state['depth'] -= 1
        return self.inner.__exit__(*a)


CB.alloc_qubits = alloc
CB.free = free
CB.mark = lambda self, name: _Mark(self, name)
b = km.CircuitBuilder()
S = [b.alloc_qubits(12) for _ in range(70)]
E = [b.alloc_bits(12) for _ in range(4095)]
mod.decode(b, S, E)
for st in state['maxidx']:
    print('%-16s max live %4d   max index+1 %4d' % (st, state['maxlive'][st], state['maxidx'][st]))
for st, (m, n, hi) in state['snap'].items():
    print('--- %s: index peak %d reached with %d allocated; registers above 840 by call site:' % (st, m, n))
    for s, c in hi:
        print('   %4d  %s' % (c, s))
for st, (n, hi) in state.get('livesnap', {}).items():
    print('=== %s: LIVE peak %d; registers above 840 by call site:' % (st, n))
    for s, c in hi:
        print('   %4d  %s' % (c, s))
for st, sets in state.get('atsets', {}).items():
    print('### %s: %d distinct live sets at or above %s' % (st, len(sets), os.environ.get('DQI_X12D_TRACE_AT')))
    for key, n in sorted(sets.items(), key=lambda kv: -kv[1]):
        print('  %d: %s' % (n, '; '.join('%s %d' % (s, c) for s, c in key)))
