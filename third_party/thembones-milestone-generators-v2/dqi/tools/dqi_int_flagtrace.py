#!/usr/bin/env python3
"""dqi_int_flagtrace.py -- workstream dqi_int: which one-qubit flags are live at the converter's peak moments.

Builds the flow's converter alone (convert, then unconvert with the forward's phase token, as the decode calls them;
the method of tools/dqi_x12d_convtest24.py) in the environment of tools/sh/dqi_int_env.sh, with the builder's
alloc_qubits / free wrapped (the method of tools/dqi_x12d_alloctrace.py). For every allocation that brings the live
count to >= AT it records the live set: multi-qubit registers by allocating site, one-qubit flags by
(level m, allocating call chain). Prints the distinct live sets per half, and the life of every flag live there.
Usage: python dqi_int_flagtrace.py <AT> [chain depth=7]
"""
import collections
import importlib.util
import os
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
import kickmix as km                                                       # noqa: E402
sys.modules['circ'] = km

AT = int(sys.argv[1])
DEPTH = int(sys.argv[2]) if len(sys.argv) > 2 else 7

spec = importlib.util.spec_from_file_location('dqi_int_entry', os.path.join(TOOLS, 'dqi_int_entry.py'))
mod = importlib.util.module_from_spec(spec)
sys.modules['dqi_int_entry'] = mod
spec.loader.exec_module(mod)
conv = mod.decode.__globals__['g']['convert_safe']

CB = km.CircuitBuilder
_alloc, _free = CB.alloc_qubits, CB.free
live = {}                     # index -> descriptor
born = {}                     # index -> (seq, half)
state = dict(half=None, seq=0, maxlive=collections.Counter(), sets={})


def describe(n):
    f = sys._getframe(2)
    names, m = [], None
    while f is not None and len(names) < DEPTH:
        fn = f.f_code.co_filename
        if 'kickmix' not in fn and 'flagtrace' not in fn:
            names.append('%s:%d %s' % (os.path.basename(fn), f.f_lineno, f.f_code.co_name))
            if m is None:
                v = f.f_locals.get('m')
                if isinstance(v, int):
                    m = v
        f = f.f_back
    if n == 1:
        return ('flag', m, ' < '.join(names))
    return ('reg', n, names[0] if names else '?')


def idx(q):
    if isinstance(q, int):
        return q
    s = str(q)
    return int(s[1:]) if s.startswith('q') else None


def alloc(self, n):
    regs = _alloc(self, n)
    d = describe(n)
    state['seq'] += 1
    for q in regs:
        live[idx(q)] = d
        born[idx(q)] = (state['seq'], state['half'])
    h = state['half']
    if h is not None:
        state['maxlive'][h] = max(state['maxlive'][h], len(live))
        if len(live) >= AT:
            regs_c = collections.Counter('%s x%d' % (d[2], d[1]) for k, d in live.items() if k >= 840 and d[0] == 'reg')
            flags = tuple(sorted((d[1] if d[1] is not None else -1, d[2], born[k][0], born[k][1])
                                 for k, d in live.items() if k >= 840 and d[0] == 'flag'))
            key = (tuple(sorted(regs_c.items())), tuple((f[0], f[1]) for f in flags))
            s = state['sets'].setdefault(h, {})
            ent = s.setdefault(key, [0, len(live), flags])
            ent[0] += 1
            ent[1] = max(ent[1], len(live))
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


CB.alloc_qubits = alloc
CB.free = free
b = km.CircuitBuilder()
A = [b.alloc_qubits(12) for _ in range(35)]
B = [b.alloc_qubits(12) for _ in range(35)]
state['half'] = 'convert'
ct = conv(b, A, B, monic_root=True)
state['half'] = 'unconvert'
conv(b, A, B, True, phase_token=ct)
state['half'] = None
c = b.finish_circuit()
print('qubits %d, static Toffoli-like %d, live peaks %s' % (c.num_qubits, c.operation_statistics().num_toffoli_like,
                                                          dict(state['maxlive'])), flush=True)
for h, sets in state['sets'].items():
    print('### %s: %d distinct live sets at >= %d' % (h, len(sets), AT))
    flagsseen = {}
    for key, (cnt, n, flags) in sorted(sets.items(), key=lambda kv: -kv[1][1]):
        regs = '; '.join('%s %d' % (s, k) for s, k in key[0])
        print('  %d (x%d): %s | flags m=%s' % (n, cnt, regs, [f[0] for f in flags]))
        for f in flags:
            flagsseen[(f[0], f[1])] = (f[2], f[3])
    print('  flags live at these moments (level m, born seq / half, chain):')
    for (m, ch), (seq, half) in sorted(flagsseen.items(), key=lambda kv: kv[1][0]):
        print('    m=%s seq %d (%s): %s' % (m, seq, half, ch[:600]))
