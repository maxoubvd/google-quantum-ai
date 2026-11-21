#!/usr/bin/env python3
"""dqi_int_flagidle.py -- workstream dqi_int (B5, round 3): for the converter's peak moments, which one-qubit flags
are IDLE there (no gate touches them) and over which interval, so that one of them can be parked (X-measured) and
rebuilt around the moment.

Same harness as tools/dqi_int_flagtrace.py (converter alone: convert, then unconvert with the forward's token, in the
environment of tools/sh/dqi_int_env.sh). Every builder method that takes qubits is wrapped: each call is one "op";
for every qubit lifetime (index, birth) the ops that touch it are recorded with their call site. At every allocation
that brings the live count to >= AT a "moment" is recorded. For each moment and each live flag the idle interval
(last touch before, first touch after) is printed, with the live count at the interval's end (room for a rebuild).
Usage: python dqi_int_flagidle.py <AT> [depth=5]
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
DEPTH = int(sys.argv[2]) if len(sys.argv) > 2 else 5

spec = importlib.util.spec_from_file_location('dqi_int_entry', os.path.join(TOOLS, 'dqi_int_entry.py'))
mod = importlib.util.module_from_spec(spec)
sys.modules['dqi_int_entry'] = mod
spec.loader.exec_module(mod)
conv = mod.decode.__globals__['g']['convert_safe']

CB = km.CircuitBuilder
ORIG = {n: getattr(CB, n) for n in dir(CB) if not n.startswith('__')}
st = dict(op=0, half=None, nlive=0)
cur = {}            # qubit index -> lifetime id
lives = []          # lifetime id -> dict(idx, n, desc, born_op, half, touches=[(op, name, site, nlive)], end)
moments = []        # (op, half, nlive, [lifetime ids of live flags], regs descriptor)


def site(depth=DEPTH):
    f = sys._getframe(2)
    names, m = [], None
    while f is not None and len(names) < depth:
        fn = f.f_code.co_filename
        if 'kickmix' not in fn and 'flagidle' not in fn:
            names.append('%s:%d %s' % (os.path.basename(fn), f.f_lineno, f.f_code.co_name))
            if m is None:
                v = f.f_locals.get('m')
                if isinstance(v, int):
                    m = v
        f = f.f_back
    return m, ' < '.join(names)


def qidx(x, out):
    if isinstance(x, (list, tuple)):
        for y in x:
            qidx(y, out)
        return
    s = str(x)
    if len(s) > 1 and s[0] == 'q' and s[1:].isdigit():
        out.append(int(s[1:]))


def alloc(self, n):
    regs = ORIG['alloc_qubits'](self, n)
    m, s = site()
    for q in regs:
        i = int(str(q)[1:])
        lid = len(lives)
        lives.append(dict(idx=i, n=n, m=m, site=s, born_op=st['op'], half=st['half'], touches=[], end=None))
        cur[i] = lid
    st['nlive'] += n
    if st['half'] is not None and st['nlive'] >= AT:
        fl = [cur[i] for i in cur if lives[cur[i]]['n'] == 1 and i >= 840]
        regs_c = collections.Counter('%s x%d' % (lives[cur[i]]['site'].split(' < ')[0], lives[cur[i]]['n'])
                                     for i in cur if lives[cur[i]]['n'] > 1 and i >= 840)
        moments.append((st['op'], st['half'], st['nlive'], sorted(fl), tuple(sorted(regs_c.items()))))
    return regs


def free(self, regs):
    out = []
    qidx(regs, out)
    for i in out:
        lid = cur.pop(i, None)
        if lid is not None:
            lives[lid]['end'] = st['op']
            st['nlive'] -= 1
    return ORIG['free'](self, regs)


def wrap(name):
    fn = ORIG[name]

    def call(self, *a, **k):
        out = []
        qidx(list(a) + list(k.values()), out)
        st['op'] += 1
        if out:
            s = None
            for i in out:
                lid = cur.get(i)
                if lid is not None and lives[lid]['n'] == 1 and i >= 840:
                    if s is None:
                        s = site(3)[1]
                    lives[lid]['touches'].append((st['op'], name, s, st['nlive']))
        return fn(self, *a, **k)
    return call


SKIP = {'alloc_qubits', 'free', 'alloc_bits', 'mark', 'push_condition', 'pop_condition', 'finish_circuit',
        'max_qubits', 'set_max_qubits', 'num_allocated_qubits', 'num_free_qubits', 'debug_print', 'info_tag',
        'bit_store', 'create_classical_register', 'create_quantum_register', 'reset'}
for n in ORIG:
    if n not in SKIP and callable(ORIG[n]):
        setattr(CB, n, wrap(n))
CB.alloc_qubits = alloc
CB.free = free

b = km.CircuitBuilder()
A = [b.alloc_qubits(12) for _ in range(35)]
B = [b.alloc_qubits(12) for _ in range(35)]
st['half'] = 'convert'
ct = conv(b, A, B, monic_root=True)
st['half'] = 'unconvert'
conv(b, A, B, True, phase_token=ct)
st['half'] = None
print('moments at >= %d: %d (convert %d, unconvert %d), ops %d' % (
    AT, len(moments), sum(1 for x in moments if x[1] == 'convert'), sum(1 for x in moments if x[1] == 'unconvert'),
    st['op']), flush=True)


def idle(lid, op):
    t = lives[lid]['touches']
    prev = [x for x in t if x[0] <= op]
    nxt = [x for x in t if x[0] > op]
    return (prev[-1] if prev else None), (nxt[0] if nxt else None)


import pickle
OUT = os.environ.get('DQI_INT_IDLE_OUT', os.path.expanduser('~/qcc/ws/dqi_int/logs/flagidle.pkl'))
keep = set(l for mo in moments for l in mo[3])
pickle.dump(dict(moments=moments, lives={l: lives[l] for l in keep}, ops=st['op']), open(OUT, 'wb'))
print('saved', OUT, 'flags', len(keep))
