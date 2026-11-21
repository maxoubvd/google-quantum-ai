#!/usr/bin/env python3
"""dqi_zn_flagidle.py -- builder Z1 (Nov 14, fourth holder): which one-qubit flags are IDLE at the converter's
peak moments, in the lct environment (the line's converter alone as r12hz's decode calls it, plus my ZN_LCT_P and
Z2's Z2_P installers), so that one of them can be parked (X-measured) and rebuilt around the moments.
Mechanism of B5's tools/dqi_int_flagidle.py (read-only, not imported): every builder method that takes qubits is
wrapped; for every one-qubit lifetime (index >= 840) the ops that touch it are recorded with their call site; every
allocation that brings the live count to >= AT is a "moment".
Usage (in the line's env, e.g. bash sh/dqi_zn_lenv.sh <cores> python -B dqi_zn_flagidle.py <AT> [out.txt]):
  output: per group of moments (half, allocation site), the live flags with their birth site and the idle
  interval around the group (last touch before the first moment, first touch after the last one).
"""
import collections
import importlib
import importlib.util
import os
import sys
import types

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import dqi_int_convtest as CT  # noqa: E402  (kickmix as 'circ', live tracking of the harness)

km = CT.km
AT = int(sys.argv[1])
OUT = sys.argv[2] if len(sys.argv) > 2 else None
DEPTH = 6


def site(depth=DEPTH, skip=2):
    f = sys._getframe(skip)
    names, m = [], None
    while f is not None and len(names) < depth:
        fn = f.f_code.co_filename
        if 'kickmix' not in fn and 'flagidle' not in fn and 'convtest' not in fn:
            names.append('%s:%d:%s' % (os.path.basename(fn), f.f_lineno, f.f_code.co_name))
            if m is None:
                v = f.f_locals.get('m')
                if isinstance(v, int):
                    m = v
        f = f.f_back
    return m, ' < '.join(names)


def qidx(x, out):
    if isinstance(x, (list, tuple)) or type(x).__name__ == 'array':
        for y in x:
            qidx(y, out)
        return
    s = str(x)
    if len(s) > 1 and s[0] == 'q' and s[1:].isdigit():
        out.append(int(s[1:]))


st = dict(op=0, half=None, nlive=0)
cur = {}
lives = []
moments = []


def main():
    spec = importlib.util.spec_from_file_location('dqi_sw_entry', os.path.join(TOOLS, 'dqi_sw_entry.py'))
    mod = importlib.util.module_from_spec(spec)
    sys.modules['dqi_sw_entry'] = mod
    spec.loader.exec_module(mod)
    for env, pre in (('ZN_LCT_P', 'dqi_zn_'), ('Z2_P', 'dqi_z2_')):
        for name in [p for p in os.environ.get(env, '').split(',') if p]:
            m = importlib.import_module(pre + name)
            print('installed', name, m.install(km, types.SimpleNamespace(decode=mod.decode)), flush=True)
    conv = mod.decode.__globals__['g']['convert_safe']
    saved = CT.hankel_saved(mod.decode)
    CB = km.CircuitBuilder
    ORIG = {n: getattr(CB, n) for n in dir(CB) if not n.startswith('__')}

    def alloc(self, n):
        regs = ORIG['alloc_qubits'](self, n)
        idx = []
        qidx(regs, idx)
        m, s = site()
        for i in idx:
            lid = len(lives)
            lives.append(dict(idx=i, n=n, m=m, site=s, born=st['op'], half=st['half'], touches=[], end=None))
            cur[i] = lid
        st['nlive'] += len(idx)
        if st['half'] is not None and st['nlive'] >= AT:
            fl = sorted(cur[i] for i in cur if lives[cur[i]]['n'] == 1 and i >= 840)
            words = collections.Counter('%s x%d' % (lives[cur[i]]['site'].split(' < ')[0], lives[cur[i]]['n'])
                                        for i in cur if lives[cur[i]]['n'] > 1 and i >= 840)
            moments.append(dict(op=st['op'], half=st['half'], live=st['nlive'], m=m, site=s, flags=fl,
                                words=dict(words)))
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
            s = None
            for i in out:
                lid = cur.get(i)
                if lid is not None and lives[lid]['n'] == 1 and i >= 840:
                    if s is None:
                        s = site(4, 2)[1]
                    lives[lid]['touches'].append((st['op'], name, s, st['nlive']))
            return fn(self, *a, **k)
        return call

    skip = {'alloc_qubits', 'free', 'alloc_bits', 'mark', 'push_condition', 'pop_condition', 'finish_circuit',
            'max_qubits', 'set_max_qubits', 'num_allocated_qubits', 'num_free_qubits', 'debug_print', 'info_tag',
            'bit_store', 'create_classical_register', 'create_quantum_register', 'reset'}
    for n in ORIG:
        if n not in skip and callable(ORIG[n]):
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
    saved.clear()
    report()


def touch_around(lid, op0, op1):
    t = lives[lid]['touches']
    prev = [x for x in t if x[0] <= op0]
    nxt = [x for x in t if x[0] > op1]
    inside = [x for x in t if op0 < x[0] <= op1]
    return (prev[-1] if prev else None), (nxt[0] if nxt else None), inside


def report():
    w = open(OUT, 'w') if OUT else sys.stdout
    pr = lambda *a: print(*a, file=w, flush=True)  # noqa: E731
    pr('moments at >= %d: %d (convert %d, unconvert %d), ops %d, peak flags tracked %d' % (
        AT, len(moments), sum(1 for x in moments if x['half'] == 'convert'),
        sum(1 for x in moments if x['half'] == 'unconvert'), st['op'], len(lives)))
    groups = collections.OrderedDict()
    for mo in moments:
        key = (mo['half'], mo['m'], ' < '.join(mo['site'].split(' < ')[:3]))
        groups.setdefault(key, []).append(mo)
    for key, mos in groups.items():
        op0, op1 = mos[0]['op'], mos[-1]['op']
        pr('\n== %s level %s at %s: %d moments, live max %d, ops %d..%d' % (
            key[0], key[1], key[2], len(mos), max(x['live'] for x in mos), op0, op1))
        pr('   words live (first moment): %s' % ', '.join('%s:%d' % kv for kv in sorted(mos[0]['words'].items())))
        common = set(mos[0]['flags'])
        for x in mos[1:]:
            common &= set(x['flags'])
        for lid in sorted(common):
            L = lives[lid]
            p, nx, ins = touch_around(lid, op0, op1)
            pr('   flag q%d born op %d (%s, level %s) %s' % (L['idx'], L['born'], L['half'], L['m'],
                                                             L['site'].split(' < ')[0]))
            pr('      last touch before: %s' % (('op %d %s %s live %d' % p) if p else '-'))
            pr('      next touch after:  %s' % (('op %d %s %s live %d' % nx) if nx else '-'))
            if ins:
                pr('      touched INSIDE the group %d times, first: op %d %s %s' % (len(ins), *ins[0][:3]))
    if OUT:
        w.close()


if __name__ == '__main__':
    main()
