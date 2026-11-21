#!/usr/bin/env python
"""dqi_k1_divprof.py -- K1 round 3: part-by-part anatomy of the un-synthesis division forms (x -> x / y in
GF(2^24), y kept), built from the kit exactly as tools/dqi_kern_divtest.py builds them (copied logic, read-only
imports), in ONE process.

For each form: executed Toffolis per division (simulator, random shots, checked x_new * y == x_old, y kept, no
phase flip, helpers clean), the peak of helpers beyond x and y, and a table of top-level parts: calls, static
Toffolis (unconditioned / classically conditioned), helpers live at entry, peak inside, live at exit.
Usage (through tools/sh/dqi_k1_jtest.sh): dqi_k1_divprof.py <form,form,...> [shots] [seed]
Forms: lateu latef (the line's late-M forms: kit + r6n + latem; with DQI_JUNC_R8P=latem,div6 the div6 inverse)
       and any name in dqi_kern_div.VARIANTS (tw3 twu twf dr dru drf dt42 dt tw2 ...).
"""
import os
import sys
import types

import numpy as np

sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
sys.dont_write_bytecode = True
import kickmix as km                                                   # noqa: E402

sys.modules['circ'] = km
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dqi_junc_r8_entry as entry                                      # noqa: E402

forms = (sys.argv[1] if len(sys.argv) > 1 else 'lateu').split(',')
shots = int(sys.argv[2]) if len(sys.argv) > 2 else 256
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1

walk = sys.modules['walk_bound']
owner = next(f for f in walk.functions(entry.decode)
             if f.__code__.co_filename.endswith('decode_inverse_retained.py') and 'produce_retained' in f.__globals__)
G = owner.__globals__
ar = G['ar']
K = ar.extension_columns()[2]
mod = sys.modules['norm_deferred_flag_return']
r6n = sys.modules['dqi_build_patch_r6n']
latem = sys.modules.get('dqi_junc_r8p_latem')
mul99 = next(m for m in list(sys.modules.values()) if getattr(m, 'tensor99_extension', None) is not None)
TENSOR = mul99.tensor99_extension(K)
import dqi_kern_div as D                                               # noqa: E402

# ---------------------------------------------------------------- part attribution (outermost wrapped call wins)
STACK = []
REC = {}          # name -> dict(calls, unc, cond, entry_max, peak, exit_max)
CUR = dict(unc=0, cond=0, live=0, peak=0, cdepth=0)


def _wrap(name, fn):
    if getattr(fn, '_k1_wrapped', False):
        return fn

    def w(*a, **k):
        top = not STACK
        STACK.append(name)
        if top:
            u0, c0, l0 = CUR['unc'], CUR['cond'], CUR['live']
            CUR['peak'] = l0
        try:
            return fn(*a, **k)
        finally:
            STACK.pop()
            if top:
                r = REC.setdefault(name, dict(calls=0, unc=0, cond=0, entry=0, peak=0, exit=0))
                r['calls'] += 1
                r['unc'] += CUR['unc'] - u0
                r['cond'] += CUR['cond'] - c0
                r['entry'] = max(r['entry'], l0)
                r['peak'] = max(r['peak'], CUR['peak'])
                r['exit'] = max(r['exit'], CUR['live'])
    w._k1_wrapped = True
    w.__name__ = getattr(fn, '__name__', name)
    return w


def wrap_dict(d, names, prefix=''):
    for n in names:
        if n in d and callable(d[n]):
            d[n] = _wrap(prefix + n, d[n])


def wrap_obj(o, names, prefix=''):
    for n in names:
        if hasattr(o, n):
            setattr(o, n, _wrap(prefix + n, getattr(o, n)))


wrap_dict(G, ['multiply_forward99', 'norm_xor_noscratch', 'mul_xor_noscratch', 'produce_retained', 'rebuild_retained',
              'inverse6_xor20', 'phase_inverse6_classical6', 'phase_product24', 'inverse_stream_noscratch'])
wrap_obj(ar, ['mul_xor0', 'phase_product', 'erase_product', '_erase_norm', 'erase_inverse'])
wrap_dict(vars(mod), ['zero_flag'])
for mn, names in (('exp.arithmetic.division', ['zero_flag']), ('exp.jet24.pro_jet', ['cached_inverse']),
                  ('exp.mul99.phase36', ['phase_inverse36_pivot'])):
    try:
        m = __import__(mn, fromlist=['x'])
        wrap_obj(m, names)
    except Exception as e:                                               # noqa: BLE001
        print('wrap skip', mn, e)
wrap_obj(D, ['inv6_iso20', 'inv6_oop'])


def kit_variant(flagged, streamed):
    """The line's division source (r6n edits + latem's late M), as tools/dqi_kern_divtest.py builds it."""
    text = open(mod.__file__).read()
    start = text.index('def divide_from_norms(')
    src = text[start:text.index('\ndef sweep(', start)].rstrip() + '\n'
    edits = (r6n.FLAGGED if flagged else r6n.UNFLAGGED) + (latem.LATE_M if latem else ())
    fn = r6n._variant(mod, src, edits, 'kern_div_' + ('f' if flagged else 'u'))
    own = owner
    if not streamed:
        g2 = {k: v for k, v in G.items() if k not in ('produce_retained', 'norm_keep2')}
        own = types.FunctionType(owner.__code__, g2, owner.__name__)

    def run(b, x, y):
        fm = b.alloc_bits(1)
        b.bit_store(fm, False)
        fn(b, x, y, fm, D.test_before(b, x), lambda: None, own, K)
    return run


class Proxy:
    """Counts Toffoli-like gates (unconditioned / under a classical control or condition) and live helpers."""

    def __init__(self, real):
        self.real = real
        self.qset = set()

    def __getattr__(self, n):
        return getattr(self.real, n)

    def _isbit(self, a):
        if isinstance(a, km.b):
            return True
        if isinstance(a, km.array) and len(a) and isinstance(a[0], km.b):
            return True
        return False

    def _tof(self, args, kw):
        c = CUR['cdepth'] > 0 or any(self._isbit(a) for a in args) or any(self._isbit(v) for v in kw.values())
        CUR['cond' if c else 'unc'] += 1

    def ccx(self, *a, **k):
        self._tof(a, k)
        return self.real.ccx(*a, **k)

    def ccz(self, *a, **k):
        self._tof(a, k)
        return self.real.ccz(*a, **k)

    def cswap(self, *a, **k):
        if not self._isbit(a[0]):
            self._tof(a, k)
        return self.real.cswap(*a, **k)

    def push_condition(self, *a, **k):
        CUR['cdepth'] += 1
        return self.real.push_condition(*a, **k)

    def pop_condition(self, *a, **k):
        CUR['cdepth'] = max(0, CUR['cdepth'] - 1)
        return self.real.pop_condition(*a, **k)

    def alloc_qubits(self, n):
        r = self.real.alloc_qubits(n)
        for q in r:
            self.qset.add(q.id)
        CUR['live'] = len(self.qset)
        CUR['peak'] = max(CUR['peak'], CUR['live'])
        GPEAK[0] = max(GPEAK[0], CUR['live'])
        return r

    def free(self, r):
        items = [r] if isinstance(r, (km.q, km.b)) else list(r)
        for q in items:
            if isinstance(q, km.q):
                self.qset.discard(q.id)
        CUR['live'] = len(self.qset)
        return self.real.free(r)


GPEAK = [0]


def bits(v, n):
    return ((v[None, :] >> np.arange(n)[:, None]) & 1).astype(bool)


def mul(a, c):
    z = np.zeros_like(a)
    for mask, v in TENSOR:
        pa = np.zeros_like(a)
        pc = np.zeros_like(a)
        for i in range(24):
            if (mask >> i) & 1:
                pa ^= (a >> i) & 1
                pc ^= (c >> i) & 1
        z ^= np.where((pa & pc) == 1, v, 0)
    return z


def run_form(form):
    REC.clear()
    for k in ('unc', 'cond', 'live', 'peak', 'cdepth'):
        CUR[k] = 0
    GPEAK[0] = 0
    if form in ('lateu', 'latef'):
        flagged = form == 'latef'
        body = kit_variant(flagged, True)
    else:
        if form in D.VARIANTS:
            fn = D.VARIANTS[form]
        else:
            import dqi_k1_divu as DU                        # K1's door-2 forms
            fn = DU.VARIANTS[form]
        flagged = getattr(fn, 'flagged', False)
    nb = 0 if form in ('lateu', 'latef') else getattr(fn, 'borrowed', 0)
    real = km.CircuitBuilder()
    x = real.alloc_qubits(24)
    y = real.alloc_qubits(24)
    Gw = real.alloc_qubits(nb) if nb else []          # K1: BORROWED wires (random data, must come back unchanged)
    if form not in ('lateu', 'latef'):
        def body(b, x, y):
            if nb:
                fn(b, x, y, K, G, Gw)
            else:
                fn(b, x, y, K, G)
    P = Proxy(real)
    body(P, x, y)
    circ = real.finish_circuit()
    nq = circ.num_qubits
    rng = np.random.default_rng(seed)
    xv = rng.integers(0, 1 << 24, size=shots, dtype=np.int64)
    yv = rng.integers(1, 1 << 24, size=shots, dtype=np.int64)
    xv[3::17] = 0
    yv[5::19] = 1
    yv[7::23] = 1 << 12
    if flagged:
        yv[::13] = 0
        xv[::13] = 0
    sim = km.Simulator(batch_size=shots, debug_print=False, counting_operations=True, seed=seed + 100)
    sim.write_across_shots(km.array(list(x)), bits(xv, 24))
    sim.write_across_shots(km.array(list(y)), bits(yv, 24))
    if nb:
        gv = rng.integers(0, 1 << nb, size=shots, dtype=np.int64)
        sim.write_across_shots(km.array(list(Gw)), bits(gv, nb))
    sim.do(circ)
    state = np.asarray(sim.read_across_shots(km.array([km.q(i) for i in range(nq)]))).astype(np.int64)
    flips = np.asarray(sim.read_phase_flipped_across_shots())
    tof = int(sim.operation_statistics().num_toffoli_like) / shots
    xn = sum(state[q.id] << k for k, q in enumerate(x))
    yn = sum(state[q.id] << k for k, q in enumerate(y))
    okv = (mul(xn, yv) == xv) & (yn == yv)
    if flagged:
        okv = np.where(yv == 0, (xn == int(os.environ.get('DQI_KERN_BEFORE', '0'))) & (yn == yv), okv)
    others = [q for q in range(nq) if q not in set(q.id for q in x) | set(q.id for q in y) | set(q.id for q in Gw)]
    dirty = state[others].any(axis=0) if others else np.zeros(shots, bool)
    gkept = (sum(state[q.id] << k for k, q in enumerate(Gw)) == gv) if nb else np.ones(shots, bool)
    if nb:
        print(f'   borrowed wires: {nb}, restored {int(gkept.sum())}/{shots}', flush=True)
    su = sum(r['unc'] for r in REC.values())
    sc = sum(r['cond'] for r in REC.values())
    print(f'== {form}: executed {tof:.2f} per division; helper peak {GPEAK[0]}; values {int(okv.sum())}/{shots}, '
          f'phase flips {int(flips.sum())}, dirty {int(dirty.sum())}; static in parts {su} unc + {sc} cond, '
          f'outside parts {CUR["unc"] - su} unc + {CUR["cond"] - sc} cond', flush=True)
    for name, r in sorted(REC.items(), key=lambda kv: -kv[1]['unc'] - kv[1]['cond']):
        print(f'   {name:30s} calls {r["calls"]:2d}  tof {r["unc"]:4d} + {r["cond"]:3d}c  live in {r["entry"]:2d} '
              f'peak {r["peak"]:2d} out {r["exit"]:2d}', flush=True)


for f in forms:
    try:
        run_form(f)
    except Exception as e:                                               # noqa: BLE001
        import traceback
        traceback.print_exc()
        print('== FAIL', f, e, flush=True)
