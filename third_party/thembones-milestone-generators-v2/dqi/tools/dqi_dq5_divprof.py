#!/usr/bin/env python
"""dqi_dq5_divprof.py -- DQ5 (Nov 18): EXECUTED Toffolis per top-level part of the un-synthesis division forms
(x -> x / y in GF(2^24), y kept), and variant forms written as source edits of the kit's divide_from_norms.

Built exactly as tools/dqi_k1_divprof.py builds the line's late-M forms (K1's logic, copied; read-only imports).
New here: (1) every top-level part ends with an INFO_TAG, and the executed count of each part is read by simulating
the circuit's prefixes with the same seed and inputs (count(prefix k) - count(prefix k-1)); (2) init_and is counted
(K1's proxy counted only ccx / ccz / cswap); (3) variants: DQ5 edits applied on top of r6n + latem.
Usage (through tools/sh/dqi_k1_jtest.sh, DQI_JUNC_R8P=latem): dqi_dq5_divprof.py <form,form,...> [shots] [seed] [noparts]
Forms: lateu latef (the line), and lateu+<variant> / latef+<variant> for the variants in dqi_dq5_divvar.VARIANTS.
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
PARTS = not (len(sys.argv) > 4 and sys.argv[4] == 'noparts')

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
import dqi_dq5_divvar as V                                             # noqa: E402

# ---------------------------------------------------------------- part markers (outermost wrapped call only)
STACK = []
BUILDER = [None]
NAMES = []


def _wrap(name, fn):
    if getattr(fn, '_dq5_wrapped', False):
        return fn

    def w(*a, **k):
        STACK.append(name)
        try:
            return fn(*a, **k)
        finally:
            STACK.pop()
            if not STACK and BUILDER[0] is not None and PARTS:
                NAMES.append(name)
                BUILDER[0].info_tag('DQ5PART %d %s' % (len(NAMES) - 1, name))
    w._dq5_wrapped = True
    w.__wrapped__ = fn
    w.__name__ = getattr(fn, '__name__', name)
    return w


def wrap_dict(d, names):
    for n in names:
        if n in d and callable(d[n]):
            d[n] = _wrap(n, d[n])


def wrap_obj(o, names):
    for n in names:
        if hasattr(o, n):
            setattr(o, n, _wrap(n, getattr(o, n)))


PART_NAMES_G = ['multiply_forward99', 'norm_xor_noscratch', 'mul_xor_noscratch', 'produce_retained', 'rebuild_retained',
                'inverse6_xor20', 'phase_inverse6_classical6', 'phase_product24', 'inverse_stream_noscratch']
wrap_dict(G, PART_NAMES_G)
wrap_obj(ar, ['mul_xor0', 'phase_product', 'erase_product', '_erase_norm', 'erase_inverse'])
wrap_dict(vars(mod), ['zero_flag'])
for mn, names in (('exp.arithmetic.division', ['zero_flag']), ('exp.jet24.pro_jet', ['cached_inverse'])):
    try:
        m = __import__(mn, fromlist=['x'])
        wrap_obj(m, names)
    except Exception as e:                                               # noqa: BLE001
        print('wrap skip', mn, e)


def kit_variant(flagged, extra):
    """The line's division source (r6n edits + latem's late M) + DQ5's extra edits."""
    text = open(mod.__file__).read()
    start = text.index('def divide_from_norms(')
    src = text[start:text.index('\ndef sweep(', start)].rstrip() + '\n'
    edits = (r6n.FLAGGED if flagged else r6n.UNFLAGGED) + (latem.LATE_M if latem else ()) + tuple(extra)
    fn = r6n._variant(mod, src, edits, 'kern_div_' + ('f' if flagged else 'u'))
    fn.__globals__.update(V.helpers(G, ar, K))

    def run(b, x, y):
        fm = b.alloc_bits(1)
        b.bit_store(fm, False)
        before = lambda F: None                                          # noqa: E731
        fn(b, x, y, fm, before, lambda: None, owner, K)
    return run


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


def circuit_lines(circ):
    t = repr(circ)
    body = t[t.index("'''") + 3:t.rindex("'''")]
    return [ln.strip() for ln in body.split('\n') if ln.strip()]


def simulate(text_lines, nq, x, y, xv, yv):
    c = km.Circuit('\n'.join(text_lines))
    sim = km.Simulator(batch_size=shots, debug_print=False, counting_operations=True, seed=seed + 100)
    sim.write_across_shots(km.array(list(x)), bits(xv, 24))
    sim.write_across_shots(km.array(list(y)), bits(yv, 24))
    sim.do(c)
    return sim, int(sim.operation_statistics().num_toffoli_like) / shots


def run_form(form):
    base, _, var = form.partition('+')
    flagged = base == 'latef'
    spec = V.VARIANTS[var] if var else dict(edits=lambda f: (), hooks=None)
    extra = spec['edits'](flagged)
    saved = G['produce_retained'], G['rebuild_retained']
    if spec.get('hooks'):
        G['produce_retained'], G['rebuild_retained'] = spec['hooks'](G)
        wrap_dict(G, ['produce_retained', 'rebuild_retained'])
    try:
        _run_form(form, flagged, extra)
    finally:
        G['produce_retained'], G['rebuild_retained'] = saved


def _run_form(form, flagged, extra):
    NAMES.clear()
    real = km.CircuitBuilder()
    BUILDER[0] = real
    x = real.alloc_qubits(24)
    y = real.alloc_qubits(24)
    peak = [0]
    base_live = real.num_allocated_qubits

    class P:                                    # live-qubit peak tracker around the real builder
        def __getattr__(self, n):
            return getattr(real, n)

        def alloc_qubits(self, n):
            r = real.alloc_qubits(n)
            peak[0] = max(peak[0], real.num_allocated_qubits - base_live)
            return r
    body = kit_variant(flagged, extra)
    body(P(), x, y)
    BUILDER[0] = None
    circ = real.finish_circuit()
    nq = circ.num_qubits
    lines = circuit_lines(circ)
    rng = np.random.default_rng(seed)
    xv = rng.integers(0, 1 << 24, size=shots, dtype=np.int64)
    yv = rng.integers(1, 1 << 24, size=shots, dtype=np.int64)
    xv[3::17] = 0
    yv[5::19] = 1
    yv[7::23] = 1 << 12
    if flagged:
        yv[::13] = 0
        xv[::13] = 0
    sim, tof = simulate(lines, nq, x, y, xv, yv)
    state = np.asarray(sim.read_across_shots(km.array([km.q(i) for i in range(nq)]))).astype(np.int64)
    flips = np.asarray(sim.read_phase_flipped_across_shots())
    xn = sum(state[q.id] << k for k, q in enumerate(x))
    yn = sum(state[q.id] << k for k, q in enumerate(y))
    okv = (mul(xn, yv) == xv) & (yn == yv)
    if flagged:
        okv = np.where(yv == 0, (xn == 0) & (yn == yv), okv)
    others = [q for q in range(nq) if q not in set(q.id for q in x) | set(q.id for q in y)]
    dirty = state[others].any(axis=0) if others else np.zeros(shots, bool)
    static = sum(1 for ln in lines if ln.split()[0] in ('CCX', 'CCZ', 'CSWAP'))
    print(f'== {form}: executed {tof:.2f} per division; helper peak {peak[0]}; values {int(okv.sum())}/{shots}, '
          f'phase flips {int(flips.sum())}, dirty {int(dirty.sum())}; static Toffoli-like lines {static}', flush=True)
    if not PARTS:
        return
    cuts = [i for i, ln in enumerate(lines) if ln.startswith('INFO_TAG "DQ5PART')]
    prev = 0.0
    agg = {}
    for i, ci in enumerate(cuts):
        _, t = simulate(lines[:ci + 1], nq, x, y, xv, yv)
        name = NAMES[i]
        a = agg.setdefault(name, [0, 0.0])
        a[0] += 1
        a[1] += t - prev
        prev = t
    a = agg.setdefault('(after last part)', [0, 0.0])
    a[1] += tof - prev
    for name, (n, t) in sorted(agg.items(), key=lambda kv: -kv[1][1]):
        print(f'   {name:30s} calls {n:2d}  executed {t:7.2f}', flush=True)


for f in forms:
    try:
        run_form(f)
    except Exception as e:                                               # noqa: BLE001
        import traceback
        traceback.print_exc()
        print('== FAIL', f, e, flush=True)
