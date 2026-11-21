#!/usr/bin/env python3
"""dqi_dq3_npair.py -- DQ3 (Nov 18): ROW 'npair' of the census: the converter's norm-pair 12-bit inverse in the kit's
CACHED form.

The line's unconvert runs research/agents/a12/known_zero_inverse.inverse_normpair_known_zero(b, a, Z) 15 times a
decode (entry 866-875 live, +24) [M census cen3]: N = Norm(a) (15), t = 1/N by the known-zero GF(64) leaf (17),
ai = conj(a) t by two 6-bit products (30): 62 Toffolis; it returns (ai, Z, (N, t)) and the pair (N, t) is kept for
the later erase_normpair.  This module computes the same three values with the kit's cached 12-bit inverse
(exp/jet24/pro_jet.cached_inverse with PLAN12: the norm's 15 products kept as raw wires, the same known-zero leaf on
D = Norm(a), a 25-product final stage that reads the kept raw wires) with two changes: D is COPIED (6 CNOTs) before
the cached form folds it into the output, and k = 1/D is KEPT instead of erased: 15 + 17 + 25 = 57 [P], peak +24
(raw 10 + k 6 + the copy 6 after the leaf; out 12 + k 6 + copy 6 at the products), the original's width.
Exactness: `python dqi_dq3_npair.py test` compares ai, N, t with the line's kernel on ALL 4,096 inputs a (Z = [a == 0])
with B5's bit-parallel trace builder, and prices both forms with kickmix's own statistics.
Install: token dqi_dq3_npair LAST in DQI_INT_EXT.  Switches: DQI_DQ3_NPAIR_CEIL (901), DQI_DQ3_NPAIR_EXTRA (24),
DQI_DQ3_NPAIR_OFF=1, DQI_DQ3_NPAIR_OUT=<json> (per call: live at entry, form).
"""
import atexit
import gc
import json
import os
import sys
import types

CEIL = int(os.environ.get('DQI_DQ3_NPAIR_CEIL', '901'))
EXTRA = int(os.environ.get('DQI_DQ3_NPAIR_EXTRA', '24'))
STATS = dict(cached=0, orig=0, patched=0)
CALLS = []
sys.modules.setdefault('dqi_dq3_npair', sys.modules[__name__])


def cached(b, a, Z, G):
    """(ai, Z, (N, t)) as inverse_normpair_known_zero, by the cached form.  G: the a12 module's globals."""
    known_leaf = G['inverse6_known_zero']
    out, N, k = cached_pair(b, a, G['jet'], lambda b_, D_: known_leaf(b_, D_, Z))
    return out, Z, (N, k)


def cached_pair(b, a, jet, leaf):
    """The kit's cached 12-bit inverse (PLAN12) with the GF(64) leaf `leaf(b, D) -> k` (D preserved), D = Norm(a)
    copied out before the fold and k kept: returns (out = 1/a, N = Norm(a) copy, k = 1/N) in the kit's bases."""
    ref, ar = jet.ref, jet.ar
    plan = jet.PLAN12
    n, h = 12, 6
    to, frm, K = jet.tower(n)
    jet.linear_inplace(b, a, to)
    terms = plan['terms']
    raw = jet.raw_cache(b, a, plan)
    D = raw[:h]
    jet.linear_xor(b, a, D, plan['affine'])            # D = Norm(a)
    k = leaf(b, D)                                      # k = 1/D, D preserved (the known-zero leaf: 17)
    N = b.alloc_qubits(h)
    b.cx(D, N)                                          # the kept norm
    jet.linear_xor(b, a, D, plan['affine'])
    jet.linear_inplace(b, raw, plan['reduce'])
    kept = plan['size'] - len(plan['retire'])
    for q, row in zip(raw[kept:], plan['retire']):
        mask = b.alloc_bits(1)[0]
        b.hmr(q, mask)
        for j in jet.bits(row):
            left, right, _ = terms[j]
            with ref.expose_pair(b, a, left << 1, right << 1) as (p, r):
                b.ccz(p, r, mask)
        b.free(mask)
        b.free(q)
    out = list(raw[:kept]) + list(b.alloc_qubits(n - kept))
    jet.linear_inplace(b, out, plan['out'])
    source = list(a) + list(k)
    for forms, vectors, leaf in plan['leaves']:
        outputs, gates, low, offset = leaf

        def remap(form, forms=forms):
            result = form & 1
            for j in jet.bits(form >> 1):
                result ^= forms[j] << 1
            return result
        for e, v in zip(outputs, vectors):
            if e & 1:
                for j in jet.bits(v):
                    b.x(out[j])
            for a_ in jet.bits((e & (1 << low) - 1) >> 1):
                for i in jet.bits(forms[a_]):
                    for j in jet.bits(v):
                        b.cx(source[i], out[j])
        for i, (left, right) in enumerate(gates):
            vec = 0
            for e, v in zip(outputs, vectors):
                if e >> offset + i & 1:
                    vec ^= v
            with ref.expose_pair(b, source, remap(left), remap(right)) as (p, q):
                jet.rankone_xor(b, p, q, out, vec)
    jet.linear_inplace(b, a, frm)
    jet.linear_inplace(b, out, frm)
    return out, N, list(k)


_STUB = '''
def inverse_normpair_known_zero(b, a, Z):
    _m = __import__('sys').modules['dqi_dq3_npair']
    if b.num_allocated_qubits + _m.EXTRA <= _m.CEIL:
        _m.STATS['cached'] += 1
        _m.CALLS.append((b.num_allocated_qubits, 1))
        return _m.cached(b, a, Z, globals())
    _m.STATS['orig'] += 1
    _m.CALLS.append((b.num_allocated_qubits, 0))
    return _dq3_npair_orig(b, a, Z)
'''


def _target(fn):
    co = getattr(fn, '__code__', None)
    return (isinstance(fn, types.FunctionType) and co is not None and co.co_name == 'inverse_normpair_known_zero'
            and co.co_filename.endswith('a12/known_zero_inverse.py') and 'dqi_dq3_npair' not in co.co_consts)


def _patch(fn):
    gl = fn.__globals__
    if '_dq3_npair_orig' not in gl:
        gl['_dq3_npair_orig'] = types.FunctionType(fn.__code__, gl, 'inverse_normpair_known_zero_line',
                                                   fn.__defaults__, None)
    ns = {}
    pad = '\n' * (fn.__code__.co_firstlineno - 1)
    exec(compile(pad + _STUB, fn.__code__.co_filename, 'exec'), gl, ns)
    fn.__code__ = ns['inverse_normpair_known_zero'].__code__
    STATS['patched'] += 1


def _dump():
    path = os.environ.get('DQI_DQ3_NPAIR_OUT')
    if path:
        try:
            json.dump(dict(STATS, ceil=CEIL, extra=EXTRA, calls=CALLS), open(path, 'w'))
        except Exception:
            pass
    sys.stderr.write('dqi_dq3_npair: %s\n' % json.dumps(dict(STATS, ceil=CEIL, extra=EXTRA)))


def install(km, ns):
    if os.environ.get('DQI_DQ3_NPAIR_OFF', '0') == '1':
        return dict(off=True)
    for obj in gc.get_objects():
        if _target(obj):
            _patch(obj)
    atexit.register(_dump)
    return dict(STATS)


def test(root):
    """Exhaustive value test against the line's kernel + Toffoli / width prices (kickmix statistics)."""
    import importlib.util
    import numpy as np
    sys.path[:0] = [os.path.expanduser('~/qcc/data/kickmix'), root, os.path.join(root, 'playbook/kit/decoder'),
                    os.path.join(root, 'exp/r-main/narrow240'), os.path.join(root, 'research/agents/a08'),
                    os.path.join(root, 'research/agents/a01'), os.path.join(root, 'research/production'),
                    os.path.join(root, 'exp/round7/root'), os.path.join(root, 'exp/next/carry'),
                    os.path.dirname(os.path.abspath(__file__))]
    import kickmix as km
    sys.modules.setdefault('circ', km)
    spec = importlib.util.spec_from_file_location('kz', os.path.join(root, 'research/agents/a12/known_zero_inverse.py'))
    kz = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(kz)
    os.environ['DQI_INT_TRACE_N'] = '4096'
    import dqi_int_inv6trace as T
    T.TB.num_allocated_qubits = property(lambda self: 0)     # TB.free keeps its assertion: no dirty qubit is freed

    def run_tb(fn):
        b = T.TB()
        a = b.alloc_qubits(12)
        for i in range(12):
            a[i].v[:] = (np.arange(4096) >> i) & 1
        Z = b.alloc_qubits(1)[0]
        Z.v[:] = np.arange(4096) == 0
        ai, Z2, (N, t) = fn(b, a, Z)
        val = lambda ws: sum(w.v.astype(np.int64) << i for i, w in enumerate(ws))
        ok_a = all((a[i].v == ((np.arange(4096) >> i) & 1).astype(bool)).all() for i in range(12))
        return val(ai), val(N), val(t), ok_a

    o = run_tb(kz.inverse_normpair_known_zero)
    c = run_tb(lambda b, a, Z: cached(b, a, Z, kz.__dict__))
    print('ai equal on 4096: %s; N equal: %s; t equal: %s; a restored (line, cached): %s %s' % (
        bool((o[0] == c[0]).all()), bool((o[1] == c[1]).all()), bool((o[2] == c[2]).all()), o[3], c[3]))
    gf = km.GF2Field(12)
    want = np.array([0] + [gf.pow(v, 4094) for v in range(1, 4096)])
    print('cached ai == a^-1 (field) on 4096: %s' % bool((c[0] == want).all()))

    def price(fn):
        b = km.CircuitBuilder()
        a = b.alloc_qubits(12)
        Z = b.alloc_qubits(1)[0]
        base = b.num_allocated_qubits
        st = dict(peak=base)
        orig = km.CircuitBuilder.alloc_qubits

        def alloc(self, nq):
            r = orig(self, nq)
            st['peak'] = max(st['peak'], self.num_allocated_qubits)
            return r
        km.CircuitBuilder.alloc_qubits = alloc
        try:
            fn(b, a, Z)
        finally:
            km.CircuitBuilder.alloc_qubits = orig
        s = b.finish_circuit().operation_statistics()
        return s.num_toffoli_like, st['peak'] - base
    print('line kernel: Toffolis %d, peak extra %d' % price(kz.inverse_normpair_known_zero))
    print('cached form: Toffolis %d, peak extra %d' % price(lambda b, a, Z: cached(b, a, Z, kz.__dict__)))


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'test':
    test(sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser(
        '~/qcc/ws/dqi_dq3/gen/zc2p/home/qcc/ws/dqi_zn/src/kv_unfused5'))
