"""dqi_junc_r8p_rinv24b.py -- the kit's reverse 24-bit in-place inverse (playbook/.../inverse12_escape/joint.py,
inverse24_reverse209, the fallback the prefix block's un-extraction always takes) with its 6-bit cache u not held
through the half scalings.

The fallback keeps N = Norm(x) (12), t = 1/N (12), u = Norm6(t) (6, kept only so that N's erasure at the end costs
19 Toffolis: erase_norm_via_inverse_cache19) and a 12-qubit half target together: 42, un-extraction's peak (907).
Here u is erased free as soon as t is built (u = Norm6(t), a quadratic form of t) and rebuilt by one 6-bit norm
(15 Toffolis) right before that erasure: the same gates otherwise, 36 wide, +15 per call (rinv24's form: +163).
Same contract as the fallback. Installed in place of the fallback in every live clone of decode_partial_many.
Usage (test against the fallback): python dqi_junc_r8p_rinv24b.py test
"""
import sys

SRC = '''
def inverse24_reverse_nou(b, x, mask, K=None):
    cols, back, _ = extension_columns(K); linear_inplace(b, x, cols); to, frm, K24 = tower(24); linear_inplace(b, x, to)
    N = b.alloc_qubits(12); norm_xor_noscratch(b, x, N, K24)
    # t = inv(N), as inverse12_keep_inverse_norm64, but u = norm(t) erased at once
    to12, frm12, K12 = tower(12); linear_inplace(b, N, to12)
    D = b.alloc_qubits(6); norm_xor_noscratch(b, N, D, K12); u = inverse6_init19(b, D); ar._erase_norm(b, N, D, K12); b.free(D)
    t = b.alloc_qubits(12); b.cx(N[6:], N[:6]); mul_xor_noscratch(b, N[:6], u, t[:6]); mul_xor_noscratch(b, N[6:], u, t[6:]); b.cx(N[6:], N[:6])
    ar._erase_norm(b, t, u, K12); b.free(u)
    linear_inplace(b, N, frm12); linear_inplace(b, t, frm12)
    out = b.alloc_qubits(12)
    mul_xor_noscratch(b, x[12:], t, out); erase_product(b, out, N, x[12:]); mul_xor_noscratch(b, x[:12], t, x[12:]); b.cx(out, x[12:]); b.cx(out, x[12:]); erase_product(b, x[12:], N, x[:12]); b.cx(out, x[12:]); b.swap(x[:12], x[12:]); b.swap(x[12:], out); b.free(out)
    # u rebuilt for the cheap erasure of N
    linear_inplace(b, t, to12); u = b.alloc_qubits(6); norm_xor_noscratch(b, t, u, K12); linear_inplace(b, t, frm12)
    erase_norm_via_inverse_cache19(b, t, N, u); ar._erase_norm(b, x, t, K24); b.free(t); linear_inplace(b, x, frm); b.cz(x, mask); b.free(mask); linear_inplace(b, x, back)
'''


def make():
    joint = sys.modules['exp.inverse12_escape.joint']
    scope = dict(vars(joint))
    exec(compile(SRC, __file__ + ':nou', 'exec'), scope)
    return scope['inverse24_reverse_nou']


def _targets(module):
    walk = sys.modules['walk_bound']
    out = {}
    for fn in walk.functions(module.decode):
        if fn.__name__ == 'inverse24_reverse_partial' and fn.__code__.co_filename.endswith('/a16/decode_partial_many.py'):
            out[id(fn.__globals__)] = fn.__globals__
    return list(out.values())


def install(km, module):
    fn = make()
    n = 0
    for G in _targets(module):
        fb = G.get('fallback')
        if fb is not None and getattr(fb, '__name__', '') == 'inverse24_reverse209':
            G['fallback'] = fn
            n += 1
    assert n, 'fallback inverse24_reverse209 not found'
    return dict(rinv24b=n)


def _test():
    import os
    import numpy as np
    sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
    import kickmix as km
    sys.modules['circ'] = km
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import dqi_junc_r8_entry as entry                      # noqa: F401
    fallback = sys.modules['exp.inverse12_escape.joint'].inverse24_reverse209
    narrow = make()
    shots = int(os.environ.get('SHOTS', '512'))
    rng = np.random.default_rng(5)
    xv = rng.integers(1, 1 << 24, size=shots, dtype=np.int64)
    mv = rng.integers(0, 1 << 24, size=shots, dtype=np.int64)
    res = {}
    for name, fn in (('fallback', fallback), ('nou', narrow)):
        b = km.CircuitBuilder()
        x = b.alloc_qubits(24)
        mask = b.alloc_bits(24)
        base = b.num_allocated_qubits
        pk = [base]
        orig = km.CircuitBuilder.alloc_qubits

        def alloc(self, n, *a, **kw):
            r = orig(self, n, *a, **kw)
            pk[0] = max(pk[0], self.num_allocated_qubits)
            return r
        km.CircuitBuilder.alloc_qubits = alloc
        try:
            fn(b, x, mask)
        finally:
            km.CircuitBuilder.alloc_qubits = orig
        c = b.finish_circuit()
        sim = km.Simulator(batch_size=shots, debug_print=False, counting_operations=True, seed=3)
        sim.write_across_shots(km.array(list(x)), ((xv[None, :] >> np.arange(24)[:, None]) & 1) == 1)
        sim.write_across_shots(km.array(list(mask)), ((mv[None, :] >> np.arange(24)[:, None]) & 1) == 1)
        sim.do(c)
        out = np.asarray(sim.read_across_shots(km.array(list(x)))).astype(np.int64)
        val = (out * (1 << np.arange(24))[:, None]).sum(axis=0)
        others = km.array([km.q(i) for i in range(c.num_qubits) if i >= 24])
        dirty = np.asarray(sim.read_across_shots(others)).any(axis=0) if c.num_qubits > 24 else np.zeros(shots, bool)
        res[name] = (val, np.asarray(sim.read_phase_flipped_across_shots()), dirty,
                     sim.operation_statistics().num_toffoli_like / shots, pk[0] - base)
    f, n = res['fallback'], res['nou']
    print('fallback: executed %.1f, peak ancillas %d; nou: executed %.1f, peak ancillas %d' % (f[3], f[4], n[3], n[4]))
    print('RESULT values equal %d of %d, phases equal %d of %d, dirty ancillas: fallback %d, nou %d'
          % (int((f[0] == n[0]).sum()), shots, int((f[1] == n[1]).sum()), shots, int(f[2].sum()), int(n[2].sum())))


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'test':
    _test()
