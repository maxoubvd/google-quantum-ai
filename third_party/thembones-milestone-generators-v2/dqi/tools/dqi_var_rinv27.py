"""dqi_var_rinv27.py -- B7 (dqi_var), fourth holder, Nov 14: the line's reverse 24-bit in-place inverse (J's rinv24b
form of the kit's inverse24_reverse209, tools/dqi_junc_r8p_rinv24b.py, read only) with its LAST erasure changed, as
R35-DQI-BLANK-G's re-pricer writes it (schedule S2, step 'p'):

  rinv24b:  N = inv0(t) is erased by rebuilding u = Norm6(t) (15 Toffolis) and the kit's cache erasure
            erase_norm_via_inverse_cache19 (19): 34 per call.
  here:     HMR(N) -> m, then the phase (-1)^(m . inv0(t)) by the kit's 12-bit one-operand inverse phase
            exact_subs.phase_inverse12_virtual (D1's 27-AND kernel once pub_inv27 is installed; the kit's own is
            <= 32): 27 per call, scratch 13 with N already free (25 beyond x at that moment; the kernel's peak
            stays the half scalings' 36).

Everything else is rinv24b's gate list. Same contract as the fallback (x in, 1/x out in place, mask phase).
Token var_rinv27: install after J's rinv24b (DQI_JUNC_R8P) and anywhere among the DQI_INT_P tokens (the phase
function is looked up when the kernel runs, so pub_inv27's patch applies). It replaces the fallback in every live
clone of decode_partial_many's inverse24_reverse_partial.
Usage (values, phases and Toffolis against the fallback and rinv24b, x = 0 included):
  python dqi_var_rinv27.py test        (SHOTS=512; DQI_VAR_RINV27_INV27=1 installs D1's pub_inv27 first)
"""
import os
import sys

SRC = '''
def inverse24_reverse_i27(b, x, mask, K=None):
    cols, back, _ = extension_columns(K); linear_inplace(b, x, cols); to, frm, K24 = tower(24); linear_inplace(b, x, to)
    N = b.alloc_qubits(12); norm_xor_noscratch(b, x, N, K24)
    # t = inv(N), as inverse12_keep_inverse_norm64, but u = norm(t) erased at once (rinv24b)
    to12, frm12, K12 = tower(12); linear_inplace(b, N, to12)
    D = b.alloc_qubits(6); norm_xor_noscratch(b, N, D, K12); u = inverse6_init19(b, D); ar._erase_norm(b, N, D, K12); b.free(D)
    t = b.alloc_qubits(12); b.cx(N[6:], N[:6]); mul_xor_noscratch(b, N[:6], u, t[:6]); mul_xor_noscratch(b, N[6:], u, t[6:]); b.cx(N[6:], N[:6])
    ar._erase_norm(b, t, u, K12); b.free(u)
    linear_inplace(b, N, frm12); linear_inplace(b, t, frm12)
    out = b.alloc_qubits(12)
    mul_xor_noscratch(b, x[12:], t, out); erase_product(b, out, N, x[12:]); mul_xor_noscratch(b, x[:12], t, x[12:]); b.cx(out, x[12:]); b.cx(out, x[12:]); erase_product(b, x[12:], N, x[:12]); b.cx(out, x[12:]); b.swap(x[:12], x[12:]); b.swap(x[12:], out); b.free(out)
    # N = inv0(t) (both in the kit basis): HMR, then (-1)^(m . inv0(t)) by the one-operand inverse phase
    m = b.alloc_bits(12); b.hmr(N, m); b.free(N); _B7_PHASE_INV12()(b, t, m); b.free(m)
    ar._erase_norm(b, x, t, K24); b.free(t); linear_inplace(b, x, frm); b.cz(x, mask); b.free(mask); linear_inplace(b, x, back)
'''
STATS = dict(calls=0)


def _phase_fn():
    """the kit's exact_subs.phase_inverse12_virtual, looked up when the kernel runs (after every install)"""
    mod = sys.modules.get('exact_subs')
    if mod is None:
        import exact_subs as mod                           # noqa: F811  (narrow240 is on the kit's sys.path)
    STATS['calls'] += 1
    return mod.phase_inverse12_virtual


def make():
    joint = sys.modules['exp.inverse12_escape.joint']
    scope = dict(vars(joint))
    scope['_B7_PHASE_INV12'] = _phase_fn
    exec(compile(SRC, __file__ + ':i27', 'exec'), scope)
    return scope['inverse24_reverse_i27']


def _targets(module):
    walk = sys.modules['walk_bound']
    out = {}
    for fn in walk.functions(module.decode):
        if fn.__name__ == 'inverse24_reverse_partial' and fn.__code__.co_filename.endswith('/a16/decode_partial_many.py'):
            out[id(fn.__globals__)] = fn.__globals__
    return list(out.values())


def install(km, module):
    fn = make()
    n, was = 0, []
    for G in _targets(module):
        fb = G.get('fallback')
        nm = getattr(fb, '__name__', '')
        if fb is not None and nm in ('inverse24_reverse_nou', 'inverse24_reverse209'):
            was.append(nm)
            G['fallback'] = fn
            n += 1
    assert n, 'var_rinv27: no reverse 24-bit fallback found'
    import atexit
    atexit.register(lambda: print('dqi_var_rinv27 kernel calls', STATS, flush=True))
    return dict(rinv27=n, replaced=sorted(set(was)))


def _test():
    import numpy as np
    sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
    import kickmix as km
    sys.modules['circ'] = km
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import dqi_junc_r8_entry as entry                      # noqa: F401
    import dqi_junc_r8p_rinv24b as rb
    import exact_subs                                      # noqa: F401  (the kit's; narrow240 is on sys.path)
    if os.environ.get('DQI_VAR_RINV27_INV27', '0') == '1':
        import dqi_pub_inv27 as P
        print('pub_inv27:', P.install(km, None), flush=True)
    fallback = sys.modules['exp.inverse12_escape.joint'].inverse24_reverse209
    shots = int(os.environ.get('SHOTS', '512'))
    rng = np.random.default_rng(7)
    xv = rng.integers(1, 1 << 24, size=shots, dtype=np.int64)
    xv[:8] = 0                                             # x = 0 included (inv0 convention)
    mv = rng.integers(0, 1 << 24, size=shots, dtype=np.int64)
    res = {}
    for name, fn in (('fallback', fallback), ('rinv24b', rb.make()), ('rinv27', make())):
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
    f = res['fallback']
    for name in ('fallback', 'rinv24b', 'rinv27'):
        r = res[name]
        print('%-9s executed %.2f per call, peak ancillas %d; values equal to the fallback %d of %d, phases equal %d '
              'of %d, dirty ancilla shots %d' % (name, r[3], r[4], int((r[0] == f[0]).sum()), shots,
                                                 int((r[1] == f[1]).sum()), shots, int(r[2].sum())))
    print('kernel phase lookups', STATS)


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'test':
    _test()
