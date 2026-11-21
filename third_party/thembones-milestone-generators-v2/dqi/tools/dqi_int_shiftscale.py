"""dqi_int_shiftscale.py -- workstream dqi_int (B5): the pivot word joins D1's shifted pair at the reverse aligned levels.

a12's reverse `handed` (research/agents/a12/known_zero_inverse.py, aligned levels 6..16) runs, per level,
    for v in payloads: dk.shifted_product(b, g[:-1], v[-1], v[:-1], F)        # D1's pub_pair pairs the two calls
    b.cx(Z0, a[0]); debts = [scale_before_inverse(b, v[-1], a) for v in payloads]; b.cx(Z0, a[0])
scale_before_inverse is out = head * a' (one product into a fresh word), the old head X-measured (debt), swap. Both
heads are the SAME two scalars as in the shifted pair, and a' = a ^ Z0 (bit 0) is untouched by the shifted products
(they read g[:-1] only), so a' is one more word of the same pair program: 67 more CCX instead of 2 x 42 (-17 per level).
Width: the points of D1's pair run first (unchanged), then the two fresh words are allocated, the pivot word is emitted,
and the one shared cache (D1 emits it first; XORs commute and the scalars do not change) last: the two fresh words are
live together where the original had one (+12 wires at that moment).
Install AFTER pub_pair (DQI_PUB_PAIR_SHIFT=1): token `shiftscale` in DQI_INT_P. DQI_INT_SS_OFF=1: nothing (control).
a12's handed is patched by replacing its __code__ (every reference sees it), its two lines rewritten.
"""
import inspect
import os
import sys
import textwrap

OLD = ("        for v in payloads:scope['dk'].shifted_product(b,g[:-1],v[-1],v[:-1],F)\n"
       "        b.cx(Z0,a[0]);debts=[scalar.scale_before_inverse(b,v[-1],a)for v in payloads];b.cx(Z0,a[0])\n")
NEW = "        b.cx(Z0,a[0]);debts=_B5_SHIFTSCALE(b,g[:-1],payloads,F,a,scope,scalar);b.cx(Z0,a[0])\n"
LOG = dict(merged=0, fallback=0)


def _pair_with_pivot(pp, b, gw, x1, out1, x2, out2, F, a, o1, o2):
    """D1's shifted_pair (tools/dqi_pub_pair.py, read-only: its kernels are called, its order changed): the points first,
    then the pivot word a into (o1, o2) = (a x1, a x2), then the shared cache over every instance."""
    prog, K = pp.PAIRSP['prog'], pp.PAIRSP['K']
    N = len(gw)
    V, W = K.vandermonde(N)
    gf_mul = K.mul
    pts = [t for t in range(N) if any(W[j][t] for j in range(max(len(out1), len(out2))))]

    def tset_for(t):
        ts = []
        for out in (out1, out2):
            for e in range(12):
                qs = []
                for j in range(len(out)):
                    c = W[j][t]
                    if c:
                        v = gf_mul(c, 1 << e)
                        qs.extend(out[j][i] for i in range(12) if (v >> i) & 1)
                ts.append(pp._xorset([qs]))
        return ts
    tsets = {t: tset_for(t) for t in pts}
    z0 = b.alloc_qubits(12)           # D1's allocation parity with the two kit calls
    m0 = b.alloc_bits(12)
    b.free(z0)
    b.free(m0)
    z = b.alloc_qubits(12)
    mask = b.alloc_bits(12)
    for t in pts:
        av, c = V[t][0], (V[t][1] if N > 1 else 0)
        xg, p = K.expose(b, gw, V[t][:N])
        s = av ^ c
        if s == 0:
            if av != 1:
                K.scale_const(b, xg, av)
            pp.emit_word_map(b, prog, xg, x1, x2, tsets[t])
            if av != 1:
                K.scale_const(b, xg, K.inv(av))
        else:
            K.scale_const(b, xg, s)
            for j in range(12):
                b.init_and(F, xg[j], z[j])
            K.scale_const(b, xg, K.inv(s))
            if av:
                if av != 1:
                    K.scale_const(b, xg, av)
                b.cx(xg, z)
                if av != 1:
                    K.scale_const(b, xg, K.inv(av))
            pp.emit_word_map(b, prog, z, x1, x2, tsets[t])
            b.hmr(z, mask)
            if av:
                if av != 1:
                    K.scale_const(b, xg, av)
                b.cz(xg, mask)
                if av != 1:
                    K.scale_const(b, xg, K.inv(av))
            K.scale_const(b, xg, s)
            for j in range(12):
                b.ccz(F, xg[j], mask[j])
            K.scale_const(b, xg, K.inv(s))
        K.unexpose(b, gw, V[t][:N], p)
    b.free(z)
    b.free(mask)
    o1[:] = b.alloc_qubits(12)
    o2[:] = b.alloc_qubits(12)
    tp = pp._plain_tset(o1, o2)
    pp.emit_word_map(b, prog, a, x1, x2, tp)
    pp.emit_cache_map(b, prog, x1, x2, [tsets[t] for t in pts] + [tp])


MAXM = int(os.environ.get('DQI_INT_SS_MAXM', '11'))     # a14's restored frames hold levels 17..12: never there
# (Nov 14 16:10Z) DQI_INT_SS_LEVELS=6,7,9,11: an explicit level set in place of m <= MAXM (unset: as before)
SSL = {int(x) for x in os.environ.get('DQI_INT_SS_LEVELS', '').split(',') if x}


def _shiftscale(b, gw, payloads, F, a, scope, scalar):
    pp = sys.modules.get('dqi_pub_pair')
    m = len(gw) + 1                                     # handed's level: g has m words, gw = g[:-1]
    inlv = (m in SSL) if SSL else (m <= MAXM)
    ok = (pp is not None and pp.PAIRSP.get('prog') is not None and pp.PAIRSP.get('pending') is None
          and len(payloads) == 2 and inlv and os.environ.get('DQI_INT_SS_OFF', '0') != '1')
    if not ok:
        LOG['fallback'] += 1
        for v in payloads:
            scope['dk'].shifted_product(b, gw, v[-1], v[:-1], F)
        return [scalar.scale_before_inverse(b, v[-1], a) for v in payloads]
    v1, v2 = payloads
    o1, o2 = [], []
    _pair_with_pivot(pp, b, list(gw), v1[-1], v1[:-1], v2[-1], v2[:-1], F, a, o1, o2)
    debts = []
    for v, o in ((v1, o1), (v2, o2)):
        m = b.alloc_bits(12)
        b.hmr(v[-1], m)
        b.swap(v[-1], o)
        debts.append(m)
    b.free(o2)
    b.free(o1)
    pp.PAIRSP['pairs'] += 1
    LOG['merged'] += 1
    LOG.setdefault('levels', []).append(m)
    return debts


def install(km, ns):
    if os.environ.get('DQI_INT_SS_OFF', '0') == '1':
        return {'off': True}
    walk = sys.modules['walk_bound']
    done = 0
    seen = set()
    for f in walk.functions(ns.decode):
        if id(f) in seen or getattr(f, '__name__', '') != 'handed':
            continue
        if not f.__code__.co_filename.endswith('known_zero_inverse.py'):
            continue
        seen.add(id(f))
        raw = inspect.getsource(f)
        assert raw.count(OLD) == 1, 'shiftscale: handed lines not found exactly once'
        body = textwrap.dedent(raw.replace(OLD, NEW))
        fv = f.__code__.co_freevars
        outer = 'def _b5_outer():\n    ' + ' = '.join(fv) + ' = None\n' + textwrap.indent(body, '    ') + \
                '    return handed\n'
        nsx = {}
        exec(compile(outer, f.__code__.co_filename, 'exec'), {}, nsx)
        code = [c for c in nsx['_b5_outer'].__code__.co_consts if hasattr(c, 'co_name') and c.co_name == 'handed'][0]
        assert code.co_freevars == fv, (code.co_freevars, fv)
        f.__globals__['_B5_SHIFTSCALE'] = _shiftscale
        f.__code__ = code
        done += 1
    assert done, 'shiftscale: no a12 handed found'
    return dict(patched=done)
