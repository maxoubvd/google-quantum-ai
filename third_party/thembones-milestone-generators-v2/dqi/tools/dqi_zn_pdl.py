#!/usr/bin/env python3
"""dqi_zn_pdl.py -- builder Z1 (Nov 14): paired payload digits installed into the LINE's converter (r12hz), in the
direct region (blocks (11, 10) and (9, 8): the levels where the line runs narrow240's direct converter).

Milestone 1 (this file): the DESCENT is paired, the line's return is unchanged. After the descent the digits are the
line's own (conversion and carries at the bottom), so the converter's output is the line's word for word.
Hooks, all in the scopes whose functions run at the paired levels:
  scope['direct_pre']  wrapper: records each level's (f, g, payloads, F, direction) in CTX
  scope['dk']          proxy: shifted_product at an upper level U = no-op; at a lower level L = the paired kernel
                       (D1 = X q1' into a temp, the shifted digit pair, the N-point product); everything else passes
  scope['scale']       wrapper: the digit scaling at a lower level L also runs the kernel's P part (two free fan-outs
                       of the digit word around s ^= X CG, CG = q1' e2, e2 = the divisor's second top word)
  scope['base_case']   wrapper: conversions Y = s ^= [not B] q0' X and the carries (bottom-up) at the bottom;
                       the reverse runs their inverse at the reverse bottom
q1' is taken early as f_L[-1] * ai (one product, erased by measurement) so the line's locator order is untouched.
Usage: ZN_LCT_P=pdl python dqi_zn_lct.py run ...    (env ZN_PDL_BLOCKS=11,9 = upper levels)
"""
import os
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)

UPPER = tuple(int(x) for x in os.environ.get('ZN_PDL_BLOCKS', '11,9').split(','))
LOWER = tuple(u - 1 for u in UPPER)
LOG = dict(scopes=0, pre_wraps=0, base_wraps=0, kernel_calls=0, scale_hooks=0, bottom=0)
EVENTS = []       # (event, level, detail) in the order they happen (diagnostics)


def ev(*a):
    if len(EVENTS) < 400:
        EVENTS.append(a)
CTX = {}          # level -> dict(f, g, payloads, F, undo); reset at every converter build (top level 17 forward)


class DKProxy:
    """The line's direct_kernels module with shifted_product dispatched by level."""

    def __init__(self, dk):
        self._dk = dk

    def __getattr__(self, name):
        return getattr(self._dk, name)

    def shifted_product(self, b, g, x, out, F):
        m = len(g) + 1
        if m in UPPER and m in CTX:
            return None                          # upper level of a paired block: no rank updates
        if m in LOWER and m in CTX:
            LOG['kernel_calls'] += 1
            return lower_conv(b, m, x, out, F)
        return self._dk.shifted_product(b, g, x, out, F)


def _payload_index(m, x):
    for k, v in enumerate(CTX[m]['payloads']):
        if v[-1] is x:
            return k
    return None


def _ar():
    from exp.safe_scalar.converter import ar
    return ar


def lower_setup(b, m, ai):
    """Q1t = f_L[-1] ai (= q1'), B = [Q1t == 0], Q1t clipped; e2 = the divisor's effective second-top word moved
    into g[-2] by one conditional swap; CG = Q1t e2."""
    from exp.safe_scalar import converter as SS
    ar = _ar()
    c = CTX[m]
    f, g, F = c['f'], c['g'], c['F']
    Q1t = b.alloc_qubits(12)
    ar.mul_xor0(b, f[-1], ai, Q1t)
    B = SS.fresh_zero_flag(b, Q1t)
    b.cx(B, Q1t[0])
    c.update(Q1t=Q1t, B=B, ai=ai)


def lower_cleanup(b, m):
    from exp.safe_scalar import converter as SS
    ar = _ar()
    c = CTX[m]
    b.cx(c['B'], c['Q1t'][0])
    SS.clear_flag(b, c['Q1t'], c['B'])
    ar.erase_product(b, c['f'][-1], c['ai'], c['Q1t'])
    b.free(c['Q1t'])
    for k in ('Q1t', 'B', 'ai'):
        c.pop(k)


def cg_open(b, m):
    ar = _ar()
    c = CTX[m]
    g, F = c['g'], c['F']
    b.cswap(F, g[-2], g[-3])
    CG = b.alloc_qubits(12)
    ar.mul_xor0(b, c['Q1t'], g[-2], CG)
    c['CG'] = CG


def cg_close(b, m):
    ar = _ar()
    c = CTX[m]
    g, F = c['g'], c['F']
    ar.erase_product(b, c['Q1t'], g[-2], c['CG'])
    b.free(c['CG'])
    c.pop('CG')
    b.cswap(F, g[-2], g[-3])


def p_part(b, m, k):
    import dqi_zn_pdk as PK
    ar = _ar()
    c = CTX[m]
    v = c['payloads'][k]
    t = c['tops'][k]
    n = len(v) - 1
    r = PK.r_poly(n)
    PK.fan_out(b, v[-1], v[:-1], r)
    ar.mul_xor0(b, t, c['CG'], v[-1])
    PK.fan_out(b, v[-1], v[:-1], r)


def lower_conv(b, m, x, out, F):
    """Forward: the point products of the paired step for the payload whose digit word is x (rest words out).
    Reverse (the line's handed pre undoes the payload before the pivot inverse exists): no-op here; the whole
    payload undo runs in the consume hook, right after the inverse is formed."""
    ar = _ar()
    c = CTX[m]
    k = _payload_index(m, x)
    assert k is not None
    if c['undo']:
        ev('conv_rev', m, k)
        c['pending_undo'] = True
        LOG['conv_rev_%d' % m] = LOG.get('conv_rev_%d' % m, 0) + 1
        return None
    t = c['tops'][k]
    D1 = b.alloc_qubits(12)
    ar.mul_xor0(b, t, c['Q1t'], D1)
    with c['dk'].shifted_pair(b, [x, D1], F) as coeff:
        c['dk'].conv(b, c['g'][:-1], coeff, out, len(out))
    ar.erase_product(b, t, c['Q1t'], D1)
    b.free(D1)
    if k == len(c['payloads']) - 1:
        lower_cleanup(b, m)


def lower_undo(b, m, ai, kit_scale):
    """The exact inverse of the forward paired step at lower level m (inside the pivot inverse context, clipped)."""
    ar = _ar()
    c = CTX[m]
    a = c['g'][-1]
    lower_setup(b, m, ai)
    for k in reversed(range(len(c['payloads']))):
        v = c['payloads'][k]
        t = c['tops'][k]
        D1 = b.alloc_qubits(12)
        ar.mul_xor0(b, t, c['Q1t'], D1)
        with c['dk'].shifted_pair(b, [v[-1], D1], c['F']) as coeff:
            c['dk'].conv(b, c['g'][:-1], coeff, v[:-1], len(v) - 1, True)
        ar.erase_product(b, t, c['Q1t'], D1)
        b.free(D1)
    for k in reversed(range(len(c['payloads']))):
        kit_scale(b, c['payloads'][k][-1], a, ai, False)
    cg_open(b, m)
    for k in reversed(range(len(c['payloads']))):
        p_part(b, m, k)
    cg_close(b, m)
    lower_cleanup(b, m)


def _find_ai():
    f = sys._getframe(2)
    while f is not None:
        if 'ai' in f.f_locals:
            return f.f_locals['ai']
        f = f.f_back
    raise AssertionError('no ai in the calling frames')


def make_scale(orig):
    def scale(b, x, a, ai, divide=True):
        for m in LOWER:
            if m in CTX and CTX[m].get('active', 0) > 0:
                k = _payload_index(m, x)
                if k is not None:
                    LOG['scale_hooks'] += 1
                    c = CTX[m]
                    key = 'scale_%s_%d' % ('rev' if c['undo'] else 'fwd', m)
                    ev('scale', m, k, c['undo'], c.get('pending_undo', False), c.get('skip_scales', 0))
                    LOG[key] = LOG.get(key, 0) + 1
                    if not c['undo']:
                        if k == 0:                         # both P parts first: CG never meets the scaling temp
                            lower_setup(b, m, ai)
                            cg_open(b, m)
                            for kk in range(len(c['payloads'])):
                                p_part(b, m, kk)
                            cg_close(b, m)
                        orig(b, x, a, ai, divide)
                    else:                                  # reverse through the plain unscale (ai exists)
                        if c.pop('pending_undo', False):
                            lower_undo(b, m, ai, orig)
                            c['skip_scales'] = len(c['payloads'])
                        if c.get('skip_scales', 0) > 0:
                            c['skip_scales'] -= 1
                            return None
                        return orig(b, x, a, ai, divide)
                    return None
        return orig(b, x, a, ai, divide)
    scale._zn_pdl = True
    return scale


def make_direct_pre(orig, scope):
    def direct_pre(b, f, g, payloads, F, undo=False, incoming=None, defer=False, **kw):
        m = len(f) - 1
        if m == 17 and not undo:
            CTX.clear()
        mine = m in UPPER or m in LOWER
        if mine:
            ev('pre_enter', m, undo, kw.get('handoff_zero') is not None)
            c = CTX.setdefault(m, {})
            if not c.get('active'):
                c.update(f=f, g=g, payloads=payloads, F=F, undo=undo,
                         dk=scope['dk']._dk if isinstance(scope['dk'], DKProxy) else scope['dk'])
                if m in UPPER:
                    CTX.setdefault(m - 1, {})['tops'] = [v[-1] for v in payloads]
            c['active'] = c.get('active', 0) + 1
        try:
            return orig(b, f, g, payloads, F, undo, incoming, defer, **kw)
        finally:
            if mine:
                CTX[m]['active'] -= 1
                if CTX[m]['active'] == 0:
                    assert not CTX[m].get('pending_undo'), 'paired undo never ran at level %d' % m
                    CTX[m].pop('skip_scales', None)
    direct_pre._zn_pdl = True
    return direct_pre


REC = {}          # level -> dict(f, payloads): every level 1..17 (the bottom needs levels down to L - 3)


def make_recorder(orig, name):
    """Wrap a level function (pre / direct_pre of any variant) to record (f, payloads) per level."""
    def rec(b, f, g, payloads, *a, **kw):
        REC[len(f) - 1] = dict(f=f, payloads=payloads)
        return orig(b, f, g, payloads, *a, **kw)
    rec._zn_pdl = True
    rec.__name__ = name
    return rec


def tops(m):
    return [v[-1] for v in REC[m]['payloads']]


def route_block(b, U, transposed=False):
    """Carry of block (U, U-1) at the bottom (the pd3 routes of tools/dqi_zn_pdc.py, read with the line's level
    records): flags = zero tests of the quotient words q1 of levels U-1 .. U-4, held while the block is routed."""
    import dqi_zn_pdc as PC
    from exp.safe_scalar import converter as SS
    ar = _ar()
    q1 = {m: REC[m]['f'][-1] for m in (U - 1, U - 2, U - 3, U - 4)}
    fl = {m: SS.fresh_zero_flag(b, q1[m]) for m in sorted(q1)}
    Bq, Cq, Fq, C2q = fl[U - 1], fl[U - 2], fl[U - 3], fl[U - 4]
    q0u = REC[U - 2]['f'][-2]
    h1, l1 = PC.flag_and(b, [(Bq, True), (Cq, False)])
    ZB = PC.masked_word(b, h1, q0u)
    h2, l2 = PC.flag_and(b, [(Bq, True), (Cq, True), (Fq, True)])
    b.cx(h2, ZB[0])
    sA, lA = PC.flag_and(b, [(Bq, True), (Cq, False), (C2q, True)])
    sB, lB = PC.flag_and(b, [(Bq, True), (Cq, True), (Fq, False), (C2q, True)])
    sU2, lU2 = PC.flag_and(b, [(Bq, True), (Cq, False), (C2q, False), (Fq, True)])
    X, tu, tl, t2 = tops(U), tops(U - 2), tops(U - 3), tops(U - 4)
    for k in range(len(X)):
        if not transposed:
            ar.mul_xor0(b, ZB, X[k], tu[k])
            for sel in (sA, sB):
                for j in range(12):
                    b.ccx(sel, X[k][j], tl[k][j])
            for j in range(12):
                b.ccx(sU2, X[k][j], t2[k][j])
        else:
            ar.mul_xor0(b, ZB, tu[k], X[k])
            for sel in (sA, sB):
                for j in range(12):
                    b.ccx(sel, tl[k][j], X[k][j])
            for j in range(12):
                b.ccx(sU2, t2[k][j], X[k][j])
    PC.erase_flag_and(b, sU2, lU2)
    PC.erase_flag_and(b, sB, lB)
    PC.erase_flag_and(b, sA, lA)
    b.cx(h2, ZB[0])
    PC.erase_flag_and(b, h2, l2)
    PC.erase_masked_word(b, h1, q0u, ZB)
    PC.erase_flag_and(b, h1, l1)
    for m in sorted(fl, reverse=True):
        SS.clear_flag(b, q1[m], fl[m])


def y_block(b, U):
    """Conversion Y of block (U, L): s_L ^= [not B] q0' X, B = [q1' == 0] (q1', q0' = level L's quotient words)."""
    import dqi_zn_pdc as PC
    from exp.safe_scalar import converter as SS
    ar = _ar()
    L = U - 1
    fL = REC[L]['f']
    B = SS.fresh_zero_flag(b, fL[-1])
    Z0 = PC.masked_word(b, B, fL[-2], neg=True)
    for x, s in zip(tops(U), tops(L)):
        ar.mul_xor0(b, Z0, x, s)
    PC.erase_masked_word(b, B, fL[-2], Z0, neg=True)
    SS.clear_flag(b, fL[-1], B)


def bottom(b, reverse):
    LOG['bottom'] += 1
    ups = sorted(UPPER)
    if not reverse:
        for U in ups:
            y_block(b, U)
        for U in ups:                                     # carries bottom-up
            route_block(b, U)
    else:
        for U in reversed(ups):
            route_block(b, U)
        for U in ups:
            y_block(b, U)


def make_base(orig):
    def base_case(b, f, g, payloads, H, parent, reverse=False, phase=None):
        REC[1] = dict(f=f, payloads=payloads)
        outer = not CTX.get('_in_base')
        CTX['_in_base'] = True
        try:
            if not reverse:
                if outer:
                    bottom(b, False)
                return orig(b, f, g, payloads, H, parent, reverse, phase)
            out = orig(b, f, g, payloads, H, parent, reverse, phase)
            if outer:
                bottom(b, True)
            return out
        finally:
            if outer:
                CTX['_in_base'] = False
    base_case._zn_pdl = True
    return base_case


def install(km, ns):
    from exp.safe_scalar import converter as SSmod
    walk = sys.modules['walk_bound']
    done = set()
    fns = list(walk.functions(ns.decode))
    dicts = []
    for fn in fns:
        dicts.append(fn.__globals__)
        for cell in (fn.__closure__ or ()):          # scope dicts captured in closures (private rebinds)
            try:
                v = cell.cell_contents
            except ValueError:
                continue
            if isinstance(v, dict) and any(k in v for k in ('consume', 'dk', 'scale', 'direct_pre', 'base_case')):
                dicts.append(v)
    for gl in dicts:
        if id(gl) in done:
            continue
        done.add(id(gl))
        hit = False
        dkm = gl.get('dk')
        if type(dkm).__name__ == 'module' and dkm.__name__.endswith('direct_kernels'):
            gl['dk'] = DKProxy(gl['dk'])
            hit = True
        sc = gl.get('scale')
        co = getattr(sc, '__code__', None)
        if (callable(sc) and co is not None and not getattr(sc, '_zn_pdl', False)
                and co.co_varnames[:5] == ('b', 'x', 'a', 'ai', 'divide')):
            gl['scale'] = make_scale(sc)
            LOG.setdefault('scale_files', set()).add('%s:%d' % (os.path.basename(co.co_filename), co.co_firstlineno))
            hit = True
        for name in ('direct_pre', 'pre'):
            f0 = gl.get(name)
            if callable(f0) and not getattr(f0, '_zn_pdl', False) and hasattr(f0, '__code__'):
                w = make_recorder(f0, name)
                if name == 'direct_pre':
                    w = make_direct_pre(w, gl)
                gl[name] = w
                LOG['pre_wraps'] += 1
        cs = gl.get('consume')
        if callable(cs) and not getattr(cs, '_zn_pdl', False) and getattr(cs, '__name__', '') == 'consume':
            gl['consume'] = make_consume(cs)
            LOG['consume_wraps'] = LOG.get('consume_wraps', 0) + 1
        bc = gl.get('base_case')
        if callable(bc) and not getattr(bc, '_zn_pdl', False) and hasattr(bc, '__code__'):
            gl['base_case'] = make_base(bc)
            LOG['base_wraps'] += 1
        LOG['scopes'] += int(hit)
    sc3 = sys.modules.get('scalar_compact_v3')
    if sc3 is not None and not getattr(sc3.scale_before_inverse, '_zn_pdl', False):
        sc3.scale_before_inverse = make_sbi(sc3.scale_before_inverse)
        LOG['sbi'] = 1
    out = dict(LOG)
    out['scale_files'] = sorted(LOG.get('scale_files', ()))
    return out


def _active_lower_undo():
    for m in LOWER:
        c = CTX.get(m)
        if c and c.get('active', 0) > 0 and c.get('undo'):
            return m
    return None


def make_sbi(orig):
    """scale_before_inverse at a paired lower level in the reverse: deferred to the consume hook (zero mask)."""
    def scale_before_inverse(b, x, a):
        m = _active_lower_undo()
        if m is not None and _payload_index(m, x) is not None:
            LOG['sbi_%d' % m] = LOG.get('sbi_%d' % m, 0) + 1
            mask = b.alloc_bits(12)
            b.bit_store(mask, False)
            return mask
        return orig(b, x, a)
    scale_before_inverse._zn_pdl = True
    return scale_before_inverse


def make_consume(orig):
    def consume(b, ai, Z, incoming):
        out = orig(b, ai, Z, incoming)
        m = _active_lower_undo()
        if m is not None:
            ev('consume', m, CTX[m].get('pending_undo', False))
        if m is not None and CTX[m].pop('pending_undo', False):
            from exp.safe_scalar import converter as SS
            a = CTX[m]['g'][-1]
            b.cx(Z, a[0])
            b.cx(Z, ai[0])
            lower_undo(b, m, ai, SS.scale)
            b.cx(Z, ai[0])
            b.cx(Z, a[0])
            CTX[m]['skip_scales'] = len(CTX[m]['payloads'])  # the line's own unscales that may follow are done
            LOG['undo_runs'] = LOG.get('undo_runs', 0) + 1
        return out
    consume._zn_pdl = True
    return consume
