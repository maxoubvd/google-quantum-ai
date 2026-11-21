#!/usr/bin/env python3
"""dqi_zn_alias.py -- builder Z1 (sixth holder, Nov 14): R36-DQI-LOC's patch A ("cross-word alias") in the line's
converter. Design: research/panel/R36-DQI-LOC.md, "For the builders" (designer's kill test tools/r36_dqi_loc_alias.py).

In a descent the kit forms the quotient in place (scale q1, cross, scale q0) and then runs
dk.shifted_conv(b, g[:-1], [q0, q1], f[:-2], F) = conv(..., N = m): m evaluation products for m - 1 output words;
the product's head coefficient (index m - 1) is X = q1 * sel(g[-2], g[-3], F), which `cross` has just XORed into the
q0 slot. With N = m - 1 points that head folds into the low words as X * corr[j], corr = x^(m-1) mod prod (x - t);
the repair is free: f[j] ^= corr[j] * (q0 slot) before AND after the cross (net corr[j] * X).
Here: every scope's `cross` is wrapped (alias before and after, the q0 word marked); the kernel module's
`shifted_conv` uses N = len(g) (m - 1) when its q0 word is marked, else the kit's N = len(g) + 1. Levels in
ZN_ALIAS_LEVELS (default 6..16; the root's cross is inline and stays as it is).
Install after the entry and after dqi_zn_pdl2 (ZN_LCT_P=pdl2,kpark,alias or DQI_INT_EXT=...,dqi_zn_alias).
"""
import os
import sys

LEVELS = {int(x) for x in os.environ.get('ZN_ALIAS_LEVELS', ','.join(str(m) for m in range(6, 17))).split(',') if x}
LOG = dict(cross_wraps=0, aliased=0, conv_short=0, conv_full=0, cross_skip=0, conv_unmarked_pending=0)
MARK = {}            # id(q0 word) -> m, set by the wrapped cross, consumed by the next shifted_conv on that word
CORR = {}


def _corr(dk, m):
    """corr[j] = sum_t W[j][t] t^N for N = m - 1 points (x^N mod prod_(t < N) (x - t)), j < m - 1."""
    if m not in CORR:
        f = dk.ar.field(12)
        N = m - 1
        V, W = dk.k.vandermonde(N)
        out = []
        for j in range(m - 1):
            c = 0
            for t in range(N):
                c ^= f.mul(W[j][t], f.pow(t, N))
            out.append(c)
        CORR[m] = out
    return CORR[m]


def _find_f(q0):
    """The level's dividend list f (f[-2] is q0) from the calling frames (cross is called inside the action
    closures of direct_pre / direct_post, where f is a free variable)."""
    fr = sys._getframe(2)
    while fr is not None:
        f = fr.f_locals.get('f')
        if isinstance(f, list) and len(f) >= 3 and f[-2] is q0:
            return f
        fr = fr.f_back
    return None


def make_cross(orig, dkmod):
    def cross(b, g, q1, q0, F):
        f = _find_f(q0)
        m = len(f) - 1 if f is not None else None
        if f is None or m not in LEVELS:
            LOG['cross_skip'] += 1
            return orig(b, g, q1, q0, F)
        if sys._getframe(1).f_locals.get('undo') is True:   # an undo order (conv BEFORE the cross): not patched
            LOG['cross_undo_skip'] = LOG.get('cross_undo_skip', 0) + 1
            return orig(b, g, q1, q0, F)
        if id(q0) in MARK:
            LOG['stale_mark'] = LOG.get('stale_mark', 0) + 1
        corr = _corr(dkmod, m)
        low = f[:-2]
        assert len(low) == len(corr) == m - 1
        for j, w in enumerate(low):
            if corr[j]:
                dkmod.k.cx(b, q0, w, corr[j])
        out = orig(b, g, q1, q0, F)
        for j, w in enumerate(low):
            if corr[j]:
                dkmod.k.cx(b, q0, w, corr[j])
        MARK[id(q0)] = m
        LOG['aliased'] += 1
        return out
    cross._zn_alias = True
    return cross


def make_shifted_conv(orig, dkmod):
    def shifted_conv(b, g, q, out, F, reverse=False):
        m = MARK.pop(id(q[0]), None)
        if m is None:
            LOG['conv_full'] += 1
            return orig(b, g, q, out, F, reverse)
        assert len(g) == m - 1 and len(out) == m - 1, (len(g), len(out), m)
        with dkmod.shifted_pair(b, q, F) as coeff:
            dkmod.conv(b, g, coeff, out, len(g), reverse)
        LOG['conv_short'] += 1
    shifted_conv._zn_alias = True
    return shifted_conv


# ---- patch B (R36-DQI-LOC): the same identity on the low levels (measured_heads' pre and post_look reverse) ----
BLEVELS = {int(x) for x in os.environ.get('ZN_ALIAS_B', '').split(',') if x}     # e.g. 3,4,5


def make_runops(orig, dkmod):
    """The cross op ('m', f[-1], g[-2], f[-2]) of a low level gets the alias before and after it; its q0 word is
    marked so that the level's convolution (inside this ops list or called right after it) runs at N = m - 1."""
    def runops(b, ops, a, ai, undo=False):
        if BLEVELS:
            f = sys._getframe(1).f_locals.get('f')
            if isinstance(f, list) and len(f) >= 4 and len(f) - 1 in BLEVELS:
                m = len(f) - 1
                idx = [i for i, op in enumerate(ops) if op[0] == 'm' and op[1] is f[-1] and op[3] is f[-2]]
                if len(idx) == 1:
                    i = idx[0]
                    q0, low, corr = f[-2], f[:-2], _corr(dkmod, m)
                    MARK[id(q0)] = m
                    LOG['b_aliased'] = LOG.get('b_aliased', 0) + 1
                    segs = [ops[:i], ops[i:i + 1], ops[i + 1:]]
                    if undo:
                        segs = segs[::-1]
                    for k, seg in enumerate(segs):
                        if k == 1:
                            for j, w in enumerate(low):
                                if corr[j]:
                                    dkmod.k.cx(b, q0, w, corr[j])
                        if seg:
                            orig(b, seg, a, ai, undo)
                        if k == 1:
                            for j, w in enumerate(low):
                                if corr[j]:
                                    dkmod.k.cx(b, q0, w, corr[j])
                    return None
        return orig(b, ops, a, ai, undo)
    runops._zn_alias = True
    return runops


def make_convolution(orig, dkmod):
    def convolution_xor(b, x, y, z, reverse=False):
        if len(y) == 2:
            m = MARK.pop(id(y[0]), None)
            if m is not None:
                assert len(x) == len(z) == m - 1, (len(x), len(z), m)
                dkmod.conv(b, x, y, z, len(x), reverse)
                LOG['b_conv_short'] = LOG.get('b_conv_short', 0) + 1
                return None
        return orig(b, x, y, z, reverse)
    convolution_xor._zn_alias = True
    return convolution_xor


# ---- patch A at the root (level 17): the ordinary (non-restore) root_action runs ... x, s q0, c with the cross
# inline; nothing touches q0 or f[:-2] before the x op in either direction, so the first alias sits at the start of
# the call, the second right before the q0 scaling (a wrapper of old.scale), and the conv takes the mark ----
ROOT = os.environ.get('ZN_ALIAS_ROOT') == '1'
PENDING2 = {}


def make_root(orig, dkmod):
    def ordinary_root_action(b, f, g, payloads, ai, Z, F, post, undo):
        m = len(f) - 1
        q0, low = f[-2], f[:-2]
        corr = _corr(dkmod, m)
        for j, w in enumerate(low):
            if corr[j]:
                dkmod.k.cx(b, q0, w, corr[j])
        MARK[id(q0)] = m
        PENDING2[id(q0)] = (low, corr)
        out = orig(b, f, g, payloads, ai, Z, F, post, undo)
        assert id(q0) not in PENDING2 and id(q0) not in MARK, 'root alias: q0 scaling or conv not reached'
        LOG['root_aliased'] = LOG.get('root_aliased', 0) + 1
        return out
    ordinary_root_action._zn_alias = True
    return ordinary_root_action


def make_old_scale(orig, dkmod):
    def scale(b, x, a, ai, divide=True):
        p = PENDING2.pop(id(x), None)
        if p is not None:
            low, corr = p
            for j, w in enumerate(low):
                if corr[j]:
                    dkmod.k.cx(b, x, w, corr[j])
        return orig(b, x, a, ai, divide)
    scale._zn_alias = True
    return scale


# ---- patch C (R36-DQI-LOC): the reverse return at chosen levels in the plain in-place form with the alias, in
# place of a12's handed path (high_full_restore + debts + close_heads): pivot; inverse with the handed known-zero flag
# (a12's inverse_normpair_known_zero); consume; clip; descent conv undone at N = m - 1; q0 unscaled; alias, cross,
# alias; q1 unscaled; the payload's undo; unclip; a08's erase_normpair; pivot. Levels: ZN_ALIAS_C (e.g. 6,7,8,9,11,13,15:
# the levels where Z1's reverse kernels do not hook the handed path). ----
CLEVELS = {int(x) for x in os.environ.get('ZN_ALIAS_C', '').split(',') if x}
CL = os.environ.get('ZN_ALIAS_CL') == '1'          # the paired lower levels too, with pdl2's kernel undo inside
NPARK2 = os.environ.get('ZN_ALIAS_NPARK', '1') == '1'   # at those levels, park the inverse's norm N (6 qubits)


def _alias(b, dkmod, q0, low, corr):
    for j, w in enumerate(low):
        if corr[j]:
            dkmod.k.cx(b, q0, w, corr[j])


def make_handed_c(prev, scope, dkmod):
    def handed(b, f, g, payloads, F, incoming, defer, Z0):
        m = len(f) - 1
        if m not in CLEVELS:
            return prev(b, f, g, payloads, F, incoming, defer, Z0)
        assert not defer
        undo = True                                   # noqa: F841  (read by the cross wrapper: no alias of its own)
        kzi = sys.modules['known_zero_inverse']
        scalar = sys.modules['scalar_compact_v3']
        q0, q1, low = f[-2], f[-1], f[:-2]
        corr = _corr(dkmod, m)
        scope['pivot'](b, g, F)
        a = g[-1]
        pd = sys.modules.get('dqi_zn_pdl2')
        kern = CL and pd is not None and m in pd.LOWER and pd._active(m) and pd.CTX[m].get('undo')
        if kern:
            # Z1's paired kernel, undone in the exact reverse of the forward with LATECONV: its late part (convs, Y,
            # unclip ran after the forward's inverse was erased) first, BEFORE this level's inverse exists (the
            # forward's late width); the in-action part (P, D0 scaling, clip) inside the inverse context below
            B = pd.clip_open(b, m)
            pd.y_parts(b, m, B)
            pd.convs(b, m, True)
        ai, Z, saved = kzi.inverse_normpair_known_zero(b, a, Z0)
        scope['consume'](b, ai, Z, incoming)
        nmask = None
        if kern and NPARK2 and isinstance(saved, tuple) and len(saved) >= 2:
            # the norm N = norm(a) (6 qubits, read only by the erase) parked over the in-context kernel part and the
            # locator undo, recomputed from the unclipped pivot right before the erase (pdl2's NPARK, measured exact)
            N = saved[0]
            nmask = b.alloc_bits(len(N))
            b.hmr(N, nmask)
            b.free(N)
        b.cx(Z, a[0])
        b.cx(Z, ai[0])
        if kern:
            c = pd.CTX[m]
            for v in c['payloads']:
                c['scale'](b, v[-1], a, ai, False)
            pd.p_parts(b, m)
            pd.clip_close(b, m, B)
            c['skip_scales'] = len(c['payloads'])
            c['inv_done'] = True
            LOG['cl_kernels'] = LOG.get('cl_kernels', 0) + 1
        MARK[id(q0)] = m
        scope['dk'].shifted_conv(b, g[:-1], [q0, q1], low, F, True)
        assert id(q0) not in MARK, 'patch C: the conv did not take the mark'
        scope['scale'](b, q0, a, ai, False)
        _alias(b, dkmod, q0, low, corr)
        scope['cross'](b, g, q1, q0, F)
        _alias(b, dkmod, q0, low, corr)
        scope['scale'](b, q1, a, ai, False)
        for v in payloads:
            scope['dk'].shifted_product(b, g[:-1], v[-1], v[:-1], F)
        for v in payloads:
            scope['scale'](b, v[-1], a, ai, False)
        b.cx(Z, ai[0])
        b.cx(Z, a[0])
        if nmask is not None:
            to, frm, K = scalar.ar.tower(12)
            scalar.ar.linear_inplace(b, a, to)
            N2 = b.alloc_qubits(len(nmask))
            scalar.norm_xor_noscratch(b, a, N2, K)
            scalar.ar.linear_inplace(b, a, frm)
            b.cz(N2, nmask)
            b.free(nmask)
            saved = (N2,) + tuple(saved[1:])
            LOG['c_npark'] = LOG.get('c_npark', 0) + 1
        scalar.erase_normpair(b, a, ai, Z, saved)
        scope['pivot'](b, g, F)
        LOG['c_levels'] = LOG.get('c_levels', 0) + 1
    handed._zn_alias = True
    return handed


def install_c(fns, dkm):
    """Every reachable function with a free variable `handed` (the line's direct_pre, possibly inside other seats'
    wrappers' closures): its cell gets the patch-C dispatcher; the helpers come from that function's globals (the
    converter scope: pivot, dk, scale, cross, consume)."""
    n = 0
    for fn in fns:
        co = getattr(fn, '__code__', None)
        if co is None or 'handed' not in co.co_freevars or not fn.__closure__:
            continue
        cells = dict(zip(co.co_freevars, fn.__closure__))
        try:
            cur = cells['handed'].cell_contents
        except ValueError:
            continue
        if getattr(cur, '_zn_alias', False) or not callable(cur):
            continue
        scope = fn.__globals__
        hco = getattr(cur, '__code__', None)
        if hco is not None and cur.__closure__ and 'scope' in hco.co_freevars:    # a12's handed: scope is a cell
            sc = dict(zip(hco.co_freevars, cur.__closure__))['scope'].cell_contents
            if isinstance(sc, dict):
                scope = sc
        if not all(k in scope for k in ('pivot', 'dk', 'scale', 'cross', 'consume')):
            LOG['c_missing'] = sorted(k for k in ('pivot', 'dk', 'scale', 'cross', 'consume') if k not in scope)
            LOG['c_noscope'] = LOG.get('c_noscope', 0) + 1
            continue
        cells['handed'].cell_contents = make_handed_c(cur, scope, dkm)
        n += 1
    LOG['c_wraps'] = n


def install(km, ns):
    walk = sys.modules['walk_bound']
    dicts = []
    for fn in walk.functions(ns.decode):
        dicts.append(fn.__globals__)
        for cell in (fn.__closure__ or ()):
            try:
                v = cell.cell_contents
            except ValueError:
                continue
            if isinstance(v, dict) and 'cross' in v:
                dicts.append(v)
    dkmods = {}
    done = set()
    for gl in dicts:
        if id(gl) in done:
            continue
        done.add(id(gl))
        dk = gl.get('dk')
        dkm = getattr(dk, '_dk', dk)                     # pdl2's DKProxy keeps the module in _dk
        if type(dkm).__name__ != 'module' or not dkm.__name__.endswith('direct_kernels'):
            continue
        dkmods[id(dkm)] = dkm
        cr = gl.get('cross')
        if callable(cr) and not getattr(cr, '_zn_alias', False) and hasattr(cr, '__code__') \
                and cr.__code__.co_varnames[:5] == ('b', 'g', 'q1', 'q0', 'F'):
            gl['cross'] = make_cross(cr, dkm)
            LOG['cross_wraps'] += 1
        if BLEVELS:
            ro = gl.get('runops')
            if callable(ro) and not getattr(ro, '_zn_alias', False):
                gl['runops'] = make_runops(ro, dkm)
                LOG['runops_wraps'] = LOG.get('runops_wraps', 0) + 1
            cx = gl.get('convolution_xor')
            if callable(cx) and not getattr(cx, '_zn_alias', False):
                gl['convolution_xor'] = make_convolution(cx, dkm)
                LOG['conv_wraps'] = LOG.get('conv_wraps', 0) + 1
    for dkm in dkmods.values():
        if not getattr(dkm.shifted_conv, '_zn_alias', False):
            dkm.shifted_conv = make_shifted_conv(dkm.shifted_conv, dkm)
    if ROOT and dkmods:
        dkm = list(dkmods.values())[0]
        seen = set()
        for fn in walk.functions(ns.decode):
            gl = fn.__globals__
            if id(gl) in seen or 'ordinary_root_action' not in gl:
                continue
            seen.add(id(gl))
            ra = gl['ordinary_root_action']
            if callable(ra) and not getattr(ra, '_zn_alias', False):
                gl['ordinary_root_action'] = make_root(ra, dkm)
                LOG['root_wraps'] = LOG.get('root_wraps', 0) + 1
            old = gl.get('old')
            if old is not None and callable(getattr(old, 'scale', None)) \
                    and not getattr(old.scale, '_zn_alias', False):
                old.scale = make_old_scale(old.scale, dkm)
                LOG['old_scale_wraps'] = LOG.get('old_scale_wraps', 0) + 1
    if CLEVELS and dkmods:
        install_c(list(walk.functions(ns.decode)), list(dkmods.values())[0])
    LOG['dk_modules'] = len(dkmods)
    return dict(LOG)
