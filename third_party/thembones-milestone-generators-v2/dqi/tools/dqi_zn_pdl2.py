#!/usr/bin/env python3
"""dqi_zn_pdl2.py -- builder Z1 (Nov 14): paired payload digits in the LINE's converter, direct region, version 2.
Paired DESCENT (the line's return unchanged); blocks (U, L) = (11, 10), (9, 8) by default (ZN_PDL_BLOCKS = uppers).

Where the kernels sit (chosen so that the level's quotient words q1' = f_L[-1], q0' = f_L[-2] are in place and at
most one temporary word + one flag exist beside the line's own moment: measured room 878 - 854 = 24 there):
  forward, level U:  the line's digit scaling (X = t / a_U) kept; dk.shifted_product -> no-op.
  forward, level L:  the line's digit scalings and dk.shifted_product -> no-ops (ai recorded); after the line's
                     locator step (its dk.shifted_conv) the paired kernel: B = [q1' == 0] (held), q1' clipped in
                     place; CG = q1' e2 (e2 = sel(g[-2], g[-3], F) by one cswap), P = X CG between two free fan-outs
                     of the digit word; D0 = the line's scale; g[:-1] rotated by F (the kit's frame, no shift word);
                     D1 = X q1' (temp), N-point product, D1 erased free; Y = D0 + [not B] q0' X; unclip.
  reverse, level L:  the line's handed pre restores the locator first: the exact inverse runs inside a
                     high_full_restore hook, BEFORE the restore (quotient words still in place), unscaling with the
                     line's own deferred-phase trick; the masks are handed to the line's scale_before_inverse calls,
                     whose debts the line repays after its inverse. (Plain undo path: inverse at the line's
                     shifted_conv undo, with ai.)
  bottom:            the carries (pd3 routes, zero-region terms), bottom-up; reverse: the same XORs top-down.
Usage: ZN_LCT_P=pdl2 python dqi_zn_lct.py run ...
"""
import os
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)

UPPER = tuple(int(x) for x in os.environ.get('ZN_PDL_BLOCKS', '11,9').split(','))
FWDONLY = os.environ.get('ZN_PDL2_FWDONLY') == '1'   # forward converter paired; the reverse is the line's own
LOWER = tuple(u - 1 for u in UPPER)
LOG = dict(scopes=0, pre_wraps=0, base_wraps=0, fwd_kernels=0, inv_sbi=0, inv_ai=0, bottom=0)
EVENTS = []
CTX = {}
REC = {}


def ev(*a):
    if len(EVENTS) < 300:
        EVENTS.append(a)


def _ar():
    from exp.safe_scalar.converter import ar
    return ar


def _ss():
    from exp.safe_scalar import converter as SS
    return SS


def _rotate():
    from repair import rotate
    return rotate


def _payload_index(m, x):
    for k, v in enumerate(CTX[m]['payloads']):
        if v[-1] is x:
            return k
    return None


def _active(m):
    c = CTX.get(m)
    return c is not None and c.get('active', 0) > 0


def _r(n):
    import dqi_zn_pdk as PK
    return PK.r_poly(n)


def _fan(b, s, low):
    import dqi_zn_pdk as PK
    PK.fan_out(b, s, low, _r(len(low)))


def p_parts(b, m):
    """CG = q1' e2 (clipped q1' in place, e2 moved into g[-2] by one cswap); s ^= X CG between two fan-outs."""
    ar = _ar()
    c = CTX[m]
    f, g, F = c['f'], c['g'], c['F']
    b.cswap(F, g[-2], g[-3])
    CG = b.alloc_qubits(12)
    ar.mul_xor0(b, f[-1], g[-2], CG)
    for v, t in zip(c['payloads'], c['tops']):
        _fan(b, v[-1], v[:-1])
        ar.mul_xor0(b, t, CG, v[-1])
        _fan(b, v[-1], v[:-1])
    ar.erase_product(b, f[-1], g[-2], CG)
    b.free(CG)
    b.cswap(F, g[-2], g[-3])


def convs_shift(b, m, undo):
    """convs with the F shift on the DIGIT (the line's shifted_pair trick) instead of rotating g[:-1]: slots
    (X', D1, 0) or (0, X', D1); with g[-2] = 0 whenever F = 1 the polynomial is the same, so the N-point output and
    the P correction are unchanged. 36 CSWAP per payload instead of 24 (L - 2) per block; one more word (z)."""
    ar = _ar()
    SS = _ss()
    c = CTX[m]
    f, g, F = c['f'], c['g'], c['F']
    sp = os.environ.get('ZN_PDL2_SHIFTPARK', '')   # '1' = both directions; 'fwd' = only when not undoing
    park = (sp == '1' or (sp == 'fwd' and not undo)) and c.get('B') is not None
    order = list(zip(c['payloads'], c['tops']))
    for v, t in (reversed(order) if undo else order):
        D1 = b.alloc_qubits(12)
        ar.mul_xor0(b, t, f[-1], D1)
        if park:                                    # B off the wire during the conv (one zero test to rebuild it)
            B = c['B']
            b.cx(B, f[-1][0])                       # unclip: B = [q1' == 0] is again a zero test of f[-1]
            mb = b.alloc_bits(1)
            b.hmr(B, mb[0])
            b.free(B)
        z = b.alloc_qubits(12)
        b.cswap(F, D1, z)
        b.cswap(F, v[-1], D1)
        c['dk'].conv(b, g[:-1], [v[-1], D1, z], v[:-1], len(v) - 1, undo)
        b.cswap(F, v[-1], D1)
        b.cx(z, D1)                                 # z = F D1 now: erased by measurement, phase repaid free
        mask = b.alloc_bits(12)
        b.hmr(z, mask)
        b.free(z)
        for j in range(12):
            b.ccz(F, D1[j], mask[j])
        b.free(mask)
        if park:
            B = SS.fresh_zero_flag(b, f[-1])
            b.cz(B, mb[0])                          # the parking measurement's phase, repaid on the rebuilt flag
            b.free(mb)
            b.cx(B, f[-1][0])                       # clip again
            c['B'] = B
            LOG['b_parks'] = LOG.get('b_parks', 0) + 1
        ar.erase_product(b, t, f[-1], D1)
        b.free(D1)


def _shift_here(m):
    sw = os.environ.get('ZN_PDL2_SHIFT', '')       # '1' = every lower level; else a list of lower levels
    return sw == '1' or str(m) in sw.split(',')


def convs(b, m, undo):
    if _shift_here(m):
        return convs_shift(b, m, undo)
    ar = _ar()
    c = CTX[m]
    f, g, F = c['f'], c['g'], c['F']
    rot = _rotate()
    rot(b, g[:-1], F)
    order = list(zip(c['payloads'], c['tops']))
    for v, t in (reversed(order) if undo else order):
        D1 = b.alloc_qubits(12)
        ar.mul_xor0(b, t, f[-1], D1)
        c['dk'].conv(b, g[:-1], [v[-1], D1], v[:-1], len(v) - 1, undo)
        ar.erase_product(b, t, f[-1], D1)
        b.free(D1)
    rot(b, g[:-1], F, True)


def y_parts(b, m, B):
    import dqi_zn_pdc as PC
    ar = _ar()
    c = CTX[m]
    q0 = c['f'][-2]
    if os.environ.get('ZN_PDL2_YSHARE') == '1':        # fourth holder: one masked word [not B] q0' for both payloads
        Z0 = PC.masked_word(b, B, q0, neg=True)
        for v, t in zip(c['payloads'], c['tops']):
            ar.mul_xor0(b, Z0, t, v[-1])
        PC.erase_masked_word(b, B, q0, Z0, neg=True)
        return
    for v, t in zip(c['payloads'], c['tops']):
        Z0 = PC.masked_word(b, B, q0, neg=True)
        ar.mul_xor0(b, Z0, t, v[-1])
        PC.erase_masked_word(b, B, q0, Z0, neg=True)


def clip_open(b, m):
    SS = _ss()
    q1 = CTX[m]['f'][-1]
    B = SS.fresh_zero_flag(b, q1)
    b.cx(B, q1[0])
    return B


def clip_close(b, m, B):
    SS = _ss()
    q1 = CTX[m]['f'][-1]
    b.cx(B, q1[0])
    SS.clear_flag(b, q1, B)


def lower_forward(b, m):
    """At the end of level L's forward action (after the line's locator step; g pivoted; ai recorded)."""
    LOG['fwd_kernels'] += 1
    c = CTX[m]
    ev('fwd', m)
    B = clip_open(b, m)
    p_parts(b, m)
    for v in c['payloads']:
        c['scale'](b, v[-1], c['a'], c['ai'], True)
    c['B'] = B
    if LATECONV:                                     # fourth holder: the rest after the line has erased ai
        LATE.append(m)
        return
    convs(b, m, False)
    B = c.pop('B')
    y_parts(b, m, B)
    clip_close(b, m, B)


LATECONV = os.environ.get('ZN_PDL2_LATECONV') == '1'
LATE = []


def late_forward(b, m):
    """The paired kernel's convs, Y and unclip, run when the line's inverse_apply has returned (its pivot inverse ai
    erased: 12 qubits fewer than inside the action), before the line's final unpivot (g still pivoted)."""
    c = CTX[m]
    convs(b, m, False)
    B = c.pop('B')
    y_parts(b, m, B)
    clip_close(b, m, B)
    LOG['late'] = LOG.get('late', 0) + 1


def make_inverse_apply(orig):
    def inverse_apply(b, a, action, *rest, **kw):
        n0 = len(LATE)
        out = orig(b, a, action, *rest, **kw)
        while len(LATE) > n0:
            late_forward(b, LATE.pop())
        return out
    inverse_apply._zn_pdl2 = True
    return inverse_apply


def lower_inverse_sbi(b, m):
    """Reverse, handed path: inside high_full_restore's hook, before the restore (g not pivoted; the line has
    flipped bit 0 of g[-2] by the handed zero Z0 for the restore: undone here and redone after)."""
    LOG['inv_sbi'] += 1
    ar = _ar()
    c = CTX[m]
    g, F, Z0 = c['g'], c['F'], c['handoff_zero']
    ev('inv_sbi', m)
    b.cx(Z0, g[-2][0])
    b.cswap(F, g[-1], g[-2])                       # pivot, as in the forward kernel's frame
    a = g[-1]
    B = clip_open(b, m)
    y_parts(b, m, B)
    convs(b, m, True)
    c['masks'] = {}
    b.cx(Z0, a[0])                                 # the line clips a by Z0 around its own deferred unscale
    for k, v in enumerate(c['payloads']):
        s = v[-1]
        out = b.alloc_qubits(12)
        ar.mul_xor0(b, s, a, out)
        mask = b.alloc_bits(12)
        b.hmr(s, mask)
        b.swap(s, out)
        b.free(out)
        c['masks'][k] = mask
    b.cx(Z0, a[0])
    clip_close(b, m, B)
    b.cswap(F, g[-1], g[-2])                       # unpivot
    b.cx(Z0, g[-2][0])
    c['need_p_undo'] = True                        # after the line repays the deferred masks (erase_normpair hook)


def _find_zero_flag():
    """The zero flag [a == 0] of the inverse_apply that runs the current action (kit: local F; a14 copy: Z)."""
    f = sys._getframe(2)
    while f is not None:
        if f.f_code.co_name.startswith('inverse_apply'):
            for name in ('F', 'Z'):
                v = f.f_locals.get(name)
                if v is not None and not isinstance(v, (list, tuple, dict)):
                    return v
        f = f.f_back
    raise AssertionError('no inverse_apply frame above the action')


def lower_inverse_ai(b, m, ai):
    """Reverse, plain undo path: at the line's shifted_conv undo (inside the action, g pivoted, ai given).
    If the level was handed a known-zero flag Z0 (= [a == 0], built by receive_reverse_gap) and switched to the
    plain path (ZN_PDL2_PLAIN), Z0 is erased here by the plain inverse's own zero flag (one CNOT, as the line's
    pre-a12 handed path did)."""
    LOG['inv_ai'] += 1
    c = CTX[m]
    ev('inv_ai', m)
    Z0 = c.pop('plain_Z0', None)
    if Z0 is not None:
        b.cx(_find_zero_flag(), Z0)
        b.free(Z0)
        LOG['plain_z0'] = LOG.get('plain_z0', 0) + 1
    a = c['g'][-1]
    B = clip_open(b, m)
    y_parts(b, m, B)
    convs(b, m, True)
    for v in c['payloads']:
        c['scale'](b, v[-1], a, ai, False)
    p_parts(b, m)
    clip_close(b, m, B)
    c['skip_scales'] = len(c['payloads'])


def _find_ai():
    f = sys._getframe(2)
    while f is not None:
        if 'ai' in f.f_locals:
            return f.f_locals['ai']
        f = f.f_back
    raise AssertionError('no ai in the calling frames')


class DKProxy:
    def __init__(self, dk):
        self._dk = dk

    def __getattr__(self, name):
        return getattr(self._dk, name)

    def shifted_product(self, b, g, x, out, F):
        m = len(g) + 1
        if (m in UPPER or m in LOWER) and _active(m):
            ev('noprod', m, CTX[m]['undo'])
            return None
        return self._dk.shifted_product(b, g, x, out, F)

    def shifted_conv(self, b, g, q, out, F, reverse=False):
        m = len(g) + 1
        if m in LOWER and _active(m):
            c = CTX[m]
            if not reverse and not c['undo']:
                r = self._dk.shifted_conv(b, g, q, out, F, reverse)
                lower_forward(b, m)
                return r
            if reverse and c['undo'] and not c.get('inv_done'):
                lower_inverse_ai(b, m, _find_ai())
                c['inv_done'] = True
        return self._dk.shifted_conv(b, g, q, out, F, reverse)


def make_scale(orig):
    def scale(b, x, a, ai, divide=True):
        for m in LOWER:
            if _active(m):
                k = _payload_index(m, x)
                if k is not None:
                    c = CTX[m]
                    c['scale'] = orig
                    if not c['undo']:
                        c['a'], c['ai'] = a, ai
                        ev('noscale', m, k)
                        return None
                    if c.get('skip_scales', 0) > 0:
                        c['skip_scales'] -= 1
                        return None
        return orig(b, x, a, ai, divide)
    scale._zn_pdl2 = True
    return scale


def make_restore(orig):
    def high_full_restore(b, f, g, F, *a, **kw):
        m = len(f) - 1
        if m in LOWER and _active(m) and CTX[m]['undo'] and not CTX[m].get('inv_done') \
                and CTX[m].get('handoff_zero') is not None:
            lower_inverse_sbi(b, m)
            CTX[m]['inv_done'] = True
        return orig(b, f, g, F, *a, **kw)
    high_full_restore._zn_pdl2 = True
    return high_full_restore


def make_sbi(orig):
    def scale_before_inverse(b, x, a):
        for m in LOWER:
            if _active(m) and CTX[m]['undo'] and CTX[m].get('masks'):
                k = _payload_index(m, x)
                if k is not None and k in CTX[m]['masks']:
                    ev('sbi_mask', m, k)
                    return CTX[m]['masks'].pop(k)
        return orig(b, x, a)
    scale_before_inverse._zn_pdl2 = True
    return scale_before_inverse


def p_undo_late(b, m, a, ai, Z):
    """Reverse, handed path: the P parts' inverse, after the line has repaid the deferred unscale debts and before
    it erases its inverse (g pivoted, f restored: f[-1] = the dividend's top; q1' re-formed in place by one
    scaling and restored after)."""
    c = CTX[m]
    f = c['f']
    ev('p_undo', m)
    b.cx(Z, a[0])
    b.cx(Z, ai[0])
    c['scale'](b, f[-1], a, ai, True)
    B = clip_open(b, m)
    p_parts(b, m)
    clip_close(b, m, B)
    c['scale'](b, f[-1], a, ai, False)
    b.cx(Z, ai[0])
    b.cx(Z, a[0])
    c['need_p_undo'] = False


def make_erase_normpair(orig):
    def erase_normpair(b, a, ai, Z, saved, *rest, **kw):
        for m in LOWER:
            if _active(m) and CTX[m]['undo'] and CTX[m].get('need_p_undo'):
                if os.environ.get('ZN_PDL2_NPARK') == '1' and isinstance(saved, tuple) and len(saved) >= 2:
                    # fourth holder: the line's norm N = norm(a) (6 qubits, read only by its own erase) is parked
                    # across the late P undo (its 879 moment) and recomputed after it from the unclipped pivot
                    SC = sys.modules['scalar_compact_v3']
                    N = saved[0]
                    mask = b.alloc_bits(len(N))
                    b.hmr(N, mask)
                    b.free(N)
                    p_undo_late(b, m, a, ai, Z)
                    to, frm, K = SC.ar.tower(12)
                    SC.ar.linear_inplace(b, a, to)
                    N2 = b.alloc_qubits(len(mask))
                    SC.norm_xor_noscratch(b, a, N2, K)
                    SC.ar.linear_inplace(b, a, frm)
                    b.cz(N2, mask)
                    b.free(mask)
                    saved = (N2,) + tuple(saved[1:])
                    LOG['npark'] = LOG.get('npark', 0) + 1
                else:
                    p_undo_late(b, m, a, ai, Z)
        return orig(b, a, ai, Z, saved, *rest, **kw)
    erase_normpair._zn_pdl2 = True
    return erase_normpair


def make_recorder(orig, name):
    def rec(b, f, g, payloads, *a, **kw):
        REC[len(f) - 1] = dict(f=f, payloads=payloads)
        return orig(b, f, g, payloads, *a, **kw)
    rec._zn_pdl2 = True
    rec.__name__ = name
    return rec


def make_direct_pre(orig, scope):
    def direct_pre(b, f, g, payloads, F, undo=False, incoming=None, defer=False, **kw):
        m = len(f) - 1
        mine = (m in UPPER or m in LOWER) and not (FWDONLY and undo)
        if mine:
            c = CTX.setdefault(m, {})
            if not c.get('active'):
                hz = kw.get('handoff_zero')
                if hz is not None and undo and m in LOWER and os.environ.get('ZN_PDL2_PLAIN', '0') == '1':
                    kw.pop('handoff_zero')                 # plain undo path at the paired lower level
                    c['plain_Z0'] = hz
                    hz = None
                c.update(f=f, g=g, payloads=payloads, F=F, undo=undo, handoff_zero=hz,
                         inv_done=False,
                         dk=scope['dk']._dk if isinstance(scope.get('dk'), DKProxy) else scope.get('dk'))
                c.pop('masks', None)
                if m in UPPER:
                    CTX.setdefault(m - 1, {})['tops'] = [v[-1] for v in payloads]
            c['active'] = c.get('active', 0) + 1
            ev('pre', m, undo, kw.get('handoff_zero') is not None)
        try:
            return orig(b, f, g, payloads, F, undo, incoming, defer, **kw)
        finally:
            if mine:
                CTX[m]['active'] -= 1
                if CTX[m]['active'] == 0 and m in LOWER:
                    assert not undo or CTX[m].get('inv_done'), 'paired inverse never ran at level %d' % m
                    assert not CTX[m].get('masks'), 'unused deferred masks at level %d' % m
                    assert not CTX[m].get('need_p_undo'), 'P undo never ran at level %d' % m
                    assert CTX[m].get('plain_Z0') is None, 'handed zero never erased at level %d' % m
                    CTX[m].pop('skip_scales', None)
    direct_pre._zn_pdl2 = True
    return direct_pre


def tops(m):
    return [v[-1] for v in REC[m]['payloads']]


def route_pass(b, ups, reverse):
    """Every block's carry in one bottom pass: the algebra of dqi_zn_pdl.route_block (pd3 routes, zero-region
    terms), cheaper: the zero flags of q1 (levels U-1 .. U-4) are held in a sliding window shared by adjacent
    blocks, and the selections are merged:
        t = not C and not F;  h1 = not B and C;  h2 = not B and t;  ZB = h1 q0_(U-2) + h2
        selL = not B and not C2 and not t   (= route_block's sA or sB: disjoint, so one CCX row instead of two)
        selU2 = h1 and C2 and not F
    Bottom-up forward, top-down reverse (each block's carry is an involution for fixed flags)."""
    import dqi_zn_pdc as PC
    SS = _ss()
    ar = _ar()
    fl = {}

    def q1(m):
        return REC[m]['f'][-1]
    def drop(m):                                      # fourth holder: ZN_PDL2_HAND levels go to the moved q_posts
        if not reverse and m in HANDSET:
            HAND[m] = (fl.pop(m), q1(m))
            LOG['handed'] = LOG.get('handed', 0) + 1
        else:
            SS.clear_flag(b, q1(m), fl.pop(m))
    for U in (list(reversed(ups)) if reverse else list(ups)):
        need = {U - 1, U - 2, U - 3, U - 4}
        for m in sorted(set(fl) - need, reverse=True):
            drop(m)
        for m in sorted(need - set(fl)):
            fl[m] = SS.fresh_zero_flag(b, q1(m))
        Bq, Cq, Fq, C2q = fl[U - 1], fl[U - 2], fl[U - 3], fl[U - 4]
        q0u = REC[U - 2]['f'][-2]
        t, lt = PC.flag_and(b, [(Cq, True), (Fq, True)])
        h1, l1 = PC.flag_and(b, [(Bq, True), (Cq, False)])
        ZB = PC.masked_word(b, h1, q0u)
        h2, l2 = PC.flag_and(b, [(Bq, True), (t, False)])
        b.cx(h2, ZB[0])
        sL, lL = PC.flag_and(b, [(Bq, True), (C2q, True), (t, True)])
        sU2, lU2 = PC.flag_and(b, [(h1, False), (C2q, False), (Fq, True)])
        X, tu, tl, t2 = tops(U), tops(U - 2), tops(U - 3), tops(U - 4)
        for k in range(len(X)):
            ar.mul_xor0(b, ZB, X[k], tu[k])
            for j in range(12):
                b.ccx(sL, X[k][j], tl[k][j])
            for j in range(12):
                b.ccx(sU2, X[k][j], t2[k][j])
        PC.erase_flag_and(b, sU2, lU2)
        PC.erase_flag_and(b, sL, lL)
        b.cx(h2, ZB[0])
        PC.erase_flag_and(b, h2, l2)
        PC.erase_masked_word(b, h1, q0u, ZB)
        PC.erase_flag_and(b, h1, l1)
        PC.erase_flag_and(b, t, lt)
    for m in sorted(fl, reverse=True):
        drop(m)


HANDSET = {int(x) for x in os.environ.get('ZN_PDL2_HAND', '').split(',') if x}
HAND = {}


def install_hand(decode):
    """The kit q_post that Z2's bottom pass calls (scalar_compact_v3's `original` = decode_backhalf_q7's q_post) takes
    its gap flag F = [q == 0] from fresh_zero_flag(b, q): when route_pass holds that very predicate for the same word
    (HAND), it is handed over instead of a fresh zero test (route_pass's clear and the q_post's test both saved)."""
    walk = sys.modules['walk_bound']
    n = 0
    for fn in walk.functions(decode):
        if fn.__name__ == 'q_post' and fn.__code__.co_filename.endswith('scalar_compact_v3.py') and fn.__closure__:
            cells = dict(zip(fn.__code__.co_freevars, fn.__closure__))
            if 'original' not in cells:
                continue
            gl = cells['original'].cell_contents.__globals__
            orig = gl.get('fresh_zero_flag')
            if orig is None or getattr(orig, '_zn_hand', False):
                continue

            def fresh_zero_flag(b, w, _o=orig):
                for m, (F, word) in list(HAND.items()):
                    if word is w:
                        del HAND[m]
                        LOG['hand_used'] = LOG.get('hand_used', 0) + 1
                        return F
                return _o(b, w)
            fresh_zero_flag._zn_hand = True
            gl['fresh_zero_flag'] = fresh_zero_flag
            n += 1
    LOG['hand_scopes'] = n
    return n


def make_root_action(orig):
    """The un-fused root (level 17) runs root_action, not direct_pre: activate it as an upper level for its PRE
    action (post=False: the forward descent and its undo at the reverse return), so that DKProxy skips its
    shifted_products ('u' ops); record its payloads; hand its tops to level 16."""
    def root_action(b, f, g, payloads, ai, Z, F, post, undo, *rest, **kw):
        m = len(f) - 1
        mine = (m in UPPER) and not post and not (FWDONLY and undo)
        if mine:
            REC[m] = dict(f=f, payloads=payloads)
            c = CTX.setdefault(m, {})
            if not c.get('active'):
                c.update(f=f, g=g, payloads=payloads, F=F, undo=undo, handoff_zero=None, inv_done=False)
                CTX.setdefault(m - 1, {})['tops'] = [v[-1] for v in payloads]
            c['active'] = c.get('active', 0) + 1
            ev('root', m, undo)
        try:
            return orig(b, f, g, payloads, ai, Z, F, post, undo, *rest, **kw)
        finally:
            if mine:
                CTX[m]['active'] -= 1
    root_action._zn_pdl2 = True
    return root_action


def bottom(b, reverse):
    import dqi_zn_pdl as V1
    V1.REC = REC                                   # route_block reads V1.REC through its module globals
    LOG['bottom'] += 1
    if os.environ.get('ZN_PDL2_NOCARRY') == '1':      # diagnostics: forward words then differ, round trip must not
        return
    if FWDONLY and reverse:
        return
    ups = sorted(UPPER)
    if os.environ.get('ZN_PDL2_ROUTE', 'pass') == 'v1':
        for U in (reversed(ups) if reverse else ups):
            V1.route_block(b, U)
        return
    route_pass(b, ups, reverse)


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
    base_case._zn_pdl2 = True
    return base_case


def install(km, ns):
    walk = sys.modules['walk_bound']
    fns = list(walk.functions(ns.decode))
    dicts = []
    for fn in fns:
        dicts.append(fn.__globals__)
        for cell in (fn.__closure__ or ()):
            try:
                v = cell.cell_contents
            except ValueError:
                continue
            if isinstance(v, dict) and any(k in v for k in ('dk', 'scale', 'direct_pre', 'base_case',
                                                            'high_full_restore')):
                dicts.append(v)
    done = set()
    for gl in dicts:
        if id(gl) in done:
            continue
        done.add(id(gl))
        hit = False
        dkm = gl.get('dk')
        if type(dkm).__name__ == 'module' and dkm.__name__.endswith('direct_kernels'):
            gl['dk'] = DKProxy(dkm)
            hit = True
        sc = gl.get('scale')
        co = getattr(sc, '__code__', None)
        if callable(sc) and co is not None and not getattr(sc, '_zn_pdl2', False) \
                and co.co_varnames[:5] == ('b', 'x', 'a', 'ai', 'divide'):
            gl['scale'] = make_scale(sc)
            hit = True
        hr = gl.get('high_full_restore')
        if callable(hr) and not getattr(hr, '_zn_pdl2', False):
            gl['high_full_restore'] = make_restore(hr)
            LOG['restore_wraps'] = LOG.get('restore_wraps', 0) + 1
        for name in ('direct_pre', 'pre'):
            f0 = gl.get(name)
            if callable(f0) and not getattr(f0, '_zn_pdl2', False) and hasattr(f0, '__code__'):
                w = make_recorder(f0, name)
                if name == 'direct_pre':
                    w = make_direct_pre(w, gl)
                gl[name] = w
                LOG['pre_wraps'] += 1
        ia = gl.get('inverse_apply')
        if LATECONV and callable(ia) and not getattr(ia, '_zn_pdl2', False) and 'direct_pre' in gl:
            gl['inverse_apply'] = make_inverse_apply(ia)
            LOG['ia_wraps'] = LOG.get('ia_wraps', 0) + 1
        bc = gl.get('base_case')
        if callable(bc) and not getattr(bc, '_zn_pdl2', False) and hasattr(bc, '__code__'):
            gl['base_case'] = make_base(bc)
            LOG['base_wraps'] += 1
        ra = gl.get('root_action')
        if 17 in UPPER and callable(ra) and not getattr(ra, '_zn_pdl2', False) and hasattr(ra, '__code__') \
                and ra.__code__.co_varnames[:9] == ('b', 'f', 'g', 'payloads', 'ai', 'Z', 'F', 'post', 'undo'):
            gl['root_action'] = make_root_action(ra)
            LOG['root_wraps'] = LOG.get('root_wraps', 0) + 1
        LOG['scopes'] += int(hit)
    nsbi = 0
    for mod in list(sys.modules.values()):           # every module object that carries scale_before_inverse
        fn = getattr(mod, 'scale_before_inverse', None)
        if callable(fn) and not getattr(fn, '_zn_pdl2', False):
            try:
                setattr(mod, 'scale_before_inverse', make_sbi(fn))
                nsbi += 1
            except (AttributeError, TypeError):
                pass
    for gl in dicts:                                  # and namespaces in scopes that carry it
        for key, val in list(gl.items()):
            fn = getattr(val, 'scale_before_inverse', None) if not isinstance(val, (int, str, float)) else None
            if callable(fn) and not getattr(fn, '_zn_pdl2', False):
                try:
                    setattr(val, 'scale_before_inverse', make_sbi(fn))
                    nsbi += 1
                except (AttributeError, TypeError):
                    pass
    LOG['sbi_wraps'] = nsbi
    nen = 0
    for holder in list(sys.modules.values()) + [v for gl in dicts for v in list(gl.values())]:
        if isinstance(holder, (int, str, float, dict, list, tuple)) or holder is None:
            continue
        fn = getattr(holder, 'erase_normpair', None)
        if callable(fn) and not getattr(fn, '_zn_pdl2', False):
            try:
                setattr(holder, 'erase_normpair', make_erase_normpair(fn))
                nen += 1
            except (AttributeError, TypeError):
                pass
    for gl in dicts:                                  # cloned scopes that hold the name directly
        fn = gl.get('erase_normpair')
        if callable(fn) and not getattr(fn, '_zn_pdl2', False):
            gl['erase_normpair'] = make_erase_normpair(fn)
            nen += 1
    LOG['erase_normpair_wraps'] = nen
    if HANDSET:
        install_hand(ns.decode)
    return dict(LOG)
