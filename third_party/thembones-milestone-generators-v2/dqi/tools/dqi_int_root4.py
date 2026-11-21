#!/usr/bin/env python3
"""dqi_int_root4.py -- B5 (14th holder, Nov 16): R48-DQI-CONV2 ledger row 4, the root's return without its Gram pair.

The root's post (forward: root_fusion.root_action, post=True) and its undo (reverse: direct.root_action, post=True,
undo=True) both run, per payload top X:  scale(X, a, ai, False) ; X ^= sel_F(D0, D1) ; scale(X, a, ai, True),
i.e. X <- (X a + D) ai.  With a ai = 1 (the line's scale already needs it) this is X ^= ai D (identity C, the
designer's model 600 / 600).  Here: the two scalings are skipped and the dot becomes
    d0, d1 = the kit's two transposed outputs (direct_kernels.selected_dot: 24 fresh wires, one product per point)
    cswap(F, d0, d1) ; X ^= ai * d0 ; cswap(F, d0, d1) ; d0, d1 erased exactly as the kit erases them.
Per payload and direction: 16 points + 12 + 42 + 12 against the line's 16 points + 12 (a16's one-word select) + 84:
-30 [P]; 2 payloads x 2 directions = -120 per decode.  Width: +12 against a16's one-word dot at that moment
(root dot sites 866 fwd / 852 rev in the converter harness [M h17 s16ref]: 878 / 864).
Install as a DQI_INT_EXT module after the modules that wrap the root (pdl2, kpark, alias, pdr) and before
r46_dqi_conv_zpark.  Switch: B5_ROOT4=0 installs nothing.
"""
import os
import sys

ON = os.environ.get('B5_ROOT4', '1') == '1'
# 15th holder: B5_ROOT4_DIR = both (default) | fwd (convert's root post only) | rev (unconvert's root post-undo only)
DIR = os.environ.get('B5_ROOT4_DIR', 'both')
SINGLE = os.environ.get('B5_ROOT4_SINGLE', '0') == '1'   # allocate the two dot words one qubit at a time
# 15th holder: B5_ROOT4_EVICT=1 frees one wire for the REVERSE root post-undo (there the line runs at its peak - 12 and
# root4 needs + 24: one over [M r4w]). The root's clip flag Z = [a == 0] is not read by the two dots: before the first
# dot a is unclipped and Z measured out (bit zb); Z's wire is lent to the dot word d0; after the second dot Z is
# recomputed into its own wire by the kit's zero_toggle (11 ANDs), the measurement's phase repaired on it (CZ by zb,
# free) and a clipped again, so root_action's own unclip and erase_product_zero see the line's state. Cost +11.
EVICT = os.environ.get('B5_ROOT4_EVICT', '0') == '1'
LOG = dict(scopes=0, old_wraps=0, posts=0, scales_skipped=0, dots=0, fallbacks=0)
ST = dict(tops=None, ai={}, skipped=set())


def _same(x, w):
    return x is w or (len(x) == len(w) and len(x) > 0 and x[0] is w[0])


def _top_index(x):
    tops = ST['tops']
    if tops is None:
        return None
    for k, t in enumerate(tops):
        if _same(x, t):
            return k
    return None


def scaled_dot(b, g, v, target, F, ai):
    """target ^= ai * (F ? D1 : D0), D0 / D1 = the kit's two transposed outputs of <g, v> (direct_kernels)."""
    dkm = sys.modules['direct_kernels']
    k, ar, rank_update = dkm.k, dkm.ar, dkm.rank_update
    ev = ST.get('ev')
    lent = None
    if ev is not None:                                   # reverse root post-undo with B5_ROOT4_EVICT=1
        lent = _evict_start(b, ev)
        d = [[lent] + list(b.alloc_qubits(11)), b.alloc_qubits(12)]
    elif SINGLE:                                         # 15th holder: one qubit at a time (fills index holes)
        d = [[b.alloc_qubits(1)[0] for _ in range(12)] for _ in range(2)]
    else:
        d = [b.alloc_qubits(12), b.alloc_qubits(12)]
    N = len(g)
    V, W = k.vandermonde(N)
    for t in range(N):
        rank_update(b, g, V[t], v, [W[i][t] for i in range(N)], d, [1, t])
    b.cswap(F, d[0], d[1])
    ar.mul_xor0(b, ai, d[0], target)
    b.cswap(F, d[0], d[1])
    masks = [b.alloc_bits(12), b.alloc_bits(12)]
    for x, m in zip(d, masks):
        b.hmr(x, m)
        if lent is not None and x[0] is lent:
            b.free(x[1:])                                # Z's wire stays allocated (|0> after the hmr)
        else:
            b.free(x)
    f = ar.field(12)
    for t in range(N):                                   # the kit's erasure, unchanged
        mask = b.alloc_bits(12)
        b.bit_store(mask, False)
        b.cx(masks[0], mask)
        for j in range(12):
            value = f.mul(t, 1 << j)
            for i in ar.bits(value):
                b.cx(masks[1][i], mask[j])
        x, p = k.expose(b, g, V[t])
        y, r = k.expose(b, v, [W[i][t] for i in range(N)])
        ar.phase_product(b, x, y, mask)
        k.unexpose(b, v, [W[i][t] for i in range(N)], r)
        k.unexpose(b, g, V[t], p)
        b.free(mask)
    for m in masks:
        b.free(m)
    if ev is not None:
        _evict_end(b, ev)


def _one(z):
    """The root's Z arrives as a one-qubit span: the single-bit hmr / cz calls need the qubit itself."""
    try:
        return z[0] if len(z) == 1 else z
    except TypeError:
        return z


def _evict_start(b, ev):
    """Reverse root post-undo, first dot: the wire to lend. If dqi_int_flagpark's part h holds its clean placeholder
    (h_hold: H's wire kept |0> from the reverse root's start to h_rebuild, after this post-undo), that wire is lent
    as it is (0 Toffolis: the dot's own hmr leaves it |0>, still allocated). Otherwise a unclipped, Z measured out
    (bit kept), and Z's wire lent (rebuilt after the last dot, +11)."""
    if ev.get('lend_only'):
        return ev['Z']
    if ev['n'] == 0:
        b.cx(ev['Z'], ev['a'][0])
        zb = b.alloc_bits(1)
        b.hmr(ev['Z'], zb[0])
        ev['zb'] = zb
        LOG['evicted'] = LOG.get('evicted', 0) + 1
    return ev['Z']


def _evict_end(b, ev):
    """After each dot; after the last one Z is rebuilt into its own wire (zero_toggle, 11 ANDs), its measurement
    phase repaired on the rebuilt copy (CZ by the kept bit) and a clipped again."""
    ev['n'] += 1
    if ev.get('lend_only'):
        if ev['n'] == ev['total']:
            ev['done'] = True
            LOG['lent_placeholder'] = LOG.get('lent_placeholder', 0) + 1
        return
    if ev['n'] == ev['total']:
        from exp.safe_scalar import converter as SS
        SS.zero_toggle(b, ev['a'], ev['Z'])
        b.cz(ev['Z'], ev['zb'][0])
        b.free(ev['zb'])
        b.cx(ev['Z'], ev['a'][0])
        ev['done'] = True
        LOG['restored'] = LOG.get('restored', 0) + 1


def scaled_dot_d(b, g, v, target, F, ai):
    """R48-DQI-CONV2 re-pricer's form d (copied from tools/r48_dqi_conv2_rp_root4d.py, run 5, harness r6d exact):
    target ^= ai * (D0 + F (D0 + D1)): e ^= D0, d ^= D0 + D1 (a16's pair [1, 1 ^ t]), e ^= F d (12 CCX), X ^= ai e (42),
    e X-measured, its F part repaid by ccz(F, d_j, re_j) under the classical bit (0 Toffolis) while d is live, d
    X-measured, the bilinear rest repaid point by point by the kit's erasure.  16 points + 54 per payload and
    direction.  With B5_ROOT4_EVICT (reverse): e borrows the root clip flag's wire (see _evict_*)."""
    dkm = sys.modules['direct_kernels']
    k, ar, rank_update = dkm.k, dkm.ar, dkm.rank_update
    ev = ST.get('ev')
    lent = _evict_start(b, ev) if ev is not None else None
    e = ([lent] + list(b.alloc_qubits(11))) if lent is not None else b.alloc_qubits(12)
    d = b.alloc_qubits(12)
    N = len(g)
    V, W = k.vandermonde(N)
    for t in range(N):                                   # e ^= P_t, d ^= (1 + t) P_t
        rank_update(b, g, V[t], v, [W[i][t] for i in range(N)], [e, d], [1, 1 ^ t])
    for j in range(12):                                  # e = sel_F(D0, D1)
        b.ccx(F, d[j], e[j])
    ar.mul_xor0(b, ai, e, target)                        # target ^= ai * sel
    re = b.alloc_bits(12)
    b.hmr(e, re)
    b.free(e[1:] if lent is not None else e)             # a lent wire stays allocated (|0> after the hmr)
    for j in range(12):
        b.ccz(F, d[j], re[j])
    rd = b.alloc_bits(12)
    b.hmr(d, rd)
    b.free(d)
    f = ar.field(12)
    for t in range(N):                                   # mask_t = re + (1 + t)^T rd (the kit's erasure)
        mask = b.alloc_bits(12)
        b.bit_store(mask, False)
        b.cx(re, mask)
        for j in range(12):
            value = f.mul(1 ^ t, 1 << j)
            for i in ar.bits(value):
                b.cx(rd[i], mask[j])
        x, p = k.expose(b, g, V[t])
        y, r = k.expose(b, v, [W[i][t] for i in range(N)])
        ar.phase_product(b, x, y, mask)
        k.unexpose(b, v, [W[i][t] for i in range(N)], r)
        k.unexpose(b, g, V[t], p)
        b.free(mask)
    b.free(re)
    b.free(rd)
    if ev is not None:
        _evict_end(b, ev)


FORM = os.environ.get('B5_ROOT4_FORM', 'b5')             # b5 (scaled_dot, this file's first form) | d


class Root4DK:
    def __init__(self, dk):
        self._dk = dk

    def __getattr__(self, name):
        return getattr(self._dk, name)

    def selected_dot(self, b, g, v, target, F):
        kk = _top_index(target)
        if kk is not None and kk in ST['skipped']:
            (scaled_dot_d if FORM == 'd' else scaled_dot)(b, g, v, target, F, ST['ai'][kk])
            LOG['dots'] += 1
            return None
        if kk is not None:
            LOG['fallbacks'] += 1
        return self._dk.selected_dot(b, g, v, target, F)


def make_old_scale(orig):
    def scale(b, x, a, ai, divide=True):
        kk = _top_index(x)
        if kk is not None:
            if not divide and kk not in ST['skipped']:     # the first scaling (X a): skipped, ai remembered
                ST['skipped'].add(kk)
                ST['ai'][kk] = ai
                LOG['scales_skipped'] += 1
                return None
            if divide and kk in ST['skipped']:             # the last scaling (X ai): skipped
                ST['skipped'].discard(kk)
                LOG['scales_skipped'] += 1
                return None
        return orig(b, x, a, ai, divide)
    scale._b5_root4 = True
    return scale


def make_root_action(orig):
    def root_action(b, f, g, payloads, ai, Z, F, post, undo, *rest, **kw):
        if not post or (undo and DIR == 'fwd') or (not undo and DIR == 'rev'):
            return orig(b, f, g, payloads, ai, Z, F, post, undo, *rest, **kw)
        prev = (ST['tops'], dict(ST['ai']), set(ST['skipped']), ST.get('ev'))
        ST.update(tops=[v[-1] for v in payloads], ai={}, skipped=set(), ev=None)
        if EVICT and undo:                               # the reverse root post-undo: one wire lent to the dots
            fp = sys.modules.get('dqi_int_flagpark')
            hold = getattr(fp, 'H_HOLD', None) if fp is not None else None
            if hold:                                     # flagpark part h's clean placeholder (Z is then False)
                ST['ev'] = dict(Z=hold[-1], lend_only=True, n=0, total=len(payloads), done=False)
            elif Z is not False and Z is not None:       # no part h: the root's clip flag Z itself
                ST['ev'] = dict(Z=_one(Z), a=g[-1], n=0, zb=None, total=len(payloads), done=False)
            else:
                LOG['no_lend'] = LOG.get('no_lend', 0) + 1
        try:
            out = orig(b, f, g, payloads, ai, Z, F, post, undo, *rest, **kw)
            assert not ST['skipped'], 'root4: a first scaling without its last'
            if ST['ev'] is not None:
                assert ST['ev']['done'] and ST['ev']['n'] == ST['ev']['total'], 'root4: eviction not closed'
            LOG['posts'] += 1
            return out
        finally:
            ST['tops'], ST['ai'], ST['skipped'], ST['ev'] = prev[0], prev[1], prev[2], prev[3]
    root_action._b5_root4 = True
    return root_action


def install(km, ns):
    if not ON:
        return dict(LOG, on=False)
    walk = sys.modules['walk_bound']
    seen = set()
    olds = set()
    for fn in walk.functions(ns.decode):
        cands = [fn.__globals__]
        for cell in (fn.__closure__ or ()):
            try:
                v = cell.cell_contents
            except ValueError:
                continue
            if isinstance(v, dict):
                cands.append(v)
        for gl in cands:
            if id(gl) in seen:
                continue
            seen.add(id(gl))
            if not all(key in gl for key in ('ordinary_root_action', 'old', 'dk', 'root_restore', 'root_action')):
                continue
            if gl.get('_b5_root4_done'):
                continue
            gl['_b5_root4_done'] = True
            gl['dk'] = Root4DK(gl['dk'])
            o = gl['old']
            if id(o) not in olds and not getattr(getattr(o, 'scale', None), '_b5_root4', False):
                o.scale = make_old_scale(o.scale)
                olds.add(id(o))
                LOG['old_wraps'] += 1
            if not getattr(gl['root_action'], '_b5_root4', False):
                gl['root_action'] = make_root_action(gl['root_action'])
            LOG['scopes'] += 1
    assert LOG['scopes'] >= 1, 'root4: root_fusion scope not found'
    import atexit
    atexit.register(lambda: print('dqi_int_root4', LOG, flush=True))
    return dict(LOG)
