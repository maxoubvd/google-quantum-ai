#!/usr/bin/env python3
"""r48_dqi_conv2_rp_root4d.py -- R48-DQI-CONV2 re-pricer (run 5): ledger row 4 (the root's return without its
Gram pair) in a cheaper form than B5's dqi_int_root4.py (ec1dea59, r4a: -60 per direction [M]), plus a width probe.

COPY of B5's tools/dqi_int_root4.py (install logic, scale skipping and root_action wrapper verbatim); only the dot
that replaces the root's (scale ; one-word dot ; scale) changes.  The line, per payload top X and direction:
    scale(X, a) [42] ; a16's one-word dot X ^= D0 + F (D0 + D1) [16 points + 12] ; scale(X, ai) [42]
    = X ^= ai sel_F(D0, D1)   (identity C; a ai = 1 is what the line's scale already needs).
B5's form (R48RP_R4FORM=b5, verbatim): the kit's two outputs d0 = D0, d1 = D1, cswap(F) [12], X ^= ai d0 [42],
cswap back [12], free erasures: 16 points + 66 per payload and direction (-30).
THIS form (R48RP_R4FORM=d): e ^= D0, d ^= D0 + D1 (a16's coefficient pair [1, 1 ^ t] into two fresh words),
e ^= F d [12 CCX] -> e = sel_F(D0, D1); X ^= ai e [42]; e X-measured (mask re): its phase is
(-1)^(re . D0) (-1)^(F re . d), and with d still live the second factor is repaid by ccz(F, d_j, re_j) (a CZ under a
classical control: 0 Toffolis, the same gate phase_product uses); then d X-measured (mask rd); the bilinear rest
(-1)^(re . D0 + rd . (D0 + D1)) is repaid point by point by the kit's own erasure with mask_t = re + (1 + t)^T rd.
16 points + 54 per payload and direction: -42, i.e. -84 per direction, -168 per decode [P].  Width: e and d are
24 wires, the same as B5's d0, d1 (+12 over a16's one-word dot at that moment).
R48RP_R4FORM=off installs no gate change (probe only).  R48RP_DIRS=fwd changes only the forward root post (undo
False; the reverse keeps the line's scale ; dot ; scale, the same map X ^= ai sel, so the round trip stays exact):
the decode timeline wd1 puts unconvert's root post-undo at 879 + 12 = 891 with a16's one-word dot [M rows 237-239],
i.e. 879 + 24 = 903 with two words, while convert's root post sits at 878 + 24 = 902 [M rows 712-714].
Probe (any form): R48RP_CAP=<n> sets the builder's own cap (CircuitBuilder.set_max_qubits) at its first
allocation, so ANY allocation, also inside native kernels, beyond n raises at the exact call site; without the cap
the module records the builder's own live count (num_allocated_qubits) high-water per harness mode at Python-level
allocations and prints it at exit with its call site.  R48RP_PROBE=0 installs no probe (no alloc wrapper: decode
builds).
Install as a DQI_INT_EXT module where B5 installs dqi_int_root4 (after r46_dqi_ext_fgret, before r46_dqi_conv_zpark).
"""
import os
import sys
import traceback

FORM = os.environ.get('R48RP_R4FORM', 'd')
CAP = int(os.environ.get('R48RP_CAP', '0') or 0)
DIRS = os.environ.get('R48RP_DIRS', 'both')                  # both | fwd (convert's root post only) | rev
LOG = dict(form=FORM, cap=CAP, dirs=DIRS, scopes=0, old_wraps=0, posts=0, scales_skipped=0, dots=0, fallbacks=0)
ST = dict(tops=None, ai={}, skipped=set())
HW = {}


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


def _erase_points(b, g, v, k, ar, N, V, W, m0, m1, c1):
    """Repay (-1)^(m0 . D0 + m1 . sum_t c1(t) P_t) point by point (the kit's erasure; P_t = <g, V_t> <v, W_t>)."""
    f = ar.field(12)
    for t in range(N):
        mask = b.alloc_bits(12)
        b.bit_store(mask, False)
        if m0 is not None:
            b.cx(m0, mask)
        for j in range(12):
            value = f.mul(c1(t), 1 << j)
            for i in ar.bits(value):
                b.cx(m1[i], mask[j])
        x, p = k.expose(b, g, V[t])
        y, r = k.expose(b, v, [W[i][t] for i in range(N)])
        ar.phase_product(b, x, y, mask)
        k.unexpose(b, v, [W[i][t] for i in range(N)], r)
        k.unexpose(b, g, V[t], p)
        b.free(mask)


def scaled_dot_b5(b, g, v, target, F, ai):
    """B5's form, verbatim: target ^= ai * (F ? D1 : D0), the kit's two transposed outputs."""
    dkm = sys.modules['direct_kernels']
    k, ar, rank_update = dkm.k, dkm.ar, dkm.rank_update
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


def scaled_dot_d(b, g, v, target, F, ai):
    """target ^= ai * (D0 + F (D0 + D1)): one 12-CCX select into a fresh word, one product, free erasures."""
    dkm = sys.modules['direct_kernels']
    k, ar, rank_update = dkm.k, dkm.ar, dkm.rank_update
    e = b.alloc_qubits(12)
    d = b.alloc_qubits(12)
    N = len(g)
    V, W = k.vandermonde(N)
    for t in range(N):                                   # e ^= P_t, d ^= (1 + t) P_t  (a16's coefficient pair)
        rank_update(b, g, V[t], v, [W[i][t] for i in range(N)], [e, d], [1, 1 ^ t])
    for j in range(12):                                  # e = D0 + F (D0 + D1) = sel_F(D0, D1)
        b.ccx(F, d[j], e[j])
    ar.mul_xor0(b, ai, e, target)                        # target ^= ai * sel
    re = b.alloc_bits(12)
    b.hmr(e, re)
    b.free(e)
    for j in range(12):                                  # the F part of (-1)^(re . e), d live: classical-control CZ
        b.ccz(F, d[j], re[j])
    rd = b.alloc_bits(12)
    b.hmr(d, rd)
    b.free(d)
    _erase_points(b, g, v, k, ar, N, V, W, re, rd, lambda t: 1 ^ t)
    b.free(re)
    b.free(rd)


SCALED = dict(b5=scaled_dot_b5, d=scaled_dot_d)


class Root4DK:
    def __init__(self, dk):
        self._dk = dk

    def __getattr__(self, name):
        return getattr(self._dk, name)

    def selected_dot(self, b, g, v, target, F):
        kk = _top_index(target)
        if kk is not None and kk in ST['skipped']:
            SCALED[FORM](b, g, v, target, F, ST['ai'][kk])
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
        if not post or (undo and DIRS == 'fwd') or (not undo and DIRS == 'rev'):
            return orig(b, f, g, payloads, ai, Z, F, post, undo, *rest, **kw)
        prev = (ST['tops'], dict(ST['ai']), set(ST['skipped']))
        ST.update(tops=[v[-1] for v in payloads], ai={}, skipped=set())
        try:
            out = orig(b, f, g, payloads, ai, Z, F, post, undo, *rest, **kw)
            assert not ST['skipped'], 'root4d: a first scaling without its last'
            LOG['posts'] += 1
            return out
        finally:
            ST['tops'], ST['ai'], ST['skipped'] = prev[0], prev[1], prev[2]
    root_action._b5_root4 = True
    return root_action


def _site():
    st = traceback.extract_stack()[:-2]
    keep = [f for f in st if 'kickmix' not in f.filename and 'root4d' not in f.filename and 'convtest' not in f.filename]
    return ' < '.join('%s:%d:%s' % (os.path.basename(f.filename), f.lineno, f.name) for f in reversed(keep[-6:]))


def _install_probe(km):
    if km is None:                                           # dqi_int_ext passes km = None: take kickmix itself
        km = sys.modules.get('kickmix') or sys.modules.get('circ')
    if km is None:
        LOG['probe'] = 'no kickmix module'
        return
    CB = km.CircuitBuilder
    inner = CB.alloc_qubits

    def alloc_probe(self, n, *a, **kw):
        if CAP and self.max_qubits > CAP:
            self.set_max_qubits(CAP)
        r = inner(self, n, *a, **kw)
        ct = sys.modules.get('dqi_int_convtest')
        mode = ct.TR.get('mode', '?') if ct is not None else '?'
        tr_live = len(ct.TR['live']) if ct is not None else -1
        own = self.num_allocated_qubits
        h = HW.get(mode)
        if h is None or own > h[0]:
            HW[mode] = (own, tr_live, _site())
        return r
    CB.alloc_qubits = alloc_probe


def install(km, ns):
    if os.environ.get('R48RP_PROBE', '1') == '1':          # R48RP_PROBE=0: no alloc wrapper (decode builds)
        _install_probe(km)
    import atexit
    atexit.register(lambda: print('r48_dqi_conv2_rp_root4d', LOG, 'own-live high-water by mode:',
                                  {m: (h[0], h[1]) for m, h in HW.items()}, '\n  sites:',
                                  {m: h[2] for m, h in HW.items()}, flush=True))
    if FORM == 'off':
        return dict(LOG)
    assert FORM in SCALED, FORM
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
    assert LOG['scopes'] >= 1, 'root4d: root_fusion scope not found'
    return dict(LOG)
