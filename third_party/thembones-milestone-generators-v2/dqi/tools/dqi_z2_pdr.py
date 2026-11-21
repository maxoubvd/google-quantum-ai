#!/usr/bin/env python3
"""dqi_z2_pdr.py -- builder Z2 (Nov 14): the paired RETURN of the converter, installed into the LINE (r12hz's
converter, narrow240 direct levels), forward direction. Paired blocks (U, L = U - 1), U in Z2_PDR_BLOCKS (default 13).

Identity used (line's direct convention, derived in research/ws_dqi_z2.md): U's return dot over L's final payload is
  <g_U, v_L> = q0_L' s_L,pre + q1_L' <w g_L, v_L> + s_(L-1),pre        (generic; gaps: Z1's routes, kit pd3)
so U's dots are dropped and its top X receives:
  at the bottom (after every q_post of levels 2..max U, moved there in ascending order):
      X ^= [not B] q0_L' s_L,pre  +  ZB s_(U-2) + (sA + sB) s_(U-3) + sU2 s_(U-4)    (Z1's route_block, transposed)
  at L's post (the three-output transposed product replaces L's dot; F = L's gap flag):
      s ^= o0 + F (o0 + o1);   X ^= c (o1 + F (o1 + o2));   then after L's scaling  X ^= CG (s + <r, low>)
  with o_k = sum_t p_t^k g(p_t) y_t (N = len(low) points), c = q1_L' + B, B = [q1_L' == 0], CG = c (F ? g[-3] : g[-2]).
Temporaries at L's point products: three words (d = o0 + o1, e1 = o1, e2 = o1 + o2), all F-independent, erased by
measurement with free phase products.
Usage: Z2_P=pdr python dqi_z2_lct.py run ...        (env Z2_PDR_BLOCKS=13,11,9)
"""
import os
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)

UPPER_F = tuple(int(x) for x in os.environ.get('Z2_PDR_BLOCKS', '13').split(',') if x)
UPPER_R = tuple(int(x) for x in os.environ.get('Z2_PDR_RBLOCKS', os.environ.get('Z2_PDR_BLOCKS', '13')).split(',') if x)
SETS = {}
# Z2_PDR_LOCAL=1 (second holder): only the q_posts of levels min(U) - 4 .. max(U) move (the routes read slots
# U-1 .. U-4 after their Q and before their ops; lower levels keep the line's own order), so the bottom pass runs at
# level min(U) - 4's post and the reverse never needs the kit levels' q_post inverses below level 6.
LOCAL = os.environ.get('Z2_PDR_LOCAL', '0') == '1'
for _d, _u in (('fwd', UPPER_F), ('rev', UPPER_R)):
    _lo = max(2, min(_u) - 4) if (_u and LOCAL) else 2
    SETS[_d] = dict(upper=_u, lower=tuple(u - 1 for u in _u), moved=tuple(range(_lo, max(_u) + 1)) if _u else ())
UPPER, LOWER, MOVED = UPPER_F, SETS['fwd']['lower'], SETS['fwd']['moved']     # forward sets (install log)
ALL_MOVED = tuple(sorted(set(SETS['fwd']['moved']) | set(SETS['rev']['moved'])))


def S(name):
    """The block sets of the current direction."""
    return SETS['rev' if ST['mode'] == 'rev' else 'fwd'][name]
LOG = dict(scopes=0, kit_scopes=0, bottom=0, upper_skips=0, lower_ops=0, cg=0, qpost_skips=0, ctx=0)
ST = dict(mode='-', ctx={}, bottom_done=False, post=None, origq={}, dk=None, scope=None, cur=None, kitq=None,
          pend={}, rpend={})
REV = os.environ.get('Z2_PDR_REV', '0') == '1'
QONLY = os.environ.get('Z2_PDR_MODE', '') == 'qonly'      # debug: only the q_post move
DEBUG = os.environ.get('Z2_PDR_MODE', '') == 'debug'      # debug: route / L terms into kept words, U's dot kept
DBG = {}                                                  # name -> word (never freed), read by the harness


# ------------------------------------------------------------------------------------------- small helpers
def _ar():
    return ST['dk'].ar


def _field():
    return _ar().field(12)


def masked_word(b, F, x, neg=False):
    z = b.alloc_qubits(12)
    if neg:
        b.x(F)
    for j in range(12):
        b.init_and(F, x[j], z[j])
    if neg:
        b.x(F)
    return z


def erase_masked_word(b, F, x, z, neg=False):
    m = b.alloc_bits(12)
    b.hmr(z, m)
    b.free(z)
    if neg:
        b.x(F)
    for j in range(12):
        b.ccz(F, x[j], m[j])
    if neg:
        b.x(F)
    b.free(m)


def flag_and(b, lits):
    out = b.alloc_qubits(1)[0]
    for q, neg in lits:
        if neg:
            b.x(q)
    if len(lits) == 1:
        b.cx(lits[0][0], out)
    else:
        anc, left = [], lits[0][0]
        for q, _ in lits[1:-1]:
            a = b.alloc_qubits(1)[0]
            b.init_and(left, q, a)
            anc.append((left, q, a))
            left = a
        b.init_and(left, lits[-1][0], out)
        for x, y, z in reversed(anc):
            b.del_and(x, y, z)
            b.free(z)
    for q, neg in lits:
        if neg:
            b.x(q)
    return out


def erase_flag_and(b, out, lits):
    for q, neg in lits:
        if neg:
            b.x(q)
    if len(lits) == 1:
        b.cx(lits[0][0], out)
        b.free(out)
    else:
        m = b.alloc_bits(1)
        b.hmr(out, m[0])
        b.free(out)
        if len(lits) == 2:
            b.ccz(lits[0][0], lits[1][0], m[0])
        else:
            anc, left = [], lits[0][0]
            for q, _ in lits[1:-1]:
                a = b.alloc_qubits(1)[0]
                b.init_and(left, q, a)
                anc.append((left, q, a))
                left = a
            b.ccz(left, lits[-1][0], m[0])
            for x, y, z in reversed(anc):
                b.del_and(x, y, z)
                b.free(z)
        b.free(m)
    for q, neg in lits:
        if neg:
            b.x(q)


def zero_flag(b, w):
    return ST['scope']['fresh_zero_flag'](b, w)


LIVEB = os.environ.get('Z2_PDR_LIVEB', '0') == '1'


def get_B(b, L, q1):
    """B = [q1_L' == 0] = the gap flag of level U = L + 1. Z2_PDR_LIVEB=1 (second holder): take the live qubit the
    line already holds for that predicate (forward: returned_gap[U] when U is in KEEP_CHILD, i.e. not parked by the
    bottom pass; reverse: given_gap[U], live through L's whole reverse post), else a fresh zero test.
    Returns (flag, owned): an owned flag is cleared by the caller."""
    ph = ST.get('phase')
    if LIVEB and ph is not None:
        F = ph['returned_gap' if ST['mode'] == 'fwd' else 'given_gap'].get(L + 1)
        if F is not None:
            LOG['liveB'] = LOG.get('liveB', 0) + 1
            return F, False
    return zero_flag(b, q1), True


def clear_zero_flag(b, w, F):
    ST['scope']['clear_flag'](b, w, F)


_R = {}


def r_poly(n):
    """r(w) = w^n mod prod_{t<n} (w + t): the low coefficients of prod (w + t) (points 0..n-1 of k.vandermonde)."""
    if n not in _R:
        f = _field()
        p = [1]
        for t in range(n):
            q = [0] * (len(p) + 1)
            for i, c in enumerate(p):
                q[i + 1] ^= c
                q[i] ^= f.mul(c, t)
            p = q
        _R[n] = p[:n]
    return _R[n]


def fan_in(b, low, s, r):
    k = ST['dk'].k
    for j, rj in enumerate(r):
        k.cx(b, low[j], s, rj)


# ------------------------------------------------------------------------------------------- the L-level step
def point_products(b, g, low, outs, coefs):
    """outs[i] ^= sum_t coefs[i](p_t) g(p_t) y_t, y_t = <W[:, t], low>: one product per point (kit rank_update)."""
    dk = ST['dk']
    N = len(g)
    V, W = dk.k.vandermonde(N)
    f = _field()
    for t in range(N):
        yc = [W[j][t] for j in range(len(low))]
        if any(yc):
            dk.rank_update(b, g, V[t][:len(g)], low, yc, outs, [c(t, f) for c in coefs])


def erase_acc(b, g, low, acc, coef):
    """Erase acc = sum_t coef(p_t) g(p_t) y_t by measurement; free classically controlled phase products."""
    dk = ST['dk']
    ar = dk.ar
    N = len(g)
    V, W = dk.k.vandermonde(N)
    f = _field()
    m = b.alloc_bits(12)
    b.hmr(acc, m)
    b.free(acc)
    for t in range(N):
        yc = [W[j][t] for j in range(len(low))]
        ct = coef(t, f)
        if not any(yc) or ct == 0:
            continue
        mk = b.alloc_bits(12)
        b.bit_store(mk, False)
        for j in range(12):                       # mk = transpose(multiply by ct) applied to m
            for i in ar.bits(f.mul(ct, 1 << j)):
                b.cx(m[i], mk[j])
        x, px = dk.k.expose(b, g, V[t][:len(g)])
        y, py = dk.k.expose(b, low, yc)
        ar.phase_product(b, x, y, mk)
        dk.k.unexpose(b, low, yc, py)
        dk.k.unexpose(b, g, V[t][:len(g)], px)
        b.free(mk)
    b.free(m)


C_D = lambda t, f: 1 ^ t                         # noqa: E731  d  = o0 + o1
C_E1 = lambda t, f: t                            # noqa: E731  e1 = o1
C_E2 = lambda t, f: t ^ f.mul(t, t)              # noqa: E731  e2 = o1 + o2


def lower_payload(b, L, k, g, v, F):
    """L's return dot for payload k replaced: s ^= S0_eff; X ^= c S1_eff (+ aliasing, closed by cg_part)."""
    ar = _ar()
    s, low = v[-1], v[:-1]
    X = dbg_word('L%d_%d' % (L, k)) if DEBUG else lv(L + 1)['payloads'][k][-1]
    q1 = lv(L)['f'][-1]
    d, e1, e2 = b.alloc_qubits(12), b.alloc_qubits(12), b.alloc_qubits(12)
    point_products(b, g, low, [s, d, e1, e2], [lambda t, f: 1, C_D, C_E1, C_E2])
    for j in range(12):
        b.ccx(F, d[j], s[j])
    erase_acc(b, g, low, d, C_D)
    for j in range(12):
        b.ccx(F, e2[j], e1[j])
    ar.mul_xor0(b, q1, e1, X)                    # X ^= q1' S1sel
    B, own = get_B(b, L, q1)
    for j in range(12):
        b.ccx(B, e1[j], X[j])                    # X ^= B S1sel  (c = q1' + B)
    if own:
        clear_zero_flag(b, q1, B)
    for j in range(12):
        b.ccx(F, e2[j], e1[j])
    erase_acc(b, g, low, e1, C_E1)
    erase_acc(b, g, low, e2, C_E2)
    LOG['lower_ops'] += 1


def cg_part(b, L, g, F):
    """After L's payload scalings: X ^= CG (s + <r, low>), CG = (q1' + B) (F ? g[-3] : g[-2])."""
    ar = _ar()
    q1 = lv(L)['f'][-1]
    b.cswap(F, g[-2], g[-3])
    CG = b.alloc_qubits(12)
    ar.mul_xor0(b, q1, g[-2], CG)
    B, own = get_B(b, L, q1)
    for j in range(12):
        b.ccx(B, g[-2][j], CG[j])
    for k, v in enumerate(lv(L)['payloads']):
        s, low = v[-1], v[:-1]
        X = dbg_word('L%d_%d' % (L, k)) if DEBUG else lv(L + 1)['payloads'][k][-1]
        r = r_poly(len(low))
        fan_in(b, low, s, r)
        ar.mul_xor0(b, CG, s, X)
        fan_in(b, low, s, r)
    for j in range(12):                          # erase CG = q1' g2 + B g2 (B kept live: 869 here)
        b.ccx(B, g[-2][j], CG[j])
    if own:
        clear_zero_flag(b, q1, B)
    ar.erase_product(b, q1, g[-2], CG)
    b.free(CG)
    b.cswap(F, g[-2], g[-3])
    LOG['cg'] += 1


# ------------------------------------------------------------------------------------------- the bottom pass
def lv(m):
    """Level m's (f, g, payloads): recorded at its descent, or derived from level m + 1 (the child's f IS the
    parent's g list, its g the parent's f[:-2], its payloads the parent's v[:-1]: the same qubit words)."""
    c = ST['ctx']
    if m not in c:
        p = lv(m + 1)
        c[m] = dict(f=p['g'], g=p['f'][:-2], payloads=[v[:-1] for v in p['payloads']])
    return c[m]


def tops(m):
    return [v[-1] for v in lv(m)['payloads']]


def dbg_word(name):
    if name not in DBG:
        DBG[name] = ST['b'].alloc_qubits(12)
    return DBG[name]


def route_block_T(b, U, given=None):
    """Transposed carry of block (U, L) plus the moved conversion term: X ^= [not B] q0_L' s_L + ZB s_(U-2)
    + (sA + sB) s_(U-3) + sU2 s_(U-4) (Z1's route_block, kit pd3, transposed).
    given: {m: live flag = [q1_m == 0]} (Z2_PDR_SHARE): used as they are and NOT cleared here."""
    ar = _ar()
    given = given or {}
    c = {m: lv(m) for m in (U, U - 1, U - 2, U - 3, U - 4)}
    L = U - 1
    q = {m: c[m]['f'][-1] for m in (L, U - 2, U - 3, U - 4)}
    fl = {m: given[m] if m in given else zero_flag(b, q[m]) for m in (L, U - 2, U - 3, U - 4)}
    Bq, Cq, Fq, C2q = fl[L], fl[U - 2], fl[U - 3], fl[U - 4]
    X = [dbg_word('R%d_%d' % (U, k)) for k in range(len(tops(U)))] if DEBUG else tops(U)
    Z0 = masked_word(b, Bq, c[L]['f'][-2], neg=True)          # [not B] q0_L'
    for x, s in zip(X, tops(L)):
        ar.mul_xor0(b, Z0, s, x)
    erase_masked_word(b, Bq, c[L]['f'][-2], Z0, neg=True)
    l1 = [(Bq, True), (Cq, False)]
    h1 = flag_and(b, l1)
    ZB = masked_word(b, h1, c[U - 2]['f'][-2])
    l2 = [(Bq, True), (Cq, True), (Fq, True)]
    h2 = flag_and(b, l2)
    b.cx(h2, ZB[0])
    lo = [(Cq, True), (Fq, True)]                             # o = (not C)(not F)
    o = flag_and(b, lo)
    lL = [(Bq, True), (C2q, True), (o, True)]                 # selL = not B, not C2, (C or F)
    lU2 = [(Bq, True), (Cq, False), (C2q, False), (Fq, True)]
    sL, sU2 = flag_and(b, lL), flag_and(b, lU2)
    tu, tl, t2 = tops(U - 2), tops(U - 3), tops(U - 4)
    for k in range(len(X)):
        if 'zb' not in RDROP:                    # Z2_PDR_DROP (diagnostic): leave a route term out
            ar.mul_xor0(b, ZB, tu[k], X[k])
        if 'sl' not in RDROP:
            for j in range(12):
                b.ccx(sL, tl[k][j], X[k][j])
        if 'su2' not in RDROP:
            for j in range(12):
                b.ccx(sU2, t2[k][j], X[k][j])
    erase_flag_and(b, sU2, lU2)
    erase_flag_and(b, sL, lL)
    erase_flag_and(b, o, lo)
    b.cx(h2, ZB[0])
    erase_flag_and(b, h2, l2)
    erase_masked_word(b, h1, c[U - 2]['f'][-2], ZB)
    erase_flag_and(b, h1, l1)
    for m in sorted(fl, reverse=True):
        if m not in given:
            clear_zero_flag(b, q[m], fl[m])


SHARE = os.environ.get('Z2_PDR_SHARE', '0') == '1'
RDROP = set(x for x in os.environ.get('Z2_PDR_DROP', '').split(',') if x)


def route_levels(U):
    return (U - 1, U - 2, U - 3, U - 4)


def bottom_pass_share(b, phase):
    """Z2_PDR_SHARE=1 (second holder): the moved q_posts' handed-up flags F_(m+1) = [q1_m == 0] are exactly the
    route's flags. The first (top) route takes them live; a flag is parked after the last route that needs it, and
    a parked flag that a later route needs is rebuilt once from its kept bit (zero test + CZ, as the pickup does)."""
    sc = ST['scope']
    blocks = sorted(S('upper'), reverse=True)
    first = set(route_levels(blocks[0])) if blocks else set()
    live = {}
    for m in S('moved'):
        d, up = lv(m), lv(m + 1)
        parent = (tops(m + 1), up['f'][-1], up['f'][-2])
        ST['kitq'](b, d['f'][-1], d['payloads'], None, m, parent, False, phase)
        if m + 1 not in sc['KEEP_CHILD']:
            F = phase['returned_gap'].pop(m + 1)
            if m in first:
                live[m] = F
            else:
                ST['pend'][m + 1] = sc['hold_flag'](b, F)
    for i, U in enumerate(blocks):
        for m in route_levels(U):
            if m not in live and m + 1 in ST['pend']:
                F = zero_flag(b, lv(m)['f'][-1])
                sc['consume_flag'](b, F, ST['pend'].pop(m + 1))
                live[m] = F
                LOG['share_rebuilds'] = LOG.get('share_rebuilds', 0) + 1
        route_block_T(b, U, {m: live[m] for m in route_levels(U) if m in live})
        LOG['share_used'] = LOG.get('share_used', 0) + sum(1 for m in route_levels(U) if m in live)
        later = set()
        for V in blocks[i + 1:]:
            later |= set(route_levels(V))
        for m in sorted(live):
            if m not in later:
                ST['pend'][m + 1] = sc['hold_flag'](b, live.pop(m))
    assert not live
    ST['bottom_done'] = True
    LOG['bottom'] += 1


def bottom_pass(b, phase):
    ST['b'] = b
    if SHARE and not QONLY:
        return bottom_pass_share(b, phase)
    for m in S('moved'):
        d, up = lv(m), lv(m + 1)
        parent = (tops(m + 1), up['f'][-1], up['f'][-2])
        ST['kitq'](b, d['f'][-1], d['payloads'], None, m, parent, False, phase)   # the chain's end (exact)
        if m + 1 not in ST['scope']['KEEP_CHILD']:   # the flag handed up to level m + 1: parked until its pickup
            F = phase['returned_gap'].pop(m + 1)
            ST['pend'][m + 1] = ST['scope']['hold_flag'](b, F)
    if not QONLY:
        for U in sorted(S('upper'), reverse=True):
            route_block_T(b, U)
    ST['bottom_done'] = True
    LOG['bottom'] += 1


def bottom_pass_rev_share(b, phase):
    """Z2_PDR_SHARE=1, reverse: each route flag F_(m+1) = [q1_m' == 0] rebuilt once from its parked bit (as the q_post
    inverse would), shared by consecutive routes, parked after the last route that needs it; the last route's flags
    go live into the q_post inverses (descending) instead of a second zero test."""
    sc = ST['scope']
    blocks = sorted(S('upper'))
    live = {}
    for i, U in enumerate(blocks):
        for m in route_levels(U):
            if m not in live:
                F = zero_flag(b, lv(m)['f'][-1])
                if m + 1 in ST['rpend']:
                    sc['consume_flag'](b, F, ST['rpend'].pop(m + 1))
                live[m] = F
        route_block_T(b, U, {m: live[m] for m in route_levels(U)})
        LOG['share_used'] = LOG.get('share_used', 0) + 4
        if i < len(blocks) - 1:
            later = set()
            for V in blocks[i + 1:]:
                later |= set(route_levels(V))
            for m in sorted(live):
                if m not in later:
                    ST['rpend'][m + 1] = sc['hold_flag'](b, live.pop(m))
    for m in reversed(S('moved')):
        d, up = lv(m), lv(m + 1)
        parent = (tops(m + 1), up['f'][-1], up['f'][-2])
        if m in live:
            Fr = live.pop(m)
            LOG['share_rev_live'] = LOG.get('share_rev_live', 0) + 1
        else:
            Fr = zero_flag(b, d['f'][-1])
            if m + 1 in ST['rpend']:
                sc['consume_flag'](b, Fr, ST['rpend'].pop(m + 1))
        phase['given_gap'][m + 1] = Fr
        line_rev = os.environ.get('Z2_PDR_RQ', 'line') == 'line' and m >= sc.get('THRESHOLD', 6)
        qp = ST['origq'][m] if line_rev else ST['kitq']
        qp(b, d['f'][-1], d['payloads'], None, m, parent, True, phase)
    assert not live
    ST['bottom_done'] = True
    LOG['bottom_rev'] = LOG.get('bottom_rev', 0) + 1


def bottom_pass_rev(b, phase):
    """Inverse of bottom_pass: the transposed routes (XORs), then the kit's q_post inverses, descending."""
    if SHARE:
        return bottom_pass_rev_share(b, phase)
    for U in sorted(S('upper')):
        route_block_T(b, U)
    for m in reversed(S('moved')):
        d, up = lv(m), lv(m + 1)
        parent = (tops(m + 1), up['f'][-1], up['f'][-2])
        Fr = zero_flag(b, d['f'][-1])                            # = [q1_m' == 0], the parent's gap flag
        if m + 1 in ST['rpend']:
            ST['scope']['consume_flag'](b, Fr, ST['rpend'].pop(m + 1))
        phase['given_gap'][m + 1] = Fr
        line_rev = os.environ.get('Z2_PDR_RQ', 'line') == 'line' and m >= ST['scope'].get('THRESHOLD', 6)
        qp = ST['origq'][m] if line_rev else ST['kitq']     # direct levels: the line's (cheaper) reverse chain
        qp(b, d['f'][-1], d['payloads'], None, m, parent, True, phase)
    ST['bottom_done'] = True
    LOG['bottom_rev'] = LOG.get('bottom_rev', 0) + 1


# ------------------------------------------------------------------------------------------- hooks
class DKProxy:
    def __init__(self, dk):
        self._dk = dk

    def __getattr__(self, name):
        return getattr(self._dk, name)

    def selected_dot(self, b, g, v, target, F):
        m = len(g) + 1
        if QONLY:
            return self._dk.selected_dot(b, g, v, target, F)
        if (ST['mode'] == 'fwd' and ST['bottom_done']) or (REV and ST['mode'] == 'rev' and not ST['bottom_done']):
            if m in S('upper'):
                LOG['upper_skips'] += 1
                if DEBUG:                            # keep the line's dot; a second copy into a kept word
                    k = [i for i, w in enumerate(lv(m)['payloads']) if w[-1] is target][0]
                    self._dk.selected_dot(b, g, v, dbg_word('D%d_%d' % (m, k)), F)
                    return self._dk.selected_dot(b, g, v, target, F)
                return None
            if m in S('lower'):
                for k, w in enumerate(lv(m)['payloads']):
                    if w[-1] is target:
                        lower_payload(b, m, k, g, list(v) + [target], F)   # v here is the payload's low part
                        if k == len(lv(m)['payloads']) - 1 and ST['mode'] == 'fwd':
                            ST['post'] = dict(L=m, g=list(g) + [None], F=F, last=target)
                        return None
                raise AssertionError('pdr: payload not found at level %d' % m)
        return self._dk.selected_dot(b, g, v, target, F)


def make_scale(orig):
    def scale(b, x, a, ai, divide=True):
        cur = ST['cur']
        if (REV and ST['mode'] == 'rev' and not ST['bottom_done'] and cur is not None and cur['reverse']
                and cur['m'] in S('lower') and not divide and not cur['cg_done']
                and x is lv(cur['m'])['payloads'][0][-1]):
            cur['cg_done'] = True
            cg_part(b, cur['m'], lv(cur['m'])['g'], cur['F'])
        out = orig(b, x, a, ai, divide)
        p = ST['post']
        if p is not None and ST['mode'] == 'fwd' and x is p['last'] and divide:
            ST['post'] = None
            L = p['L']
            gfull = lv(L)['g']
            cg_part(b, L, gfull, p['F'])
        return out
    scale._z2_pdr = True
    return scale


def make_qpost(orig, kind):
    def q_post(b, q, payloads, H, m, parent, reverse=False, phase=None):
        if SY['mode'] == 'capture' and not reverse and phase is not None:     # record the line's handed-up wire
            out = orig(b, q, payloads, H, m, parent, reverse, phase)
            F = phase['returned_gap'].get(m + 1)
            if F is not None:
                HOME_WIRE[m] = int(str(F)[1:])
            return out
        if ST['mode'] == 'fwd' and not reverse and m in S('moved'):
            if not ST['bottom_done']:
                bottom_pass(b, phase)
            LOG['qpost_skips'] += 1
            return None
        if REV and ST['mode'] == 'rev' and reverse and m in S('moved'):
            F = phase['given_gap'].pop(m + 1)       # parked here (free), rebuilt once at the bottom
            if os.environ.get('Z2_PDR_RPARK', '1') == '1':
                ST['rpend'][m + 1] = ST['scope']['hold_flag'](b, F)
            else:
                clear_zero_flag(b, lv(m)['f'][-1], F)
            LOG['rev_book'] = LOG.get('rev_book', 0) + 1
            if m == min(S('moved')):
                bottom_pass_rev(b, phase)
            return None
        return orig(b, q, payloads, H, m, parent, reverse, phase)
    q_post._z2_pdr = kind
    return q_post


def make_pre(orig, kind):
    def pre(b, f, g, payloads, *a, **kw):
        if kind == 'direct':
            undo = kw.get('undo', a[1] if len(a) > 1 else False)
        elif kind == 'kit':
            undo = kw.get('undo', a[0] if a else False)
        else:
            undo = False                                       # pre_partial: forward only
        if ST['mode'] == 'fwd' and not undo:
            m = len(f) - 1
            ST['ctx'][m] = dict(f=f, g=g, payloads=payloads)
            LOG['ctx'] += 1
        return orig(b, f, g, payloads, *a, **kw)
    pre._z2_pdr = kind
    return pre


HOME_WIRE = {}      # Z2_PDR_SYNC capture: level m -> wire of the flag the line's own q_post(m) hands up (forward)


def pickup(b, m):
    """After level m's forward post: rebuild the parked flag of level m + 1 (= [f_m[-1] == 0]) for the parent.
    With Z2_PDR_SYNC=1 it is rebuilt on the wire the line's own q_post(m) gives it (recorded in the capture), so a
    flag that the line then holds across the scan sits where the line holds it."""
    if ST['mode'] == 'fwd' and ST['bottom_done'] and m + 1 in ST['pend']:
        sc = ST['scope']
        F = None
        if SYNC and m in HOME_WIRE:
            pp2 = sys.modules.get('dqi_int_polepark2')
            fpk = sys.modules.get('dqi_int_flagpark')
            if pp2 is not None and fpk is not None:
                F = pp2.alloc_pref(b, HOME_WIRE[m])
                fpk.zero_toggle(b, lv(m)['f'][-1], F)
                LOG['home'] = LOG.get('home', 0) + int(int(str(F)[1:]) == HOME_WIRE[m])
        if F is None:
            F = zero_flag(b, lv(m)['f'][-1])
        sc['consume_flag'](b, F, ST['pend'].pop(m + 1))
        LOG['pickups'] = LOG.get('pickups', 0) + 1
        return F
    return None


def make_post(orig):
    def direct_post(b, f, g, payloads, F, *a, **kw):
        reverse = kw.get('reverse', a[2] if len(a) > 2 else False)
        phase = kw.get('phase', a[5] if len(a) > 5 else None)
        if phase is not None:
            ST['phase'] = phase
        prev = ST['cur']
        ST['cur'] = dict(m=len(f) - 1, F=F, reverse=bool(reverse), cg_done=False)
        try:
            out = orig(b, f, g, payloads, F, *a, **kw)
        finally:
            ST['cur'] = prev
        Fr = pickup(b, len(f) - 1)
        if Fr is not None:
            phase['returned_gap'][len(f)] = Fr
        return out
    direct_post._z2_pdr = True
    return direct_post


def make_post_look(orig):
    def post_look(b, f, g, payloads, H, parent, *a, **kw):
        reverse = kw.get('reverse', a[0] if a else False)
        phase = kw.get('phase', a[3] if len(a) > 3 else None)
        out = orig(b, f, g, payloads, H, parent, *a, **kw)
        if not reverse:
            Fr = pickup(b, len(f) - 1)
            if Fr is not None:
                phase['returned_gap'][len(f)] = Fr
        return out
    post_look._z2_pdr = True
    return post_look


CANON = os.environ.get('Z2_PDR_CANON', '0') == '1'
CANON_LOW = int(os.environ.get('Z2_PDR_CANON_LOW', '840'))


def canon_free(b):
    """Z2_PDR_CANON=1 (second holder): the circuit builder hands out the LAST freed wire first, so a patched converter's
    free order moves every later allocation, and the decode's final route_words needs the syndrome wires (< 840) free
    at the end (the layout trap). At the converter's exit the free wires below the high-water mark are re-stacked:
    syndrome wires at the bottom, every other wire above them, each group in its own order. No gate, no new wire."""
    spare = b.finish_circuit().num_qubits - b.num_allocated_qubits
    if spare <= 0:
        return
    free_now = list(b.alloc_qubits(spare))           # top of the stack first
    low = [q for q in free_now if int(str(q)[1:]) < CANON_LOW]
    high = [q for q in free_now if int(str(q)[1:]) >= CANON_LOW]
    for q in reversed(low):
        b.free(q)
    for q in reversed(high):
        b.free(q)
    LOG['canon'] = LOG.get('canon', 0) + 1


SYNC = os.environ.get('Z2_PDR_SYNC', '0') == '1'
SY = dict(mode='off', rec=[], i=0, fixes=0, misses=0)


def sync_point(b):
    """Z2_PDR_SYNC=1 (second holder): the free stack at every converter exit as it would be WITHOUT pdr. Capture mode
    (one scratch decode at install, pdr's blocks emptied) records the order; replay re-stacks to it when the free set
    is the same (B7's sync3 method, recorded after every other module). kickmix hands out the last freed wire first,
    so the rest of the decode then allocates exactly as without pdr: pdr cannot move a later stage's wires."""
    n = b.finish_circuit().num_qubits - b.num_allocated_qubits
    got = list(b.alloc_qubits(n)) if n > 0 else []
    if SY['mode'] == 'capture':
        SY['rec'].append([int(str(q)[1:]) for q in got])
    elif SY['mode'] == 'replay' and SY['i'] < len(SY['rec']):
        rec = SY['rec'][SY['i']]
        SY['i'] += 1
        have = {int(str(q)[1:]): q for q in got}
        if set(have) == set(rec):
            for i in reversed(rec):
                b.free(have[i])
            SY['fixes'] += 1
            print('pdr sync restored exit %d (%d free wires)' % (SY['i'] - 1, len(rec)), flush=True)
            return
        SY['misses'] += 1
        print('pdr sync MISS at exit %d: free now, not recorded %s; recorded, not free %s' % (
            SY['i'] - 1, sorted(set(have) - set(rec))[:8], sorted(set(rec) - set(have))[:8]), flush=True)
    for q in reversed(got):
        b.free(q)


def sync_capture(ns):
    """One scratch decode with pdr's block sets emptied (every hook passes through), recording the exits' free order."""
    km_ = sys.modules.get('kickmix') or __import__('kickmix')
    CB = km_.CircuitBuilder
    saved_sets = {d: dict(SETS[d]) for d in SETS}
    for d in SETS:
        SETS[d].update(upper=(), lower=(), moved=())
    others = {}
    a3 = sys.modules.get('dqi_pub_a3')
    if a3 is not None and hasattr(a3, 'STACK'):
        others['a3'] = dict(a3.STACK)
    s3 = sys.modules.get('dqi_var_sync3')
    if s3 is not None:
        others['s3'] = dict(s3.ST)
    tagged = CB.mark
    fv = getattr(getattr(tagged, '__code__', None), 'co_freevars', ())
    if 'done' in fv and 'original' in fv:                    # do not spend the re-roll tag (as D1 and B7 do)
        CB.mark = dict(zip(fv, tagged.__closure__))['original'].cell_contents
    SY.update(mode='capture', rec=[], i=0)
    try:
        b0 = CB()
        S0 = [b0.alloc_qubits(12) for _ in range(70)]
        E0 = [b0.alloc_bits(12) for _ in range(4095)]
        ns.decode(b0, S0, E0)
    finally:
        CB.mark = tagged
        for d in SETS:
            SETS[d].update(saved_sets[d])
        if 'a3' in others:
            a3.STACK.clear()
            a3.STACK.update(others['a3'])
        if 's3' in others:
            s3.ST.clear()
            s3.ST.update(others['s3'])
        for k in ('pend', 'rpend', 'ctx', 'origq_dummy'):
            if isinstance(ST.get(k), dict):
                ST[k].clear()
        ST.update(mode='-', bottom_done=False, post=None, cur=None)
    SY.update(mode='replay', i=0)
    return len(SY['rec'])


def make_conv(orig):
    def convert_safe(b, A, B, reverse=False, *a, **kw):
        ST.update(mode='rev' if reverse else 'fwd', bottom_done=False, post=None)
        assert not ST['pend'], 'pdr: parked flags left over %s' % sorted(ST['pend'])
        if not reverse:
            ST['ctx'] = {}
            DBG.clear()
        try:
            out = orig(b, A, B, reverse, *a, **kw)
            if CANON:
                canon_free(b)
            if SYNC:
                sync_point(b)
            return out
        finally:
            ST['mode'] = '-'
    convert_safe._z2_pdr = True
    return convert_safe


def install(km, ns):
    walk = sys.modules['walk_bound']
    done = set()
    for fn in walk.functions(ns.decode):
        gl = fn.__globals__
        if id(gl) in done:
            continue
        if fn.__name__ == 'direct_post' and fn.__code__.co_filename.endswith('high_stream.py'):
            done.add(id(gl))
            need = ('dk', 'scale', 'q_post', 'direct_pre', 'fresh_zero_flag', 'clear_flag', 'hold_flag',
                    'consume_flag', 'KEEP_CHILD')
            if not all(k in gl for k in need):
                LOG.setdefault('skipped', []).append([k for k in need if k not in gl])
                continue
            ST['scope'] = ST['scope'] or gl
            if not isinstance(gl['dk'], DKProxy):
                gl['dk'] = DKProxy(gl['dk'])
            ST['dk'] = ST['dk'] or gl['dk']._dk
            if not getattr(gl['scale'], '_z2_pdr', False):
                gl['scale'] = make_scale(gl['scale'])
            if not getattr(gl['q_post'], '_z2_pdr', False):
                orig = gl['q_post']
                LOG.setdefault('qchain', []).append(('direct', orig.__code__.co_filename.split('/')[-1],
                                                     orig.__code__.co_firstlineno))
                gl['q_post'] = make_qpost(orig, 'direct')
                for m in ALL_MOVED:
                    if m >= gl.get('THRESHOLD', 6):
                        ST['origq'][m] = orig
            if not getattr(gl['direct_pre'], '_z2_pdr', False):
                gl['direct_pre'] = make_pre(gl['direct_pre'], 'direct')
            if 'pre_partial' in gl and not getattr(gl['pre_partial'], '_z2_pdr', False):
                gl['pre_partial'] = make_pre(gl['pre_partial'], 'partial')
            if not getattr(gl['direct_post'], '_z2_pdr', False):
                gl['direct_post'] = make_post(gl['direct_post'])
            LOG['scopes'] += 1
    for fn in walk.functions(ns.decode):
        gl = fn.__globals__
        if fn.__name__ == 'post_look' and 'q_post' in gl and id(gl) not in done:
            done.add(id(gl))
            if not getattr(gl['q_post'], '_z2_pdr', False):
                orig = gl['q_post']
                LOG.setdefault('qchain', []).append(('kit', orig.__code__.co_filename.split('/')[-1],
                                                     orig.__code__.co_firstlineno))
                gl['q_post'] = make_qpost(orig, 'kit')
                for m in ALL_MOVED:
                    if m < 6:
                        ST['origq'].setdefault(m, orig)
            if 'pre' in gl and not getattr(gl['pre'], '_z2_pdr', False):
                gl['pre'] = make_pre(gl['pre'], 'kit')
            if not getattr(gl['post_look'], '_z2_pdr', False):
                gl['post_look'] = make_post_look(gl['post_look'])
            LOG['kit_scopes'] += 1
    for fn in walk.functions(ns.decode):
        if fn.__name__ == 'q_post' and fn.__code__.co_filename.endswith('scalar_compact_v3.py') and fn.__closure__:
            cells = dict(zip(fn.__code__.co_freevars, fn.__closure__))
            if 'original' in cells:
                ST['kitq'] = cells['original'].cell_contents
    assert ST['kitq'] is not None, 'pdr: kit q_post not found'
    LOG['kitq'] = (ST['kitq'].__code__.co_filename.split('/')[-1], ST['kitq'].__code__.co_firstlineno)
    g = ns.decode.__globals__['g']
    if not getattr(g['convert_safe'], '_z2_pdr', False):
        g['convert_safe'] = make_conv(g['convert_safe'])
    missing = [m for m in ALL_MOVED if m not in ST['origq']]
    assert not missing, 'pdr: no q_post owner for levels %s' % missing
    if SYNC:
        LOG['sync_rec'] = sync_capture(ns)
    sc = ST['scope']
    return dict(LOG, upper_fwd=UPPER_F, upper_rev=UPPER_R if REV else (), moved=ALL_MOVED,
                keep_child=sorted(sc.get('KEEP_CHILD', ())), keep_full=sorted(sc.get('KEEP_FULL', ())))
