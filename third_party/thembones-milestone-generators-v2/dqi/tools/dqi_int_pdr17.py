#!/usr/bin/env python3
"""dqi_int_pdr17.py -- B5 (13th holder, Nov 16): the line's pdr (snap_b5r46n dqi_z2_pdr.py 690ce6dd, byte
for byte below except the marked B5 edits) + THE ROOT BLOCK (U = 17, L = 16) paired return (R48-DQI-CONV2 attack 1).
Forward: the root inverse lives on through level 16's descent (its retirement deferred), the root digit is scaled
back to its input value right after level 16's descent (T = clip(q1_16) * level 16's pivot = the root pivot), the
root's return dots and its first scaling are dropped, term 1 ([not B] q0_16 s16pre) is read at level 16's post after
root_child_post, term 3 and the gap routes at the bottom (route_block_T(17) without its Z0 row), term 2 by lower_payload.
Original header follows.
dqi_z2_pdr.py -- builder Z2 (Nov 14): the paired RETURN of the converter, installed into the LINE (r12hz's
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
    _hi = min(max(_u), 15) if _u else 0                     # B5: levels 16, 17 never move
    SETS[_d] = dict(upper=_u, lower=tuple(u - 1 for u in _u), moved=tuple(range(_lo, _hi + 1)) if _u else ())
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
    if LIVEB and ph is not None and L != 16:                 # B5: kpark parks Froot through level 16's post
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
    if L == 16 and k == 0 and ST['mode'] == 'fwd':           # B5: term 1 of the root block, s16 still pre here
        term1_17(b, q1)
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
    if L == 16 and ST['mode'] == 'rev' and k == len(lv(16)['payloads']) - 1:   # B5: term 1's undo
        term1_17(b, q1)


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
    if L == 16:
        check(b, 'after cg_part(16)')
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
    if U != 17:                                               # B5: block 17's Z0 row runs at level 16's post
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
    NOSU2 = os.environ.get('R46_NOSU2') == '1'
    sL = flag_and(b, lL)
    sU2 = None if NOSU2 else flag_and(b, lU2)
    tu, tl, t2 = tops(U - 2), tops(U - 3), tops(U - 4)
    for k in range(len(X)):
        if 'zb' not in RDROP:                    # Z2_PDR_DROP (diagnostic): leave a route term out
            ar.mul_xor0(b, ZB, tu[k], X[k])
        if 'sl' not in RDROP:
            for j in range(12):
                b.ccx(sL, tl[k][j], X[k][j])
        if 'su2' not in RDROP and not NOSU2:
            for j in range(12):
                b.ccx(sU2, t2[k][j], X[k][j])
    if not NOSU2:
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
    if U == 17:
        check(b, 'after route_block_T(17)')


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
    if blocks and blocks[0] == 17 and len(blocks) > 1:        # B5
        first |= set(route_levels(blocks[1]))
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
    for m in sorted([x for x in live if x not in S('moved')], reverse=True):   # B5: level 16 (block 17)
        clear_zero_flag(b, lv(m)['f'][-1], live.pop(m))
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
        if ST['mode'] == 'fwd' and not reverse and ST.get('desc') in ('l16', 'below'):
            ST['desc'] = 'done'                                  # B5 diagnostics: the descent ends here
        if m == 16 and ST['mode'] == 'fwd' and not reverse and P17F:
            check(b, 'before q_post(16)')
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
        if kind == 'direct' and len(f) - 1 == 16 and ST['mode'] == 'rev' and undo and P17R:
            unearly17(b, f, g, a[0])                             # B5
        out = orig(b, f, g, payloads, *a, **kw)
        if kind == 'direct' and len(f) - 1 == 16 and ST['mode'] == 'fwd' and not undo and P17F:
            early17(b, f, g, a[0])                               # B5: a[0] = level 16's gap flag F
        return out
    pre._z2_pdr = kind
    return pre


def pickup(b, m):
    """After level m's forward post: rebuild the parked flag of level m + 1 (= [f_m[-1] == 0]) for the parent."""
    if ST['mode'] == 'fwd' and ST['bottom_done'] and m + 1 in ST['pend']:
        sc = ST['scope']
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
        if len(f) - 1 == 16 and P17F and not reverse:
            check(b, 'level 16 post start')
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


def make_conv(orig):
    def convert_safe(b, A, B, reverse=False, *a, **kw):
        ST.update(mode='rev' if reverse else 'fwd', bottom_done=False, post=None)
        assert not ST['pend'], 'pdr: parked flags left over %s' % sorted(ST['pend'])
        if not reverse:
            ST['ctx'] = {}
            DBG.clear()
        try:
            return orig(b, A, B, reverse, *a, **kw)
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
    if P17F or P17R:
        install_root(km, ns)                                      # B5: after pdr's own convert_safe wrapper
    missing = [m for m in ALL_MOVED if m not in ST['origq']]
    assert not missing, 'pdr: no q_post owner for levels %s' % missing
    sc = ST['scope']
    return dict(LOG, upper_fwd=UPPER_F, upper_rev=UPPER_R if REV else (), moved=ALL_MOVED,
                keep_child=sorted(sc.get('KEEP_CHILD', ())), keep_full=sorted(sc.get('KEEP_FULL', ())))


# ============================================================================ B5: the root block (U = 17, L = 16)
# Forward only in this version (Z2_PDR_BLOCKS holds 17); the reverse stays the line's unless 17 is in Z2_PDR_RBLOCKS.
P17F = 17 in UPPER_F
P17R = REV and 17 in UPPER_R
LOG17 = dict(early=0, term1=0, root_dots_skipped=0, root_scales_skipped=0, deferred=0, retired=0)
_CB = {}
DESC = {'l16': (0, ''), 'below': (0, '')}


CHECK = os.environ.get('B5_P17_CHECK', '0') == '1'      # diagnostics: validate the circuit after each piece
NODEFER = os.environ.get('B5_P17_NODEFER', '0') == '1'  # diagnostics: no deferral, no early scaling (wrong words)
# 14th holder, diagnostics with Z2_PDR_MODE=debug: every pdr term goes to a kept word (T17_k term 1, L16_k lower terms,
# R17_k routes) and the root keeps the line's dots with a copy in D17_k; by default the line's root scalings too
# (early17 off). B5_P17_DBGSCALE=1 keeps early17 and the first-scale skip (tests T against the root's clipped pivot).
DBGSCALE = os.environ.get('B5_P17_DBGSCALE', '0') == '1'
SKIP_FIRST = (not DEBUG) or DBGSCALE                    # the root post's first scaling is dropped


def _p17_scaling():
    """True when the deferral of the root inverse and early17 run."""
    return (not NODEFER) and SKIP_FIRST


def check(b, where):
    if not CHECK:
        return
    try:
        b.finish_circuit()
        print('pdr17 check ok', where, ST['mode'], flush=True)
    except Exception as e:                                   # noqa: BLE001
        print('pdr17 check FAILED', where, ST['mode'], repr(e)[:120], flush=True)
        raise


def _same(x, w):
    return x is w or (len(x) == len(w) and len(x) > 0 and x[0] is w[0])


def _kit_scale():
    from exp.safe_scalar import converter as SS
    return SS.scale


def term1_17(b, q1):
    """X17 ^= [not B] q0_16' s16pre for every payload: Z2's Z0 row of route_block_T, read at level 16's post right
    after root_child_post (the only moment where digit 16 is in its pre state); B = root_child_post's flag (live)."""
    ar = _ar()
    check(b, 'term1 start')
    B, own = get_B(b, 16, q1)
    q0 = lv(16)['f'][-2]
    if CHECK:
        print('pdr17 ids B', B, 'own', own, 'q0', q0[0], q0[-1], 'q1', q1[0],
              [(v16[-1][0], v16[-1][-1], v17[-1][0], v17[-1][-1]) for v16, v17 in
               zip(lv(16)['payloads'], lv(17)['payloads'])], flush=True)
    Z0 = masked_word(b, B, q0, neg=True)
    if CHECK:
        print('pdr17 ids Z0', Z0[0], Z0[-1], flush=True)
    check(b, 'term1 after masked word')
    for k, (v16, v17) in enumerate(zip(lv(16)['payloads'], lv(17)['payloads'])):
        ar.mul_xor0(b, Z0, v16[-1], dbg_word('T17_%d' % k) if DEBUG else v17[-1])
        check(b, 'term1 after a product')
    erase_masked_word(b, B, q0, Z0, neg=True)
    if own:
        clear_zero_flag(b, q1, B)
    check(b, 'after term1')
    LOG17['term1'] += 1


def early17(b, f, g, F):
    """After level 16's forward descent (Z1's paired kernel has read the root digit d17 = X ai17): X17 <- d17 T with
    T = clip(q1_16) * (level 16's pivot after its own pivot swap) = the root's pivot (generic: q1_16 = a17 / a16;
    root gap: q1_16 = 0 and level 16's pivot is the root's pivoted word). The old X17 is erased with the root inverse
    ai17 (free phase); ai17 is then retired here (the line retires it before level 16: that hmr was deferred)."""
    check(b, 'before early17')
    if not _p17_scaling():
        return
    d = ST.get('defer')
    assert d is not None and d['state'] == 'held', 'pdr17: root inverse not held at level 16 (%s)' % (d and d['state'])
    ar = _ar()
    sc = ST['scope']
    q1 = f[-1]
    B = sc['fresh_zero_flag'](b, q1)
    b.cx(B, q1[0])                                          # clip: q1' = q1 + [q1 == 0]
    b.cswap(F, g[-1], g[-2])                                # level 16's pivot
    T = b.alloc_qubits(12)
    ar.mul_xor0(b, q1, g[-1], T)
    ai = d['ai']
    scale = _kit_scale()
    for v in lv(17)['payloads']:
        scale(b, v[-1], T, ai, False)                       # X17 = d17 T; the old X17 = new * ai17, erased free
    ar.erase_product(b, q1, g[-1], T)
    b.free(T)
    b.cswap(F, g[-1], g[-2])
    b.cx(B, q1[0])
    sc['clear_flag'](b, q1, B)
    _CB['hmr'](b, ai, d['bits'])                            # the deferred retirement of the root inverse
    _CB['free'](b, ai)
    d['state'] = 'done'
    if ST.get('desc') == 'l16':
        ST['desc'] = 'below'
    check(b, 'after early17')
    LOG17['early'] += 1
    LOG17['retired'] += 1


def unearly17(b, f, g, F):
    """Reverse, at the START of level 16's descent undo (before Z1's reverse kernel reads d17): the root inverse,
    retired by a12 after level 16's reverse post, is rebuilt as 1 / T (T as in early17) into the same word, X17 <- X17
    ai17 = d17 (the old X17 = new T erased free), T erased; ai17 then stays live to the root's end, where a12's own
    rebuild is replaced by a hand-over (OldProxy._inverse_zero)."""
    ar = _ar()
    sc = ST['scope']
    held, _mask = ST['phase']['_a12_reverse_root']
    q1 = f[-1]
    B = sc['fresh_zero_flag'](b, q1)
    b.cx(B, q1[0])
    b.cswap(F, g[-1], g[-2])
    T = b.alloc_qubits(12)
    ar.mul_xor0(b, q1, g[-1], T)
    inv = sc.get('inverse_initialized')
    if inv is None:
        from exp.inverse19_joint.kernel import inverse_initialized as inv
    new = inv(b, T)
    held[:] = new
    scale = _kit_scale()
    for v in lv(17)['payloads']:
        scale(b, v[-1], T, held, True)                      # X17 = X ai17 = d17; the old X17 = new T, erased free
    ar.erase_product(b, q1, g[-1], T)
    b.free(T)
    b.cswap(F, g[-1], g[-2])
    b.cx(B, q1[0])
    sc['clear_flag'](b, q1, B)
    ST['rev_ai17'] = held
    LOG17['unearly'] = LOG17.get('unearly', 0) + 1


class RootDK:
    """The root scope's dk: the root's return dots are dropped in the forward when block 17 is paired."""
    def __init__(self, dk):
        self._dk = dk

    def __getattr__(self, name):
        return getattr(self._dk, name)

    def selected_dot(self, b, g, v, target, F):
        if len(g) + 1 == 17 and ((P17F and ST['mode'] == 'fwd' and ST['bottom_done'])
                                 or (P17R and ST['mode'] == 'rev' and not ST['bottom_done'])):
            if DEBUG:                                        # keep the line's dot; a copy into a kept word
                k = [i for i, w in enumerate(lv(17)['payloads']) if _same(w[-1], target)][0]
                self._dk.selected_dot(b, g, v, dbg_word('D17_%d' % k), F)
                LOG17['root_dots_debug'] = LOG17.get('root_dots_debug', 0) + 1
                return self._dk.selected_dot(b, g, v, target, F)
            LOG17['root_dots_skipped'] += 1
            return None
        return self._dk.selected_dot(b, g, v, target, F)


def make_old_scale(orig):
    """The root scope's old.scale (wrapped in place on the namespace: the harness and walk_bound reach a11's hankel
    through it): the forward drops the root post's first scaling (X = d a); the reverse drops its twin (X / a)."""
    def scale(b, x, a, ai, divide=True):
        top = 17 in ST['ctx'] and any(_same(x, v[-1]) for v in ST['ctx'][17]['payloads'])
        if top and SKIP_FIRST and ((P17F and ST['mode'] == 'fwd' and ST['bottom_done'] and not divide)
                                   or (P17R and ST['mode'] == 'rev' and not ST['bottom_done'] and divide)):
            LOG17['root_scales_skipped'] += 1
            return None
        return orig(b, x, a, ai, divide)
    scale._b5_p17 = True
    return scale


def make_old_inverse_zero(orig):
    """Reverse, the root's end (a12's rebuild of the root inverse): when it is already live (rebuilt at level 16's
    reverse descent), hand it over with a fresh zero flag of a (erased by the line's erase_product_zero)."""
    def inverse_zero(b, a, *r, **kw):
        held = ST.get('rev_ai17')
        if held is None:
            return orig(b, a, *r, **kw)
        ST['rev_ai17'] = None
        Z = ST['scope']['fresh_zero_flag'](b, a)
        LOG17['rev_handover'] = LOG17.get('rev_handover', 0) + 1
        return list(held), Z
    inverse_zero._b5_p17 = True
    return inverse_zero


def make_ordinary(orig):
    def ordinary_root_action(b, f, g, payloads, ai, Z, F, post, undo, *rest, **kw):
        if P17F and ST['mode'] == 'fwd' and not post and not undo:
            ST['ctx'][17] = dict(f=f, g=g, payloads=payloads)
            if _p17_scaling():
                ST['defer'] = dict(ai=f[-1], state='armed', bits=None, freed=False)
            ST['desc'] = 'l16'
            LOG17['deferred'] += 1
            check(b, 'root pre start')
        out = orig(b, f, g, payloads, ai, Z, F, post, undo, *rest, **kw)
        check(b, 'after ordinary root action post=%s undo=%s' % (post, undo))
        return out
    ordinary_root_action._b5_p17 = True
    return ordinary_root_action


def _builder_class():
    return sys.modules['kickmix'].CircuitBuilder          # as r46_dqi_conv_zpark finds it


def install_root(km, ns):
    """Wrap root_fusion's scope (dk, old, ordinary_root_action) and the builder's hmr / free (deferred retirement)."""
    walk = sys.modules['walk_bound']
    n = 0
    seen = set()
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
            # root_fusion's cloned scope ('short'): it holds the root's restore and ordinary actions, old and dk
            if not all(k in gl for k in ('ordinary_root_action', 'old', 'dk', 'root_restore', 'root_action')):
                continue
            if gl.get('_b5_p17_done'):
                continue
            gl['_b5_p17_done'] = True
            gl['dk'] = RootDK(gl['dk'])
            o = gl['old']
            if not getattr(getattr(o, 'scale', None), '_b5_p17', False):
                o.scale = make_old_scale(o.scale)
                if hasattr(o, 'inverse_zero'):
                    o.inverse_zero = make_old_inverse_zero(o.inverse_zero)
                LOG17['old_wraps'] = LOG17.get('old_wraps', 0) + 1
            if not getattr(gl['ordinary_root_action'], '_b5_p17', False):
                gl['ordinary_root_action'] = make_ordinary(gl['ordinary_root_action'])
            n += 1
    LOG17['root_scopes'] = n
    assert n >= 1, 'pdr17: root_fusion scope not found'
    CB = _builder_class()
    if not getattr(CB.hmr, '_b5_p17', False):
        o_hmr, o_free = CB.hmr, CB.free
        _CB.update(hmr=o_hmr, free=o_free)

        def hmr(self, q, bits, *a, **kw):
            d = ST.get('defer')
            if d is not None and d['state'] == 'armed' and q is d['ai']:
                d['state'] = 'held'
                d['bits'] = bits
                return None
            return o_hmr(self, q, bits, *a, **kw)

        def free(self, x, *a, **kw):
            d = ST.get('defer')
            if d is not None and d['state'] == 'held' and x is d['ai'] and not d['freed']:
                d['freed'] = True
                return None
            return o_free(self, x, *a, **kw)
        hmr._b5_p17 = True
        free._b5_p17 = True
        CB.hmr, CB.free = hmr, free
    if os.environ.get('B5_P17_DESC') == '1' and not getattr(CB.alloc_qubits, '_b5_desc', False):
        o_alloc3 = CB.alloc_qubits                   # diagnostics: live high-water of the forward descent

        def alloc_qubits3(self, n, *a, **kw):
            r = o_alloc3(self, n, *a, **kw)
            ph = ST.get('desc')
            if ST['mode'] == 'fwd' and ph in ('l16', 'below'):
                live = self.num_allocated_qubits
                if live > DESC[ph][0]:
                    fr, out = sys._getframe(1), []
                    while fr is not None and len(out) < 6:
                        out.append(fr.f_code.co_name)
                        fr = fr.f_back
                    DESC[ph] = (live, '<'.join(out))
            return r
        alloc_qubits3._b5_desc = True
        CB.alloc_qubits = alloc_qubits3
        atexit_desc = True
        import atexit as _ax
        _ax.register(lambda: print('pdr17 descent high-water', DESC, flush=True))
    W = os.environ.get('B5_P17_WATCH')
    if W and not getattr(CB.alloc_qubits, '_b5_watch', False):       # diagnostics: who frees / allocates qubit W
        wid = int(W)
        o_alloc2, o_free2 = CB.alloc_qubits, CB.free

        def _stack():
            fr, out = sys._getframe(2), []
            while fr is not None and len(out) < 10:
                out.append('%s:%d' % (fr.f_code.co_name, fr.f_lineno))
                fr = fr.f_back
            return ' < '.join(out)

        def _has(x):
            for q in (x if isinstance(x, (list, tuple)) else [x]):
                if getattr(q, 'id', None) == wid:
                    return True
            return False

        def alloc_qubits(self, n, *a, **kw):
            r = o_alloc2(self, n, *a, **kw)
            if ST['mode'] == 'fwd' and _has(r):
                print('pdr17 watch ALLOC q%d n=%d %s' % (wid, n, _stack()), flush=True)
            return r

        def free2(self, x, *a, **kw):
            if ST['mode'] == 'fwd' and _has(x):
                print('pdr17 watch FREE q%d %s' % (wid, _stack()), flush=True)
            return o_free2(self, x, *a, **kw)
        alloc_qubits._b5_watch = True
        CB.alloc_qubits, CB.free = alloc_qubits, free2
    g = ns.decode.__globals__['g']
    inner = g['convert_safe']

    def convert_safe(b, A, B, reverse=False, *a, **kw):
        ST['defer'] = None
        ST['rev_ai17'] = None
        out = inner(b, A, B, reverse, *a, **kw)
        if P17R and reverse:
            assert ST.get('rev_ai17') is None and LOG17.get('unearly', 0) == LOG17.get('rev_handover', 0),                 'pdr17: reverse root inverse not handed over'
        if P17F and not reverse and _p17_scaling():
            d = ST.get('defer')
            assert d is not None and d['state'] == 'done', 'pdr17: root inverse retirement not done (%s)' % (d,)
        return out
    convert_safe._b5_p17 = True
    g['convert_safe'] = convert_safe
    import atexit
    atexit.register(lambda: print('dqi_int_pdr17 root', LOG17, flush=True))
    return n
