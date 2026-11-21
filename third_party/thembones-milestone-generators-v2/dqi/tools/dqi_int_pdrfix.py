"""dqi_int_pdrfix.py -- workstream dqi_int (B5, Nov 14): repair of the DOUBLE-GAP class of Z2's paired return (pdr).

The class [M, my predecessor, harness 13:45Z]: when two consecutive quotient heads are zero (B = [q1_L == 0] and
C = [q1_(U-2) == 0] both 1, U = L + 1; the remainder degree drops by 3), U's return dot equals s_L, which pdr's
route never adds; and pdr's L-level terms (lower_payload's X ^= B S1sel, cg_part's CG = (q1' + B) g2) must then carry
B.not(C) instead of B. Law measured at block 14 on 16 / 16 synthetic charts; this module is the repair:
  route_block_T:  hD = B C (one AND) is XORed into bit 0 of Z0 = [not B] q0_L' around Z0's products, so Z0 = 1 on the
                  double gap and X ^= s_L rides the existing product (Z0 is 0 whenever B = 1, so nothing else moves).
  lower_payload:  X ^= (B not F) S1sel, F = L's own gap flag (= [q1_(L-1) == 0] = C, live: it is the dot's F).
  cg_part:        CG = (q1' + B not F) (F ? g[-3] : g[-2]), the same flag for the term and its erasure.
Price: 1 AND per route, 1 AND per payload in lower_payload, 1 AND in cg_part: 4 per block-direction with two payloads
(the erasures are by measurement). Width: one more flag qubit at each of those moments.
v2 (16:50Z, after the harness h16 measured that a double gap breaks TWO blocks): the double gap with zero heads at
(a, a - 1) meets a block as (B, C), (C, F) or (F, C2). Laws fitted on 32 charts each [M, tools/dqi_int_lawfit.py]:
(F, C2) with B = C = 0 adds s_(U-4); (C, F) with B = C2 = 0 needs q0_(U-3)' s_(U-3) + s_(U-4) in place of s_(U-3)
(a third product: NOT repaired here). PDRFIX_ROUTE=1 (default) rewrites Z2's selectors exactly and cheaper
(h2 = not B, o; sU2 = not B, C2, C xor F: -4 ANDs per route) and adds the (F, C2) law; PDRFIX_BC=1 (default) the
(B, C) repair above. Both: about 0 per block-direction, 4 of the 7 zero-head pairs (14,13 / 12,11 / 10,9 / 8,7) exact.
Install right AFTER dqi_z2_pdr (DQI_INT_EXT=...,dqi_z2_pdr,dqi_int_pdrfix): it replaces the three functions in the
loaded pdr module (pdr's hooks look them up as module globals at call time). DQI_INT_PDRFIX=0 installs nothing.
The function bodies are copies of Z1's pdr copy (sha dc6e2b22, snap_rv6) with the repair lines marked NEW.
"""
import os
import sys

P = None          # the loaded dqi_z2_pdr module
LOG = dict(route=0, lower=0, cg=0)
BC = os.environ.get('PDRFIX_BC', '1') == '1'          # the (B, C) repair: route hD + lower / CG with B not F
ROUTE = os.environ.get('PDRFIX_ROUTE', '1') == '1'    # v2 selectors (cheaper) + the (F, C2) law


def _bnf(b, B, F):
    """The lower / CG flag: B not F with the (B, C) repair, else B itself (nothing allocated)."""
    if not BC:
        return B, None
    lits = [(B, False), (F, True)]
    return P.flag_and(b, lits), lits


def _bnf_erase(b, h, lits):
    if lits is not None:
        P.erase_flag_and(b, h, lits)


XDBG = os.environ.get('PDRFIX_XDBG') == '1'     # debug runs only (with Z2_PDR_MODE=debug): candidate words kept


def _xdbg(b, U, fl):
    """Copies, at the route's moment, of the words a gap law may read: payload tops s_m, q0_m' (f[-2]) and q1_m'
    (f[-1]) of levels U-1 .. U-6, and the four route flags (bit 0); names x<kind><U>_<m>[_<k>] (printed by pdrdbg)."""
    def copy(name, src):
        w = P.dbg_word(name)
        for j, s in enumerate(src):
            b.cx(s, w[j])
    for m in range(U - 1, max(U - 7, 1), -1):
        d = P.lv(m)
        for k, s in enumerate(P.tops(m)):
            copy('xS%d_%d_%d' % (U, m, k), s)
        copy('xQ%d_%d' % (U, m), d['f'][-2])
        copy('xH%d_%d' % (U, m), d['f'][-1])
    for nm, m in (('B', U - 1), ('C', U - 2), ('F', U - 3), ('C2', U - 4)):
        copy('xG%s%d' % (nm, U), [fl[m]])


def route_block_T(b, U, given=None):
    ar = P._ar()
    given = given or {}
    c = {m: P.lv(m) for m in (U, U - 1, U - 2, U - 3, U - 4)}
    L = U - 1
    q = {m: c[m]['f'][-1] for m in (L, U - 2, U - 3, U - 4)}
    fl = {m: given[m] if m in given else P.zero_flag(b, q[m]) for m in (L, U - 2, U - 3, U - 4)}
    Bq, Cq, Fq, C2q = fl[L], fl[U - 2], fl[U - 3], fl[U - 4]
    if XDBG and P.DEBUG:
        _xdbg(b, U, fl)
    X = [P.dbg_word('R%d_%d' % (U, k)) for k in range(len(P.tops(U)))] if P.DEBUG else P.tops(U)
    Z0 = P.masked_word(b, Bq, c[L]['f'][-2], neg=True)          # [not B] q0_L'
    if BC:
        lD = [(Bq, False), (Cq, False)]                         # NEW: the double gap B C
        hD = P.flag_and(b, lD)                                  # NEW
        b.cx(hD, Z0[0])                                         # NEW: Z0 = 1 on the double gap
    for x, s in zip(X, P.tops(L)):
        ar.mul_xor0(b, Z0, s, x)
    if BC:
        b.cx(hD, Z0[0])                                         # NEW
        P.erase_flag_and(b, hD, lD)                             # NEW
    P.erase_masked_word(b, Bq, c[L]['f'][-2], Z0, neg=True)
    if ROUTE:
        _route_rest(b, ar, U, c, X, Bq, Cq, Fq, C2q)
    else:
        _route_rest_orig(b, ar, U, c, X, Bq, Cq, Fq, C2q)
    for m in sorted(fl, reverse=True):
        if m not in given:
            P.clear_zero_flag(b, q[m], fl[m])
    LOG['route'] += 1


def _route_rest(b, ar, U, c, X, Bq, Cq, Fq, C2q):
    """v2 (16:50Z): the same terms with cheaper selectors, and the (F, C2) law [M harness xd13: with B = C = 0,
    F = C2 = 1 the dot needs s_(U-4)]: h2 = not B, o (o = not C, not F, computed first: 1 AND for h2 instead of 2 + 1);
    sU2 = not B, C2, (C xor F) (C xor F formed in F's qubit by two CNOTs: 3 ANDs in all instead of 5) -- equal to
    Z2's not B, C, C2, not F except at C = 0, F = 1, C2 = 1, where it adds s_(U-4)."""
    l1 = [(Bq, True), (Cq, False)]
    h1 = P.flag_and(b, l1)
    ZB = P.masked_word(b, h1, c[U - 2]['f'][-2])
    lo = [(Cq, True), (Fq, True)]
    o = P.flag_and(b, lo)
    l2 = [(Bq, True), (o, False)]
    h2 = P.flag_and(b, l2)
    b.cx(h2, ZB[0])
    lL = [(Bq, True), (C2q, True), (o, True)]
    sL = P.flag_and(b, lL)
    b.cx(Cq, Fq)                                                # F's qubit holds C xor F
    lU2 = [(Bq, True), (Fq, False), (C2q, False)]
    sU2 = P.flag_and(b, lU2)
    b.cx(Cq, Fq)
    tu, tl, t2 = P.tops(U - 2), P.tops(U - 3), P.tops(U - 4)
    for k in range(len(X)):
        ar.mul_xor0(b, ZB, tu[k], X[k])
        for j in range(12):
            b.ccx(sL, tl[k][j], X[k][j])
        for j in range(12):
            b.ccx(sU2, t2[k][j], X[k][j])
    b.cx(Cq, Fq)
    P.erase_flag_and(b, sU2, lU2)
    b.cx(Cq, Fq)
    P.erase_flag_and(b, sL, lL)
    b.cx(h2, ZB[0])
    P.erase_flag_and(b, h2, l2)
    P.erase_flag_and(b, o, lo)
    P.erase_masked_word(b, h1, c[U - 2]['f'][-2], ZB)
    P.erase_flag_and(b, h1, l1)


def _route_rest_orig(b, ar, U, c, X, Bq, Cq, Fq, C2q):
    """Z2's selectors, unchanged (Z1's pdr copy dc6e2b22)."""
    l1 = [(Bq, True), (Cq, False)]
    h1 = P.flag_and(b, l1)
    ZB = P.masked_word(b, h1, c[U - 2]['f'][-2])
    l2 = [(Bq, True), (Cq, True), (Fq, True)]
    h2 = P.flag_and(b, l2)
    b.cx(h2, ZB[0])
    lo = [(Cq, True), (Fq, True)]
    o = P.flag_and(b, lo)
    lL = [(Bq, True), (C2q, True), (o, True)]
    lU2 = [(Bq, True), (Cq, False), (C2q, False), (Fq, True)]
    sL, sU2 = P.flag_and(b, lL), P.flag_and(b, lU2)
    tu, tl, t2 = P.tops(U - 2), P.tops(U - 3), P.tops(U - 4)
    for k in range(len(X)):
        if 'zb' not in P.RDROP:
            ar.mul_xor0(b, ZB, tu[k], X[k])
        if 'sl' not in P.RDROP:
            for j in range(12):
                b.ccx(sL, tl[k][j], X[k][j])
        if 'su2' not in P.RDROP:
            for j in range(12):
                b.ccx(sU2, t2[k][j], X[k][j])
    P.erase_flag_and(b, sU2, lU2)
    P.erase_flag_and(b, sL, lL)
    P.erase_flag_and(b, o, lo)
    b.cx(h2, ZB[0])
    P.erase_flag_and(b, h2, l2)
    P.erase_masked_word(b, h1, c[U - 2]['f'][-2], ZB)
    P.erase_flag_and(b, h1, l1)


def lower_payload(b, L, k, g, v, F):
    ar = P._ar()
    s, low = v[-1], v[:-1]
    X = P.dbg_word('L%d_%d' % (L, k)) if P.DEBUG else P.lv(L + 1)['payloads'][k][-1]
    q1 = P.lv(L)['f'][-1]
    d, e1, e2 = b.alloc_qubits(12), b.alloc_qubits(12), b.alloc_qubits(12)
    P.point_products(b, g, low, [s, d, e1, e2], [lambda t, f: 1, P.C_D, P.C_E1, P.C_E2])
    for j in range(12):
        b.ccx(F, d[j], s[j])
    P.erase_acc(b, g, low, d, P.C_D)
    for j in range(12):
        b.ccx(F, e2[j], e1[j])
    ar.mul_xor0(b, q1, e1, X)                    # X ^= q1' S1sel
    B, own = P.get_B(b, L, q1)
    h, lBn = _bnf(b, B, F)                       # NEW: B not C, C = F (L's own gap flag)
    for j in range(12):
        b.ccx(h, e1[j], X[j])                    # X ^= (B not C) S1sel
    _bnf_erase(b, h, lBn)                        # NEW
    if own:
        P.clear_zero_flag(b, q1, B)
    for j in range(12):
        b.ccx(F, e2[j], e1[j])
    P.erase_acc(b, g, low, e1, P.C_E1)
    P.erase_acc(b, g, low, e2, P.C_E2)
    P.LOG['lower_ops'] += 1
    LOG['lower'] += 1


def cg_part(b, L, g, F):
    ar = P._ar()
    q1 = P.lv(L)['f'][-1]
    b.cswap(F, g[-2], g[-3])
    CG = b.alloc_qubits(12)
    ar.mul_xor0(b, q1, g[-2], CG)
    B, own = P.get_B(b, L, q1)
    h, lBn = _bnf(b, B, F)                       # NEW: B not C, C = F
    for j in range(12):
        b.ccx(h, g[-2][j], CG[j])
    for k, v in enumerate(P.lv(L)['payloads']):
        s, low = v[-1], v[:-1]
        X = P.dbg_word('L%d_%d' % (L, k)) if P.DEBUG else P.lv(L + 1)['payloads'][k][-1]
        r = P.r_poly(len(low))
        P.fan_in(b, low, s, r)
        ar.mul_xor0(b, CG, s, X)
        P.fan_in(b, low, s, r)
    for j in range(12):
        b.ccx(h, g[-2][j], CG[j])
    _bnf_erase(b, h, lBn)                        # NEW
    if own:
        P.clear_zero_flag(b, q1, B)
    ar.erase_product(b, q1, g[-2], CG)
    b.free(CG)
    b.cswap(F, g[-2], g[-3])
    P.LOG['cg'] += 1
    LOG['cg'] += 1


def install(km, ns):
    global P
    if os.environ.get('DQI_INT_PDRFIX', '1') == '0':
        return {'off': True}
    P = sys.modules.get('dqi_z2_pdr')
    assert P is not None, 'pdrfix: install after dqi_z2_pdr'
    P.route_block_T = route_block_T
    P.lower_payload = lower_payload
    P.cg_part = cg_part
    return dict(installed=True, BC=BC, ROUTE=ROUTE, XDBG=XDBG)
