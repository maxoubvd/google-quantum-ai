#!/usr/bin/env python3
"""r49_dqi_fold_rp_merge.py -- R49-DQI-FOLD re-pricer: ledger row 7 of R49-DQI-FOLD built FORWARD ONLY in the line's
converter (EXT module; list it after dqi_zn_pdl2 and dqi_z2_pdr, before r46_dqi_conv_zpark).

Row 7: Z1's carry product by ZB and the kit q_post's Gram masked product of level U - 3 land in the same word
tops(U-2) with the same multiplier q0_(U-2) under exclusive flags; made as ONE product.
Forward, per Z1 block U in R49_M_SITES (default: every block of ZN_PDL_BLOCKS):
  route_pass (base case): h1, the ZB word, h2 and the two products tops(U-2) ^= ZB X are left out (t and the sL rows
     stay; the line's sU2 drop R46_NOSU2=1 is required);
  q_post(U-3) of Z2's bottom pass (the kit q_post ST['kitq'] = decode_backhalf_q7.q_post): its masked_products step
     becomes  h1 = not B and C and not F,  h2 = not B and not C and not F  (B, C = [q1 == 0] at levels U-1, U-2:
     Z1's HAND flag when HAND holds it, else a fresh zero test of the word as it stands, cleared after; inside the
     q_post the word p holds q1_(U-3) under F, so a test of p reads C or F, which the not-F literal absorbs),
     w = h1 xor F (one CX: exclusive),  u = w c (12 ANDs),  u_0 ^= h2,
     per payload  cswap(F, X, v[-1]); t ^= u X; cswap(F, X, v[-1])      (X = tops(U): under F the q_post's own operand)
     u_0 ^= h2; u erased by HMR, phase repaid by classically controlled CCZ(w, c_i, mask_i); w -> h1; flags erased.
Reverse (R49_M_DIRS=fwd, default): unchanged (the line's own route_pass and reverse q_posts). The forward map must
equal the line's, so the round trip stays clean and the forward words equal the reference's exactly when the merge
is exact.  R49_M_DIRS=both: the reverse route_pass leaves the same ZB products out and the line's reverse q_post
(a08's form, ST['origq'][U-3], which calls the kit's masked_products clone in the direct scope) takes the merged
product, with B, C made fresh (the reverse has no HAND).
Switches: R49_M_MODE=on|off (default on), R49_M_SITES=17,15,13,11 (levels U), R49_M_FRESH=1 (never read HAND),
R49_M_DIRS=fwd|both, R49_M_RHOLD=1 (with both: the reverse q_posts run descending and each turns its flag
F = [q1_m == 0] into a held bit at its end; for m = U-1, U-2 of a merged block that hold is delayed until reverse
q_post(U-3) has used the live flag as B or C; then the line's own hold_flag runs and its mask goes where the
q_post put it: same gates, later; +2 wires at most).
R49_M_HFWD=1 (forward): a flag the merged site must make fresh for a level whose own kit q_post runs LATER in the
same bottom pass (site 12's B = [q1_14 == 0], site 14's C = [q1_15 == 0]) is not cleared but handed to that q_post
through Z1's HAND (its zero test and the clear are both saved: -22 each); C is tested BEFORE the site's q_post starts
(inside it the word p is swapped under F). +1 wire across one or two q_posts.
"""
import atexit
import os
import sys

MODE = os.environ.get('R49_M_MODE', 'on')
FRESH = os.environ.get('R49_M_FRESH', '0') == '1'
BOTH = os.environ.get('R49_M_DIRS', 'fwd') == 'both'
RHOLD = BOTH and os.environ.get('R49_M_RHOLD', '0') == '1'
HFWD = os.environ.get('R49_M_HFWD', '0') == '1'
RSHARE = os.environ.get('R49_M_RSHARE', '0') == '1'   # reverse moved q_posts: one AND bank for both payloads
PRE = {}            # forward: level L -> (flag, word) tested before the site's q_post, handed on after it
LOG = dict(route_blocks_merged=0, route_blocks_plain=0, sites=0, sites_rev=0, flags_hand=0, flags_fresh=0,
           flags_delayed=0, holds_delayed=0, holds_released=0, mp_plain=0, kitq_calls=0, revq_calls=0,
           skipped_reverse=0)
CUR = []
DELAY = {}          # level L -> (flag F_L live, phase dict, key m+1) while a later reverse merge site needs it
SENT = object()
RSC = []            # the dict the line's reverse q_post reads masked_products / hold_flag from
RET = {}            # reverse: level m -> the flag its q_post left live in phase['returned_gap'][m + 1]


def _need():
    """level L -> set of merge sites (U - 3) whose B or C is [q1_L == 0]."""
    need = {}
    for U in _sites():
        for L in (U - 1, U - 2):
            need.setdefault(L, set()).add(U - 3)
    return need


NEED = {}


def _sites():
    import dqi_zn_pdl2 as Z1
    s = os.environ.get('R49_M_SITES', '')
    return set(int(x) for x in s.split(',') if x) if s else set(Z1.UPPER)


def make_route_pass(Z1, PC, orig):
    def route_pass(b, ups, reverse):
        if (reverse and not BOTH) or MODE != 'on':
            return orig(b, ups, reverse)
        assert os.environ.get('R46_NOSU2') == '1', 'r49 merge: written for the line (R46_NOSU2=1)'
        sites = _sites()
        if not reverse:                                     # a new converter build: reset the reverse bookkeeping
            assert not DELAY, 'r49 merge: delayed holds left from the last build %s' % sorted(DELAY)
            NEED.clear()
            NEED.update(_need())
            RET.clear()
        SS = Z1._ss()
        ar = Z1._ar()
        fl = {}
        REC = Z1.REC

        def q1(m):
            return REC[m]['f'][-1]

        def drop(m):
            if not reverse and m in Z1.HANDSET:
                Z1.HAND[m] = (fl.pop(m), q1(m))
                Z1.LOG['handed'] = Z1.LOG.get('handed', 0) + 1
            else:
                SS.clear_flag(b, q1(m), fl.pop(m))
        for U in (list(reversed(ups)) if reverse else list(ups)):
            need = {U - 1, U - 2, U - 3, U - 4}
            for m in sorted(set(fl) - need, reverse=True):
                drop(m)
            for m in sorted(need - set(fl)):
                fl[m] = SS.fresh_zero_flag(b, q1(m))
            Bq, Cq, Fq, C2q = fl[U - 1], fl[U - 2], fl[U - 3], fl[U - 4]
            X, tu, tl = Z1.tops(U), Z1.tops(U - 2), Z1.tops(U - 3)
            if U in sites:                                  # merged: the ZB product moves to q_post(U - 3)
                t, lt = PC.flag_and(b, [(Cq, True), (Fq, True)])
                sL, lL = PC.flag_and(b, [(Bq, True), (C2q, True), (t, True)])
                for k in range(len(X)):
                    for j in range(12):
                        b.ccx(sL, X[k][j], tl[k][j])
                PC.erase_flag_and(b, sL, lL)
                PC.erase_flag_and(b, t, lt)
                LOG['route_blocks_merged'] += 1
            else:                                           # the line's block (R46_NOSU2=1 form)
                q0u = REC[U - 2]['f'][-2]
                t, lt = PC.flag_and(b, [(Cq, True), (Fq, True)])
                h1, l1 = PC.flag_and(b, [(Bq, True), (Cq, False)])
                ZB = PC.masked_word(b, h1, q0u)
                h2, l2 = PC.flag_and(b, [(Bq, True), (t, False)])
                b.cx(h2, ZB[0])
                sL, lL = PC.flag_and(b, [(Bq, True), (C2q, True), (t, True)])
                for k in range(len(X)):
                    ar.mul_xor0(b, ZB, X[k], tu[k])
                    for j in range(12):
                        b.ccx(sL, X[k][j], tl[k][j])
                PC.erase_flag_and(b, sL, lL)
                b.cx(h2, ZB[0])
                PC.erase_flag_and(b, h2, l2)
                PC.erase_masked_word(b, h1, q0u, ZB)
                PC.erase_flag_and(b, h1, l1)
                PC.erase_flag_and(b, t, lt)
                LOG['route_blocks_plain'] += 1
        for m in sorted(fl, reverse=True):
            drop(m)
    route_pass._r49m = True
    return route_pass


def same_word(x, y):
    return x is y or list(x) == list(y)


def _moved():
    import dqi_z2_pdr as Z2
    return set(Z2.SETS['fwd']['moved'])


def make_kitq(orig, key='kitq_calls', scope=None):
    def kitq(b, q, payloads, H, m, parent, reverse=False, phase=None):
        LOG[key] += 1
        CUR.append(dict(m=m, q=q, parent=parent, payloads=payloads, reverse=reverse, phase=phase))
        stub = None
        pre = None
        if HFWD and MODE == 'on' and not reverse and parent is not None and (m + 3) in _sites():
            import dqi_zn_pdl2 as Z1
            L = m + 1                                # C's level; its own q_post is the next one of the pass
            if L not in Z1.HAND and L in _moved():
                pre = (L, Z1._ss().fresh_zero_flag(b, parent[1]), parent[1])   # p before this q_post swaps it
                PRE[L] = (pre[1], pre[2])
                LOG['flags_pre'] = LOG.get('flags_pre', 0) + 1
        if RHOLD and reverse and scope is not None and NEED.get(m) and phase is not None:
            Fm = phase['given_gap'].get(m + 1)       # the flag this reverse q_post pops first
            real = scope['hold_flag']
            box = {}

            def hold_flag(b_, F_, *a, **kw):
                if F_ is Fm and not box:
                    box['F'] = F_
                    return SENT
                return real(b_, F_, *a, **kw)
            stub = (real, box)
            scope['hold_flag'] = hold_flag
        try:
            return orig(b, q, payloads, H, m, parent, reverse, phase)
        finally:
            CUR.pop()
            if pre is not None:                      # hand C on to the next q_post (its own test is saved)
                import dqi_zn_pdl2 as Z1
                L, Fl, word = pre
                PRE.pop(L, None)
                Z1.HAND[L] = (Fl, word)
                LOG['flags_handed_on'] = LOG.get('flags_handed_on', 0) + 1
            if stub is not None:
                real, box = stub
                scope['hold_flag'] = real
                if box and phase['return_mask'].get(m + 1) is SENT:
                    DELAY[m] = (box['F'], phase, m + 1)
                    LOG['holds_delayed'] += 1
            if reverse and phase is not None and phase.get('returned_gap', {}).get(m + 1) is not None:
                RET[m] = phase['returned_gap'][m + 1]          # kept live by the kit (m + 1 in KEEP_CHILD)
    kitq._r49m = True
    kitq._r49m_orig = orig
    return kitq


def release(b, scope, site):
    for L in list(DELAY):
        NEED.get(L, set()).discard(site)
        if not NEED.get(L):
            F, phase, key = DELAY.pop(L)
            phase['return_mask'][key] = scope['hold_flag'](b, F)
            LOG['holds_released'] += 1


def shared_step(b, F, payloads, c, targets, ar):
    """The kit's masked_products (decode_backhalf_q7): ONE AND bank u = F c for both payloads (the reverse q_post of
    the line, conversion.py's form, builds one bank per payload)."""
    u = b.alloc_qubits(12)
    for i in range(12):
        b.init_and(F, c[i], u[i])
    for v, t in zip(payloads, targets):
        ar.mul_xor0(b, u, v[-1], t)
    mask = b.alloc_bits(12)
    b.hmr(u, mask)
    b.free(u)
    for i in range(12):
        b.ccz(F, c[i], mask[i])
    b.free(mask)


def make_mp1(Z1, PC, orig, ar):
    """masked_product(b, F, x, c, t) of the line's reverse q_post (conversion.py form: called once per payload).
    At a reverse merge site, or with R49_M_RSHARE=1 at any wrapped reverse q_post, the first call does the whole
    step for both payloads (merged product, or one shared AND bank) and the second call does nothing."""
    def masked_product(b, F, x, c, t):
        cur = CUR[-1] if CUR else None
        if MODE != 'on' or cur is None or not cur['reverse'] or cur['parent'] is None:
            return orig(b, F, x, c, t)
        merge = BOTH and (cur['m'] + 3) in _sites()
        if not (merge or RSHARE):
            return orig(b, F, x, c, t)
        k = cur.get('mp1', 0)
        cur['mp1'] = k + 1
        payloads, tops = cur['payloads'], cur['parent'][0]
        assert same_word(x, payloads[k][-1]) and same_word(t, tops[k]), 'r49 merge: masked_product call order'
        if k > 0:
            return
        if merge:
            merged_step(b, F, payloads, c, tops, cur, Z1, PC, ar)
        else:
            shared_step(b, F, payloads, c, tops, ar)
            LOG['rev_shared'] = LOG.get('rev_shared', 0) + 1
    masked_product._r49m = True
    return masked_product


def make_mp(Z1, PC, orig, ar):
    def merged_products(b, F, payloads, c, targets):
        cur = CUR[-1] if CUR else None
        if MODE != 'on' or cur is None or cur['parent'] is None or (cur['m'] + 3) not in _sites():
            LOG['mp_plain'] += 1
            return orig(b, F, payloads, c, targets)
        if cur['reverse'] and not BOTH:
            LOG['skipped_reverse'] += 1
            return orig(b, F, payloads, c, targets)
        merged_step(b, F, payloads, c, targets, cur, Z1, PC, ar)
    merged_products._r49m = True
    return merged_products


def merged_step(b, F, payloads, c, targets, cur, Z1, PC, ar):
    if True:
        rev = cur['reverse']
        m = cur['m']
        U = m + 3
        tops, p, cc = cur['parent']
        assert cc is c and all(same_word(x, y) for x, y in zip(tops, targets))
        assert all(same_word(x, y) for x, y in zip(targets, Z1.tops(U - 2))), 'r49 merge: targets are not tops(U-2)'
        assert same_word(c, Z1.REC[U - 2]['f'][-2]), 'r49 merge: c is not q0_(U-2)'
        X = Z1.tops(U)
        SS = Z1._ss()
        fresh = []
        phase = cur.get('phase') or {}

        def flag(level, word, hand_on=False):
            if rev and level in DELAY:
                LOG['flags_delayed'] += 1
                return DELAY[level][0]
            if rev and level in RET and phase.get('returned_gap', {}).get(level + 1) is RET[level]:
                LOG['flags_returned'] = LOG.get('flags_returned', 0) + 1     # the kit keeps it live (KEEP_CHILD)
                return RET[level]
            if not rev and level in PRE:
                return PRE[level][0]                 # tested before this q_post; handed on by the wrapper
            if hand_on:                              # fresh, then handed to the later q_post of that level
                Fl = SS.fresh_zero_flag(b, word)
                handon.append((level, Fl, word))
                LOG['flags_fresh'] += 1
                return Fl
            h = None if (FRESH or rev) else Z1.HAND.get(level)
            if h is not None and same_word(h[1], word):
                LOG['flags_hand'] += 1
                return h[0]
            Fl = SS.fresh_zero_flag(b, word)
            fresh.append((word, Fl))
            LOG['flags_fresh'] += 1
            return Fl
        wB = Z1.REC[U - 1]['f'][-1]
        assert same_word(p, Z1.REC[U - 2]['f'][-1]), 'r49 merge: p is not q1_(U-2)'
        handon = []
        B = flag(U - 1, wB, hand_on=(HFWD and not rev and (U - 1) in _moved() and (U - 1) not in Z1.HAND))
        C = flag(U - 2, p)                      # p holds q1_(U-3) under F here: a fresh test reads C or F
        h1, l1 = PC.flag_and(b, [(B, True), (C, False), (F, True)])
        h2, l2 = PC.flag_and(b, [(B, True), (C, True), (F, True)])
        b.cx(F, h1)                              # w = h1 xor F
        u = b.alloc_qubits(12)
        for i in range(12):
            b.init_and(h1, c[i], u[i])
        b.cx(h2, u[0])
        for v, t, x in zip(payloads, targets, X):
            b.cswap(F, x, v[-1])
            ar.mul_xor0(b, u, x, t)
            b.cswap(F, x, v[-1])
        b.cx(h2, u[0])
        mask = b.alloc_bits(12)
        b.hmr(u, mask)
        b.free(u)
        for i in range(12):
            b.ccz(h1, c[i], mask[i])
        b.free(mask)
        b.cx(F, h1)
        PC.erase_flag_and(b, h2, l2)
        PC.erase_flag_and(b, h1, l1)
        for word, Fl in reversed(fresh):
            SS.clear_flag(b, word, Fl)
        for level, Fl, word in handon:               # B of level U-1: its own q_post comes two calls later
            Z1.HAND[level] = (Fl, word)
            LOG['flags_handed_on'] = LOG.get('flags_handed_on', 0) + 1
        LOG['sites_rev' if rev else 'sites'] += 1
        if rev and RHOLD:
            release(b, RSC[0], m)


def install(km, ns):
    import dqi_zn_pdl2 as Z1
    import dqi_z2_pdr as Z2
    import dqi_zn_pdc as PC
    if not getattr(Z1.route_pass, '_r49m', False):
        Z1.route_pass = make_route_pass(Z1, PC, Z1.route_pass)
    kq = Z2.ST['kitq']
    assert kq is not None, 'r49 merge: install after dqi_z2_pdr'
    base = getattr(kq, '_r49m_orig', None) or kq
    gl = base.__globals__
    if not getattr(gl['masked_products'], '_r49m', False):
        gl['masked_products'] = make_mp(Z1, PC, gl['masked_products'], gl['ar'])
    LOG['mp_dicts'] = 1
    if not getattr(kq, '_r49m', False):
        Z2.ST['kitq'] = make_kitq(kq)
    if BOTH or RSHARE:
        NEED.clear()
        NEED.update(_need())
        levels = set()
        if BOTH:
            levels |= set(U - 3 for U in _sites())
            if RHOLD:
                levels |= set(L for L in NEED if L in Z2.ST['origq'])
        if RSHARE:
            levels |= set(m for m in Z2.SETS['rev']['moved'] if m in Z2.ST['origq'])
        files = set()
        for m in sorted(levels):
            oq = Z2.ST['origq'].get(m)
            if oq is None or getattr(oq, '_r49m', False):
                continue
            # the line's reverse q_post at the moved levels is conversion.py's (it reads masked_product, one bank per
            # payload, and hold_flag from its globals); a08's form reads masked_products from its closure `scope`
            cells = dict(zip(oq.__code__.co_freevars, oq.__closure__ or ()))
            sc = cells['scope'].cell_contents if 'scope' in cells else oq.__globals__
            ar = sc.get('ar', gl['ar'])
            if 'masked_product' in sc and not getattr(sc['masked_product'], '_r49m', False):
                sc['masked_product'] = make_mp1(Z1, PC, sc['masked_product'], ar)
                LOG['mp_dicts'] += 1
            if 'masked_products' in sc and not getattr(sc['masked_products'], '_r49m', False):
                sc['masked_products'] = make_mp(Z1, PC, sc['masked_products'], ar)
                LOG['mp_dicts'] += 1
            RSC[:] = [sc]
            files.add('%s:%d' % (oq.__code__.co_filename.split('/')[-1], oq.__code__.co_firstlineno))
            Z2.ST['origq'][m] = make_kitq(oq, 'revq_calls', sc if RHOLD else None)
        LOG['rev_wrapped'] = sorted(levels)
        LOG['rev_scope_is_fwd'] = bool(RSC) and RSC[0] is gl
        LOG['revq_file'] = sorted(files)
    LOG['sites_set'] = sorted(_sites())
    LOG['mode'] = MODE
    LOG['kitq_file'] = base.__code__.co_filename.split('/')[-1]
    return dict(LOG)


@atexit.register
def _report():
    print('r49_dqi_fold_rp_merge', dict(LOG), flush=True)
