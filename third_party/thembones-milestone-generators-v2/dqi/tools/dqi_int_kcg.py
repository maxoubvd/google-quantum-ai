"""dqi_int_kcg.py -- B5 (20th holder, Nov 17): R62-DQI-CONVN's KCG in its REORDERED form, form B only (convert down).

R62-DQI-CONVN measured [M 102,400 shots] that at the paired lower levels L in Z1's LOWER (16, 14, 12, 10) the kit's
`cross` (q0 ^= q1 x sel(g[-2], g[-3], F)) and Z1's `p_parts` (CG = q1' x g[-2] after cswap(F, g[-2], g[-3]); q1' = q1
clipped by B = [q1 == 0]) make the same product up to B e2: CG = q1 e2 + B e2. In the forward (convert down) nothing
between the two touches q1, g or the payload words (only q0 and f[:-2]: the locator step), so this patch does at the
cross what both did: the clip flag B's zero test (moved from Z1's clip_open, held to it: 1 wire), cswap(F, g[-2],
g[-3]), CG = q1 g[-2] (one product), q0 ^= CG (the cross), CG ^= B g[-2] (12 CCX: CG = q1' e2), the payload products
X ^= t CG between their fan-outs (moved from p_parts), CG erased by HMR with its phase (phase_product on (q1, e2) and
CZ(B, e2) on the mask bits), cswap back; Z1's clip_open then takes the held B and its p_parts does nothing.
Per site: the kit's selected word and one 42-product go [P count]. Reverse paths are untouched: the final state of the
level step is the parent's on every register (the same products, reordered across ops that commute with them).
Install LAST in DQI_INT_EXT (after dqi_zn_pdl2). DQI_INT_KCG_OFF=1: install nothing. DQI_INT_KCG_LEVELS=16,14 limits.
Form C (DQI_INT_KCG_C=1; unconvert up, the mirror of form B): in Z1's patch-C handed path (tools/dqi_zn_alias.py, read
only) the level's p_parts runs BEFORE the locator undo and the cross; there it does nothing, clip_close only unclips
and HOLDS the flag B = [q1 == 0] (1 wire through the locator undo, the q0 unscaling and the alias), and the cross does
the whole block (CG = q1 e2 once, q0 ^= CG, CG ^= B e2, the payload products, CG erased with its phase), then clears B.
Between p_parts and the cross only q0 and f[:-2] change (conv undo, q0 scaling, alias): q1, g, F and the payload
words are those p_parts read.
"""
import os
import sys

LOG = dict(scopes=0, applied=0, levels=[])


def _z1():
    return sys.modules['dqi_zn_pdl2']


FORM_A = os.environ.get('DQI_INT_KCG_A', '0') == '1'
FORM_C = os.environ.get('DQI_INT_KCG_C', '0') == '1'


def _from_handed_c():
    """True when the caller of the wrapper that calls this is Z1's patch-C handed path (dqi_zn_alias.handed)."""
    fr = sys._getframe(2)
    return fr.f_code.co_name == 'handed' and fr.f_globals.get('__name__') == 'dqi_zn_alias'


def _form_c(b, g, q1, q0, F, m):
    """Form C at the cross of the patch-C handed path (reverse, level m): runs the deferred p_parts block here."""
    Z1 = _z1()
    c = Z1.CTX.get(m)
    if not (FORM_C and c is not None and c.get('undo') and c.get('kcgC_pending')):
        return False
    assert 'kcgC_B' in c, 'kcg C: clip_close did not hold the flag at level %d' % m
    assert c['f'][-1] is q1 and c['f'][-2] is q0 and c['g'][-2] is g[-2] and c['g'][-3] is g[-3] and c['F'] is F, \
        'kcg C: the cross at level %d is not the deferred p_parts frame' % m
    ar = Z1._ar()
    SS = Z1._ss()
    B = c.pop('kcgC_B')
    c.pop('kcgC_pending')
    b.cswap(F, g[-2], g[-3])                          # e2 into g[-2]
    CG = b.alloc_qubits(12)
    ar.mul_xor0(b, q1, g[-2], CG)                     # CG = q1 e2
    b.cx(CG, q0)                                      # the cross
    for j in range(12):
        b.ccx(B, g[-2][j], CG[j])                     # CG = q1' e2
    for v, t in zip(c['payloads'], c['tops']):        # p_parts' payload products, moved here
        Z1._fan(b, v[-1], v[:-1])
        ar.mul_xor0(b, t, CG, v[-1])
        Z1._fan(b, v[-1], v[:-1])
    mask = b.alloc_bits(12)
    b.hmr(CG, mask)
    b.free(CG)
    ar.phase_product(b, q1, g[-2], mask)
    for k in range(12):
        b.ccz(B, g[-2][k], mask[k])
    b.free(mask)
    b.cswap(F, g[-2], g[-3])
    SS.clear_flag(b, q1, B)                           # clip_close's clear, moved here (q1 unclipped since then)
    LOG['applied_C'] = LOG.get('applied_C', 0) + 1
    LOG.setdefault('levels_C', []).append(m)
    return True


def _form_a(b, g, q1, q0, F, m):
    """Form A (unconvert down; DQI_INT_KCG_A=1): Z2's cg_part (tools/dqi_z2_pdr.py, read only) runs at the first payload
    scaling of a reverse direct_post at its lower levels and makes CG = q1 e2 again; here its whole body runs AT the
    cross (which then costs nothing more than q0 ^= CG) and Z2's hook is told it is done. Between the two only q0 and
    f[:-2] change (the reverse locator step). Returns True when it ran."""
    Z2 = sys.modules.get('dqi_z2_pdr')
    if Z2 is None or not FORM_A:
        return False
    cur = Z2.ST.get('cur')
    if not (Z2.REV and Z2.ST.get('mode') == 'rev' and not Z2.ST.get('bottom_done') and cur is not None
            and cur.get('reverse') and cur.get('m') == m and m in Z2.S('lower') and not cur.get('cg_done')):
        return False
    lvL, lvU = Z2.lv(m), Z2.lv(m + 1)
    if not (lvL['f'][-1] is q1 and lvL['g'][-2] is g[-2] and lvL['g'][-3] is g[-3] and cur['F'] is F):
        LOG['a_mismatch'] = LOG.get('a_mismatch', 0) + 1
        return False
    ar = Z2._ar()
    b.cswap(F, g[-2], g[-3])
    CG = b.alloc_qubits(12)
    ar.mul_xor0(b, q1, g[-2], CG)                     # CG = q1 e2
    b.cx(CG, q0)                                      # the cross
    B, own = Z2.get_B(b, m, q1)
    for j in range(12):
        b.ccx(B, g[-2][j], CG[j])                     # CG = (q1 + B) e2
    for k, v in enumerate(lvL['payloads']):
        s, low = v[-1], v[:-1]
        X = lvU['payloads'][k][-1]
        r = Z2.r_poly(len(low))
        Z2.fan_in(b, low, s, r)
        ar.mul_xor0(b, CG, s, X)
        Z2.fan_in(b, low, s, r)
    for j in range(12):
        b.ccx(B, g[-2][j], CG[j])                     # back to q1 e2
    if own:
        Z2.clear_zero_flag(b, q1, B)
    ar.erase_product(b, q1, g[-2], CG)
    b.free(CG)
    b.cswap(F, g[-2], g[-3])
    cur['cg_done'] = True
    LOG['applied_A'] = LOG.get('applied_A', 0) + 1
    return True


def make_cross(orig, levels):
    def cross(b, g, q1, q0, F):
        Z1 = _z1()
        m = len(g)
        if _form_a(b, g, q1, q0, F, m):
            return None
        if _form_c(b, g, q1, q0, F, m):
            return None
        c = Z1.CTX.get(m)
        if (m in levels and c is not None and Z1._active(m) and not c.get('undo')
                and c.get('tops') is not None and c.get('payloads') is not None and not c.get('kcg_pending')):
            ar = Z1._ar()
            SS = Z1._ss()
            B = SS.fresh_zero_flag(b, q1)                 # the clip flag (Z1's clip_open makes it later otherwise)
            b.cswap(F, g[-2], g[-3])                      # e2 into g[-2]
            CG = b.alloc_qubits(12)
            ar.mul_xor0(b, q1, g[-2], CG)                 # CG = q1 e2
            b.cx(CG, q0)                                  # the cross: q0 ^= q1 e2
            for j in range(12):
                b.ccx(B, g[-2][j], CG[j])                 # CG = (q1 + B) e2 = q1' e2
            for v, t in zip(c['payloads'], c['tops']):    # p_parts' payload products, moved here
                Z1._fan(b, v[-1], v[:-1])
                ar.mul_xor0(b, t, CG, v[-1])
                Z1._fan(b, v[-1], v[:-1])
            mask = b.alloc_bits(12)
            b.hmr(CG, mask)
            b.free(CG)
            ar.phase_product(b, q1, g[-2], mask)          # (-1)^(mask . q1 e2)
            for k in range(12):
                b.ccz(B, g[-2][k], mask[k])               # (-1)^(mask . B e2)
            b.free(mask)
            b.cswap(F, g[-2], g[-3])
            c['kcg_B'] = B
            c['kcg_pending'] = True
            LOG['applied'] += 1
            LOG['levels'].append(m)
            return None
        return orig(b, g, q1, q0, F)
    cross._b5_kcg = True
    return cross


def make_clip_open(orig):
    def clip_open(b, m):
        c = _z1().CTX.get(m)
        if c is not None and c.get('kcg_pending') and 'kcg_B' in c:
            B = c.pop('kcg_B')
            b.cx(B, c['f'][-1][0])
            return B
        return orig(b, m)
    clip_open._b5_kcg = True
    return clip_open


def make_p_parts(orig, levels=()):
    def p_parts(b, m):
        c = _z1().CTX.get(m)
        if c is not None and c.pop('kcg_pending', False):
            assert 'kcg_B' not in c, 'kcg: clip_open did not take the held flag at level %d' % m
            return None
        if FORM_C and c is not None and m in levels and c.get('undo') and not c.get('kcgC_pending') \
                and _from_handed_c():
            c['kcgC_pending'] = True                  # form C: the block runs at this level's cross
            LOG['deferred_C'] = LOG.get('deferred_C', 0) + 1
            return None
        return orig(b, m)
    p_parts._b5_kcg = True
    return p_parts


def make_clip_close(orig):
    def clip_close(b, m, B):
        c = _z1().CTX.get(m)
        if c is not None and c.get('kcgC_pending') and 'kcgC_B' not in c:
            b.cx(B, c['f'][-1][0])                    # unclip only; B = [q1 == 0] held to the cross
            c['kcgC_B'] = B
            return None
        return orig(b, m, B)
    clip_close._b5_kcg = True
    return clip_close


def install(km, ns):
    if os.environ.get('DQI_INT_KCG_OFF', '0') == '1':
        return dict(LOG, off=True)
    Z1 = _z1()
    levels = tuple(int(x) for x in os.environ.get('DQI_INT_KCG_LEVELS', ','.join(map(str, Z1.LOWER))).split(',') if x)
    walk = sys.modules['walk_bound']
    seen = set()
    for fn in walk.functions(ns.decode):
        ds = [fn.__globals__] + [cl.cell_contents for cl in (fn.__closure__ or ())
                                 if isinstance(getattr(cl, 'cell_contents', None), dict)]
        for d in ds:
            if id(d) in seen:
                continue
            seen.add(id(d))
            cr = d.get('cross')
            if not (callable(cr) and 'direct_pre' in d and 'dk' in d):
                continue
            # Z1's alias (tools/dqi_zn_alias.py, read only) wraps the kit's cross at the same levels and brackets it by
            # q0-multiples XORed into f[:-2]: the reordered form must sit INSIDE that bracket, so the kit's own cross
            # held in the alias wrapper's closure ('orig') is the one replaced; a scope without the alias is wrapped
            # directly.
            if getattr(cr, '_zn_alias', False) and cr.__closure__:
                cells = dict(zip(cr.__code__.co_freevars, cr.__closure__))
                oc = cells.get('orig')
                if oc is not None and not getattr(oc.cell_contents, '_b5_kcg', False):
                    oc.cell_contents = make_cross(oc.cell_contents, levels)
                    LOG['alias_inner'] = LOG.get('alias_inner', 0) + 1
            elif not getattr(cr, '_b5_kcg', False):
                d['cross'] = make_cross(cr, levels)
                LOG['scopes'] += 1
    if not getattr(Z1.clip_open, '_b5_kcg', False):
        Z1.clip_open = make_clip_open(Z1.clip_open)
    if not getattr(Z1.p_parts, '_b5_kcg', False):
        Z1.p_parts = make_p_parts(Z1.p_parts, levels)
    if FORM_C and not getattr(Z1.clip_close, '_b5_kcg', False):
        Z1.clip_close = make_clip_close(Z1.clip_close)
    LOG['want_levels'] = list(levels)
    LOG['forms'] = 'B' + ('A' if FORM_A else '') + ('C' if FORM_C else '')
    return LOG
