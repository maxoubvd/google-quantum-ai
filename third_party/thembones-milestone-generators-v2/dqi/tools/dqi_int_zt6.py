"""dqi_pub_zt.py -- builder D1: R32-DQI-EXT's construction Z-T in the line's zero node (research/panel/R32-DQI-EXT.md).

Z-T: the node z = 0 runs as a telescoped group on a 12-bit unit of the carried T1 scale G = (low, high) instead of
starting from the constant 1 and measuring its 12-bit gauge. split: f = [low = 0], CSWAP(f, low, high) makes the low
half a unit nu; the zero node's chain starts on nu and ends on nu' = nu x (the nonzero one of N_f, D_f), which stays
in G's low half; the return's reversed chain maps nu' back to nu. No gauge is measured, so the return runs no 12-bit
scalar chain (1,344 [M]) and no gauge payment. Designer's kernels (tools/r32_dqi-ext_ztkern.py, read-only) measured
-1,099 for the extraction pair at Q 902 [M kernels, +-40] and window 26 beyond chart + G.

What it replaces: the line's zero-node code (`extraction_companion` of exp/frontier/a12/decode_trailing_companion.py,
whose forward is the kit's `extraction` through replay_erase), both trips, by the recipe of R32-DQI-EXT section 3:
  forward: external_extraction; split; zj_in (D_1 = nu X, L_1 = kap_0 nu^2 (L0 + b0), N_1 = kap_0 nu, the two old
           words X-measured with masks kept for the return); the kit's 33 `step`s; pole toggle, pole swap,
           `normal_form`, `controlled_inverse` as before; NO gauge measurement; unsplit (f erased by its zero test);
           free SWAPs put the head back on the chart words (u in words[-2], v in words[-1]) and nu' in G's low wires.
  return:  `controlled_inverse` undo as before; the SWAPs undone; split; zret_a (1 / nu' rebuilt, normal-form mask
           paid, u, v scaled back, 1 / nu' erased by the kit's inverse phase, pole swap and pole flag undone); the
           kit's reverse `step`s; zret_b (1 / nu, the two words divided back, zj_in's masks paid, 1 / nu erased);
           unsplit; external_extraction (return).
The pole flag stays in the token (`pole`), so B5's polepark2 (parks it after the synthesis tail, rebuilds it before
the reverse tail, and uses it at unextraction's start) works unchanged. G is idle between the two trips (only carried).
Install (tools/dqi_pub_entry.py): DQI_PUB_P=...,zt (any order after B5's decode is built).
"""
import sys

import os

STATS = dict(forward=0, reverse=0)
# DQI_PUB_ZJOIN=1 (with the pfd module's matching switch): heads 1..4 of the zero node's chain without L products.
ZJOIN = os.environ.get('DQI_PUB_ZJOIN', '0') == '1'
# B5 copy (tools/dqi_int_zt6.py, Nov 15): under DQI_INT_J6=1 the prefix block runs at J = 6, so heads 1..5 of the
# zero node's chain go without their L products (step_n); otherwise this file behaves exactly as dqi_pub_zt.py.
_ZJ_LAST = 5 if os.environ.get('DQI_INT_J6', '') == '1' else 4


def _make(G, banks):
    """G: globals of a live extraction_companion clone (provides ref, external_extraction); banks: the T1 bank dict
    of the live node chain (its `_t1_banks` global: fix round 1, 13:02Z, sys.modules' copy of the patch is not the
    instance the build line's loader installed)."""
    ref = G['ref']
    ar, scale, inv = ref.ar, ref.scale, ref.inv
    ks, corr = ref.ks, ref.corrections
    step, normal_form, scale_word = ref.step, ref.normal_form, ref.scale_word
    square12, inverse_initialized, phase_inverse34 = ref.square12, ref.inverse_initialized, ref.phase_inverse34
    zero_toggle, zero_flag, controlled_inverse = ref.guarded.zero_toggle, ref.zero_flag, ref.controlled_inverse
    # Both trips use the companion's own external_extraction: the patched clone (prefix block, T1 node chains) that
    # the line's return calls. Fix round 4 (13:20Z): ref.extraction's global of that name is an unpatched copy (no
    # prefix block, no T1 scale was ever made), which is why rounds 1-3 found no scale.
    external_extraction = G['external_extraction']
    external_forward = external_extraction
    assert '_dqb_block' in external_extraction.__code__.co_names, 'companion external_extraction is not the patched one'
    phase_zero = sys.modules['exp.arithmetic.division'].phase_zero
    assert corr[0] == 0

    held = {}

    def bank_of(b):
        """The T1 bank entry of this decode: the one entry, among the bank dicts of every live node chain, that
        holds a scale. (Fix rounds 2-3, 13:07Z / 13:15Z: the node chains run on a wrapped builder, and the clones
        that run are not all reachable from the decode's globals, so all live functions are searched.)"""
        import gc
        dicts = {id(banks): banks}
        for o in gc.get_objects():
            if isinstance(o, type(bank_of)) and o.__name__ == 'external_node_chain':
                d = o.__globals__.get('_t1_banks')
                if isinstance(d, dict):
                    dicts[id(d)] = d
        live = [v for d in dicts.values() for v in d.values()
                if isinstance(v, dict) and v.get('scale') is not None]
        print(f'dqi_pub_zt: {len(dicts)} T1 bank dicts, {len(live)} entries hold a scale', file=sys.stderr, flush=True)
        assert len(live) == 1, f'{len(live)} T1 bank entries hold a scale'
        return live[0]

    def step_n(b, head, N, D, L, kap, reverse=False):
        """The kit's zero-node `step` without its L product (heads 1..4 under DQI_PUB_ZJOIN: the prefix block has
        added their L contribution to L_0); L is still scaled by kap at every level."""
        a = head[0]
        if reverse:
            scale(b, N, inv(kap))
            b.swap(N, D)
            ar.mul_xor0(b, a, D, N)
            scale(b, L, inv(kap))
        else:
            scale(b, L, kap)
            ar.mul_xor0(b, a, D, N)
            b.swap(N, D)
            scale(b, N, kap)

    def split(b, low, high):
        f = zero_flag(b, low)
        b.cswap(f, low, high)
        return f

    def unsplit(b, f, low, high):
        b.cswap(f, low, high)
        bit = b.alloc_bits(1)
        b.hmr(f, bit[0])
        b.free(f)
        phase_zero(b, low, bit[0])
        b.free(bit)

    def zt_extraction(b, words, reverse=False, token=None):
        heads = [words[2 * i:2 * i + 2] for i in range(34)]
        N, L = words[-2:]
        a0, h0 = heads[0]
        if not reverse:
            STATS['forward'] += 1
            ext = external_forward(b, words)
            entry = bank_of(b)
            held[id(b)] = entry
            Gs = entry['scale']                              # T1's carried scale, left by group 16
            C, MU = Gs[:12], Gs[12:]
            f = split(b, C, MU)
            # zj_in: first step with the unit
            b.cx(a0, N)
            b.cx(h0, L)
            T = b.alloc_qubits(12)
            ar.mul_xor0(b, C, N, T)                          # T = nu (u0 + a0)
            mN = b.alloc_bits(12)
            b.hmr(N, mN)
            square12(b, C)
            ar.mul_xor0(b, C, L, N)                          # N = nu^2 (v0 + b0)
            square12(b, C, True)
            mL = b.alloc_bits(12)
            b.hmr(L, mL)
            scale(b, C, ks[0])                               # C = kap_0 nu = N_1
            scale(b, N, ks[0])                               # N-wires = L_1
            b.swap(T, L)                                     # L-wires = D_1; T clean again
            b.free(T)
            Nc, Dc, Lc = C, L, N
            for i in range(1, 34):
                if ZJOIN and i <= _ZJ_LAST:
                    step_n(b, heads[i], Nc, Dc, Lc, ks[i])
                else:
                    step(b, heads[i], Nc, Dc, Lc, ks[i], corr[i])
            pole = b.alloc_qubits(1)[0]
            zero_toggle(b, Dc, pole)
            b.cswap(pole, Nc, Dc)
            norm = normal_form(b, Nc, Dc, Lc)
            ci = controlled_inverse(b, Lc, pole)
            unsplit(b, f, Dc, MU)                            # Dc (L-wires) holds nu' (or 0 with nu' in MU)
            # head back on the chart words, nu' back on G's low wires: C=u, N=v, L=nu'  ->  N=u, L=v, C=nu'
            b.swap(C, N)
            b.swap(C, L)
            a3 = _a3()
            if a3 is not None and REC['pole'] is not None and a3._qidx(pole) != REC['pole']:
                pole = _move_pole(b, pole, REC['pole'])
            _restore_at(b, 'F')
            return dict(external=ext, pole=pole, normalization=norm, chart_inverse=ci, zt_masks=(mN, mL))
        STATS['reverse'] += 1
        pole = token['pole']
        controlled_inverse(b, L, pole, True, token['chart_inverse'])
        entry = held.pop(id(b), None) or bank_of(b)
        Gs = entry['scale']
        C, MU = Gs[:12], Gs[12:]
        b.swap(C, L)                                         # C = v, L = nu'
        b.swap(C, N)                                         # C = u, N = v
        Nc, Dc, Lc = C, L, N
        f = split(b, Dc, MU)
        # zret_a
        di = inverse_initialized(b, Dc)
        b.cz(di, token['normalization'])
        b.free(token['normalization'])
        scale_word(b, Nc, Dc, di)
        square12(b, Dc)
        square12(b, di)
        scale_word(b, Lc, Dc, di)
        square12(b, di, True)
        square12(b, Dc, True)
        m2 = b.alloc_bits(12)
        b.hmr(di, m2)
        b.free(di)
        phase_inverse34(b, Dc, m2)
        b.free(m2)
        b.cswap(pole, Nc, Dc)
        zero_toggle(b, Dc, pole)
        b.free(pole)
        for i in reversed(range(1, 34)):
            if ZJOIN and i <= _ZJ_LAST:
                step_n(b, heads[i], Nc, Dc, Lc, ks[i], True)
            else:
                step(b, heads[i], Nc, Dc, Lc, ks[i], corr[i], True)
        # zret_b
        scale(b, Nc, inv(ks[0]))
        scale(b, Lc, inv(ks[0]))
        di = inverse_initialized(b, Nc)
        scale_word(b, Dc, di, Nc)
        square12(b, Nc)
        square12(b, di)
        scale_word(b, Lc, di, Nc)
        square12(b, di, True)
        square12(b, Nc, True)
        mN, mL = token['zt_masks']
        b.cz(Dc, mN)
        b.free(mN)
        b.cz(Lc, mL)
        b.free(mL)
        m3 = b.alloc_bits(12)
        b.hmr(di, m3)
        b.free(di)
        phase_inverse34(b, Nc, m3)
        b.free(m3)
        b.cx(a0, Dc)                                         # L-wires = u0
        b.cx(h0, Lc)                                         # N-wires = v0
        b.swap(N, L)                                         # u0 back on N's wires, v0 on L's
        unsplit(b, f, Nc, MU)                                # G restored on its own wires
        _restore_at(b, 'R')
        external_extraction(b, words, True, token['external'])

    return zt_extraction


_TRAMPOLINE = '''
def extraction_companion(b, words, reverse=False, token=None):
    return _d1_zt_extraction(b, words, reverse, token)
'''

# ---- free-stack parity with the line before Z-T (fix round 5, 13:28Z) -------------------------------------------
# kickmix reuses freed qubits last-in first-out, so Z-T's different allocations in the zero node move every later
# placement; a3's own sync heals the forward at its converter points, but Z-T's return runs in unextraction, after
# them, and the final routing then failed on q21 (p2zte). A scratch decode of the line WITHOUT Z-T records the free
# stack (and the zero node's pole wire) at F = exit of the forward zero node and R = just before the return's
# external_extraction; the Z-T build reproduces both there (pole moved to the recorded wire by a free SWAP).
REC = dict(F=None, R=None, pole=None, fixes=0, misses=0)


def _a3():
    return sys.modules.get('dqi_pub_a3')


def _move_pole(b, pole, want):
    """Put the pole flag on wire `want` (free in this build): take every spare wire, keep `want`, SWAP, free."""
    a3 = _a3()
    got = a3._spare(b)
    have = {a3._qidx(q): q for q in got}
    if want not in have:
        for q in reversed(got):
            b.free(q)
        return pole
    target = have[want]
    for q in reversed(got):
        if a3._qidx(q) != want:
            b.free(q)
    b.swap(pole, target)
    b.free(pole)
    return target


def _restore_at(b, key):
    a3 = _a3()
    order = REC[key]
    if order is None or a3 is None:
        REC['misses'] += 1
        return
    got = a3._spare(b)
    have = {a3._qidx(q): q for q in got}
    if set(have) != set(order):
        REC['misses'] += 1
        for q in reversed(got):
            b.free(q)
        print(f'dqi_pub_zt: stack {key} MISS ({len(set(have) ^ set(order))} wires differ)', file=sys.stderr, flush=True)
        return
    for i in reversed(order):
        b.free(have[i])
    REC['fixes'] += 1


def _capture(km, f, G):
    """Scratch decode of the line as it is (Z-T not installed), recording F, R and the pole wire."""
    a3 = _a3()
    import types as _t
    orig = _t.FunctionType(f.__code__, f.__globals__, f.__name__, f.__defaults__, f.__closure__)
    real_ee = G['external_extraction']

    def ee_rec(b, words, reverse=False, token=None):
        if reverse:
            REC['R'] = a3._snapshot(b)
        return real_ee(b, words, reverse, token)

    def comp_rec(b, words, reverse=False, token=None):
        out = orig(b, words, reverse, token)
        if not reverse:
            REC['F'] = a3._snapshot(b)
            REC['pole'] = a3._qidx(out['pole'])
        return out

    scope = {}
    exec(compile(_TRAMPOLINE.replace('_d1_zt_extraction', '_d1_zt_capture'), __file__ + ':cap', 'exec'), scope)
    saved_code = f.__code__
    G['_d1_zt_capture'] = comp_rec
    G['external_extraction'] = ee_rec
    f.__code__ = scope['extraction_companion'].__code__
    st = dict(a3.STACK)                                  # a3's replay counters must not be spent on this scratch
    CB = km.CircuitBuilder
    tagged = CB.mark
    fv = getattr(getattr(tagged, '__code__', None), 'co_freevars', ())
    if 'done' in fv and 'original' in fv:                # r6tag's once-per-process tag: not on the scratch build
        CB.mark = dict(zip(fv, tagged.__closure__))['original'].cell_contents
    try:
        b0 = CB()
        S0 = [b0.alloc_qubits(12) for _ in range(70)]
        E0 = [b0.alloc_bits(12) for _ in range(4095)]
        _CAPTURE_DECODE[0](b0, S0, E0)
    finally:
        CB.mark = tagged
        f.__code__ = saved_code
        G['external_extraction'] = real_ee
        a3.STACK.update(st)
    del b0, S0, E0
    import gc
    gc.collect()


_CAPTURE_DECODE = [None]


def install(km, ns):
    walk = sys.modules['walk_bound']
    scope = {}
    exec(compile(_TRAMPOLINE, __file__ + ':zt', 'exec'), scope)
    code = scope['extraction_companion'].__code__
    fns = walk.functions(ns.decode)
    bank_dicts = {id(f.__globals__['_t1_banks']): f.__globals__['_t1_banks'] for f in fns
                  if f.__name__ == 'external_node_chain' and '_t1_banks' in f.__globals__}
    assert len(bank_dicts) == 1, f'expected one T1 bank dict, found {len(bank_dicts)}'
    banks = next(iter(bank_dicts.values()))
    comps = [f for f in fns if f.__name__ == 'extraction_companion' and 'forward_original' in f.__globals__
             and 'external_extraction' in f.__globals__ and 'ref' in f.__globals__]
    assert comps, 'no extraction_companion clone found'
    if _a3() is not None:                               # record the line's free stack at F and R (scratch decode)
        _CAPTURE_DECODE[0] = ns.decode
        _capture(km or sys.modules['circ'], comps[0], comps[0].__globals__)
    made = {}
    for f in comps:
        k = id(f.__globals__)
        if k not in made:
            made[k] = _make(f.__globals__, banks)
        f.__globals__['_d1_zt_extraction'] = made[k]
        f.__code__ = code
    return dict(companions=len(comps), scopes=len(made), banks=len(bank_dicts),
                captured=(REC['F'] is not None, REC['R'] is not None, REC['pole']))
