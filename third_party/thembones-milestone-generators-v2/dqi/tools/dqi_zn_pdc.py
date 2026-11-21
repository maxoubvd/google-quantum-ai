#!/usr/bin/env python3
"""dqi_zn_pdc.py -- builder Z1 (Nov 14): the kit's bez step with PAIRED payload digits, built in kickmix
(forward = convert's half-Euclid, reverse = unconvert's), on top of the kit's exp/safe_scalar converter.
Construction pd3 of tools/dqi_zn_pd.py (class PairedR2), word-exact there against the kit on every gap pattern:
  blocks (U, L = U - 1), U in UPPER (default 17, 15, 13, 11);
  descent: at U only the kit's digit scaling (X = t / a_U); at L the locator step first, then B = [q1' == 0],
           q1' clipped, Z0 = [not B] q0', CG = q1' g[-2]; per payload P = X CG with two free fan-outs, the kit's
           scaling (D0), D1 = X q1' (temp), the kit's N-point product, D1 erased, Y = s ^= Z0 X;
  bottom:  carries (bottom-up), every q_post of levels 1..17 (moved down, the kit's own, with its tokens),
           carries transposed (top-down);
  return:  at L the transposed kernel (Y^T, transposed product into s and a temp, t ^= temp q1', temp erased,
           scaling, t ^= CG (s + <r, low>)); at U only the kit's last scaling; elsewhere the kit's dots + last.
Every temporary word is a plain product or a masked word, erased by measurement with free phase products.
Usage: python dqi_zn_pdc.py <kinds=0,1,2,3,4> <seed> [upper=17,15,13,11]
"""
import os
import sys

import numpy as np

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
import dqi_zn_pdk as PK  # noqa: E402  (sets up kickmix and the kit's paths; kernels pieces)
km = PK.km
from exp.safe_scalar import converter as SS  # noqa: E402
from repair import rotate  # noqa: E402

ar = SS.ar
K = PK.K
DK = PK.DK
CFG = dict(upper=(17, 15, 13, 11))


def upper():
    return set(CFG['upper'])


def lower():
    return {u - 1 for u in CFG['upper']}


# ---------------------------------------------------------------------------------------------- small helpers
def masked_word(b, F, x, neg=False):
    """New word = [F] x (or [not F] x): 12 ANDs."""
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


def flag_and(b, terms):
    """One new qubit = AND of literals (qubit, negated); returns (qubit, terms) for erase_flag_and."""
    out = b.alloc_qubits(1)[0]
    lits = list(terms)
    for q, neg in lits:
        if neg:
            b.x(q)
    if len(lits) == 1:
        b.cx(lits[0][0], out)
        anc = []
    else:
        anc = []
        left = lits[0][0]
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
    return out, lits


def erase_flag_and(b, out, lits):
    """Erase an AND of literals by measurement (phase repaired with free classically controlled gates)."""
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
        else:                                   # rebuild the AND of all but the last, then a classically ctl. CCZ
            anc = []
            left = lits[0][0]
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


# ---------------------------------------------------------------------------------------------- level steps
def loc_ops(f, g):
    q0, q1 = f[-2], f[-1]
    return [('s', q1, True), ('m', q1, g[-2], q0), ('s', q0, True), ('c', g[:-1], [q0, q1], f[:-2])]


def loc_restore_ops(f, g):
    q0, q1 = f[-2], f[-1]
    return [('c', g[:-1], [q0, q1], f[:-2]), ('s', q0, False), ('m', q1, g[-2], q0), ('s', q1, False)]


def lower_down(b, f, g, payloads, tops, a, ai, undo):
    """Paired step at the lower level (inside the pivot's inverse context). undo = the exact inverse."""
    q0, q1 = f[-2], f[-1]
    n = len(g) - 1
    r = PK.r_poly(n)
    if not undo:
        SS.runops(b, loc_ops(f, g), a, ai)
    B = SS.fresh_zero_flag(b, q1)
    b.cx(B, q1[0])

    def yop(t, s):                                    # s ^= [not B] q0' X with a transient masked word
        Z0 = masked_word(b, B, q0, neg=True)
        ar.mul_xor0(b, Z0, t, s)
        erase_masked_word(b, B, q0, Z0, neg=True)

    def p_part():
        CG = b.alloc_qubits(12)
        ar.mul_xor0(b, q1, g[-2], CG)
        for v, t in zip(payloads, tops):
            PK.fan_out(b, v[-1], v[:-1], r)
            ar.mul_xor0(b, t, CG, v[-1])
            PK.fan_out(b, v[-1], v[:-1], r)
        ar.erase_product(b, q1, g[-2], CG)
        b.free(CG)

    def per_payload(v, t):
        s, low = v[-1], v[:-1]
        if not undo:
            SS.scale(b, s, a, ai, True)
        else:
            yop(t, s)
        D1 = b.alloc_qubits(12)
        ar.mul_xor0(b, t, q1, D1)
        DK.conv(b, g[:-1], [s, D1], low, n, undo)
        ar.erase_product(b, t, q1, D1)
        b.free(D1)
        if not undo:
            yop(t, s)
        else:
            SS.scale(b, s, a, ai, False)

    if not undo:
        p_part()
        for v, t in zip(payloads, tops):
            per_payload(v, t)
    else:
        for v, t in zip(payloads, tops):
            per_payload(v, t)
        p_part()
    b.cx(B, q1[0])
    SS.clear_flag(b, q1, B)
    if undo:
        SS.runops(b, loc_ops(f, g), a, ai, True)


def lower_up(b, f, g, payloads, tops, a, ai, undo):
    """Transposed paired step at the lower level (return; inside last's inverse context). undo = inverse."""
    q0, q1 = f[-2], f[-1]
    n = len(g) - 1
    r = PK.r_poly(n)
    if undo:
        SS.runops(b, loc_restore_ops(f, g), a, ai, True)
    B = SS.fresh_zero_flag(b, q1)
    b.cx(B, q1[0])

    def yop(s, t):                                    # t ^= [not B] q0' s with a transient masked word
        Z0 = masked_word(b, B, q0, neg=True)
        ar.mul_xor0(b, Z0, s, t)
        erase_masked_word(b, B, q0, Z0, neg=True)

    def cg_part():
        CG = b.alloc_qubits(12)
        ar.mul_xor0(b, q1, g[-2], CG)
        for v, t in zip(payloads, tops):
            PK.fan_in(b, v[:-1], v[-1], r)
            ar.mul_xor0(b, CG, v[-1], t)
            PK.fan_in(b, v[:-1], v[-1], r)
        ar.erase_product(b, q1, g[-2], CG)
        b.free(CG)

    def per_payload(v, t):
        s, low = v[-1], v[:-1]
        if not undo:
            yop(s, t)
        else:
            SS.scale(b, s, a, ai, False)
        D1 = b.alloc_qubits(12)
        PK.conv_transpose(b, g[:-1], low, [s, D1], n)
        ar.mul_xor0(b, D1, q1, t)
        PK.erase_conv_transpose(b, g[:-1], low, D1, n)
        if not undo:
            SS.scale(b, s, a, ai, True)
        else:
            yop(s, t)

    if not undo:
        for v, t in zip(payloads, tops):
            per_payload(v, t)
        cg_part()
    else:
        cg_part()
        for v, t in zip(payloads, tops):
            per_payload(v, t)
    b.cx(B, q1[0])
    SS.clear_flag(b, q1, B)
    if not undo:
        SS.runops(b, loc_restore_ops(f, g), a, ai)


def pre_pd(b, f, g, payloads, parent, m, undo=False, incoming=None, defer=False):
    if m in upper():
        def action(ai):
            SS.runops(b, [('s', v[-1], True) for v in payloads] + loc_ops(f, g), g[-1], ai, undo)
    elif m in lower():
        def action(ai):
            lower_down(b, f, g, payloads, parent[0], g[-1], ai, undo)
    else:
        def action(ai):
            n = m - 1
            ops = [('s', v[-1], True) for v in payloads]
            ops += [('m', g[i], v[-1], v[i]) for v in payloads for i in range(n)]
            SS.runops(b, ops + loc_ops(f, g), g[-1], ai, undo)
    return SS.inverse_apply(b, g[-1], action, incoming, defer)


def post_pd(b, f, g, payloads, parent, m, undo=False, incoming=None, defer=False):
    """The return at level m without its q_post (done at the bottom): kit dots + last, or the paired forms."""
    n = m - 1
    if m in upper():
        def action(ai):
            SS.runops(b, [('s', v[-1], True) for v in payloads] + loc_restore_ops(f, g), g[-1], ai, undo)
        return SS.inverse_apply(b, g[-1], action, incoming, defer)
    if m in lower():
        def action(ai):
            lower_up(b, f, g, payloads, parent[0], g[-1], ai, undo)
        return SS.inverse_apply(b, g[-1], action, incoming, defer)

    def dots():
        for v in payloads:
            for i in range(n):
                ar.mul_xor0(b, g[i], v[i], v[-1])

    def last(u, inc, delay):
        def action(ai):
            SS.runops(b, [('s', v[-1], True) for v in payloads] + loc_restore_ops(f, g), g[-1], ai, u)
        return SS.inverse_apply(b, g[-1], action, inc, delay)
    if not undo:
        dots()
        return last(False, incoming, defer)
    carry = last(True, incoming, defer)
    dots()
    return carry


# ---------------------------------------------------------------------------------------------- the bottom
def flags_at_bottom(b, ctx):
    """Zero tests of the quotient words q1 of levels 16..7 (those the routes read), held during the bottom pass."""
    need = set()
    for U in CFG['upper']:
        need |= {U - 1, U - 2, U - 3, U - 4}
    fl = {}
    for m in sorted(need):
        fl[m] = SS.fresh_zero_flag(b, ctx['lv'][m]['f'][-1])
    return fl


def clear_flags(b, ctx, fl):
    for m in sorted(fl, reverse=True):
        SS.clear_flag(b, ctx['lv'][m]['f'][-1], fl[m])


def tops_of(ctx, m, base_payloads):
    return [v[0] for v in base_payloads] if m == 1 else [v[-1] for v in ctx['lv'][m]['payloads']]


def route_block(b, ctx, fl, U, base_payloads, transposed):
    """Carry of block (U, U-1) (or its transpose): d_U' ^= ZB X, d_L' ^= selL X, d_U'' ^= selU2 X, where
    ZB = [not B and C] q0_U' + [not B and not C and not F], selL = not B and ((C and not C2) or (not C and F)),
    selU2 = not B and C and C2 and not F (zero region: F(U') and F(L') both set, no carry)."""
    L = U - 1
    Bq, Cq, Fq, C2q = fl[L], fl[U - 2], fl[U - 3], fl[U - 4]
    q0u = ctx['lv'][U - 2]['f'][-2]
    h1, l1 = flag_and(b, [(Bq, True), (Cq, False)])
    ZB = masked_word(b, h1, q0u)
    h2, l2 = flag_and(b, [(Bq, True), (Cq, True), (Fq, True)])
    b.cx(h2, ZB[0])
    # zero region of padded charts: F(U') and F(L') together mean the divisor at U' is 0: no carry there
    sA, lA = flag_and(b, [(Bq, True), (Cq, False), (C2q, True)])
    sB, lB = flag_and(b, [(Bq, True), (Cq, True), (Fq, False), (C2q, True)])
    sU2, lU2 = flag_and(b, [(Bq, True), (Cq, False), (C2q, False), (Fq, True)])
    X = [v[-1] for v in ctx['lv'][U]['payloads']]
    tu, tl, t2 = (tops_of(ctx, U - 2, base_payloads), tops_of(ctx, U - 3, base_payloads),
                  tops_of(ctx, U - 4, base_payloads))
    for k in range(len(X)):
        if not transposed:
            ar.mul_xor0(b, ZB, X[k], tu[k])
            for sel in (sA, sB):
                for j in range(12):
                    b.ccx(sel, X[k][j], tl[k][j])
            for j in range(12):
                b.ccx(sU2, X[k][j], t2[k][j])
        else:
            ar.mul_xor0(b, ZB, tu[k], X[k])
            for sel in (sA, sB):
                for j in range(12):
                    b.ccx(sel, tl[k][j], X[k][j])
            for j in range(12):
                b.ccx(sU2, t2[k][j], X[k][j])
    erase_flag_and(b, sU2, lU2)
    erase_flag_and(b, sB, lB)
    erase_flag_and(b, sA, lA)
    b.cx(h2, ZB[0])
    erase_flag_and(b, h2, l2)
    erase_masked_word(b, h1, q0u, ZB)
    erase_flag_and(b, h1, l1)


def bottom_pass(b, ctx, base_payloads, phase, reverse):
    """Forward: carries (bottom-up), every q_post 1..17, carries transposed (top-down). Reverse: the inverse."""
    ups = sorted(CFG['upper'])
    lv = ctx['lv']

    def carries(transposed, order):
        fl = {}
        for k, U in enumerate(order):                 # sliding window: at most four flags held
            need = {U - 1, U - 2, U - 3, U - 4}
            for m in sorted(set(fl) - need, reverse=True):
                SS.clear_flag(b, lv[m]['f'][-1], fl.pop(m))
            for m in sorted(need - set(fl)):
                fl[m] = SS.fresh_zero_flag(b, lv[m]['f'][-1])
            route_block(b, ctx, fl, U, base_payloads, transposed)
        for m in sorted(fl, reverse=True):
            SS.clear_flag(b, lv[m]['f'][-1], fl.pop(m))

    def qposts(rev):
        for m in (range(17, 0, -1) if rev else range(1, 18)):
            d = lv[m]
            SS.q_post(b, d['f'][-1], d['payloads'], None, m, d['parent'], rev, phase)

    if not reverse:
        carries(False, ups)
        qposts(False)
        carries(True, list(reversed(ups)))
    else:
        carries(True, ups)
        qposts(True)
        carries(False, list(reversed(ups)))


# ---------------------------------------------------------------------------------------------- recursion
def base_pd(b, f, g, payloads, parent, reverse, phase, ctx):
    a = g[0]
    q = f[-1]
    ctx['lv'][1] = dict(f=f, g=g, payloads=payloads, parent=parent)

    def first(undo, inc=None, delay=False):
        def action(ai):
            SS.runops(b, [('s', v[0], True) for v in payloads] + [('s', q, True)], a, ai, undo)
        return SS.inverse_apply(b, a, action, inc, delay)

    def last(undo, inc=None, delay=False):
        def action(ai):
            SS.runops(b, [('s', v[0], True) for v in payloads] + [('s', q, False)], a, ai, undo)
        return SS.inverse_apply(b, a, action, inc, delay)
    if not reverse:
        carry = first(False, delay=True)
        bottom_pass(b, ctx, payloads, phase, False)
        phase['a'][1] = last(False, carry, True)
    else:
        carry = last(True, phase['a'].pop(1), True)
        bottom_pass(b, ctx, payloads, phase, True)
        first(True, carry)


def look_pd(b, f, g, payloads, parent, reverse, phase, ctx):
    m = len(f) - 1
    if m == 1:
        return base_pd(b, f, g, payloads, parent, reverse, phase, ctx)
    ctx['lv'][m] = dict(f=f, g=g, payloads=payloads, parent=parent)
    F = SS.fresh_zero_flag(b, g[-1])
    rotate(b, g, F)
    if not reverse:
        carry = pre_pd(b, f, g, payloads, parent, m, defer=True)
    else:
        carry = post_pd(b, f, g, payloads, parent, m, True, incoming=phase['a'].pop(m), defer=True)
    SS.canonical_permute(b, f, g, F)
    SS.clear_flag(b, g[-1], F)
    child = ([v[-1] for v in payloads], f[-1], f[-2])
    look_pd(b, g, f[:-2], [v[:-1] for v in payloads], child, reverse, phase, ctx)
    F = SS.fresh_zero_flag(b, g[-1])
    SS.canonical_permute(b, f, g, F, True)
    if not reverse:
        phase['a'][m] = post_pd(b, f, g, payloads, parent, m, incoming=carry, defer=True)
    else:
        pre_pd(b, f, g, payloads, parent, m, True, incoming=carry)
    rotate(b, g, F, True)
    SS.clear_flag(b, g[-1], F)


def bez_pd(b, O, E, payloads, reverse=False, phase=None):
    """Drop-in for SS.bez_lookahead (history word allocated here)."""
    phase = {'a': {}, 'q': {}} if phase is None else phase
    ctx = {'lv': {}}

    def reduce():
        for i in range(len(O)):
            ar.mul_xor0(b, E[-1], O[i], E[i])
    reduce()
    lead = list(b.alloc_qubits(12))
    b.x(lead[0])
    look_pd(b, O + [lead], E[:-1], payloads, None, reverse, phase, ctx)
    b.x(lead[0])
    b.free(lead)
    reduce()
    return phase


# ---------------------------------------------------------------------------------------------- test
def build(which, mode):
    b = km.CircuitBuilder()
    regs = dict(O=[b.alloc_qubits(12) for _ in range(17)], E=[b.alloc_qubits(12) for _ in range(18)],
                u=[b.alloc_qubits(12) for _ in range(17)], w=[b.alloc_qubits(12) for _ in range(17)])
    P = [regs['u'], regs['w']]
    if which == 'kit':
        ph = {'a': {}, 'q': {}}
        SS.bez_lookahead(b, regs['O'], regs['E'], P, False, None, ph)
        if mode == 'round':
            SS.bez_lookahead(b, regs['O'], regs['E'], P, True, None, ph)
    else:
        ph = bez_pd(b, regs['O'], regs['E'], P, False)
        if mode == 'round':
            bez_pd(b, regs['O'], regs['E'], P, True, ph)
    return b.finish_circuit(), regs


def clean_check(cr, vals, seed):
    """Phase flips and nonzero ancillas (every qubit outside the 69 register words) after a round trip."""
    import dqi_int_convtest as CT
    c, regs = cr
    n = vals['O'].shape[0]
    sim = km.Simulator(batch_size=n, debug_print=False, counting_operations=False, seed=seed)
    sim.use_same_registers_as(c)
    inreg = set()
    for k in ('O', 'E', 'u', 'w'):
        qs = [q for w in regs[k] for q in w]
        inreg |= {str(q) for q in qs}
        sim.write_across_shots(km.array(qs), CT.to_bits(vals[k]))
    sim.do(c)
    anc = [km.q(i) for i in range(c.num_qubits) if str(km.q(i)) not in inreg]
    dirty = sim.read_across_shots(km.array(anc)).any(axis=0)
    return np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool), dirty


def main():
    import dqi_int_convtest as CT
    import dqi_zn_pd as PD
    import dqi_zn_pd_kit as PKIT
    kinds = (sys.argv[1] if len(sys.argv) > 1 else '0,1,2,3,4').split(',')
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    if len(sys.argv) > 3:
        CFG['upper'] = tuple(int(x) for x in sys.argv[3].split(','))
    circ = {}
    for which in ('kit', 'pd'):
        for mode in ('fwd', 'round'):
            circ[(which, mode)] = build(which, mode)
            c = circ[(which, mode)][0]
            print('%-3s %-5s qubits %d, static Toffoli-like %d' % (which, mode, c.num_qubits,
                                                                    c.operation_statistics().num_toffoli_like),
                  flush=True)
    rng = np.random.default_rng(seed)
    tot = {}

    def sets():
        for kd in kinds:
            if kd != 'g':
                yield from PD.chart_sets([kd], rng)
                continue
            CT.N = 64                                     # every forced gap pattern of tools/dqi_zn_pd_gaps.py
            pats = [{m} for m in range(17, 1, -1)]
            pats += [{m, m - d} for d in (2, 3, 4) for m in range(17, 2 + d, -1)]
            A = np.concatenate([CT.gap_charts(rng, 64, p)[0] for p in pats])
            B = np.concatenate([CT.gap_charts(rng, 64, p)[1] for p in pats])
            yield '52 gap patterns x 64', A, B, None
    for label, A, B, wt in sets():
        O, E, u, w = PD.bez_inputs(A, B)
        vals = dict(O=O, E=E, u=u, w=w)
        out = {}
        for key, (c, regs) in circ.items():
            out[key] = PKIT.run(c, regs, vals, seed)
            tot[key] = tot.get(key, 0.0) + out[key][1]
        kf, pf = out[('kit', 'fwd')][0], out[('pd', 'fwd')][0]
        eq = np.ones(O.shape[0], dtype=bool)
        for k in kf:
            eq &= (kf[k] == pf[k]).all(axis=1)
        pr = out[('pd', 'round')][0]
        back = np.ones(O.shape[0], dtype=bool)
        for k, v in vals.items():
            back &= (pr[k] == v).all(axis=1)
        ph, dirty = clean_check(circ[('pd', 'round')], vals, seed)
        print('%-22s pd forward == kit forward %d/%d; pd round trip values back %d/%d, phase flips %d, dirty '
              'ancillas %d' % (label, eq.sum(), len(eq), back.sum(), len(back), ph.sum(), dirty.sum()), flush=True)
    nb = len(kinds)
    for key in circ:
        print('executed per shot %-3s %-5s %.1f' % (key[0], key[1], tot[key] / nb), flush=True)


if __name__ == '__main__':
    main()
