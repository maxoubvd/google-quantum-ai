#!/usr/bin/env python
"""r37_dqi_ext_rp_phase.py -- re-pricer of R37-DQI-EXT: word-level PHASE simulation of the prefix block with the
designer's patches C (compute, both calls), P (group loop, scheduled as cmd_blockloop) and T (return erasure with Pa
kept through the loop): the part of the design that no builder has built.

The designer's model checks word values only. Here every kit primitive is mocked on 12-bit words with a tracked phase
bit, with kickmix's semantics (research/panel/BRIEF.md): HMR = X-measure with a random outcome m and the phase kick
(-1)^(m . w); cz(w, m) = (-1)^(m . w); phase_product(x, y, m) = (-1)^(m . xy) (CZs only); erase_product(x, y, z) = HMR
of z = xy then phase_product(x, y, m); mul_xor0(x, y, z): z ^= xy, one 42-AND product. Head words, L / N words of the
14 groups and the zero node's words are random per input. The masks of one "mask sample" are the same for all inputs
(the schedule is data-independent), so the block is exact iff (a) every value is as specified at every check and
(b) for each mask sample the final phase bit is the same on every input (a data-dependent phase is a failed shot).

Checks: U words after C against the block's own point chain (prefix5's formulas); L_g after the forward call against
the block identity L0 += U(w_g) X_g^2 / sc_g + G(w_g) / Q1g^2 (U, G evaluated directly in E from prefix5's docstring
polynomials, independent of the point values); every word restored after the return; peak clean words; products.
Modes (python -I -B r37_dqi_ext_rp_phase.py <mode> [inputs] [masks] [seed]):
  full    loop scheduled (cmd_blockloop's ops, zero node joined: 32 products per call)
  fast    loop applied as its specification (same counts), for long sweeps
  neg     negative controls: each run drops ONE phase payment (every k-th phase_product); must FAIL each time
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r37_dqi_ext_chain as ch  # noqa: E402  (the designer's field, nodes and constants, read-only)

mul, inv, emul, escale, esq, econj, fsqrt = ch.mul, ch.inv, ch.emul, ch.escale, ch.esq, ch.econj, ch.fsqrt
NF, MASK, M, KAP, CS, HOM = ch.NF, ch.MASK, ch.M, ch.KAP, ch.CS, ch.HOM
k0, k1 = KAP[0], KAP[1]
k0s, k1s = mul(k0, k0), mul(k1, k1)
BETA = {0: 0, 1: 1, 4: 2}
NODES = [(g, esq(CS[g] | NF)) for g in HOM]


def sq(v):
    return mul(v, v)


def par(x):
    return bin(x).count('1') & 1


class Mock:
    def __init__(self, regs, mrng, drop=0):
        self.R = dict(regs)
        self.tmp = set()
        self.peak = 0
        self.phase = 0
        self.prods = 0
        self.mrng = mrng
        self.drop = drop
        self.npay = 0

    def alloc(self, k):
        assert k not in self.R, k
        self.R[k] = 0
        self.tmp.add(k)
        self.peak = max(self.peak, len(self.tmp))

    def free(self, k):
        assert self.R[k] == 0, ('free dirty', k)
        del self.R[k]
        self.tmp.discard(k)

    def cx(self, s, d, c=1):
        if c:
            self.R[d] ^= mul(self.R[s], c)

    def scale(self, w, c):
        assert c
        self.R[w] = mul(self.R[w], c)

    def sq(self, w, undo=False):
        self.R[w] = fsqrt(self.R[w]) if undo else sq(self.R[w])

    def mulx(self, x, y, z):
        self.R[z] ^= mul(self.R[x], self.R[y])
        self.prods += 1

    def hmr(self, w):
        m = self.mrng.randrange(NF)
        self.phase ^= par(m & self.R[w])
        self.R[w] = 0
        self.free(w)
        return m

    def cz(self, w, m):
        self.phase ^= par(m & self.R[w])

    def pp(self, x, y, m):
        self.npay += 1
        if self.drop and self.npay % self.drop == 0:
            return
        self.phase ^= par(m & mul(self.R[x], self.R[y]))

    def erase(self, x, y, z):
        assert self.R[z] == mul(self.R[x], self.R[y]), ('erase: word is not the product', z)
        m = self.hmr(z)
        self.pp(x, y, m)

    def form(self, a, bb, beta, body):
        """body(reg) while reg holds (a + beta bb)^2 (beta None: bb^2), as prefix5's form()."""
        reg = bb if beta is None else a
        if beta:
            self.cx(bb, a, beta)
        self.sq(reg)
        body(reg)
        self.sq(reg, True)
        if beta:
            self.cx(bb, a, beta)

    def fsq(self, w, body):
        self.sq(w)
        body(w)
        self.sq(w, True)


# ---------------------------------------------------------------- reference formulas (prefix5's docstring, in E) ---
def ref_YGU(h, w):
    """Y, G, U at w in E (or F), directly from the polynomials: Y = a3 + (k1 + w) b3 + b4 (a3^2 + w b3^2),
    G = q1(w)(b2 + b4) + (a2^2 + w b2^2) Y, U = q1(w)(a1 + (k0 + w) b1) + q0(w) Y + (a1^2 + w b1^2) G."""
    a1, b1, a2, b2, a3, b3, b4 = (h[k] for k in ('a1', 'b1', 'a2', 'b2', 'a3', 'b3', 'b4'))
    w2 = esq(w)
    q0, q1 = w2 ^ w ^ k0s, w2 ^ w ^ k1s
    Y = a3 ^ escale(k1 ^ w, b3) ^ emul(b4, sq(a3) ^ escale(w, sq(b3)))
    G = emul(q1, b2 ^ b4) ^ emul(sq(a2) ^ escale(w, sq(b2)), Y)
    U = emul(q1, a1 ^ emul(k0 ^ w, b1)) ^ emul(q0, Y) ^ emul(sq(a1) ^ escale(w, sq(b1)), G)
    return Y, G, U


def ref_point(h, p):
    """The block's point chain (prefix5 temp_point / point0), alpha = p in F or None (leading coefficient)."""
    a1, b1, a2, b2, a3, b3, b4 = (h[k] for k in ('a1', 'b1', 'a2', 'b2', 'a3', 'b3', 'b4'))
    if p is None:
        Y = b3 ^ mul(b4, sq(b3))
        Gv = b2 ^ b4 ^ mul(Y, sq(b2))
        return b1 ^ Y ^ mul(Gv, sq(b1))
    beta = BETA[p]
    x, y, s = a1 ^ mul(beta, b1), a2 ^ mul(beta, b2), a3 ^ mul(beta, b3)
    q0, q1 = sq(p) ^ p ^ k0s, sq(p) ^ p ^ k1s
    Y = a3 ^ mul(k1 ^ p, b3) ^ mul(b4, sq(s))
    Gv = mul(q1, b2 ^ b4) ^ mul(sq(y), Y)
    return mul(q1, a1 ^ mul(k0 ^ p, b1)) ^ mul(q0, Y) ^ mul(sq(x), Gv)


# ---------------------------------------------------------------- the loop (patch P), constants cached ------------
def loop_plan():
    allnodes = [w for _, w in NODES] + [econj(w) for _, w in NODES] + [0]          # zero node joined
    pts = list(range(len(allnodes) + 2))
    fin = [0, 1, 4]

    def lam_x(node, t):
        num = den = 1
        for nu in allnodes:
            if nu != node:
                if t is not None:
                    num = emul(num, t ^ nu)
                den = emul(den, node ^ nu)
        return emul(num, ch.einv(den))

    def lam_u(t):
        out, allp = {}, 1
        for p in fin:
            num, den = 1, 1
            for r in fin:
                if r != p:
                    num = mul(num, t ^ r)
                    den = mul(den, p ^ r)
            out[p] = mul(num, inv(den))
            allp = mul(allp, t ^ p)
        out[None] = allp
        return out

    def spread_c(t, w):
        if t is None:
            out = 1
            for s in pts:
                out = emul(out, w ^ s)
            return out
        num, den = 1, 1
        for s in pts:
            if s != t:
                num = emul(num, w ^ s)
                den = mul(den, t ^ s)
        return escale(num, inv(den))

    sc = {g: mul(sq(k0 ^ KAP[g]), sq(k1 ^ KAP[g])) for g, _ in NODES}
    zc = mul(k0s, k1s)
    plan = []
    for t in pts + [None]:
        if t is None:
            ukey, lu = ('U', None), None
        elif t in fin:
            ukey, lu = ('U', t), None
        else:
            ukey, lu = ('U', 1), lam_u(t)
        cf = {}
        for g, w in NODES:
            lam = lam_x(w, t)
            cf[('N', g, 0)] = lam >> M
            cf[('N', g, 1)] = (lam & MASK) ^ (lam >> M)
        cf[('N', 'z', 0)] = lam_x(0, t) & MASK
        ykey = next(k for k in cf if cf[k])
        cc = {}
        for g, w in NODES:
            c = escale(spread_c(t, w), inv(sc[g]))
            cc[('L', g, 0)], cc[('L', g, 1)] = c & MASK, c >> M
        cc[('L', 'z', 0)] = mul(spread_c(t, 0) & MASK, inv(zc))
        tkey = next(k for k in cc if cc[k])
        plan.append((ukey, lu, cf, ykey, cc, tkey))
    return plan, sc, zc


PLAN, SC, ZC = loop_plan()


def loop_sched(mk):
    """Patch P exactly as cmd_blockloop: expose U(t) on a U word, XX(t) on one N word, spread, one product, undo."""
    for ukey, lu, cf, ykey, cc, tkey in PLAN:
        if lu is not None:
            mk.scale(ukey, lu[1])
            for p in (0, 4, None):
                mk.cx(('U', p), ukey, lu[p])
        mk.scale(ykey, cf[ykey])
        for k in cf:
            if k != ykey:
                mk.cx(k, ykey, cf[k])
        mk.scale(tkey, inv(cc[tkey]))
        for k in cc:
            if k != tkey:
                mk.cx(tkey, k, cc[k])
        mk.mulx(ukey, ykey, tkey)
        for k in cc:
            if k != tkey:
                mk.cx(tkey, k, cc[k])
        mk.scale(tkey, cc[tkey])
        for k in cf:
            if k != ykey:
                mk.cx(k, ykey, cf[k])
        mk.scale(ykey, inv(cf[ykey]))
        if lu is not None:
            for p in (0, 4, None):
                mk.cx(('U', p), ukey, lu[p])
            mk.scale(ukey, inv(lu[1]))


def loop_spec(mk):
    """The loop's specification (the same 32 products counted): L_g ^= U(w_g) X_g^2 / sc_g, zL ^= U(0) zX / zc."""
    R = mk.R
    Upt = {p: R[('U', p)] for p in (0, 1, 4, None)}
    for g, w in NODES:
        X2 = R[('N', g, 0)] | (R[('N', g, 1)] << M)
        v = escale(emul(u_of(Upt, w), X2), inv(SC[g]))
        R[('L', g, 0)] ^= v & MASK
        R[('L', g, 1)] ^= v >> M
    R[('L', 'z', 0)] ^= mul(mul(Upt[0], R[('N', 'z', 0)]), inv(ZC))
    mk.prods += len(PLAN)


def u_of(Upt, w):
    """The cubic through U(0), U(1), U(4) with leading coefficient U(inf), at w in E."""
    out = 0
    allp = 1
    for p in (0, 1, 4):
        num, den = 1, 1
        for r in (0, 1, 4):
            if r != p:
                num = emul(num, w ^ r)
                den = mul(den, p ^ r)
        out ^= emul(escale(num, inv(den)), Upt[p])
        allp = emul(allp, w ^ p)
    return out ^ emul(allp, Upt[None])


# ---------------------------------------------------------------- patch C (the designer's 9 steps) ----------------
FANIN = {}
for _g, _w in NODES:
    _c0 = inv(sq(k1 ^ KAP[_g]))
    FANIN[_g] = (_c0, escale(_w, _c0), escale(esq(_w), _c0))


def g_fanin(mk):
    """The block's G fan-ins: L_g ^= G(w_g) / Q1g^2 = (g0 + w_g g1 + w_g^2 g2) / (k1 + kap_g)^2; zL ^= g0 / k1^2."""
    for g, _ in NODES:
        c0, c1, c2 = FANIN[g]
        mk.cx('g0', ('L', g, 0), c0)
        mk.cx('g1', ('L', g, 0), c1 & MASK)
        mk.cx('g1', ('L', g, 1), c1 >> M)
        mk.cx('g2', ('L', g, 0), c2 & MASK)
        mk.cx('g2', ('L', g, 1), c2 >> M)
    mk.cx('g0', ('L', 'z', 0), inv(k1s))


def coef_words(mk):
    """Pa, Pb (y0, y1 exposed) and V0, V2, V1 (g0, g1, g2 with their linear parts): 5 products."""
    mk.alloc('Pa')
    mk.fsq('a3', lambda r: mk.mulx('b4', r, 'Pa'))
    mk.alloc('Pb')
    mk.fsq('b3', lambda r: mk.mulx('b4', r, 'Pb'))
    mk.cx('a3', 'Pa')
    mk.cx('b3', 'Pa', k1)                      # y0
    mk.cx('b3', 'Pb')                          # y1
    mk.alloc('g0')
    mk.fsq('a2', lambda r: mk.mulx('Pa', r, 'g0'))
    mk.alloc('g2')
    mk.fsq('b2', lambda r: mk.mulx('Pb', r, 'g2'))
    mk.alloc('g1')
    mk.cx('b2', 'a2')
    mk.cx('Pb', 'Pa')
    mk.fsq('a2', lambda r: mk.mulx('Pa', r, 'g1'))
    mk.cx('Pb', 'Pa')
    mk.cx('b2', 'a2')
    mk.cx('g0', 'g1')
    mk.cx('g2', 'g1')                          # V1
    mk.cx('b2', 'g0', k1s)
    mk.cx('b4', 'g0', k1s)
    for gw in ('g1', 'g2'):
        mk.cx('b2', gw)
        mk.cx('b4', gw)


def c_compute(mk, keep_pa):
    """Patch C. keep_pa (return call): step 9 runs before step 8 and Pa stays live as y0 for patch T."""
    coef_words(mk)                                               # steps 1, 2
    g_fanin(mk)
    mk.cx('b3', 'Pb')                                            # step 3: erase Pb, Pa (free)
    mk.fsq('b3', lambda r: mk.erase('b4', r, 'Pb'))
    mk.cx('b3', 'Pa', k1)
    mk.cx('a3', 'Pa')
    mk.fsq('a3', lambda r: mk.erase('b4', r, 'Pa'))
    for p in (1, 4):                                             # step 4
        U = ('U', p)
        mk.alloc(U)
        q1 = sq(p) ^ p ^ k1s
        mk.cx('a1', U, q1)
        mk.cx('b1', U, mul(q1, k0 ^ p))
        mk.scale('g1', p)
        mk.cx('g0', 'g1')
        mk.cx('g2', 'g1', sq(p))                                 # G(p) on g1
        mk.form('a1', 'b1', BETA[p], lambda r: mk.mulx('g1', r, U))
        mk.cx('g2', 'g1', sq(p))
        mk.cx('g0', 'g1')
        mk.scale('g1', inv(p))
    m1 = mk.hmr('g1')                                            # step 5: g1 = b2 + b4 + a2^2 y1 + b2^2 y0
    mk.cz('b2', m1)
    mk.cz('b4', m1)
    mk.alloc(('U', 0))                                           # step 6
    mk.cx('a1', ('U', 0), k1s)
    mk.cx('b1', ('U', 0), mul(k1s, k0))
    mk.fsq('a1', lambda r: mk.mulx('g0', r, ('U', 0)))
    m0 = mk.hmr('g0')                                            # g0 = k1^2 (b2 + b4) + a2^2 y0
    mk.cx('b4', 'b2')
    mk.scale('b2', k1s)
    mk.cz('b2', m0)
    mk.scale('b2', inv(k1s))
    mk.cx('b4', 'b2')
    mk.alloc(('U', None))                                        # step 7
    mk.cx('b1', ('U', None))
    mk.fsq('b1', lambda r: mk.mulx('g2', r, ('U', None)))
    m2 = mk.hmr('g2')                                            # g2 = b2 + b4 + b2^2 y1
    mk.cz('b2', m2)
    mk.cz('b4', m2)

    def step8(keep):                                             # Pa again: y0 live
        mk.alloc('Pa')
        mk.fsq('a3', lambda r: mk.mulx('b4', r, 'Pa'))
        mk.cx('a3', 'Pa')
        mk.cx('b3', 'Pa', k1)
        for p in (0, 1, 4):
            mk.cx('Pa', ('U', p), sq(p) ^ p ^ k0s)
        mk.fsq('a2', lambda r: mk.pp(r, 'Pa', m0))               # m0's a2^2 y0
        mk.fsq('b2', lambda r: mk.pp(r, 'Pa', m1))               # m1's b2^2 y0
        if not keep:
            mk.cx('b3', 'Pa', k1)
            mk.cx('a3', 'Pa')
            mk.fsq('a3', lambda r: mk.erase('b4', r, 'Pa'))

    def step9():                                                 # Pb again: y1 live
        mk.alloc('Pb')
        mk.fsq('b3', lambda r: mk.mulx('b4', r, 'Pb'))
        mk.cx('b3', 'Pb')
        for p in (1, 4):
            mk.cx('Pb', ('U', p), mul(sq(p) ^ p ^ k0s, p))
        mk.cx('Pb', ('U', None))
        mk.fsq('b2', lambda r: mk.pp(r, 'Pb', m2))               # m2's b2^2 y1
        mk.fsq('a2', lambda r: mk.pp(r, 'Pb', m1))               # m1's a2^2 y1
        mk.cx('b3', 'Pb')
        mk.fsq('b3', lambda r: mk.erase('b4', r, 'Pb'))

    if keep_pa:
        step9()
        step8(True)
    else:
        step8(False)
        step9()


def t_erase(mk):
    """Patch T with Pa (y0) kept through the loop: X-measure the four U words, Pb + V0, V2, V1 (4 products), pay
    every mask (linear part by CZ on an exposed word, h_1(p)^2 G(p) by one phase_product), erase all (free)."""
    masks = {p: mk.hmr(('U', p)) for p in (None, 4, 1, 0)}
    mk.alloc('Pb')
    mk.fsq('b3', lambda r: mk.mulx('b4', r, 'Pb'))
    mk.cx('b3', 'Pb')                                            # y1
    mk.alloc('g0')
    mk.fsq('a2', lambda r: mk.mulx('Pa', r, 'g0'))
    mk.alloc('g2')
    mk.fsq('b2', lambda r: mk.mulx('Pb', r, 'g2'))
    mk.alloc('g1')
    mk.cx('b2', 'a2')
    mk.cx('Pb', 'Pa')
    mk.fsq('a2', lambda r: mk.mulx('Pa', r, 'g1'))
    mk.cx('Pb', 'Pa')
    mk.cx('b2', 'a2')
    mk.cx('g0', 'g1')
    mk.cx('g2', 'g1')
    mk.cx('b2', 'g0', k1s)
    mk.cx('b4', 'g0', k1s)
    for gw in ('g1', 'g2'):
        mk.cx('b2', gw)
        mk.cx('b4', gw)
    for p in (0, 1, 4, None):
        m = masks[p]
        if p is None:                                            # U(inf) = b1 + y1 + b1^2 g2
            mk.cx('b1', 'Pb')
            mk.cz('Pb', m)
            mk.cx('b1', 'Pb')
            mk.fsq('b1', lambda r: mk.pp(r, 'g2', m))
            continue
        q0, q1 = sq(p) ^ p ^ k0s, sq(p) ^ p ^ k1s                # U(p) = q1 (a1 + (k0+p) b1) + q0 Y(p) + x^2 G(p)
        mk.scale('Pa', q0)
        mk.cx('a1', 'Pa', q1)
        mk.cx('b1', 'Pa', mul(q1, k0 ^ p))
        mk.cx('Pb', 'Pa', mul(q0, p))
        mk.cz('Pa', m)
        mk.cx('Pb', 'Pa', mul(q0, p))
        mk.cx('b1', 'Pa', mul(q1, k0 ^ p))
        mk.cx('a1', 'Pa', q1)
        mk.scale('Pa', inv(q0))
        mk.cx('g1', 'g0', p)
        mk.cx('g2', 'g0', sq(p))
        mk.form('a1', 'b1', BETA[p], lambda r: mk.pp(r, 'g0', m))
        mk.cx('g2', 'g0', sq(p))
        mk.cx('g1', 'g0', p)
    for gw in ('g1', 'g2'):
        mk.cx('b4', gw)
        mk.cx('b2', gw)
    mk.cx('b4', 'g0', k1s)
    mk.cx('b2', 'g0', k1s)
    mk.cx('g0', 'g1')
    mk.cx('g2', 'g1')
    mk.cx('b2', 'a2')
    mk.cx('Pb', 'Pa')
    mk.fsq('a2', lambda r: mk.erase('Pa', r, 'g1'))
    mk.cx('Pb', 'Pa')
    mk.cx('b2', 'a2')
    mk.fsq('a2', lambda r: mk.erase('Pa', r, 'g0'))
    mk.fsq('b2', lambda r: mk.erase('Pb', r, 'g2'))
    mk.cx('b3', 'Pa', k1)
    mk.cx('a3', 'Pa')
    mk.fsq('a3', lambda r: mk.erase('b4', r, 'Pa'))
    mk.cx('b3', 'Pb')
    mk.fsq('b3', lambda r: mk.erase('b4', r, 'Pb'))


# ---------------------------------------------------------------- one decode: forward block, return block ---------
HEADS = ('a1', 'b1', 'a2', 'b2', 'a3', 'b3', 'b4')


def instance_regs(rng):
    """Block inputs from a weight-35 decoder instance (the designer's instance(): jets + reference heads): heads 1..3
    and b_4, N_g = (x_g + H_0(gamma_g))^2, L_g = the jet's v word; the zero node's words random."""
    J, heads = ch.instance(rng)
    R = {}
    (R['a1'], R['b1']), (R['a2'], R['b2']), (R['a3'], R['b3']) = heads[1], heads[2], heads[3]
    R['b4'] = heads[4][1]
    a0, b0 = heads[0]
    for g, _ in NODES:
        X2 = esq(J[g][0] ^ ch.head_value(a0, b0, CS[g]))
        R[('N', g, 0)], R[('N', g, 1)] = X2 & MASK, X2 >> M
        R[('L', g, 0)], R[('L', g, 1)] = J[g][1] & MASK, J[g][1] >> M
    R[('N', 'z', 0)] = rng.randrange(NF)
    R[('L', 'z', 0)] = rng.randrange(NF)
    return R


def random_regs(rng, special):
    if special == 'inst':
        return instance_regs(rng)
    R = {k: rng.randrange(NF) for k in HEADS}
    if special:                                   # structured inputs: some head words zero
        for k in HEADS:
            if rng.random() < 0.3:
                R[k] = 0
    for g, _ in NODES:
        for hh in (0, 1):
            R[('N', g, hh)] = rng.randrange(NF)
            R[('L', g, hh)] = rng.randrange(NF)
    R[('N', 'z', 0)] = rng.randrange(NF)
    R[('L', 'z', 0)] = rng.randrange(NF)
    return R


def one_decode(R0, mseed, loop, drop=0):
    mk = Mock(R0, random.Random(mseed), drop)
    h = {k: R0[k] for k in HEADS}
    stats = {}
    # forward call
    p0 = mk.prods
    c_compute(mk, keep_pa=False)
    for p in (0, 1, 4, None):
        assert mk.R[('U', p)] == ref_point(h, p), ('U after C (forward)', p)
    loop(mk)
    for p in (0, 1, 4, None):
        assert mk.R[('U', p)] == ref_point(h, p), ('U after loop', p)
    fwd = {p: mk.hmr(('U', p)) for p in (None, 4, 1, 0)}       # the forward's deferral (masks to the return)
    stats['fwd_products'] = mk.prods - p0
    stats['fwd_peak'] = mk.peak
    for g, w in NODES:                                           # block identity on the L words
        Y, G, U = ref_YGU(h, w)
        X2 = R0[('N', g, 0)] | (R0[('N', g, 1)] << M)
        want = (R0[('L', g, 0)] | (R0[('L', g, 1)] << M)) ^ escale(emul(U, X2), inv(SC[g])) \
            ^ escale(G, inv(sq(k1 ^ KAP[g])))
        assert (mk.R[('L', g, 0)] | (mk.R[('L', g, 1)] << M)) == want, ('L after forward', g)
    Y0, G0, U0 = ref_YGU(h, 0)
    assert mk.R[('L', 'z', 0)] == R0[('L', 'z', 0)] ^ mul(mul(U0, R0[('N', 'z', 0)]), inv(ZC)) ^ mul(G0, inv(k1s))
    for k in HEADS:
        assert mk.R[k] == R0[k], ('head after forward', k)
    assert not mk.tmp, mk.tmp
    # return call
    mk.peak = 0
    p0 = mk.prods
    c_compute(mk, keep_pa=True)
    for p in (0, 1, 4, None):
        assert mk.R[('U', p)] == ref_point(h, p), ('U after C (return)', p)
        mk.cz(('U', p), fwd[p])                                  # pay the forward's measured U words
    loop(mk)
    t_erase(mk)
    stats['ret_products'] = mk.prods - p0
    stats['ret_peak'] = mk.peak
    assert not mk.tmp, mk.tmp
    for k in R0:
        assert mk.R[k] == R0[k], ('not restored', k)
    return mk.phase, stats


def run(mode, inputs, masks, seed):
    loop = loop_sched if mode in ('full', 'inst') else loop_spec
    rng = random.Random(seed)
    bad_phase = 0
    stats = None
    for s in range(masks):
        mseed = rng.randrange(1 << 30)
        ph = set()
        for i in range(inputs):
            R0 = random_regs(rng, special='inst' if mode == 'inst' else (i % 4 == 3))
            phase, stats = one_decode(R0, mseed, loop)
            ph.add(phase)
        bad_phase += len(ph) > 1
    print(f'{mode}: {masks} mask samples x {inputs} inputs (every 4th with zeroed head words); value checks all '
          f'passed (U after C both calls, L block identity, full restoration); data-dependent phase in {bad_phase} '
          f'of {masks} mask samples; per call: forward {stats["fwd_products"]} products, peak {stats["fwd_peak"]} '
          f'clean words; return {stats["ret_products"]} products, peak {stats["ret_peak"]} clean words')
    print(f'{mode.upper()} PASS' if not bad_phase else f'{mode.upper()} FAIL')


def run_neg(inputs, masks, seed):
    """Negative controls: drop every k-th phase payment (k = 1..npay); each must show a data-dependent phase."""
    rng = random.Random(seed)
    R0 = random_regs(rng, False)
    cnt = Mock(R0, random.Random(1))                             # count the payments of one decode
    c_compute(cnt, False)
    loop_spec(cnt)
    fwdm = {p: cnt.hmr(('U', p)) for p in (None, 4, 1, 0)}
    c_compute(cnt, True)
    for p in (0, 1, 4, None):
        cnt.cz(('U', p), fwdm[p])
    loop_spec(cnt)
    t_erase(cnt)
    total = cnt.npay
    caught = 0
    for k in range(1, total + 1):
        hit = False
        for s in range(masks):
            mseed = rng.randrange(1 << 30)
            ph = set()
            for i in range(inputs):
                Rr = random_regs(rng, False)
                ph.add(one_decode_drop(Rr, mseed, k))
            if len(ph) > 1:
                hit = True
                break
        caught += hit
    print(f'neg: {total} phase payments per decode; dropping payment k alone (k = 1..{total}) gives a '
          f'data-dependent phase for {caught} of {total}')
    print('NEG PASS (harness detects every dropped payment)' if caught == total else 'NEG FAIL (harness blind)')


def one_decode_drop(R0, mseed, k):
    mk = Mock(R0, random.Random(mseed))

    def pp(x, y, m):
        mk.npay += 1
        if mk.npay == k:
            return
        mk.phase ^= par(m & mul(mk.R[x], mk.R[y]))
    mk.pp = pp
    c_compute(mk, False)
    loop_spec(mk)
    fwd = {p: mk.hmr(('U', p)) for p in (None, 4, 1, 0)}
    c_compute(mk, True)
    for p in (0, 1, 4, None):
        mk.cz(('U', p), fwd[p])
    loop_spec(mk)
    t_erase(mk)
    return mk.phase


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'fast'
    inputs = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    masks = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    if mode == 'neg':
        run_neg(inputs, masks, seed)
    else:
        run(mode, inputs, masks, seed)
