#!/usr/bin/env python3
"""dqi_int_j6sched.py -- workstream dqi_int (B5): the prefix block at J = 6, closed forms and schedule model.

Part 1 (mode `forms`): the J = 6 correction as an adjoint (Horner) recurrence over polynomials in w = gamma_g^2:
    Z  = Q2^2 b4 + h4 mu5                      (degree 2)
    Y' = Q2^2 mu3 + Q1^2 mu5 + h3 Z            (degree 3)
    G' = Q1^2 (Q2^2 b2 + Z) + h2 Y'            (degree 4)
    U' = Q1^2 Q2^2 mu1 + Q0^2 Y' + h1 G'       (degree 5)
with h_i = a_i^2 + w b_i^2, mu_i = a_i + (k + w) b_i (k = kappa of the head's node), Qk^2 = w^2 + w + k_k^2; the block
adds U'(w_g) X^2 / (Q0^2 Q1^2 Q2^2) + G'(w_g) / (Q1^2 Q2^2) to L0. Checked against dqi_int_j6model.correction (the
affine tracking that is itself checked against the true chain). J = 5 (Z -> b4, no head 5) is checked the same way.
Usage: python dqi_int_j6sched.py forms [trials] [seed]
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def sq(v):
    from dqi_research_model import mul
    return mul(v, v)


def polys(h, w, j6=True):
    """(U', G') (J = 6) or (U, G) (J = 5) at w (an E element); h = {name: F word}."""
    from dqi_research_model import mul, emul, esq, escale, KAP
    k0, k1, k2 = KAP[0], KAP[1], KAP[2]
    w2 = esq(w)
    Q0s, Q1s, Q2s = w2 ^ w ^ sq(k0), w2 ^ w ^ sq(k1), w2 ^ w ^ sq(k2)

    def hh(i):
        return sq(h['a%d' % i]) ^ escale(w, sq(h['b%d' % i]))

    def mu(i, k):
        return h['a%d' % i] ^ mul(k, h['b%d' % i]) ^ escale(w, h['b%d' % i])

    if j6:
        Z = escale(Q2s, h['b4']) ^ emul(hh(4), mu(5, k2))
        Y = emul(Q2s, mu(3, k1)) ^ emul(Q1s, mu(5, k2)) ^ emul(hh(3), Z)
        Gp = emul(Q1s, escale(Q2s, h['b2']) ^ Z) ^ emul(hh(2), Y)
        U = emul(emul(Q1s, Q2s), mu(1, k0)) ^ emul(Q0s, Y) ^ emul(hh(1), Gp)
    else:
        Y = mu(3, k1) ^ escale(hh(3), h['b4'])
        Gp = escale(Q1s, h['b2'] ^ h['b4']) ^ emul(hh(2), Y)
        U = emul(Q1s, mu(1, k0)) ^ emul(Q0s, Y) ^ emul(hh(1), Gp)
    return U, Gp


def cmd_forms(trials, rng):
    from dqi_research_model import mul, esq, escale, KAP, GAMMA, G
    from dqi_int_j6model import correction
    bad = {True: 0, False: 0}
    n = 0
    for _ in range(trials):
        heads = [(rng.randrange(4096), rng.randrange(4096)) for _ in range(2 * G)]
        h = {}
        for i in range(1, 6):
            h['a%d' % i], h['b%d' % i] = heads[i]
        for g in range(3, G):
            w = esq(GAMMA[g])
            q = [sq(KAP[k] ^ KAP[g]) for k in range(3)]
            for j6 in (True, False):
                S = {1, 2, 3, 4, 5} if j6 else {1, 2, 3, 4}
                M, V = correction(heads[:2 * g], g, S)
                U, Gp = polys(h, w, j6)
                if j6:
                    okU = U == escale(M, mul(mul(q[0], q[1]), q[2]))
                    okG = Gp == escale(V, mul(q[1], q[2]))
                else:
                    okU = U == escale(M, mul(q[0], q[1]))
                    okG = Gp == escale(V, q[1])
                bad[j6] += (not okU) + (not okG)
            n += 1
    print(f'FORMS: {trials} head sets x 14 groups = {n}; mismatches J=6 {bad[True]}, J=5 {bad[False]}')
    return bad[True] == 0 and bad[False] == 0


# ===================================================================================================================
# Part 2 (mode `sched`): schedule B5-1 on mock kit primitives (kickmix phase semantics, as tools/r37_dqi_ext_rp_phase.py)
# Field and nodes of the designer's model (r37_dqi_ext_chain, read-only), the same field as dqi_research_model.
# ===================================================================================================================
import r37_dqi_ext_chain as ch  # noqa: E402

cmul, cinv, cemul, cescale, cesq, cfsqrt = ch.mul, ch.inv, ch.emul, ch.escale, ch.esq, ch.fsqrt
NF, MASK, MB, CKAP, CCS, HOM = ch.NF, ch.MASK, ch.M, ch.KAP, ch.CS, ch.HOM
K0, K1, K2 = CKAP[0], CKAP[1], CKAP[2]
NODES = [(g, cesq(CCS[g] | NF)) for g in HOM]
HEADS = ('a1', 'b1', 'a2', 'b2', 'a3', 'b3', 'a4', 'b4', 'a5', 'b5')
BETAS = (1, 2, 3, 5, 7)                         # clean points p = beta^2: p1..p4 by 2-word chains, p5 in-place chain


def csq(v):
    return cmul(v, v)


def qs(k, p):                                   # Qk^2 at an F point p
    return csq(p) ^ p ^ csq(k)


def par(x):
    return bin(x).count('1') & 1


def tmask(m, c):
    """m' with m' . x = m . (c x) for every x (the transpose of multiplication by c)."""
    out = 0
    for j in range(MB):
        out |= par(m & cmul(c, 1 << j)) << j
    return out


class Mock:
    def __init__(self, regs, mrng):
        self.R = dict(regs)
        self.tmp = set()
        self.peak = 0
        self.phase = 0
        self.prods = 0
        self.mrng = mrng

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
            self.R[d] ^= cmul(self.R[s], c)

    def scale(self, w, c):
        assert c
        self.R[w] = cmul(self.R[w], c)

    def sq(self, w, undo=False):
        self.R[w] = cfsqrt(self.R[w]) if undo else csq(self.R[w])

    def mulx(self, x, y, z):
        assert len({x, y, z}) == 3
        self.R[z] ^= cmul(self.R[x], self.R[y])
        self.prods += 1

    def hmr(self, w):
        m = self.mrng.randrange(NF)
        self.phase ^= par(m & self.R[w])
        self.R[w] = 0
        self.free(w)
        return m

    def cz(self, w, m):
        self.phase ^= par(m & self.R[w])

    def pp(self, x, y, m):                       # phase_product: CZs only
        self.phase ^= par(m & cmul(self.R[x], self.R[y]))

    def czc(self, w, m, c):                      # (-1)^(m . c w): in the circuit, scale w by c, CZ, scale back
        if c:
            self.cz(w, tmask(m, c))

    def ppc(self, x, y, m, c):                   # (-1)^(m . c x y): in the circuit, scale y by c around the product
        if c:
            self.pp(x, y, tmask(m, c))

    def mfree(self, m):                          # a mask register no longer needed (the circuit frees its bits)
        pass

    def loop(self):
        loop_spec(self)

    def form(self, a, bb, beta, body):
        """body(a) while a holds (a + beta bb)^2."""
        self.cx(bb, a, beta)
        self.sq(a)
        body(a)
        self.sq(a, True)
        self.cx(bb, a, beta)

    def cz_lin(self, terms, m):
        """(-1)^(m . sum c_r r): the linear form exposed in place on its first register (free), CZ, undone."""
        terms = [(r, c) for r, c in terms if c]
        if not terms or not m:
            return
        (h, c0), rest = terms[0], terms[1:]
        self.scale(h, c0)
        for r, c in rest:
            self.cx(r, h, c)
        self.cz(h, m)
        for r, c in rest:
            self.cx(r, h, c)
        self.scale(h, cinv(c0))


# ---- linear parts of the chain at an F point p (lists of (register, constant)) -------------------------------------
def mu_terms(i, k, p, c=1):                     # c * mu_i(p) = c (a_i + (k + p) b_i)
    return [('a%d' % i, c), ('b%d' % i, cmul(c, k ^ p))]


def c4_terms(p):
    return [('b4', qs(K2, p))]


def c3_terms(p):
    return mu_terms(3, K1, p, qs(K2, p)) + mu_terms(5, K2, p, qs(K1, p))


def c2_terms(p):
    return [('b2', cmul(qs(K1, p), qs(K2, p)))]


def c1_terms(p):
    return mu_terms(1, K0, p, cmul(qs(K1, p), qs(K2, p)))


def add_lin(mk, terms, dst):
    for r, c in terms:
        mk.cx(r, dst, c)


def with_mu5(mk, p, body):
    """body() while a5 holds mu5(p) = a5 + (k2 + p) b5."""
    mk.cx('b5', 'a5', K2 ^ p)
    body()
    mk.cx('b5', 'a5', K2 ^ p)


def lagr(points, x, emode=False):
    """Lagrange weights {t: L_t(x)} over F points; x in F (emode False) or E (emode True)."""
    out = {}
    for t in points:
        num, den = 1, 1
        for r in points:
            if r != t:
                num = cemul(num, x ^ r) if emode else cmul(num, x ^ r)
                den = cmul(den, t ^ r)
        out[t] = cescale(num, cinv(den)) if emode else cmul(num, cinv(den))
    return out


PTS = [csq(b) for b in BETAS]
BETA_OF = {csq(b): b for b in BETAS}
BETA_OF[0] = 0
GSET = PTS[:5]                                  # G' interpolation set for the fan-in (degree 4)
UPTS = PTS[:5] + [0]                            # the six U' words: five clean points and U'(0) in b1
SC = {g: cmul(cmul(csq(K0 ^ CKAP[g]), csq(K1 ^ CKAP[g])), csq(K2 ^ CKAP[g])) for g, _ in NODES}
GC = {g: cinv(cmul(csq(K1 ^ CKAP[g]), csq(K2 ^ CKAP[g]))) for g, _ in NODES}
ZC = cmul(cmul(csq(K0), csq(K1)), csq(K2))
ZG = cinv(cmul(csq(K1), csq(K2)))
FAN = {g: {t: cescale(v, GC[g]) for t, v in lagr(GSET, w, True).items()} for g, w in NODES}
FANZ = {t: cmul(v, ZG) for t, v in lagr(GSET, 0).items()}


def g_fanin(mk, t, word):
    """L_g ^= (weight of G'(t) in G'(w_g)) / (Q1g^2 Q2g^2) for every g; zL ^= (weight at 0) / (k1^2 k2^2)."""
    for g, _ in NODES:
        c = FAN[g][t]
        mk.cx(word, ('L', g, 0), c & MASK)
        mk.cx(word, ('L', g, 1), c >> MB)
    mk.cx(word, ('L', 'z', 0), FANZ[t])


# ---- the chains ------------------------------------------------------------------------------------------------------
def chain2(mk, p, tag, fwd_gmask=None, fanin=True):
    """2-word chain at p: W1 = Z(p) -> G'(p), W2 = Y'(p) -> U'(p); W1 measured (mask returned). 4 products."""
    beta = BETA_OF[p]
    W1, W2 = ('W1', tag), ('U', p)
    mk.alloc(W1)
    add_lin(mk, c4_terms(p), W1)
    with_mu5(mk, p, lambda: mk.form('a4', 'b4', beta, lambda r: mk.mulx(r, 'a5', W1)))      # Z(p)
    mk.alloc(W2)
    add_lin(mk, c3_terms(p), W2)
    mk.form('a3', 'b3', beta, lambda r: mk.mulx(r, W1, W2))                                 # Y'(p)
    mk.scale(W1, qs(K1, p))
    add_lin(mk, c2_terms(p), W1)
    mk.form('a2', 'b2', beta, lambda r: mk.mulx(r, W2, W1))                                 # G'(p)
    if fanin and p in GSET:
        g_fanin(mk, p, W1)
    mk.scale(W2, qs(K0, p))
    add_lin(mk, c1_terms(p), W2)
    mk.form('a1', 'b1', beta, lambda r: mk.mulx(r, W1, W2))                                 # U'(p)
    if fwd_gmask is not None:
        mk.cz(W1, fwd_gmask)                     # the forward's measured G'(p), paid on the same value
        mk.mfree(fwd_gmask)
    return mk.hmr(W1)


def chain_inplace_clean(mk, p):
    """The 5th clean point: Z(p) -> b4, Y'(p) -> b3 in place (roots kept in a4, a3), G'(p) over Z in b4, U'(p) into the
    clean word, then the three in-place stages undone. 7 products."""
    beta = BETA_OF[p]
    q0, q1, q2 = qs(K0, p), qs(K1, p), qs(K2, p)
    cy = cmul(q2, beta ^ K1 ^ p)
    assert cy and q2
    mk.cx('b4', 'a4', beta)                      # t4 = a4 + beta b4 kept in a4
    mk.scale('b4', q2)
    with_mu5(mk, p, lambda: mk.form('a4', 'b4', 0, lambda r: mk.mulx(r, 'a5', 'b4')))       # b4 = Z(p)
    mk.cx('b3', 'a3', beta)                      # s = a3 + beta b3 kept in a3
    mk.scale('b3', cy)
    mk.cx('a3', 'b3', q2)
    add_lin(mk, mu_terms(5, K2, p, q1), 'b3')
    mk.form('a3', 'b3', 0, lambda r: mk.mulx(r, 'b4', 'b3'))                                # b3 = Y'(p)
    mk.scale('b4', q1)
    add_lin(mk, c2_terms(p), 'b4')
    mk.form('a2', 'b2', beta, lambda r: mk.mulx(r, 'b3', 'b4'))                             # b4 = G'(p)
    if p in GSET:
        g_fanin(mk, p, 'b4')
    W = ('U', p)
    mk.alloc(W)
    mk.cx('b3', W, q0)
    add_lin(mk, c1_terms(p), W)
    mk.form('a1', 'b1', beta, lambda r: mk.mulx(r, 'b4', W))                                # W = U'(p)
    mk.form('a2', 'b2', beta, lambda r: mk.mulx(r, 'b3', 'b4'))                             # undo G' -> Z
    add_lin(mk, c2_terms(p), 'b4')
    mk.scale('b4', cinv(q1))
    mk.form('a3', 'b3', 0, lambda r: mk.mulx(r, 'b4', 'b3'))                                # undo Y' -> b3
    add_lin(mk, mu_terms(5, K2, p, q1), 'b3')
    mk.cx('a3', 'b3', q2)
    mk.scale('b3', cinv(cy))
    mk.cx('b3', 'a3', beta)
    with_mu5(mk, p, lambda: mk.form('a4', 'b4', 0, lambda r: mk.mulx(r, 'a5', 'b4')))       # undo Z -> b4
    mk.scale('b4', cinv(q2))
    mk.cx('b4', 'a4', beta)


def point0(mk, compute, part='all'):
    """U'(0) in place in b1 with Z(0), Y'(0), G'(0) in place in b4, b3, b2 (beta = 0: no roots). 4 products.
    Undo in two parts if asked: 'b1' (U'(0) -> b1, while G'(0), Y'(0) stay in b2, b3) then 'rest'."""
    k0s, k1s, k2s = csq(K0), csq(K1), csq(K2)
    if not compute and part == 'b1':
        mk.form('a1', 'b1', 0, lambda r: mk.mulx(r, 'b2', 'b1'))
        mk.cx('b3', 'b1', k0s)
        mk.cx('a1', 'b1', cmul(k1s, k2s))
        mk.scale('b1', cinv(cmul(cmul(k1s, k2s), K0)))
        return
    if not compute and part in ('rest', 'b2b3'):
        mk.form('a2', 'b2', 0, lambda r: mk.mulx(r, 'b3', 'b2'))
        mk.cx('b4', 'b2', k1s)
        mk.scale('b2', cinv(cmul(k1s, k2s)))
        mk.form('a3', 'b3', 0, lambda r: mk.mulx(r, 'b4', 'b3'))
        add_lin(mk, mu_terms(5, K2, 0, k1s), 'b3')
        mk.cx('a3', 'b3', k2s)
        mk.scale('b3', cinv(cmul(k2s, K1)))
        if part == 'b2b3':
            return
        with_mu5(mk, 0, lambda: mk.form('a4', 'b4', 0, lambda r: mk.mulx(r, 'a5', 'b4')))
        mk.scale('b4', cinv(k2s))
        return
    if not compute and part == 'b4z0':
        # fix round 2: P40 = a4^2 mu5(0) computed ONCE into the clean word z0; b4 ^= z0 undoes Z(0) -> k2^2 b4, and
        # z0 + k2^2 b4 is the coefficient word z0 = Z(0) the shared-Z erasure needs (one product for both)
        mk.alloc('z0')
        with_mu5(mk, 0, lambda: mk.form('a4', 'b4', 0, lambda r: mk.mulx(r, 'a5', 'z0')))
        mk.cx('z0', 'b4')
        mk.scale('b4', cinv(k2s))
        mk.cx('b4', 'z0', k2s)
        return
    if compute:
        mk.scale('b4', k2s)
        with_mu5(mk, 0, lambda: mk.form('a4', 'b4', 0, lambda r: mk.mulx(r, 'a5', 'b4')))
        mk.scale('b3', cmul(k2s, K1))
        mk.cx('a3', 'b3', k2s)
        add_lin(mk, mu_terms(5, K2, 0, k1s), 'b3')
        mk.form('a3', 'b3', 0, lambda r: mk.mulx(r, 'b4', 'b3'))
        mk.scale('b2', cmul(k1s, k2s))
        mk.cx('b4', 'b2', k1s)
        mk.form('a2', 'b2', 0, lambda r: mk.mulx(r, 'b3', 'b2'))
        mk.scale('b1', cmul(cmul(k1s, k2s), K0))
        mk.cx('a1', 'b1', cmul(k1s, k2s))
        mk.cx('b3', 'b1', k0s)
        mk.form('a1', 'b1', 0, lambda r: mk.mulx(r, 'b2', 'b1'))
    else:
        mk.form('a1', 'b1', 0, lambda r: mk.mulx(r, 'b2', 'b1'))
        mk.cx('b3', 'b1', k0s)
        mk.cx('a1', 'b1', cmul(k1s, k2s))
        mk.scale('b1', cinv(cmul(cmul(k1s, k2s), K0)))
        mk.form('a2', 'b2', 0, lambda r: mk.mulx(r, 'b3', 'b2'))
        mk.cx('b4', 'b2', k1s)
        mk.scale('b2', cinv(cmul(k1s, k2s)))
        mk.form('a3', 'b3', 0, lambda r: mk.mulx(r, 'b4', 'b3'))
        add_lin(mk, mu_terms(5, K2, 0, k1s), 'b3')
        mk.cx('a3', 'b3', k2s)
        mk.scale('b3', cinv(cmul(k2s, K1)))
        with_mu5(mk, 0, lambda: mk.form('a4', 'b4', 0, lambda r: mk.mulx(r, 'a5', 'b4')))
        mk.scale('b4', cinv(k2s))


def loop_spec(mk, n_points=33):
    """The loop's specification (its products counted): L_g ^= U'(w_g) X_g^2 / sc_g; zL ^= U'(0) zX / zc."""
    R = mk.R
    vals = {p: R[('U', p)] for p in UPTS[:5]}
    vals[0] = R['b1']                            # U'(0) in place
    for g, w in NODES:
        Uw = 0
        for t, c in lagr(UPTS, w, True).items():
            Uw ^= cescale(c, vals[t])
        X2 = R[('N', g, 0)] | (R[('N', g, 1)] << MB)
        v = cescale(cemul(Uw, X2), cinv(SC[g]))
        R[('L', g, 0)] ^= v & MASK
        R[('L', g, 1)] ^= v >> MB
    R[('L', 'z', 0)] ^= cmul(cmul(vals[0], R[('N', 'z', 0)]), cinv(ZC))
    mk.prods += n_points + 1


def compute(mk, fwd_g=None):
    """The block's compute on either call; fwd_g: the forward's G'(p) masks (return call). Returns its own G' masks."""
    gm = {}
    for i, p in enumerate(PTS[:4]):
        gm[p] = chain2(mk, p, i, None if fwd_g is None else fwd_g[p])
    chain_inplace_clean(mk, PTS[4])
    point0(mk, True)
    mk.cx('b2', ('L', 'z', 0), ZG) if False else None      # zero node's G'(0) comes through FANZ (GSET)
    return gm


def erase_return(mk, gm):
    """Return call after its loop: the five clean U' words measured; their masks and the return's own G'(p) masks
    paid with G' at 5 points and Y' at 4 (Lagrange): point 0's in-place values first, then point 0 undone, then
    3-word chains (Z, Y', G') at PTS[0..3]. 4 + 12 products."""
    um = {p: mk.hmr(('U', p)) for p in PTS[:5]}
    gpts = [0] + PTS[:4]
    ypts = [0] + PTS[:3]
    LG = {p: lagr(gpts, p) for p in PTS[:5]}
    LY = {p: lagr(ypts, p) for p in PTS[:5]}

    def pay(t, gword, yword):
        for p in PTS[:5]:
            beta = BETA_OF[p]
            if gword is not None:
                mk.form('a1', 'b1', beta, lambda r, _c=LG[p][t], _m=um[p]: mk.ppc(r, gword, _m, _c))
                if p in gm:
                    mk.czc(gword, gm[p], LG[p][t])
            if yword is not None and t in LY[p]:
                mk.czc(yword, um[p], cmul(qs(K0, p), LY[p][t]))
    point0(mk, False, 'b1')                      # b1 back first: h1(p) and c1(p) need it
    for p in PTS[:5]:                            # linear parts of U'(p)
        mk.cz_lin(c1_terms(p), um[p])
    pay(0, 'b2', 'b3')                           # G'(0), Y'(0) still in place in b2, b3
    if ZSHARE and Z0SHARE:
        point0(mk, False, 'b2b3')
        point0(mk, False, 'b4z0')
    else:
        point0(mk, False, 'rest')
    if ZSHARE:
        erase_chains_zshare(mk, pay, ypts, z0_ready=Z0SHARE)
        for m in list(um.values()) + list(gm.values()):
            mk.mfree(m)
        return
    for i, t in enumerate(PTS[:4]):
        beta = BETA_OF[t]
        Wz, Wy, Wg = ('E1', i), ('E2', i), ('E3', i)
        mk.alloc(Wz)
        add_lin(mk, c4_terms(t), Wz)
        with_mu5(mk, t, lambda: mk.form('a4', 'b4', beta, lambda r: mk.mulx(r, 'a5', Wz)))
        mk.alloc(Wy)
        add_lin(mk, c3_terms(t), Wy)
        mk.form('a3', 'b3', beta, lambda r: mk.mulx(r, Wz, Wy))
        mk.alloc(Wg)
        mk.cx(Wz, Wg, qs(K1, t))
        add_lin(mk, c2_terms(t), Wg)
        mk.form('a2', 'b2', beta, lambda r: mk.mulx(r, Wy, Wg))
        pay(t, Wg, Wy if t in ypts else None)
        m = mk.hmr(Wg)                           # G'(t) = Q1^2 Z + c2 + h2 Y'
        mk.czc(Wz, m, qs(K1, t))
        mk.cz_lin(c2_terms(t), m)
        mk.form('a2', 'b2', beta, lambda r: mk.pp(r, Wy, m))
        mk.mfree(m)
        m = mk.hmr(Wy)                           # Y'(t) = c3 + h3 Z
        mk.cz_lin(c3_terms(t), m)
        mk.form('a3', 'b3', beta, lambda r: mk.pp(r, Wz, m))
        mk.mfree(m)
        m = mk.hmr(Wz)                           # Z(t) = c4 + h4 mu5
        mk.cz_lin(c4_terms(t), m)
        with_mu5(mk, t, lambda: mk.form('a4', 'b4', beta, lambda r: mk.pp(r, 'a5', m)))
        mk.mfree(m)
    for m in list(um.values()) + list(gm.values()):
        mk.mfree(m)


ZSHARE = os.environ.get('DQI_INT_J6_ZSHARE', '') == '1'    # fix round 1 (B5-2): Z as coefficient words in the erasure
Z0SHARE = os.environ.get('DQI_INT_J6_Z0SHARE', '') == '1'  # fix round 2 (B5-3): one product for b4's undo and z0


def erase_chains_zshare(mk, pay, ypts, z0_ready=False):
    """The erasure's G' / Y' values at PTS[0..3] from Z held as coefficient words: z0 = k2^2 b4 + a4^2 m50,
    z2 = b4 + b4^2 b5, z1 = P41 + z0 + z2 + k2^2 b4 with P41 = (a4 + b4)^2 (m50 + b5) (3 products); per point
    Y'(t) and G'(t) (2 products) with Z(t) exposed in place on z0. 11 products against 12; 5 words.
    z0_ready: the word z0 already holds Z(0) (fix round 2)."""
    k2s = csq(K2)
    for n in (('z1', 'z2') if z0_ready else ('z0', 'z1', 'z2')):
        mk.alloc(n)
    if not z0_ready:
        mk.cx('b4', 'z0', k2s)                   # z0
        with_mu5(mk, 0, lambda: mk.form('a4', 'b4', 0, lambda r: mk.mulx(r, 'a5', 'z0')))
    mk.cx('b4', 'z2')                            # z2
    mk.form('b4', 'a4', 0, lambda r: mk.mulx(r, 'b5', 'z2'))
    mk.cx('b4', 'a4')                            # z1: (a4 + b4)^2 (m50 + b5) + z0 + z2 + k2^2 b4
    mk.cx('b5', 'a5', K2 ^ 1)
    mk.form('a4', 'b4', 0, lambda r: mk.mulx(r, 'a5', 'z1'))
    mk.cx('b5', 'a5', K2 ^ 1)
    mk.cx('b4', 'a4')
    mk.cx('z0', 'z1')
    mk.cx('z2', 'z1')
    mk.cx('b4', 'z1', k2s)

    def expose(t, undo):                         # z0 <-> Z(t) = z0 + t z1 + t^2 z2
        mk.cx('z1', 'z0', t)
        mk.cx('z2', 'z0', csq(t))
    for i, t in enumerate(PTS[:4]):
        beta = BETA_OF[t]
        Wy, Wg = ('E2', i), ('E3', i)
        expose(t, False)
        mk.alloc(Wy)
        add_lin(mk, c3_terms(t), Wy)
        mk.form('a3', 'b3', beta, lambda r: mk.mulx(r, 'z0', Wy))                            # Y'(t)
        mk.alloc(Wg)
        mk.cx('z0', Wg, qs(K1, t))
        add_lin(mk, c2_terms(t), Wg)
        mk.form('a2', 'b2', beta, lambda r: mk.mulx(r, Wy, Wg))                              # G'(t)
        pay(t, Wg, Wy if t in ypts else None)
        m = mk.hmr(Wg)
        mk.czc('z0', m, qs(K1, t))
        mk.cz_lin(c2_terms(t), m)
        mk.form('a2', 'b2', beta, lambda r: mk.pp(r, Wy, m))
        mk.mfree(m)
        m = mk.hmr(Wy)
        mk.cz_lin(c3_terms(t), m)
        mk.form('a3', 'b3', beta, lambda r: mk.pp(r, 'z0', m))
        mk.mfree(m)
        expose(t, True)
    m = mk.hmr('z1')                             # z1 = b4 + a4^2 b5 + b4^2 m50
    mk.cz_lin([('b4', 1)], m)
    mk.form('a4', 'b4', 0, lambda r: mk.pp(r, 'b5', m))
    with_mu5(mk, 0, lambda: mk.form('b4', 'a4', 0, lambda r: mk.pp(r, 'a5', m)))
    mk.mfree(m)
    m = mk.hmr('z2')                             # z2 = b4 + b4^2 b5
    mk.cz_lin([('b4', 1)], m)
    mk.form('b4', 'a4', 0, lambda r: mk.pp(r, 'b5', m))
    mk.mfree(m)
    m = mk.hmr('z0')                             # z0 = k2^2 b4 + a4^2 m50
    mk.cz_lin([('b4', k2s)], m)
    with_mu5(mk, 0, lambda: mk.form('a4', 'b4', 0, lambda r: mk.pp(r, 'a5', m)))
    mk.mfree(m)


def block_forward(mk, check=None):
    """The forward call of the J = 6 block: compute, loop, five U' words measured (masks to the return), point 0
    undone. Returns the masks the return pays."""
    gm = compute(mk)
    if check:
        check('forward compute')
    mk.loop()
    um = {p: mk.hmr(('U', p)) for p in PTS[:5]}
    point0(mk, False)
    return {'g': gm, 'u': um}


def block_return(mk, fwd, check=None):
    """The return call: compute (paying the forward's G' masks on the same words), the forward's U' masks paid on
    the recomputed words, loop, erasure."""
    gm = compute(mk, fwd['g'])
    if check:
        check('return compute')
    for p in PTS[:5]:
        mk.cz(('U', p), fwd['u'][p])
        mk.mfree(fwd['u'][p])
    mk.loop()
    erase_return(mk, gm)


def polys_c(h, w):
    """U', G' at w (E) from the closed forms (ch's field)."""
    w2 = cesq(w)
    Q0s, Q1s, Q2s = w2 ^ w ^ csq(K0), w2 ^ w ^ csq(K1), w2 ^ w ^ csq(K2)

    def hh(i):
        return csq(h['a%d' % i]) ^ cescale(w, csq(h['b%d' % i]))

    def mu(i, k):
        return h['a%d' % i] ^ cmul(k, h['b%d' % i]) ^ cescale(w, h['b%d' % i])
    Z = cescale(Q2s, h['b4']) ^ cemul(hh(4), mu(5, K2))
    Y = cemul(Q2s, mu(3, K1)) ^ cemul(Q1s, mu(5, K2)) ^ cemul(hh(3), Z)
    Gp = cemul(Q1s, cescale(Q2s, h['b2']) ^ Z) ^ cemul(hh(2), Y)
    U = cemul(cemul(Q1s, Q2s), mu(1, K0)) ^ cemul(Q0s, Y) ^ cemul(hh(1), Gp)
    return U, Gp


def random_regs(rng, special):
    R = {k: rng.randrange(NF) for k in HEADS}
    if special:
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


def one_decode(R0, mseed):
    mk = Mock(R0, random.Random(mseed))
    h = {k: R0[k] for k in HEADS}
    st = {}

    def check(where):
        for p in PTS[:5]:
            assert mk.R[('U', p)] == polys_c(h, p)[0], ('U after ' + where, p)
        assert mk.R['b1'] == polys_c(h, 0)[0], 'U(0) in b1 after ' + where
    fwd = block_forward(mk, check)
    st['fwd_products'], st['fwd_peak'] = mk.prods, mk.peak
    for g, w in NODES:
        U, Gp = polys_c(h, w)
        X2 = R0[('N', g, 0)] | (R0[('N', g, 1)] << MB)
        want = (R0[('L', g, 0)] | (R0[('L', g, 1)] << MB)) ^ cescale(cemul(U, X2), cinv(SC[g])) ^ cescale(Gp, GC[g])
        assert (mk.R[('L', g, 0)] | (mk.R[('L', g, 1)] << MB)) == want, ('L after forward', g)
    U0, G0 = polys_c(h, 0)
    assert mk.R[('L', 'z', 0)] == R0[('L', 'z', 0)] ^ cmul(cmul(U0, R0[('N', 'z', 0)]), cinv(ZC)) ^ cmul(G0, ZG)
    for k in HEADS:
        assert mk.R[k] == R0[k], ('head after forward', k)
    assert not mk.tmp, mk.tmp
    mk.peak, p0 = 0, mk.prods
    block_return(mk, fwd, check)
    st['ret_products'], st['ret_peak'] = mk.prods - p0, mk.peak
    assert not mk.tmp, mk.tmp
    for k in R0:
        assert mk.R[k] == R0[k], ('not restored', k)
    return mk.phase, st


def cmd_sched(inputs, masks, rng):
    bad = 0
    st = None
    for _ in range(masks):
        mseed = rng.randrange(1 << 30)
        ph = set()
        for i in range(inputs):
            phase, st = one_decode(random_regs(rng, i % 4 == 3), mseed)
            ph.add(phase)
        bad += len(ph) > 1
    print(f'SCHED B5-1: {masks} mask samples x {inputs} inputs; values pass (U words both calls, L block identity, '
          f'restoration); data-dependent phase in {bad} of {masks}; forward {st["fwd_products"]} products peak '
          f'{st["fwd_peak"]} words; return {st["ret_products"]} products peak {st["ret_peak"]} words')
    print('SCHED PASS' if not bad else 'SCHED FAIL')
    return not bad


# ===================================================================================================================
# Part 3 (mode `sched6`): the count of fix 1 -- the same block if the block had a SIXTH clean word (12 more qubits).
# Five two-word chains, the sixth clean word by the in-place-transient chain, no in-place U' word; the return pays its
# six measured U' words with G' at 5 points and Y' at 4 from Z coefficient words. A counting model of a layout that
# does not fit today's 62 wires (it peaks at 6 clean words), not a build.
# ===================================================================================================================
PT6 = csq(11)
BETA_OF[PT6] = 11
UPTS6 = PTS[:5] + [PT6]


def loop_spec6(mk, n_points=33):
    R = mk.R
    vals = {p: R[('U', p)] for p in UPTS6}
    for g, w in NODES:
        Uw = 0
        for t, c in lagr(UPTS6, w, True).items():
            Uw ^= cescale(c, vals[t])
        X2 = R[('N', g, 0)] | (R[('N', g, 1)] << MB)
        v = cescale(cemul(Uw, X2), cinv(SC[g]))
        R[('L', g, 0)] ^= v & MASK
        R[('L', g, 1)] ^= v >> MB
    U0 = 0
    for t, c in lagr(UPTS6, 0).items():          # U'(0) exposed from the six words (free) for the zero node
        U0 ^= cmul(c, vals[t])
    R[('L', 'z', 0)] ^= cmul(cmul(U0, R[('N', 'z', 0)]), cinv(ZC))
    mk.prods += n_points + 1


def compute6(mk, fwd_g=None):
    gm = {}
    for i, p in enumerate(PTS[:5]):
        gm[p] = chain2(mk, p, i, None if fwd_g is None else fwd_g[p])
    chain_inplace_clean(mk, PT6)
    return gm


def erase_return6(mk, gm):
    um = {p: mk.hmr(('U', p)) for p in UPTS6}
    gpts, ypts = PTS[:5], PTS[:4]
    LG = {p: lagr(gpts, p) for p in UPTS6}
    LY = {p: lagr(ypts, p) for p in UPTS6}

    def pay(t, gword, yword):
        for p in UPTS6:
            mk.form('a1', 'b1', BETA_OF[p], lambda r, _c=LG[p][t], _m=um[p]: mk.ppc(r, gword, _m, _c))
            if p in gm:
                mk.czc(gword, gm[p], LG[p][t])
            if yword is not None:
                mk.czc(yword, um[p], cmul(qs(K0, p), LY[p][t]))
    for p in UPTS6:
        mk.cz_lin(c1_terms(p), um[p])
    k2s = csq(K2)
    for n in ('z0', 'z1', 'z2'):
        mk.alloc(n)
    mk.cx('b4', 'z0', k2s)
    with_mu5(mk, 0, lambda: mk.form('a4', 'b4', 0, lambda r: mk.mulx(r, 'a5', 'z0')))
    mk.cx('b4', 'z2')
    mk.form('b4', 'a4', 0, lambda r: mk.mulx(r, 'b5', 'z2'))
    mk.cx('b4', 'a4')
    mk.cx('b5', 'a5', K2 ^ 1)
    mk.form('a4', 'b4', 0, lambda r: mk.mulx(r, 'a5', 'z1'))
    mk.cx('b5', 'a5', K2 ^ 1)
    mk.cx('b4', 'a4')
    mk.cx('z0', 'z1')
    mk.cx('z2', 'z1')
    mk.cx('b4', 'z1', k2s)
    for i, t in enumerate(gpts):
        beta = BETA_OF[t]
        Wy, Wg = ('E2', i), ('E3', i)
        mk.cx('z1', 'z0', t)
        mk.cx('z2', 'z0', csq(t))
        mk.alloc(Wy)
        add_lin(mk, c3_terms(t), Wy)
        mk.form('a3', 'b3', beta, lambda r: mk.mulx(r, 'z0', Wy))
        mk.alloc(Wg)
        mk.cx('z0', Wg, qs(K1, t))
        add_lin(mk, c2_terms(t), Wg)
        mk.form('a2', 'b2', beta, lambda r: mk.mulx(r, Wy, Wg))
        pay(t, Wg, Wy if t in ypts else None)
        m = mk.hmr(Wg)
        mk.czc('z0', m, qs(K1, t))
        mk.cz_lin(c2_terms(t), m)
        mk.form('a2', 'b2', beta, lambda r: mk.pp(r, Wy, m))
        m = mk.hmr(Wy)
        mk.cz_lin(c3_terms(t), m)
        mk.form('a3', 'b3', beta, lambda r: mk.pp(r, 'z0', m))
        mk.cx('z2', 'z0', csq(t))
        mk.cx('z1', 'z0', t)
    m = mk.hmr('z1')
    mk.cz_lin([('b4', 1)], m)
    mk.form('a4', 'b4', 0, lambda r: mk.pp(r, 'b5', m))
    with_mu5(mk, 0, lambda: mk.form('b4', 'a4', 0, lambda r: mk.pp(r, 'a5', m)))
    m = mk.hmr('z2')
    mk.cz_lin([('b4', 1)], m)
    mk.form('b4', 'a4', 0, lambda r: mk.pp(r, 'b5', m))
    m = mk.hmr('z0')
    mk.cz_lin([('b4', k2s)], m)
    with_mu5(mk, 0, lambda: mk.form('a4', 'b4', 0, lambda r: mk.pp(r, 'a5', m)))


def one_decode6(R0, mseed):
    mk = Mock(R0, random.Random(mseed))
    h = {k: R0[k] for k in HEADS}
    fwd_g = compute6(mk)
    for p in UPTS6:
        assert mk.R[('U', p)] == polys_c(h, p)[0], ('U6 after compute', p)
    loop_spec6(mk)
    fwd_u = {p: mk.hmr(('U', p)) for p in UPTS6}
    st = dict(fwd_products=mk.prods, fwd_peak=mk.peak)
    for g, w in NODES:
        U, Gp = polys_c(h, w)
        X2 = R0[('N', g, 0)] | (R0[('N', g, 1)] << MB)
        want = (R0[('L', g, 0)] | (R0[('L', g, 1)] << MB)) ^ cescale(cemul(U, X2), cinv(SC[g])) ^ cescale(Gp, GC[g])
        assert (mk.R[('L', g, 0)] | (mk.R[('L', g, 1)] << MB)) == want, ('L6 after forward', g)
    for k in HEADS:
        assert mk.R[k] == R0[k], ('head6 after forward', k)
    mk.peak, p0 = 0, mk.prods
    gm = compute6(mk, fwd_g)
    for p in UPTS6:
        mk.cz(('U', p), fwd_u[p])
    loop_spec6(mk)
    erase_return6(mk, gm)
    st.update(ret_products=mk.prods - p0, ret_peak=mk.peak)
    assert not mk.tmp, mk.tmp
    for k in R0:
        assert mk.R[k] == R0[k], ('not restored (6)', k)
    return mk.phase, st


def cmd_sched6(inputs, masks, rng):
    bad, st = 0, None
    for _ in range(masks):
        mseed = rng.randrange(1 << 30)
        ph = set()
        for i in range(inputs):
            phase, st = one_decode6(random_regs(rng, i % 4 == 3), mseed)
            ph.add(phase)
        bad += len(ph) > 1
    print(f'SCHED6 (six clean words, fix 1 count): {masks} x {inputs}; data-dependent phase in {bad}; forward '
          f'{st["fwd_products"]} products peak {st["fwd_peak"]} words; return {st["ret_products"]} peak {st["ret_peak"]}')
    print('SCHED6 PASS' if not bad else 'SCHED6 FAIL')
    return not bad


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'forms'
    if mode == 'sched6':
        sys.exit(0 if cmd_sched6(int(sys.argv[2]) if len(sys.argv) > 2 else 8,
                                 int(sys.argv[4]) if len(sys.argv) > 4 else 4,
                                 random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)) else 1)
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    if mode == 'forms':
        sys.exit(0 if cmd_forms(trials, rng) else 1)
    if mode == 'sched':
        masks = int(sys.argv[4]) if len(sys.argv) > 4 else 8
        sys.exit(0 if cmd_sched(trials, masks, rng) else 1)
