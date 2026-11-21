"""dqi_int_pfdp.py -- B5's copy of D1's tools/dqi_pub_pfd.py (sha 495877d9, LF) + R37-DQI-EXT's patch P.

Used in place of dqi_pub_pfd.py through an overlay snapshot (the file is copied in under that name). Unset
DQI_INT_PFDP: exactly D1's module. DQI_INT_PFDP=obs: D1's module, and the free stack before / after the block's group
loop is compared and printed (no gate change). DQI_INT_PFDP=1: patch P (R37-DQI-EXT.md, "the block's group loop as
one polynomial product"; schedule = cmd_blockloop of tools/r37_dqi_ext_chain.py): all 14 N_g hold X_g^2 at once, and
for each of 31 points t (F elements 0..29 and infinity) U(t) is exposed in place on one U word, XX(t) = sum_g
Tr(X_g^2 lam_g(t)) in place on one N word, one 42-AND product goes into one L word inside a free spread that adds
(1 / sc_g) W_g,t p_t to every L_g; 31 x 42 = 1,302 against 14 x 99 = 1,386 per block call. After the loop the free
stack is set back to its order before the loop (no gates), as the original loop leaves it when it is stack-neutral
(DQI_INT_PFDP=obs measures that).

D1's docstring follows.
dqi_pub_pfd.py -- builder D1: the prefix block (J = 5) with its forward erasure deferred to the return.

The build line's block (tools/dqi_build_patch_prefix5.py, read-only) is called identically on both trips (it is an
involution on the later groups' L words). Each call computes U at the points alpha = 1, 4, infinity into three clean
words and at alpha = 0 in place, applies the 14 group products, and then erases the three clean U words by
recomputing Y and G at each point (2 + 2 + 3 = 7 twelve-bit products, 294 Toffolis).
Deferral: on the FORWARD call the three clean U words are X-measured instead (masks m_p kept as classical bits).
The head words that U is made of (heads 1..4) hold the same values when the RETURN call runs, so the return call,
right after it has recomputed U(alpha_p) into the same kind of word, pays the forward's phase (-1)^(m_p . U_p) with
a classically controlled Z on that word (free), and then erases its own U words as before. Point 0 (in place in the
head words) is still undone on both trips. The skipped recomputation allocated and freed two 12-qubit transients
per point; kickmix reuses freed qubits last-in first-out per qubit, so each such cycle reverses that chunk of the
free stack: the forward replays the same allocate / free sequence without gates, which keeps every later placement
identical (without the replay the final routing failed on q21, build p2pfd).
Measured: p2pfdb = p1a1f + pfd, 360,702.8 executed (quick), -293.9 against p1a1f; static -294; Q 902.
Implementation: the source of tools/dqi_build_patch_prefix5.py is executed in a private module with two exact text
replacements (asserted to match once each); nothing in the original file changes.
Install (from tools/dqi_pub_entry.py): DQI_PUB_P=...,pfd  (after B5's decode has installed prefix5).
"""
import os
import sys
import types

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, 'dqi_build_patch_prefix5.py')

_PAY_OLD = ("            point0(True)\n"
            "            R1, R2, R3 = U[1], U[4], U[None]\n")
_PAY_NEW = ("            if not _p0:\n"
            "                point0(True)\n"
            "                _B1 = b1\n"
            "            _pfd_pay(b, U, b1)\n"
            "            if _PFD_ZJOIN:\n"
            "                _zX, _zL, _za0 = words[68], words[69], words[0]\n"
            "                b.cx(_za0, _zX); square12(b, _zX)\n"
            "                _zc = mul(sq(k0), sq(k1))\n"
            "                scale(b, _zL, _zc); ar.mul_xor0(b, _B1, _zX, _zL); scale(b, _zL, inv(_zc))\n"
            "                square12(b, _zX, True); b.cx(_za0, _zX)\n"
            "                if not _p0:\n"
            "                    cx(b, b4, _zL, inv(sq(k1)))\n"
            "            R1, R2, R3 = U[1], U[4], U[None]\n")
# DQI_PUB_PFD_P0=1: on the FORWARD call only, the point alpha = 0 runs out of place and first: Y(0), G(0), U(0) into
# three clean words (3 products), Y(0) and G(0) erased by measurement at once (free: their repair is bilinear in live
# words), U(0) kept for the group loop and X-measured at the end like the other U words; the return, whose point 0
# stays in place (b_1 holds U(0) there), pays the mask on b_1. Forward point 0: 6 products -> 3 (-126 [P count]).
# The forward block peaks at 60 clean qubits instead of 48 (12 wires created early, at high-water 888 -> 900); at the
# block's exit the free stack is set to the one the unchanged line has there (recorded by a scratch decode at
# install), with the early wires at the bottom, smallest index on top: kickmix pops before it creates and creates in
# index order, so each early wire is taken exactly where the unchanged line creates it (wires 888..901 are created
# later in the same forward extract [M, S3's q902s trace]).
_START_OLD = ("            U = {}\n"
              "            for p in TEMP_POINTS:\n")
_START_NEW = ("            U = {}\n"
              "            _p0 = _pfd_p0(b)\n"
              "            if _p0:\n"
              "                _k0s, _k1s = sq(k0), sq(k1)\n"
              "                _Y0 = b.alloc_qubits(12)\n"
              "                b.cx(a3, _Y0); cx(b, b3, _Y0, k1); sqmul(a3, b4, _Y0)\n"
              "                _G0 = b.alloc_qubits(12)\n"
              "                cx(b, b2, _G0, _k1s); cx(b, b4, _G0, _k1s); sqmul(a2, _Y0, _G0)\n"
              "                _B1 = b.alloc_qubits(12)\n"
              "                cx(b, b1, _B1, mul(_k1s, k0)); cx(b, a1, _B1, _k1s); cx(b, _Y0, _B1, _k0s)\n"
              "                sqmul(a1, _G0, _B1)\n"
              "                if _PFD_ZJOIN:\n"
              "                    cx(b, _G0, words[69], inv(_k1s))\n"
              "                cx(b, b2, _G0, _k1s); cx(b, b4, _G0, _k1s)\n"
              "                square12(b, a2); ar.erase_product(b, _Y0, a2, _G0); square12(b, a2, True)\n"
              "                b.free(_G0)\n"
              "                b.cx(a3, _Y0); cx(b, b3, _Y0, k1)\n"
              "                square12(b, a3); ar.erase_product(b, b4, a3, _Y0); square12(b, a3, True)\n"
              "                b.free(_Y0)\n"
              "            for p in TEMP_POINTS:\n")
_LOOP_OLD = [("                cx(b, b1, R1, lam[0][0]); cx(b, b1, R2, lam[0][1])\n",
              "                cx(b, _B1, R1, lam[0][0]); cx(b, _B1, R2, lam[0][1])\n"),
             ("                cx(b, b1, R2, lam[0][1]); cx(b, b1, R1, lam[0][0])\n",
              "                cx(b, _B1, R2, lam[0][1]); cx(b, _B1, R1, lam[0][0])\n")]
P0 = os.environ.get('DQI_PUB_PFD_P0', '0') == '1'
P0_STATE = dict(active=P0, recording=False, order=None, fixes=0, misses=0, extras=0)
# DQI_PUB_ZJOIN=1 (needs the zt module's matching switch): the node z = 0 joins the block. At point 0 the block
# holds U(0) in b_1 and G(0) in b_4, and z = 0 is a later node with w = 0, kappa_g = 0, so the zero node's L word
# gets U(0) (x_0 + a_0)^2 / (k0^2 k1^2) + G(0) / k1^2 (one 12-bit product) and its chain runs heads 1..4 without
# their L products (4 x 42 per trip). Identity: tools/dqi_pub_zcheck.py, 500 / 500 [M classical].
ZJOIN = os.environ.get('DQI_PUB_ZJOIN', '0') == '1'
_END_OLD = ("            point0(False)\n"
            "            for p in reversed(TEMP_POINTS):\n"
            "                temp_point(p, U[p], False, p == TEMP_POINTS[-1])\n"
            "                b.free(U[p])\n")
_END_NEW = ("            if not _p0:\n"
            "                point0(False)\n"
            "            _pfd_end(b, U, temp_point, p0u=_B1 if _p0 else None)\n")
_TP_OLD = "                temp_point(p, U[p], True, p == TEMP_POINTS[-1])\n"
_TP_NEW = "                _pfd_tp(b, temp_point, p, U[p], True, p == TEMP_POINTS[-1])\n"

# DQI_PUB_PFD_RINF=1: on the RETURN call only, the point at infinity uses a clean G word instead of G in place in
# the b_4 register (saves the in-place undo product in its compute and in its erasure: -2 products). The return
# block runs when the circuit's high-water mark is already 901 (S3's q902s trace), so its 60-qubit peak creates no
# wire; the extra 12-qubit G cycle is undone on the free stack by an exact two-chunk reorder (_fix24).
RINF = os.environ.get('DQI_PUB_PFD_RINF', '0') == '1'

STATS = dict(forward=0, reverse=0, rinf=0)

# ---- B5: patch P (R37-DQI-EXT) -------------------------------------------------------------------------------------
PFDP = os.environ.get('DQI_INT_PFDP', '')
PFDP_STATS = dict(mode=PFDP, loops=0, products=0, stack_same=[], restored=0, restore_miss=0)
_LOOPP_OLD = ("            m = combine([R1, R2])\n"
              "            for gi in hom:\n")
_LOOPP_NEW = ("            m = combine([R1, R2])\n"
              "            _pfdp_pre(b)\n"
              "            if _PFDP_ON:\n"
              "                _pfdp_loop(b, words, hom, {0: _B1, 1: R1, 4: R2, None: R3}, a0, h0, combine, square24,\n"
              "                           k0, k1, cx, scale, ar, mul, inv, K, groups)\n"
              "            for gi in (() if _PFDP_ON else hom):\n")
_LOOPE_OLD = "                mix(R1, R2, lam[1], lam[4], True)\n"
_LOOPE_NEW = ("                mix(R1, R2, lam[1], lam[4], True)\n"
              "            _pfdp_post(b)\n")
_PRE = {}


def _stack_order(b):
    """The free stack as a list of qubit indices, top first (pops everything and pushes it back: no gates)."""
    got = _spare(b)
    order = [_qidx(q) for q in got]
    for q in reversed(got):
        b.free(q)
    return order


def _stack_set(b, order):
    """Set the free stack to `order` (top first) if it holds exactly those qubits; False (unchanged) otherwise."""
    got = _spare(b)
    have = {_qidx(q): q for q in got}
    if set(have) != set(order) or len(order) != len(got):
        for q in reversed(got):
            b.free(q)
        return False
    for i in reversed(order):
        b.free(have[i])
    return True


def _pfdp_pre(b):
    if PFDP:
        _PRE[id(b)] = _stack_order(b)


def _pfdp_post(b):
    if not PFDP:
        return
    pre = _PRE.pop(id(b))
    now = _stack_order(b)
    PFDP_STATS['stack_same'].append(now == pre)
    if PFDP == '1' and now != pre:
        if _stack_set(b, pre):
            PFDP_STATS['restored'] += 1
        else:
            PFDP_STATS['restore_miss'] += 1
    print(f'dqi_int_pfdp: group loop exit ({PFDP}): free stack {"unchanged" if now == pre else "CHANGED"} '
          f'({len(pre)} -> {len(now)} wires)', file=sys.stderr, flush=True)


def _pfdp_loop(b, words, hom, U, a0, h0, combine, square24, k0, k1, cx, scale, ar, mul, inv, K, groups):
    """Patch P: L_g ^= U(w_g) X_g^2 / sc_g for every g in hom as 31 rank-one updates (see the module docstring)."""
    def emul(u, v):
        t = mul(u[1], v[1])
        return (mul(u[0], v[0]) ^ mul(K, t), mul(u[0], v[1]) ^ mul(u[1], v[0]) ^ t)

    def einv(u):
        a, c = u
        n = mul(a, a) ^ mul(a, c) ^ mul(K, mul(c, c))
        ni = inv(n)
        return (mul(a ^ c, ni), mul(c, ni))

    def eadd(u, v):
        return (u[0] ^ v[0], u[1] ^ v[1])

    def escale(u, s):
        return (mul(u[0], s), mul(u[1], s))

    def sq(v):
        return mul(v, v)

    hom = list(hom)
    W = {gi: (groups[gi]['c'] ^ groups[gi]['kap'], 1) for gi in hom}           # w_g = gamma_g^2 as (lo, hi)
    allnodes = [W[gi] for gi in hom] + [(W[gi][0] ^ W[gi][1], W[gi][1]) for gi in hom]   # + conjugates
    assert len(set(allnodes)) == len(allnodes)
    pts = list(range(len(allnodes) + 2))                                         # deg U XX = 30: 30 points + inf
    fin = (0, 1, 4)

    def lam_x(node, t):                         # prod_{nu != node} (t - nu) / (node - nu); t None: 1 / prod(...)
        num, den = (1, 0), (1, 0)
        for nu in allnodes:
            if nu != node:
                if t is not None:
                    num = emul(num, eadd((t, 0), nu))
                den = emul(den, eadd(node, nu))
        return emul(num, einv(den))

    def lam_u(t):                               # U(t) = sum_p lam_u[p] U(p), t in F; U(None) = leading coefficient
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

    def spread_c(t, w):                         # coefficient of p_t in (U XX)(w)
        if t is None:
            out = (1, 0)
            for s in pts:
                out = emul(out, eadd(w, (s, 0)))
            return out
        num, den = (1, 0), 1
        for s in pts:
            if s != t:
                num = emul(num, eadd(w, (s, 0)))
                den = mul(den, t ^ s)
        return escale(num, inv(den))

    sc = {gi: mul(sq(k0 ^ groups[gi]['kap']), sq(k1 ^ groups[gi]['kap'])) for gi in hom}
    base = {gi: groups[gi]['word'] for gi in hom}
    Nw = {gi: (words[base[gi]], words[base[gi] + 1]) for gi in hom}
    Lw = {gi: (words[base[gi] + 2], words[base[gi] + 3]) for gi in hom}

    def n_to_x2(gi, undo):                      # N_g <-> N_g + H_0(gamma_g), squared (all linear)
        N = combine(words[base[gi]:base[gi] + 2])
        ct = groups[gi]['c']
        if undo:
            square24(b, N, True)
        cx(b, h0, a0, ct); b.cx(combine([a0, h0]), N); cx(b, h0, a0, ct)
        if not undo:
            square24(b, N)

    for gi in hom:
        n_to_x2(gi, False)
    for t in pts + [None]:
        lu = None                               # expose x = U(t) on one U word
        if t is None:
            ux = U[None]
        elif t in fin:
            ux = U[t]
        else:
            lu = lam_u(t)
            ux = U[1]
            scale(b, ux, lu[1])
            for p in (0, 4, None):
                cx(b, U[p], ux, lu[p])
        cf = []                                 # expose y = XX(t) on one N word: HI(lam N) = lam1 n0 + (lam0+lam1) n1
        for gi in hom:
            lam = lam_x(W[gi], t)
            cf.append((Nw[gi][0], lam[1]))
            cf.append((Nw[gi][1], lam[0] ^ lam[1]))
        yi = next(i for i, (_, c) in enumerate(cf) if c)
        yw, yc = cf[yi]
        scale(b, yw, yc)
        for i, (w_, c) in enumerate(cf):
            if i != yi and c:
                cx(b, w_, yw, c)
        cc = []                                 # spread: L_g ^= (1 / sc_g) W_g,t p
        for gi in hom:
            c = escale(spread_c(t, W[gi]), inv(sc[gi]))
            cc.append((Lw[gi][0], c[0]))
            cc.append((Lw[gi][1], c[1]))
        ti = next(i for i, (_, c) in enumerate(cc) if c)
        tw, tc = cc[ti]
        scale(b, tw, inv(tc))
        for i, (w_, c) in enumerate(cc):
            if i != ti and c:
                cx(b, tw, w_, c)
        ar.mul_xor0(b, ux, yw, tw)              # the one 42-AND product of this point
        PFDP_STATS['products'] += 1
        for i, (w_, c) in enumerate(cc):
            if i != ti and c:
                cx(b, tw, w_, c)
        scale(b, tw, tc)
        for i, (w_, c) in enumerate(cf):        # un-expose y, then x
            if i != yi and c:
                cx(b, w_, yw, c)
        scale(b, yw, inv(yc))
        if lu is not None:
            for p in (0, 4, None):
                cx(b, U[p], ux, lu[p])
            scale(b, ux, inv(lu[1]))
    for gi in hom:
        n_to_x2(gi, True)
    PFDP_STATS['loops'] += 1


# ---- B5: patch T (R37-DQI-EXT, form ii), DQI_INT_PFDT=1, RETURN call only -----------------------------------------
# The return takes the forward's P0 layout (U(0) in a clean word), keeps G in b_4 at its last temp point (RINF off
# there: width), and after the group loop X-measures its four U words and pays their masks with Y, G held as
# COEFFICIENT words: Pa = b4 a3^2, Pb = b4 b3^2 (y0 = a3 + k1 b3 + Pa, y1 = b3 + Pb), V0 = a2^2 y0, V2 = b2^2 y1,
# V1 = (a2^2 + b2^2)(y0 + y1) + V0 + V2 (g0 = k1^2 (b2 + b4) + V0, g1 = b2 + b4 + V1, g2 = b2 + b4 + V2): 5 products;
# U(p) = q1 (a1 + (k0 + p) b1) + q0 (y0 + p y1) + (a1 + beta b1)^2 (g0 + p g1 + p^2 g2), U(inf) = b1 + y1 + b1^2 g2;
# each mask phase = CZs on an exposed linear form + one free phase_product. Return block: 18 products against 21.
PFDT = os.environ.get('DQI_INT_PFDT', '') == '1'
PFDT_STATS = dict(erasures=0, products=0)
_PAYT_OLD = "            _pfd_pay(b, U, b1)\n"
_PAYT_NEW = "            _pfd_pay(b, U, _B1)\n"
_ENDT_OLD = "            _pfd_end(b, U, temp_point, p0u=_B1 if _p0 else None)\n"
_ENDT_NEW = ("            _pfdt_end(b, U, temp_point, _B1 if _p0 else None, a1, b1, a2, b2, a3, b3, b4, k0, k1, cx, scale,\n"
             "                      square12, ar, mul, inv, form)\n")


def _pfdt_end(b, U, temp_point, u0, a1, b1, a2, b2, a3, b3, b4, k0, k1, cx, scale, square12, ar, mul, inv, form):
    if _DIR[id(b)] or u0 is None:
        _pfd_end(b, U, temp_point, p0u=u0)
        return
    _pfdt_erase(b, U, u0, a1, b1, a2, b2, a3, b3, b4, k0, k1, cx, scale, square12, ar, mul, inv, form)


def _pfdt_erase(b, U, u0, a1, b1, a2, b2, a3, b3, b4, k0, k1, cx, scale, square12, ar, mul, inv, form):
    def sq(v):
        return mul(v, v)

    def cxc(src, dst, c):                       # dst ^= c src, skipped for c = 0
        if c:
            cx(b, src, dst, c)

    def sqmul(src, fac, dst):                   # dst ^= fac src^2: one 42-AND product
        square12(b, src)
        ar.mul_xor0(b, fac, src, dst)
        square12(b, src, True)
        PFDT_STATS['products'] += 1

    beta_of = {0: 0, 1: 1, 4: 2}
    Uw = {1: U[1], 4: U[4], None: U[None], 0: u0}
    masks = {}
    for p in (None, 4, 1, 0):                   # the four U words out by X-measurement (masks kept)
        m = b.alloc_bits(12)
        b.hmr(Uw[p], m)
        masks[p] = m
        b.free(Uw[p])
    y0 = _PFDC_PA.pop(id(b), None)              # patch C keeps Pa = y0 live through the loop
    if y0 is None:
        y0 = b.alloc_qubits(12)
        sqmul(a3, b4, y0)                       # Pa
        b.cx(a3, y0)
        cxc(b3, y0, k1)                         # y0 = a3 + k1 b3 + Pa
    y1 = b.alloc_qubits(12)
    sqmul(b3, b4, y1)                           # Pb
    b.cx(b3, y1)                                # y1 = b3 + Pb
    g0 = b.alloc_qubits(12)
    sqmul(a2, y0, g0)                           # V0 = a2^2 y0
    g2 = b.alloc_qubits(12)
    sqmul(b2, y1, g2)                           # V2 = b2^2 y1
    g1 = b.alloc_qubits(12)
    b.cx(b2, a2)
    b.cx(y1, y0)
    sqmul(a2, y0, g1)                           # (a2 + b2)^2 (y0 + y1)
    b.cx(y1, y0)
    b.cx(b2, a2)
    b.cx(g0, g1)
    b.cx(g2, g1)                                # V1
    k1s = sq(k1)
    cxc(b2, g0, k1s)
    cxc(b4, g0, k1s)                            # g0 = k1^2 (b2 + b4) + V0
    for g in (g1, g2):
        b.cx(b2, g)
        b.cx(b4, g)                             # g1, g2 = b2 + b4 + V1, V2
    for p in (0, 1, 4, None):                   # pay (-1)^(m_p . U(p))
        m = masks[p]
        if p is None:
            b.cx(b1, y1)
            b.cz(y1, m)
            b.cx(b1, y1)
            square12(b, b1)
            ar.phase_product(b, b1, g2, m)
            square12(b, b1, True)
        else:
            q0, q1 = sq(p) ^ p ^ sq(k0), sq(p) ^ p ^ sq(k1)
            assert q0 and q1
            scale(b, y0, q0)
            cxc(a1, y0, q1)
            cxc(b1, y0, mul(q1, k0 ^ p))
            cxc(y1, y0, mul(q0, p))
            b.cz(y0, m)
            cxc(y1, y0, mul(q0, p))
            cxc(b1, y0, mul(q1, k0 ^ p))
            cxc(a1, y0, q1)
            scale(b, y0, inv(q0))
            cxc(g1, g0, p)
            cxc(g2, g0, sq(p))
            form(a1, b1, beta_of[p], lambda r: ar.phase_product(b, r, g0, m))
            cxc(g2, g0, sq(p))
            cxc(g1, g0, p)
        b.free(m)
    for g in (g1, g2):                          # coefficient words out: linear parts off, then free erasures
        b.cx(b4, g)
        b.cx(b2, g)
    cxc(b4, g0, k1s)
    cxc(b2, g0, k1s)
    b.cx(g0, g1)
    b.cx(g2, g1)                                # g1 = (a2 + b2)^2 (y0 + y1)
    b.cx(b2, a2)
    b.cx(y1, y0)
    square12(b, a2)
    ar.erase_product(b, y0, a2, g1)
    square12(b, a2, True)
    b.cx(y1, y0)
    b.cx(b2, a2)
    square12(b, a2)
    ar.erase_product(b, y0, a2, g0)             # V0
    square12(b, a2, True)
    square12(b, b2)
    ar.erase_product(b, y1, b2, g2)             # V2
    square12(b, b2, True)
    b.cx(a3, y0)
    cxc(b3, y0, k1)
    square12(b, a3)
    ar.erase_product(b, b4, a3, y0)             # Pa
    square12(b, a3, True)
    b.cx(b3, y1)
    square12(b, b3)
    ar.erase_product(b, b4, b3, y1)             # Pb
    square12(b, b3, True)
    for w in (g1, g2, g0, y1, y0):
        b.free(w)
    PFDT_STATS['erasures'] += 1


# ---- B5: patch C (R37-DQI-EXT-reprice "first build step"), DQI_INT_PFDC=1 (needs DQI_INT_PFDT=1), both calls -------
# The block's compute phase as C's 11-product schedule (op for op as c_compute in tools/r37_dqi_ext_rp_phase.py): the
# four U words from coefficient words (Pa, Pb, V0, V2, V1; U~(1), U~(4), U~(0), U~(inf); Pa, Pb again), the G fan-ins
# into the L words and the zero node's L word done from g0, g1, g2; every measured coefficient word paid by phases
# bilinear in live words. Forward: step 8 then 9, Pa, Pb erased. Return: step 9 then 8 with Pa KEPT as y0 for T's
# erasure (5 -> 4 products). Off while P0 is inactive (D1's capture decode records the unchanged line's exit stack).
PFDC = os.environ.get('DQI_INT_PFDC', '') == '1'
# B5 (Nov 15): DQI_INT_J6=1 replaces the whole block (both calls) by the J = 6 block of tools/dqi_int_j6blk.py and
# makes head 5 N-only in the chain; the zero node's chain needs the matching zt copy (tools/dqi_int_zt6.py).
J6 = os.environ.get('DQI_INT_J6', '') == '1'
PFDC_STATS = dict(computes=0, products=0)
_PFDC_PA = {}
_C1_OLD = ("            _p0 = _pfd_p0(b)\n"
           "            if _p0:\n")
_C1_NEW = ("            _p0 = _pfd_p0(b)\n"
           "            _c = _pfdc_on(b)\n"
           "            if _c:\n"
           "                _B1 = _pfdc_compute(b, U, words, hom, a1, b1, a2, b2, a3, b3, b4, k0, k1, cx, scale,\n"
           "                                    square12, ar, mul, inv, K, groups, form)\n"
           "                _p0 = True\n"
           "            elif _p0:\n")
_C2_OLD = ("            for p in TEMP_POINTS:\n"
           "                U[p] = b.alloc_qubits(12)\n")
_C2_NEW = ("            for p in (() if _c else TEMP_POINTS):\n"
           "                U[p] = b.alloc_qubits(12)\n")


def _pfdc_on(b):
    return PFDC and P0_STATE['active']


def _pfdc_compute(b, U, words, hom, a1, b1, a2, b2, a3, b3, b4, k0, k1, cx, scale, square12, ar, mul, inv, K,
                  groups, form):
    keep_pa = not _DIR[id(b)]                    # return call: Pa stays live as y0 for patch T

    def sq(v):
        return mul(v, v)

    def cxc(s, d, c=1):                          # d ^= c s (skipped for c = 0)
        if c == 1:
            b.cx(s, d)
        elif c:
            cx(b, s, d, c)

    def fsq(w, body):
        square12(b, w)
        body(w)
        square12(b, w, True)

    def mulx(x, y, z):
        ar.mul_xor0(b, x, y, z)
        PFDC_STATS['products'] += 1

    def hmr(w):
        m = b.alloc_bits(12)
        b.hmr(w, m)
        b.free(w)
        return m

    def erase(x, y, z):
        ar.erase_product(b, x, y, z)
        b.free(z)

    def pp(x, y, m):
        ar.phase_product(b, x, y, m)

    def emul(u, v):
        t = mul(u[1], v[1])
        return (mul(u[0], v[0]) ^ mul(K, t), mul(u[0], v[1]) ^ mul(u[1], v[0]) ^ t)

    k0s, k1s = sq(k0), sq(k1)
    beta_of = {0: 0, 1: 1, 4: 2}
    # steps 1, 2: coefficient words (coef_words)
    Pa = b.alloc_qubits(12)
    fsq(a3, lambda r: mulx(b4, r, Pa))
    Pb = b.alloc_qubits(12)
    fsq(b3, lambda r: mulx(b4, r, Pb))
    cxc(a3, Pa)
    cxc(b3, Pa, k1)                              # y0
    cxc(b3, Pb)                                  # y1
    g0 = b.alloc_qubits(12)
    fsq(a2, lambda r: mulx(Pa, r, g0))
    g2 = b.alloc_qubits(12)
    fsq(b2, lambda r: mulx(Pb, r, g2))
    g1 = b.alloc_qubits(12)
    cxc(b2, a2)
    cxc(Pb, Pa)
    fsq(a2, lambda r: mulx(Pa, r, g1))
    cxc(Pb, Pa)
    cxc(b2, a2)
    cxc(g0, g1)
    cxc(g2, g1)                                  # V1
    cxc(b2, g0, k1s)
    cxc(b4, g0, k1s)
    for gw in (g1, g2):
        cxc(b2, gw)
        cxc(b4, gw)
    # G fan-ins (g_fanin): L_g ^= (g0 + w_g g1 + w_g^2 g2) / (k1 + kap_g)^2; zL ^= g0 / k1^2
    for gi in hom:
        g = groups[gi]
        c0 = inv(sq(k1 ^ g['kap']))
        w = (g['c'] ^ g['kap'], 1)
        c1 = (mul(w[0], c0), mul(w[1], c0))
        w2 = emul(w, w)
        c2 = (mul(w2[0], c0), mul(w2[1], c0))
        L0, L1 = words[g['word'] + 2], words[g['word'] + 3]
        cxc(g0, L0, c0)
        cxc(g1, L0, c1[0])
        cxc(g1, L1, c1[1])
        cxc(g2, L0, c2[0])
        cxc(g2, L1, c2[1])
    cxc(g0, words[69], inv(k1s))
    # step 3: erase Pb, Pa (free)
    cxc(b3, Pb)
    fsq(b3, lambda r: erase(b4, r, Pb))
    cxc(b3, Pa, k1)
    cxc(a3, Pa)
    fsq(a3, lambda r: erase(b4, r, Pa))
    # step 4: U~(1), U~(4)
    for p in (1, 4):
        U[p] = b.alloc_qubits(12)
        q1 = sq(p) ^ p ^ k1s
        cxc(a1, U[p], q1)
        cxc(b1, U[p], mul(q1, k0 ^ p))
        scale(b, g1, p)
        cxc(g0, g1)
        cxc(g2, g1, sq(p))                       # G(p) on g1
        form(a1, b1, beta_of[p], lambda r, _u=U[p]: mulx(g1, r, _u))
        cxc(g2, g1, sq(p))
        cxc(g0, g1)
        scale(b, g1, inv(p))
    m1 = hmr(g1)                                 # step 5: g1 = b2 + b4 + a2^2 y1 + b2^2 y0
    b.cz(b2, m1)
    b.cz(b4, m1)
    U0 = b.alloc_qubits(12)                      # step 6
    cxc(a1, U0, k1s)
    cxc(b1, U0, mul(k1s, k0))
    fsq(a1, lambda r: mulx(g0, r, U0))
    m0 = hmr(g0)                                 # g0 = k1^2 (b2 + b4) + a2^2 y0
    cxc(b4, b2)
    scale(b, b2, k1s)
    b.cz(b2, m0)
    scale(b, b2, inv(k1s))
    cxc(b4, b2)
    U[None] = b.alloc_qubits(12)                 # step 7
    cxc(b1, U[None])
    fsq(b1, lambda r: mulx(g2, r, U[None]))
    m2 = hmr(g2)                                 # g2 = b2 + b4 + b2^2 y1
    b.cz(b2, m2)
    b.cz(b4, m2)
    Uw = {0: U0, 1: U[1], 4: U[4], None: U[None]}

    def step8(keep):                             # Pa again: y0 live
        pa = b.alloc_qubits(12)
        fsq(a3, lambda r: mulx(b4, r, pa))
        cxc(a3, pa)
        cxc(b3, pa, k1)
        for p in (0, 1, 4):
            cxc(pa, Uw[p], sq(p) ^ p ^ k0s)
        fsq(a2, lambda r: pp(r, pa, m0))         # m0's a2^2 y0
        fsq(b2, lambda r: pp(r, pa, m1))         # m1's b2^2 y0
        if keep:
            return pa
        cxc(b3, pa, k1)
        cxc(a3, pa)
        fsq(a3, lambda r: erase(b4, r, pa))
        return None

    def step9():                                 # Pb again: y1 live
        pb = b.alloc_qubits(12)
        fsq(b3, lambda r: mulx(b4, r, pb))
        cxc(b3, pb)
        for p in (1, 4):
            cxc(pb, Uw[p], mul(sq(p) ^ p ^ k0s, p))
        cxc(pb, Uw[None])
        fsq(b2, lambda r: pp(r, pb, m2))         # m2's b2^2 y1
        fsq(a2, lambda r: pp(r, pb, m1))         # m1's a2^2 y1
        cxc(b3, pb)
        fsq(b3, lambda r: erase(b4, r, pb))

    if keep_pa:
        step9()
        _PFDC_PA[id(b)] = step8(True)
    else:
        step8(False)
        step9()
    for m in (m0, m1, m2):
        b.free(m)
    PFDC_STATS['computes'] += 1
    return U0


def _pfdp_report():
    s = dict(PFDP_STATS)
    s['stack_same'] = ''.join('1' if x else '0' for x in s['stack_same'])
    print(f'dqi_int_pfdp stats {s} T {PFDT_STATS} C {PFDC_STATS}', file=sys.stderr, flush=True)


if PFDP or PFDT or PFDC:
    import atexit
    atexit.register(_pfdp_report)
# ---- end of B5's patch P part --------------------------------------------------------------------------------------
_CALLS = {}       # id(builder) -> number of block calls so far (even: the next call is a forward call)
_MASKS = {}       # id(builder) -> stack of {point: classical mask bits} from forward calls
_DIR = {}         # id(builder) -> direction of the call in progress


def _pfd_start(b):
    """Block entry: the first call of a builder is the forward trip, the second the return."""
    k = id(b)
    n = _CALLS.get(k, 0)
    _CALLS[k] = n + 1
    _DIR[k] = n % 2 == 0


def _fix24(b):
    """After temp_point with a clean G where the original kept G in place: the original popped Y (x1..x12) and
    pushed it back, leaving the next chunk y untouched; the clean version also popped and pushed G (y1..y12), which
    reverses y. Pop both chunks and push them back in the original order (no gates)."""
    reg = b.alloc_qubits(24)                  # pops x12..x1 then y12..y1
    for i in range(12, 24):                   # push y12, ..., y1  (y back in its original order)
        b.free(reg[i])
    for i in reversed(range(12)):             # push x1, ..., x12  (x12 on top, as the original left it)
        b.free(reg[i])


def _pfd_tp(b, temp_point, p, U, compute, inplace):
    if inplace and RINF and not _DIR[id(b)] and not (PFDT and P0_STATE['active']):   # B5 T: G stays in b_4
        temp_point(p, U, compute, False)
        _fix24(b)
        STATS['rinf'] += 1
    else:
        temp_point(p, U, compute, inplace)


def _pfd_p0(b):
    """True on a forward call when the forward point 0 runs out of place (DQI_PUB_PFD_P0=1); with B5's patch T
    (DQI_INT_PFDT=1) on the return call as well."""
    return P0_STATE['active'] and (_DIR[id(b)] or PFDT)


def _pfd_pay(b, U, b1):
    """Start of the group loop: on a return call pay the forward's measured U words (key 0: U(0), which the return
    holds in place in b_1)."""
    k = id(b)
    if _DIR[k]:
        STATS['forward'] += 1
        return
    STATS['reverse'] += 1
    masks = _MASKS[k].pop()
    for p, m in masks.items():
        b.cz(b1 if p == 0 else U[p], m)   # (-1)^(m_p . U(alpha_p)): the forward's phase, U recomputed here
        b.free(m)


def _qidx(q):
    return int(str(q)[1:])


def _spare(b):
    n = b.finish_circuit().num_qubits - b.num_allocated_qubits
    return list(b.alloc_qubits(n)) if n > 0 else []


def _p0_snapshot(b):
    got = _spare(b)
    order = [_qidx(q) for q in got]                  # stack top first
    for q in reversed(got):
        b.free(q)
    return order


def _p0_restore(b):
    """Forward block exit under P0: the unchanged line's free stack, the wires P0 created early below it."""
    order = P0_STATE['order']
    got = _spare(b)
    have = {_qidx(q): q for q in got}
    extra = sorted(set(have) - set(order or ()))
    if order is None or not set(order) <= set(have) or (extra and order and extra[0] <= max(order)):
        P0_STATE['misses'] += 1
        for q in reversed(got):
            b.free(q)
        print(f'dqi_pub_pfd: P0 stack MISS (recorded {None if order is None else len(order)}, spare {len(have)})',
              file=sys.stderr, flush=True)
        return
    for i in reversed(extra):                        # bottom: the early wires, smallest index on top of them
        b.free(have[i])
    for i in reversed(order):
        b.free(have[i])
    P0_STATE['fixes'] += 1
    P0_STATE['extras'] = len(extra)
    print(f'dqi_pub_pfd: P0 stack restored ({len(order)} wires, {len(extra)} early wires at the bottom)',
          file=sys.stderr, flush=True)


def _pfd_end(b, U, temp_point, points=None, p0u=None):
    """End of the block: forward call measures the three clean U words; return call erases them as before."""
    k = id(b)
    pts = points if points is not None else _MOD.TEMP_POINTS
    if _DIR[k]:
        masks = {}
        for p in reversed(pts):
            # keep the free stack identical to the original erasure: temp_point(p, U, False, inplace) allocates
            # Y (and G when not in place) and frees them last-in first-out; kickmix reuses freed qubits in that
            # order, so the same allocate / free sequence is replayed here without gates (fix round 1, 12:01Z:
            # without it the final routing failed on q21, as a3 v5 did).
            inplace = p == pts[-1]
            Y = b.alloc_qubits(12)
            if not inplace:
                Gr = b.alloc_qubits(12)
                b.free(Gr)
            b.free(Y)
            m = b.alloc_bits(12)
            b.hmr(U[p], m)
            masks[p] = m
            b.free(U[p])
        if p0u is not None:                          # P0: U(0) out of place, measured like the others
            m = b.alloc_bits(12)
            b.hmr(p0u, m)
            masks[0] = m
            b.free(p0u)
        _MASKS.setdefault(k, []).append(masks)
    else:
        for p in reversed(pts):
            _pfd_tp(b, temp_point, p, U[p], False, p == pts[-1])
            b.free(U[p])


def _load_module():
    src = open(_SRC).read()
    pairs = [(_PAY_OLD, _PAY_NEW), (_END_OLD, _END_NEW), (_TP_OLD, _TP_NEW), (_START_OLD, _START_NEW)] + _LOOP_OLD
    if PFDP:                                    # B5: hooks around the group loop (off: D1's source unchanged)
        pairs += [(_LOOPP_OLD, _LOOPP_NEW), (_LOOPE_OLD, _LOOPE_NEW)]
    if PFDT:                                    # B5: patch T (applied after D1's own replacements)
        pairs += [(_PAYT_OLD, _PAYT_NEW), (_ENDT_OLD, _ENDT_NEW)]
    if PFDC:                                    # B5: patch C (needs T: the return pays U(0) on _B1, Pa kept)
        assert PFDT, 'DQI_INT_PFDC=1 needs DQI_INT_PFDT=1'
        pairs += [(_C1_OLD, _C1_NEW), (_C2_OLD, _C2_NEW)]
    for old, new in pairs:
        assert src.count(old) == 1, 'prefix5 source changed: deferral not applied'
        src = src.replace(old, new)
    mod = types.ModuleType('dqi_pub_pfd_prefix5')
    mod.__file__ = _SRC
    mod._pfd_pay = _pfd_pay
    mod._pfd_end = _pfd_end
    mod._pfd_tp = _pfd_tp
    mod._pfd_p0 = _pfd_p0
    mod._PFD_ZJOIN = ZJOIN
    mod._PFDP_ON = PFDP == '1'
    mod._pfdp_pre = _pfdp_pre
    mod._pfdp_post = _pfdp_post
    mod._pfdp_loop = _pfdp_loop
    mod._pfdt_end = _pfdt_end
    mod._pfdc_on = _pfdc_on
    mod._pfdc_compute = _pfdc_compute
    exec(compile(src, _SRC + ':pfd', 'exec'), mod.__dict__)
    return mod


_MOD = _load_module()


def install(km, ns):
    """Replace the prefix block of every live external_extraction clone with the deferring copy."""
    walk = sys.modules['walk_bound']
    fns = walk.functions(ns.decode)
    need = ('ar', 'cx', 'scale', 'scalepair', 'combine', 'groups', 'ks', 'K', 'mul', 'inv', 'square12',
            'external_square24', 'mul24_xor99')
    mul_head = next(f for f in fns if f.__name__ == 'external_mul_head' and all(k in f.__globals__ for k in need))
    inner = _MOD.make_block(mul_head.__globals__)
    n6 = 0
    if J6:                                       # B5: the prefix block at J = 6 (tools/dqi_int_j6blk.py, schedule B5-1)
        import dqi_int_j6blk
        j6 = dqi_int_j6blk.make_j6_block(mul_head.__globals__, lambda bb: _DIR[id(bb)], _MASKS)
        inner_j5 = inner

        def inner(b, words, homogeneous):        # the P0 capture decode records the unchanged line's exit stack
            if P0_STATE['recording']:
                return inner_j5(b, words, homogeneous)
            return j6(b, words, homogeneous)
        for f in fns:                            # head 5 becomes an N-only chain step
            if f.__name__ == 'chain_step' and '_dqb_nonly' in f.__globals__:
                f.__globals__['_dqb_nonly'] = frozenset({1, 2, 3, 4, 5})
                n6 += 1
        assert n6, 'no prefix chain_step with _dqb_nonly found'
        assert P0 and ZJOIN, 'DQI_INT_J6=1 needs DQI_PUB_PFD_P0=1 and DQI_PUB_ZJOIN=1'
        print(f'dqi_int_pfdp: J = 6 block installed, head 5 N-only in {n6} chain_step clones', file=sys.stderr,
              flush=True)

    def block(b, words, homogeneous):
        _pfd_start(b)
        out = inner(b, words, homogeneous)
        if _DIR[id(b)]:
            if P0_STATE['recording']:
                P0_STATE['order'] = _p0_snapshot(b)
            elif P0_STATE['active']:
                _p0_restore(b)
        return out

    n = 0
    for f in fns:
        if f.__name__ == 'external_extraction' and '_dqb_block' in f.__code__.co_names:
            f.__globals__['_dqb_block'] = block
            n += 1
    assert n, 'no external_extraction with the prefix block found (is prefix5 installed?)'
    if P0:
        _p0_capture(km or sys.modules['circ'], ns.decode)
    return dict(extraction_clones=n, rinf=RINF, zjoin=ZJOIN, p0=P0, pfdp=PFDP, pfdt=PFDT, pfdc=PFDC, j6=n6,
                p0_recorded=None if P0_STATE['order'] is None else len(P0_STATE['order']))


def _p0_capture(km, decode):
    """Scratch decode of the line with the forward point 0 in place (as it is without P0), recording the free stack
    at the forward block's exit (as tools/dqi_pub_zt.py records its own two points)."""
    a3 = sys.modules.get('dqi_pub_a3')
    st = dict(a3.STACK) if a3 is not None else None   # a3's replay counters must not be spent on this scratch
    CB = km.CircuitBuilder
    tagged = CB.mark
    fv = getattr(getattr(tagged, '__code__', None), 'co_freevars', ())
    if 'done' in fv and 'original' in fv:                # r6tag's once-per-process tag: not on the scratch build
        CB.mark = dict(zip(fv, tagged.__closure__))['original'].cell_contents
    P0_STATE.update(active=False, recording=True)
    b0 = None
    try:
        b0 = CB()
        S0 = [b0.alloc_qubits(12) for _ in range(70)]
        E0 = [b0.alloc_bits(12) for _ in range(4095)]
        decode(b0, S0, E0)
    finally:
        CB.mark = tagged
        if a3 is not None:
            a3.STACK.update(st)
        P0_STATE.update(active=True, recording=False)
        if b0 is not None:
            for d in (_CALLS, _DIR, _MASKS):
                d.pop(id(b0), None)
    del b0, S0, E0
    import gc
    gc.collect()
    print(f'dqi_pub_pfd: P0 capture: forward block exit stack of {len(P0_STATE["order"] or ())} wires',
          file=sys.stderr, flush=True)
