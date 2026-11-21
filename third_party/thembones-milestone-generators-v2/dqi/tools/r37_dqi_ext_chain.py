#!/usr/bin/env python
"""r37_dqi_ext_chain.py -- designer R37-DQI-EXT: exact word-level model of the extraction pair's chain and junction.

Own file. The field, the nodes, the jets and the affine reference extraction are COPIED from the builders' classical
model tools/dqi_research_model.py (read-only there); everything below "R37" is this seat's.

Modes (one core, seconds unless noted):
  ident [shots] [seed]   identities of the anatomy (design file section 1), on weight-35 instances:
                         (a) every chain step is  N ^= H D, Lam ^= l D  with l_even = sqrt(b) (an F scalar) and
                             l_odd = sqrt(a + k_j b) + sqrt(b) gamma_g  (node enters only through gamma_g);
                         (b) the pair form  L_2 = q^2 L + q^2 b D^2 + (H' + q b') D_1^2;
                         (c) T1's rider flow through groups 3..16 (carried gauge 1/n) gives the affine heads;
                         (d) the return identities used by ledger rows 12, 13.
  count                  products per direction by kind (to compare with K1's census of best_q3b).
  census [shots] [seed] [groups]   rank census at every 24-bit product site of the rider flow (minutes; box):
                         T1 operand ranks over GF(2); T2 number of bilinear forms sum c_ij x_i y_j of the two
                         operands that are GF(2)-linear in the live state before the product (an "alias": a
                         coefficient of the product that is already live, R36's device). 0 expected if none.
"""
import random
import sys

M = 12
MOD = 0x1053
NF = 1 << M
MASK = NF - 1
EXP = []
LOG = [-1] * NF
_x = 1
for _i in range(NF - 1):
    LOG[_x] = _i
    EXP.append(_x)
    _x <<= 1
    if _x & NF:
        _x ^= MOD
assert _x == 1
EXP = EXP * 2
K = 2048


def mul(a, b):
    return EXP[LOG[a] + LOG[b]] if a and b else 0


def inv(a):
    assert a
    return EXP[NF - 1 - LOG[a]]


def emul(u, v):
    a, b, c, d = u & MASK, u >> M, v & MASK, v >> M
    bd = mul(b, d)
    return (mul(a, c) ^ mul(K, bd)) | ((mul(a, d) ^ mul(b, c) ^ bd) << M)


def esq(u):
    return emul(u, u)


def einv(u):
    a, b = u & MASK, u >> M
    n = mul(a, a) ^ mul(a, b) ^ mul(K, mul(b, b))
    assert n
    ni = inv(n)
    return mul(a ^ b, ni) | (mul(b, ni) << M)


def escale(u, s):
    return mul(u & MASK, s) | (mul(u >> M, s) << M)


G = 17
CS = [2 * g for g in range(G)]
KAP = [mul(c, c) ^ c ^ K for c in CS]


def head_value(a, b, c):
    return (a ^ mul(b, c)) | (b << M)


def coords(u, c):
    b = u >> M
    return (u & MASK) ^ mul(b, c), b


def jets(roots, weights):
    out = []
    for g in range(G):
        u = v = 0
        for beta, w in zip(roots, weights):
            t = einv((CS[g] | NF) ^ beta)
            u ^= escale(t, w)
            v ^= escale(esq(t), w)
        out.append((u, v))
    return out


class Pole(Exception):
    pass


def reference_extract(J):
    heads = []
    for g in range(G):
        u, v = J[g]
        for i in range(2 * g):
            a, b = heads[i]
            q = KAP[i // 2] ^ KAP[g]
            u ^= head_value(a, b, CS[g])
            v ^= b
            if not u:
                raise Pole(g, i)
            ui = einv(u)
            v = emul(v, esq(ui))
            u = ui
            v = u ^ escale(v, q)
            u = escale(u, q)
        a, b = coords(u, CS[g])
        heads.append((a, b))
        v ^= b
        if not v:
            raise Pole(g, 'local')
        a2, b2 = coords(einv(v), CS[g])
        heads.append((a2, b2))
    return heads


def chain_params(g):
    out = []
    c = 0
    for i in range(2 * g):
        q = KAP[i // 2] ^ KAP[g]
        out.append((q, c))
        c ^= inv(q)
    assert c == 0
    return out


def kit_step(H, b, N, D, L, q, c):
    """The kit's step (external_chain_step): L <- q (L + (c H + b) D^2), (N, D) <- (q D, N + H D)."""
    L = escale(L ^ emul(escale(H, c) ^ b, esq(D)), q)
    return escale(D, q), N ^ emul(H, D), L


# ---------------------------------------------------------------- R37 ------------------------------------------
def fsqrt(a):
    return EXP[(LOG[a] + (NF - 1) * (LOG[a] & 1)) // 2] if a else 0


SK = fsqrt(K)


def esqrt(u):
    """sqrt in E: sqrt(y) = y + sqrt(K)."""
    a, b = fsqrt(u & MASK), fsqrt(u >> M)
    return (a ^ mul(b, SK)) | (b << M)


def ell(i, g, a, b):
    """The second linear polynomial of step i at node g, as an E value: l_even = sqrt(b); l_odd = A + B gamma_g."""
    if i % 2 == 0:
        return fsqrt(b)
    A, B = fsqrt(a ^ mul(KAP[i // 2], b)), fsqrt(b)
    return (A ^ mul(B, CS[g])) | (B << M)


def lam_step(i, g, a, b, N, D, Lam, q):
    """One step in the square-root form: N ^= H D, Lam ^= l D, then the free scalings."""
    H = head_value(a, b, CS[g])
    P = emul(H, D)
    R = emul(ell(i, g, a, b), D)
    sq = fsqrt(q)
    if i % 2 == 0:
        Lam = escale(Lam ^ R, sq)
    else:
        Lam = escale(Lam, sq) ^ R
    return escale(D, q), N ^ P, Lam


def instance(rng, weight=35):
    while True:
        roots = rng.sample(range(1, NF), weight)
        weights = [rng.randrange(1, NF) for _ in roots]
        J = jets(roots, weights)
        try:
            return J, reference_extract(J)
        except Pole:
            continue


HOM = range(3, 17)


def rider_flow(J, heads, sites=None):
    """T1 with the in-place junction (carried gauge 1/n). Returns the heads it produces (as E values per group).
    sites: optional callback(kind, g, i, x, y, state) called before every 24-bit product (x, y operands)."""
    out = {}
    rider = None
    for g in HOM:
        x, L0 = J[g]
        a0, b0 = heads[0]
        q0 = KAP[0] ^ KAP[g]
        X = x ^ head_value(a0, b0, CS[g])
        Lam0 = esqrt(L0 ^ b0)

        def state(extra):
            if sites is None:
                return None
            st = []
            for j in range(G):
                if j in out:
                    st.extend(out[j])
                elif j < 3:
                    st.extend((head_value(*heads[2 * j], CS[j]), head_value(*heads[2 * j + 1], CS[j])))
                elif j != g:
                    st.extend(J[j])
            st.extend(extra)
            return st

        if rider is None:
            N, D, Lam = q0, X, escale(Lam0, fsqrt(q0))
        else:
            if sites:
                sites('j_in_X', g, 0, rider, X, state([rider, X, Lam0]))
            D = emul(rider, X)
            if sites:
                sites('j_in_L', g, 0, rider, Lam0, state([rider, D, Lam0]))
            Lam = escale(emul(rider, Lam0), fsqrt(q0))
            N = escale(rider, q0)
        for i in range(1, 2 * g):
            a, b = heads[i]
            q = KAP[i // 2] ^ KAP[g]
            if sites:
                H = head_value(a, b, CS[g])
                sites('N_only' if i <= 4 else ('odd_N' if i % 2 else 'even_N'), g, i, H, D, state([N, D, Lam]))
                if i > 4:
                    sites('odd_L' if i % 2 else 'even_L', g, i, ell(i, g, a, b), D, state([N, D, Lam]))
            N, D, Lam = lam_step(i, g, a, b, N, D, Lam, q)
        Gi = einv(D)
        if sites:
            sites('j_out_N', g, 0, N, Gi, state([N, Gi, Lam]))
        u = emul(N, Gi)
        if sites:
            sites('j_out_L', g, 0, Lam, Gi, state([u, Gi, Lam]))
        sv = emul(Lam, Gi)
        v = esq(sv)
        w = einv(v ^ (u >> M))
        out[g] = (u, w)
        rider = Gi
    return out


def cmd_ident(shots, seed):
    rng = random.Random(seed)
    bad = dict(a=0, b=0, c=0, d=0)
    steps = 0
    for _ in range(shots):
        J, heads = instance(rng)
        # (a), (b): the kit's chain against the square-root form, every group, started from the raw jet
        for g in range(1, G):
            N, D, L = J[g][0], 1, J[g][1]
            N2, D2, Lam = N, D, esqrt(L)
            params = chain_params(g)
            for i, (q, c) in enumerate(params):
                a, b = heads[i]
                H = head_value(a, b, CS[g])
                if i % 2 == 0:
                    Hn = head_value(*heads[i + 1], CS[g])
                    bn = heads[i + 1][1]
                    D1 = N ^ emul(H, D)
                    q2 = mul(q, q)
                    pair = escale(L, q2) ^ escale(esq(D), mul(q2, b)) ^ emul(Hn ^ mul(q, bn), esq(D1))
                N, D, L = kit_step(H, b, N, D, L, q, c)
                N2, D2, Lam = lam_step(i, g, a, b, N2, D2, Lam, q)
                steps += 1
                if (N, D) != (N2, D2) or esq(Lam) != L:
                    bad['a'] += 1
                if i % 2 == 1 and L != pair:
                    bad['b'] += 1
        # (c) rider flow
        out = rider_flow(J, heads)
        for g in HOM:
            u, w = out[g]
            if u != head_value(*heads[2 * g], CS[g]) or w != head_value(*heads[2 * g + 1], CS[g]):
                bad['c'] += 1
        # (d) return identities: with rider G = 1/n, y = v + d, w = 1/y, Z = Lam' + sqrt(d) n:
        #     n = Z sqrt(w),  Lam' = 1/(G sqrt(w)) + sqrt(d) / G,  N_prev G = 1 + H' u / q  (last step)
        for g in HOM:
            N, D, Lam = J[g][0], 1, esqrt(J[g][1])
            prev = None
            for i in range(2 * g):
                a, b = heads[i]
                q = KAP[i // 2] ^ KAP[g]
                prev = (N, D, a, b, q)
                N, D, Lam = lam_step(i, g, a, b, N, D, Lam, q)
            Gi = einv(D)
            u = emul(N, Gi)
            d = u >> M
            w = einv(esq(emul(Lam, Gi)) ^ d)
            Z = Lam ^ escale(D, fsqrt(d))
            Hp = head_value(prev[2], prev[3], CS[g])
            ok = emul(Z, esqrt(w)) == D
            ok &= Lam == einv(emul(Gi, esqrt(w))) ^ escale(D, fsqrt(d))
            ok &= emul(prev[0], Gi) == 1 ^ escale(emul(Hp, u), inv(prev[4]))
            if not ok:
                bad['d'] += 1
    print(f'ident: shots {shots}, chain steps checked {steps}, mismatches {bad}')
    print('IDENT PASS' if not any(bad.values()) else 'IDENT FAIL')


def cmd_count():
    kinds = {}
    def site(kind, g, i, x, y, st):
        kinds[kind] = kinds.get(kind, 0) + 1
    rng = random.Random(1)
    J, heads = instance(rng)
    rider_flow(J, heads, site)
    price = dict(N_only=99, odd_N=99, odd_L=87, even_N=99, even_L=79, j_in_X=99, j_in_L=99, j_out_N=99, j_out_L=99)
    tot = 0
    for k in sorted(kinds):
        tot += kinds[k] * price[k]
        print(f'{k:8s} {kinds[k]:4d} x {price[k]:3d} = {kinds[k] * price[k]:6d}')
    print(f'per direction {tot} (fused odd = odd_N + odd_L = 186, even = 178); '
          f'steps: odd {kinds["odd_N"]}, even {kinds["even_N"]}, N-only {kinds["N_only"]}')


def rank(rows):
    piv = {}
    r = 0
    for v in rows:
        while v:
            h = v.bit_length() - 1
            p = piv.get(h)
            if p is None:
                piv[h] = v
                r += 1
                break
            v ^= p
    return r


def cmd_census(shots, seed, groups):
    rng = random.Random(seed)
    data = {}
    def site(kind, g, i, x, y, st):
        if g not in groups:
            return
        s = 0
        for wd in st:
            s = (s << 24) | wd
        data.setdefault((g, kind, i), []).append((x, y, (s << 1) | 1))   # the constant 1: affine functions
    for _ in range(shots):
        J, heads = instance(rng)
        rider_flow(J, heads, site)
    tot = dict(sites=0, low_operand=0, alias=0)
    for key in sorted(data):
        rec = data[key]
        rx = rank([x for x, _, _ in rec])
        ry = rank([y for _, y, _ in rec])
        rows_s, rows_all = [], []
        for x, y, s in rec:
            mono = 0
            for bi in range(24):
                if (x >> bi) & 1:
                    mono |= y << (24 * bi)
            rows_s.append(s)
            rows_all.append((s << 576) | mono)
        rs = rank(rows_s)
        ra = rank(rows_all)
        alias = rs + 576 - ra
        tot['sites'] += 1
        tot['low_operand'] += (rx < 24) + (ry < 24)
        tot['alias'] += alias
        flag = '' if (rx == 24 and ry == 24 and alias == 0) else '   <----'
        print(f'g {key[0]:2d} {key[1]:8s} step {key[2]:2d}: operand ranks {rx} {ry}; state rank {rs}; '
              f'aliases {alias}{flag}', flush=True)
    print(f'census: shots {shots}, {tot}')


def econj(u):
    return ((u & MASK) ^ (u >> M)) | (u & ~MASK)


def cmd_blockpoly(trials, seed, zero):
    """Patch P: the prefix block's group loop, L0_g ^= U(w_g) X_g^2 for the 14 homogeneous groups (14 products of
    99, + one 42 for the zero node when it joins), as ONE polynomial product U(w) * XX(w) over F, XX the CRT image
    of the X_g^2 (w_g = gamma_g^2): deg U = 3, deg XX = 27 (28 with the zero node) -> 31 (32) evaluation products
    of 42. Checks, on random words, that the 31 (32) rank-one updates give the same 14 (15) values."""
    rng = random.Random(seed)
    nodes = [(g, esq(CS[g] | NF)) for g in HOM]                     # w_g in E
    allnodes = [w for _, w in nodes] + [econj(w) for _, w in nodes] + ([0] if zero else [])
    assert len(set(allnodes)) == len(allnodes)
    npts = len(allnodes) + 2                                         # finite points (deg U XX = npts); + infinity
    pts = list(range(npts))                                          # t = 0 .. npts-1 in F (0, 1, 4 among them)
    fin = [0, 1, 4]

    def lag_x(node, t):                                              # prod_{nu != node} (t - nu) / (node - nu), in E
        num = den = 1
        for nu in allnodes:
            if nu != node:
                num = emul(num, t ^ nu)
                den = emul(den, node ^ nu)
        return emul(num, einv(den)), einv(den)

    LX = {(g, t): lag_x(w, t)[0] for g, w in nodes for t in pts}
    MU = {g: lag_x(w, 0)[1] for g, w in nodes}
    if zero:
        LZ = {t: lag_x(0, t)[0] for t in pts}
        MUZ = lag_x(0, 0)[1]

    def u_at(U, t):                                                  # the cubic through U_0, U_1, U_4 with lead U_inf
        out = 0
        allp = 1
        for p in fin:
            num, den = 1, 1
            for r in fin:
                if r != p:
                    num = emul(num, t ^ r)
                    den = mul(den, p ^ r)
            out ^= emul(escale(num, inv(den)), U[p])
            allp = emul(allp, t ^ p)
        return out ^ emul(allp, U[None])

    def lag_c(t, w):                                                 # prod_{s != t} (w - s) / (t - s)
        num, den = 1, 1
        for s in pts:
            if s != t:
                num = emul(num, w ^ s)
                den = mul(den, t ^ s)
        return escale(num, inv(den))

    def prod_all(w):
        out = 1
        for s in pts:
            out = emul(out, w ^ s)
        return out

    targets = nodes + ([('z', 0)] if zero else [])
    LC = {(g, t): lag_c(t, w) for g, w in targets for t in pts}
    PA = {g: prod_all(w) for g, w in targets}
    # exposure / spread constants are usable: every Lagrange value at an F point off the nodes is nonzero, so XX(t)
    # can be exposed on any group's word and the product spread from any target word
    # (with the zero node, t = 0 is that node: XX(0) is its own word and its target reads only the product at 0)
    assert all(LX[k] for k in LX if not (zero and k[1] == 0)) and all(MU[g] for g in MU)
    assert all(LC[k] for k in LC if k[0] != 'z') and all(PA[g] for g in PA if g != 'z')
    if zero:
        assert all(LC[('z', t)] == (1 if t == 0 else 0) for t in pts) and PA['z'] == 0 and LZ[0] == 1
    for t in pts:                                                    # U(t) on a U word: its own coefficient != 0
        coef = {p: u_at({q: int(q == p) for q in (0, 1, 4, None)}, t) & MASK for p in (0, 1, 4, None)}
        assert (t in fin and coef[t] == 1 and sum(map(bool, coef.values())) == 1) or all(coef.values()), t
    bad = 0
    for _ in range(trials):
        U = {p: rng.randrange(NF) for p in (0, 1, 4, None)}
        X = {g: rng.randrange(NF * NF) for g, _ in nodes}            # X_g^2, any E value
        xz = rng.randrange(NF)
        want = {g: emul(u_at(U, w), X[g]) for g, w in nodes}         # the line: 14 products of 99
        if zero:
            want['z'] = mul(U[0], xz)                                #           + one of 42
        prods = 0
        got = {g: 0 for g, _ in targets}
        for t in pts:
            xt = 0
            for g, _ in nodes:
                xt ^= emul(X[g], LX[(g, t)]) >> M                    # Tr(X_g lambda_g(t)): F-linear in X_g's two words
            if zero:
                xt ^= mul(xz, LZ[t] & MASK)
            p = mul(u_at(U, t) & MASK, xt)                            # ONE 42-AND product
            prods += 1
            for g, _ in targets:
                got[g] ^= escale(LC[(g, t)], p)
        xinf = 0
        for g, _ in nodes:
            xinf ^= emul(X[g], MU[g]) >> M
        if zero:
            xinf ^= mul(xz, MUZ & MASK)
        p = mul(U[None], xinf)
        prods += 1
        for g, _ in targets:
            got[g] ^= escale(PA[g], p)
        if zero:
            got['z'] &= MASK
        bad += got != want
    line = 14 * 99 + (42 if zero else 0)
    print(f'blockpoly: zero node {"joined" if zero else "out"}, trials {trials}, mismatches {bad}; '
          f'products {prods} x 42 = {prods * 42} against the line {line}: {prods * 42 - line:+d} per block run')
    print('BLOCKPOLY PASS' if not bad else 'BLOCKPOLY FAIL')


def cmd_blockloop(trials, seed, zero):
    """Patch P as a SCHEDULE of the kit's primitives on word registers: scale(word, c), cx(src, dst, c) (dst ^= c src)
    and one mul_xor0 per point. Registers: the four U words (keys 0, 1, 4, None), per group the two F words of
    N_g = X_g^2 and of L_g, and with the zero node its F words zX = (x_0 + a_0)^2 and zL. Checks that every
    register is restored and L_g ^= U(w_g) X_g^2 / sc_g (zL ^= U(0) zX / (k0^2 k1^2)), and counts the linear ops."""
    rng = random.Random(seed)
    k0, k1 = KAP[0], KAP[1]
    sq = lambda v: mul(v, v)
    nodes = [(g, esq(CS[g] | NF)) for g in HOM]
    allnodes = [w for _, w in nodes] + [econj(w) for _, w in nodes] + ([0] if zero else [])
    npts = len(allnodes) + 2
    pts = list(range(npts))
    fin = [0, 1, 4]

    def lam_x(node, t):                                             # t None: the leading-coefficient functional
        num = den = 1
        for nu in allnodes:
            if nu != node:
                if t is not None:
                    num = emul(num, t ^ nu)
                den = emul(den, node ^ nu)
        return emul(num, einv(den))

    def lam_u(t):                                                    # U(t) = sum_p lam_u[p] U_p, t in F
        out = {}
        allp = 1
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

    def spread_c(t, w):                                              # coefficient of p_t in (U XX)(w), w in E
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

    sc = {g: mul(sq(k0 ^ KAP[g]), sq(k1 ^ KAP[g])) for g, _ in nodes}
    zc = mul(sq(k0), sq(k1))
    ops = dict(scale=0, cx=0, mul=0)
    bad = 0
    for _ in range(trials):
        R = {('U', p): rng.randrange(NF) for p in (0, 1, 4, None)}
        for g, _ in nodes:
            for h in (0, 1):
                R[('N', g, h)] = rng.randrange(NF)
                R[('L', g, h)] = rng.randrange(NF)
        if zero:
            R[('N', 'z', 0)] = rng.randrange(NF)
            R[('L', 'z', 0)] = rng.randrange(NF)
        R0 = dict(R)

        def scale(k, c):
            assert c
            R[k] = mul(R[k], c)
            ops['scale'] += 1

        def cx(src, dst, c):
            if c:
                R[dst] ^= mul(R[src], c)
                ops['cx'] += 1

        for t in pts + [None]:
            # --- expose x = U(t) on one U word
            if t is None:
                ukey, lu = ('U', None), None
            elif t in fin:
                ukey, lu = ('U', t), None
            else:
                lu = lam_u(t)
                ukey = ('U', 1)
                scale(ukey, lu[1])
                for p in (0, 4, None):
                    cx(('U', p), ukey, lu[p])
            # --- expose y = XX(t) on one N word: coefficients (A, B) of (n0, n1) are (lam1, lam0 + lam1)
            cf = {}
            for g, w in nodes:
                lam = lam_x(w, t)
                cf[('N', g, 0)] = lam >> M
                cf[('N', g, 1)] = (lam & MASK) ^ (lam >> M)
            if zero:
                cf[('N', 'z', 0)] = lam_x(0, t) & MASK
            ykey = next(k for k in cf if cf[k])
            scale(ykey, cf[ykey])
            for k in cf:
                if k != ykey:
                    cx(k, ykey, cf[k])
            # --- spread, one product, spread back
            cc = {}
            for g, w in nodes:
                c = escale(spread_c(t, w), inv(sc[g]))
                cc[('L', g, 0)], cc[('L', g, 1)] = c & MASK, c >> M
            if zero:
                cc[('L', 'z', 0)] = mul(spread_c(t, 0) & MASK, inv(zc))
            tkey = next(k for k in cc if cc[k])
            scale(tkey, inv(cc[tkey]))
            for k in cc:
                if k != tkey:
                    cx(tkey, k, cc[k])
            R[tkey] ^= mul(R[ukey], R[ykey])                         # the one 42-AND product
            ops['mul'] += 1
            for k in cc:
                if k != tkey:
                    cx(tkey, k, cc[k])
            scale(tkey, cc[tkey])
            # --- un-expose
            for k in cf:
                if k != ykey:
                    cx(k, ykey, cf[k])
            scale(ykey, inv(cf[ykey]))
            if lu is not None:
                for p in (0, 4, None):
                    cx(('U', p), ukey, lu[p])
                scale(ukey, inv(lu[1]))
        # --- check
        Uw = lambda w: (emul(escale(emul(w ^ 1, w ^ 4), inv(mul(1, 4))), R0[('U', 0)])
                        ^ emul(escale(emul(w, w ^ 4), inv(mul(1, 1 ^ 4))), R0[('U', 1)])
                        ^ emul(escale(emul(w, w ^ 1), inv(mul(4, 4 ^ 1))), R0[('U', 4)])
                        ^ emul(emul(emul(w, w ^ 1), w ^ 4), R0[('U', None)]))
        ok = all(R[k] == R0[k] for k in R0 if k[0] in ('U', 'N'))
        for g, w in nodes:
            X = R0[('N', g, 0)] | (R0[('N', g, 1)] << M)
            want = (R0[('L', g, 0)] | (R0[('L', g, 1)] << M)) ^ escale(emul(Uw(w), X), inv(sc[g]))
            ok &= (R[('L', g, 0)] | (R[('L', g, 1)] << M)) == want
        if zero:
            ok &= R[('L', 'z', 0)] == R0[('L', 'z', 0)] ^ mul(mul(R0[('U', 0)], R0[('N', 'z', 0)]), inv(zc))
        bad += not ok
    per = {k: v // trials for k, v in ops.items()}
    print(f'blockloop: zero node {"joined" if zero else "out"}, trials {trials}, mismatches {bad}; per block call: '
          f'{per["mul"]} products of 42, {per["cx"]} constant-multiply XORs and {per["scale"]} scalings of 12-bit words')
    print('BLOCKLOOP PASS' if not bad else 'BLOCKLOOP FAIL')


def cmd_blockcoef(trials, seed):
    """Patch C: the block's four U words in 11 products at 60 clean qubits (the line: 13). Schedule, with the clean
    words live after each step (12 qubits each):
      1 Pa, Pb                        (2 products)  [Pa Pb]
      2 V0, V2, V1 (Karatsuba)        (3)           [Pa Pb V0 V1 V2]   g_i = V_i + lin; G fan-ins from here
      3 erase Pa, Pb (free: b_4 a_3^2, b_4 b_3^2)   [g0 g1 g2]
      4 U1~, U4~ = lin1(p) + h_1(p)^2 G(p)  (2)     [g0 g1 g2 U1 U4]
      5 measure g1 (mask m1)                        [g0 g2 U1 U4]
      6 U0~ = lin1(0) + a_1^2 g0  (1); measure g0   [g2 U1 U4 U0]
      7 Uinf~ = b_1 + b_1^2 g2    (1); measure g2   [U1 U4 U0 Uinf]
      8 Pa again (1): y0 live: U(p) ^= Q0^2(p) y0; pay m0 and the b_2^2 y0 part of m1; erase Pa (free)
      9 Pb again (1): y1 live: U(p) ^= p Q0^2(p) y1, Uinf ^= y1; pay m2 and the a_2^2 y1 part of m1; erase Pb
    Checks the four words against the block's point chain and tracks the peak number of clean words."""
    rng = random.Random(seed)
    k0, k1 = KAP[0], KAP[1]
    sq = lambda v: mul(v, v)
    bad = 0
    for _ in range(trials):
        a1, b1, a2, b2, a3, b3, b4 = (rng.randrange(NF) for _ in range(7))

        def chain(p):
            if p is None:
                Y = b3 ^ mul(b4, sq(b3))
                Gv = b2 ^ b4 ^ mul(Y, sq(b2))
                return b1 ^ Y ^ mul(Gv, sq(b1))
            beta = fsqrt(p)
            x, y, s = a1 ^ mul(beta, b1), a2 ^ mul(beta, b2), a3 ^ mul(beta, b3)
            q0, q1 = sq(p) ^ p ^ sq(k0), sq(p) ^ p ^ sq(k1)
            Y = a3 ^ mul(k1 ^ p, b3) ^ mul(b4, sq(s))
            Gv = mul(q1, b2 ^ b4) ^ mul(sq(y), Y)
            return mul(q1, a1 ^ mul(k0 ^ p, b1)) ^ mul(q0, Y) ^ mul(sq(x), Gv)

        prods, live, peak = 0, set(), 0

        def alloc(*names):
            nonlocal peak
            live.update(names)
            peak = max(peak, len(live))

        alloc('Pa', 'Pb')
        Pa, Pb = mul(b4, sq(a3)), mul(b4, sq(b3)); prods += 2
        y0, y1 = a3 ^ mul(k1, b3) ^ Pa, b3 ^ Pb
        alloc('V0', 'V2')
        V0, V2 = mul(sq(a2), y0), mul(sq(b2), y1); prods += 2
        alloc('V1')
        V1 = mul(sq(a2) ^ sq(b2), y0 ^ y1) ^ V0 ^ V2; prods += 1
        g0, g1, g2 = mul(sq(k1), b2 ^ b4) ^ V0, b2 ^ b4 ^ V1, b2 ^ b4 ^ V2
        live -= {'Pa', 'Pb'}
        U = {}
        for p in (1, 4):
            alloc(('U', p))
            q1 = sq(p) ^ p ^ sq(k1)
            Gp = g0 ^ mul(p, g1) ^ mul(sq(p), g2)
            U[p] = mul(q1, a1 ^ mul(k0 ^ p, b1)) ^ mul(sq(a1) ^ mul(p, sq(b1)), Gp); prods += 1
        live -= {'V1'}
        alloc(('U', 0))
        U[0] = mul(sq(k1), a1 ^ mul(k0, b1)) ^ mul(sq(a1), g0); prods += 1
        live -= {'V0'}
        alloc(('U', None))
        U[None] = b1 ^ mul(sq(b1), g2); prods += 1
        live -= {'V2'}
        alloc('Pa')
        y0 = a3 ^ mul(k1, b3) ^ mul(b4, sq(a3)); prods += 1
        for p in (0, 1, 4):
            U[p] ^= mul(sq(p) ^ p ^ sq(k0), y0)
        live -= {'Pa'}
        alloc('Pb')
        y1 = b3 ^ mul(b4, sq(b3)); prods += 1
        for p in (1, 4):
            U[p] ^= mul(mul(sq(p) ^ p ^ sq(k0), p), y1)
        U[None] ^= y1
        live -= {'Pb'}
        for p in (0, 1, 4, None):
            bad += U[p] != chain(p)
    print(f'blockcoef: trials {trials}, mismatches {bad}; {prods} products (the line: 13), peak {peak} clean words '
          f'= {12 * peak} qubits')
    print('BLOCKCOEF PASS' if not bad else 'BLOCKCOEF FAIL')


def cmd_blockerase(trials, seed):
    """Patch T: on the RETURN block call the clean U words are X-measured first and their masks paid with Y and G
    held as COEFFICIENT words (y_0, y_1: two products; g_0, g_1, g_2: three products), every U(p) then being
    lin(p) + Q0^2(p) Y(p) + h_1(p)^2 G(p) with Y(p), G(p) free linear forms: 5 products for all clean points
    against 2 per point. Checks the coefficient forms against the block's own point chain (prefix5's docstring)."""
    rng = random.Random(seed)
    k0, k1 = KAP[0], KAP[1]
    sq = lambda v: mul(v, v)
    bad = 0
    for _ in range(trials):
        a1, b1, a2, b2, a3, b3, b4 = (rng.randrange(NF) for _ in range(7))
        # the block's point chain (three chained products per point)
        def chain(p):
            if p is None:
                Y = b3 ^ mul(b4, sq(b3))
                Gv = b2 ^ b4 ^ mul(Y, sq(b2))
                return b1 ^ Y ^ mul(Gv, sq(b1))
            beta = fsqrt(p)
            x, y, s = a1 ^ mul(beta, b1), a2 ^ mul(beta, b2), a3 ^ mul(beta, b3)
            q0, q1 = sq(p) ^ p ^ sq(k0), sq(p) ^ p ^ sq(k1)
            Y = a3 ^ mul(k1 ^ p, b3) ^ mul(b4, sq(s))
            Gv = mul(q1, b2 ^ b4) ^ mul(sq(y), Y)
            return mul(q1, a1 ^ mul(k0 ^ p, b1)) ^ mul(q0, Y) ^ mul(sq(x), Gv)
        # coefficient words: 2 + 3 products
        prods = 0
        Pa, Pb = mul(b4, sq(a3)), mul(b4, sq(b3)); prods += 2
        y0, y1 = a3 ^ mul(k1, b3) ^ Pa, b3 ^ Pb
        V0, V2 = mul(sq(a2), y0), mul(sq(b2), y1); prods += 2
        V1 = mul(sq(a2) ^ sq(b2), y0 ^ y1) ^ V0 ^ V2; prods += 1
        g0, g1, g2 = mul(sq(k1), b2 ^ b4) ^ V0, b2 ^ b4 ^ V1, b2 ^ b4 ^ V2
        def coef(p):
            # U(p) = lin(p) + Q0^2(p) Y(p) + h_1(p)^2 G(p); the last term is the only product, and it is a
            # bilinear phase between two exposed live words when U(p) has been measured
            if p is None:
                return b1 ^ y1 ^ mul(g2, sq(b1))
            q0, q1 = sq(p) ^ p ^ sq(k0), sq(p) ^ p ^ sq(k1)
            Yp = y0 ^ mul(p, y1)
            Gp = g0 ^ mul(p, g1) ^ mul(sq(p), g2)
            return mul(q1, a1 ^ mul(k0 ^ p, b1)) ^ mul(q0, Yp) ^ mul(sq(a1) ^ mul(p, sq(b1)), Gp)
        for p in (0, 1, 4, None):
            bad += chain(p) != coef(p)
    print(f'blockerase: trials {trials}, mismatches {bad}; coefficient words cost {prods} products for all clean '
          f'U words; the line recomputes Y and G per clean point (2 each)')
    print('BLOCKERASE PASS' if not bad else 'BLOCKERASE FAIL')


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'ident'
    if mode == 'blockcoef':
        cmd_blockcoef(int(sys.argv[2]) if len(sys.argv) > 2 else 2000, int(sys.argv[3]) if len(sys.argv) > 3 else 1)
        sys.exit(0)
    if mode == 'blockloop':
        cmd_blockloop(int(sys.argv[2]) if len(sys.argv) > 2 else 20, int(sys.argv[3]) if len(sys.argv) > 3 else 1,
                      len(sys.argv) > 4 and sys.argv[4] == 'zero')
        sys.exit(0)
    if mode == 'blockerase':
        cmd_blockerase(int(sys.argv[2]) if len(sys.argv) > 2 else 2000, int(sys.argv[3]) if len(sys.argv) > 3 else 1)
        sys.exit(0)
    if mode == 'blockpoly':
        cmd_blockpoly(int(sys.argv[2]) if len(sys.argv) > 2 else 200, int(sys.argv[3]) if len(sys.argv) > 3 else 1,
                      len(sys.argv) > 4 and sys.argv[4] == 'zero')
        sys.exit(0)
    if mode == 'ident':
        cmd_ident(int(sys.argv[2]) if len(sys.argv) > 2 else 20, int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    elif mode == 'count':
        cmd_count()
    elif mode == 'census':
        shots = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
        seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
        groups = [int(t) for t in sys.argv[4].split(',')] if len(sys.argv) > 4 else list(HOM)
        cmd_census(shots, seed, groups)
