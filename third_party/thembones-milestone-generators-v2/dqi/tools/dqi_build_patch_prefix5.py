"""Patch for the q26 generator: hoisted-L prefix, J = 5, point-value form (workstream dqi_build).

Identity (tools/dqi_build_prefix5_check.py): for a later node g the L-products of chain steps 1..4 equal
    L0 += U(w_g) / (Q0g^2 Q1g^2) * (x + H_0)^2 + G(w_g) / Q1g^2,   w_g = gamma_g^2,
with the polynomials Y = mu_3 + b_4 h_3^(2), G = Q1^2 (b_2 + b_4) + Y h_2^(2), U = Q1^2 mu_1 + Q0^2 Y + G h_1^(2).
At an F' point alpha = beta^2 the values are a chain of three 12-bit products:
    Y = a_3 + (k1 + alpha) b_3 + b_4 s^2,  G = Q1^2(alpha)(b_2 + b_4) + y^2 Y,
    U = Q1^2(alpha)(a_1 + (k0 + alpha) b_1) + Q0^2(alpha) Y + x^2 G,   (x, y, s) = h_1, h_2, h_3 at beta.
U is held as its values at alpha = 1, 4, infinity (three clean words) and alpha = 0 (in place in the b_1 word).
Chain steps with heads 1..4 become N-only products. Peak clean use of the block: 48 qubits.
Replaces dqi_build_patch_prefix3 (same hooks); do not install both.
"""
import sys

import dqi_build_patch_prefix3 as p3

NONLY = frozenset({1, 2, 3, 4})
TEMP_POINTS = (1, 4, None)        # alpha of the clean U words in allocation order; the last keeps G in b_4
BETA = {0: 0, 1: 1, 4: 2}         # beta = sqrt(alpha)


def make_block(G):
    ar, cx, scale, scalepair, combine = G['ar'], G['cx'], G['scale'], G['scalepair'], G['combine']
    groups, K, mul, inv = G['groups'], G['K'], G['mul'], G['inv']
    square12, square24, mul24 = G['square12'], G['external_square24'], G['mul24_xor99']
    k0, k1 = groups[0]['kap'], groups[1]['kap']
    assert G['ks'][0] == k0 and G['ks'][2] == k1 and k0 and k1

    def sq(v):
        return mul(v, v)

    def fmul(u, v):
        t = mul(u[1], v[1])
        return (mul(u[0], v[0]) ^ mul(K, t), mul(u[0], v[1]) ^ mul(u[1], v[0]) ^ t)

    def fscale(u, c):
        return (mul(u[0], c), mul(u[1], c))

    def lagrange(points, w):
        fin = [p for p in points if p is not None]
        lam = {}
        for p in fin:
            num, den = (1, 0), 1
            for r in fin:
                if r != p:
                    num = fmul(num, (w[0] ^ r, w[1]))
                    den = mul(den, p ^ r)
            lam[p] = fscale(num, inv(den))
        if None in points:
            num = (1, 0)
            for r in fin:
                num = fmul(num, (w[0] ^ r, w[1]))
            lam[None] = num
        return lam

    def consts(p):
        if p is None:
            return None, 1, 1
        return BETA[p], sq(p) ^ p ^ sq(k0), sq(p) ^ p ^ sq(k1)

    def block(b, words, homogeneous):
        a0, h0, a1, b1, a2, b2, a3, b3, b4 = (words[i] for i in (0, 1, 2, 3, 4, 5, 6, 7, 9))
        hom = sorted(homogeneous)
        print(f'dqb_prefix5 block entry: allocated {b.num_allocated_qubits}', file=sys.stderr)

        def form(a, bb, beta, body):
            """Run body(reg) while reg holds (a + beta bb)^2 (beta None: bb^2)."""
            reg = bb if beta is None else a
            if beta:
                cx(b, bb, a, beta)
            square12(b, reg)
            body(reg)
            square12(b, reg, True)
            if beta:
                cx(b, bb, a, beta)

        def lin_y(p, t):
            if p is None:
                b.cx(b3, t)
            else:
                b.cx(a3, t); cx(b, b3, t, k1 ^ p)

        def lin_u(p, t, q1):
            if p is None:
                b.cx(b1, t)
            else:
                cx(b, a1, t, q1); cx(b, b1, t, mul(q1, k0 ^ p))

        def temp_point(p, U, compute, inplace):
            beta, q0, q1 = consts(p)
            Y = b.alloc_qubits(12)
            lin_y(p, Y)
            form(a3, b3, beta, lambda r: ar.mul_xor0(b, b4, r, Y))
            if inplace:
                b.cx(b2, b4); scale(b, b4, q1)
                Gr = b4
            else:
                Gr = b.alloc_qubits(12)
                cx(b, b2, Gr, q1); cx(b, b4, Gr, q1)
            form(a2, b2, beta, lambda r: ar.mul_xor0(b, Y, r, Gr))
            lin_u(p, U, q1); cx(b, Y, U, q0)
            if compute:
                form(a1, b1, beta, lambda r: ar.mul_xor0(b, Gr, r, U))
                for gi in hom:
                    g = groups[gi]
                    L = combine(words[g['word'] + 2:g['word'] + 4])
                    lam = lagrange(TEMP_POINTS, (g['c'] ^ g['kap'], 1))[p]
                    c = fscale(lam, inv(sq(k1 ^ g['kap'])))
                    cx(b, Gr, L[:12], c[0]); cx(b, Gr, L[12:], c[1])
            else:
                form(a1, b1, beta, lambda r: ar.erase_product(b, Gr, r, U))
            if inplace:
                form(a2, b2, beta, lambda r: ar.mul_xor0(b, Y, r, Gr))
                scale(b, b4, inv(q1)); b.cx(b2, b4)
            else:
                cx(b, b2, Gr, q1); cx(b, b4, Gr, q1)
                form(a2, b2, beta, lambda r: ar.erase_product(b, Y, r, Gr))
                b.free(Gr)
            lin_y(p, Y)
            form(a3, b3, beta, lambda r: ar.erase_product(b, b4, r, Y))
            b.free(Y)

        def sqmul(src, fac, dst):
            square12(b, src); ar.mul_xor0(b, fac, src, dst); square12(b, src, True)

        def point0(compute):
            """alpha = 0 in place: b_3 <- Y(0), b_4 <- G(0), b_1 <- U(0)."""
            k0s, k1s = sq(k0), sq(k1)
            if compute:
                scale(b, b3, k1); b.cx(a3, b3); sqmul(a3, b4, b3)
                b.cx(b2, b4); scale(b, b4, k1s); sqmul(a2, b3, b4)
                scale(b, b1, k0); b.cx(a1, b1); scale(b, b1, k1s); cx(b, b3, b1, k0s); sqmul(a1, b4, b1)
            else:
                sqmul(a1, b4, b1); cx(b, b3, b1, k0s); scale(b, b1, inv(k1s)); b.cx(a1, b1); scale(b, b1, inv(k0))
                sqmul(a2, b3, b4); scale(b, b4, inv(k1s)); b.cx(b2, b4)
                sqmul(a3, b4, b3); b.cx(a3, b3); scale(b, b3, inv(k1))

        def mix(R1, R2, u, v, undo):
            """(R1, R2) <- F element u R1 + v R2 as (lo, hi) = (u0 R1 + v0 R2, u1 R1 + v1 R2)."""
            p, q, r, s = u[0], v[0], u[1], v[1]
            det = mul(p, s) ^ mul(q, r)
            assert p and det
            if not undo:
                scale(b, R1, p); cx(b, R2, R1, q); scale(b, R2, mul(det, inv(p))); cx(b, R1, R2, mul(r, inv(p)))
            else:
                cx(b, R1, R2, mul(r, inv(p))); scale(b, R2, mul(p, inv(det))); cx(b, R2, R1, q); scale(b, R1, inv(p))

        with b.mark('dqb_prefix5'):
            U = {}
            for p in TEMP_POINTS:
                U[p] = b.alloc_qubits(12)
                temp_point(p, U[p], True, p == TEMP_POINTS[-1])
            point0(True)
            R1, R2, R3 = U[1], U[4], U[None]
            m = combine([R1, R2])
            for gi in hom:
                g = groups[gi]
                ct, kap = g['c'], g['kap']
                N = combine(words[g['word']:g['word'] + 2])
                L = combine(words[g['word'] + 2:g['word'] + 4])
                lam = lagrange((0, 1, 4, None), (ct ^ kap, 1))
                sc = mul(sq(k0 ^ kap), sq(k1 ^ kap))
                mix(R1, R2, lam[1], lam[4], False)
                cx(b, R3, R1, lam[None][0]); cx(b, R3, R2, lam[None][1])
                cx(b, b1, R1, lam[0][0]); cx(b, b1, R2, lam[0][1])
                cx(b, h0, a0, ct); b.cx(combine([a0, h0]), N); cx(b, h0, a0, ct)
                square24(b, N)
                scalepair(b, L, sc); mul24(b, m, N, L, K); scalepair(b, L, inv(sc))
                square24(b, N, True)
                cx(b, h0, a0, ct); b.cx(combine([a0, h0]), N); cx(b, h0, a0, ct)
                cx(b, b1, R2, lam[0][1]); cx(b, b1, R1, lam[0][0])
                cx(b, R3, R2, lam[None][1]); cx(b, R3, R1, lam[None][0])
                mix(R1, R2, lam[1], lam[4], True)
            point0(False)
            for p in reversed(TEMP_POINTS):
                temp_point(p, U[p], False, p == TEMP_POINTS[-1])
                b.free(U[p])
    return block


def install(km, module):
    walk = sys.modules['walk_bound']
    fns = walk.functions(module.decode)
    scope = {}
    exec(compile(p3.CHAIN_SRC, __file__ + ':chain', 'exec'), scope)
    chain_code = scope['chain_step'].__code__
    scope = {}
    exec(compile(p3.EXTRACTION_SRC, __file__ + ':extraction', 'exec'), scope)
    extraction_code = scope['external_extraction'].__code__
    need = ('ar', 'cx', 'scale', 'scalepair', 'combine', 'groups', 'ks', 'K', 'mul', 'inv', 'square12',
            'external_square24', 'mul24_xor99')
    mul_head = next(f for f in fns if f.__name__ == 'external_mul_head' and all(k in f.__globals__ for k in need))
    block = make_block(mul_head.__globals__)
    log = dict(chain=0, extraction=0)
    for f in fns:
        if f.__name__ == 'chain_step' and 'joint_tensor' in f.__code__.co_names:
            assert f.__closure__ is None
            f.__globals__.update(_dqb_nonly=NONLY, _dqb_mul_head=mul_head)
            f.__code__ = chain_code
            log['chain'] += 1
        elif f.__name__ == 'external_extraction' and 'inverse24_forward190' in f.__code__.co_names:
            assert f.__closure__ is None
            f.__globals__['_dqb_block'] = block
            f.__code__ = extraction_code
            log['extraction'] += 1
    assert log['chain'] and log['extraction'], log
    return log
