#!/usr/bin/env python
"""Classical check of the hoisted-L prefix J = 5 in point-value form (workstream dqi_build).

L0 += M_g (x + H_0)^2 + V_g, then the chain without the L-products of steps 1..4 equals the true chain, with
  M_g = U(w_g) / (Q0g^2 Q1g^2),  V_g = G(w_g) / Q1g^2,  w_g = gamma_g^2, and the polynomials in w
  Y = mu_3 + b_4 h_3^(2),  G = Q1^2 (b_2 + b_4) + Y h_2^(2),  U = Q1^2 mu_1 + Q0^2 Y + G h_1^(2),
  mu_i(w) = (a_i + kappa b_i) + b_i w,  h^(2)(w) = a^2 + b^2 w,  Qk^2(w) = w^2 + w + kappa_k^2.
Each polynomial is represented by its values at the F' points alpha in POINTS (None = leading coefficient):
at alpha = beta^2 the chain is three 12-bit products,
  Y = a_3 + (k1 + alpha) b_3 + b_4 s^2,  G = Q1^2(alpha)(b_2 + b_4) + y^2 Y,
  U = Q1^2(alpha)(a_1 + (k0 + alpha) b_1) + Q0^2(alpha) Y + x^2 G,   (x, y, s) = h_1, h_2, h_3 at beta.
"""
import random
from dqi_research_model import *
from dqi_build_prefix_check import chain_skip

POINTS = (0, 1, 4, None)          # alpha; beta = sqrt(alpha) = 0, 1, 2; None is the point at infinity
BETA = {0: 0, 1: 1, 4: 2}


def sq(x):
    return mul(x, x)


def point_values(heads, alpha):
    (a1, b1), (a2, b2), (a3, b3), b4 = heads[1], heads[2], heads[3], heads[4][1]
    k0, k1 = KAP[0], KAP[1]
    if alpha is None:
        Y = b3 ^ mul(b4, sq(b3))
        Gv = b2 ^ b4 ^ mul(sq(b2), Y)
        return Y, Gv, b1 ^ Y ^ mul(sq(b1), Gv)
    beta = BETA[alpha]
    x, y, s = a1 ^ mul(beta, b1), a2 ^ mul(beta, b2), a3 ^ mul(beta, b3)
    q0, q1 = sq(alpha) ^ alpha ^ sq(k0), sq(alpha) ^ alpha ^ sq(k1)
    Y = a3 ^ mul(k1 ^ alpha, b3) ^ mul(b4, sq(s))
    Gv = mul(q1, b2 ^ b4) ^ mul(sq(y), Y)
    return Y, Gv, mul(q1, a1 ^ mul(k0 ^ alpha, b1)) ^ mul(q0, Y) ^ mul(sq(x), Gv)


def lagrange(points, w):
    """F-valued weights lam[p] with P(w) = sum lam[p] * P(p) for deg P < len(points); None = leading coeff."""
    fin = [p for p in points if p is not None]
    lam = {}
    for p in fin:
        num, den = 1, 1
        for r in fin:
            if r != p:
                num = emul(num, w ^ r)
                den = mul(den, p ^ r)
        lam[p] = escale(num, inv(den))
    if None in points:
        num = 1
        for r in fin:
            num = emul(num, w ^ r)
        lam[None] = num
    return lam


def correction(heads, g, upoints=POINTS, gpoints=(1, 4, None)):
    w = esq(GAMMA[g])
    Q0, Q1 = KAP[0] ^ KAP[g], KAP[1] ^ KAP[g]
    vals = {p: point_values(heads, p) for p in POINTS}
    lu, lg = lagrange(upoints, w), lagrange(gpoints, w)
    U = 0
    for p in upoints:
        U ^= escale(lu[p], vals[p][2])
    Gw = 0
    for p in gpoints:
        Gw ^= escale(lg[p], vals[p][1])
    return escale(U, inv(mul(sq(Q0), sq(Q1)))), escale(Gw, inv(sq(Q1)))


if __name__ == '__main__':
    rng = random.Random(5)
    n = 0
    for trial in range(400):
        g = rng.randrange(3, G)
        heads = [(rng.randrange(4096), rng.randrange(4096)) if rng.random() < 0.85 else (0, 0) for _ in range(2 * g)]
        x, L0 = rng.randrange(1 << 24), rng.randrange(1 << 24)
        true = chain_forward(heads, g, x, 1, L0)
        M, V = correction(heads, g)
        D1 = x ^ head_value(heads[0][0], heads[0][1], CS[g])
        got = chain_skip(heads, g, x, 1, L0 ^ emul(M, esq(D1)) ^ V, {1, 2, 3, 4})
        assert got == true, (trial, g)
        n += 1
    print(f'hoisted-L prefix J=5 (point values at alpha = 0, 1, 4, inf): {n} random cases, groups 3..16: PASS')
