#!/usr/bin/env python
"""Classical check of the hoisted-L prefix (J = 3) for the homogeneous chain (workstream dqi_build).

Claim: for a later node g, running the chain WITHOUT the L-products of steps 1 and 2 on the jet whose
derivative was first replaced by  L0 + U_g * (x + H_0)^2 + b_2  gives exactly the true (N, D, L), where
    U_g = [a_1 + b_1 (gamma^2 + kappa_0) + b_2 (a_1^2 + b_1^2 gamma^2)] / q_0(gamma)^2
is GF(2)-linear in the four words (a_1, b_1, b_2 a_1^2, b_2 b_1^2).
Checked on random heads (including zero heads) and all groups g >= 2.
"""
import random
from dqi_research_model import *


def chain_skip(heads, g, N, D, L, skip):
    for i, (q, c) in enumerate(chain_params(g)):
        a, b = heads[i]
        H = head_value(a, b, CS[g])
        if i in skip:
            L = escale(L, q)
        else:
            L = escale(L ^ emul(escale(H, c) ^ b, esq(D)), q)
        N, D = escale(D, q), N ^ emul(H, D)
    return N, D, L


rng = random.Random(3)
n = 0
for trial in range(300):
    g = rng.randrange(2, G)
    heads = [(rng.randrange(4096), rng.randrange(4096)) if rng.random() < 0.85 else (0, 0) for _ in range(2 * g)]
    x, L0 = rng.randrange(1 << 24), rng.randrange(1 << 24)
    true = chain_forward(heads, g, x, 1, L0)
    (a0, b0), (a1, b1), (a2, b2) = heads[0], heads[1], heads[2]
    gam = GAMMA[g]
    g2 = esq(gam)
    q0 = KAP[0] ^ KAP[g]
    U = (a1 ^ mul(b2, mul(a1, a1))) ^ escale(g2 ^ KAP[0], b1) ^ escale(g2, mul(b2, mul(b1, b1)))
    U = escale(U, inv(mul(q0, q0)))
    D1 = x ^ head_value(a0, b0, CS[g])
    L0c = L0 ^ emul(U, esq(D1)) ^ b2
    got = chain_skip(heads, g, x, 1, L0c, {1, 2})
    assert got == true, (trial, g)
    n += 1
print(f'hoisted-L prefix J=3: {n} random cases (groups 2..16, zero heads included): PASS')
