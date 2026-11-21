#!/usr/bin/env python
"""Builder D1: classical check that the prefix block's point-0 values serve the node z = 0 (the zero node).

The zero node's chain is the extraction chain at gamma = 0 in GF(2^12): level i has q_i = kappa_(i//2) (z^2 + z +
kappa at z = 0), c_i = sum_(j<i) 1 / q_j, head value h_i(0) = a_i:
    L <- q (L + (c a + b) D^2),   (N, D) <- (q D, N + a D);   level 0 is the kit's first step.
Claim: skipping the L-products of levels 1..4 and adding to L_0, before the chain,
    U(0) / (kappa_0^2 kappa_1^2) * (x + a_0)^2 + G(0) / kappa_1^2,
with (Y, G, U) at alpha = 0 from the block's own formulas (tools/dqi_build_prefix5_check.point_values), gives the
same (N, D, L) at the end of the chain. Usage: python dqi_pub_zcheck.py [trials=500]
"""
import random
import sys

from dqi_research_model import mul, inv, KAP
from dqi_build_prefix5_check import point_values


def chain0(heads, x, L0, skip=()):
    """The zero node's chain over heads 0..33 at z = 0, homogeneous start (N, D, L) = (x, 1, L0)."""
    N, D, L = x, 1, L0
    c = 0
    for i, (a, b) in enumerate(heads):
        q = KAP[i // 2]
        if i not in skip:
            L = mul(L ^ mul(mul(c, a) ^ b, mul(D, D)), q)
        else:
            L = mul(L, q)
        N, D = mul(D, q), N ^ mul(a, D)
        c ^= inv(q)
    return N, D, L


def main():
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    rng = random.Random(13)
    ok = 0
    for t in range(trials):
        heads = [(rng.randrange(4096), rng.randrange(4096)) if rng.random() < 0.9 else (0, 0) for _ in range(34)]
        x, L0 = rng.randrange(4096), rng.randrange(4096)
        true = chain0(heads, x, L0)
        Y, G, U = point_values(heads, 0)
        k0, k1 = KAP[0], KAP[1]
        X = x ^ heads[0][0]
        corr = mul(mul(U, inv(mul(mul(k0, k0), mul(k1, k1)))), mul(X, X)) ^ mul(G, inv(mul(k1, k1)))
        got = chain0(heads, x, L0 ^ corr, skip={1, 2, 3, 4})
        ok += got == true
    print(f'zero node joins the prefix block (point 0): {ok} of {trials} {"PASS" if ok == trials else "FAIL"}')


if __name__ == '__main__':
    main()
