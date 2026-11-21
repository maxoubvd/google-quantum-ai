#!/usr/bin/env python3
"""r32_dqi-ext_zt.py -- construction Z-T of R32-DQI-EXT on every error weight (classical, exact).

Z-T: the node z = 0 (35th head) is run as a telescoped group on a 12-bit unit nu of the carried scale G
(G = nu + mu y; if nu = 0 a flag swaps the halves). The chain starts from (nu u0, nu, nu^2 v0) instead of
(u0, 1, v0); nothing is measured at its end: the unit's wires hold nu' = nu * (the nonzero one of N_f, D_f).
On the return the reversed chain starts from nu' and ends at (nu u0, nu, nu^2 v0): the unit is back, so the
12-bit scalar chain that re-creates the measured gauge today (33 products up and down) is not run.

Conventions are the kit's (final_decoder/decode.py, functions step / first_step / normal_form / extraction):
  (N, D) <- (kap D, N + a D),  L <- kap (L + (c a + b) D^2),  c_i = sum_{j<i} 1/kap_j,  kap_i = q_i(0),
  pole = [D_f = 0] swaps N and D, then N and L are divided by D and D^2.
Field arithmetic, the instance sampler and the zero-safe heads of the 17 quadratic groups are imported
(read only) from tools/dqi_research_model.py and tools/dqi_research_term.py (workstream dqi_research).

Checks, for every weight 0..35 (reps instances each):
  Z1  the normalised pair (u, v) and the pole flag do not depend on nu (nu = 1 against a random unit)
  Z2  the unit after the group is never 0 (both of N_f, D_f zero would be the kit's own failure)
  Z3  the reversed chain from nu' ends at (nu u0, nu, nu^2 v0); the unit returns exactly
  Z4  with nu = 0 the swapped half does the same (the flag f = [low half = 0] is recomputable after the group:
      the low half is still 0)
Usage: python r32_dqi-ext_zt.py [reps=12] [seed=5]
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dqi_research_model import mul, inv, KAP, G, NF, instance, jets          # noqa: E402
from dqi_research_term import affine_zero_safe                              # noqa: E402

KS = [KAP[i // 2] for i in range(2 * G)]
CORR = []
_c = 0
for _k in KS:
    CORR.append(_c)
    _c ^= inv(_k)
assert _c == 0


def sq(a):
    return mul(a, a)


def zero_jet(roots, weights):
    u = v = 0
    for beta, w in zip(roots, weights):
        t = inv(beta)
        u ^= mul(t, w)
        v ^= mul(sq(t), w)
    return u, v


def forward(heads, N, D, L):
    for i in range(2 * G):
        a, b = heads[i]
        L = mul(L ^ mul(mul(a, CORR[i]) ^ b, sq(D)), KS[i])
        N, D = mul(D, KS[i]), N ^ mul(a, D)
    return N, D, L


def backward(heads, N, D, L):
    for i in reversed(range(2 * G)):
        a, b = heads[i]
        qi = inv(KS[i])
        N, D = D ^ mul(a, mul(N, qi)), mul(N, qi)
        L = mul(L, qi) ^ mul(mul(a, CORR[i]) ^ b, sq(D))
    return N, D, L


def group(heads, u0, v0, nu):
    """Forward Z-T group on the unit nu. Returns (u, v, pole, unit after) or None if both N_f, D_f vanish."""
    N, D, L = forward(heads, mul(nu, u0), nu, mul(sq(nu), v0))
    pole = int(D == 0)
    if pole:
        N, D = D, N
    if D == 0:
        return None
    di = inv(D)
    return mul(N, di), mul(L, sq(di)), pole, D


def ungroup(heads, u, v, pole, unit):
    N, D, L = mul(u, unit), unit, mul(v, sq(unit))
    if pole:
        N, D = D, N
    return backward(heads, N, D, L)


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    rng = random.Random(seed)
    ok = {k: True for k in ('Z1', 'Z2', 'Z3', 'Z4')}
    n = 0
    poles_by_weight = {}
    both_zero = 0
    for weight in range(0, 36):
        for _ in range(reps):
            roots, weights = instance(rng, weight)
            heads = affine_zero_safe(jets(roots, weights))
            u0, v0 = zero_jet(roots, weights)
            ref = group(heads, u0, v0, 1)
            if ref is None:
                both_zero += 1
                ok['Z2'] = False
                continue
            poles_by_weight.setdefault(weight, set()).add(ref[2])
            nu = rng.randrange(1, NF)
            got = group(heads, u0, v0, nu)
            ok['Z1'] &= (got is not None and got[:3] == ref[:3])
            ok['Z2'] &= (got is not None and got[3] != 0 and got[3] == mul(nu, ref[3]))
            back = ungroup(heads, got[0], got[1], got[2], got[3])
            ok['Z3'] &= (back == (mul(nu, u0), nu, mul(sq(nu), v0)))
            # Z4: carried scale with a zero low half: the flag swaps the halves, the group runs on the other one
            low, high = 0, rng.randrange(1, NF)
            f = int(low == 0)
            unit, other = (high, low) if f else (low, high)
            got4 = group(heads, u0, v0, unit)
            unit2 = got4[3]
            low2, high2 = (other, unit2) if f else (unit2, other)
            ok['Z4'] &= (got4[:3] == ref[:3] and int(low2 == 0) == f and high2 != 0)
            back4 = ungroup(heads, got4[0], got4[1], got4[2], high2 if f else low2)
            ok['Z4'] &= (back4 == (mul(unit, u0), unit, mul(sq(unit), v0)))
            n += 1
    print('instances: %d (weights 0..35, %d each), seed %d; both N_f and D_f zero: %d' % (n, reps, seed, both_zero))
    for k in ('Z1', 'Z2', 'Z3', 'Z4'):
        print('%s %s' % (k, 'PASS' if ok[k] else 'FAIL'))
    pat = ''.join('p' if poles_by_weight.get(w) == {1} else ('.' if poles_by_weight.get(w) == {0} else '?')
                  for w in range(36))
    print('pole at the node z = 0 by weight 0..35 (p = always, . = never, ? = mixed): %s' % pat)
    print('products per direction: chain 33 x (42 + 42) unchanged; scale-in 2 x 42; no scalar chain')
    print('ALL PASS' if all(ok.values()) else 'SOME FAILED')


if __name__ == '__main__':
    main()
