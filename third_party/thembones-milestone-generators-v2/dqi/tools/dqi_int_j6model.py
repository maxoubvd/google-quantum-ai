#!/usr/bin/env python3
"""dqi_int_j6model.py -- workstream dqi_int (B5): the hoisted-L prefix at J = 6 (R37-DQI-EXT ledger row 22), model only.

For the skip set S of chain steps whose L-products are hoisted (J = 5: {1, 2, 3, 4}, the line; J = 6: {1, ..., 5}),
the correction is L0 += M_g X^2 + V_g (X = x + H_0, D_i affine in X), M_g = sum_i w_i alpha_i^2, V_g = sum_i w_i
beta_i^2, w_i = (c_i H_i + b_i) / prod_{j<i} q_j. Checked against the true chain (chain_forward / chain_skip of the
research model). Then, for candidate products of Qk^2 = (kappa_k + kappa_g)^2 (k = 0, 1, 2), the lowest degree d such
that M_g Q-product (and V_g Q-product) is ONE polynomial in w_g = gamma_g^2 with GF(2^12) coefficients for all 14 groups
(interpolation through conjugate node pairs, verified on the other groups): the number of U (G) point values the block
would hold.
Usage: python dqi_int_j6model.py [trials=20] [seed=1]
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dqi_research_model import (mul, inv, emul, esq, escale, head_value, chain_params, chain_forward, CS, KAP,  # noqa
                                GAMMA, G)

M12 = 12
MASK = (1 << M12) - 1


def econj(u):
    return ((u & MASK) ^ (u >> M12)) | (u & ~MASK)


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


def correction(heads, g, S):
    """(M_g, V_g) for skip set S, by the affine tracking of D_i in X."""
    params = chain_params(g)
    q0, _ = params[0]
    Na, Nb = 0, q0          # N_1 = q0 (alpha 0, beta q0)
    Da, Db = 1, 0           # D_1 = X
    pref = q0               # prod_{j < i} q_j for i = 1
    M = V = 0
    for i in range(1, max(S) + 1):
        q, c = params[i]
        a, b = heads[i]
        H = head_value(a, b, CS[g])
        if i in S:
            w = escale(escale(H, c) ^ b, inv(pref)) if False else emul(escale(H, c) ^ b, einv_f(pref))
            M ^= emul(w, esq(Da))
            V ^= emul(w, esq(Db))
        Na, Nb, Da, Db = escale(Da, q), escale(Db, q), Na ^ emul(H, Da), Nb ^ emul(H, Db)
        pref = mul(pref, q)
    return M, V


def einv_f(s):
    return inv(s)           # an F scalar: q_j are F constants (kappa_k + kappa_g)


def lagrange_eval(nodes, vals, t):
    out = 0
    for i, (ni, vi) in enumerate(zip(nodes, vals)):
        num, den = 1, 1
        for j, nj in enumerate(nodes):
            if j != i:
                num = emul(num, t ^ nj)
                den = emul(den, ni ^ nj)
        out ^= emul(emul(vi, num), einv(den))
    return out


def einv(u):
    a, b = u & MASK, u >> M12
    K = 2048
    n = mul(a, a) ^ mul(a, b) ^ mul(K, mul(b, b))
    ni = inv(n)
    return mul(a ^ b, ni) | (mul(b, ni) << M12)


def min_degree(groups, w, v, maxd=12):
    nodes, vals = [], []
    for g in groups:
        nodes += [w[g], econj(w[g])]
        vals += [v[g], econj(v[g])]
    for d in range(0, maxd + 1):
        k = d + 1
        if k > len(nodes) - 2:
            break
        ok = all(lagrange_eval(nodes[:k], vals[:k], nodes[j]) == vals[j] for j in range(k, len(nodes)))
        if ok:
            return d
    return None


def main():
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    groups = list(range(3, G))
    bad = 0
    degs = {}
    for t in range(trials):
        heads = [(rng.randrange(4096), rng.randrange(4096)) for _ in range(2 * G)]
        for S, name in (({1, 2, 3, 4}, 'J5'), ({1, 2, 3, 4, 5}, 'J6')):
            w, Mv, Vv = {}, {}, {}
            for g in groups:
                hg = heads[:2 * g]
                x, L0 = rng.randrange(1 << 24), rng.randrange(1 << 24)
                M, V = correction(hg, g, S)
                X = x ^ head_value(hg[0][0], hg[0][1], CS[g])
                if chain_skip(hg, g, x, 1, L0 ^ emul(M, esq(X)) ^ V, S) != chain_forward(hg, g, x, 1, L0):
                    bad += 1
                w[g] = esq(GAMMA[g])
                Mv[g], Vv[g] = M, V
            for qs in ((0, 1), (0, 1, 2), (1,), (1, 2), (2,), (0, 2)):
                for label, val in (('M', Mv), ('V', Vv)):
                    scaled = {}
                    for g in groups:
                        f = 1
                        for k in qs:
                            f = mul(f, mul(KAP[k] ^ KAP[g], KAP[k] ^ KAP[g]))
                        scaled[g] = escale(val[g], f)
                    d = min_degree(groups, w, scaled)
                    key = (name, label, qs)
                    degs.setdefault(key, set()).add(d)
    print(f'J=5 / J=6 corrections checked against the true chain: {trials} x 14 groups each, mismatches {bad}')
    for key in sorted(degs, key=str):
        print(key, 'degree', sorted(degs[key], key=lambda d: (d is None, d)))


if __name__ == '__main__':
    main()
