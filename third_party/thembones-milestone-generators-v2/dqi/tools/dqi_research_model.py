#!/usr/bin/env python
"""Exact classical model of the DQI extraction used by workstream dqi_research.

Conventions follow the kit's final decoder (final_decoder/decode.py):
  F' = GF(2^12) mod 0x1053, F = F'[y]/(y^2 + y + K), K = 2048, elements of F are a | (b << 12) = a + b y.
  17 double nodes gamma_g = 2g + y (g = 0..16), kappa_g = c^2 + c + K, q_i(z) = z^2 + z + kappa_{i//2},
  heads h_i(z) = a_i + b_i z, continued fraction f_i = h_i + q_i / f_{i+1}, f_0 = r = sum_k w_k / (z - beta_k).
  Reference (affine) jet step and the homogeneous chain step
      L <- q (L + (c H + b) D^2),   (N, D) <- (q D, N + H D),   c_i = sum_{j<i} 1/q_j(gamma_g)
  are transcribed from external_ordinary_node / external_chain_step.
The zero node (35th head) is not modelled here; it is a GF(2^12) copy of the same chain.
"""
import random

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
GAMMA = [c | NF for c in CS]


def head_value(a, b, c):
    """h(gamma) for h = a + b z at gamma = c + y."""
    return (a ^ mul(b, c)) | (b << M)


def coords(u, c):
    b = u >> M
    return (u & MASK) ^ mul(b, c), b


def jets(roots, weights):
    out = []
    for g in range(G):
        u = v = 0
        for beta, w in zip(roots, weights):
            t = einv(GAMMA[g] ^ beta)
            u ^= escale(t, w)
            v ^= escale(esq(t), w)
        out.append((u, v))
    return out


class Pole(Exception):
    pass


def reference_extract(J):
    """Affine (ordinary) extraction. Returns the 34 heads [(a, b)] or raises Pole."""
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


def chain_forward(heads, g, N, D, L):
    """Homogeneous chain over heads 0..2g-1 at node g, started from (N, D, L)."""
    for i, (q, c) in enumerate(chain_params(g)):
        a, b = heads[i]
        H = head_value(a, b, CS[g])
        L = escale(L ^ emul(escale(H, c) ^ b, esq(D)), q)
        N, D = escale(D, q), N ^ emul(H, D)
    return N, D, L


def chain_backward(heads, g, N, D, L):
    """Exact inverse of chain_forward."""
    for i, (q, c) in reversed(list(enumerate(chain_params(g)))):
        a, b = heads[i]
        H = head_value(a, b, CS[g])
        qi = inv(q)
        N, D = D ^ emul(H, escale(N, qi)), escale(N, qi)
        L = escale(L, qi) ^ emul(escale(H, c) ^ b, esq(D))
    return N, D, L


def instance(rng, weight):
    roots = rng.sample(range(1, NF), weight)
    weights = [rng.randrange(1, NF) for _ in roots]
    return roots, weights


if __name__ == '__main__':
    rng = random.Random(20261110)
    ok = poles = 0
    for weight in range(0, 36):
        for _ in range(3):
            roots, weights = instance(rng, weight)
            J = jets(roots, weights)
            try:
                heads = reference_extract(J)
            except Pole:
                poles += 1
                continue
            for g in range(1, G):
                N, D, L = chain_forward(heads, g, J[g][0], 1, J[g][1])
                assert D, 'unexpected pole in chain'
                di = einv(D)
                u = emul(N, di)
                v = emul(L, esq(di))
                a, b = heads[2 * g]
                assert u == head_value(a, b, CS[g])
                a2, b2 = heads[2 * g + 1]
                assert einv(v ^ b) == head_value(a2, b2, CS[g])
                assert chain_backward(heads, g, N, D, L) == (J[g][0], 1, J[g][1])
            ok += 1
    print(f'instances checked {ok}, affine poles {poles}')
    print('PASS: homogeneous chain reproduces the affine extraction (value and derivative) and inverts exactly')
