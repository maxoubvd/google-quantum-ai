#!/usr/bin/env python3
"""dqi_pub_inv27model.py -- builder D1 (holder 4): classical model of a one-operand inverse phase over the cubic
tower GF(4096) = GF(16)[theta], theta in GF(8) (theta^3 = theta + 1), GF(16) = F2[beta] (beta^4 = beta + 1).

Mechanism (re-derived here from the field; credit: the public DQI decoder's 27-Toffoli kernel, Syed Asad Ali,
Apache-2.0, read only for the idea): (-1)^(m . inv0(z)) = (-1)^Tr(alpha / z) with alpha the trace-dual of the
classical mask m; z' = z / alpha (a free classical scaling) leaves the FIXED phase (-1)^Tr(1 / z'). With
z' = a + b theta + c theta^2 (a, b, c in GF(16)), the multiplication matrix of z' gives
    e2 = a^2 + b^2 + c^2 + b c,    N = a e2 + (b^3 + b c^2 + c^3),    Tr(1 / z') = Tr16(e2 / N).
This script checks those identities on every z and prints the 8 quadratic target forms (over the 8 bits of b, c)
of the joint map (b, c) -> (b c, b^3 + b c^2 + c^3).
"""
import json
import sys

MOD = 0x1053


def mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & 4096:
            a ^= MOD
    return r


def pw(a, e):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a)
        a = mul(a, a)
        e >>= 1
    return r


def inv0(a):
    return pw(a, 4094) if a else 0


def tr(a, k=12):          # absolute trace of an element of GF(2^k) inside GF(4096)
    s, v = 0, a
    for _ in range(k):
        s ^= v
        v = mul(v, v)
    return s


def setup():
    beta = next(v for v in range(2, 4096) if pw(v, 4) ^ v ^ 1 == 0)
    theta = next(v for v in range(2, 4096) if pw(v, 3) ^ v ^ 1 == 0)
    g16 = []
    for i in range(16):
        v = 0
        for j in range(4):
            if (i >> j) & 1:
                v ^= pw(beta, j)
        g16.append(v)
    c16 = {v: i for i, v in enumerate(g16)}
    basis = [mul(pw(theta, i), pw(beta, j)) for i in range(3) for j in range(4)]   # index 4 i + j
    return beta, theta, g16, c16, basis


def main():
    beta, theta, g16, c16, basis = setup()
    # coordinates of z in the basis
    piv = {}
    for i, v in enumerate(basis):
        c = 1 << i
        for pb in sorted(piv, reverse=True):
            if (v >> pb) & 1:
                v ^= piv[pb][0]
                c ^= piv[pb][1]
        assert v, 'basis dependent'
        piv[v.bit_length() - 1] = (v, c)

    def coords(z):
        w, c = z, 0
        for pb in sorted(piv, reverse=True):
            if (w >> pb) & 1:
                w ^= piv[pb][0]
                c ^= piv[pb][1]
        assert w == 0
        return [g16[(c >> (4 * i)) & 15] for i in range(3)]      # a, b, c as GF(16) elements inside GF(4096)
    bad = 0
    for z in range(1, 4096):
        a, b, c = coords(z)
        e2 = mul(a, a) ^ mul(b, b) ^ mul(c, c) ^ mul(b, c)
        K = mul(mul(b, b), b) ^ mul(b, mul(c, c)) ^ mul(mul(c, c), c)
        N = mul(a, e2) ^ K
        lhs = tr(inv0(z))
        rhs = tr(mul(e2, inv0(N)), 12)        # e2 / N lies in GF(16); its trace to GF(2) taken in GF(4096)
        # Tr_{4096/2} of a GF(16) element = 3 x Tr_{16/2} = Tr_{16/2} (odd degree): the same value
        if lhs != rhs or N == 0:
            bad += 1
    print(json.dumps(dict(identity_wrong=bad, checked=4095)))
    # the joint map's 8 quadratic targets over the 8 bits (b bits 0-3, c bits 4-7), coordinates in (1, beta, ..)
    targets = []
    for out in range(8):
        T = set()
        lin = 0
        for u in range(8):
            for v in range(u, 8):
                pass
        targets.append(None)

    def joint(bb, cc):
        b_, c_ = g16[bb], g16[cc]
        return c16[mul(b_, c_)] | (c16[mul(mul(b_, b_), b_) ^ mul(b_, mul(c_, c_)) ^ mul(mul(c_, c_), c_)] << 4)
    # algebraic normal form of each output bit (degree <= 2 expected)
    anf = []
    for out in range(8):
        f = [(joint(x & 15, x >> 4) >> out) & 1 for x in range(256)]
        coef = f[:]
        for i in range(8):
            for x in range(256):
                if (x >> i) & 1:
                    coef[x] ^= coef[x ^ (1 << i)]
        mons = [x for x in range(256) if coef[x]]
        assert all(bin(x).count('1') <= 2 for x in mons), ('degree > 2', out)
        anf.append(mons)
    json.dump(dict(anf=anf), open(sys.argv[1], 'w')) if len(sys.argv) > 1 else None
    print(json.dumps(dict(quad_terms=[sum(1 for x in m if bin(x).count('1') == 2) for m in anf],
                          lin_terms=[sum(1 for x in m if bin(x).count('1') == 1) for m in anf])))


if __name__ == '__main__':
    main()
