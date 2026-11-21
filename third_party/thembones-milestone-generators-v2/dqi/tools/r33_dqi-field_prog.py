#!/usr/bin/env python3
"""r33_dqi-field_prog.py -- the 12-bit program "one word times R scalars" built from a per-point GF(8) program
of r33_dqi-field_quad.py (rule yfree1 or yfree), with its exact check (seat DQI-FIELD, round R33).

Field: GF(2^12) with the competition modulus x^12 + x^6 + x^4 + x + 1 (4179).  K = GF(8) inside it; an element
is a cubic over K in theta; a product is the product of two cubics over K: 7 conditions = 6 finite points of K
and the point at infinity (the kit's 42 = 7 x 6 scheme).  At each point the word's value X (3 bits) meets the
R scalars' values Y_1..Y_R (3 bits each) in the SAT program (L gates, mixed forms):
    X * Y_r = sum_i W[r][i] g_i  +  linear terms  +  Q_r(Y)        (Q_r: quadratic in the scalars alone)
The Q_r part depends on the scalars only: it is paid ONCE per scalar tuple (the "cache": one AND per independent
garbage form) and XORed into every target word by CNOT fan-out; nothing is held on wires.

Usage: r33_dqi-field_prog.py PROGRAM.json [trials=20000] [seed=1]
Prints: gates per word, cache ANDs per scalar tuple, and the exact check (random triples and all unit vectors).
"""
import json
import random
import sys

MOD = 4179
K3 = 0b1011  # x^3 + x + 1, the SAT script's GF(8)


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


def power(a, e):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a)
        a = mul(a, a)
        e >>= 1
    return r


def inv(a):
    return power(a, 4094)


def solve_gf2(rows, nvars):
    """rows: list of (mask over nvars, rhs int payload); Gauss-Jordan; returns pivot dict col -> (mask, payload)."""
    piv = {}
    for mask, pay in rows:
        for c, (pm, pp) in piv.items():
            if (mask >> c) & 1:
                mask ^= pm
                pay ^= pp
        if mask:
            c = mask.bit_length() - 1
            for c2 in list(piv):
                pm, pp = piv[c2]
                if (pm >> c) & 1:
                    piv[c2] = (pm ^ mask, pp ^ pay)
            piv[c] = (mask, pay)
    return piv


def build_scheme():
    """Returns E (list over 7 points of 3 masks: bit j of the point value as a parity mask on the 12 input bits)
    and J (list over 7 points of 3 field elements: the product contribution of bit j of that point's product)."""
    sub = [v for v in range(4096) if power(v, 8) == v]
    assert len(sub) == 8
    om = next(v for v in sub if v > 1 and (mul(mul(v, v), v) ^ v ^ 1) == 0)   # omega^3 = omega + 1
    kb = [1, om, mul(om, om)]

    def kel(bits):   # 3 bits -> element of K inside L
        r = 0
        for j in range(3):
            if (bits >> j) & 1:
                r ^= kb[j]
        return r
    kelem = [kel(b) for b in range(8)]
    kbits = {kelem[b]: b for b in range(8)}
    # theta: an element of degree 4 over K
    theta = None
    for cand in range(2, 4096):
        basis = [mul(kb[j], power(cand, d)) for d in range(4) for j in range(3)]
        if len(solve_gf2([(b, 0) for b in basis], 12)) == 12:
            theta = cand
            break
    assert theta is not None
    tp = [power(theta, d) for d in range(7)]
    basis = [mul(kb[j], tp[d]) for d in range(4) for j in range(3)]   # index 3 d + j
    # coordinates: for each input bit e (unit vector 1 << e) its 12 coordinates in the basis
    piv = solve_gf2([(basis[i], 1 << i) for i in range(12)], 12)     # express unit vectors by the basis
    coord_of_unit = []
    for e in range(12):
        # reduce the unit vector by the pivots to read its combination
        mask, pay = 1 << e, 0
        for c in sorted(piv, reverse=True):
            pm, pp = piv[c]
            if (mask >> c) & 1:
                mask ^= pm
                pay ^= pp
        assert mask == 0
        coord_of_unit.append(pay)        # bit (3 d + j) set: K-coordinate j of the theta^d digit

    def digits(x):                       # x -> list of 4 K-elements (as 3-bit values)
        c = 0
        for e in range(12):
            if (x >> e) & 1:
                c ^= coord_of_unit[e]
        return [(c >> (3 * d)) & 7 for d in range(4)]

    def kmul_bits(a, b):
        return kbits[mul(kelem[a], kelem[b])]
    pts = [0, 1, 2, 3, 4, 5]             # six finite points of K (3-bit labels), plus infinity
    # value at a point as a function of the input: linear over GF(2); build masks on the 12 input bits
    E = []
    for t in pts + ['inf']:
        masks = [0, 0, 0]
        for e in range(12):
            dg = digits(1 << e)
            if t == 'inf':
                val = dg[3]
            else:
                val, tpow = 0, 1
                for d in range(4):
                    val ^= kmul_bits(dg[d], tpow)
                    tpow = kmul_bits(tpow, t)
            for j in range(3):
                if (val >> j) & 1:
                    masks[j] |= 1 << e
        E.append(masks)
    # interpolation: c(t) = P_t for finite t, c_6 = P_inf; solve c_0..c_5 over K (work inside L)
    J = []
    V = [[power(kelem[t], d) if (t or d) else 1 for d in range(6)] for t in pts]
    for t in pts:
        V[pts.index(t)][0] = 1
    # invert V over K by Gauss-Jordan in L arithmetic
    n = 6
    M = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(V)]
    for col in range(n):
        pr = next(r for r in range(col, n) if M[r][col])
        M[col], M[pr] = M[pr], M[col]
        ic = inv(M[col][col])
        M[col] = [mul(ic, v) for v in M[col]]
        for r in range(n):
            if r != col and M[r][col]:
                f = M[r][col]
                M[r] = [v ^ mul(f, w) for v, w in zip(M[r], M[col])]
    Vinv = [row[n:] for row in M]        # c_d = sum_t Vinv[d][t] * (P_t - c_6 t^6)
    for ti, t in enumerate(pts + ['inf']):
        contrib = []
        for j in range(3):
            unit = kb[j]
            if t == 'inf':
                # P_inf = unit: c_6 = unit; finite right-hand sides become - c_6 t^6
                res = mul(unit, tp[6])
                for d in range(6):
                    cd = 0
                    for si, s in enumerate(pts):
                        cd ^= mul(Vinv[d][si], mul(unit, power(kelem[s], 6)) if s else 0)
                    res ^= mul(cd, tp[d])
            else:
                res = 0
                for d in range(6):
                    res ^= mul(mul(Vinv[d][ti], unit), tp[d])
            contrib.append(res)
        J.append(contrib)
    return E, J


def parity(x):
    return bin(x).count('1') & 1


def main():
    prog = json.load(open(sys.argv[1]))
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    assert prog['k'] == 3 and prog['result'] == 'SAT'
    R, L = prog['n'], prog['L']
    nv = 3 + 3 * R
    l = [[int(ch) for ch in row] for row in prog['program']['l']]
    lp = [[int(ch) for ch in row] for row in prog['program']['lp']]
    W = [[int(ch) for ch in row] for row in prog['program']['W']]
    nout = 3 * R
    # per-output corrections of the point program
    lin = [[0] * nv for _ in range(nout)]
    quad = [dict() for _ in range(nout)]
    for c in range(nout):
        for i in range(L):
            if not W[c][i]:
                continue
            for a in range(nv):
                if l[i][a] & lp[i][a]:
                    lin[c][a] ^= 1
                for b in range(a + 1, nv):
                    if (l[i][a] & lp[i][b]) ^ (l[i][b] & lp[i][a]):
                        quad[c][(a, b)] = quad[c].get((a, b), 0) ^ 1
    # the quadratic part must be: cross terms = the GF(8) products, x-x nothing, y-y = garbage
    garbage = [sorted(p for p, v in quad[c].items() if v and p[0] >= 3) for c in range(nout)]
    xx = [p for c in range(nout) for p, v in quad[c].items() if v and p[1] < 3]
    assert not xx, 'word-side garbage does not cancel'
    cross_scalar = sorted({p for g in garbage for p in g if (p[0] - 3) // 3 != (p[1] - 3) // 3})
    # cache: per scalar, the span of its garbage forms (inside Lambda^2 of 3 bits every form is one AND)
    cache = 0
    per_scalar = []
    for s in range(R):
        idx = {(3 + 3 * s + a, 3 + 3 * s + b): n for n, (a, b) in enumerate(((0, 1), (0, 2), (1, 2)))}
        vecs = []
        for c in range(nout):
            v = 0
            for p in garbage[c]:
                if p in idx:
                    v |= 1 << idx[p]
            vecs.append(v)
        d = len(solve_gf2([(v, 0) for v in vecs if v], 3))
        per_scalar.append(d)
        cache += d
    ncross = len(cross_scalar)
    E, J = build_scheme()
    # sanity of the scheme: plain products by the 7 point products
    def kmul3(a, b):
        r = 0
        for i in range(3):
            if (b >> i) & 1:
                r ^= a << i
        for i in (4, 3):
            if (r >> i) & 1:
                r ^= K3 << (i - 3)
        return r

    def point_values(x):
        return [sum(parity(x & E[p][j]) << j for j in range(3)) for p in range(7)]

    def inject(p, bits):
        r = 0
        for j in range(3):
            if (bits >> j) & 1:
                r ^= J[p][j]
        return r

    def run(g, ds):
        """the program: returns the R products and the number of ANDs spent on the word (gates)."""
        X = point_values(g)
        Ys = [point_values(d) for d in ds]
        out = [0] * R
        gates = 0
        for p in range(7):
            v = [(X[p] >> j) & 1 for j in range(3)]
            for r in range(R):
                v += [(Ys[r][p] >> j) & 1 for j in range(3)]
            gv = []
            for i in range(L):
                a = 0
                b = 0
                for q in range(nv):
                    a ^= l[i][q] & v[q]
                    b ^= lp[i][q] & v[q]
                gv.append(a & b)
                gates += 1
            for r in range(R):
                bits = 0
                for j in range(3):
                    c = 3 * r + j
                    val = 0
                    for i in range(L):
                        val ^= W[c][i] & gv[i]
                    for q in range(nv):
                        val ^= lin[c][q] & v[q]
                    for (a, b) in garbage[c]:          # the scalar-only part (from the cache)
                        val ^= v[a] & v[b]
                    bits |= val << j
                out[r] ^= inject(p, bits)
        return out, gates
    # scheme sanity
    for _ in range(2000):
        x, y = rng.randrange(4096), rng.randrange(4096)
        X, Y = point_values(x), point_values(y)
        z = 0
        for p in range(7):
            z ^= inject(p, kmul3(X[p], Y[p]))
        assert z == mul(x, y), 'evaluation scheme wrong'
    bad = 0
    tests = [(1 << a, [1 << b] * R) for a in range(12) for b in range(12)]
    tests += [(rng.randrange(4096), [rng.randrange(4096) for _ in range(R)]) for _ in range(trials)]
    gates = None
    for g, ds in tests:
        out, gates = run(g, ds)
        if any(out[r] != mul(g, ds[r]) for r in range(R)):
            bad += 1
    res = {'R': R, 'gates_per_point': L, 'gates_per_word': gates, 'bilinear_per_word': 42 * R,
           'cache_forms_per_scalar_per_point': per_scalar, 'cache_ands_per_scalar_tuple': 7 * cache,
           'cross_scalar_garbage_monomials_per_point': ncross,
           'tests': len(tests), 'wrong': bad, 'verdict': 'PASS' if bad == 0 else 'FAIL'}
    print(json.dumps(res))


if __name__ == '__main__':
    main()
