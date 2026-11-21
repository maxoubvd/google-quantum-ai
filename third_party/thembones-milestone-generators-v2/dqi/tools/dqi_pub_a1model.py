#!/usr/bin/env python3
"""dqi_pub_a1model.py -- builder D1, item a1: classical check of the CRT body for our un-synthesis (numpy-free).

The kit's body (quotient_add): with S, T in q-adic digit words (qadic_matrix(m, k), q = x^2 + x + k), l the head
after the division (l = U / W, W = S mod q, U the residue before it), it adds Qt . floor(l S_c / q) to T.
The CRT body adds the same polynomial F = (l S_c - U) / q (exact) as
    premap(U)  : T += Qt . sum_j  U(t_j) / q(t_j) . L_j        (free, before the division)
    updates    : T += Qt . sum_j  l(t_j) S_c(t_j) / q(t_j) . L_j   (len(T) rank-one products, after it)
with L_j the Lagrange basis at t_j = 0..len(T)-1. This script builds those forms and checks them against the
floor map on random inputs. Usage: python dqi_pub_a1model.py [trials=200] [seed=1]
"""
import random
import sys

MOD = 0x1053
EXP = [0] * 8190
LOG = [0] * 4096
x = 1
for i in range(4095):
    EXP[i] = x
    LOG[x] = i
    x <<= 1
    if x & 4096:
        x ^= MOD
for i in range(4095, 8190):
    EXP[i] = EXP[i - 4095]


def mul(a, b):
    if a == 0 or b == 0:
        return 0
    return EXP[LOG[a] + LOG[b]]


def inv(a):
    return EXP[(4095 - LOG[a]) % 4095]


def power(a, n):
    r = 1
    for _ in range(n):
        r = mul(r, a)
    return r


def divmodp(a, b):
    """Polynomial division over GF(2^12), coefficients ascending."""
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    if len(a) < len(b):
        return [0], a
    out = [0] * (len(a) - len(b) + 1)
    bi = inv(b[-1])
    for k in range(len(a) - len(b), -1, -1):
        t = mul(a[k + len(b) - 1], bi)
        out[k] = t
        for j, v in enumerate(b):
            a[k + j] ^= mul(t, v)
    return out, a[:len(b) - 1]


def qadic_matrix(n, k):
    """The kit's paired_coefficient.qadic_matrix: coefficients -> q-adic digits (rows = digits)."""
    cols = []
    for j in range(n):
        p = [0] * j + [1]
        out = []
        for _ in range((n + 1) // 2):
            p, r = divmodp(p, [k, 1, 1])
            out += (r + [0, 0])[:2]
        cols.append(out[:n])
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def matvec(M, v):
    out = []
    for row in M:
        s = 0
        for a, b in zip(row, v):
            s ^= mul(a, b)
        out.append(s)
    return out


def matinv(M):
    n = len(M)
    A = [list(r) + [int(i == j) for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c])
        A[c], A[p] = A[p], A[c]
        iv = inv(A[c][c])
        A[c] = [mul(v, iv) for v in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                f = A[r][c]
                A[r] = [a ^ mul(f, b) for a, b in zip(A[r], A[c])]
    return [r[n:] for r in A]


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x_ in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] ^= mul(x_, y)
    return out


def peval(p, t):
    r = 0
    for c in reversed(p):
        r = mul(r, t) ^ c
    return r


def crt_forms(m, n, k):
    """Forms for a body with source length m and target length n (n = m - 1).
    Returns (pre, updates): pre[i][w] = coefficient of U[w] in T[i]; updates = list of (ac, bc, tc)."""
    assert n == m - 1
    Qs, Qt = qadic_matrix(m, k), qadic_matrix(n, k)
    Qs_inv = matinv(Qs)
    pts = list(range(n))
    V = [[power(t, j) for j in range(n)] for t in pts]
    W = matinv(V)                               # coefficients of L_j: column j of W
    q = [k, 1, 1]
    tcs = []
    for j, t in enumerate(pts):
        qt = peval(q, t)
        Lj = [W[i][j] for i in range(n)]
        tcs.append([mul(c, inv(qt)) for c in matvec(Qt, Lj)])
    pre = [[0, 0] for _ in range(n)]
    for j, t in enumerate(pts):
        for i in range(n):
            pre[i][0] ^= tcs[j][i]
            pre[i][1] ^= mul(tcs[j][i], t)
    updates = []
    for j, t in enumerate(pts):
        row = [power(t, i) for i in range(m)]           # S_c(t) as a form of the coefficients
        bc = [0] * m                                    # ... of the digit words: row . Qs^-1
        for i in range(m):
            for c in range(m):
                bc[c] ^= mul(row[i], Qs_inv[i][c])
        updates.append(([1, t], bc, tcs[j]))
    return pre, updates


def gf24_mul(a, b, k):
    """(a0 + a1 x)(b0 + b1 x) mod x^2 + x + k."""
    p = pmul(a, b) + [0, 0, 0]
    _, r = divmodp(p[:3], [k, 1, 1])
    return (r + [0, 0])[:2]


def check(trials, seed):
    rnd = random.Random(seed)
    bad = 0
    for trial in range(trials):
        n = rnd.choice([3, 5, 7, 9, 15, 33, 35]) - (trial % 2)   # targets odd and even
        m = n + 1
        c = 2 * rnd.randrange(17)
        k = 2048 ^ mul(c, c) ^ c
        S = [rnd.randrange(4096) for _ in range(m)]
        T = [rnd.randrange(4096) for _ in range(n)]
        l = [rnd.randrange(4096), rnd.randrange(1, 4096)]
        Qs, Qt = qadic_matrix(m, k), qadic_matrix(n, k)
        Sc = matvec(matinv(Qs), S)
        # the kit: T += Qt . floor(l Sc / q)
        fl, rem = divmodp(pmul(l, Sc), [k, 1, 1])
        fl = (fl + [0] * n)[:n]
        kit = [a ^ b for a, b in zip(T, matvec(Qt, fl))]
        # the residue before the division: U = l * (Sc mod q)
        _, Wr = divmodp(Sc, [k, 1, 1])
        U = gf24_mul(l, (Wr + [0, 0])[:2], k)
        assert U == (rem + [0, 0])[:2]
        pre, ups = crt_forms(m, n, k)
        mine = list(T)
        for i in range(n):
            mine[i] ^= mul(pre[i][0], U[0]) ^ mul(pre[i][1], U[1])
        for ac, bc, tc in ups:
            lv = mul(ac[0], l[0]) ^ mul(ac[1], l[1])
            sv = 0
            for a, s in zip(bc, S):
                sv ^= mul(a, s)
            pr = mul(lv, sv)
            for i in range(n):
                mine[i] ^= mul(tc[i], pr)
        bad += mine != kit
    return bad


if __name__ == '__main__':
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    print('mismatches', check(trials, seed), 'of', trials)
