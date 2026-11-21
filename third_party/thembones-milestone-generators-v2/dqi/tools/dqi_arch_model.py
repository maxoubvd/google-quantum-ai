#!/usr/bin/env python3
"""dqi_arch_model.py -- classical model for the DQI architecture workstream (numpy only).

Checks, on random errors of weight 18..35 drawn like the official verifier, the algebra used in
research/ws_dqi_arch.md:
  C0  key equation: Berlekamp-Massey on the 70 syndromes gives the locator of the sampled error
  C1  A_j = sum_i B_(i+j+1) S_i is the numerator of sum_k Y'_k/(z - X_k), Y'_k = Y_k X_k
  C2  baseline scan: C = A/B' mod B, phase = sum_j [B(x_j)=0] <M_j, C(x_j)/x_j> equals the verifier phase
  C3  B = E^2 + z O^2, B' = O^2
  C4  D = A_e O + A_o E:  D^2 = (A B)',  D(X_k) = sqrt(Y'_k) O(X_k)^2
  C5  c = interpolant of sqrt(Y'_k): c^2 = C mod B;  N = O C_e + E C_o: N^2 = A mod B, N(X_k) = c(X_k) O(X_k)
  C6  design W (Bezout pair, bilinear scan): ca E + cb O = 1, Psi = cb^2 + z ca^2 = 1/B' mod B,
      sqrt(Y'_k) = a(X_k) (cb(X_k) + sqrt(X_k) ca(X_k)); phase by this scan equals the verifier phase
  C7  T = E/O mod B takes the classical value sqrt(X_k) at every root; T^2 + z = Psi B exactly
  C8  three more free phase terms: s_0, s_1, s_2 (degree <= 2 in the words of A and B) equal E_71, E_72, E_73
  C9  in-place triangular Toeplitz map by middle-product recursion (the kit's tri): T(35) = 147 products,
      and reverse(tri(reversed(B[1:]), X)) = C for X = Laurent coefficients of C/B
  C10 design R (reciprocal scan, no convert): e_j = (a(x_j)/O(x_j))^2 / x_j at the roots; V = D/A mod B = 1/c
Usage: python dqi_arch_model.py [instances_per_weight=3] [seed=1]
"""
import sys
import numpy as np

MOD = 4179          # x^12 + x^6 + x^4 + x + 1, km.GF2Field(12).modulus
ORD = 4095
EXP = np.zeros(2 * ORD, dtype=np.int64)
LOG = np.zeros(4096, dtype=np.int64)
_x = 1
for _i in range(ORD):
    EXP[_i] = _x
    LOG[_x] = _i
    _x <<= 1
    if _x & 4096:
        _x ^= MOD
EXP[ORD:] = EXP[:ORD]
I64 = np.int64


def mul(a, b):
    a = np.asarray(a, dtype=I64)
    b = np.asarray(b, dtype=I64)
    return np.where((a == 0) | (b == 0), 0, EXP[(LOG[a] + LOG[b]) % ORD])


def inv(a):
    a = np.asarray(a, dtype=I64)
    assert np.all(a != 0), 'inverse of zero'
    return EXP[(ORD - LOG[a]) % ORD]


def fpow(a, e):
    a = np.asarray(a, dtype=I64)
    return np.where(a == 0, 0, EXP[(LOG[a] * (e % ORD)) % ORD])


def fsqrt(a):
    return fpow(a, 2048)


def fsq(a):
    return fpow(a, 2)


def sm(a, b):
    return int(mul(a, b))


def trim(p):
    p = np.asarray(p, dtype=I64)
    n = len(p)
    while n > 1 and p[n - 1] == 0:
        n -= 1
    return p[:n].copy()


def deg(p):
    p = trim(p)
    return -1 if (len(p) == 1 and p[0] == 0) else len(p) - 1


def padd(a, b):
    r = np.zeros(max(len(a), len(b)), dtype=I64)
    r[:len(a)] ^= np.asarray(a, dtype=I64)
    r[:len(b)] ^= np.asarray(b, dtype=I64)
    return trim(r)


def pmul(a, b):
    a = trim(a)
    b = trim(b)
    r = np.zeros(len(a) + len(b) - 1, dtype=I64)
    for i, ai in enumerate(a):
        if ai:
            r[i:i + len(b)] ^= mul(ai, b)
    return trim(r)


def pdivmod(a, b):
    a = trim(a).copy()
    b = trim(b)
    db = len(b) - 1
    assert deg(b) >= 0
    if len(a) - 1 < db:
        return np.zeros(1, dtype=I64), a
    q = np.zeros(len(a) - db, dtype=I64)
    ilc = inv(b[-1])
    for i in range(len(a) - 1, db - 1, -1):
        if a[i]:
            c = mul(a[i], ilc)
            q[i - db] = c
            a[i - db:i + 1] ^= mul(c, b)
    return trim(q), trim(a[:max(db, 1)])


def pmod(a, b):
    return pdivmod(a, b)[1]


def peval(p, xs):
    xs = np.asarray(xs, dtype=I64)
    r = np.zeros(xs.shape, dtype=I64)
    for c in np.asarray(p, dtype=I64)[::-1]:
        r = mul(r, xs) ^ c
    return r


def pderiv(p):
    p = trim(p)
    r = np.zeros(max(len(p) - 1, 1), dtype=I64)
    for i in range(1, len(p), 2):
        r[i - 1] = p[i]
    return trim(r)


def psquare(p):
    p = trim(p)
    r = np.zeros(2 * len(p) - 1, dtype=I64)
    r[0::2] = fsq(p)
    return trim(r)


def pshift(p, k=1):
    return trim(np.concatenate([np.zeros(k, dtype=I64), np.asarray(p, dtype=I64)]))


def split(p):
    """p = pe^2 + z po^2 (coefficients of pe, po are square roots of those of p)."""
    p = trim(p)
    pe = fsqrt(p[0::2])
    po = fsqrt(p[1::2]) if len(p) > 1 else np.zeros(1, dtype=I64)
    return trim(pe), trim(po)


def xgcd(a, b):
    """g, s, t with s a + t b = g, g monic."""
    r0, r1 = trim(a), trim(b)
    s0, s1 = np.array([1], dtype=I64), np.array([0], dtype=I64)
    t0, t1 = np.array([0], dtype=I64), np.array([1], dtype=I64)
    while deg(r1) >= 0:
        q, r = pdivmod(r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, padd(s0, pmul(q, s1))
        t0, t1 = t1, padd(t0, pmul(q, t1))
    c = inv(r0[-1])
    return trim(mul(c, r0)), trim(mul(c, s0)), trim(mul(c, t0))


def pinvmod(a, m):
    g, s, _ = xgcd(pmod(a, m), m)
    assert deg(g) == 0, 'not invertible'
    return pmod(s, m)


def from_roots(xs):
    p = np.array([1], dtype=I64)
    for x in xs:
        p = pmul(p, np.array([x, 1], dtype=I64))
    return p


def numerator(xs, ws, B):
    """sum_k ws[k] B(z)/(z - xs[k])."""
    r = np.zeros(1, dtype=I64)
    for x, w in zip(xs, ws):
        q, rem = pdivmod(B, np.array([x, 1], dtype=I64))
        assert deg(rem) == -1
        r = padd(r, mul(w, q))
    return r


def interp(xs, ys, B):
    dB = pderiv(B)
    return numerator(xs, mul(ys, inv(peval(dB, xs))), B)


def peq(a, b):
    a, b = trim(a), trim(b)
    return len(a) == len(b) and bool(np.all(a == b))


def berlekamp_massey(S):
    C, Bp, L, m, b = [1], [1], 0, 1, 1
    for n in range(len(S)):
        d = int(S[n])
        for i in range(1, L + 1):
            d ^= sm(C[i], S[n - i])
        if d == 0:
            m += 1
            continue
        coef = sm(d, int(inv(b)))
        T = C[:]
        C = C + [0] * max(0, len(Bp) + m - len(C))
        for i, v in enumerate(Bp):
            C[i + m] ^= sm(coef, v)
        if 2 * L <= n:
            L, Bp, b, m = n + 1 - L, T, d, 1
        else:
            m += 1
    return C[:L + 1], L


# --- in-place triangular Toeplitz map by middle-product recursion (transcribed from the kit's tri) ---
_VW = {}


def vandermonde(n):
    if n not in _VW:
        V = [[int(fpow(t, j)) if (t or j) else 1 for j in range(n)] for t in range(n)]
        for t in range(n):
            V[t][0] = 1
        M = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(V)]
        for col in range(n):
            piv = next(r for r in range(col, n) if M[r][col])
            M[col], M[piv] = M[piv], M[col]
            ic = int(inv(M[col][col]))
            M[col] = [sm(ic, v) for v in M[col]]
            for r in range(n):
                if r != col and M[r][col]:
                    f = M[r][col]
                    M[r] = [v ^ sm(f, w) for v, w in zip(M[r], M[col])]
        _VW[n] = (V, [row[n:] for row in M])
    return _VW[n]


def tri(q, x, counter):
    """x <- L(q) x in place, L unit lower triangular Toeplitz: x_i += sum_{j<i} q_{i-j} x_j, q = (q_1..q_{n-1})."""
    n = len(x)
    if n == 1:
        return
    k = n // 2
    m = n - k
    hi = x[k:]
    tri(q[:m - 1], hi, counter)
    x[k:] = hi
    npts = k + m - 1
    V, W = vandermonde(npts)
    for t in range(npts):
        ua = 0
        for i in range(k):
            ua ^= sm(V[t][i], x[k - 1 - i])
        cb = 0
        for s in range(npts):
            cb ^= sm(W[s][t], q[s])
        pr = sm(ua, cb)
        counter[0] += 1
        for i2 in range(m):
            x[k + i2] ^= sm(V[t][i2], pr)
    lo = x[:k]
    tri(q[:k - 1], lo, counter)
    x[:k] = lo


def tri_direct(q, x):
    n = len(x)
    y = list(x)
    for i in range(n):
        for j in range(i):
            y[i] ^= sm(q[i - j - 1], x[j])
    return y


def weight_law():
    m, l, b, r = 4095, 35, 12, 2016
    q = 1 << b
    k = np.arange(l + 1, dtype=np.float64)
    off = np.sqrt(k[1:] * (m - k[1:] + 1))
    d = (q - 2 * r) / np.sqrt(r * (q - r))
    _, vecs = np.linalg.eigh(np.diag(k * d) + np.diag(off, 1) + np.diag(off, -1))
    w = np.abs(vecs[:, -1])
    return (w / np.linalg.norm(w)) ** 2


def run_instance(w, rng, stats):
    pos = rng.choice(ORD, size=w, replace=False)
    vals = rng.integers(1, ORD, size=w, endpoint=True, dtype=I64)
    masks = rng.integers(0, ORD, size=ORD, endpoint=True, dtype=I64)
    X = EXP[pos]
    Y = vals
    S = np.array([int(np.bitwise_xor.reduce(mul(Y, fpow(X, i + 1)))) for i in range(70)], dtype=I64)
    err = np.zeros(ORD, dtype=I64)
    err[pos] = vals
    phase_true = int(sum(int(v).bit_count() for v in (err & masks))) & 1
    xs_all = EXP[:ORD]
    ok = {}

    # C0
    Cc, L = berlekamp_massey(S)
    B = from_roots(X)
    sigma_rev = np.array([Cc[w - j] if w - j < len(Cc) else 0 for j in range(w + 1)], dtype=I64) if L == w else None
    ok['C0'] = (L == w) and peq(sigma_rev, B)

    # C1
    Yp = mul(Y, X)
    A = numerator(X, Yp, B)
    A2 = np.zeros(w, dtype=I64)
    for j in range(w):
        acc = 0
        for i in range(0, w - j):
            acc ^= sm(B[i + j + 1], S[i])
        A2[j] = acc
    ok['C1'] = peq(A, A2)

    # C2 baseline scan
    dB = pderiv(B)
    Cpol = pmod(pmul(A, pinvmod(dB, B)), B)
    flags = peval(B, xs_all) == 0
    e_scan = np.where(flags, mul(peval(Cpol, xs_all), inv(xs_all)), 0)
    ph = int(sum(int(v).bit_count() for v in (e_scan & masks))) & 1
    ok['C2'] = bool(np.all(e_scan == err)) and ph == phase_true and bool(np.all(peval(Cpol, X) == Yp))

    # C3
    E, O = split(B)
    ok['C3'] = peq(padd(psquare(E), pshift(psquare(O))), B) and peq(dB, psquare(O))

    # C4
    Ae, Ao = split(A)
    D = padd(pmul(Ae, O), pmul(Ao, E))
    sY = fsqrt(Yp)
    OX = peval(O, X)
    ok['C4'] = peq(psquare(D), pderiv(pmul(A, B))) and bool(np.all(peval(D, X) == mul(sY, fsq(OX))))

    # C5
    c = interp(X, sY, B)
    Ce, Co = split(Cpol)
    Nn = padd(pmul(O, Ce), pmul(E, Co))
    ok['C5'] = (peq(pmod(psquare(c), B), Cpol) and deg(Nn) <= w - 1 and peq(pmod(psquare(Nn), B), A)
                and bool(np.all(peval(Nn, X) == mul(peval(c, X), OX))))

    # C6 design W
    g, ca, cb = xgcd(E, O)
    Psi = padd(psquare(cb), pshift(psquare(ca)))
    wall = fsqrt(xs_all)
    a_all = peval(Ae, xs_all) ^ mul(wall, peval(Ao, xs_all))
    om_all = peval(cb, xs_all) ^ mul(wall, peval(ca, xs_all))
    e_w = np.where(flags, mul(fsq(mul(a_all, om_all)), inv(xs_all)), 0)
    ph_w = int(sum(int(v).bit_count() for v in (e_w & masks))) & 1
    ok['C6'] = (deg(g) == 0 and peq(pmod(pmul(Psi, dB), B), np.array([1])) and bool(np.all(e_w == err))
                and ph_w == phase_true and deg(ca) < max(deg(O), 1) and deg(cb) < max(deg(E), 1))

    # C7
    T = pmod(pmul(E, pinvmod(O, B)), B)
    qT, rT = pdivmod(padd(psquare(T), np.array([0, 1], dtype=I64)), B)
    ok['C7'] = bool(np.all(peval(T, X) == fsqrt(X))) and deg(rT) == -1 and peq(qT, pmod(Psi, B))

    # C8
    Y71 = mul(Y, fpow(X, 71))
    A71 = np.zeros(w, dtype=I64)
    t71 = numerator(X, Y71, B)
    A71[:len(t71)] = t71
    Ei = [int(np.bitwise_xor.reduce(mul(Y, fpow(X, i)))) for i in (71, 72, 73)]
    s0 = int(A71[w - 1])
    s1 = int(A71[w - 2]) ^ sm(B[w - 1], A71[w - 1])
    s2 = int(A71[w - 3]) ^ sm(B[w - 1], A71[w - 2]) ^ sm(int(fsq(B[w - 1])) ^ int(B[w - 2]), A71[w - 1])
    ok['C8'] = [s0, s1, s2] == Ei

    # C9
    u = mul(peval(Cpol, X), inv(peval(dB, X)))
    Xl = [int(np.bitwise_xor.reduce(mul(u, fpow(X, j)))) for j in range(w)]
    qv = [int(v) for v in B[1:w][::-1]]
    cnt = [0]
    y = list(Xl)
    tri(qv, y, cnt)
    Cfull = np.zeros(w, dtype=I64)
    Cfull[:len(Cpol)] = Cpol
    ok['C9'] = (y == tri_direct(qv, Xl)) and bool(np.all(np.array(y[::-1], dtype=I64) == Cfull))
    stats.setdefault('tri', {})[w] = cnt[0]

    # C10 design R
    O_all = peval(O, xs_all)
    safe = np.where(O_all == 0, 1, O_all)
    e_r = np.where(flags, mul(fsq(mul(a_all, inv(safe))), inv(xs_all)), 0)
    V = pmod(pmul(D, pinvmod(A, B)), B) if deg(xgcd(A, B)[0]) == 0 else None
    ok['C10'] = bool(np.all(e_r == err)) and (V is None or bool(np.all(mul(peval(V, X), peval(c, X)) == 1)))
    if V is None:
        stats['A_not_coprime'] = stats.get('A_not_coprime', 0) + 1
    return ok


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    rng = np.random.default_rng(seed)
    p = weight_law()
    print('weight law: mean %.3f  P(35) %.4f  P(<=25) %.2e  P(<=22) %.2e  P(<=20) %.2e  P(<=17) %.2e'
          % (float(np.arange(36) @ p), p[35], p[:26].sum(), p[:23].sum(), p[:21].sum(), p[:18].sum()))
    names = ['C%d' % i for i in range(11)]
    tot = {n: 0 for n in names}
    fail = {n: [] for n in names}
    stats = {}
    count = 0
    for w in range(18, 36):
        for _ in range(reps):
            ok = run_instance(w, rng, stats)
            count += 1
            for n in names:
                tot[n] += bool(ok[n])
                if not ok[n]:
                    fail[n].append(w)
    print('instances: %d (weights 18..35, %d each)' % (count, reps))
    for n in names:
        print('%-4s %d/%d %s' % (n, tot[n], count, '' if tot[n] == count else 'FAILED at weights %s' % fail[n]))
    print('tri product counts by size:', stats['tri'])
    print('A not coprime to B (V check skipped):', stats.get('A_not_coprime', 0))
    print('ALL PASS' if all(tot[n] == count for n in names) else 'SOME FAILED')


if __name__ == '__main__':
    main()
