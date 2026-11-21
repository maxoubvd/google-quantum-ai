#!/usr/bin/env python3
"""r81_dqi_fold_klift.py -- lift a GF(8) scheme of a bipartite product kernel K(p,n) to K = GF(2^12) (seat R81-DQI-FOLD).

PUBLIC MATHEMATICS: K = GF(8)[y] / f(y), deg f = 4. A product of two elements of K is a product of two polynomials of
degree <= 3 over GF(8), fixed by its values at 7 points of P^1(GF(8)) (Toom-Cook); every map around the point values
(evaluation, interpolation, reduction mod f) is GF(2)-linear. So a bilinear scheme of r terms for the GF(8) kernel
"p heads x n words, every product its own output word" gives one of 7 r terms for the same kernel over K: at each of
the 7 points the GF(8) scheme is applied to the point values of the heads and the words, and each output word of the
scheme is spread into the K-word of its product by that point's interpolation map.

usage: r81_dqi_fold_klift.py <GF(8) scheme of grid p,n (or 'outer' p = 2)> <p> <n> <out> [modulus hex = 1053]
       r81_dqi_fold_klift.py selftest     (lifts one 6-term GF(8) product: must be a 42-term K-product scheme)
       r81_dqi_fold_klift.py dots <K-level scheme of K(p,n)> <p> <n> <out> [modulus hex = 1053]
           the SAME scheme with its legs permuted (the K-product's trilinear form Tr(x y z) is symmetric): p dot
           products of ONE n-word vector x with p n-word vectors z_0 .. z_(p-1): out_k ^= sum_j x_j z_(k,j).
           Operand a = x (n words), operand b = (z_0, .., z_(p-1)) (p n words), outputs p words. Same term count.
           Checked as a tensor identity against the dots computed from their definition.
The output is checked here as a tensor identity against the kernel over K computed from its definition.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r81_dqi_fold_seed as fs  # noqa: E402

S8, P8 = 3, 0xB


def kmul(x, y, poly):
    return fs.gfmul(x, y, 12, poly)


def kpow(x, e, poly):
    r = 1
    while e:
        if e & 1:
            r = kmul(r, x, poly)
        x = kmul(x, x, poly)
        e >>= 1
    return r


def setup(poly):
    # GF(8) inside K: beta with beta^3 + beta + 1 = 0
    beta = next(b for b in range(2, 4096) if kpow(b, 8, poly) == b and (kpow(b, 3, poly) ^ b ^ 1) == 0)

    def phi(c):                       # GF(8) (basis 1, g, g^2; g^3 = g + 1) -> K
        r = 0
        pw = 1
        for i in range(3):
            if c >> i & 1:
                r ^= pw
            pw = kmul(pw, beta, poly)
        return r
    y = 2
    ypow = [kpow(y, i, poly) for i in range(7)]
    # x in K <-> (x_0..x_3) in GF(8)^4, x = sum phi(x_i) y^i: columns of the 12 x 12 matrix over GF(2)
    cols = [kmul(phi(1 << k), ypow[i], poly) for i in range(4) for k in range(3)]   # column index 3 i + k
    # invert: for each K basis vector e_j find its GF(8)^4 coordinates
    rows = [(cols[c], 1 << c) for c in range(12)]
    piv = {}
    for v, comb in rows:
        while v:
            t = v.bit_length() - 1
            if t in piv:
                v ^= piv[t][0]
                comb ^= piv[t][1]
            else:
                piv[t] = (v, comb)
                break
        else:
            raise SystemExit("1, y, y^2, y^3 is not a GF(8)-basis of K")

    def coords(x):                    # 12-bit packed (x_0..x_3), 3 bits each
        comb = 0
        while x:
            t = x.bit_length() - 1
            x ^= piv[t][0]
            comb ^= piv[t][1]
        return comb
    return phi, ypow, coords


def g8mul(a, b):
    return fs.gfmul(a, b, S8, P8)


def g8inv(a):
    return next(b for b in range(1, 8) if g8mul(a, b) == 1)


def interp_matrix(points):
    """W[k][t] over GF(8): z_k = sum_t W[k][t] Z(p_t) for deg Z <= 6 and the 7 points (a point may be 'inf')."""
    n = len(points)
    V = []
    for p in points:
        if p == "inf":
            V.append([0] * (n - 1) + [1])
        else:
            row, pw = [], 1
            for _ in range(n):
                row.append(pw)
                pw = g8mul(pw, p)
            V.append(row)
    # invert V over GF(8) (Gauss-Jordan)
    A = [V[i][:] + [1 if i == j else 0 for j in range(n)] for i in range(n)]
    for c in range(n):
        r = next(i for i in range(c, n) if A[i][c])
        A[c], A[r] = A[r], A[c]
        iv = g8inv(A[c][c])
        A[c] = [g8mul(x, iv) for x in A[c]]
        for i in range(n):
            if i != c and A[i][c]:
                f = A[i][c]
                A[i] = [x ^ g8mul(f, yv) for x, yv in zip(A[i], A[c])]
    return [row[n:] for row in A]


def lift(terms, p_, n_, poly, points=(0, 1, 2, 3, 4, 5, "inf")):
    phi, ypow, coords = setup(poly)
    W = interp_matrix(list(points))
    # ev_t: K (12 bits) -> GF(8) (3 bits), as the 3 masks over K bits of its output bits
    evmask = []
    for p in points:
        masks = [0, 0, 0]
        for j in range(12):
            c = coords(1 << j)
            xi = [c >> (3 * i) & 7 for i in range(4)]
            if p == "inf":
                v = xi[3]
            else:
                v, pw = 0, 1
                for i in range(4):
                    v ^= g8mul(xi[i], pw)
                    pw = g8mul(pw, p)
            for b in range(3):
                if v >> b & 1:
                    masks[b] |= 1 << j
        evmask.append(masks)
    # interp_t: GF(8) value w at point t -> K element sum_k phi(W[k][t] w) y^k
    def interp(t, w):
        r = 0
        for k in range(7):
            c = g8mul(W[k][t], w)
            if c:
                r ^= kmul(phi(c), ypow[k], poly)
        return r

    def form(mask3, t):               # GF(8) form (3-bit mask) composed with ev_t: a 12-bit mask over K bits
        r = 0
        for b in range(3):
            if mask3 >> b & 1:
                r ^= evmask[t][b]
        return r
    out = []
    for t in range(len(points)):
        for a, b, c in terms:
            am = 0
            for i in range(p_):
                am |= form(a >> (3 * i) & 7, t) << (12 * i)
            bm = 0
            for j in range(n_):
                bm |= form(b >> (3 * j) & 7, t) << (12 * j)
            cm = 0
            for w in range(p_ * n_):
                g = c >> (3 * w) & 7
                if g:
                    cm |= interp(t, g) << (12 * w)
            out.append((am, bm, cm))
    return out


def trace_tables(poly):
    def tr(x):
        t, y = 0, x
        for _ in range(12):
            t ^= y
            y = kmul(y, y, poly)
        assert t in (0, 1)
        return t
    # vec2form[c] = mask m with bit i = Tr(c x^i); form2vec = its inverse
    v2f = {}
    for c in range(4096):
        m = 0
        for i in range(12):
            if tr(kmul(c, 1 << i, poly)):
                m |= 1 << i
        v2f[c] = m
    f2v = {m: c for c, m in v2f.items()}
    assert len(f2v) == 4096
    return v2f, f2v


def dots_tensor(p_, n_, poly):
    """T[x bit][z bit] = output mask of the kernel out_k ^= sum_j x_j z_(k,j)."""
    na, nb = 12 * n_, 12 * p_ * n_
    T = [[0] * nb for _ in range(na)]
    for j in range(n_):
        for k in range(p_):
            for i in range(12):
                for l in range(12):
                    T[12 * j + i][12 * (k * n_ + j) + l] = kmul(1 << i, 1 << l, poly) << (12 * k)
    return T


def main():
    if sys.argv[1] == "dots":
        src, p_, n_, out = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
        poly = int(sys.argv[6], 16) if len(sys.argv) > 6 else 0x1053
        na, nb, nc, terms = fs.read_scheme(src)
        assert (na, nb, nc) == (12 * p_, 12 * n_, 12 * p_ * n_)
        v2f, f2v = trace_tables(poly)
        new = []
        for a, b, c in terms:
            bm = 0
            for w in range(p_ * n_):
                bm |= v2f[c >> (12 * w) & 4095] << (12 * w)
            cm = 0
            for k in range(p_):
                cm |= f2v[a >> (12 * k) & 4095] << (12 * k)
            new.append((b, bm, cm))
        ok = fs.tensor_of(new, 12 * n_, 12 * p_ * n_) == dots_tensor(p_, n_, poly)
        print("%s -> %d dot products of one %d-word vector, %d terms: tensor identity %s (separate K-products: %d)"
              % (src, p_, n_, len(new), "PASS" if ok else "FAIL", 42 * p_ * n_))
        if ok:
            fs.write_scheme(out, new, 12 * n_, 12 * p_ * n_, 12 * p_)
            print("wrote", out)
        sys.exit(0 if ok else 1)
    if sys.argv[1] == "selftest":
        poly = 0x1053
        base = [(a, b) for a, b in fs.base_pairs(3, P8, "opt", [], None)]
        T8 = fs.target_tensor("grid", 3, P8, 101)
        cs = fs.solve_c(T8, base, 3, 3, 3)
        t8 = [(a, b, c) for (a, b), c in zip(base, cs)]
        up = lift(t8, 1, 1, poly)
        ok = fs.tensor_of(up, 12, 12) == fs.target_tensor("grid", 12, poly, 101)
        print("selftest: one GF(8) product of %d terms lifts to a K-product scheme of %d terms: tensor %s"
              % (len(t8), len(up), "PASS" if ok else "FAIL"))
        sys.exit(0 if ok else 1)
    src, p_, n_, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    poly = int(sys.argv[5], 16) if len(sys.argv) > 5 else 0x1053
    na, nb, nc, terms = fs.read_scheme(src)
    assert (na, nb, nc) == (3 * p_, 3 * n_, 3 * p_ * n_), "the scheme is not a GF(8) grid of these sizes"
    assert fs.tensor_of(terms, na, nb) == fs.target_tensor("grid", 3, P8, 100 * p_ + n_), "GF(8) scheme: tensor"
    up = lift(terms, p_, n_, poly)
    m = 100 * p_ + n_
    ok = fs.tensor_of(up, 12 * p_, 12 * n_) == fs.target_tensor("grid", 12, poly, m)
    print("%s: K(%d,%d) over GF(8) at %d terms -> over K (modulus %x) at %d terms: tensor identity %s"
          % (src, p_, n_, len(terms), poly, len(up), "PASS" if ok else "FAIL"))
    if ok:
        fs.write_scheme(out, up, 12 * p_, 12 * n_, 12 * p_ * n_)
        print("wrote", out, "(against %d separate K-products at 42: %d)" % (p_ * n_, 42 * p_ * n_))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
