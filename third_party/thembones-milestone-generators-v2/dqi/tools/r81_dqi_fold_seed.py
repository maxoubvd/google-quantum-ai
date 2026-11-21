#!/usr/bin/env python3
"""r81_dqi_fold_seed.py -- fold and shared-scalar kernel tensors over GF(2^s): seeds and exact checks (seat R81-DQI-FOLD).

PUBLIC MATHEMATICS ONLY (products of polynomials over finite fields; Toom-Cook, Karatsuba, schoolbook seeds).
Scheme file format = tools/dqi_dq6_flip.c: "na nb nc R" then R lines "a b c" (hex masks).

Kernels (kind):
  fold  s poly m : (h0 + z h1) x (A_0 + .. + z^(m-1) A_(m-1)) over GF(2^s): na = 2s, nb = m s, nc = (m+1) s
  star  s poly m : one scalar k times m words: (k, A) -> (k A_0, .., k A_(m-1)): na = s, nb = m s, nc = m s
  modm  s poly m : diagonal head: outputs (h0 + t_i h1) A_i at m points t_i of P^1(GF(2^s)): na = 2s, nb = nc = m s
  outer s poly m : all 2 m word products h_i A_j, each its own output word (the complete bipartite pattern K_{2,m};
                   the fold of 2 m - 1 words restricts onto it): na = 2s, nb = m s, nc = 2 m s

Commands:
  seed <kind> <s> <polyhex> <m> <style> <out> [rngseed] [basefile ...]
      styles: toom[:p0,p1,..]  Toom at the listed points (field elements as ints, 'inf' allowed; default 0,1,..)
              school           every word product h_i A_j (star: k A_j)
              kblock           blocks of two words by Karatsuba (3 products), a last single word by 2
              toomx:<extra>    Toom with <extra> redundant points (loose)
      field-product scheme a product uses: env BASE = opt (default: Karatsuba, c solved) | naive | file (basefile list,
      cycled); env GAUGE = 0 (same scheme everywhere) | 1 (a random gauge (alpha x)(beta y)/(alpha beta) per product)
  check <kind> <s> <polyhex> <m> <scheme>     tensor identity + all inputs when na + nb <= 24 (bilinear or 'alt' files)
  embed <s> <polyhex> <m> <outer-2 scheme> <out> <j0,j1[;j0,j1..]> [singles j,j,..]
      a fold seed of m words from K(2,2) schemes on the word pairs (j0, j1) (j1 - j0 != 1 not required: overlapping
      output words are XORed) plus 12-term stars on the single words: an IRREGULAR loose seed
  embedo <s> <polyhex> <m> <outer-2 scheme> <out> <j0,j1[;..]> [singles]   the same for the K(2,m) kernel 'outer'
  graphseed <s> <polyhex> <nh> <nd> <out> <rngseed> <i:j> [<i:j> ..]
      a product graph: nh heads, nd words, one product h_i d_j per listed edge, each its own output word (in the
      order listed), every product a 6-term scheme in a random gauge (BASE=naive: 9 terms). The walker takes the
      seed's own tensor as its target, so the seed file IS the kernel's definition (built product by product).
  gridseed <s> <polyhex> <p> <n> <out> <rngseed> [<K(2,k) scheme> <i0,i1,j0,..,j(k-1)[;..]>] ...
      kernel 'grid' (m = 100 p + n): K(p,n), p heads x n words, every product its own output word; the listed blocks
      (heads i0, i1; words j0 ..) take the given K(2,k) scheme, every other product a 6-term scheme in a random gauge
  tensor <kind> <s> <polyhex> <m>             print dims
"""
import os
import random
import sys


def gfmul(x, y, s, poly):
    r = 0
    for i in range(s):
        if y >> i & 1:
            r ^= x << i
    for k in range(2 * s - 2, s - 1, -1):
        if r >> k & 1:
            r ^= poly << (k - s)
    return r


def gfinv(x, s, poly):
    r = 1
    e = (1 << s) - 2
    b = x
    while e:
        if e & 1:
            r = gfmul(r, b, s, poly)
        b = gfmul(b, b, s, poly)
        e >>= 1
    return r


def par(x):
    return bin(x).count("1") & 1


def pull(a, alpha, s, poly):
    """functional a' with a'(u) = a(alpha u)."""
    r = 0
    for i in range(s):
        if par(a & gfmul(alpha, 1 << i, s, poly)):
            r |= 1 << i
    return r


def dims(kind, s, m):
    if kind == "fold":
        return 2 * s, m * s, (m + 1) * s
    if kind == "star":
        return s, m * s, m * s
    if kind == "modm":
        return 2 * s, m * s, m * s
    if kind == "outer":
        return 2 * s, m * s, 2 * m * s
    if kind == "grid":
        p, n = divmod(m, 100)
        return p * s, n * s, p * n * s
    if kind == "fold2":
        return 4 * s, m * s, 2 * (m + 1) * s
    raise SystemExit("kind")


def modm_points(s, m):
    pts = list(range(1 << s)) + ["inf"]
    return pts[:m]


def apply_map(kind, s, poly, m, x, y):
    """the bilinear map on integer-packed operands."""
    msk = (1 << s) - 1
    if kind == "fold":
        h = [x & msk, x >> s & msk]
        out = 0
        for j in range(m):
            aj = y >> (j * s) & msk
            if not aj:
                continue
            for i in range(2):
                if h[i]:
                    out ^= gfmul(h[i], aj, s, poly) << ((i + j) * s)
        return out
    if kind == "star":
        out = 0
        for j in range(m):
            out ^= gfmul(x & msk, y >> (j * s) & msk, s, poly) << (j * s)
        return out
    if kind == "outer":
        out = 0
        for i in range(2):
            for j in range(m):
                out ^= gfmul(x >> (i * s) & msk, y >> (j * s) & msk, s, poly) << ((i * m + j) * s)
        return out
    if kind == "fold2":
        lo = apply_map("fold", s, poly, m, x & ((1 << (2 * s)) - 1), y)
        hi = apply_map("fold", s, poly, m, x >> (2 * s), y)
        return lo | hi << ((m + 1) * s)
    if kind == "grid":
        p, n = divmod(m, 100)
        out = 0
        for i in range(p):
            for j in range(n):
                out ^= gfmul(x >> (i * s) & msk, y >> (j * s) & msk, s, poly) << ((i * n + j) * s)
        return out
    if kind == "modm":
        h0, h1 = x & msk, x >> s & msk
        out = 0
        for j, t in enumerate(modm_points(s, m)):
            xt = h1 if t == "inf" else h0 ^ gfmul(t, h1, s, poly)
            out ^= gfmul(xt, y >> (j * s) & msk, s, poly) << (j * s)
        return out
    raise SystemExit("kind")


def target_tensor(kind, s, poly, m):
    na, nb, nc = dims(kind, s, m)
    return [[apply_map(kind, s, poly, m, 1 << i, 1 << j) for j in range(nb)] for i in range(na)]


def tensor_of(terms, na, nb):
    T = [[0] * nb for _ in range(na)]
    for a, b, c in terms:
        for i in range(na):
            if a >> i & 1:
                row = T[i]
                for j in range(nb):
                    if b >> j & 1:
                        row[j] ^= c
    return T


def solve_c(T, ab, na, nb, nc):
    """c_k with sum_k a_k (x) b_k (x) c_k = T, or None."""
    vecs = []
    for a, b in ab:
        v = 0
        for i in range(na):
            if a >> i & 1:
                v ^= b << (i * nb)
        vecs.append(v)
    piv = {}
    for k, v in enumerate(vecs):
        comb = 1 << k
        while v:
            t = v.bit_length() - 1
            if t in piv:
                v ^= piv[t][0]
                comb ^= piv[t][1]
            else:
                piv[t] = (v, comb)
                break
    cs = [0] * len(ab)
    for o in range(nc):
        v = 0
        for i in range(na):
            for j in range(nb):
                if T[i][j] >> o & 1:
                    v ^= 1 << (i * nb + j)
        comb = 0
        while v:
            t = v.bit_length() - 1
            if t not in piv:
                return None
            v ^= piv[t][0]
            comb ^= piv[t][1]
        for k in range(len(ab)):
            if comb >> k & 1:
                cs[k] |= 1 << o
    return cs


def read_scheme(path):
    with open(path) as f:
        head = f.readline().split()
        if head[0] == "alt":
            D, nc, r = int(head[1]), int(head[2]), int(head[3])
            terms = [tuple(int(t, 16) for t in f.readline().split()) for _ in range(r)]
            return ("alt", D, nc, terms)
        na, nb, nc, r = map(int, head)
        terms = [tuple(int(t, 16) for t in f.readline().split()) for _ in range(r)]
    return (na, nb, nc, terms)


def write_scheme(path, terms, na, nb, nc):
    with open(path, "w") as f:
        f.write("%d %d %d %d\n" % (na, nb, nc, len(terms)))
        for a, b, c in terms:
            f.write("%x %x %x\n" % (a, b, c))


def base_pairs(s, poly, style, files, rng):
    """(a, b) form pairs of a field-product scheme on s bits (c is solved later for the whole kernel)."""
    if style == "naive":
        return [(1 << i, 1 << j) for i in range(s) for j in range(s)]
    if style == "file":
        fn = files[rng.randrange(len(files))]
        na, nb, nc, terms = read_scheme(fn)
        assert na == s and nb == s
        return [(a, b) for a, b, c in terms]
    # opt: recursive Karatsuba pairs (3 for s = 2, 6 for s = 3)
    def kp(bits):
        if len(bits) == 1:
            return [(1 << bits[0],)]
        if len(bits) == 3:
            x0, x1, x2 = (1 << b for b in bits)
            return [(x0,), (x1,), (x2,), (x0 | x1,), (x0 | x2,), (x1 | x2,)]
        h = len(bits) // 2
        lo, hi = kp(bits[:h]), kp(bits[h:])
        mid = []
        if len(bits) == 2:
            mid = [(lo[0][0] | hi[0][0],)]
        else:
            raise SystemExit("opt base: s = 2 or 3 only; use BASE=file")
        return lo + hi + mid
    return [(p[0], p[0]) for p in kp(list(range(s)))]


def gauge(pairs, s, poly, rng, on):
    if not on:
        return pairs
    al = rng.randrange(1, 1 << s)
    be = rng.randrange(1, 1 << s)
    return [(pull(a, al, s, poly), pull(b, be, s, poly)) for a, b in pairs]


def lift(pairs, s, poly, xc, yc):
    """forms of the product (sum_i xc[i] X_i)(sum_j yc[j] Y_j): masks over the word-packed operands."""
    out = []
    for a, b in pairs:
        am = 0
        for i, c in enumerate(xc):
            if c:
                am |= pull(a, c, s, poly) << (i * s)
        bm = 0
        for j, c in enumerate(yc):
            if c:
                bm |= pull(b, c, s, poly) << (j * s)
        out.append((am, bm))
    return out


def products(kind, s, poly, m, style):
    """list of (xc, yc): coefficient vectors over GF(2^s) of the two factors of each word product."""
    q = 1 << s
    name, _, arg = style.partition(":")
    P = []
    if kind == "star":
        return [([1], [1 if j == w else 0 for j in range(m)]) for w in range(m)]
    if kind == "fold2":
        one = products("fold", s, poly, m, style)
        return [(xc + [0, 0], yc) for xc, yc in one] + [([0, 0] + xc, yc) for xc, yc in one]
    if kind == "modm":
        for j, t in enumerate(modm_points(s, m)):
            xc = [0, 1] if t == "inf" else [1, t]
            P.append((xc, [1 if i == j else 0 for i in range(m)]))
        return P
    if name == "school":
        for i in range(2):
            for j in range(m):
                P.append(([1 if k == i else 0 for k in range(2)], [1 if k == j else 0 for k in range(m)]))
        return P
    if name == "kblock":
        j = 0
        while j + 1 < m:
            e0 = [1 if k == j else 0 for k in range(m)]
            e1 = [1 if k == j + 1 else 0 for k in range(m)]
            P += [([1, 0], e0), ([0, 1], e1), ([1, 1], [x | y for x, y in zip(e0, e1)])]
            j += 2
        if j < m:
            e0 = [1 if k == j else 0 for k in range(m)]
            P += [([1, 0], e0), ([0, 1], e0)]
        return P
    if name in ("toom", "toomx"):
        if name == "toom" and arg:
            pts = [p if p == "inf" else int(p) for p in arg.split(",")]
        else:
            extra = int(arg) if (name == "toomx" and arg) else 0
            if m + extra > q:
                raise SystemExit("not enough points")
            pts = list(range(m + extra)) + ["inf"]
        for t in pts:
            if t == "inf":
                P.append(([0, 1], [0] * (m - 1) + [1]))
            else:
                yc, pw = [], 1
                for j in range(m):
                    yc.append(pw)
                    pw = gfmul(pw, t, s, poly)
                P.append(([1, t], yc))
        return P
    raise SystemExit("style")


def build(kind, s, poly, m, style, rng, base, gauge_on, files):
    na, nb, nc = dims(kind, s, m)
    ab = []
    for xc, yc in products(kind, s, poly, m, style):
        pairs = gauge(base_pairs(s, poly, base, files, rng), s, poly, rng, gauge_on)
        ab += lift(pairs, s, poly, xc, yc)
    T = target_tensor(kind, s, poly, m)
    cs = solve_c(T, ab, na, nb, nc)
    if cs is None:
        raise SystemExit("the products do not span the target")
    terms = [(a, b, c) for (a, b), c in zip(ab, cs) if c]
    assert tensor_of(terms, na, nb) == T
    return terms, (na, nb, nc)


def check(kind, s, poly, m, path):
    sch = read_scheme(path)
    na, nb, nc = dims(kind, s, m)
    T = target_tensor(kind, s, poly, m)
    if sch[0] == "alt":
        _, D, nc2, terms = sch
        assert D == na + nb and nc2 == nc, "dimensions"
        # sum_k c_k (x) (l_k ^ m_k) = Omega
        W = {}
        for l, mm, c in terms:
            for i in range(D):
                for j in range(i + 1, D):
                    if ((l >> i & 1) & (mm >> j & 1)) ^ ((l >> j & 1) & (mm >> i & 1)):
                        W[(i, j)] = W.get((i, j), 0) ^ c
        ok_t = True
        for i in range(D):
            for j in range(i + 1, D):
                want = T[i][j - na] if (i < na <= j) else 0
                if W.get((i, j), 0) != want:
                    ok_t = False
        ok_all = None
        if D <= 24:
            import numpy as np
            N = 1 << D
            v = np.arange(N, dtype=np.uint64)

            def parv(mask):
                w = v & np.uint64(mask)
                p = np.zeros(N, dtype=np.uint8)
                while w.any():
                    p ^= (w & np.uint64(1)).astype(np.uint8)
                    w = w >> np.uint64(1)
                return p.astype(bool)
            f = np.zeros(N, dtype=np.uint64)
            for l, mm, c in terms:
                f[parv(l) & parv(mm)] ^= np.uint64(c)
            xs = v & np.uint64((1 << na) - 1)
            ys = v >> np.uint64(na)
            ref = ref_table(kind, s, poly, m, xs, ys)
            g = f ^ ref                       # must be GF(2)-linear in v
            lin = np.zeros(N, dtype=np.uint64)
            for i in range(D):
                gi = int(g[1 << i])
                lin[(v >> np.uint64(i) & np.uint64(1)).astype(bool)] ^= np.uint64(gi)
            ok_all = bool((g == lin).all())
        return ok_t, ok_all, len(terms), "alt"
    na2, nb2, nc2, terms = sch
    assert (na2, nb2, nc2) == (na, nb, nc), "dimensions"
    ok_t = tensor_of(terms, na, nb) == T and all(a and b and c for a, b, c in terms)
    ok_all = None
    if na + nb <= 24:
        import numpy as np
        N = 1 << (na + nb)
        v = np.arange(N, dtype=np.uint64)
        xs = v & np.uint64((1 << na) - 1)
        ys = v >> np.uint64(na)

        def parv(w):
            p = np.zeros(N, dtype=np.uint8)
            while w.any():
                p ^= (w & np.uint64(1)).astype(np.uint8)
                w = w >> np.uint64(1)
            return p.astype(bool)
        f = np.zeros(N, dtype=np.uint64)
        for a, b, c in terms:
            f[parv(xs & np.uint64(a)) & parv(ys & np.uint64(b))] ^= np.uint64(c)
        ok_all = bool((f == ref_table(kind, s, poly, m, xs, ys)).all())
    return ok_t, ok_all, len(terms), "bilinear"


def ref_table(kind, s, poly, m, xs, ys):
    """the map on numpy arrays of packed operands, from the definition (carry-less products + reduction)."""
    import numpy as np
    msk = np.uint64((1 << s) - 1)

    def fm(u, w):
        r = np.zeros(u.shape, dtype=np.uint64)
        for i in range(s):
            r ^= np.where((w >> np.uint64(i) & np.uint64(1)).astype(bool), u << np.uint64(i), np.uint64(0))
        for k in range(2 * s - 2, s - 1, -1):
            r ^= np.where((r >> np.uint64(k) & np.uint64(1)).astype(bool), np.uint64(poly << (k - s)), np.uint64(0))
        return r
    out = np.zeros(xs.shape, dtype=np.uint64)
    if kind == "fold":
        h = [xs & msk, xs >> np.uint64(s) & msk]
        for j in range(m):
            aj = ys >> np.uint64(j * s) & msk
            for i in range(2):
                out ^= fm(h[i], aj) << np.uint64((i + j) * s)
    elif kind == "star":
        for j in range(m):
            out ^= fm(xs & msk, ys >> np.uint64(j * s) & msk) << np.uint64(j * s)
    elif kind == "outer":
        for i in range(2):
            for j in range(m):
                out ^= fm(xs >> np.uint64(i * s) & msk, ys >> np.uint64(j * s) & msk) << np.uint64((i * m + j) * s)
    elif kind == "grid":
        p, n = divmod(m, 100)
        for i in range(p):
            for j in range(n):
                out ^= fm(xs >> np.uint64(i * s) & msk, ys >> np.uint64(j * s) & msk) << np.uint64((i * n + j) * s)
    elif kind == "fold2":
        for hd in range(2):
            h = [xs >> np.uint64((2 * hd) * s) & msk, xs >> np.uint64((2 * hd + 1) * s) & msk]
            for j in range(m):
                aj = ys >> np.uint64(j * s) & msk
                for i in range(2):
                    out ^= fm(h[i], aj) << np.uint64((hd * (m + 1) + i + j) * s)
    elif kind == "modm":
        h0, h1 = xs & msk, xs >> np.uint64(s) & msk
        for j, t in enumerate(modm_points(s, m)):
            xt = h1 if t == "inf" else h0 ^ fm(h1, np.full(xs.shape, t, dtype=np.uint64))
            out ^= fm(xt, ys >> np.uint64(j * s) & msk) << np.uint64(j * s)
    return out


def main():
    cmd = sys.argv[1]
    if cmd == "seed":
        kind, s, poly, m, style, out = sys.argv[2], int(sys.argv[3]), int(sys.argv[4], 16), int(sys.argv[5]), sys.argv[6], sys.argv[7]
        rng = random.Random(int(sys.argv[8]) if len(sys.argv) > 8 else 1)
        files = sys.argv[9:]
        base = os.environ.get("BASE", "file" if files else "opt")
        terms, (na, nb, nc) = build(kind, s, poly, m, style, rng, base, os.environ.get("GAUGE", "0") == "1", files)
        write_scheme(out, terms, na, nb, nc)
        print("wrote %s: %s s=%d m=%d %s rank %d dims %d %d %d" % (out, kind, s, m, style, len(terms), na, nb, nc))
    elif cmd == "check":
        kind, s, poly, m, path = sys.argv[2], int(sys.argv[3]), int(sys.argv[4], 16), int(sys.argv[5]), sys.argv[6]
        ok_t, ok_all, r, model = check(kind, s, poly, m, path)
        print("%s: %s rank %d tensor %s all-inputs %s" % (path, model, r, "PASS" if ok_t else "FAIL",
                                                        "n/a" if ok_all is None else ("PASS" if ok_all else "FAIL")))
        sys.exit(0 if ok_t and ok_all in (None, True) else 1)
    elif cmd in ("embed", "embedo"):
        s_, poly, m, src, out = int(sys.argv[2]), int(sys.argv[3], 16), int(sys.argv[4]), sys.argv[5], sys.argv[6]
        pairs = [tuple(int(x) for x in p.split(",")) for p in sys.argv[7].split(";") if p]
        singles = [int(x) for x in sys.argv[8].split(",")] if len(sys.argv) > 8 and sys.argv[8] else []
        kind = "fold" if cmd == "embed" else "outer"
        na, nb, nc = dims(kind, s_, m)
        _, _, _, ot = read_scheme(src)
        msk = (1 << s_) - 1

        def slot(i, j):                      # output word of the product h_i A_j
            return (i + j) if kind == "fold" else (i * m + j)
        terms = []
        for j0, j1 in pairs:
            for a, b, c in ot:
                bm = (b & msk) << (j0 * s_) | (b >> s_ & msk) << (j1 * s_)
                cm = 0
                for i in range(2):
                    for jj, j in enumerate((j0, j1)):
                        cm ^= (c >> ((i * 2 + jj) * s_) & msk) << (slot(i, j) * s_)
                terms.append((a, bm, cm))
        rng = random.Random(7)
        for j in singles:
            for i in range(2):
                for a, b in gauge(base_pairs(s_, poly, "opt", [], rng), s_, poly, rng, True):
                    terms.append((a << (i * s_), b << (j * s_), None, slot(i, j)))
        # the stars' c: solve word by word (a single GF(2^s) product into one output word)
        fixed = [t for t in terms if len(t) == 3]
        T = target_tensor(kind, s_, poly, m)
        Tf = tensor_of(fixed, na, nb)
        rest = [[T[i][j] ^ Tf[i][j] for j in range(nb)] for i in range(na)]
        star = [t for t in terms if len(t) == 4]
        cs = solve_c(rest, [(a, b) for a, b, _, _ in star], na, nb, nc) if star else []
        if cs is None:
            raise SystemExit("embed: stars do not span the rest")
        allt = [t for t in fixed if t[2]] + [(a, b, c) for (a, b, _, _), c in zip(star, cs) if c]
        assert tensor_of(allt, na, nb) == T, "embed: tensor"
        write_scheme(out, allt, na, nb, nc)
        print("wrote %s: %s s=%d m=%d embed rank %d dims %d %d %d" % (out, kind, s_, m, len(allt), na, nb, nc))
    elif cmd == "graphseed":
        s_, poly, nh, nd, out = int(sys.argv[2]), int(sys.argv[3], 16), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
        rng = random.Random(int(sys.argv[7]))
        edges = []
        for k, e in enumerate(sys.argv[8:]):          # "i:j" (own output word), "i:j:e" (output word e: sums),
            q = e.split(":")                          # "i:j:e*c+e2*c2" (the product times the constants c, c2
            outs = []                                 #  of GF(2^s), as ints, into the words e, e2: moments)
            if len(q) > 2:
                for item in q[2].split("+"):
                    w, _, cst = item.partition("*")
                    outs.append((int(w), int(cst) if cst else 1))
            else:
                outs.append((k, 1))
            edges.append((int(q[0]), int(q[1]), outs))
        nout = 1 + max(w for _, _, outs in edges for w, _ in outs)
        T1 = target_tensor("grid", s_, poly, 101)
        terms = []
        for i, j, outs in edges:
            pairs = gauge(base_pairs(s_, poly, os.environ.get("BASE", "opt"), [], rng), s_, poly, rng, True)
            cs = solve_c(T1, pairs, s_, s_, s_)
            for (a, b), c in zip(pairs, cs):
                cm = 0
                for w, cst in outs:
                    cm ^= gfmul(cst, c, s_, poly) << (w * s_)
                if cm:
                    terms.append((a << (i * s_), b << (j * s_), cm))
        write_scheme(out, terms, nh * s_, nd * s_, nout * s_)
        print("wrote %s: graph of %d products rank %d dims %d %d %d" % (out, len(edges), len(terms), nh * s_, nd * s_,
                                                                    nout * s_))
    elif cmd == "gridseed":
        s_, poly, p_, n_, out = int(sys.argv[2]), int(sys.argv[3], 16), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
        rng = random.Random(int(sys.argv[7]))
        m = 100 * p_ + n_
        na, nb, nc = dims("grid", s_, m)
        msk = (1 << s_) - 1
        fixed, done = [], set()
        for q in range(8, len(sys.argv) - 1, 2):        # pairs: <K(2,k) scheme> <i0,i1,j0,..,j(k-1)[;..]>
            ona, onb, onc, ot = read_scheme(sys.argv[q])
            k_ = onb // s_
            assert ona == 2 * s_ and onc == 2 * k_ * s_, "block scheme must be a K(2,k) grid"
            for blk in sys.argv[q + 1].split(";"):
                idx = [int(x) for x in blk.split(",")]
                hs, ws = idx[:2], idx[2:]
                assert len(ws) == k_
                for a, b, c in ot:
                    am = (a & msk) << (hs[0] * s_) | (a >> s_ & msk) << (hs[1] * s_)
                    bm = 0
                    for jj, j in enumerate(ws):
                        bm |= (b >> (jj * s_) & msk) << (j * s_)
                    cm = 0
                    for ii, i in enumerate(hs):
                        for jj, j in enumerate(ws):
                            cm ^= (c >> ((ii * k_ + jj) * s_) & msk) << ((i * n_ + j) * s_)
                    fixed.append((am, bm, cm))
                for i in hs:
                    for j in ws:
                        assert (i, j) not in done
                        done.add((i, j))
        ab = []
        for i in range(p_):
            for j in range(n_):
                if (i, j) in done:
                    continue
                for a, b in gauge(base_pairs(s_, poly, os.environ.get("BASE", "opt"), [], rng), s_, poly, rng, True):
                    ab.append((a << (i * s_), b << (j * s_)))
        T = target_tensor("grid", s_, poly, m)
        Tf = tensor_of(fixed, na, nb)
        rest = [[T[i][j] ^ Tf[i][j] for j in range(nb)] for i in range(na)]
        cs = solve_c(rest, ab, na, nb, nc) if ab else []
        if cs is None:
            raise SystemExit("gridseed: span")
        allt = fixed + [(a, b, c) for (a, b), c in zip(ab, cs) if c]
        assert tensor_of(allt, na, nb) == T, "gridseed: tensor"
        write_scheme(out, allt, na, nb, nc)
        print("wrote %s: grid s=%d K(%d,%d) rank %d dims %d %d %d" % (out, s_, p_, n_, len(allt), na, nb, nc))
    elif cmd == "tensor":
        kind, s, m = sys.argv[2], int(sys.argv[3]), int(sys.argv[5])
        print(dims(kind, s, m))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
