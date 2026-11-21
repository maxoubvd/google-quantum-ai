"""R56-CHEM-ENGINE3 kill test (b): the middle-size lab (random instance n outputs x r high rows x a low coordinates).

Cell index = t * n + o (t low coordinate, o output). W = the r rows of the tensor in the cell space.
  1. LOCAL ENGINE: the blocks family (3-coordinate blocks): relations per block, rooms V(b, y), lines, features
     paired greedily over lines (cross-block and same-block), restarts. Count = pairs + singles; law = sum of
     block ranks / 2.
  2. b1, ALTERNATING EXACT SOLVES: at the plan, with the fan-outs u and one factor pair (alpha, a) of every
     product fixed, the linear system for the other pairs (beta, b): its nullity beyond the gauge (one per
     product), and for every product whether the system without it is solvable (the removal test).
  3. b2, NON-LOCAL WALK: one product removed; all other products' planes and fan-outs are free to move one
     at a time (random plane / random fan-out / swap of a point); objective = uncovered dimensions
     dim(W + S) - dim S (exact, by ranks); annealing; reports the best value reached.

  python -I r56_chem_engine3_mid.py --n 10 --r 26 --a 15 --seed 1 [--restarts 40] [--walk 20000]
"""
import argparse
import math
import random
import time


def bits(x):
    out = []
    while x:
        l = x & -x
        out.append(l.bit_length() - 1)
        x ^= l
    return out


class Span:
    """Row-reduced span of integers with combination tracking."""

    def __init__(self):
        self.piv = {}

    def reduce(self, v, c=0):
        while v:
            h = v.bit_length() - 1
            e = self.piv.get(h)
            if e is None:
                return v, c
            v ^= e[0]
            c ^= e[1]
        return 0, c

    def add(self, v, tag=0):
        v, c = self.reduce(v, tag)
        if v:
            self.piv[v.bit_length() - 1] = (v, c)
            return True
        return False

    def __len__(self):
        return len(self.piv)


def rank(vs):
    sp = Span()
    return sum(1 for v in vs if sp.add(v))


def nullspace(rows, ncols):
    """Basis of {x in F2^ncols : <row, x> = 0 for every row} (rows as ints over ncols bits)."""
    # eliminate on the transposed system: columns as ints over the rows
    cols = []
    for j in range(ncols):
        c = 0
        for i, rw in enumerate(rows):
            if rw >> j & 1:
                c |= 1 << i
        cols.append(c)
    sp = Span()
    out = []
    for j, c in enumerate(cols):
        v, comb = sp.reduce(c, 1 << j)
        if v:
            sp.piv[v.bit_length() - 1] = (v, comb)
        else:
            out.append(comb)
    return out


def span_basis(vs):
    sp = Span()
    out = []
    for v in vs:
        if sp.add(v):
            out.append(v)
    return out


def intersect(A, B, n):
    """Basis of span(A) cap span(B) in F2^n."""
    if not A or not B:
        return []
    # x = sum a_i A_i = sum b_j B_j : nullspace of [A; B] as combination
    sp = Span()
    out = []
    k = len(A)
    for i, v in enumerate(A + B):
        w, comb = sp.reduce(v, 1 << i)
        if w:
            sp.piv[w.bit_length() - 1] = (w, comb)
        else:
            x = 0
            for q in bits(comb & ((1 << k) - 1)):
                x ^= A[q]
            if x:
                out.append(x)
    return span_basis(out)


def feat_vec(y, s, n):
    v = 0
    for t in bits(y):
        v |= s << (t * n)
    return v


def load_instance(path):
    import json
    d = json.load(open(path))
    return d["n"], d["r"], d["a"], [int(x, 16) for x in d["rows"]]


def make_instance(n, r, a, seed):
    rng = random.Random(seed)
    while True:
        rows = [rng.getrandbits(n * a) for _ in range(r)]
        if rank(rows) == min(r, n * a):
            return rows


def counting(n, r, a):
    b = n * r * a
    per = n + 2 * (r + a) - 4
    P = 1
    while P * per - math.lgamma(P + 1) / math.log(2) < b:
        P += 1
    return P


# ------------------------------------------------------------------ the blocks family
def blocks_plan(rows, n, r, a, rng, verbose=False):
    nb = a // 3
    full = (1 << n) - 1
    rel = []
    rooms = {}
    brank = []
    for b in range(nb):
        sh = 3 * b * n
        proj = [(w >> sh) & ((1 << (3 * n)) - 1) for w in rows]
        brank.append(rank(proj))
        phi = nullspace(proj, 3 * n)             # relations of the block: functionals on its 3n cells
        rel.append(phi)
        for y in range(1, 8):
            fun = []
            for p in phi:
                f = 0
                for j in bits(y):
                    f ^= (p >> (j * n)) & full
                fun.append(f)
            rooms[(b, y)] = nullspace(fun, n)
    # candidate lines
    lines = []
    for b in range(nb):
        for y in range(1, 8):
            for b2 in range(b, nb):
                for y2 in range(1, 8):
                    if b2 == b and y2 <= y:
                        continue
                    I = intersect(rooms[(b, y)], rooms[(b2, y2)], n)
                    if I:
                        lines.append((b, y, b2, y2, I))
    sp = [Span() for _ in range(nb)]
    need = list(brank)
    prods = []                                    # (y1 (full low mask), y2, s)
    order = list(range(len(lines)))
    for _round in range(6):
        rng.shuffle(order)
        for li in order:
            b, y, b2, y2, I = lines[li]
            if len(sp[b]) >= need[b] or len(sp[b2]) >= need[b2]:
                continue
            if b == b2 and len(sp[b]) + 2 > need[b]:
                continue
            # a random vector of the line
            s = 0
            for v in I:
                if rng.getrandbits(1):
                    s ^= v
            if not s:
                s = I[0]
            f1 = feat_vec(y, s, n)
            f2 = feat_vec(y2, s, n)
            if b == b2:
                t = Span()
                t.piv = dict(sp[b].piv)
                if t.add(f1) and t.add(f2):
                    sp[b] = t
                    prods.append((y << (3 * b), y2 << (3 * b), s))
            else:
                if sp[b].reduce(f1)[0] and sp[b2].reduce(f2)[0]:
                    sp[b].add(f1)
                    sp[b2].add(f2)
                    prods.append((y << (3 * b), y2 << (3 * b2), s))
    pairs = len(prods)
    singles = 0
    for b in range(nb):
        ys = list(range(1, 8))
        guard = 0
        while len(sp[b]) < need[b] and guard < 2000:
            guard += 1
            y = rng.choice(ys)
            R = rooms[(b, y)]
            if not R:
                continue
            s = 0
            for v in R:
                if rng.getrandbits(1):
                    s ^= v
            if s and sp[b].add(feat_vec(y, s, n)):
                prods.append((y << (3 * b), 0, s))
                singles += 1
        if len(sp[b]) < need[b]:
            return None
    return dict(prods=prods, pairs=pairs, singles=singles, brank=brank, nlines=len(lines),
                relations=[len(p) for p in rel], rooms=[len(rooms[(b, 1)]) for b in range(nb)])


def features_of(prods, n):
    """List of (product index, slot, cell vector)."""
    out = []
    for i, (y1, y2, s) in enumerate(prods):
        out.append((i, 0, feat_vec(y1, s, n)))
        if y2:
            out.append((i, 1, feat_vec(y2, s, n)))
    return out


def solve_H(rows, feats):
    """H_f (ints over the rows) with rows = sum_f H_f (x) feature_f, or None."""
    sp = Span()
    for q, (_, _, v) in enumerate(feats):
        sp.add(v, 1 << q)
    H = [0] * len(feats)
    for M, w in enumerate(rows):
        v, c = sp.reduce(w)
        if v:
            return None
        for q in bits(c):
            H[q] |= 1 << M
    return H


def uncovered(rows, vecs):
    sp = Span()
    for v in vecs:
        sp.add(v)
    d = len(sp)
    for w in rows:
        sp.add(w)
    return len(sp) - d, d


# ------------------------------------------------------------------ b1: the alternating system
def alternating_system(rows, prods, n, r, a):
    """With u_p, alpha_p (= H of slot 0) and a_p (= y of slot 1) fixed: unknowns b_p (a bits) and beta_p (r bits).
    Equation space: (M, t, o) -> bit M * (a n) + t * n + o. Returns (unknowns, rank, nullity, removable list)."""
    feats = features_of(prods, n)
    H = solve_H(rows, feats)
    assert H is not None
    Hof = {}
    for q, (i, slot, _) in enumerate(feats):
        Hof[(i, slot)] = H[q]
    an = a * n
    cols = []                                     # (product, column vector)
    for i, (y1, y2, s) in enumerate(prods):
        alpha = Hof[(i, 0)]
        # b_p: unit e_t  -> alpha (x) e_t (x) s
        for t in range(a):
            v = 0
            for M in bits(alpha):
                v |= (s << (t * n)) << (M * an)
            cols.append((i, v))
        # beta_p: unit e_M -> e_M (x) a_p (x) s   (a_p = y2; if the product has one feature, a_p = 0: no column)
        if y2:
            fv = feat_vec(y2, s, n)
            for M in range(r):
                cols.append((i, fv << (M * an)))
    target = 0
    for M, w in enumerate(rows):
        target |= w << (M * an)
    sp = Span()
    for q, (_, v) in enumerate(cols):
        sp.add(v, 1 << q)
    rk = len(sp)
    resid, _ = sp.reduce(target)
    assert resid == 0
    nonzero = sum(1 for _, v in cols if v)
    # removal test: is the target in the span of the other products' columns?
    removable = []
    for i in range(len(prods)):
        sp2 = Span()
        for (j, v) in cols:
            if j != i:
                sp2.add(v)
        if sp2.reduce(target)[0] == 0:
            removable.append(i)
    return nonzero, rk, nonzero - rk, removable


# ------------------------------------------------------------------ b2: the non-local walk
def random_plane(a, rng):
    while True:
        y1 = rng.getrandbits(a)
        y2 = rng.getrandbits(a)
        if y1 and y2 and y1 != y2:
            return y1, y2


def walk(rows, prods, n, a, steps, rng, t0=0.6):
    P = len(prods)
    cur = [list(p) for p in prods]

    def vecs(ps):
        out = []
        for y1, y2, s in ps:
            out.append(feat_vec(y1, s, n))
            if y2:
                out.append(feat_vec(y2, s, n))
        return out
    e, d = uncovered(rows, vecs(cur))
    best = e
    hist = {e: 1}
    for it in range(steps):
        i = rng.randrange(P)
        old = list(cur[i])
        kind = rng.randrange(4)
        if kind == 0:
            cur[i][0], cur[i][1] = random_plane(a, rng)
        elif kind == 1:
            cur[i][2] = rng.getrandbits(n) or 1
        elif kind == 2:
            cur[i][rng.randrange(2)] ^= 1 << rng.randrange(a)
            if not cur[i][0] or cur[i][0] == cur[i][1]:
                cur[i] = old
                continue
        else:
            cur[i][2] ^= 1 << rng.randrange(n)
            if not cur[i][2]:
                cur[i] = old
                continue
        e2, d2 = uncovered(rows, vecs(cur))
        temp = t0 * (1 - it / steps) + 1e-9
        if e2 <= e or rng.random() < math.exp((e - e2) / temp):
            e = e2
            hist[e] = hist.get(e, 0) + 1
            best = min(best, e)
            if e == 0:
                return 0, it + 1, hist
        else:
            cur[i] = old
    return best, steps, hist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--r", type=int, default=26)
    ap.add_argument("--a", type=int, default=15)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--restarts", type=int, default=40)
    ap.add_argument("--walk", type=int, default=0)
    ap.add_argument("--drop", type=int, default=1)
    ap.add_argument("--noalt", action="store_true")
    ap.add_argument("--load", default="")
    A = ap.parse_args()
    n, r, a = A.n, A.r, A.a
    rng = random.Random(1000 + A.seed)
    if A.load:
        n, r, a, rows = load_instance(A.load)
    else:
        rows = make_instance(n, r, a, A.seed)
    cnt = counting(n, r, a)
    fixy = n * r * a / (n + 2 * r)
    t0 = time.time()
    best = None
    for _ in range(A.restarts):
        pl = blocks_plan(rows, n, r, a, rng)
        if pl is None:
            continue
        c = pl["pairs"] + pl["singles"]
        if best is None or c < best["pairs"] + best["singles"]:
            best = pl
    law = sum(best["brank"])
    P = best["pairs"] + best["singles"]
    feats = features_of(best["prods"], n)
    assert solve_H(rows, feats) is not None, "plan does not cover"
    print(f"instance n {n} r {r} a {a} seed {A.seed}: cells {n * a}, block ranks {best['brank']}, relations "
          f"{best['relations']}, rooms {best['rooms']}, lines {best['nlines']}", flush=True)
    print(f"LOCAL blocks engine: {P} products ({best['pairs']} pairs + {best['singles']} singles), features "
          f"{len(feats)}; law {law} features = {-(-law // 2)} products; fixed-y bound {fixy:.1f}; counting {cnt} "
          f"({time.time() - t0:.1f} s)", flush=True)
    if not A.noalt:
        t0 = time.time()
        unk, rk, nul, rem = alternating_system(rows, best["prods"], n, r, a)
        gauge = best["pairs"]
        print(f"b1 ALTERNATING SYSTEM (u, alpha, a fixed; beta, b free): unknowns {unk}, equations {n * r * a}, "
              f"rank {rk}, nullity {nul} = gauge {gauge} + {nul - gauge}; removable products {len(rem)} of {P} "
              f"({time.time() - t0:.1f} s)", flush=True)
    if A.walk:
        t0 = time.time()
        ps = list(best["prods"])
        rng.shuffle(ps)
        ps = ps[A.drop:]
        vec0 = [v for (_, _, v) in features_of(ps, n)]
        e0, d0 = uncovered(rows, vec0)
        bestv, used, hist = walk(rows, ps, n, a, A.walk, rng)
        print(f"b2 NON-LOCAL WALK: {A.drop} product(s) dropped ({len(ps)} left), uncovered at start {e0}; after "
              f"{used} moves best uncovered {bestv}; visits by uncovered {dict(sorted(hist.items()))} "
              f"({time.time() - t0:.1f} s)", flush=True)
    print(f"RESULT n {n} r {r} a {a} seed {A.seed}: local {P}, law {-(-law // 2)}, counting {cnt}", flush=True)


if __name__ == "__main__":
    main()
