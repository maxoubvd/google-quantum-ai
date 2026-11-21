"""R56-CHEM-ENGINE3: the dual search with moves that KEEP THE LINES (the reopener of section 4, tested here).

From the blocks family's plan keep a random set K of its products. Every relation space Phi that vanishes on
K keeps those products (their planes stay in the kernels of their fan-outs). The relations that vanish on W
and on K form a linear space of dimension D = cells - rank(W + K). K is grown until D = c + gap, c = cells
- 2 P the number of relations a plan of P products needs. Then
  gap 1: EVERY hyperplane of that space is tried (2^D - 1 relation spaces, exhaustive);
  gap 2+: random subspaces of codimension gap.
Each Phi is priced exactly (supply against the annihilator, pairs at one fan-out, singles) by the dual
search's evaluator; a plan at P products is rebuilt and checked (W inside the span of its products).

  python -I r56_chem_engine3_lns.py --n 12 --r 31 --a 6 --seed 1 --P 30 --rounds 20 [--gap 1] [--samples 4000]
"""
import argparse
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r56_chem_engine3_mid as mid
import r56_chem_engine3_dual as dual

bits = mid.bits


class LineEvaluator:
    """Exact pricing of a relation space when the rooms are too large to enumerate their fan-outs (real
    output counts): rooms V(y) for every low point (2^a - 1 points, a <= 9), the supply from the rooms,
    products = two points with a common fan-out (the meet of two rooms), greedy over the lines, singles."""

    def __init__(self, n, a, c):
        self.n, self.a, self.c = n, a, c
        self.N = n * a
        self.mask = (1 << n) - 1
        self.Ypts = list(range(1, 1 << a))

    def score(self, phis, rng=None, orders=2):
        n, a = self.n, self.a
        cols = [[(p >> (t * n)) & self.mask for t in range(a)] for p in phis]
        rooms = {}
        for y in self.Ypts:
            ts = bits(y)
            fun = []
            for pc in cols:
                f = 0
                for t in ts:
                    f ^= pc[t]
                fun.append(f)
            V = mid.nullspace(fun, n)
            if V:
                rooms[y] = V
        need = self.N - mid.rank(phis)
        sp = mid.Span()
        tot = 0
        for y, V in rooms.items():
            tot += len(V)
            for sv in V:
                sp.add(mid.feat_vec(y, sv, n))
        miss = need - len(sp)
        if miss > 0:
            return miss, need, tot, None
        pts = sorted(rooms, key=lambda y: -len(rooms[y]))
        lines = []
        for i in range(len(pts)):
            di = len(rooms[pts[i]])
            for j in range(i + 1, len(pts)):
                if di + len(rooms[pts[j]]) < n - 6:
                    break
                I = mid.intersect(rooms[pts[i]], rooms[pts[j]], n)
                if I:
                    lines.append((pts[i], pts[j], I))
        best = None
        rr = rng or random.Random(1)
        for o in range(orders):
            sp2 = mid.Span()
            asg = []
            order = list(range(len(lines)))
            for rep in range(5):
                rr.shuffle(order)
                for li in order:
                    if len(sp2) + 2 > need:
                        break
                    y1, y2, I = lines[li]
                    sv = 0
                    for v in I:
                        if rr.getrandbits(1):
                            sv ^= v
                    if not sv:
                        sv = I[0]
                    save = dict(sp2.piv)
                    if sp2.add(mid.feat_vec(y1, sv, n)) and sp2.add(mid.feat_vec(y2, sv, n)):
                        asg.append((sv, [y1, y2]))
                    else:
                        sp2.piv = save
            for y in pts:
                if len(sp2) >= need:
                    break
                for sv in rooms[y]:
                    if sp2.add(mid.feat_vec(y, sv, n)):
                        asg.append((sv, [y]))
            cand = (need - len(sp2), len(asg))
            if best is None or cand < best[0]:
                best = (cand, asg)
        return best[0][0], best[0][1], tot, best[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--r", type=int, default=31)
    ap.add_argument("--a", type=int, default=6)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--P", type=int, default=30)
    ap.add_argument("--rounds", type=int, default=10)
    ap.add_argument("--gap", type=int, default=1)
    ap.add_argument("--samples", type=int, default=4000)
    ap.add_argument("--rseed", type=int, default=0)
    ap.add_argument("--load", default="")
    ap.add_argument("--lines", action="store_true", help="price by rooms and lines (real output counts)")
    A = ap.parse_args()
    n, r, a = A.n, A.r, A.a
    if A.load:
        n, r, a, rows = mid.load_instance(A.load)
    else:
        rows = mid.make_instance(n, r, a, A.seed)
    N = n * a
    rng = random.Random(4000 + 17 * A.seed + A.rseed)
    best = None
    for _ in range(20):
        pl = mid.blocks_plan(rows, n, r, a, rng)
        if pl and (best is None or pl["pairs"] + pl["singles"] < best["pairs"] + best["singles"]):
            best = pl
    law = best["pairs"] + best["singles"]
    c = N - 2 * A.P
    ev = LineEvaluator(n, a, c) if A.lines else dual.Evaluator(n, a, c)
    print(f"instance n {n} r {r} a {a} seed {A.seed}{' ' + A.load if A.load else ''}: blocks plan {law} products; target {A.P} products = {c} "
          f"relations; counting {mid.counting(n, r, a)}; gap {A.gap}", flush=True)
    pairs = [p for p in best["prods"] if p[1]]
    t0 = time.time()
    overall = None
    hist = {}
    for rd in range(A.rounds):
        rng.shuffle(pairs)
        sp = mid.Span()
        for w in rows:
            sp.add(w)
        kept = []
        kv = []
        fails = 0
        for (y1, y2, s) in pairs:
            f1, f2 = mid.feat_vec(y1, s, n), mid.feat_vec(y2, s, n)
            save = dict(sp.piv)
            sp.add(f1)
            sp.add(f2)
            if N - len(sp) < c + A.gap:
                sp.piv = save
                fails += 1
                if fails > 6:
                    break
                continue
            kept.append((y1, y2, s))
            kv += [f1, f2]
        Wp = mid.nullspace(rows + kv, N)
        D = len(Wp)
        if D < c:
            continue
        nbest = None
        ntried = 0
        hyper = D - c == 1
        it = iter(range(1, 1 << D)) if hyper and D <= 13 else None
        while True:
            if hyper:
                if it is not None:
                    lam = next(it, None)
                else:
                    lam = rng.getrandbits(D) if ntried < A.samples else None
                if not lam:
                    if lam is None:
                        break
                    continue
                i0 = (lam & -lam).bit_length() - 1
                basis = []
                for i in range(D):
                    if i != i0:
                        basis.append((1 << i) | ((1 << i0) if lam >> i & 1 else 0))
            else:
                if ntried >= A.samples:
                    break
                spc = mid.Span()
                basis = []
                while len(basis) < c:
                    cf = rng.getrandbits(D)
                    if cf and spc.add(cf):
                        basis.append(cf)
            ntried += 1
            phis = [dual.combine(Wp, cf) for cf in basis]
            miss, cost, tot, asg = ev.score(phis, rng, 2)
            if miss == 0 and A.P < cost <= A.P + 2:
                m2, c2, t2, a2 = ev.score(phis, rng, 24)
                if (m2, c2) < (miss, cost):
                    miss, cost, asg = m2, c2, a2
            key = (miss, cost)
            hist[key] = hist.get(key, 0) + 1
            if nbest is None or key < nbest:
                nbest = key
            if miss == 0 and cost <= A.P:
                ok, prods, dim, nfan, npts, nblk = dual.verify(rows, n, a, asg)
                assert ok and prods <= A.P
                print(f"  round {rd}: kept {len(kept)} products, relations on offer {D}: FOUND {prods} products, "
                      f"verified (W inside the span, dim {dim}, {nfan} fan-outs, {npts} low points of which "
                      f"{nblk} on block points); {ntried} relation spaces tried, {time.time() - t0:.0f} s",
                      flush=True)
                print(f"RESULT n {n} r {r} a {a} seed {A.seed}: blocks {law}, line-keeping dual search {prods}, "
                      f"counting {mid.counting(n, r, a)}", flush=True)
                return
        if overall is None or (nbest is not None and nbest < overall):
            overall = nbest
        print(f"  round {rd}: kept {len(kept)} products, relations on offer {D}, {ntried} relation spaces tried: "
              f"best (missing, products) {nbest} ({time.time() - t0:.0f} s)", flush=True)
    top = sorted(hist.items())[:6]
    print(f"RESULT n {n} r {r} a {a} seed {A.seed}: blocks {law}, line-keeping dual search: no plan at {A.P} in "
          f"{A.rounds} rounds; best (missing, products) {overall}; lowest states {top}", flush=True)


if __name__ == "__main__":
    main()
