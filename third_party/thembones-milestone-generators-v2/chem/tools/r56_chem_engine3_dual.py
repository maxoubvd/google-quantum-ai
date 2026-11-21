"""R56-CHEM-ENGINE3: a THIRD non-local start, in the DUAL: search over the relation space.

A plan with P products = a space Phi of c = n a - 2 P relations (functionals on the n a cells that vanish on
every row of W) whose annihilator is spanned by products: for a fan-out s the usable low vectors are
Y(s) = kernel of the c x a matrix Phi.s; the supply is sum_s Y(s) (x) s. The plan is valid iff the supply
spans the annihilator (dimension n a - c); its cost is sum_s ceil(d_s / 2) for the dimensions d_s taken at s.
The fan-outs are ENUMERATED (2^n - 1 of them), so this works for n <= 8 or so; the search moves Phi inside
W-perp (annealing on: missing dimensions, then odd leftovers).

  python -I r56_chem_engine3_dual.py --n 4 --r 11 --a 6 --seed 2 --P 9 [--iters 200000] [--tries 4]
"""
import argparse
import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r56_chem_engine3_mid as mid

bits = mid.bits


def par(x):
    return bin(x).count("1") & 1


def combine(basis, coef):
    v = 0
    for q in bits(coef):
        v ^= basis[q]
    return v


class Evaluator:
    def __init__(self, n, a, c):
        self.n, self.a, self.c = n, a, c
        self.N = n * a
        self.mask = (1 << n) - 1
        self.mode = "s" if n <= a or n <= 8 else "y"       # enumerate the smaller side
        self.S = list(range(1, 1 << n)) if self.mode == "s" else None
        self.Ypts = list(range(1, 1 << a)) if self.mode == "y" else None

    def supply(self, phis):
        """per fan-out s: basis of Y(s) (low vectors)."""
        n, a = self.n, self.a
        cols = [[(p >> (t * n)) & self.mask for t in range(a)] for p in phis]     # phi_t per relation
        out = {}
        if self.mode == "y":
            # rooms V(y) = kernel of the functionals y^T phi; then Y(s) = {y : s in V(y)}
            acc = {}
            for y in self.Ypts:
                ts = bits(y)
                fun = []
                for pc in cols:
                    f = 0
                    for t in ts:
                        f ^= pc[t]
                    fun.append(f)
                V = mid.nullspace(fun, n)
                d = len(V)
                if d == 0 or d > 11:
                    if d > 11:
                        V = V[:11]
                        d = 11
                    else:
                        continue
                for m in range(1, 1 << d):
                    sv = 0
                    for q in bits(m):
                        sv ^= V[q]
                    acc.setdefault(sv, []).append(y)
            for sv, ys in acc.items():
                out[sv] = mid.span_basis(ys)
            return out
        for s in self.S:
            rowsm = []
            for pc in cols:
                v = 0
                for t in range(a):
                    if par(pc[t] & s):
                        v |= 1 << t
                rowsm.append(v)
            Y = mid.nullspace(rowsm, a)
            if Y:
                out[s] = Y
        return out

    def score(self, phis, rng=None, orders=3):
        """(missing dimensions, products, supply, assignment): missing = (N - c) - rank(supply);
        assignment = [(s, [low vectors taken at s])] of the best order tried. Pass 1 takes PAIRS only (two
        new dimensions at one fan-out), pass 2 the singles."""
        sup = self.supply(phis)
        need = self.N - mid.rank(phis)
        keys = sorted(sup, key=lambda s: -len(sup[s]))
        multi = [s for s in keys if len(sup[s]) >= 2]
        best = None
        n = self.n
        for o in range(orders):
            if o and rng is not None:
                rng.shuffle(multi)
            sp = mid.Span()
            took = {}
            cost = 0
            for s in multi:
                pend = None
                for y in sup[s]:
                    v = mid.feat_vec(y, s, n)
                    if pend is None:
                        if sp.reduce(v)[0]:
                            pend = (y, v)
                    else:
                        save = dict(sp.piv)
                        sp.add(pend[1])
                        if sp.add(v):
                            took.setdefault(s, []).extend((pend[0], y))
                            cost += 1
                            pend = None
                        else:
                            sp.piv = save
            for s in keys:
                for y in sup[s]:
                    if sp.add(mid.feat_vec(y, s, n)):
                        took.setdefault(s, []).append(y)
                        cost += 1                      # a single: one product for one feature
            miss = need - len(sp)
            # singles at one fan-out can still be paired two by two
            cost = sum((len(v) + 1) // 2 for v in took.values())
            cand = (miss, cost)
            if best is None or cand < best[0]:
                best = (cand, [(s, v) for s, v in took.items()])
            if rng is None:
                break
        return best[0][0], best[0][1], sum(len(v) for v in sup.values()), best[1]


def block_relations(rows, n, a):
    out = []
    for b in range(a // 3):
        sh = 3 * b * n
        proj = [(w >> sh) & ((1 << (3 * n)) - 1) for w in rows]
        for p in mid.nullspace(proj, 3 * n):
            out.append(p << sh)
    return out


def search(rows, n, r, a, P, iters, rng, t_start=1.5, verbose=False, warm=False, deep=24):
    N = n * a
    c = N - 2 * P
    Wp = mid.nullspace(rows, N)
    D = len(Wp)
    if c > D or c < 0:
        return None
    ev = Evaluator(n, a, c)
    coefs = []
    sp = mid.Span()
    if warm:                                   # start from the blocks family's relations
        tag = mid.Span()
        for q, v in enumerate(Wp):
            tag.add(v, 1 << q)
        br = block_relations(rows, n, a)
        rng.shuffle(br)
        for p in br:
            if len(coefs) >= c:
                break
            res, cf = tag.reduce(p)
            assert res == 0
            if cf and sp.add(cf):
                coefs.append(cf)
    while len(coefs) < c:
        cf = rng.getrandbits(D)
        if cf and sp.add(cf):
            coefs.append(cf)
    phis = [combine(Wp, cf) for cf in coefs]

    def energy(ph):
        miss, cost, tot, asg = ev.score(ph, rng, 2)
        if miss == 0 and P < cost <= P + 2:            # endgame: the supply spans; look harder for the pairs
            miss2, cost2, tot2, asg2 = ev.score(ph, rng, deep)
            if (miss2, cost2) < (miss, cost):
                miss, cost, asg = miss2, cost2, asg2
        return 4.0 * miss + 1.0 * max(0, cost + (miss + 1) // 2 - P) - 0.02 * tot, miss, cost, asg
    e, miss, cost, asg = energy(phis)
    best = (miss, cost)
    t0 = time.time()
    for it in range(iters):
        if miss == 0 and cost <= P:
            return dict(P=cost, iters=it, secs=time.time() - t0, phis=phis, asg=asg)
        i = rng.randrange(c) if c else 0
        old_cf, old_phi = coefs[i], phis[i]
        kind = rng.randrange(3)
        if kind == 0:
            cf = old_cf ^ (1 << rng.randrange(D))
        elif kind == 1:
            cf = old_cf ^ coefs[rng.randrange(c)] if c > 1 else old_cf
        else:
            cf = old_cf ^ rng.getrandbits(D)
        if not cf:
            continue
        coefs[i] = cf
        if mid.rank(coefs) < c:
            coefs[i] = old_cf
            continue
        phis[i] = combine(Wp, cf)
        e2, miss2, cost2, asg2 = energy(phis)
        temp = t_start * (1 - it / iters) + 0.02
        if e2 <= e or rng.random() < math.exp((e - e2) / temp):
            e, miss, cost, asg = e2, miss2, cost2, asg2
            if (miss, cost) < best:
                best = (miss, cost)
                if verbose:
                    print(f"    it {it}: missing {miss}, products {cost}", flush=True)
        else:
            coefs[i], phis[i] = old_cf, old_phi
    return dict(P=None, best=best, iters=iters, secs=time.time() - t0)


def verify(rows, n, a, asg):
    """the explicit plan: products = pairs of low vectors at one fan-out; W inside their span?"""
    vecs = []
    prods = 0
    pts = set()
    for s, ys in asg:
        prods += (len(ys) + 1) // 2
        for y in ys:
            vecs.append(mid.feat_vec(y, s, n))
            pts.add(y)
    unc, dim = mid.uncovered(rows, vecs)
    nblk = sum(1 for y in pts if sum(1 for q in range(a // 3) if (y >> (3 * q)) & 7) == 1)
    return unc == 0 and dim == len(vecs), prods, dim, len(asg), len(pts), nblk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--r", type=int, default=11)
    ap.add_argument("--a", type=int, default=6)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--P", type=int, default=0, help="first target (default: the law), then downwards")
    ap.add_argument("--down", type=int, default=6)
    ap.add_argument("--iters", type=int, default=100000)
    ap.add_argument("--tries", type=int, default=3)
    ap.add_argument("--warm", action="store_true", help="start from the blocks family's relations")
    ap.add_argument("--tstart", type=float, default=1.5)
    ap.add_argument("--deep", type=int, default=24, help="orders tried for the pairs when the supply spans")
    A = ap.parse_args()
    n, r, a = A.n, A.r, A.a
    rows = mid.make_instance(n, r, a, A.seed)
    law = 0
    for b in range(a // 3):
        law += mid.rank([(w >> (3 * b * n)) & ((1 << (3 * n)) - 1) for w in rows])
    lawp = -(-law // 2)
    cnt = mid.counting(n, r, a)
    print(f"instance n {n} r {r} a {a} seed {A.seed}: law {lawp} products ({law} features); counting {cnt}; "
          f"fixed-y bound {n * r * a / (n + 2 * r):.1f}; relations on offer {n * a - r}", flush=True)
    rng = random.Random(900 + A.seed)
    P = A.P or lawp
    smallest = None
    for target in range(P, max(0, P - A.down) - 1, -1):
        found = None
        for tr in range(A.tries):
            res = search(rows, n, r, a, target, A.iters, rng, t_start=A.tstart, warm=A.warm, deep=A.deep)
            if res is None:
                break
            if res["P"] is not None:
                ok, prods, dim, nfan, npts, nblk = verify(rows, n, a, res["asg"])
                assert ok and prods <= target, "verification failed"
                found = (res["iters"], res["secs"], dim, nfan, tr, npts, nblk)
                break
            last = res
        if found:
            smallest = target
            print(f"  P {target} (relations {n * a - 2 * target}): FOUND, verified (W inside the span of the "
                  f"products, dim {found[2]}, {found[3]} fan-outs, {found[5]} low points of which {found[6]} on "
                  f"block points; try {found[4]}, {found[0]} moves, {found[1]:.0f} s)", flush=True)
        else:
            if res is None:
                print(f"  P {target}: more relations than exist", flush=True)
            else:
                print(f"  P {target} (relations {n * a - 2 * target}): not found in {A.tries} x {A.iters} moves "
                      f"(best: missing {last['best'][0]}, products {last['best'][1]})", flush=True)
            break
    print(f"RESULT n {n} r {r} a {a} seed {A.seed}: law {lawp}, dual search smallest {smallest}, counting {cnt}",
          flush=True)


if __name__ == "__main__":
    main()
