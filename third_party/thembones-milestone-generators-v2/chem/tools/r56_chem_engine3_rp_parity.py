"""R56-CHEM-ENGINE3 re-pricer: EXACT pairing of the designer's line-keeping candidates.

The designer priced every relation space Phi with a GREEDY pairing (r56_chem_engine3_lns.LineEvaluator at real
output counts, r56_chem_engine3_dual.Evaluator at small ones) and handed the exact pairing to C2 as an open test.
Here each Phi is priced EXACTLY:
  S = annihilator of Phi, need = dim S; the supply (all rank-one cells y (x) s inside S) must span S;
  a product is a plane Y (x) s with Y inside Y(s) = ker(Phi.s); any plan spanning S has at least need - nu
  products, nu = the linear matroid parity of all such planes (the most planes whose union is independent),
  and need - nu is reached when the supply spans S [P: a line adds <= 2 dimensions, and adds 2 only when it is
  independent of the earlier ones].
  nu = rank / 2 of Lovasz's matrix  sum_s (a random alternating form on Y(s)) (x) (s s^T)  over GF(2^20)
  (random substitution: an error can only UNDER-state nu, probability <= nu / 2^20 per evaluation).

Modes (the kept-set loop is the designer's own, line by line, on the same random stream: the same kept sets
as r56_chem_engine3_lns.py with the same arguments; the exact pricing uses its own numpy generator):
  python -I r56_chem_engine3_rp_parity.py --load stageB_k1.json --lines --P 73 --gap 0 --samples 1 --rounds R --rseed S
  python -I r56_chem_engine3_rp_parity.py --n 14 --r 36 --a 6 --seed 1 --P 35 --gap 0 --samples 1 --rounds R
  ... --law      price the blocks family's own relation space (validation: must read the law)
"""
import argparse
import json
import os
import random
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r56_chem_engine3_mid as mid      # noqa: E402  (designer's library, imported read-only)
import r56_chem_engine3_dual as dual    # noqa: E402
import r56_chem_engine3_lns as lns      # noqa: E402

bits = mid.bits

# ------------------------------------------------------------------ GF(2^20)
MB = 20
POLY = (1 << 20) | (1 << 3) | 1           # x^20 + x^3 + 1
ORD = (1 << MB) - 1


def _tables():
    exp = np.zeros(2 * ORD + 2, dtype=np.int64)
    log = np.zeros(1 << MB, dtype=np.int64)
    x = 1
    e = []
    for i in range(ORD):
        e.append(x)
        x <<= 1
        if x >> MB:
            x ^= POLY
    assert x == 1
    ea = np.array(e, dtype=np.int64)
    assert np.unique(ea).size == ORD, "polynomial not primitive"
    exp[:ORD] = ea
    exp[ORD:2 * ORD] = ea
    log[ea] = np.arange(ORD, dtype=np.int64)
    return exp, log


EXP, LOG = _tables()


def gf_rank(M):
    """rank over GF(2^20) of an int64 matrix (entries < 2^20); M is overwritten."""
    nr, nc = M.shape
    r = 0
    for c in range(nc):
        if r == nr:
            break
        nz = np.flatnonzero(M[r:, c])
        if nz.size == 0:
            continue
        p = r + nz[0]
        if p != r:
            M[[r, p]] = M[[p, r]]
        row = M[r]
        cols = np.flatnonzero(row)
        inv = ORD - LOG[row[c]]
        row[cols] = EXP[LOG[row[cols]] + inv]
        below = r + 1 + np.flatnonzero(M[r + 1:, c])
        if below.size:
            f = M[below, c]
            M[np.ix_(below, cols)] ^= EXP[LOG[f][:, None] + LOG[row[cols]][None, :]]
        r += 1
    return r


# pattern of the alternating form u v^T + v u^T on the a x a low matrix, as a mask over the pairs t < t'
def _pairs(a):
    return [(t, u) for t in range(a) for u in range(t + 1, a)]


class Exact:
    def __init__(self, n, a, seed=12345):
        self.n, self.a = n, a
        self.N = n * a
        self.mask = (1 << n) - 1
        self.pr = _pairs(a)
        self.gen = np.random.default_rng(seed)
        A = 1 << a
        self.pat = [[0] * A for _ in range(A)]
        for u in range(A):
            for v in range(A):
                m = 0
                for q, (t, t2) in enumerate(self.pr):
                    b = ((u >> t) & (v >> t2) & 1) ^ ((v >> t) & (u >> t2) & 1)
                    if b:
                        m |= 1 << q
                self.pat[u][v] = m

    def rooms(self, phis):
        n, a = self.n, self.a
        cols = [[(p >> (t * n)) & self.mask for t in range(a)] for p in phis]
        rooms, funs = {}, {}
        for y in range(1, 1 << a):
            ts = bits(y)
            fun = []
            for pc in cols:
                f = 0
                for t in ts:
                    f ^= pc[t]
                fun.append(f)
            funs[y] = fun
            rooms[y] = mid.nullspace(fun, n)
        return cols, rooms, funs

    def fanouts(self, cols, rooms, funs, maxmeet=16):
        """every fan-out s with dim Y(s) >= 2, with a basis of Y(s): the meets of two rooms, one per plane."""
        a = self.a
        fans = set()
        trunc = 0
        for p in range(1, 1 << a):
            Vp = rooms[p]
            if not Vp:
                continue
            for q in range(p + 1, 1 << a):
                if (p ^ q) < q or not rooms[q]:
                    continue
                # meet = {s = Vp alpha : fun_q(s) = 0}: kernel of the c x dim(Vp) matrix fun_q . Vp
                rws = []
                for f in funs[q]:
                    m = 0
                    for i, v in enumerate(Vp):
                        if (f & v).bit_count() & 1:
                            m |= 1 << i
                    rws.append(m)
                K = mid.nullspace(rws, len(Vp))
                if not K:
                    continue
                I = []
                for al in K:
                    s = 0
                    for i in bits(al):
                        s ^= Vp[i]
                    I.append(s)
                k = len(I)
                if k > maxmeet:
                    trunc += 1
                    I = I[:maxmeet]
                    k = maxmeet
                # all non-zero s of the meet (Gray code)
                s = 0
                for g in range(1, 1 << k):
                    s ^= I[(g & -g).bit_length() - 1]
                    fans.add(s)
        out = []
        for s in fans:
            rowsm = []
            for pc in cols:
                v = 0
                for t in range(a):
                    if (pc[t] & s).bit_count() & 1:
                        v |= 1 << t
                rowsm.append(v)
            Y = mid.nullspace(rowsm, a)
            if len(Y) >= 2:
                out.append((s, Y))
        return out, trunc

    def nu(self, fans, gen=None):
        """matroid parity of all planes Y (x) s, Y inside Y(s): rank / 2 of the random Lovasz matrix."""
        gen = gen or self.gen
        n, a = self.n, self.a
        L = len(fans)
        if L == 0:
            return 0
        npr = len(self.pr)
        Av = np.zeros((L, npr), dtype=np.int64)
        Sm = np.zeros((L, n), dtype=np.float32)
        for l, (s, Y) in enumerate(fans):
            for o in bits(s):
                Sm[l, o] = 1.0
            d = len(Y)
            m = d * (d - 1) // 2
            xs = gen.integers(1, 1 << MB, size=m)
            row = [0] * npr
            k = 0
            for i in range(d):
                for j in range(i + 1, d):
                    x = int(xs[k])
                    k += 1
                    pm = self.pat[Y[i]][Y[j]]
                    while pm:
                        lo = pm & -pm
                        row[lo.bit_length() - 1] ^= x
                        pm ^= lo
            Av[l] = row
        Yp = np.zeros((n, npr, n), dtype=np.int64)
        for kb in range(MB):
            Wk = ((Av >> kb) & 1).astype(np.float32)
            Z = (Wk[:, :, None] * Sm[:, None, :]).reshape(L, npr * n)
            G = Sm.T @ Z
            Yp ^= (np.rint(G).astype(np.int64) & 1).reshape(n, npr, n) << kb
        N = self.N
        M = np.zeros((N, N), dtype=np.int64)
        for q, (t, t2) in enumerate(self.pr):
            B = Yp[:, q, :]
            M[t * n:(t + 1) * n, t2 * n:(t2 + 1) * n] = B
            M[t2 * n:(t2 + 1) * n, t * n:(t + 1) * n] = B.T
        rk = gf_rank(M)
        assert rk % 2 == 0, rk
        return rk // 2

    def price(self, phis, reps=1):
        need = self.N - mid.rank(phis)
        cols, rooms, funs = self.rooms(phis)
        sp = mid.Span()
        for y, V in rooms.items():
            for v in V:
                sp.add(mid.feat_vec(y, v, self.n))
        miss = need - len(sp)
        fans, trunc = self.fanouts(cols, rooms, funs)
        nu = max(self.nu(fans) for _ in range(reps))
        maxd = max((len(Y) for _, Y in fans), default=0)
        return dict(need=need, miss=miss, nu=nu, prod=(need - nu) if miss == 0 else None, nfan=len(fans),
                    maxd=maxd, trunc=trunc)


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
    ap.add_argument("--lines", action="store_true")
    ap.add_argument("--law", action="store_true", help="price the blocks family's relation space and stop")
    ap.add_argument("--every", type=int, default=200)
    ap.add_argument("--save", default="", help="json file for the relation spaces that read <= P + 1 exactly")
    ap.add_argument("--secs", type=float, default=0, help="stop after this many seconds")
    ap.add_argument("--gmax", type=int, default=-1, help="exact pricing only when the greedy reads <= gmax "
                    "(and on the calibration rounds); -1: every relation space")
    ap.add_argument("--calib", type=int, default=0, help="exact pricing on every relation space of the first rounds")
    ap.add_argument("--check", type=int, default=0, help="print the greedy reading of the first rounds")
    ap.add_argument("--savemax", type=int, default=-1, help="save only the spaces that read <= savemax (-1: P + 1)")
    A = ap.parse_args()
    n, r, a = A.n, A.r, A.a
    if A.load:
        n, r, a, rows = mid.load_instance(A.load)
    else:
        rows = mid.make_instance(n, r, a, A.seed)
    N = n * a
    ex = Exact(n, a, seed=777 + A.rseed + 1000 * A.seed)
    if A.law:
        br = dual.block_relations(rows, n, a)
        t0 = time.time()
        res = ex.price(br, reps=2)
        print(f"LAW n {n} r {r} a {a}{' ' + A.load if A.load else ' seed ' + str(A.seed)}: block relations "
              f"{len(br)} (rank {mid.rank(br)}), need {res['need']}, missing {res['miss']}, nu {res['nu']}, "
              f"EXACT products {res['prod']}; fan-outs {res['nfan']}, max dim Y(s) {res['maxd']}, truncated meets "
              f"{res['trunc']} ({time.time() - t0:.1f} s)", flush=True)
        return
    # ---- the designer's kept-set loop (tools/r56_chem_engine3_lns.py main), unchanged in its use of rng
    rng = random.Random(4000 + 17 * A.seed + A.rseed)
    best = None
    for _ in range(20):
        pl = mid.blocks_plan(rows, n, r, a, rng)
        if pl and (best is None or pl["pairs"] + pl["singles"] < best["pairs"] + best["singles"]):
            best = pl
    law = best["pairs"] + best["singles"]
    c = N - 2 * A.P
    ev = lns.LineEvaluator(n, a, c) if A.lines else dual.Evaluator(n, a, c)
    print(f"instance n {n} r {r} a {a} seed {A.seed}{' ' + A.load if A.load else ''}: blocks plan {law} products; "
          f"target {A.P} products = {c} relations; counting {mid.counting(n, r, a)}; gap {A.gap}; rseed {A.rseed}",
          flush=True)
    pairs = [p for p in best["prods"] if p[1]]
    t0 = time.time()
    hist_g, hist_e, joint = {}, {}, {}
    saved = []
    nk = 0
    lastp = -1
    worse = 0
    for rd in range(A.rounds):
        if A.secs and time.time() - t0 > A.secs:
            break
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
            kg = (miss, cost)
            hist_g[kg] = hist_g.get(kg, 0) + 1
            if rd < A.check:
                print(f"  check round {rd}: greedy {kg}", flush=True)
            if A.gmax >= 0 and rd >= A.calib and not (miss == 0 and cost <= A.gmax):
                continue
            # ---- the exact price of the same relation space
            res = ex.price(phis)
            nk += 1
            if A.gmax >= 0:
                print(f"  exact round {rd}: greedy {kg}, exact {res['prod']} (missing {res['miss']}, nu {res['nu']}, "
                      f"fan-outs {res['nfan']})", flush=True)
            hist_g[kg] -= 1
            ke = (res["miss"], res["prod"] if res["prod"] is not None else -1)
            hist_g[kg] = hist_g.get(kg, 0) + 1
            hist_e[ke] = hist_e.get(ke, 0) + 1
            if miss == 0 and res["miss"] == 0:
                jk = (cost, res["prod"])
                joint[jk] = joint.get(jk, 0) + 1
                if res["prod"] > cost:
                    worse += 1
            if res["miss"] == 0 and res["prod"] <= A.P + 1:
                print(f"  round {rd}: kept {len(kept)}, D {D}: greedy {kg}, EXACT products {res['prod']} (need "
                      f"{res['need']}, nu {res['nu']}, fan-outs {res['nfan']}, max dim Y(s) {res['maxd']})"
                      f"{'  <= TARGET' if res['prod'] <= A.P else ''} ({time.time() - t0:.0f} s)", flush=True)
                if A.save and res["prod"] <= (A.savemax if A.savemax >= 0 else A.P + 1):
                    saved.append(dict(round=rd, greedy=list(kg), exact=res["prod"], nu=res["nu"],
                                      need=res["need"], phis=[hex(p) for p in phis]))
                    json.dump(dict(load=A.load, n=n, r=r, a=a, seed=A.seed, P=A.P, rseed=A.rseed, spaces=saved),
                              open(A.save, "w"))
        if nk and nk % A.every == 0 and nk != lastp:
            lastp = nk
            te = sorted(k for k in hist_e if k[0] == 0)[:4]
            tg = sorted(k for k in hist_g if k[0] == 0)[:4]
            print(f"  [{nk} priced exactly, round {rd}, {time.time() - t0:.0f} s] exact lowest "
                  f"{[(k[1], hist_e[k]) for k in te]}; greedy lowest {[(k[1], hist_g[k]) for k in tg]}; "
                  f"exact > greedy {worse}", flush=True)
    te = sorted(hist_e.items())[:8]
    tg = sorted(hist_g.items())[:8]
    jt = sorted(joint.items())[:24]
    ngr = sum(hist_g.values())
    print(f"RESULT n {n} r {r} a {a} seed {A.seed}{' ' + A.load if A.load else ''} P {A.P} rseed {A.rseed}: {ngr} "
          f"relation spaces (greedy), {nk} priced exactly; blocks {law}; EXACT lowest (missing, products) {te}; "
          f"GREEDY lowest {tg}; joint (greedy, exact) {jt}; exact > greedy {worse} ({time.time() - t0:.0f} s)",
          flush=True)


if __name__ == "__main__":
    main()
