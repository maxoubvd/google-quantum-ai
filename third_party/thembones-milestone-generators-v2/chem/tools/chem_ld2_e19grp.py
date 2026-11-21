#!/usr/bin/env python3
"""chem_ld2_e19grp.py: builder C2 (nineteenth holder, Nov 17). THE REDUCED LOADERS' GROUPS UNDER THE ENGINE'S LAW.

The placed-kill group engine (tools/chem_ld2_e17kill.py, my seventeenth holder's) plans each output group of a
reduced loader on its residual after the kills: cells N = b outputs x a live low coordinates, W = the 76 rows (rank
dim_w), relations W-perp; it draws a RANDOM relation space Phi of dimension e = min(#relations, bonus) and covers W
greedily with products (s, Y), Y inside ker(Phi.s): products >= (N - e) / 2. With e = a - 2 every fan-out has a
kernel of dimension >= 2; a 7-output group with a column kill (N 98, dim_w 76, a 14) reads 43 = (98 - 12) / 2.
42 needs a STRUCTURED Phi of dimension 14 (a random one leaves about 10 of 127 fan-outs with a plane).

  dump:  run the engine on a real loader, write every generic group's W basis, a, b and its own plan (local cells
         o * a + t) to JSON.
  lk:    R56's line-keeping search on one dumped group (keep part of the group's plan; Phi = the relations vanishing
         on W and on the kept products; gap g: random c-dim subspaces of a (c + g)-dim space), priced EXACTLY:
         every fan-out s (2^b - 1) with its kernel Y(s) = ker(Phi.s); supply = sum Y(s) (x) s against the
         annihilator (need = N - c); pairs = linear matroid parity of all planes Y' (x) s, Y' inside Y(s) = rank / 2
         of Lovasz's matrix sum_s (random alternating form on Y(s)) (x) s s^T over GF(2^20) (gf_rank imported
         read-only from tools/r56_chem_engine3_rp_parity.py); products = need - nu when the supply spans.
  stats: the fan-out kernel profile of the group's own relation space family (random e-dim Phi, e = a-2, a-1, a).

  python -I chem_ld2_e19grp.py dump <table.npz> --k 4 --plan "6,6,7:1:0:0,6:1:2:1" [--pad] [--seed 0] --out g.json
  python -I chem_ld2_e19grp.py lk g.json --g 2 --P 42 [--gap 0] [--secs 600] [--rseed 1] [--save s.json]
Run in WSL with ~/qcc/ws/chem_own/venv/bin/python or on a box with ~/qcc/.venv_work/bin/python.
"""
import argparse
import json
import math
import os
import random
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

MB = 20


def bits(x):
    out = []
    while x:
        lo = x & -x
        out.append(lo.bit_length() - 1)
        x ^= lo
    return out


class Span:
    def __init__(self):
        self.piv = {}

    def reduce(self, v):
        while v:
            h = v.bit_length() - 1
            e = self.piv.get(h)
            if e is None:
                return v
            v ^= e
        return 0

    def add(self, v):
        r = self.reduce(v)
        if r:
            self.piv[r.bit_length() - 1] = r
            return True
        return False

    def __len__(self):
        return len(self.piv)


def nullspace(rows, ncols):
    """basis of {x in F2^ncols : <row, x> = 0 for every row}."""
    piv = {}                                   # pivot column -> reduced row
    for r in rows:
        for c, pr in piv.items():
            if r >> c & 1:
                r ^= pr
        if r:
            c = r.bit_length() - 1
            for c2 in list(piv):
                if piv[c2] >> c & 1:
                    piv[c2] ^= r
            piv[c] = r
    free = [j for j in range(ncols) if j not in piv]
    out = []
    for f in free:
        x = 1 << f
        for c, pr in piv.items():
            if pr >> f & 1:
                x |= 1 << c
        out.append(x)
    return out


# ------------------------------------------------------------------ dump (engine side)
def block_set(words, k):
    """the block set V and the padded words exactly as route (b) makes them (tools/chem_ld2_e17gc4.py): c4.vset, then
    vset_rel when relations among the outputs constrain V (relation loaders); returns (V, padded words)."""
    import chem_loaders_core as core
    import chem_ld2_c4 as c4
    import chem_ld2_e17gc4 as gc
    V = c4.vset(words, k)
    Wp = c4.pad_table(words, {k: V})
    an = core.analyze(Wp, k, 4)
    try:
        V2 = gc.vset_rel(words, k, an)
    except Exception:  # noqa: BLE001
        V2 = None
    if V2 is not None and V2 != V:
        V = V2
        Wp = c4.pad_table(words, {k: V})
    return V, Wp


def cmd_dump(a):
    import chem_loaders_core as core
    import chem_ld2_e17kill as kl
    core.J, core.BROT = 5, 7
    words = np.asarray(np.load(a.table)["angle_words"], np.int64)
    V, kap = None, 0
    if a.pad:
        V, words0 = block_set(words, a.k)
        words = words0
    an = core.analyze(words, a.k, 4)
    if a.pad:
        import chem_ld2_c4 as c4
        kap = c4.kappa_wires(an, V)
    if a.plan == "auto":                      # the placed-kill plan of this loader's number of generic outputs
        a.plan = kl.PLANS_BY_N[sum(1 for x in an.out_level if x > 2)][0]
    plan = kl.parse_plan(a.plan)
    groups = []

    class Probe(kl.KillDecomp):
        def apply_kills(self, outs, coords, dl, d, shift=0):
            live = super().apply_kills(outs, coords, dl, d, shift)
            if live is not None:
                self._live = live
            return live

        def _select(self, keys, wb, dim_w):
            self._wb = list(wb)
            return super()._select(keys, wb, dim_w)

        def solve_group(self, outs, coords, dl=0, d=0, tag=""):
            self._live, self._wb = coords, None
            n0 = len(self.prods)
            made = super().solve_group(outs, coords, dl, d, tag)
            if tag == "generic" and self._wb is not None:
                lc = bits(self._live)
                loc = {t: i for i, t in enumerate(lc)}
                oi = {o: i for i, o in enumerate(outs)}
                prods = []
                for p in self.prods[n0:]:
                    y1 = sum(1 << loc[t] for t in bits(p[0]))
                    y2 = sum(1 << loc[t] for t in bits(p[2])) if p[2] else 0
                    s = sum(1 << oi[o] for o in bits(p[4]) if o in oi)
                    prods.append([y1, y2, s])
                groups.append(dict(outs=list(outs), coords=lc, a=len(lc), b=len(outs), dim_w=len(self._wb),
                                   wb=[hex(v) for v in self._wb], prods=prods, made=made, log=self.log[-1]))
            return made

    kd = Probe(an, a.seed, max(2, a.effort), plan, [int(x) for x in a.avoid.split(",") if x.strip()], V, kap,
               a.tries if a.pad else 1)
    prods = kd.run()
    out = dict(table=a.table, k=a.k, plan=a.plan, seed=a.seed, pad=a.pad, products=len(prods), infra=an.infra,
               nb=int(an.X0.shape[1]), nr=int(an.X0.shape[0]), groups=groups)
    json.dump(out, open(a.out, "w"), indent=0)
    print(f"k {a.k}: {an.X0.shape[1]} generic outputs, {len(prods)} products (+{an.infra} infra); groups: "
          + " | ".join(f"b={g['b']} a={g['a']} dimW={g['dim_w']} -> {g['made']}" for g in groups), flush=True)


# ------------------------------------------------------------------ exact pricer for one group
class Group:
    def __init__(self, g, seed=1):
        self.a, self.b = g["a"], g["b"]
        self.N = self.a * self.b
        self.full_a = (1 << self.a) - 1
        self.wb = [int(v, 16) for v in g["wb"]]
        self.prods = [tuple(p) for p in g["prods"]]
        self.gen = np.random.default_rng(seed)
        self.EXP, self.LOG, self.gf_rank = None, None, None
        import r56_chem_engine3_rp_parity as rp   # read-only import: GF(2^20) tables and rank
        self.gf_rank = rp.gf_rank
        self.svec = {}
        for s in range(1, 1 << self.b):
            self.svec[s] = np.array([(s >> o) & 1 for o in range(self.b)], np.int64)

    def feat(self, y, s):
        v = 0
        for o in bits(s):
            v |= y << (o * self.a)
        return v

    def restrict(self, phis, s):
        a, fa = self.a, self.full_a
        out = []
        for ph in phis:
            v = 0
            for o in bits(s):
                v ^= (ph >> (o * a)) & fa
            if v:
                out.append(v)
        return out

    def kernels(self, phis):
        fans = []
        for s in range(1, 1 << self.b):
            Y = nullspace(self.restrict(phis, s), self.a)
            if Y:
                fans.append((s, Y))
        return fans

    def lovasz(self, fans):
        a, b, N = self.a, self.b, self.N
        M = np.zeros((N, N), np.int64)
        for s, Y in fans:
            d = len(Y)
            if d < 2:
                continue
            S = np.outer(self.svec[s], self.svec[s])
            yv = [np.array([(y >> t) & 1 for t in range(a)], np.int64) for y in Y]
            Asum = np.zeros((a, a), np.int64)
            for i in range(d):
                for j in range(i + 1, d):
                    x = int(self.gen.integers(1, 1 << MB))
                    pat = np.outer(yv[i], yv[j]) ^ np.outer(yv[j], yv[i])
                    Asum ^= pat * x
            M ^= np.kron(S, Asum)
        rk = self.gf_rank(M)
        assert rk % 2 == 0, rk
        return rk // 2

    def price(self, phis, reps=1):
        c = len(phis)
        need = self.N - c
        fans = self.kernels(phis)
        sp = Span()
        for s, Y in fans:
            for y in Y:
                sp.add(self.feat(y, s))
        miss = need - len(sp)
        planes = [(s, Y) for s, Y in fans if len(Y) >= 2]
        nu = max(self.lovasz(planes) for _ in range(reps)) if planes else 0
        prof = {}
        for s, Y in fans:
            prof[len(Y)] = prof.get(len(Y), 0) + 1
        return dict(need=need, miss=miss, nu=nu, prod=(need - nu) if miss <= 0 else None, prof=prof)


def keep_space(G, c, gap, rng, maxfail=6):
    """R56's keep step on this group's own plan: grow the kept set while N - rank(W + K) >= c + gap."""
    pairs = [p for p in G.prods if p[1]]
    rng.shuffle(pairs)
    sp = Span()
    for w in G.wb:
        sp.add(w)
    kept, fails = [], 0
    for (y1, y2, s) in pairs:
        save = dict(sp.piv)
        sp.add(G.feat(y1, s))
        sp.add(G.feat(y2, s))
        if G.N - len(sp) < c + gap:
            sp.piv = save
            fails += 1
            if fails > maxfail:
                break
            continue
        kept.append((y1, y2, s))
    basis = list(sp.piv.values())
    return nullspace(basis, G.N), kept


def sub_space(D, c, rng):
    """a random c-dim subspace of the span of D (list of ints)."""
    out, sp = [], Span()
    while len(out) < c:
        v = 0
        for d in D:
            if rng.random() < 0.5:
                v ^= d
        if v and sp.add(v):
            out.append(v)
    return out


def cmd_lk(a):
    d = json.load(open(a.dump))
    g = d["groups"][a.g]
    G = Group(g, seed=a.rseed + 7)
    rng = random.Random(a.rseed)
    c = G.N - 2 * a.P
    print(f"group {a.g}: b {G.b}, a {G.a}, N {G.N}, dim W {g['dim_w']}, own plan {len(G.prods)} products "
          f"({sum(1 for p in G.prods if p[1])} pairs); target P {a.P} -> c {c} relations; W-perp "
          f"{G.N - g['dim_w']}; gap {a.gap}", flush=True)
    t0, tried, hist, found = time.time(), 0, {}, []
    rounds = 0
    while time.time() - t0 < a.secs and rounds < a.rounds:
        rounds += 1
        D, kept = keep_space(G, c, a.gap, rng)
        if len(D) < c:
            continue
        cands = [D] if len(D) == c else [sub_space(D, c, rng) for _ in range(a.samples)]
        for phis in cands:
            res = G.price(phis)
            tried += 1
            key = res["prod"] if res["miss"] <= 0 else -res["miss"]
            hist[key] = hist.get(key, 0) + 1
            if res["miss"] <= 0 and res["prod"] <= a.P:
                found.append(dict(kept=len(kept), phis=[hex(v) for v in phis], res=res))
                print(f"FOUND {res['prod']} (kept {len(kept)}, nu {res['nu']}, profile {res['prof']}) after "
                      f"{tried} spaces ({time.time() - t0:.0f} s)", flush=True)
                if a.save:
                    json.dump(found, open(a.save, "w"))
                if len(found) >= a.stop:
                    break
        if len(found) >= a.stop:
            break
        if rounds % a.every == 0:
            top = sorted(hist.items(), key=lambda kv: (kv[0] < 0, abs(kv[0])))[:6]
            print(f"round {rounds}: {tried} spaces, kept {len(kept)}, D {len(D)}; readings {top} "
                  f"({time.time() - t0:.0f} s)", flush=True)
    top = sorted(hist.items(), key=lambda kv: (kv[0] < 0, abs(kv[0])))[:10]
    print(f"END: {tried} spaces in {rounds} rounds ({time.time() - t0:.0f} s); found {len(found)} at <= {a.P}; "
          f"readings (products: count; negative = supply short by) {top}", flush=True)


def kernel_dims_np(M, a):
    """dims of the kernels of S matrices at once over GF(2): M is (S, c) uint32, row r of matrix s = an a-bit int."""
    M = M.copy()
    S, c = M.shape
    idx = np.arange(S)
    rank = np.zeros(S, np.int32)
    active = np.ones((S, c), bool)
    for col in range(a):
        bit = ((M >> col) & 1).astype(bool)
        cand = bit & active
        has = cand.any(axis=1)
        piv = np.argmax(cand, axis=1)
        prow = M[idx, piv]
        elim = bit & has[:, None]
        elim[idx, piv] = False
        M = np.where(elim, M ^ prow[:, None], M)
        active[idx[has], piv[has]] = False
        rank += has
    return a - rank


def cmd_anneal_fast(a, G, rel, c, rng):
    """the same annealer with the surrogate vectorized: rows[s, i] = Phi_i restricted by s, updated per move."""
    R = len(rel)
    S = (1 << G.b) - 1
    RJ = np.zeros((R, S), np.uint32)
    for j, r in enumerate(rel):
        for s in range(1, 1 << G.b):
            v = 0
            for o in bits(s):
                v ^= (r >> (o * G.a)) & G.full_a
            RJ[j, s - 1] = v
    C = [rng.getrandbits(R) for _ in range(c)]
    rows = np.zeros((S, c), np.uint32)
    for i, cv in enumerate(C):
        for j in bits(cv):
            rows[:, i] ^= RJ[j]

    def crank():
        sp = Span()
        return sum(1 for v in C if sp.add(v))

    def surr():
        kd = kernel_dims_np(rows, G.a)
        m = kd >= 2
        return int(m.sum()), int(kd[m].sum())

    def phis_of():
        out = []
        for cv in C:
            v = 0
            for j in bits(cv):
                v ^= rel[j]
            out.append(v)
        return out

    cur, curd = surr() if crank() == c else (-1, 0)
    best = (cur, curd)
    priced, found = set(), []
    t0, step, T = time.time(), 0, a.T0
    print(f"group {a.g}: fast annealer, W-perp {R}, target P {a.P} -> c {c}; start planes {cur}", flush=True)
    while time.time() - t0 < a.secs:
        step += 1
        i, j = rng.randrange(c), rng.randrange(R)
        C[i] ^= 1 << j
        rows[:, i] ^= RJ[j]
        if crank() < c:
            new, newd = -1, 0
        else:
            new, newd = surr()
        dE = (new - cur) + 0.05 * (newd - curd)
        if new >= 0 and (dE >= 0 or rng.random() < math.exp(dE / max(T, 1e-6))):
            cur, curd = new, newd
            improved = (cur, curd) > best
            if improved:
                best = (cur, curd)
            if cur >= a.price and improved:
                key = tuple(sorted(C))
                if key not in priced:
                    priced.add(key)
                    phis = phis_of()
                    res = G.price(phis)
                    if res["miss"] <= 0 and res["prod"] <= a.P:
                        found.append(dict(phis=[hex(v) for v in phis], res=res))
                        print(f"step {step}: planes {cur} -> FOUND {res['prod']} (nu {res['nu']}, profile "
                              f"{res['prof']}) ({time.time() - t0:.0f} s)", flush=True)
                        if a.save:
                            json.dump(found, open(a.save + ".tmp", "w"))
                            os.replace(a.save + ".tmp", a.save)
                        if not a.multi:
                            print("FOUND", flush=True)
                            return
                        for _ in range(a.kick):             # kick away from the found space
                            ii, jj = rng.randrange(c), rng.randrange(R)
                            C[ii] ^= 1 << jj
                            rows[:, ii] ^= RJ[jj]
                        cur, curd = surr() if crank() == c else (-1, 0)
                        best = (cur, curd)                  # the next find must climb again
        else:
            C[i] ^= 1 << j
            rows[:, i] ^= RJ[j]
        T = a.T0 * (1 - (time.time() - t0) / a.secs) + 0.05
        if step % a.every == 0:
            print(f"step {step}: cur planes {cur} (dims {curd}), best {best[0]}, T {T:.2f}, found {len(found)} "
                  f"({time.time() - t0:.0f} s)", flush=True)
    print(f"END: {step} steps, best planes {best[0]}, found {len(found)}", flush=True)


def cmd_anneal(a):
    """the dual search on one group: anneal Phi (c relations, each a combination of a basis of W-perp) on the
    surrogate 'fan-outs with a kernel of dimension >= 2' (each gives a plane; 42 needed for 42 products), every
    candidate that reaches --price planes priced exactly (supply + Lovasz). Moves: flip one coefficient."""
    d = json.load(open(a.dump))
    g = d["groups"][a.g]
    G = Group(g, seed=a.rseed + 7)
    rng = random.Random(a.rseed)
    rel = nullspace(G.wb, G.N)
    R = len(rel)
    c = G.N - 2 * a.P
    if not a.slow:
        return cmd_anneal_fast(a, G, rel, c, rng)

    def phis_of(C):
        out = []
        for row in C:
            v = 0
            for j in bits(row):
                v ^= rel[j]
            out.append(v)
        return out

    def surrogate(phis):
        if Span_rank(phis) < c:
            return -1, 0
        pl, dsum = 0, 0
        for s in range(1, 1 << G.b):
            k = len(nullspace(G.restrict(phis, s), G.a))
            if k >= 2:
                pl += 1
                dsum += k
        return pl, dsum

    def Span_rank(vs):
        sp = Span()
        return sum(1 for v in vs if sp.add(v))

    start = None
    if a.start:
        st = json.load(open(a.start))
        cand = st[0]["phis"] if isinstance(st, list) else st["phis"]
        phis0 = [int(x, 16) for x in cand]
        # coefficients over rel: solve each relation in the basis rel
        sp = {}
        basis = []
        for j, r in enumerate(rel):
            basis.append((r, 1 << j))
        C = []
        for ph in phis0:
            v, cv = ph, 0
            piv = {}
            for r, tag in basis:
                rr, tt = r, tag
                for h, (pr, pt) in piv.items():
                    if rr >> h & 1:
                        rr ^= pr
                        tt ^= pt
                if rr:
                    piv[rr.bit_length() - 1] = (rr, tt)
            for h in sorted(piv, reverse=True):
                if v >> h & 1:
                    v ^= piv[h][0]
                    cv ^= piv[h][1]
            assert v == 0
            C.append(cv)
        start = C
    C = start or [rng.getrandbits(R) for _ in range(c)]
    cur, curd = surrogate(phis_of(C))
    best = (cur, curd, list(C))
    t0 = time.time()
    T = a.T0
    step = 0
    priced = {}
    print(f"group {a.g}: b {G.b}, a {G.a}, N {G.N}, W-perp {R}, target P {a.P} -> c {c}; start planes {cur}",
          flush=True)
    while time.time() - t0 < a.secs:
        step += 1
        i, j = rng.randrange(c), rng.randrange(R)
        C[i] ^= 1 << j
        new, newd = surrogate(phis_of(C))
        dE = (new - cur) + 0.05 * (newd - curd)
        if new >= 0 and (dE >= 0 or rng.random() < math.exp(dE / max(T, 1e-6))):
            cur, curd = new, newd
            if (cur, curd) > best[:2]:
                best = (cur, curd, list(C))
                if cur >= a.price:
                    key = tuple(sorted(C))
                    if key not in priced:
                        res = G.price(phis_of(C))
                        priced[key] = res
                        print(f"step {step}: planes {cur} -> exact {res['prod'] if res['miss'] <= 0 else 'short ' + str(res['miss'])} "
                              f"(nu {res['nu']}, profile {res['prof']}) ({time.time() - t0:.0f} s)", flush=True)
                        if res["miss"] <= 0 and res["prod"] <= a.P:
                            if a.save:
                                json.dump([dict(phis=[hex(v) for v in phis_of(C)], res=res)], open(a.save, "w"))
                            print(f"FOUND {res['prod']}", flush=True)
                            return
        else:
            C[i] ^= 1 << j
        T = a.T0 * (1 - (time.time() - t0) / a.secs) + 0.05
        if step % a.every == 0:
            print(f"step {step}: cur planes {cur} (dims {curd}), best {best[0]} ({best[1]}), T {T:.2f} "
                  f"({time.time() - t0:.0f} s)", flush=True)
    res = G.price(phis_of(best[2]))
    print(f"END: best planes {best[0]}; exact {res['prod'] if res['miss'] <= 0 else 'short ' + str(res['miss'])} "
          f"(nu {res['nu']}, profile {res['prof']}); {step} steps", flush=True)


def extract_plan(G, phis, rng):
    """an explicit plan from a relation space that prices at P = need - nu: every plane of every fan-out kernel is a
    line; lines are deleted (random order) while the Lovasz rank keeps nu; then the kept lines' features must span
    the annihilator (need = 2 nu here) and hold W. Returns [(y1, y2, s)]."""
    res = G.price(phis, reps=2)
    fans = [(s, Y) for s, Y in G.kernels(phis) if len(Y) >= 2]
    lines = []
    for s, Y in fans:
        d = len(Y)
        pts = []
        for m in range(1, 1 << d):
            y = 0
            for q in bits(m):
                y ^= Y[q]
            pts.append(y)
        seen = set()
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                key = frozenset((pts[i], pts[j], pts[i] ^ pts[j]))
                if key not in seen:
                    seen.add(key)
                    lines.append((s, pts[i], pts[j]))
    nu = res["nu"]
    order = list(range(len(lines)))
    rng.shuffle(order)
    keep = set(order)
    for i in order:
        if len(keep) == nu:
            break
        trial = keep - {i}
        if G.lovasz([(lines[j][0], [lines[j][1], lines[j][2]]) for j in sorted(trial)]) == nu:
            keep = trial
    plan = [lines[j] for j in sorted(keep)]
    sp = Span()
    ok = all(sp.add(G.feat(y1, s)) and sp.add(G.feat(y2, s)) for s, y1, y2 in plan)
    holds = ok and all(sp.reduce(w) == 0 for w in G.wb)
    return [(y1, y2, s) for s, y1, y2 in plan], dict(nu=nu, need=res["need"], lines=len(lines), independent=ok,
                                                       holds_W=holds)


def enum_plans(G, phis, P, rng, limit):
    """every plan of P products from one relation space whose fan-outs give at most one plane each: drop
    (#fan-outs - P) fan-outs, pick one plane in each kernel of dimension >= 3; keep the candidates whose 2P features
    are independent and hold W. Random order, at most `limit` tries."""
    fans = [(s, Y) for s, Y in G.kernels(phis) if len(Y) >= 2]
    opts = []
    for s, Y in fans:
        d = len(Y)
        pts = []
        for m in range(1, 1 << d):
            y = 0
            for q in bits(m):
                y ^= Y[q]
            pts.append(y)
        pl, seen = [], set()
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                key = frozenset((pts[i], pts[j], pts[i] ^ pts[j]))
                if key not in seen:
                    seen.add(key)
                    pl.append((pts[i], pts[j]))
        opts.append((s, pl))
    out, keys = [], set()
    nf = len(opts)
    for _ in range(limit):
        use = sorted(rng.sample(range(nf), P)) if nf > P else list(range(nf))
        plan = []
        for i in use:
            s, pl = opts[i]
            y1, y2 = pl[rng.randrange(len(pl))]
            plan.append((y1, y2, s))
        key = frozenset((frozenset((y1, y2, y1 ^ y2)), s) for (y1, y2, s) in plan)
        if key in keys:
            continue
        keys.add(key)
        sp = Span()
        if not all(sp.add(G.feat(y1, s)) and sp.add(G.feat(y2, s)) for (y1, y2, s) in plan):
            continue
        if any(sp.reduce(w) for w in G.wb):
            continue
        out.append(plan)
    return out, len(opts), sum(1 for _, pl in opts if len(pl) > 1)


def cmd_plan(a):
    d = json.load(open(a.dump))
    g = d["groups"][a.g]
    G = Group(g, seed=a.rseed + 11)
    rng = random.Random(a.rseed)
    found = json.load(open(a.spaces))
    out = []
    for f in found:
        phis = [int(x, 16) for x in f["phis"]]
        P = f["res"]["prod"] if "res" in f else G.N // 2 - len(phis) // 2
        plans, nf, nmulti = enum_plans(G, phis, P, rng, a.tries)
        print(f"space: {nf} fan-outs with a plane ({nmulti} with several); {len(plans)} distinct valid plans of {P} "
              f"in {a.tries} tries", flush=True)
        out += [dict(planes=[[y1, y2, s] for (y1, y2, s) in pl]) for pl in plans]
    json.dump(out, open(a.out, "w"))
    print(f"{len(out)} plans written to {a.out}", flush=True)


def cmd_build(a):
    """the placed-kill engine with the closure, generic group number --g planned by the given plan file (planes in
    the group's local coordinates, from `plan`); the other groups by the engine; emit_c4 (68-AND, checked)."""
    import chem_loaders_core as core
    import chem_ld2_e17kill as kl
    import chem_ld2_c4 as c4
    core.J, core.BROT = 5, 7
    words = np.asarray(np.load(a.table)["angle_words"], np.int64)
    V, words = block_set(words, a.k)
    an = core.analyze(words, a.k, 4)
    kap = c4.kappa_wires(an, V)
    plans = json.load(open(a.plans))
    # injected generic groups: {index: list of plans}; a second group (--g2 / --inj2) for the relation loaders' two
    # 7-groups (its plans are made on the residual left by the first group's closed plan: --fix1)
    injs = {a.g: [[tuple(x) for x in pl["planes"]] for pl in plans[:a.maxplans]]}
    if a.fix1:
        f1 = json.load(open(a.fix1))
        f1 = f1 if "planes" in f1 else f1[str(a.g)]           # a --saveclosed file works too
        injs[a.g] = [[tuple(x) for x in f1["planes"]]]
    if a.inj2:
        injs[a.g2] = [[tuple(x) for x in pl["planes"]] for pl in json.load(open(a.inj2))[:a.maxplans]]
    dumped = []

    class PlanKill(kl.KillDecomp):
        def apply_kills(self, outs, coords, dl, d, shift=0):
            live = super().apply_kills(outs, coords, dl, d, shift)
            if live is not None:
                self._live = live
                if self._cur is not None:             # a new attempt at an injected group: its next plan
                    pi = self._pidx.get(self._cur, -1) + 1
                    self._pidx[self._cur] = pi
                    lst = injs[self._cur]
                    self._planes = lst[pi % len(lst)]
            return live

        def _tau_bad(self, wmask):
            bad = super()._tau_bad(wmask)
            if self._cur is not None:
                self._taus = getattr(self, "_taus", {})
                self._taus.setdefault((self._cur, bad), []).append(self._pidx.get(self._cur, 0))
            return bad

        def solve_group(self, outs, coords, dl=0, d=0, tag=""):
            self._live, self._cur, self._wb = coords, None, None
            self._inject = None
            n0 = len(self.prods)
            if tag == "generic":
                self._gi = getattr(self, "_gi", -1) + 1
                if self._gi in injs:
                    self._cur = self._gi
                    self._inject = (outs, coords)
            made = super().solve_group(outs, coords, dl, d, tag)
            if tag == "generic" and self._wb is not None:
                lc = bits(self._live)
                loc = {t: i for i, t in enumerate(lc)}
                oi = {o: i for i, o in enumerate(outs)}
                pr = []
                for p in self.prods[n0:]:
                    y1 = sum(1 << loc[t] for t in bits(p[0]))
                    y2 = sum(1 << loc[t] for t in bits(p[2])) if p[2] else 0
                    pr.append([y1, y2, sum(1 << oi[o] for o in bits(p[4]) if o in oi)])
                dumped.append(dict(outs=list(outs), coords=lc, a=len(lc), b=len(outs), dim_w=len(self._wb),
                                   wb=[hex(v) for v in self._wb], prods=pr, made=made, log=self.log[-1],
                                   injected=self._cur))
            self._cur = None
            return made

        def _select(self, keys, wb, dim_w):
            self._wb = list(wb)
            if self._cur is None:
                return super()._select(keys, wb, dim_w)
            outs, coords = self._inject
            lc = bits(self._live)
            A = len(lc)
            if self._cur == a.g:
                assert A == self._pa and len(outs) == self._pb, (A, len(outs))
            ylocal = []
            for ym in self.an.ymasks:
                m = sum(1 << ci for ci, t in enumerate(lc) if ym >> t & 1)
                if m and m not in ylocal:
                    ylocal.append(m)
            chosen = []
            F = Span()
            why = None
            for (y1, y2, s) in self._planes:
                ps = [p for p, ym in enumerate(ylocal) if not ((y1 | y2) & ~ym)]
                if not ps:
                    why = f"no phase holds the plane (phases {len(ylocal)})"
                    break
                key = (s, ps[0])
                for y in (y1, y2):
                    flat = sum(y << (gi * A) for gi in bits(s))
                    if not F.add(flat):
                        why = "dependent features"
                    chosen.append((key, y, flat))
            if why is None and any(F.reduce(w) for w in wb):
                why = f"W not inside the plan's span (dim W {dim_w}, {sum(1 for w in wb if F.reduce(w))} rows out)"
            self._inj_calls = getattr(self, "_inj_calls", 0) + 1
            if why is not None:
                if self._inj_calls == 1:
                    print(f"  injection refused: {why}", flush=True)
                return None
            if self._inj_calls == 1:
                print(f"  injection accepted: {len(self._planes)} products", flush=True)
            return len(self._planes), chosen

    best = None
    if a.plan == "auto":
        a.plan = kl.PLANS_BY_N[sum(1 for x in an.out_level if x > 2)][0]
    allp = injs[a.g]
    ntry = max(a.tries, max(len(v) for v in injs.values()))
    kd = PlanKill(an, a.seed, max(2, a.effort), kl.parse_plan(a.plan), [0, 1], V, kap, ntry)
    kd._planes, kd._inject, kd._cur, kd._pidx = allp[0], None, None, {}
    kd._pa, kd._pb = a.a, a.b
    try:
        prods = kd.run()
        text, st = c4.emit_c4(an, prods, V)
        line = " | ".join(x.split("generic ")[-1] for x in kd.log if "generic" in x)
        used = {g: kd._pidx.get(g, 0) % len(injs[g]) for g in injs}
        print(f"{len(allp)} plans offered; closed with plan {used[a.g]} "
              f"(attempts {kd.retries}): {len(prods)} products, emit_c4 {st['toffoli']} (closed, checked); {line}",
              flush=True)
        best = (st["toffoli"], used[a.g], text)
        if a.saveclosed:
            json.dump({str(g): dict(planes=[list(p) for p in injs[g][used[g]]]) for g in injs},
                      open(a.saveclosed, "w"))
    except (RuntimeError, AssertionError, ValueError) as exc:
        print(f"{len(allp)} plans offered: FAILED {type(exc).__name__}: {str(exc)[:160]}", flush=True)
    if a.dumpg:
        json.dump(dict(table=a.table, k=a.k, plan=a.plan, seed=a.seed, groups=dumped), open(a.dumpg, "w"))
    taus = getattr(kd, "_taus", {})
    print(f"injected group: {len(taus)} distinct mismatch values over "
          f"{len(set((g, i) for (g, _), v in taus.items() for i in v))} plans; weights "
          f"{sorted(bin(t).count('1') for (_, t) in taus)[:12]}", flush=True)
    if best and a.emit:
        os.makedirs(a.emit, exist_ok=True)
        with open(os.path.join(a.emit, f"angle_k_{a.k}.kmx"), "w") as fh:
            fh.write(best[2])
        print(f"BEST {best[0]} (plan {best[1]}) -> {a.emit}/angle_k_{a.k}.kmx", flush=True)


def cmd_stats(a):
    d = json.load(open(a.dump))
    g = d["groups"][a.g]
    G = Group(g, seed=a.rseed + 7)
    rng = random.Random(a.rseed)
    rel = nullspace(G.wb, G.N)
    print(f"group {a.g}: b {G.b}, a {G.a}, N {G.N}, dim W {g['dim_w']}, relations {len(rel)}, own plan "
          f"{len(G.prods)}", flush=True)
    for e in [G.a - 2, G.a - 1, G.a, G.a + 2]:
        if e > len(rel):
            continue
        rs = []
        for _ in range(a.samples):
            phis = sub_space(rel, e, rng)
            res = G.price(phis)
            rs.append((res["prod"] if res["miss"] <= 0 else -res["miss"], res["prof"]))
        vals = sorted(r[0] for r in rs)
        prof = {}
        for _, p in rs:
            for k, v in p.items():
                prof[k] = prof.get(k, 0) + v / len(rs)
        print(f"random Phi of dim {e}: readings {vals[:8]}...; mean kernel profile "
              f"{ {k: round(v, 1) for k, v in sorted(prof.items())} }", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("dump")
    p.add_argument("table")
    p.add_argument("--k", type=int, required=True)
    p.add_argument("--plan", required=True)
    p.add_argument("--avoid", default="")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--effort", type=int, default=3)
    p.add_argument("--pad", action="store_true")
    p.add_argument("--tries", type=int, default=400)
    p.add_argument("--out", required=True)
    p = sub.add_parser("lk")
    p.add_argument("dump")
    p.add_argument("--g", type=int, required=True)
    p.add_argument("--P", type=int, required=True)
    p.add_argument("--gap", type=int, default=0)
    p.add_argument("--samples", type=int, default=8)
    p.add_argument("--secs", type=float, default=600)
    p.add_argument("--rounds", type=int, default=10 ** 9)
    p.add_argument("--rseed", type=int, default=1)
    p.add_argument("--every", type=int, default=200)
    p.add_argument("--stop", type=int, default=3)
    p.add_argument("--save", default="")
    p = sub.add_parser("stats")
    p.add_argument("dump")
    p.add_argument("--g", type=int, required=True)
    p.add_argument("--samples", type=int, default=20)
    p.add_argument("--rseed", type=int, default=1)
    p = sub.add_parser("anneal")
    p.add_argument("dump")
    p.add_argument("--g", type=int, required=True)
    p.add_argument("--P", type=int, required=True)
    p.add_argument("--secs", type=float, default=600)
    p.add_argument("--T0", type=float, default=1.5)
    p.add_argument("--price", type=int, default=30, help="price exactly from this many planes on")
    p.add_argument("--rseed", type=int, default=1)
    p.add_argument("--every", type=int, default=2000)
    p.add_argument("--start", default="", help="a saved relation space (lk --save) to start from")
    p.add_argument("--save", default="")
    p.add_argument("--slow", action="store_true", help="the first (pure Python) surrogate")
    p.add_argument("--multi", action="store_true", help="keep annealing after a find (kick, append to --save)")
    p.add_argument("--kick", type=int, default=6)
    p = sub.add_parser("plan")
    p.add_argument("dump")
    p.add_argument("spaces", help="json list of found relation spaces (anneal / lk --save)")
    p.add_argument("--g", type=int, required=True)
    p.add_argument("--tries", type=int, default=20000)
    p.add_argument("--rseed", type=int, default=1)
    p.add_argument("--out", required=True)
    p = sub.add_parser("build")
    p.add_argument("table")
    p.add_argument("plans")
    p.add_argument("--k", type=int, required=True)
    p.add_argument("--plan", required=True)
    p.add_argument("--g", type=int, required=True, help="generic group index planned from the file (0-based)")
    p.add_argument("--a", type=int, required=True, help="the group's live low coordinates (check)")
    p.add_argument("--b", type=int, required=True, help="the group's outputs (check)")
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--effort", type=int, default=3)
    p.add_argument("--tries", type=int, default=400)
    p.add_argument("--maxplans", type=int, default=100)
    p.add_argument("--emit", default=None)
    p.add_argument("--fix1", default="", help="one plan file (dict with planes) fixing the first injected group")
    p.add_argument("--inj2", default="", help="plans for a second injected group (--g2)")
    p.add_argument("--g2", type=int, default=3)
    p.add_argument("--dumpg", default="", help="write every generic group's W basis and plan as seen in this build")
    p.add_argument("--saveclosed", default="", help="write the closing plan of every injected group")
    a = ap.parse_args()
    {"dump": cmd_dump, "lk": cmd_lk, "stats": cmd_stats, "anneal": cmd_anneal, "plan": cmd_plan,
     "build": cmd_build}[a.cmd](a)


if __name__ == "__main__":
    main()
