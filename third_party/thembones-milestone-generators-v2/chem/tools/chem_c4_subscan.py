#!/usr/bin/env python3
"""chem_c4_subscan.py: builder C4. Scan 6-output subsets of a loader's generic outputs (residual = X0: groups without
kills) for a pool-free 6-group at 38 products: every hyperplane of the subset's W-perp counted (designer R59's 'hyp'
test), the best ones priced exactly (builder C2's Group.price: kernels + Lovasz matroid parity).

  python -I chem_c4_subscan.py d.json --pool 6-24 --exclude 0-5 --samples 200 --seed 1 --out found.json
"""
import argparse
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_ld2_e19grp as eg  # noqa: E402
import chem_c4_hyp as hy  # noqa: E402


def rng_list(s):
    a, b = (int(x) for x in s.split("-"))
    return list(range(a, b + 1))


def group_of(d, outs):
    nr, nl = d["nr"], d["nl"]
    X0 = d["X0"]
    rows = []
    for M in range(nr):
        v = 0
        for gi, o in enumerate(outs):
            x = X0[M][o]
            for ci in range(nl):
                if (x >> ci) & 1:
                    v |= 1 << (gi * nl + ci)
        rows.append(v)
    sp, wb = eg.Span(), []
    for v in rows:
        if sp.add(v):
            wb.append(v)
    return dict(outs=outs, coords=list(range(nl)), a=nl, b=len(outs), dim_w=len(wb), wb=[hex(v) for v in wb],
                prods=[], made=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--pool", default="6-24")
    ap.add_argument("--samples", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--top", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = json.load(open(a.dump))
    pool = rng_list(a.pool)
    rng = random.Random(a.seed)
    found, t0 = [], time.time()
    seen = set()
    for i in range(a.samples):
        outs = tuple(sorted(rng.sample(pool, 6)))
        if outs in seen:
            continue
        seen.add(outs)
        g = group_of(d, list(outs))
        G = eg.Group(g, seed=7)
        rel = eg.nullspace(G.wb, G.N)
        if len(rel) != 15:
            continue
        hs = hy.hyperplanes(G, rel)
        best_n = hs[0][0]
        res_best = None
        if best_n >= 38:
            for n, lam in hs[:a.top]:
                if n < 38:
                    break
                res = G.price(hy.phis_of(rel, lam), reps=2)
                if res["prod"] is not None and res["prod"] <= 38:
                    res_best = (lam, res["prod"])
                    break
        if res_best:
            found.append(dict(outs=list(outs), hyper=res_best[0], prod=res_best[1], fanouts=best_n))
            print(f"FOUND {list(outs)}: 38 (hyperplane {res_best[0]:#x}, {best_n} fan-outs) [{time.time() - t0:.0f} s]",
                  flush=True)
        if (i + 1) % 25 == 0:
            print(f"  {i + 1} subsets, {len(found)} at 38 ({time.time() - t0:.0f} s)", flush=True)
    json.dump(found, open(a.out, "w"))
    print(f"DONE {len(seen)} subsets, {len(found)} at 38", flush=True)


if __name__ == "__main__":
    main()
