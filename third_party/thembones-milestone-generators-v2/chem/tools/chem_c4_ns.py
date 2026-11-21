#!/usr/bin/env python3
"""chem_c4_ns.py: builder C4. Pool-free 6-groups at 38 whose products are ALL non-silent on a chosen low coordinate t
(so that a later group's column kill on t keeps the whole pool: the kill needs 75 independent pool columns).

For random 6-output subsets of a loader (residual X0): every hyperplane of W-perp (R59's 'hyp'); for the best ones and
every coordinate t: the fan-outs whose kernel has a plane non-silent on t (2-dim kernel: the plane itself; 3-dim: kept),
priced exactly (C2's Group.lovasz: matroid parity); nu >= 38 = a 38-plan with no product silent on t. Explicit plans
(C2's enum_plans on the filtered kernels) are written for the build.

  python -I chem_c4_ns.py d.json --pool 0-24 --samples 300 --seed 1 --top 6 --out ns.json
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
import chem_c4_subscan as sc  # noqa: E402


def nonsilent_fans(G, phis, t):
    out = []
    for s, Y in G.kernels(phis):
        if len(Y) >= 3:
            out.append((s, Y))
        elif len(Y) == 2 and (t is None or any((y >> t) & 1 for y in (Y[0], Y[1], Y[0] ^ Y[1]))):
            out.append((s, Y))
    return out


def plans_from(G, fans, P, rng, tries, t):
    """explicit plans of P products from the filtered kernels, every plane non-silent on t, features independent and
    holding W."""
    opts = []
    for s, Y in fans:
        pts = []
        for m in range(1, 1 << len(Y)):
            y = 0
            for q in eg.bits(m):
                y ^= Y[q]
            pts.append(y)
        pl, seen = [], set()
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                key = frozenset((pts[i], pts[j], pts[i] ^ pts[j]))
                if key in seen:
                    continue
                seen.add(key)
                if t is None or any((y >> t) & 1 for y in key):
                    pl.append((pts[i], pts[j]))
        if pl:
            opts.append((s, pl))
    out, keys = [], set()
    if len(opts) < P:
        return out
    for _ in range(tries):
        use = sorted(rng.sample(range(len(opts)), P))
        plan = []
        for i in use:
            s, pl = opts[i]
            y1, y2 = pl[rng.randrange(len(pl))]
            plan.append((y1, y2, s))
        key = frozenset((frozenset((y1, y2, y1 ^ y2)), s) for (y1, y2, s) in plan)
        if key in keys:
            continue
        keys.add(key)
        sp = eg.Span()
        if not all(sp.add(G.feat(y1, s)) and sp.add(G.feat(y2, s)) for (y1, y2, s) in plan):
            continue
        if any(sp.reduce(w) for w in G.wb):
            continue
        out.append([list(x) for x in plan])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--pool", default="0-24")
    ap.add_argument("--samples", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--top", type=int, default=6)
    ap.add_argument("--coords", default="", help="only these coordinates t (default: all)")
    ap.add_argument("--tries", type=int, default=3000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = json.load(open(a.dump))
    pool = sc.rng_list(a.pool)
    rng = random.Random(a.seed)
    nl = d["nl"]
    coords = [int(x) for x in a.coords.split(",")] if a.coords else list(range(nl))
    found, seen, t0 = [], set(), time.time()
    for i in range(a.samples):
        outs = tuple(sorted(rng.sample(pool, 6)))
        if outs in seen:
            continue
        seen.add(outs)
        g = sc.group_of(d, list(outs))
        G = eg.Group(g, seed=7)
        rel = eg.nullspace(G.wb, G.N)
        if len(rel) != 15:
            continue
        hs = hy.hyperplanes(G, rel)
        for n, lam in hs[:a.top]:
            if n < 38:
                break
            phis = hy.phis_of(rel, lam)
            for t in coords:
                fans = nonsilent_fans(G, phis, t)
                if len(fans) < 38:
                    continue
                nu = G.lovasz(fans)
                if nu < 38:
                    continue
                plans = plans_from(G, fans, 38, rng, a.tries, t)
                if plans:
                    found.append(dict(outs=list(outs), hyper=lam, t=t, fans=len(fans), plans=plans[:400]))
                    print(f"FOUND {list(outs)} t {t}: {len(fans)} non-silent fan-outs, nu {nu}, {len(plans)} plans "
                          f"[{time.time() - t0:.0f} s]", flush=True)
        if (i + 1) % 50 == 0:
            print(f"  {i + 1} subsets, {len(found)} finds ({time.time() - t0:.0f} s)", flush=True)
    json.dump(found, open(a.out, "w"))
    print(f"DONE {len(seen)} subsets, {len(found)} finds", flush=True)


if __name__ == "__main__":
    main()
