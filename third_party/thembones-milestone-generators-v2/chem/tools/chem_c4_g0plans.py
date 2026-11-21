#!/usr/bin/env python3
"""chem_c4_g0plans.py: builder C4. Plans at 38 for group G of a C2-format dump (e.g. design seat R65's --dump), every
product non-silent on coordinate t: every hyperplane of W-perp, kernels filtered for silence, exact matroid parity,
explicit plans (chem_c4_ns.plans_from). Output: [{planes: [[y1, y2, s], ...]}, ...] for R65's / C2's --inj.

  python -I chem_c4_g0plans.py d.json --g 0 --t 14 --max 1500 --out plans.json
"""
import argparse
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_ld2_e19grp as eg  # noqa: E402
import chem_c4_hyp as hy  # noqa: E402
import chem_c4_ns as ns  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--g", type=int, default=0)
    ap.add_argument("--t", type=int, default=14)
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("--tries", type=int, default=4000)
    ap.add_argument("--max", type=int, default=1500)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = json.load(open(a.dump))
    g = d["groups"][a.g]
    G = eg.Group(g, seed=7)
    rel = eg.nullspace(G.wb, G.N)
    rng = random.Random(1)
    t = g["coords"].index(a.t) if a.t in g["coords"] else None
    relax = a.t < 0                            # holder 2: --t -1 = no silence filter (owners: dead supply pays the kill)
    out = []
    if len(rel) == 15 and (t is not None or relax):
        for n, lam in hy.hyperplanes(G, rel)[:a.top]:
            if n < 38:
                break
            fans = ns.nonsilent_fans(G, hy.phis_of(rel, lam), t)
            if len(fans) < 38 or G.lovasz(fans) < 38:
                continue
            out += [dict(planes=pl) for pl in ns.plans_from(G, fans, 38, rng, a.tries, t)]
            if len(out) >= a.max:
                break
    json.dump(out[:a.max], open(a.out, "w"))
    print(f"group {a.g} outs {g['outs']}: {len(out[:a.max])} plans at 38, none silent on coordinate {a.t}", flush=True)


if __name__ == "__main__":
    main()
