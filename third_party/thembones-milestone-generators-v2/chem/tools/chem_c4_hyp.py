#!/usr/bin/env python3
"""chem_c4_hyp.py: builder C4. A 6-group at 38 products by an EXACT search: every hyperplane Phi of the group's W-perp
(c = a - 1 = 14 relations; designer R59's 'hyp' test), priced exactly (fan-out kernels + Lovasz matroid parity, builder
C2's Group class, imported read-only); for every hyperplane that prices at 38, explicit plans (C2's enum_plans), scored
by how many products are NOT silent on the later groups' kill coordinates (the engine's kill order).

  python -I chem_c4_hyp.py g.json --g 0 --kill 3,7 --plans 40 --out plans.json
g.json: a group file in C2's e19grp dump format (tools/chem_c4_flips.py --grp). Output: list of {planes: [[y1,y2,s],..],
silent: [...], hyper: lam} sorted by the fewest silent products.
"""
import argparse
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_ld2_e19grp as eg  # noqa: E402


def hyperplanes(G, rel):
    """per fan-out s the space R_s = {F in W-perp: F s = 0} (coefficients over rel); for every hyperplane lam the
    number of fan-outs that keep a plane. Returns list of (n, lam) sorted descending."""
    R = len(rel)
    rs = []
    for s in range(1, 1 << G.b):
        rows = []
        for t in range(G.a):
            v = 0
            for j, r in enumerate(rel):
                w = 0
                for o in eg.bits(s):
                    w ^= (r >> (o * G.a)) & G.full_a
                if (w >> t) & 1:
                    v |= 1 << j
            rows.append(v)
        rs.append(eg.nullspace(rows, R))
    out = []
    for lam in range(1, 1 << R):
        n = 0
        for b in rs:
            if len(b) >= 2 or (len(b) == 1 and bin(lam & b[0]).count("1") % 2 == 0):
                n += 1
        out.append((n, lam))
    out.sort(reverse=True)
    return out


def phis_of(rel, lam):
    piv = lam.bit_length() - 1
    return [rel[j] ^ (rel[piv] if (lam >> j) & 1 else 0) for j in range(len(rel)) if j != piv]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--g", type=int, default=0)
    ap.add_argument("--kill", default="", help="global low coordinates of the later kills (silence is scored on them)")
    ap.add_argument("--top", type=int, default=60, help="hyperplanes priced (by fan-out count)")
    ap.add_argument("--plans", type=int, default=40, help="plan tries per hyperplane")
    ap.add_argument("--target", type=int, default=38)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rseed", type=int, default=1)
    a = ap.parse_args()
    d = json.load(open(a.dump))
    g = d["groups"][a.g]
    G = eg.Group(g, seed=a.rseed + 7)
    rel = eg.nullspace(G.wb, G.N)
    rng = random.Random(a.rseed)
    coords = g["coords"]
    kill = [coords.index(int(t)) for t in a.kill.split(",") if t.strip() and int(t) in coords]
    t0 = time.time()
    hs = hyperplanes(G, rel)
    need = a.target
    print(f"group {a.g}: b {G.b} a {G.a} dim_w {g['dim_w']} W-perp {len(rel)}; hyperplanes with >= {need} fan-outs: "
          f"{sum(1 for n, _ in hs if n >= need)} ({time.time() - t0:.1f} s)", flush=True)
    found = []
    for n, lam in hs[:a.top]:
        if n < need:
            break
        phis = phis_of(rel, lam)
        res = G.price(phis, reps=2)
        if res["prod"] is None or res["prod"] > need:
            continue
        plans, nf, nmulti = eg.enum_plans(G, phis, res["prod"], rng, a.plans)
        for pl in plans:
            sil = [sum(1 for (y1, y2, s) in pl if not (((y1 | y2) >> t) & 1)) for t in kill]
            found.append(dict(planes=[list(x) for x in pl], silent=sil, hyper=lam, fanouts=n, prod=res["prod"]))
        print(f"  hyperplane {lam:#x}: {n} fan-outs, priced {res['prod']}; {len(plans)} explicit plans "
              f"(silent on kills {[min(f['silent'][i] for f in found[-len(plans):]) for i in range(len(kill))] if plans and kill else '-'})",
              flush=True)
    found.sort(key=lambda f: (f["prod"], sum(f["silent"])))
    json.dump(found, open(a.out, "w"))
    print(f"{len(found)} plans at {need} -> {a.out}; best silent {found[0]['silent'] if found else None}", flush=True)


if __name__ == "__main__":
    main()
