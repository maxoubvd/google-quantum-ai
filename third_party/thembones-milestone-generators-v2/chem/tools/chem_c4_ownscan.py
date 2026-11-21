#!/usr/bin/env python3
"""chem_c4_ownscan.py: builder C4. For a design-seat-R65 owner loader (plan 6 | 6:k | 6:k | 7:k, placement by dead
supply): every first group made of 6 wires from R65's group-1 and group-4 wires (groups 2 and 3 kept as R65 placed
them), tested EXACTLY for a pool-free 38 whose products are all non-silent on group 2's kill coordinate t (hyperplanes
of W-perp + matroid parity on the non-silent kernels). W of a first group = the row space of its wires' mixed tables
(the wire basis, chem_c4_wb.py), so no engine run is needed for the test.

  python -I chem_c4_ownscan.py <table.npz> --k 10 --place "0:2,0:3,...(25 w:b in R65 group order)" --t 14 \
      --sizes 6,6,6,7 --part 0/12 --out cand.json
"""
import argparse
import itertools
import json
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_e19grp as eg  # noqa: E402
import chem_c4_hyp as hy  # noqa: E402
import chem_c4_ns as ns  # noqa: E402
import chem_c4_wb as wb  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--place", required=True)
    ap.add_argument("--sizes", default="6,6,6,7")
    ap.add_argument("--t", type=int, default=14)
    ap.add_argument("--part", default="0/1")
    ap.add_argument("--top", type=int, default=6)
    ap.add_argument("--from", dest="src", default="0,3", help="R65 groups the first group is drawn from (0-based)")
    ap.add_argument("--minkeep", type=int, default=0, help="at least this many of R65's group-1 wires")
    ap.add_argument("--cands", default="", help="JSON list of first groups (lists of 'w:b') to test instead")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    words = np.asarray(np.load(a.table)["angle_words"], np.int64)
    V, words = eg.block_set(words, a.k)
    an = core.analyze(words, a.k, 4)
    wires = wb.wire_basis(an)                      # basis output j = wire wires[j]
    pos = {tuple(an.out_wires[q]): j for j, q in enumerate(wires)}
    seq = [tuple(int(x) for x in t.split(":")) for t in a.place.split(",")]
    sizes = [int(x) for x in a.sizes.split(",")]
    cut = np.cumsum([0] + sizes)
    grp = [seq[cut[i]:cut[i + 1]] for i in range(len(sizes))]
    # holder 2: --from "0,1,2" draws the first group from R65's groups 1-3 (the last group, which carries the dead-word
    # supply of its two column kills, stays as R65 placed it); --minkeep: at least that many of R65's group-1 wires
    src = [int(x) for x in a.src.split(",")]
    pool = [w for i in src for w in grp[i]]
    nl, nr = an.nl, an.nr
    d = dict(nr=nr, nl=nl, X0=[[int(v) for v in row] for row in an.X0])
    import chem_c4_subscan as sc
    g1set = set(grp[0])
    combos = [cb for cb in itertools.combinations(range(len(pool)), 6)
              if sum(1 for i in cb if pool[i] in g1set) >= a.minkeep]
    if a.cands:                                   # holder 2: an explicit list of first groups (w:b lists)
        pool = sorted({tuple(int(x) for x in w.split(":")) for c in json.load(open(a.cands)) for w in c})
        idx = {w: i for i, w in enumerate(pool)}
        combos = [tuple(idx[tuple(int(x) for x in w.split(":"))] for w in c) for c in json.load(open(a.cands))]
    pi, pn = (int(x) for x in a.part.split("/"))
    combos = combos[pi::pn]
    rng = random.Random(pi)
    found = []
    for cb in combos:
        g1 = [pool[i] for i in cb]
        outs = sorted(pos[w] for w in g1)
        g = sc.group_of(d, outs)
        G = eg.Group(g, seed=7)
        rel = eg.nullspace(G.wb, G.N)
        if len(rel) != 15:
            continue
        hs = hy.hyperplanes(G, rel)
        if hs[0][0] < 38:
            continue
        for n, lam in hs[:a.top]:
            if n < 38:
                break
            fans = ns.nonsilent_fans(G, hy.phis_of(rel, lam), a.t)
            if len(fans) < 38 or G.lovasz(fans) < 38:
                continue
            found.append(dict(g1=[f"{w[0]}:{w[1]}" for w in g1], fans=len(fans), hyper=lam))
            print(f"FOUND g1 {found[-1]['g1']} ({len(fans)} non-silent fan-outs)", flush=True)
            break
    json.dump(found, open(a.out, "w"))
    print(f"DONE {len(combos)} candidates, {len(found)} found", flush=True)


if __name__ == "__main__":
    main()
