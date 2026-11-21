#!/usr/bin/env python3
"""chem_c4_pre.py: builder C4. Fast pre-checks for a pair build (two pool-free 38-groups injected as groups 1 and 2):
  (1) for every plan of each group: its products' high vectors (the engine's own solve: rows of the group's residual
      in the feature basis), the kill column R_p[:, t] of every product, and the closure parity (C2's tau_state);
  (2) for plan pairs: the rank of the 76 kill columns (group 3's column kill on t needs 75);
  (3) wires completed by each group under the output order (closure bits a group must meet by retries).

  python -I chem_c4_pre.py <table.npz> --k 2 --pair DIR [--pairs 40]
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402

bits = core.bits


def gf2_rank(vecs):
    piv = {}
    r = 0
    for v in vecs:
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                piv[h] = v
                r += 1
                break
    return r


def products_of(an, outs, plan):
    """the engine's solve for one injected group: rows of the residual X0 restricted to outs, coordinates in the
    plan's feature basis -> A vectors; returns products [y1, A1, y2, A2, s_global]."""
    nl, nr = an.nl, an.nr
    A_ = nl
    feats = []
    for (y1, y2, s) in plan:
        for y in (y1, y2):
            feats.append(sum(y << (gi * A_) for gi in bits(s)))
    sp = core.Span()
    for i, f in enumerate(feats):
        assert sp.add(f, 1 << i), "dependent features"
    Acoef = [0] * len(feats)
    for M in range(nr):
        v = 0
        for gi, o in enumerate(outs):
            x = int(an.X0[M, o])
            for ci in range(nl):
                if (x >> ci) & 1:
                    v |= 1 << (gi * A_ + ci)
        r, c = sp.reduce(v)
        assert r == 0, "row outside the plan's span"
        for idx in bits(c):
            Acoef[idx] |= 1 << M
    prods = []
    for k, (y1, y2, s) in enumerate(plan):
        sg = sum(1 << outs[gi] for gi in bits(s))
        prods.append([y1, Acoef[2 * k], y2, Acoef[2 * k + 1], sg, 0])
    return prods


def killcol(p, t):
    y1, A1, y2, A2 = p[:4]
    return (A1 if (y1 >> t) & 1 else 0) ^ (A2 if (y2 >> t) & 1 else 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--pair", required=True)
    ap.add_argument("--pairs", type=int, default=400)
    ap.add_argument("--wirebasis", action="store_true")
    ap.add_argument("--sizes", default="6,6,7,6", help="group sizes (24 outputs: 6,6,6,6; 23: 6,6,7,4)")
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    import chem_ld2_c4 as c4
    import chem_ld2_e19grp as eg
    words0 = np.asarray(np.load(a.table)["angle_words"], np.int64)
    V, words = eg.block_set(words0, a.k)
    an = core.analyze(words, a.k, 4)
    if a.wirebasis:
        import chem_c4_wb
        chem_c4_wb.wire_basis(an)
    kap = c4.kappa_wires(an, V)
    args = open(os.path.join(a.pair, "args.txt")).read().split()
    order = [int(x) for x in args[args.index("--order") + 1].split(",")]
    t = int(args[args.index("--kills") + 1].split(",")[0])
    S1, S2, rest = order[:6], order[6:12], order[12:]
    P1 = json.load(open(os.path.join(a.pair, "g0.json")))
    P2 = json.load(open(os.path.join(a.pair, "g1.json")))
    # wires completed by each group under this order (holder 2: group sizes from --sizes; 25 outputs: 6,6,7,6)
    sz = [int(x) for x in a.sizes.split(",")]
    groups = [S1, S2, rest[:sz[2]], rest[sz[2]:]]
    done = 0
    wm = []
    for g in groups:
        om = sum(1 << o for o in g)
        new = done | om
        w = 0
        for q, gq in enumerate(an.G):
            if gq and (gq & om) and not (gq & ~new):
                w |= 1 << q
        wm.append(w)
        done = new
    print(f"closure bits a group must meet (this order): {[bin(w).count('1') for w in wm]}", flush=True)
    cs = [0] + list(np.cumsum(sz)); nat = [list(range(cs[i], cs[i + 1])) for i in range(4)]
    done, wn = 0, []
    for g in nat:
        om = sum(1 << o for o in g)
        new = done | om
        wn.append(sum(1 << q for q, gq in enumerate(an.G) if gq and (gq & om) and not (gq & ~new)))
        done = new
    print(f"closure bits in the engine's natural order: {[bin(w).count('1') for w in wn]}", flush=True)
    # per plan: products, kill columns, closure parity on the wires completed by that group
    info = []
    for gi, (S, PL) in enumerate(((S1, P1), (S2, P2))):
        rows = []
        prev = [] if gi == 0 else None
        for pl in PL[:a.pairs]:
            pr = products_of(an, S, [tuple(x) for x in pl])
            rows.append(pr)
        info.append(rows)
    # closure: group 1's parity on wm[0], group 2's on wm[1] (with group 1's products present)
    ok1 = [i for i, pr in enumerate(info[0])
           if not ((c4.tau_state(an, c4.normalize([list(p) for p in pr]), V) ^ kap) & wm[0])]
    best = None
    tested = 0
    for i in ok1:
        for j, pr2 in enumerate(info[1]):
            allp = info[0][i] + pr2
            if (c4.tau_state(an, c4.normalize([list(p) for p in allp]), V) ^ kap) & wm[1]:
                continue
            rk = gf2_rank([killcol(p, t) for p in allp])
            tested += 1
            if best is None or rk > best[0]:
                best = (rk, i, j)
            if rk >= 75:
                break
        if best and best[0] >= 75:
            break
    print(f"plans: group 1 {len(info[0])} ({len(ok1)} pass closure), group 2 {len(info[1])}; pairs passing both "
          f"closures tested {tested}; best kill-column rank on t={t}: {best[0] if best else None} (75 needed) at plans "
          f"{best[1:] if best else None}", flush=True)
    if best and best[0] >= 75:
        json.dump([P1[best[1]]], open(os.path.join(a.pair, "g0_best.json"), "w"))
        json.dump([P2[best[2]]], open(os.path.join(a.pair, "g1_best.json"), "w"))
        print("WROTE g0_best.json / g1_best.json", flush=True)


if __name__ == "__main__":
    main()
