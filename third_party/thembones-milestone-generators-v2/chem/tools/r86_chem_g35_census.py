#!/usr/bin/env python3
"""r86_chem_g35_census.py: seat R86-CHEM-G35. For a dumped group 4 (r86_chem_g35_dump.py): the exact delivery deficit
of C (x) psi (three forced relations sharing one low vector psi) over EVERY payer with enough dead unknowns for the
rank-70 step and EVERY 3-dim C spanned by singles / pairs of the five other outputs, for many psi. Deficit = rank of the
augmented system minus the coefficient rank (0 = deliverable); also the coefficient rank against 3 (nr - 1).

  python -I r86_chem_g35_census.py k9_g4.npz --psis 30 --seed 1 [--payer-min 0]
"""
import argparse
import itertools
import json
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def par(x):
    return bin(x).count("1") & 1


def rank_aug(cols, rhs):
    piv = {}
    r = 0
    for v in cols:
        while v:
            h = v.bit_length() - 1
            e = piv.get(h)
            if e is None:
                piv[h] = v
                r += 1
                break
            v ^= e
    v = rhs
    while v:
        h = v.bit_length() - 1
        e = piv.get(h)
        if e is None:
            return r, r + 1
        v ^= e
    return r, r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--psis", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--payer-min", type=int, default=-1, help="least payer dead unknowns (default 70 - K)")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    z = np.load(a.dump)
    meta = json.load(open(a.dump.replace(".npz", ".json")))
    outs, live, nr, K = meta["outs"], meta["live"], meta["nr"], meta["K"]
    gi = {o: i for i, o in enumerate(outs)}
    R0 = [[int(v) for v in row] for row in z["R0"]]
    N = [[int(v) for v in row] for row in z["N"]]
    togg = {o: [c["Ec"] for c in meta["comp"] if c["o"] == o] for o in outs}
    pmin = a.payer_min if a.payer_min >= 0 else 70 - K
    payers = [o for o in outs if len(togg[o]) >= pmin]
    best = None
    tried = 0
    found = []
    for it in range(a.psis):
        psi = rng.getrandbits(15) & live
        if not psi:
            continue
        nv = []
        for j in range(K):
            v = 0
            for r in range(nr):
                if par(N[j][r] & psi):
                    v |= 1 << r
            nv.append(v)
        tv, rv = {}, {}
        for o in outs:
            lst = []
            for Ec in togg[o]:
                v = 0
                for r in range(nr):
                    if par(Ec[r] & psi):
                        v |= 1 << r
                lst.append(v)
            tv[o] = lst
            v = 0
            for r in range(nr):
                if par(R0[gi[o]][r] & psi):
                    v |= 1 << r
            rv[o] = v
        for p in payers:
            nonp = [o for o in outs if o != p]
            vecs = [(o,) for o in nonp] + list(itertools.combinations(nonp, 2))
            for trip in itertools.combinations(vecs, 3):
                # independence of the three combinations over GF(2)^5
                ms = [sum(1 << nonp.index(o) for o in c) for c in trip]
                if ms[0] ^ ms[1] == 0 or ms[0] ^ ms[2] == 0 or ms[1] ^ ms[2] == 0 or ms[0] ^ ms[1] ^ ms[2] == 0:
                    continue
                cols = []
                for o in nonp:
                    where = [i for i, c in enumerate(trip) if o in c]
                    if not where:
                        continue
                    for v in nv + tv[o]:
                        col = 0
                        for i in where:
                            col |= v << (i * nr)
                        cols.append(col)
                rhs = 0
                for i, c in enumerate(trip):
                    acc = 0
                    for o in c:
                        acc ^= rv[o]
                    rhs |= acc << (i * nr)
                rc, ra = rank_aug(cols, rhs)
                tried += 1
                key = (ra - rc, -rc)
                if best is None or key < best[0]:
                    best = (key, dict(psi=hex(psi), payer=p, rels=[list(c) for c in trip], coef_rank=rc, deficit=ra - rc))
                if ra == rc:
                    found.append(dict(psi=hex(psi), payer=p, rels=[list(c) for c in trip], coef_rank=rc))
    print(json.dumps(dict(dump=os.path.basename(a.dump), K=K, nr=nr,
                          toggles={str(o): len(togg[o]) for o in outs}, payers=payers, layouts_tried=tried,
                          best=best[1] if best else None, consistent_layouts=len(found), first_found=found[:5])),
          flush=True)


if __name__ == "__main__":
    main()
