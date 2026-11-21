"""R58-CHEM-KEPT: STAGE A OF THE POINT-STAGE LOADER PLANNED BY BUILDER C2's ENGINE (imported read-only).

C2's basis of the loader (15 low vectors in five blocks) with ONE block replaced by the indicators of three low
points D; C2's engine (search_protected, E.plan) plans the other four blocks on the table in that basis: a plan that
is exact on the 13 low points outside D. Output: the products in the OLD low coordinates, for
tools/r58_chem_kept_mask.py --planA.

  CHEM_BROT=7 CHEM_HB=77 python -I r58_chem_kept_c2a.py PAD.npz --k 1 --pfile k1b4.json --rep 4 --D 13,14,15
        [--seed 0] [--effort 1] --out planA.json
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_wwidth  # noqa: E402,F401
import chem_ld2_basis as lb2  # noqa: E402
import chem_ld2_engine as e2  # noqa: E402
import r55_chem_engine2_rp_c4stair as rs  # noqa: E402

bits = core.bits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pad")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--pfile", required=True)
    ap.add_argument("--rep", type=int, required=True, help="the block of the basis replaced by the point block")
    ap.add_argument("--D", default="13,14,15")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--effort", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t0 = time.time()
    core.J = 5
    words = np.asarray(np.load(a.pad)["angle_words"], np.int64)
    an = core.analyze(words, a.k, 4)
    nl = an.nl
    D = [int(x) for x in a.D.split(",")]
    LB = []
    for bs in an.basis:
        tt = 0
        for l in range(16):
            v = 0
            for m in bits(bs["anf"]):
                if (l & m) == m:
                    v ^= 1
            tt |= v << l
        LB.append(tt)
    sp = core.Span()
    for t in range(nl):
        sp.add(LB[t], 1 << t)
    mD = []
    for p in D:
        r, c = sp.reduce(1 << p)
        assert r == 0
        mD.append(c)
    p0 = json.load(open(a.pfile))
    P = [int(x, 16) for x in p0["P"]]
    five = [list(b) for b in p0["blocks"]]
    P2 = list(P)
    for i, t in enumerate(five[a.rep]):
        P2[t] = mD[i]
    try:
        lb2.mat_inv(P2, nl)
    except ValueError:
        print(f"k {a.k} rep {a.rep} D {D}: the point block is not a complement of the other four blocks")
        sys.exit(3)
    e2.install()
    an2 = lb2.transform(an, P2, free_top=False, blocks=five)
    stA = [five[b] for b in range(5) if b != a.rep]
    m = np.uint32(sum(1 << t for b in stA for t in b))
    pc, price, cod = rs.plan_stage(an2, an2.X0 & m, stA, a.seed, a.effort, 0, 0)
    prods = [[int(x) for x in q] for q in lb2.map_back(pc, P2)]
    json.dump(dict(k=a.k, D=D, rep=a.rep, pfile=os.path.basename(a.pfile), seed=a.seed, effort=a.effort,
                   products=len(prods), codim=[int(x) for x in cod], prods=[[str(x) for x in q] for q in prods]),
              open(a.out, "w"))
    print(f"k {a.k} pfile {os.path.basename(a.pfile)} rep {a.rep} D {D} seed {a.seed}: stage A (C2's engine, four "
          f"blocks) {len(prods)} products, codim {list(cod)}; {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
