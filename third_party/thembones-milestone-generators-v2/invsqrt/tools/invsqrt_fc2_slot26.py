#!/usr/bin/env python3
"""invsqrt_fc2_slot26.py -- FC2 (own code): LP margin (R37 designer's converged lp_conv, ring model, ALL inputs) of
today's m 8 arithmetic (w41e: FC2's h1b term = invsqrt_fc2_h1m2.stage_word(y, 1, koff 6, ctau [256], drops (1,8,1),
(3,6,1)), K 1533) with MAIN's word cut to 26 slots: rows 26..blast, slot 0 = zero (i.e. a ring of N = 26 tokens:
one slot rotation less) / one / rest / nrest, against the full 27-row word.
  invsqrt_fc2_slot26.py SRC WD OUT.jsonl NPROC"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r37_inv_40_model as RM            # noqa: E402  (reads sys.argv[1] = SRC)
from invsqrt_fc2_h1m2 import stage_word  # noqa: E402

CASES = [dict(name="w27_std", blast=0, slot="std"), dict(name="w26_zero", blast=1, slot="zero"),
         dict(name="w26_rest", blast=1, slot="rest"), dict(name="w26_one", blast=1, slot="one"),
         dict(name="w26_nrest", blast=1, slot="nrest")]
if len(sys.argv) > 5:                                  # variants: cubic thresholds / dropped parts / koff at 26 slots
    CASES = json.loads(sys.argv[5])


def run(case, wd):
    t0 = time.time()
    z = np.load(os.path.join(wd, "chain.npz"))
    tgt, tol = z["target"].astype(np.uint64), z["tol"].astype(np.float64)
    k, d, y, delta = RM.chain(8)
    dps = case.get("drop_parts", [[1, 8, 1], [3, 6, 1]])
    koff = int(case.get("koff", 6))
    if "Kt" in case:                                   # a target constant: koff from the natural one of this drop set
        _, _, _, K0 = stage_word(y, 1, 0, 0, 0, case.get("ctau", [256]), (), dps)
        koff = int(case["Kt"]) - K0
    U, C, kept, K = stage_word(y, 1, koff, 0, 0, case.get("ctau", [256]), (), dps)
    Y = RM.main_rows(U, d, 8, int(case.get("bt", 7)), case["blast"], case["slot"])
    fire, ang, fixed = RM.ring_rows(Y, k)
    if case["slot"] == "zero":
        assert not ((Y & 1) != 0).any()
    tot = fixed.copy()
    for a_, f_ in zip(ang, fire):
        tot += np.where(f_ == 1, np.uint64(a_), np.uint64(0))
    err = (tot - tgt).view(np.int64).astype(np.float64)
    mg, conv, bound = RM.lp_conv(err, tol, fire)
    nrot = 33 if case["slot"] == "zero" else 34
    out = dict(case, K=K, margin=mg, conv=conv, bound=bound, nrot=nrot, tper=round(RM.tper(mg, nrot), 4),
               Tsum=round(nrot * RM.tper(mg, nrot), 2), t=round(time.time() - t0, 1))
    return out


def main():
    wd, outp, nproc = sys.argv[2], sys.argv[3], int(sys.argv[4])
    import multiprocessing as mp
    with mp.Pool(nproc) as pool, open(outp, "a") as f:
        for o in pool.starmap(run, [(c_, wd) for c_ in CASES]):
            f.write(json.dumps(o) + "\n")
            print(json.dumps(o), flush=True)


if __name__ == "__main__":
    main()
