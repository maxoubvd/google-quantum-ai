#!/usr/bin/env python3
"""chem_c3_gencount.py: builder C3 (second holder). Generic output count of every loader k >= 1 of an UNPADDED table
(1,219 rows), and for the RELATION-BEARING WIRED loaders (24 or fewer generic outputs: wires + bit-1 relations, WR7's
P51 / P52) the tie-engine job: the placed-kill plan of that count (C2's PLANS_BY_N, read-only import), modes
'reg,tie3' + '-' for the later groups, stop = the plan's engine price - 1. One JSON line on stdout.
  CHEM_BROT=7 CHEM_HB=77 ~/qcc/ws/chem_own/venv/bin/python -I tools/chem_c3_gencount.py TABLE.npz
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_wwidth  # noqa: E402,F401
import chem_ld2_e17kill as kl  # noqa: E402


def main():
    core.J, core.BROT = 5, 7
    words = np.asarray(np.load(sys.argv[1])["angle_words"], np.int64)
    words = words[:1219] if words.shape[0] > 1219 else words
    import chem_ld2_e19grp as eg   # read-only: block_set = the padding R59's tie engine applies with --pad
    gen, jobs = {}, []
    for k in range(1, 15):
        try:
            # 12:50Z: count AFTER the loader's own padding (RK7 run 24's k13: 23 unpadded, 24 padded = the engine's view)
            _, wp = eg.block_set(words, k)
            an = core.analyze(wp, k, 4)
        except Exception as ex:  # noqa: BLE001
            gen[k] = f"error {ex!r}"[:80]
            continue
        n = int(an.X0.shape[1])
        gen[k] = n
        if n == 23:
            # 12:40Z [M on WR7's imposed G52 k14]: 6,6,7:1:0:0,4:1:2:1 = 39 | 38 | 43 | 28 = 148 with ties, 147 with
            # the annealed 7-group (R59's g42 driver: emitted 221 on a table with carriers); C2's PLANS_BY_N[23]
            # (one-output groups) closes at 152 on one seed in 12
            jobs.append(dict(k=k, n=n, plan="6,6,7:1:0:0,4:1:2:1", modes="reg,tie3,-,-", stop=147, driver="g42",
                             wantg="39,38,43,28"))
        elif n <= 24 and n in kl.PLANS_BY_N:
            plan, price = kl.PLANS_BY_N[n]
            ng = len(kl.parse_plan(plan))
            jobs.append(dict(k=k, n=n, plan=plan, modes=",".join(["reg", "tie3"] + ["-"] * (ng - 2)), stop=price - 1,
                             driver="tie", wantg="-"))
    print(json.dumps(dict(gen=gen, jobs=jobs)))


if __name__ == "__main__":
    main()
