#!/usr/bin/env python3
"""r86_chem_g35_hook.py: seat R86-CHEM-G35. Builder C4's chem_c4_build.main() UNCHANGED (imported, never edited), with
group R86_GI's kill step answered from a solution file of r86_chem_g35_solve.py (pool combinations per output and the
dead toggles that make C (x) psi hold and rank 70), through C4's own hook interface (chem_c4_drop.DEAD_ROW_HOOK returns
(chosen, picked)); the group's 35-product plans go in with C4's own --inj (PlanTie._select). The engine then solves the
high vectors, retries plans until C2's C-4 closure holds (tau == kappa), and C2's emit_c4 checks every wire on every
address. The same environment as the dump run (seed, C4_DEAD_ROW, C4_DEAD_BATCHES) reproduces the engine state.

  R86_SOL=k12_plan35.json R86_DUMPMETA=k12_g4b.json C4_DEAD_BATCHES=both CHEM_BROT=7 CHEM_HB=77 C4_DROP=1 \
  C4_DEAD_ROW=3:2:72 C4_DEAD_OUT=EMIT python -I r86_chem_g35_hook.py table.npz --k 12 --wirebasis \
      --plan 6,6,6:1:0:0,6:1:0:1 --inj 3=k12_plan35_plans.json --emit DIR
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_c4_dead as deadm  # noqa: E402
import chem_c4_drop as dr  # noqa: E402

GI = int(os.environ.get("R86_GI", "3"))
SOL = json.load(open(os.environ["R86_SOL"]))
META = json.load(open(os.environ["R86_DUMPMETA"]))
CALLS = {"n": 0}
_orig = deadm.install_row


def answer(outs, live, killed, R0, N, K, P, usable, sols, kern):
    assert list(outs) == META["outs"], (outs, META["outs"])
    assert [str(c) for c in sols] == META["sols"], "engine state differs from the dump"
    chosen = [int(SOL["comb"][str(o)]) for o in outs]
    comp = {}
    for c in META["comp"]:
        comp.setdefault(c["o"], []).append(c)
    picked = []
    for o in outs:
        for i in SOL["toggle_idx"][str(o)]:
            c = comp[o][i]
            cells = [(cc[0], cc[1], cc[2], (cc[3][0], cc[3][1]), cc[4], cc[5]) for cc in c["cells"]]
            picked.append((o, np.array(c["E"], np.uint32), cells))
    return chosen, picked


def install_row(cells, m, target_gi, conj, rk=None, specs=None):
    _orig(cells, m, target_gi, conj, rk, specs)
    hook = dr.DEAD_ROW_HOOK

    def wrapped(dec, outs, live, killed, R0, N, K, P, usable, sols, kern, nr, b, rows_of, base):
        if getattr(dec, "_gi", -1) == GI:
            CALLS["n"] += 1
            return answer(outs, live, killed, R0, N, K, P, usable, sols, kern)
        return hook(dec, outs, live, killed, R0, N, K, P, usable, sols, kern, nr, b, rows_of, base)
    dr.DEAD_ROW_HOOK = wrapped


deadm.install_row = install_row

if __name__ == "__main__":
    import chem_c4_build
    try:
        chem_c4_build.main()
    finally:
        print(f"R86 hook calls for group {GI}: {CALLS['n']}", flush=True)
