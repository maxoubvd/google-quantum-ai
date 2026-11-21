#!/usr/bin/env python3
"""r89_chem_cascade_strict.py (seat R89-CHEM-CASCADE): run another tool (C4's chem_c4_build.py, R86's r86_chem_g35_dump.py
or r86_chem_g35_hook.py, all read-only) with builder C4's dead-paid ROW DROP made STRICT in one group.

Why: chem_c4_dead.solve_dead_row returns None when its draws miss (measured on k10 tonight: 2 of 78 and 12 of 152 calls
solved), and the engine then keeps the group WITHOUT the drop; the closure retry accepts the first attempt whose tau
holds, so the drop usually ends unpaid (group 3 at 38, the loader at 223). R86's hook answers group 4 from a dump of the
engine state at the first group-4 call, so the dump run and the hook run must both reach group 4 with the drop PAID.
Here a miss in group R89_ROWGI (default 2 = group 3) raises RuntimeError: the engine's closure loop retries the group
(deterministic per seed, so the dump run and the hook run reach the same state). R89_ROWMAX misses in a row (default
600) fall back to C4's behaviour.
  R89_ROWGI=2 python -I r89_chem_cascade_strict.py r86_chem_g35_dump.py <its arguments>
"""
import os
import runpy
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_c4_dead as cd  # noqa: E402

GI = int(os.environ.get("R89_ROWGI", "2"))
MAXM = int(os.environ.get("R89_ROWMAX", "600"))
ST = {"calls": 0, "miss": 0, "consec": 0, "paid": 0}
_orig = cd.solve_dead_row


def strict_row(dec, outs, live, killed, R0, N, K, P, usable, sols, kern, nr, b, rows_of, base):
    res = _orig(dec, outs, live, killed, R0, N, K, P, usable, sols, kern, nr, b, rows_of, base)
    if getattr(dec, "_gi", -1) != GI:
        return res
    ST["calls"] += 1
    if res is not None:
        ST["paid"] += 1
        ST["consec"] = 0
        return res
    ST["miss"] += 1
    ST["consec"] += 1
    if ST["consec"] < MAXM:
        raise RuntimeError("r89: row drop unpaid on this attempt")
    return None


cd.solve_dead_row = strict_row                 # chem_c4_dead.install_row binds this name into chem_c4_drop.DEAD_ROW_HOOK

if __name__ == "__main__":
    target = sys.argv[1]
    sys.argv = [target] + sys.argv[2:]
    try:
        runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), target), run_name="__main__")
    finally:
        print(f"R89 strict row drop (group index {GI}): {ST}", flush=True)
