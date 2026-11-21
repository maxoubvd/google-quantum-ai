#!/usr/bin/env python3
"""r89_chem_cascade_build.py (seat R89-CHEM-CASCADE): C4's tools/chem_c4_build.py, UNCHANGED, with the dead-cell hook
replaced by r89_chem_cascade_dead.solve_struct (a full second column kill + R89_G random relations in the cells-bound
group 4). Same arguments and environment as chem_c4_build.py (set C4_DEAD to any value to enable the hook, C4_CUBIC=3 for
C4's block). Prints the hook's stats after C4's lines.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_c4_dead  # noqa: E402
import r89_chem_cascade_dead as r89  # noqa: E402

chem_c4_dead.solve_dead = r89.solve_struct          # chem_c4_dead.install() binds this name into chem_c4_drop.DEAD_HOOK
import chem_c4_build  # noqa: E402

if __name__ == "__main__":
    try:
        chem_c4_build.main()
    finally:
        print(f"R89 hook stats {r89.STATS}", flush=True)
