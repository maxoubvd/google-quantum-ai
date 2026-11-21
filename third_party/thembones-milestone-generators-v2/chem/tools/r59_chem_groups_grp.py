#!/usr/bin/env python3
"""r59_chem_groups_grp.py: designer R59-CHEM-GROUPS. Builder C2's group tool (tools/chem_ld2_e19grp.py: dump, anneal,
plan, build, lk, stats) run UNCHANGED over the tie engine: the name chem_ld2_e17kill.KillDecomp is rebound, in this
process only, to TieKill with the group modes of the environment variable R59_MODES (default "reg,tie3,-,-"), so the
builder's Probe / PlanKill subclasses (dump and build with an injected group plan) sit on top of groups 1-2 planned
with column ties. Same command line as chem_ld2_e19grp.py; no file of the builder is edited.

    R59_MODES=reg,tie3,-,- CHEM_BROT=7 CHEM_HB=77 python -I r59_chem_groups_grp.py dump <sub.npz> --k 4 \
        --plan "6,6,7:1:0:0,6:1:3:1" --avoid 0,1 --pad --seed 1 --out g.json
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_ld2_e17kill as kl  # noqa: E402
import r59_chem_groups_tie as rt  # noqa: E402  (TieKill is built on the ORIGINAL KillDecomp, bound at import)

MODES = os.environ.get("R59_MODES", "reg,tie3,-,-").split(",")
GTRIALS = int(os.environ.get("R59_GTRIALS", "40"))


class TK(rt.TieKill):
    def __init__(self, an, seed, effort, plan, avoid_pos, V=None, kap=0, tries=1):
        super().__init__(an, seed, effort, plan, avoid_pos, MODES, V, kap, tries, GTRIALS)


kl.KillDecomp = TK
import chem_ld2_e19grp as eg  # noqa: E402

if __name__ == "__main__":
    eg.main()
