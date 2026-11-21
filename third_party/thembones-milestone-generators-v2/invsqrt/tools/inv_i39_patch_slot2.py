#!/usr/bin/env python3
"""inv_i39_patch_slot2.py -- I39 (builder): STEP 2 REPAIR (a) for the 39-wire slot stage (R49-INV-39W's frame repair
"the word's forward table has no schedule within 1 temporaries"): option h1c "slot": {"wfirst": 1} runs the table word
(sw "first") BEFORE the mux add and the cubic, while the slot wire and c0 are still zero, and LENDS those two wires to
the word table's forward emission (two more temporaries: 30 live + 6 table wires + 3 temporaries = 39). The adds commute
(all into R); the word's carry hold (sw "hold") must be off (it would be held across the mux add and the cubic).
Needs tools/inv_i39_patch_slot.py first. Idempotent (marker line).   inv_i39_patch_slot2.py [SRC]   Safe: my src only."""
import os
import sys

SRC = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/qcc/ws/inv_i39/src")
P = os.path.join(SRC, "invsqrt_fc2_patch_h1c.py")
assert "/inv_i39/" in P.replace("\\", "/"), "only my own source copy"
s = open(P).read()
assert "# I39-SLOT (inv_i39_patch_slot.py)" in s, "run tools/inv_i39_patch_slot.py first"
MARK = "# I39-SLOT2 (inv_i39_patch_slot2.py)"
if MARK in s:
    print("already patched:", P)
    sys.exit(0)

A1 = "    k_f = table(bw)\n"
assert s.count(A1) == 1, "anchor 1"
R1 = ('    _lend = list(F8.__dict__.get("_i39_lend") or [])   MARK: slot and c0 lent while they are zero\n'
      "    if _lend:\n"
      "        c.free_(_lend)\n"
      "    k_f = table(bw)\n"
      "    if _lend:\n"
      "        _take(c, _lend)\n").replace("MARK", MARK)
s = s.replace(A1, R1)

A2 = "    if mux:                                            # lin 1 + lin 2 - part 2 of the top row as one mux add\n"
assert s.count(A2) == 1, "anchor 2"
R2 = ('    _sw1_pre = None                                    MARK: the table word first\n'
      '    if (o.get("slot") or {}).get("wfirst") and (o.get("sw") or {}).get("first"):\n'
      '        assert not (o.get("sw") or {}).get("hold") or o["slot"].get("wfirst_hold"), "wfirst: hold off"\n'
      '        F8.__dict__["_i39_lend"] = [R[1], cw[0]]\n'
      '        _sw1_pre = _sword_g0(F8, c, o, T, R, s, nb, lits, tog, tstate, par, first=True)\n'
      '        F8.__dict__["_i39_lend"] = []\n' + A2).replace("MARK", MARK)
s = s.replace(A2, R2)

A3 = "        _sw1 = _sword_g0(F8, c, o, T, R, s, nb, lits, tog, tstate, par, first=True)\n"
assert s.count(A3) == 1, "anchor 3"
R3 = ("        _sw1 = (_sword_g0(F8, c, o, T, R, s, nb, lits, tog, tstate, par, first=True) if _sw1_pre is None\n"
      "                else _sw1_pre)                         MARK\n").replace("MARK", MARK)
s = s.replace(A3, R3)
open(P, "w").write(s)
print("patched:", P)
