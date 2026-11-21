#!/usr/bin/env python3
"""inv_i39_patch_slot3.py -- I39 (builder): two fixes on top of tools/inv_i39_patch_slot2.py (my src copy only).
  (1) repair (a)'s lending: the table word subtracts, so its target wires are X-framed BEFORE the table is made; the
      slot (zero before the frame) then holds 1: it is X-ed back to 0 while lent and re-framed after (c0 is never
      framed: c0m's ph_real).
  (2) REPAIR (b) (R49-INV-39W: "the undos at 8 clean wires write c0's phantom top with gated CCX and conditioned CX:
      ph_fix asserts"): with the slot on, invsqrt_fc_s8m3.ph_fix also turns GATED writes into the phantom into phases
      under a derived classical bit d = m0 AND g (FC's kit extension "BAND", the official verifier's nested condition):
      XCCX a b P g -> BAND d (m0, g); CZ a b if d.   CX x P if g -> Z x if d.   X P if g -> NEG if d.
Idempotent (marker).   inv_i39_patch_slot3.py [SRC]     Safe: edits only MY src copy."""
import os
import sys

SRC = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/qcc/ws/inv_i39/src")
P = os.path.join(SRC, "invsqrt_fc2_patch_h1c.py")
assert "/inv_i39/" in P.replace("\\", "/"), "only my own source copy"
s = open(P).read()
assert "# I39-SLOT2 (inv_i39_patch_slot2.py)" in s, "run tools/inv_i39_patch_slot2.py first"
MARK = "# I39-SLOT3 (inv_i39_patch_slot3.py)"
if MARK in s:
    print("already patched:", P)
    sys.exit(0)

A1 = ("    if _lend:\n"
      "        c.free_(_lend)\n"
      "    k_f = table(bw)\n"
      "    if _lend:\n"
      "        _take(c, _lend)\n")
assert s.count(A1) == 1, "anchor 1"
R1 = ("    _lx = [w_ for w_ in _lend if w_ in _M3x.ph_real(Rt)]   MARK: framed (X) lent wires hold 1: X them to 0\n"
      "    for w_ in _lx:\n"
      "        c.x(w_)\n"
      "    if _lend:\n"
      "        c.free_(_lend)\n"
      "    k_f = table(bw)\n"
      "    if _lend:\n"
      "        _take(c, _lend)\n"
      "    for w_ in _lx:\n"
      "        c.x(w_)\n").replace("MARK", MARK)
s = s.replace(A1, R1)

A0 = "def stage(F8, c, U, cfg):\n"
assert s.count(A0) == 1, "anchor 0"
R0 = '''_I39_C = {"c": None}


def _i39_ph_fix(ops):                                  MARK: ph_fix + gated writes into the phantom
    import invsqrt_fc_s8m3 as M3_
    if not M3_.PH:
        return ops
    Pw, m0 = M3_.PH["P"], M3_.PH["m0"]
    c_ = _I39_C["c"]
    out, rest = [], []

    def dbit(g):
        d = c_.bit()
        if not hasattr(c_, "bprob"):
            c_.bprob = {}
        c_.bprob[d] = 0.5 * c_.bprob.get(g, 0.5)
        out.append(("BAND", d, ((m0, 1), (g, 1))))
        return d
    for o in ops:
        g = o[0]
        if g == "XCCX" and o[3] == Pw:
            assert Pw not in o[1:3]
            out.append(("CZ", o[1], o[2], dbit(o[4])))
        elif g == "CX" and o[2] == Pw and o[3] is not None:
            out.append(("Z", o[1], dbit(o[3])))
        elif g == "X" and o[1] == Pw and o[2] is not None:
            out.append(("NEG", dbit(o[2])))
        else:
            out.extend(M3_._i39_orig_ph_fix([o]))
    return out


'''.replace("MARK", MARK) + A0
s = s.replace(A0, R0)

A2 = '''    _slot = o.get("slot")                              # I39-SLOT (inv_i39_patch_slot.py): the rounding (scaffold: b held)\n'''
assert s.count(A2) == 1, "anchor 2"
s = s.replace(A2, A2)                                   # (unchanged; the install goes at the stage's start below)

A3 = "    L, T = list(U[:9]), list(U[9:])\n    s = T[nb]\n"
assert s.count(A3) == 1, "anchor 3"
R3 = A3 + ('''    if o.get("slot") is not None:                      MARK: gated phantom writes -> phases
        import invsqrt_fc_s8m3 as _M3i
        _I39_C["c"] = c
        if not hasattr(_M3i, "_i39_orig_ph_fix"):
            _M3i._i39_orig_ph_fix = _M3i.ph_fix
        _M3i.ph_fix = _i39_ph_fix
''').replace("MARK", MARK)
s = s.replace(A3, R3)
open(P, "w").write(s)
print("patched:", P)
