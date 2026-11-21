#!/usr/bin/env python3
"""inv_i39_patch_slot.py -- I39 (builder): STEP 2 SCAFFOLD. Patches MY source copy's stage module
(~/qcc/ws/inv_i39/src/invsqrt_fc2_patch_h1c.py, a copy of FC's k9f stage) with the PARITY SLOT, all behind the cfg
option h1c "slot": {"par": P} (absent = the file behaves exactly as before):
  (1) st_main: a fresh slot wire below DIV2's y19: U = [slot] + y19 = 2 y19 (20 wires, the stage's own width);
  (2) the ROUNDING after the passes and the guard's clear, before the carry rule: b = slot xor t (t = T0 xor P, T0 =
      Tf0 xor s while T is folded) on a fresh wire, then R[1:] += b (slot := t, + (u AND NOT t) above), as ONE block
      whose undo is the inverse with c0 replaced by the phantom top (FC8's c0m: c0 is measured before that undo runs;
      its writes become phases under m0, `invsqrt_fc_s8m3.ph_fix`);
      SCAFFOLD FORM: b stays on its wire to the mirror (the cheap forms: b as the last pass's carry-in, b erased by
      X-measurement with its phase repaid in the mirror: 'BUILD STEP 2' of research/ws_inv_i39.md);
  (3) c0 CLEARED after T += C: c0 = T'0 xor slot xor P (two CNOTs and an X), its wire freed; the mirror re-makes it
      (plain ops are mirrored) before the stage's undo blocks read it.
The model the circuit must meet: U = 2 y19 + hp(T), hp = hn + ((hn xor t) & 1) (tools/inv_i39_slotmodel.py).
Idempotent (a marker line); prints what it did.   inv_i39_patch_slot.py [SRC]     Safe: edits only MY src copy."""
import os
import sys

SRC = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/qcc/ws/inv_i39/src")
P = os.path.join(SRC, "invsqrt_fc2_patch_h1c.py")
assert "/inv_i39/" in P.replace("\\", "/"), "only my own source copy"
s = open(P).read()
MARK = "# I39-SLOT (inv_i39_patch_slot.py)"
if MARK in s:
    print("already patched:", P)
    sys.exit(0)

A0 = "def stage(F8, c, U, cfg):\n"
assert s.count(A0) == 1, "anchor 0"
R0 = '''_I39_WPOS = {"X": (1,), "CX": (1, 2), "CCX": (1, 2, 3), "SWAP": (1, 2), "CSWAP": (1, 2, 3), "Z": (1,),
             "CZ": (1, 2), "CCZ": (1, 2, 3), "ZP": (1,), "HMR": (1,), "ANDC": (1, 2, 3), "ANDE": (1, 2, 3),
             "CANDC": (1, 2, 3), "CANDE": (1, 2, 3), "NEG": (), "BBEG": (), "BEND": (), "CBEG": (), "CEND": ()}


def _i39_remap(o, a, b):                               MARK: wire a -> b in one op
    pos = _I39_WPOS.get(o[0])
    if pos is None:
        assert a not in o[1:], ("i39 remap: unknown op touching the wire", o)
        return o
    return tuple(b if (k in pos and v == a) else v for k, v in enumerate(o))


'''.replace("MARK", MARK) + A0
s = s.replace(A0, R0)

A1 = "    c0 = cw[0]\n    if nc == 1:"
assert s.count(A1) == 1, "anchor 1"
R1 = '''    _slot = o.get("slot")                              MARK: the rounding (scaffold: b held)
    if _slot is not None:
        par_ = int(_slot.get("par", 0))
        bw_ = c.alloc(1)[0]
        c.cx(L[0], bw_)                                # u = the slot after the passes
        c.cx(T[0], bw_)                                # Tf0 = T0 xor s (T is folded here)
        c.cx(s, bw_)
        if par_:
            c.x(bw_)                                   # b = u xor t, t = T0 xor par
        _bid = len(c.blocks)                           # one block: forward = the add on the real wires;
        c.blocks[_bid] = None                          # its undo = the inverse with c0 -> the phantom top
        c.ops.append(("BBEG", _bid))                   # (c0m: c0 is measured before this undo runs)
        _i0 = len(c.ops)
        A.add(c, [bw_], list(L) + [cw[0]])             # slot := t, + (u AND NOT t) above (c0 = top)
        _seg = c.ops[_i0:]
        c.ops.append(("BEND", _bid))
        import invsqrt_fc_kit as _KK
        _und = _KK.Circ.mirror_blocks(_seg, c.blocks)
        if _M3c.PH:
            _und = _M3c.ph_fix([_i39_remap(o_, cw[0], _M3c.PH["P"]) for o_ in _und])
        c.blocks[_bid] = _und
        F8.__dict__["_i39_bw"] = bw_
'''.replace("MARK", MARK) + A1
s = s.replace(A1, R1)

A2 = ("    c.free_([c1, c2])\n"
      "    _M3c.PH.clear()                                    # FC8 c0m off after the stage\n"
      "    return c0\n")
assert s.count(A2) == 1, "anchor 2"
R2 = '''    c.free_([c1, c2])
    _M3c.PH.clear()                                    # FC8 c0m off after the stage
    if _slot is not None:                              MARK: c0 = T'0 xor slot xor par
        c.cx(T[0], c0)
        c.cx(L[0], c0)
        if par_:
            c.x(c0)
        c.free_([c0])
        return None
    return c0
'''.replace("MARK", MARK)
s = s.replace(A2, R2)

A3 = ("        c.tag(\"H1\")\n"
      "        c0 = stage(F8, c, U, cfg_)\n"
      "        F8.__dict__[\"_fc2_h1c_held\"] = [c0]\n")
assert s.count(A3) == 1, "anchor 3"
R3 = '''        c.tag("H1")
        if opts(cfg_).get("slot") is not None:          MARK: U = [slot] + y19 = 2 y19
            assert len(U) == 9 + int(opts(cfg_)["nb"]), ("slot: DIV2 must leave 19 wires", len(U))
            U = [c.alloc(1)[0]] + list(U)
            F8.__dict__["_i39_U"] = list(U)
        c0 = stage(F8, c, U, cfg_)
        F8.__dict__["_fc2_h1c_held"] = [c0] if c0 is not None else [F8.__dict__["_i39_bw"]]
'''.replace("MARK", MARK)
s = s.replace(A3, R3)
open(P, "w").write(s)
print("patched:", P)
