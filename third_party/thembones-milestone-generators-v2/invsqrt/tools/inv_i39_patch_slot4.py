#!/usr/bin/env python3
"""inv_i39_patch_slot4.py -- I39 (builder): THE SLOT WITHOUT A HELD BIT ("exact clear" rounding), on top of
tools/inv_i39_patch_slot3.py, my src copy only. Option h1c "slot": {"mode": "clear", "uprog": PATH}: after the passes
the slot holds u(T) = hn(T) mod 2, a function of T alone (degree 4 on the T that occur [M tools/inv_i39_udeg.py]);
the stage computes u(T) from the folded T wires and the sign by the program PATH (FC's famg format, inputs Tf0..9, s;
`tools/inv_i39_uspec.py` makes its spec) and XORs it INTO the slot (slot := 0), then writes t (slot := T0 xor par).
The register becomes 2 y19 + hq, hq = hn - u + t: a table with the slot's parity and NO carry (bits >= 1 untouched,
so c0 and the carry rule are those of hn), no b wire, plain ops only (the circuit's mirror inverts them).
Idempotent (marker).   inv_i39_patch_slot4.py [SRC]     Safe: edits only MY src copy."""
import os
import sys

SRC = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/qcc/ws/inv_i39/src")
P = os.path.join(SRC, "invsqrt_fc2_patch_h1c.py")
assert "/inv_i39/" in P.replace("\\", "/"), "only my own source copy"
s = open(P).read()
assert "# I39-SLOT3 (inv_i39_patch_slot3.py)" in s, "run tools/inv_i39_patch_slot3.py first"
MARK = "# I39-SLOT4 (inv_i39_patch_slot4.py)"
if MARK in s:
    print("already patched:", P)
    sys.exit(0)

A1 = '''    if _slot is not None:
        par_ = int(_slot.get("par", 0))
        bw_ = c.alloc(1)[0]
'''
assert s.count(A1) == 1, "anchor 1"
R1 = '''    if _slot is not None and _slot.get("mode") == "clear":   MARK: slot := 0 by u(T), then t
        par_ = int(_slot.get("par", 0))
        import invsqrt_fw_chiprog as _CPu
        _up = json.load(open(os.path.expanduser(_slot["uprog"])))
        _up = _up.get("acc", _up.get("prog", _up))
        _upg = {"gates": _up["gates"], "outputs": _up["outputs"]}
        _ins = [T[i_] for i_ in range(nb)] + [s]          # Tf0..9 (T is folded here) and the sign
        _bud = max(1, c.cap - c.live)
        _sch = _CPu.schedule(_upg, len(_ins), _bud)
        assert _sch is not None, ("u(T) program: no schedule within %d temporaries" % _bud)
        F8.__dict__["_i39_uands"] = _CPu.emit(c, _upg, _ins, [L[0]], _bud, _sch)   # slot ^= u(T) -> 0
        c.cx(T[0], L[0])                               # slot := t = T0 xor par (T[0] = Tf0 = T0 xor s)
        c.cx(s, L[0])
        if par_:
            c.x(L[0])
        F8.__dict__["_i39_bw"] = None
    elif _slot is not None:
        par_ = int(_slot.get("par", 0))
        bw_ = c.alloc(1)[0]
'''.replace("MARK", MARK)
s = s.replace(A1, R1)

A2 = '''        F8.__dict__["_fc2_h1c_held"] = [c0] if c0 is not None else [F8.__dict__["_i39_bw"]]
'''
assert s.count(A2) == 1, "anchor 2"
R2 = '''        F8.__dict__["_fc2_h1c_held"] = ([c0] if c0 is not None else
                                        ([F8.__dict__["_i39_bw"]] if F8.__dict__.get("_i39_bw") is not None else []))
'''
s = s.replace(A2, R2)
open(P, "w").write(s)
print("patched:", P)
