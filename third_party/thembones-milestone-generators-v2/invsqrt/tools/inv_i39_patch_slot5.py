#!/usr/bin/env python3
"""inv_i39_patch_slot5.py -- I39 (builder): the table word's forward add at ONE wire less (the last 40-window of repair
(a)): with the slot's "wfirst" on, every table bit is X-measured and freed right after the add's step that reads it
(the add streams its sources: FC2's _fwd_stream calls a release after each step), so the carries of the next steps
reuse those wires (38 live + 1 carry = 39). The phases of those measurements are repaid exactly as before (the undo
recomputes the table and applies Z on each bit under its outcome); the undo takes back the same wires. Needs
tools/inv_i39_patch_slot4.py first. Idempotent (marker).   inv_i39_patch_slot5.py [SRC]   Safe: my src copy only."""
import os
import sys

SRC = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/qcc/ws/inv_i39/src")
P = os.path.join(SRC, "invsqrt_fc2_patch_h1c.py")
assert "/inv_i39/" in P.replace("\\", "/"), "only my own source copy"
s = open(P).read()
assert "# I39-SLOT4 (inv_i39_patch_slot4.py)" in s, "run tools/inv_i39_patch_slot4.py first"
MARK = "# I39-SLOT5 (inv_i39_patch_slot5.py)"
if MARK in s:
    print("already patched:", P)
    sys.exit(0)

A1 = ("            k_f, c.live, c.cap, names, None if car is None else T.index(car), bool(first)))\n"
      "    lits_ = [((lambda w_=w_: w_), (lambda: None)) for w_ in srcw]\n")
assert s.count(A1) == 1, "anchor 1"
R1 = ("            k_f, c.live, c.cap, names, None if car is None else T.index(car), bool(first)))\n"
      "    _strm = bool(F8.__dict__.get(\"_i39_stream\")) and not sw.get(\"hold\")   MARK: table bits measured as read\n"
      "    _mb = {}\n"
      "\n"
      "    def _rel(w_):\n"
      "        if _strm and w_ in bw and w_ not in _mb:\n"
      "            _mb[w_] = c.hmr(w_)\n"
      "            c.free_([w_])\n"
      "    lits_ = [((lambda w_=w_: w_), (lambda w_=w_: _rel(w_))) for w_ in srcw]\n").replace("MARK", MARK)
s = s.replace(A1, R1)

A2 = ("        debts = _fwd_stream(c, NW, lits_, Rt, NW._lit(_cw))\n"
      "    _cin(False)\n"
      "    mbits = [c.hmr(w) for w in bw]\n"
      "    c.free_(bw)\n"
      "    c.ops.append((\"BEND\", bid))\n")
assert s.count(A2) == 1, "anchor 2"
R2 = ("        debts = _fwd_stream(c, NW, lits_, Rt, NW._lit(_cw))\n"
      "    _cin(False)\n"
      "    mbits = [_mb[w] if w in _mb else c.hmr(w) for w in bw]   MARK\n"
      "    c.free_([w for w in bw if w not in _mb])\n"
      "    c.ops.append((\"BEND\", bid))\n").replace("MARK", MARK)
s = s.replace(A2, R2)

A3 = ("    Rt = _M3x.ph_target(Rt)\n"
      "    bw2 = c.alloc(nbit)\n"
      "    assert list(bw2) == list(bw), (\"table wires\", bw2, bw)\n"
      "    k_u = table(bw, _ust8)\n")
assert s.count(A3) == 1, "anchor 3"
R3 = ("    Rt = _M3x.ph_target(Rt)\n"
      "    if _mb:                                            MARK: the same wires back\n"
      "        _take(c, list(bw))\n"
      "    else:\n"
      "        bw2 = c.alloc(nbit)\n"
      "        assert list(bw2) == list(bw), (\"table wires\", bw2, bw)\n"
      "    k_u = table(bw, _ust8)\n").replace("MARK", MARK)
s = s.replace(A3, R3)

A4 = ('        F8.__dict__["_i39_lend"] = [R[1], cw[0]]\n')
assert s.count(A4) == 1, "anchor 4"
R4 = ('        F8.__dict__["_i39_lend"] = [R[1], cw[0]]\n'
      '        F8.__dict__["_i39_stream"] = bool((o.get("slot") or {}).get("wstream"))   MARK\n').replace("MARK", MARK)
s = s.replace(A4, R4)
A5 = ('        F8.__dict__["_i39_lend"] = []\n')
assert s.count(A5) == 1, "anchor 5"
R5 = ('        F8.__dict__["_i39_lend"] = []\n'
      '        F8.__dict__["_i39_stream"] = False\n')
s = s.replace(A5, R5)
open(P, "w").write(s)
print("patched:", P)
