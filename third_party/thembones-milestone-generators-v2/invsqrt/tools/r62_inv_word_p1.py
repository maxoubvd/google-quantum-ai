#!/usr/bin/env python3
"""r62_inv_word_p1.py -- R62-INV-WORD (designer): DEVICE 1 ON THE CUBIC (window W1 of the 38-wire stage).
The cubic table's five stored bits are X-measured AS THE FORWARD ADD READS THEM (option h1c "r62": {"cubstream": 1}):
each bit's wire is released right after its step of the streamed ripple, so the add holds (5 - k) table bits + 2
carries at step k instead of 5 + 2 throughout: 6 wires at its widest instead of 7. The phase of each measurement is
repaid exactly as today (Z on the re-emitted table bit in the block's explicit undo): no gate is added.
Same form as I39's "wstream" for the word (inv_i39_patch_slot5.py). Needs sw "sconj" (nothing touches the table
wires after the ripple).  Idempotent (marker).   r62_inv_word_p1.py [SRC]      Safe to run: edits only MY snapshot."""
import os
import sys

SRC = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/qcc/ws/r62_inv_word/i39/src")
P = os.path.join(SRC, "invsqrt_fc2_patch_h1c.py")
assert "/r62_inv_word/" in P.replace("\\", "/"), "only my own source copy"
s = open(P).read()
MARK = "# R62-P1 (r62_inv_word_p1.py)"
if MARK in s:
    print("already patched:", P)
    sys.exit(0)
I = " " * 20
A1 = I + "lits_ = [((lambda w_=w_: w_), (lambda: None)) for w_ in srcw]\n"
assert s.count(A1) == 1, ("anchor 1", s.count(A1))
R1 = (I + '_r62s = bool((o.get("r62") or {}).get("cubstream")) and _sc   ' + MARK + "\n"
      + I + "_r62mb = {}\n\n"
      + I + "def _r62rel(w_):\n"
      + I + "    if _r62s and w_ in bw and w_ not in _r62mb:\n"
      + I + "        _r62mb[w_] = c.hmr(w_)\n"
      + I + "        c.free_([w_])\n"
      + I + "lits_ = [((lambda w_=w_: w_), (lambda w_=w_: _r62rel(w_))) for w_ in srcw]\n")
s = s.replace(A1, R1)
A2 = (I + "mbits = [c.hmr(w) for w in bw]\n"
      + I + "c.free_(bw)\n"
      + I + 'c.ops.append(("BEND", bid))\n')
assert s.count(A2) == 1, ("anchor 2", s.count(A2))
R2 = (I + "mbits = [(_r62mb[w] if w in _r62mb else c.hmr(w)) for w in bw]   " + MARK + "\n"
      + I + "c.free_([w for w in bw if w not in _r62mb])\n"
      + I + 'c.ops.append(("BEND", bid))\n')
s = s.replace(A2, R2)
A3 = (I + "bw2 = c.alloc(nbit)\n"
      + I + 'assert list(bw2) == list(bw), ("table wires", bw2, bw)\n'
      + I + "emit_table(c, kt, vals, zin, bw)\n")
assert s.count(A3) == 1, ("anchor 3", s.count(A3))
R3 = (I + "if _r62mb:                                     " + MARK + ": the same wires back\n"
      + I + "    _take(c, list(bw))\n"
      + I + "else:\n"
      + I + "    bw2 = c.alloc(nbit)\n"
      + I + '    assert list(bw2) == list(bw), ("table wires", bw2, bw)\n'
      + I + "emit_table(c, kt, vals, zin, bw)\n")
s = s.replace(A3, R3)
open(P, "w").write(s)
print("patched:", P)
