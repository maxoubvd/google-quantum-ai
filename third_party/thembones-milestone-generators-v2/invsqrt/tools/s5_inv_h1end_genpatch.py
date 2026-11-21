#!/usr/bin/env python3
"""s5_inv_h1end_genpatch.py -- S5-INV-H1END: patch a COPY of FC's generator source (srcM) so that H1's word un-add
(invsqrt_fc2_patch_h1c.py, the non-w0slot path) LENDS T[6] (q11) and runs at one more clean wire.

Lend (env S5LEND=1 at build time; without it the patched source emits exactly what srcM emits):
  after _cin(True): X-measure T6 (outcome m) and repay (-1)^(m T6) by Cliffords from
      T6 = T5 ^ W5 ^ MAJ(T5 ^ T7, W3, W4)        [M on all 3,133,452 inputs at ops 18,028 / 18,498 of fc19_b312]
  i.e. Z T5, Z W5, CZ (T5, W3) (T7, W3) (T5, W4) (T7, W4) (W3, W4), all under m; free the wire; Cc = cap - live;
  the un-add then takes the cheapest form the generator has at that Cc (game / mixed / ...);
  before _cin(False): take T6's wire back and re-make it with ONE AND:
      T6 = T7 ^ W5 ^ (T5 ^ T7 ^ W3)(T5 ^ T7 ^ W4).
  python -I s5_inv_h1end_genpatch.py SRC_COPY     (edits SRC_COPY/invsqrt_fc2_patch_h1c.py in place; refuses srcM)"""
import os
import sys

src = os.path.abspath(os.path.expanduser(sys.argv[1]))
assert not src.rstrip("/").endswith("s18/srcM"), "never patch FC's own srcM"
fn = os.path.join(src, "invsqrt_fc2_patch_h1c.py")
txt = open(fn).read()
assert "S5LEND" not in txt, "already patched"

A1 = ("        _cin(True)\n"
      "    if _r62w0:\n"
      "        pass\n"
      "    elif NW.plan(max(0, nt - 2 - NW.free_prefix(srcw, nt, _cw)), Cc)[0] is not None and not sw.get(\"mixed\") == 2:\n")
B1 = ("        _cin(True)\n"
      "    _s5 = (not _r62w0) and os.environ.get(\"S5LEND\") == \"1\"      # S5-INV-H1END lend of T6 (q11)\n"
      "    if _s5:\n"
      "        _S5W, _t5, _t7 = T[6], T[5], T[7]\n"
      "        _w3, _w4, _w5 = srcw[3], srcw[4], srcw[5]\n"
      "        _m5 = c.hmr(_S5W)\n"
      "        c.z(_t5, _m5)\n"
      "        c.z(_w5, _m5)\n"
      "        for _a5, _b5 in ((_t5, _w3), (_t7, _w3), (_t5, _w4), (_t7, _w4), (_w3, _w4)):\n"
      "            c.cz(_a5, _b5, _m5)\n"
      "        c.free_([_S5W])\n"
      "        Cc = c.cap - c.live\n"
      "        print(\"S5LEND: lent q%d; T5 q%d T7 q%d W3 q%d W4 q%d W5 q%d; Cc %d; L %d; srcw %s; Rt %s\" % (\n"
      "            _S5W, _t5, _t7, _w3, _w4, _w5, Cc, max(0, nt - 2 - NW.free_prefix(srcw, nt, _cw)), srcw, Rt),\n"
      "            file=sys.stderr)\n"
      "    if _r62w0:\n"
      "        pass\n"
      "    elif NW.plan(max(0, nt - 2 - NW.free_prefix(srcw, nt, _cw)), Cc)[0] is not None and not sw.get(\"mixed\") == 2:\n")
A2 = ("    if not _r62w0:                                     # R62-P2 (r62_inv_word_p2.py)\n"
      "        _cin(False)\n")
B2 = ("    if _s5:                                            # S5-INV-H1END: re-make T6 by one AND\n"
      "        print(\"S5LEND: un-add form %s\" % how, file=sys.stderr)\n"
      "        _take(c, [_S5W])\n"
      "        for _a5 in (_t5, _t7):\n"
      "            c.cx(_a5, _w3)\n"
      "            c.cx(_a5, _w4)\n"
      "        c.and_compute(_w3, _w4, _S5W)\n"
      "        for _a5 in (_t5, _t7):\n"
      "            c.cx(_a5, _w3)\n"
      "            c.cx(_a5, _w4)\n"
      "        c.cx(_t7, _S5W)\n"
      "        c.cx(_w5, _S5W)\n"
      "    if not _r62w0:                                     # R62-P2 (r62_inv_word_p2.py)\n"
      "        _cin(False)\n")
assert txt.count(A1) == 1, ("anchor 1", txt.count(A1))
assert txt.count(A2) == 1, ("anchor 2", txt.count(A2))
txt = txt.replace(A1, B1).replace(A2, B2)
open(fn, "w").write(txt)
print("patched", fn)
