#!/usr/bin/env python3
"""invsqrt_fc18_patch_r25c.py -- FC eighteenth holder: apply FC18-R25C (tools/invsqrt_fc18_r25c.py) to a source COPY:
copies the module in and adds two gated hooks (cfg "_fc18r25c": 1; inert without the key):
  invsqrt_fc2_k4fuse.py  _r63_chain_top: per-bit op marks; the chain block's undo built by chain_mirror
  invsqrt_fc17_skipe.py  r24_block's undo: no re-make of w7' / w8'; skip E's m8 part paid by Cliffords
  python3 invsqrt_fc18_patch_r25c.py SRC
Refuses a source outside .../invsqrt_fc/ (FC's own copies only)."""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def sub(path, old, new):
    s = open(path).read()
    assert s.count(old) == 1, (path, old[:60], s.count(old))
    open(path, "w").write(s.replace(old, new))


def main():
    src = os.path.abspath(os.path.expanduser(sys.argv[1]))
    assert "/invsqrt_fc/" in src, "FC's own source copies only"
    shutil.copy(os.path.join(HERE, "invsqrt_fc18_r25c.py"), os.path.join(src, "invsqrt_fc18_r25c.py"))
    k4 = os.path.join(src, "invsqrt_fc2_k4fuse.py")
    sub(k4, "    for bit, lt, lf, rf in _R63_CHAIN:\n        below, tgt = wire[bit - 1], wire[bit]\n",
        "    _m18 = {}                                    # FC18-R25C: where each bit's forward ops start\n"
        "    for bit, lt, lf, rf in _R63_CHAIN:\n        _m18[bit] = len(c.ops)\n"
        "        below, tgt = wire[bit - 1], wire[bit]\n")
    sub(k4, "    c.blocks[bid] = K_.Circ.mirror(fwd)\n",
        "    c.blocks[bid] = K_.Circ.mirror(fwd)\n"
        "    if cfg.get(\"_fc18r25c\"):                     # FC18-R25C (invsqrt_fc18_r25c.py): bits 7, 8 under g only\n"
        "        import invsqrt_fc18_r25c as _R25C\n"
        "        c.blocks[bid] = _R25C.chain_mirror(c, K_, c.ops[i0:_m18[7]], c.ops[_m18[7]:_m18[8]],\n"
        "                                           c.ops[_m18[8]:_m18[9]], c.ops[_m18[9]:i0 + len(fwd)], wire)\n")
    sk = os.path.join(src, "invsqrt_fc17_skipe.py")
    sub(sk, "    c.free.remove(t8)\n    c.live += 1\n    K4._ch1_w7(c, tgt, d, v, t8, erase=False)      # w7' re-made, then the m8 part of the debt\n"
            "    c.z(t8, m8)\n    c.free.remove(top)\n    c.live += 1\n"
            "    W._r24_form(c, tgt, d, v, top, erase=False)    # w8' re-made\n    cost += 2.0 + 0.5 * nccz[0]\n",
        "    c.free.remove(t8)\n    c.live += 1\n"
        "    _c18 = (sys.modules[\"invsqrt_fk_w24ret\"]._STATE.get(\"cfg\") or {}).get(\"_fc18r25c\")   # FC18-R25C\n"
        "    if _c18:                                       # w7' / w8' not re-made: the chain mirror makes them under g\n"
        "        import invsqrt_fc18_r25c as _R25C\n"
        "        _R25C.pay_w7_phase(c, tgt, d, v, m8)\n"
        "    else:\n"
        "        K4._ch1_w7(c, tgt, d, v, t8, erase=False)      # w7' re-made, then the m8 part of the debt\n"
        "        c.z(t8, m8)\n"
        "    c.free.remove(top)\n    c.live += 1\n"
        "    if not _c18:\n"
        "        W._r24_form(c, tgt, d, v, top, erase=False)    # w8' re-made\n"
        "    cost += (0.0 if _c18 else 2.0) + 0.5 * nccz[0]\n")
    print("FC18-R25C applied to", src)


if __name__ == "__main__":
    main()
