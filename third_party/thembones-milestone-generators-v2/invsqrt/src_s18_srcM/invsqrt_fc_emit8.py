#!/usr/bin/env python3
"""invsqrt_fc_emit8.py -- FC: rebuild a g8 circuit from its published angles.json (+ cfg) and write invsqrt.kmx and
invsqrt.kmb (own code). Byte-identical to the run that produced the angles when the generator is unchanged.
  invsqrt_fc_emit8.py OUTDIR ANGLES.json [CFG.json]     (CFG defaults to the cfg stored inside ANGLES.json)
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_fc_g8 as F8   # noqa: E402
import invsqrt_fc_kit as K   # noqa: E402

sys.path.insert(0, os.path.join(os.path.expanduser("~/qcc/data"), "kickmix"))


def main():
    out, angp = sys.argv[1], sys.argv[2]
    a = json.load(open(angp))
    cfg = dict(F8.CFG)
    cfg.update(json.load(open(sys.argv[3])) if len(sys.argv) > 3 else a.get("cfg", {}))
    if "rem_ang" in a.get("cfg", {}):
        cfg["rem_ang"] = a["cfg"]["rem_ang"]          # LP angles of the extra remainder rotations
    ang = [int(t) for t in a["ang"]]
    feat = [int(t) for t in a["feat"]]
    K.KMX_CBLOCKS = bool(cfg.get("kmx_cblocks"))
    c, regs = F8.build_full(cfg, ang, feat)
    os.makedirs(out, exist_ok=True)
    p = os.path.join(out, "invsqrt.kmx")
    with open(p + ".tmp", "w") as f:
        f.write(K.kmx_text(c, regs["X"]))
    os.replace(p + ".tmp", p)
    import kickmix as km
    kc = km.Circuit.from_file(p, format="kmx")
    kc.to_file(os.path.join(out, "invsqrt.kmb.tmp"), format="kmb")
    os.replace(os.path.join(out, "invsqrt.kmb.tmp"), os.path.join(out, "invsqrt.kmb"))
    print(f"emit8 {out}: toffolis {c.toffolis()} wires {c.next}; kickmix qubits {kc.num_qubits}")


if __name__ == "__main__":
    main()
