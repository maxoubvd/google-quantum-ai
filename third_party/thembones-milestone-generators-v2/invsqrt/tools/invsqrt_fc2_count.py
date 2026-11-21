#!/usr/bin/env python3
"""invsqrt_fc2_count.py -- FC2 (own code): build ONE configuration of FC's g8 (with the FC2 patches named in the
cfg, through FC's hook) and print its counts: expected Toffolis (a conditioned AND = 1/2), static Toffolis, wires,
peak, per-stage counts. Optionally the reference check on a sample (word against FC's reference, x and ancillas).
  invsqrt_fc2_count.py FOLDER ['{cfg overrides}'] [NSAMP]      (FOLDER: cfg.json + angles.json)
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_fc_g8 as F8   # noqa: E402  (FC, read-only)


def main():
    folder = sys.argv[1]
    cfg = dict(F8.CFG)
    cfg.update(json.load(open(os.path.join(folder, "cfg.json"))))
    if len(sys.argv) > 2 and sys.argv[2]:
        for k, v in json.loads(sys.argv[2]).items():
            if v is None:
                cfg.pop(k, None)
            else:
                cfg[k] = v
    a = json.load(open(os.path.join(folder, "angles.json")))
    ang, feat = [int(t) for t in a["ang"]], [int(t) for t in a["feat"]]
    out = {}
    if len(sys.argv) > 3:
        x = F8.sample_x(int(sys.argv[3]))
        o = F8.full_check(x, dict(cfg), ang=ang, feat=feat)
        c = o["c"]
        out.update(Y_bad=o["Y_bad"], x_ok=o["x_ok"], anc_ok=o["anc_ok"], margin=o["margin"])
    else:
        c, regs = F8.build_full(dict(cfg), ang, feat)
    out.update(exp=c.toffolis_exp(), static=c.toffolis(), wires=c.next, peak=c.peak,
               stages=c.stage_counts(), stage_peaks=getattr(c, "stage_peaks", {}), fc2=cfg.get("fc2"),
               fc2_src=cfg.get("fc2_src"))
    print(json.dumps(out))


if __name__ == "__main__":
    main()
