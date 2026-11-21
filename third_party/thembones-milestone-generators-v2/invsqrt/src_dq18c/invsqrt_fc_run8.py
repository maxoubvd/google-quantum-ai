#!/usr/bin/env python3
"""invsqrt_fc_run8.py -- FC: build the g8 (SC, any m) circuit, simulate it on ALL 3,133,452 inputs, angle LP (S1's
lp_fast, read-only import), re-simulate with the LP angles, write invsqrt.kmx + angles.json + report.json (own code).
  invsqrt_fc_run8.py WD OUTDIR ['{cfg overrides}']
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_fc_g8 as F8   # noqa: E402
import invsqrt_fc_kit as K   # noqa: E402
import invsqrt_sw_model as SM  # noqa: E402  (S1, read-only: lp_fast)


def tper(m, n):
    return 0.56 * math.log2(n * 2 ** 64 / (math.pi * m)) + 5.3


def main():
    wd, outdir = sys.argv[1], sys.argv[2]
    cfg = dict(F8.CFG)
    if len(sys.argv) > 3:
        cfg.update(json.loads(sys.argv[3]))
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    z = np.load(os.path.join(wd, "chain.npz"))
    tgt = z["target"].astype(np.uint64)
    tol = z["tol"].astype(np.float64)
    x = np.arange(1, F8.MAX_X + 1, dtype=np.int64)
    K.KMX_CBLOCKS = bool(cfg.get("kmx_cblocks"))
    o = F8.full_check(x, cfg, tgt_tol=(tgt, tol))
    c = o["c"]
    rep = dict(cfg={k: v for k, v in cfg.items() if not k.startswith("_")}, ops=len(c.ops), toffolis=c.toffolis_exp(), toffolis_static=c.toffolis(),
               wires=c.next, peak=c.peak, stages=c.stage_counts(), peaks=c.stage_peaks, U_range=o["U_range"],
               U_in_range=o["U_in_range"], Y_ok=o["Y_ok"], Y_bad=o["Y_bad"], x_ok=o["x_ok"], anc_ok=o["anc_ok"],
               sign_ok=o["sign_ok"], margin_ideal=o["margin"], t_sim=round(time.time() - t0, 1))
    print("pass 1", json.dumps({k: v for k, v in rep.items() if k not in ("stages", "peaks", "cfg")}), flush=True)
    fire = o["fire"]
    regs = o["regs"]
    nfixT = 0
    if "epinfo" in regs:
        zl = regs["epinfo"]["zps"]
        assert fire.shape[0] == len(zl), (fire.shape[0], len(zl))
        prof = cfg.get("ep_profile", "generic")
        tun = [i for i, (oi, kind, idx, a) in enumerate(zl) if prof == "generic" or kind == "slot" or idx <= 6]
        nfixT = len(zl) - len(tun)
        nrot = len(tun)
        m_lp, dv = SM.lp_fast(o["err"], tol, fire[tun])
        ang = [1 << (p + 27) for p in range(27)]
        nfe = 10 if prof == "generic" else 7
        feat = [((1 << (p + 54)) - (1 << (p + 27))) % (1 << 64) for p in range(nfe)]
        rem = [0] * sum(1 for z_ in zl if z_[1] == "rem")
        for i, dd in zip(tun, dv):
            oi, kind, idx, a = zl[i]
            val = int((int(a) + int(dd)) % (1 << 64))
            if kind == "slot":
                ang[idx] = val
            elif kind == "rem":
                rem[idx] = val
            else:
                feat[idx] = val
        if rem:
            cfg["rem_ang"] = [str(a) for a in rem]
            rep["cfg"]["rem_ang"] = cfg["rem_ang"]
    else:
        nrot = fire.shape[0]
        m_lp, dv = SM.lp_fast(o["err"], tol, fire)
        base = [c.ops[i][2] for i in regs["zps"]]
        newang = [int((b + int(dd)) % (1 << 64)) for b, dd in zip(base, dv)]
        ang, feat = newang[:27], newang[27:]
    print(f"LP margin {m_lp:.6g} over {nrot} rotations ({nfixT} fixed T)", flush=True)
    o2 = F8.full_check(x, cfg, ang=ang, feat=feat, tgt_tol=(tgt, tol))
    c2 = o2["c"]
    m2 = o2["margin"]
    tp = tper(m2, nrot) if m2 > 0 else float("nan")
    cost = c2.next * (c2.toffolis_exp() + nrot * tp / 4.0 + nfixT / 4.0) if m2 > 0 else float("inf")
    rep.update(margin_lp=m_lp, margin_lp_check=m2, x_ok2=o2["x_ok"], anc_ok2=o2["anc_ok"], sign_ok2=o2["sign_ok"],
               n_rot=nrot, T_per_rot=tp, predicted_cost=cost, t_total=round(time.time() - t0, 1))
    print("pass 2", json.dumps({k: v for k, v in rep.items() if k not in ("stages", "peaks", "cfg")}), flush=True)
    with open(os.path.join(outdir, "invsqrt.kmx.tmp"), "w") as f:
        f.write(K.kmx_text(c2, regs["X"]))
    os.replace(os.path.join(outdir, "invsqrt.kmx.tmp"), os.path.join(outdir, "invsqrt.kmx"))
    json.dump(dict(ang=[str(a) for a in ang], feat=[str(a) for a in feat], cfg=rep["cfg"]),
              open(os.path.join(outdir, "angles.json"), "w"))
    json.dump(rep, open(os.path.join(outdir, "report.json"), "w"), indent=1)
    print("RESULT", json.dumps(dict(cost=cost, Q=c2.next, tof=c2.toffolis_exp(), margin=m2, tper=tp)), flush=True)


if __name__ == "__main__":
    main()
