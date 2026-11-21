#!/usr/bin/env python3
"""r54_inv_tail_dump.py -- R54-INV-TAIL (designer): THE CIRCUIT'S OWN phase function, dumped from FC's generator on
all 3,133,452 inputs, to calibrate my integer model's endpoint against it (my model, R37's ring, reads k9f's LP margin
6.8% low; FC's LP on the simulated circuit is the truth).  Loads MY copy of the source exactly as FC's runner
invsqrt_fx_xcrun.py does (same patches, same env), runs F8.full_check's steps and saves:
  word (the circuit's 27 word wires after MAIN), fire (one row per executed rotation), err, tol, the rotation list
  (kind, index, base angle), the tunable mask, FC's LP margin.
  r54_inv_tail_dump.py run SRC WD OUT.npz 'CFGJSON'        (argv order as FC's runner: mode, src, wd, out, cfg)
Safe to run: reads my copy; writes OUT.npz in my folder."""
import json
import os
import sys

import numpy as np

src = os.path.abspath(os.path.expanduser(sys.argv[2]))
sys.path.insert(0, src)
import invsqrt_fw_fcpolrun as R  # noqa: E402


def main():
    wd, outp = os.path.expanduser(sys.argv[3]), os.path.expanduser(sys.argv[4])
    F8, K, w = R.load(src)
    cfg = R.cfg_for(F8, src, json.loads(sys.argv[5]))
    if os.environ.get("FX_XC", "1") != "0":
        import invsqrt_fx_xc as XC
        XC.apply(sys.modules["invsqrt_fc_arith"], cfg)
    import invsqrt_sw_model as SM
    z = np.load(os.path.join(wd, "chain.npz"))
    tgt, tol = z["target"].astype(np.uint64), z["tol"].astype(np.float64)
    x = np.arange(1, F8.MAX_X + 1, dtype=np.int64)
    c, regs = F8.build_full(cfg, None, None)
    sim = K.Sim(len(x), c.next + 2)
    sim.fire = []
    sim.set_word(regs["X"], x)
    nf = regs["nfwd"]
    sim.run(c.ops, 0, nf)
    word = sim.get_word(regs["word"])
    extra = {}
    for nm in ("F", "rest", "held1"):
        if nm in regs:
            try:
                extra["reg_" + nm] = np.array([sim.q[w_].astype(np.uint8) for w_ in regs[nm]])
            except Exception:                                   # noqa: BLE001
                pass
    sim.run(c.ops, nf, None)
    x_ok = bool(np.array_equal(sim.get_word(regs["X"]), x))
    others = [w_ for w_ in range(c.next) if w_ not in regs["X"]]
    anc_ok = bool(not sim.q[others].any())
    tot = sim.ang + np.where(sim.sign, np.uint64(1 << 63), np.uint64(0))
    err = (tot - tgt).view(np.int64).astype(np.float64)
    fire = np.array(sim.fire, dtype=np.uint8)
    zl = regs["epinfo"]["zps"] if "epinfo" in regs else []
    prof = cfg.get("ep_profile", "generic")
    tun = [i for i, (oi, kind, idx, a) in enumerate(zl) if prof == "generic" or kind == "slot" or idx <= 6]
    m_lp, dv = SM.lp_fast(err, tol, fire[tun])
    print("toffolis %.3f peak %d x_ok %s anc_ok %s rotations %d (tunable %d) LP margin %.0f" % (
        c.toffolis_exp(), c.peak, x_ok, anc_ok, fire.shape[0], len(tun), m_lp), flush=True)
    print("rotation list:", [(kind, idx) for (oi, kind, idx, a) in zl], flush=True)
    np.savez_compressed(outp, word=word, fire=fire, err=err, tun=np.array(tun), kinds=np.array([k_ for _, k_, _, _ in zl]),
                        idxs=np.array([i_ for _, _, i_, _ in zl]), angs=np.array([str(a_) for _, _, _, a_ in zl]),
                        margin=m_lp, dv=np.array(dv, dtype=np.float64), **extra)
    print("saved", outp, flush=True)


if __name__ == "__main__":
    main()
