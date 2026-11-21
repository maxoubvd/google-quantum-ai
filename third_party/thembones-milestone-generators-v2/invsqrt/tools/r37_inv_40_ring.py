#!/usr/bin/env python3
"""r37_inv_40_ring.py -- R37 designer INV-40 (own code): the set's own arithmetic (FC's frozen w41d: FC2's h1b term
through its patched ref_U, read-only) on rings of 27, 26 and 25 slots: MAIN stopped 0, 1 or 2 rows early, the
endpoint as a ring of N slots for the word bits 27 - N .. 26 (N + 7 rotations), converged LP on all inputs.
  r37_inv_40_ring.py SRC CFG.json WD
"""
import json
import math
import os
import sys

import numpy as np

src = os.path.abspath(os.path.expanduser(sys.argv[1]))
cfgp, wd = sys.argv[2], sys.argv[3]
sys.argv = [sys.argv[0], src]
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r37_inv_40_model as M   # noqa: E402  (puts SRC on the path)
import invsqrt_fw_fcpolrun as W  # noqa: E402

F8, K, _ = W.load(src)
cfg = W.cfg_for(F8, src, json.load(open(cfgp)))
c, regs = F8.build_full(cfg)                 # applies FC2's patches (ref_U becomes the circuit's own arithmetic)
x = np.arange(1, F8.MAX_X + 1, dtype=np.int64)
# the dividend U as the circuit itself holds it when MAIN starts (FC's ref_U does not carry every FC2 option)
i_main = max(t[0] for t in c.tags if t[1] == "MAIN")
sim = K.Sim(len(x), c.next + 2)
sim.fire = []
sim.set_word(regs["X"], x)
sim.run(c.ops, 0, i_main)
U = sim.get_word(regs["U"], signed=True).astype(np.int64)
sim.run(c.ops, i_main, regs["nfwd"])
Yc = sim.get_word(regs["word"]).astype(np.int64)
n_, k, Q_, d, delta_ = F8.ref_chain(x, cfg["m"])
print("word of the circuit = ref_main(U read at MAIN's start):", bool(np.array_equal(Yc, F8.ref_main(U, d, cfg["m"],
      cfg["bt"]))), flush=True)
del sim
z = np.load(os.path.join(wd, "chain.npz"))
tgt, tol = z["target"].astype(np.uint64), z["tol"].astype(np.float64)
m, bt = cfg["m"], cfg["bt"]
print("U range", int(U.min()), int(U.max()), "m", m, "bt", bt, flush=True)
for N, bts in ((27, (bt,)), (26, (bt, bt + 1)), (25, (bt,))):
    for b_ in bts:
        bl = 27 - N
        Y = M.main_rows(U, d, m, b_, bl, "zero" if bl else "std")
        fire, ang, fixed = M.ring_rows(Y, k) if N == 27 else M.ring_rows_n(Y, k, N)
        tot = fixed.copy()
        for a_, f_ in zip(ang, fire):
            tot += np.where(f_ == 1, np.uint64(a_), np.uint64(0))
        err = (tot - tgt).view(np.int64).astype(np.float64)
        mg, conv, bound = M.lp_conv(err, tol, fire)
        nrot = fire.shape[0]
        tp = M.tper(mg, nrot)
        print("ring %d (rows 26..%d, b_t %d): LP margin %.6g converged %s, %d rotations, T per rotation %.3f, "
              "rotations %.2f Toffoli-eq per shot" % (N, bl, b_, mg, conv, nrot, tp, nrot * tp / 4.0 if mg > 0 else
                                                      float("nan")), flush=True)

for slot in ("rest", "nrest", "zero", "one"):
    Y = M.main_rows(U, d, m, bt, 1, slot)
    fire, ang, fixed = M.ring_rows(Y, k)
    tot = fixed.copy()
    for a_, f_ in zip(ang, fire):
        tot += np.where(f_ == 1, np.uint64(a_), np.uint64(0))
    err = (tot - tgt).view(np.int64).astype(np.float64)
    mg, conv, bound = M.lp_conv(err, tol, fire)
    tp = M.tper(mg, 34)
    print("ring 27, rows 26..1, slot 0 = %s: LP margin %.6g converged %s, T per rotation %.3f, rotations %.2f" % (
        slot, mg, conv, tp, 34 * tp / 4.0 if mg > 0 else float("nan")), flush=True)

# all 27 rows kept, the ring on the word bits 1 .. 26 only (bit 0 = the final sign is not a slot), U shifted by a
# constant (FC2's delta knob moves the term by half units; here whole units of U)
for dU in (0, -1, 1, -2, 2):
    Y = F8.ref_main(U + dU, d, m, bt)
    for N in (26, 25):
        fire, ang, fixed = M.ring_rows_n(Y, k, N)
        tot = fixed.copy()
        for a_, f_ in zip(ang, fire):
            tot += np.where(f_ == 1, np.uint64(a_), np.uint64(0))
        err = (tot - tgt).view(np.int64).astype(np.float64)
        mg, conv, bound = M.lp_conv(err, tol, fire)
        nrot = fire.shape[0]
        tp = M.tper(mg, nrot)
        print("rows 26..0, ring %d on bits %d..26, U %+d: LP margin %.6g converged %s, %d rotations, T per rotation "
              "%.3f, rotations %.2f" % (N, 27 - N, dU, mg, conv, nrot, tp, nrot * tp / 4.0 if mg > 0 else
                                        float("nan")), flush=True)
