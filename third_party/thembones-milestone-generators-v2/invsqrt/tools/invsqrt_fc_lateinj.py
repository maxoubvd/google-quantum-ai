#!/usr/bin/env python3
"""invsqrt_fc_lateinj.py -- FC (fourth holder): the constant U - 1 injected INSIDE MAIN, at the row whose LSB is U's LSB
(row m + 1: R -= s d 2^b + 2^b). In the circuit that row's LSB cell drops out and the carry-in becomes NOT sig (free).
Question: is the word then ref_main(U - 1) (the model's off -1, margin 7.47e7 at 26 slots), and what LP margin does
the late-injected word read on the 26-slot ring? U is read from the circuit itself at MAIN's start (all inputs), as in
the R37 designer's tools/r37_inv_40_ring.py.
  invsqrt_fc_lateinj.py SRC CFG.json WD [ENVLESS]      (SRC: frozen src with F3's runner; WD: folder with chain.npz)
"""
import json
import os
import sys
import time

import numpy as np

src = os.path.abspath(os.path.expanduser(sys.argv[1]))
cfgp, wd = sys.argv[2], sys.argv[3]
sys.argv = [sys.argv[0], src]
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r37_inv_40_model as M   # noqa: E402  (puts SRC on the path)
import invsqrt_fw_fcpolrun as W  # noqa: E402

t0 = time.time()
F8, K, _ = W.load(src)
cfg = W.cfg_for(F8, src, json.load(open(cfgp)))
c, regs = F8.build_full(cfg)
x = np.arange(1, F8.MAX_X + 1, dtype=np.int64)
i_main = max(t[0] for t in c.tags if t[1] == "MAIN")
sim = K.Sim(len(x), c.next + 2)
sim.fire = []
sim.set_word(regs["X"], x)
sim.run(c.ops, 0, i_main)
U = sim.get_word(regs["U"], signed=True).astype(np.int64)
del sim
n_, k, Q_, d, delta_ = F8.ref_chain(x, cfg["m"])
m, bt = cfg["m"], cfg["bt"]
blast = int(cfg.get("_fc2_blast", 0))
print("U read, range", int(U.min()), int(U.max()), "m", m, "bt", bt, "blast", blast, "%.0fs" % (time.time() - t0),
      flush=True)
z = np.load(os.path.join(wd, "chain.npz"))
tgt, tol = z["target"].astype(np.uint64), z["tol"].astype(np.float64)


def rows(U, inj_row=None, inj=0, bl=1):
    """FC's ref_main with rows 26..bl (final sign weight 2^bl, low bits not slots) and an optional constant inj * 2^b
    added at row b = inj_row (after that row's subtraction)."""
    N = (np.int64(1) << (28 + m)) + (U.astype(np.int64) << (m + 1))
    R = N.copy()
    Qw = np.zeros(len(d), dtype=np.int64)
    for b in range(26, bl - 1, -1):
        t = max(0, bt - b)
        dt = F8.tail_dt(d, t) if t else d
        s_ = np.where(R >= 0, 1, -1)
        R = R - s_ * (dt << b)
        if b == inj_row:
            R = R + (np.int64(inj) << b)
        Qw = Qw + s_ * (1 << b)
    Y = Qw - ((R < 0).astype(np.int64) << bl)
    if bl:
        Y = Y - (Y & ((1 << bl) - 1))
    return Y


def lp(Y, N):
    fire, ang, fixed = M.ring_rows(Y, k) if N == 27 else M.ring_rows_n(Y, k, N)
    tot = fixed.copy()
    for a_, f_ in zip(ang, fire):
        tot += np.where(f_ == 1, np.uint64(a_), np.uint64(0))
    err = (tot - tgt).view(np.int64).astype(np.float64)
    mg, conv, bound = M.lp_conv(err, tol, fire)
    nrot = fire.shape[0]
    return mg, conv, nrot, (nrot * M.tper(mg, nrot) / 4.0 if mg > 0 else float("nan"))


N = 27 - blast if blast else 26
bl = blast if blast else 1
Ye = {off: rows(U + off, bl=bl) for off in (0, -1, -2, -4)}
print("check: rows(U) == ref model main_rows:", bool(np.array_equal(Ye[0], M.main_rows(U, d, m, bt, bl, "zero"))),
      flush=True)
cases = [("exact off 0", Ye[0], 0), ("exact off -1", Ye[-1], -1)]
for r_ in (m + 1, m + 2, m + 3):                   # -1 at the LSB of row r_ = U - 2^(r_ - m - 1)
    offu = -(1 << (r_ - m - 1))
    cases.append(("late row %d (U %+d)" % (r_, offu), rows(U, inj_row=r_, inj=-1, bl=bl), offu))
for name, Y, offu in cases:
    dif = int((Y != Ye[offu]).sum())
    mg, conv, nrot, rc = lp(Y, N)
    print("%-22s diff vs exact U%+d: %7d | ring %d margin %.5g conv %s rot %d rotations/shot %.2f"
          % (name, offu, dif, N, mg, conv, nrot, rc), flush=True)
print("done %.0fs" % (time.time() - t0))
