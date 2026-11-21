#!/usr/bin/env python3
"""invsqrt_fc2_h1lp.py -- FC2 (own code): angle-LP margin of a family-C ARITHMETIC (a word Y per input) without
building the circuit: the endpoint's phase is modelled as F1's ring (27 slots, token i at slot (i + k) mod 27, slot
angle 2^(p+27); wrap features j = 0..6 tunable, j = 7, 8, 9 fixed at T, S, Z), the LP is S1's lp_fast (read-only).
Calibrate on today's SC (its circuit margin is 5.6525e7), then compare variants of the higher-order term.
  invsqrt_fc2_h1lp.py CFG.json WD
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_fc_g8 as F8   # noqa: E402
import invsqrt_fc_gen as G   # noqa: E402
import invsqrt_sw_model as SM  # noqa: E402


def lp_margin(Y, k, tgt, tol):
    N = 27
    w = [((Y >> i) & 1).astype(np.uint8) for i in range(N)]
    M = 1 << 64
    fire, ang = [], []
    for p in range(N):
        i = (p - k) % N
        bits = np.zeros(len(Y), dtype=np.uint8)
        for ii in range(N):
            sel = i == ii
            bits[sel] = w[ii][sel]
        fire.append(bits)
        ang.append(1 << (p + 27))
    fixed = np.zeros(len(Y), dtype=np.uint64)
    for j in range(10):
        on = (k >= j + 1)
        i = np.where(on, j - k + N, 0)
        bits = np.zeros(len(Y), dtype=np.uint8)
        for ii in range(N):
            sel = on & (i == ii)
            bits[sel] = w[ii][sel]
        if j <= 6:
            fire.append(bits)
            ang.append(((1 << (j + 54)) - (1 << (j + 27))) % M)
        else:
            fixed += np.where(bits == 1, np.uint64((1 << (54 + j - 7 + 7)) % M), np.uint64(0))
    fire = np.array(fire, dtype=np.uint8)
    tot = fixed.copy()
    for a_, f_ in zip(ang, fire):
        tot += np.where(f_ == 1, np.uint64(a_), np.uint64(0))
    err = (tot - tgt).view(np.int64).astype(np.float64)
    m, _ = SM.lp_fast(err, tol, fire)
    return float(m)


def ref_main_g(U, g, d, m, bt):
    """FC's ref_main (read-only copy of its recurrence) with one extra dividend bit g at position m."""
    N = (np.int64(1) << (28 + m)) + (U.astype(np.int64) << (m + 1)) + (g.astype(np.int64) << m)
    R = N.copy()
    Qw = np.zeros(len(d), dtype=np.int64)
    for b in range(26, -1, -1):
        t = max(0, bt - b)
        dt = F8.tail_dt(d, t) if t else d
        s = np.where(R >= 0, 1, -1)
        R = R - s * (dt << b)
        Qw = Qw + s * (1 << b)
    return Qw - (R < 0).astype(np.int64)


def tper(m, n=34):
    return 0.56 * math.log2(n * 2 ** 64 / (math.pi * m)) + 5.3 if m > 0 else float("nan")


def main():
    cfg = dict(F8.CFG)
    cfg.update(json.load(open(sys.argv[1])))
    z = np.load(os.path.join(sys.argv[2], "chain.npz"))
    tgt, tol = z["target"].astype(np.uint64), z["tol"].astype(np.float64)
    m = cfg["m"]
    x = np.arange(1, F8.MAX_X + 1, dtype=np.int64)
    n, k, Q, d, delta = F8.ref_chain(x, m)
    a = F8.a_top(cfg)
    D = delta.astype(np.int64) << cfg["F1"]
    sg1, q1 = G.ref_div(D, d, a, F8.sched1_of(cfg))
    sg2, P = G.ref_div(q1 + cfg["c2"] * d, d, a, cfg["sched2"])
    y = -((P + 1) >> 1)
    yf = y.astype(np.float64)
    svals = [int(t) for t in os.environ.get("FC2_S", "6,8").split(",")]
    variants = {} if os.environ.get("FC2_S") else {
        "SC_today": F8.ref_U(x, dict(cfg))[0],
        "H1_exact": np.round(2.0 ** 28 * yf / (2.0 ** 28 - 3 * yf)).astype(np.int64) - 1}
    for s in svals:
        yt = (y >> s).astype(np.float64) + 0.5
        variants[f"trunc_s{s}"] = y + np.floor(3 * yt * yt * 2.0 ** (2 * s) / 2.0 ** 28 + 0.5).astype(np.int64) - 1
        hh = np.floor(3 * yt * yt * 2.0 ** (2 * s) / 2.0 ** 28 + 0.5).astype(np.int64)
        print(f"s {s}: h in [{int(hh.min())}, {int(hh.max())}], y >> s in [{int((y >> s).min())}, "
              f"{int((y >> s).max())}]", flush=True)
    # circuit form at s = 9: T = y >> 9 (10 bits signed), w = T xor sign (9-bit magnitude, T (T + 1) = w (w + 1)),
    # h = floor((3 S + K) / 2^10) with S = w (w + 1) from the KEPT partial products only (pairs i < j at column
    # i + j + 1 >= cut, every linear term w_i (2^2i + 2^i) kept), K = 512 + 3 E[dropped] (E[w_i w_j] = 1/4)
    for cut in [int(t) for t in os.environ.get("FC2_CUTS", "").split(",") if t]:
        T = y >> 9
        w = np.where(T < 0, -1 - T, T).astype(np.int64)
        S = np.zeros_like(w)
        drop_e, kept = 0.0, 0
        for i in range(9):
            wi = (w >> i) & 1
            S += wi * ((1 << (2 * i)) + (1 << i))
            for j in range(i + 1, 9):
                if i + j + 1 >= cut:
                    S += wi * ((w >> j) & 1) * (1 << (i + j + 1))
                    kept += 1
                else:
                    drop_e += (1 << (i + j + 1)) / 4.0
        K = int(round(512 + 3 * drop_e))
        h = (3 * S + K) >> 10
        variants[f"cut{cut}_ands{kept}"] = y + h - 1
        print(f"cut {cut}: kept pairs {kept}, K {K}, h in [{int(h.min())}, {int(h.max())}]", flush=True)
        if os.environ.get("FC2_DUMP"):
            # R36-INV-CORR's DUMP: one garbage bit g left on MAIN's dividend position m (half a unit of U)
            c = ((y & 511) + h) >> 9
            variants[f"cut{cut}_dump_c0"] = (y + h - 1, (c & 1).astype(np.int64))
            variants[f"cut{cut}_dump_par"] = (y + h - 1, (y & 1).astype(np.int64))
            print(f"carry values {sorted(set(np.unique(c).tolist()))}", flush=True)
    # +-1 digit form at s = 9 (R36-INV-CORR item 2): 2T + 1 = sum s_i 2^i, s_i = 2 t'_i - 1 (t' = T + 2^9, 10 bits),
    # (2T + 1)^2 = 1023^2 - 4 sum_{i<j} (t'_i xor t'_j) 2^(i+j); keep the pairs with i + j >= cut, the dropped XOR
    # bits by their mean 1/2; h = floor((3 (2T + 1)^2 + 2048) / 4096) (= the float form when nothing is dropped)
    for cut in [int(t) for t in os.environ.get("FC2_PM_CUTS", "").split(",") if t]:
        Tp = (y >> 9) + 512
        X = np.zeros_like(Tp)
        drop_e, kept = 0.0, 0
        for i in range(10):
            for j in range(i + 1, 10):
                if i + j >= cut:
                    X += (((Tp >> i) ^ (Tp >> j)) & 1) << (i + j)
                    kept += 1
                else:
                    drop_e += (1 << (i + j)) / 2.0
        sq = 1023 * 1023 - 4 * X - int(round(4 * drop_e))
        for dk in [int(t) for t in os.environ.get("FC2_PM_DK", "0").split(",")]:   # rounding constant, 1/4096 units
            h = (3 * sq + 2048 + dk) >> 12
            variants[f"pm{cut}_pairs{kept}" + (f"_dk{dk}" if dk else "")] = y + h - 1
            if os.environ.get("FC2_DUMP") and not dk:
                c = ((y & 511) + h) >> 9
                variants[f"pm{cut}_dump_c0"] = (y + h - 1, (c & 1).astype(np.int64))
                variants[f"pm{cut}_dump_par"] = (y + h - 1, (y & 1).astype(np.int64))
        print(f"pm cut {cut}: kept pairs {kept}, h in [{int(h.min())}, {int(h.max())}]", flush=True)
    out = {}
    for bt in [int(t) for t in os.environ.get("FC2_BT", "6,7").split(",")]:
        for nm, U in variants.items():
            if nm == "SC_today" and bt == 7:
                continue
            if isinstance(U, tuple):
                Y = ref_main_g(U[0], U[1], d, m, bt)
            else:
                Y = F8.ref_main(U, d, m, bt)
            mg = lp_margin(Y, k, tgt, tol)
            out[f"{nm}_bt{bt}"] = dict(margin=f"{mg:.4g}", tper=round(tper(mg), 3))
    print(json.dumps(out))


if __name__ == "__main__":
    main()
