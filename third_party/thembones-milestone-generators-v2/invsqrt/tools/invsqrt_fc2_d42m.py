#!/usr/bin/env python3
"""invsqrt_fc2_d42m.py -- FC2 (own code): LP margin of the D42 layout (R36-INV-CORR ledger item 3; lead 16:27Z):
DIV2 one tail row shorter (y2: 18 digits, units of 2), the higher-order term at half resolution on y2, and the
carry's low bit c0 left on U's own unit wire (no held wire, no dump position: MAIN unchanged, br 8).
  y2 = -(P2 + 1) / 2 with P2 = FC's ref_div(q1, d, 8, sched2[:-1]),  T = y2 >> 8 (10 bits), L = y2 & 255,
  t' = T + 512, x_ij = t'_i xor t'_j; parts of 3 = 1 + 2 at hh-unit positions i + j - 11 and i + j - 10,
  kept iff pos >= -g (and not dropped);  acc = 2^g L + K - sum x 2^(pos + g) + 2^g c3h,
  c3h = the cubic in hh units (T_8 - T_9 halves: +-1 at hh position -1, i.e. 2^(g-1) in acc units),
  hl = acc >> g, C = hl >> 8, U2 = 256 T + hl,  U = 2 U2 + gb,  gb = C & 1 (or 1 - (C & 1)).
SPEC.json: list of {"name", "g", "koff": [...], "cub": 0|1, "gb": 1|-1, "drop_parts": [[i,j,part],...], "bt": [7]}
  invsqrt_fc2_d42m.py CFG.json WD SPEC.json OUT.jsonl
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_fc_g8 as F8                 # noqa: E402  (FC, read-only)
import invsqrt_fc_gen as G                 # noqa: E402
from invsqrt_fc2_h1lp import lp_margin, tper   # noqa: E402


def stage_d42(y2, g, koff, cub=1, gb=1, drop_parts=(), bias=()):
    """bias: list of (threshold literal, amount in acc units): acc -= amount * literal, literal = w_8 (|T| >= 256)
    or w_8 OR w_7 (|T| >= 128) etc. given as 'w8', 'w7', 'w87' (w_8 or w_7): a T-dependent constant that keeps
    the carry non-negative where the term is small."""
    T = y2 >> 8
    L = y2 & 255
    tp = T + 512
    acc = L.astype(np.int64) << g
    kept, dropped_mean = [], 0.0
    dset = {tuple(q) for q in drop_parts}
    for i in range(10):
        for j in range(i + 1, 10):
            x = ((tp >> i) ^ (tp >> j)) & 1
            for part, pos in ((1, i + j - 11), (2, i + j - 10)):
                if pos >= -g and (i, j, part) not in dset:
                    acc -= x << (pos + g)
                    kept.append((i, j, part, pos))
                else:
                    dropped_mean += 0.5 * 2.0 ** (pos + g)
    # hh = h / 2: h = 3 (2T+1)^2 / 4096 = 3 * 1023^2 / 4096 - 3 X / 1024 (U units); minus 1 as in H1 (U = q - 1)
    knat = 2.0 ** g * ((3 * 1023 ** 2 / 4096.0) / 2.0) - dropped_mean
    K = int(round(knat)) + koff
    acc += K
    if cub:
        c3 = ((T >> 8) & 1) - ((T >> 9) & 1)          # T_8 - T_9 = sign(T) [|T| >= 256] (U units: +-1)
        acc += c3 << (g - 1)                          # half of it in hh units
    w = np.where(T < 0, -1 - T, T)
    for lit, amt in bias:
        lv = {"w8": (w >> 8) & 1, "w7": (w >> 7) & 1, "w87": ((w >> 8) | (w >> 7)) & 1,
              "w876": ((w >> 8) | (w >> 7) | (w >> 6)) & 1}[lit]
        acc -= lv * int(amt)
    hl = acc >> g
    C = hl >> 8
    U2 = (T << 8) + hl
    gbit = (C & 1) if gb == 1 else (1 - (C & 1))
    return 2 * U2 + gbit, C, kept, K


def main():
    cfg = dict(F8.CFG)
    cfg.update(json.load(open(sys.argv[1])))
    z = np.load(os.path.join(sys.argv[2], "chain.npz"))
    tgt, tol = z["target"].astype(np.uint64), z["tol"].astype(np.float64)
    specs = json.load(open(sys.argv[3]))
    out = open(sys.argv[4], "a")
    m = cfg["m"]
    x = np.arange(1, F8.MAX_X + 1, dtype=np.int64)
    n, k, Q, d, delta = F8.ref_chain(x, m)
    a = F8.a_top(cfg)
    _, q1 = G.ref_div(delta.astype(np.int64) << cfg["F1"], d, a, F8.sched1_of(cfg))
    _, P2 = G.ref_div(q1 + cfg["c2"] * d, d, a, cfg["sched2"][:-1])
    y2 = -((P2 + 1) >> 1)
    print("y2 range", int(y2.min()), int(y2.max()), flush=True)
    for sp in specs:
        for koff in sp.get("koff", [0]):
            U, C, kept, K = stage_d42(y2, sp["g"], koff, sp.get("cub", 1), sp.get("gb", 1), sp.get("drop_parts", ()),
                                      [tuple(b) for b in sp.get("bias", ())])
            for bt in sp.get("bt", [7]):
                Y = F8.ref_main(U, d, m, bt)
                mg = lp_margin(Y, k, tgt, tol)
                rec = dict(name=sp["name"], g=sp["g"], koff=koff, K=K, cub=sp.get("cub", 1), gb=sp.get("gb", 1), bias=sp.get("bias"),
                           bt=bt, parts=len(kept), margin=float("%.4g" % mg), tper=round(tper(mg), 3),
                           C=[int(C.min()), int(C.max())], U=[int(U.min()), int(U.max())])
                print(json.dumps(rec), flush=True)
                out.write(json.dumps(rec) + "\n")
                out.flush()


if __name__ == "__main__":
    main()
