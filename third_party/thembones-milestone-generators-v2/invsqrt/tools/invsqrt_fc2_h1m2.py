#!/usr/bin/env python3
"""invsqrt_fc2_h1m2.py -- FC2 (own code): LP margin of the CIRCUIT arithmetic of patch h1 (the higher-order term
after DIV2 on FC's `nocorr` base), exactly as the stage computes it:
  y = U1 (19-bit two's complement: plain DIV1 + DIV2), T = y >> 9 (10-bit signed), L = y & 511, t' = T + 512.
  acc = L 2^g + K - sum over kept parts of x_ij 2^(pos + g),   x_ij = t'_i xor t'_j  (i < j <= 9),
        part 1 at pos = i + j - 10 and part 2 at pos = i + j - 9 (3 = 1 + 2), kept iff pos >= -g (and s >= smin),
        [+ cubic c3 2^g, c3 = round(5 (T + 1/2)^3 / 2^28) clipped to +-cubic]
  U = 512 T + (acc >> g); C = acc >> (g + 9) (the carry into T); dump bit = C & 1 on MAIN's dividend position 8.
  K = round(Knat) + koff, Knat = 2^g (3 1023^2 / 4096 - 1/2) - (mean of every dropped part) [= the series term
  rounded, minus 1 as in H1's U = q - 1]; koff = 2^g removes the -1.
SPEC.json: list of {"name", "g", "koff": [..], "cubic": 0|1|2, "bt": [..], "smin": 0, "dump": 1}
  invsqrt_fc2_h1m2.py CFG.json WD SPEC.json OUT.jsonl
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_fc_g8 as F8                 # noqa: E402  (FC, read-only)
import invsqrt_fc_gen as G                 # noqa: E402
from invsqrt_fc2_h1lp import lp_margin, ref_main_g, tper   # noqa: E402  (FC2's own)


def stage_word(y, g, koff, cubic=0, smin=0, ctau=None, drop_j=(), drop_parts=(), ctau_neg=()):
    """the stage's integer arithmetic (see module doc); returns U, C, and the parts kept.
    ctau: list of thresholds on w = |T| folded (w = T xor sign): c3 = sign(T) * #{tau in ctau: w >= tau} (the
    circuit's cheap cubic; [w >= 256] is the bare wire w_8). drop_j: rows j left out (by their mean);
    drop_parts: [i, j, part] (part 1 at i + j - 10, part 2 at i + j - 9) left out."""
    T = y >> 9
    L = y & 511
    tp = T + 512
    acc = L.astype(np.int64) << g
    kept, dropped_mean = [], 0.0
    dset = {tuple(q) for q in drop_parts}
    for i in range(10):
        for j in range(i + 1, 10):
            x = ((tp >> i) ^ (tp >> j)) & 1
            for part, pos in ((1, i + j - 10), (2, i + j - 9)):
                if pos >= -g and i + j >= smin and j not in drop_j and (i, j, part) not in dset:
                    acc -= x << (pos + g)
                    kept.append((i, j, pos))
                else:
                    dropped_mean += 0.5 * 2.0 ** (pos + g)
    knat = 2.0 ** g * (3 * 1023 ** 2 / 4096.0 - 0.5) - dropped_mean
    K = int(round(knat)) + koff
    acc += K
    if ctau:
        w = np.where(T < 0, -1 - T, T)
        c3 = sum((w >= t).astype(np.int64) for t in ctau)
        c3n = sum((w >= t).astype(np.int64) for t in ctau_neg) if ctau_neg else 0
        acc += np.where(T < 0, -c3 - c3n, c3) << g
    elif cubic:
        c3 = np.round(5.0 * (T + 0.5) ** 3 / 2.0 ** 28).astype(np.int64)
        acc += np.clip(c3, -cubic, cubic) << g
    hl = acc >> g
    return (T << 9) + hl, hl >> 9, kept, K


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
    D = delta.astype(np.int64) << cfg["F1"]
    ycache = {}

    def y_of(s1, s2):                                   # plain DIV1 + DIV2 with the given tail schedules
        key = (tuple(s1), tuple(s2))
        if key not in ycache:
            _, q1 = G.ref_div(D, d, a, list(s1))
            _, P = G.ref_div(q1 + cfg["c2"] * d, d, a, list(s2))
            ycache.clear()
            ycache[key] = -((P + 1) >> 1)
        return ycache[key]
    for sp in specs:
        y = y_of(sp.get("sched1") or F8.sched1_of(cfg), sp.get("sched2") or cfg["sched2"])
        for koff in sp.get("koff", [0]):
            U, C, kept, K = stage_word(y, sp["g"], koff, sp.get("cubic", 0), sp.get("smin", 0), sp.get("ctau"),
                                       sp.get("drop_j", ()), sp.get("drop_parts", ()), sp.get("ctau_neg", ()))
            gb = (C & 1) if sp.get("dump", 1) else np.zeros_like(C)
            for bt in sp.get("bt", [7]):
                Y = ref_main_g(U, gb, d, m, bt)
                mg = lp_margin(Y, k, tgt, tol)
                rec = dict(name=sp["name"], g=sp["g"], koff=koff, K=K, cubic=sp.get("cubic", 0), ctau=sp.get("ctau"), ctau_neg=sp.get("ctau_neg"),
                           smin=sp.get("smin", 0), bt=bt, parts=len(kept), margin=float("%.4g" % mg),
                           sched1=sp.get("sched1"), sched2=sp.get("sched2"),
                           tper=round(tper(mg), 3), C=[int(C.min()), int(C.max())],
                           U_err=[int((U - y).min()), int((U - y).max())])
                print(json.dumps(rec), flush=True)
                out.write(json.dumps(rec) + "\n")
                out.flush()


if __name__ == "__main__":
    main()
