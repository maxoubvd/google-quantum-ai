"""invsqrt_gen: experiment, not part of the generator: how much of the error budget does a STRUCTURED
curvature table cost? (own tool)

The public file builds its c2 from the address by a coarse table plus interpolation rows (research/
ws_invsqrt_decode.md 2b.1: 44 Toffolis against our 66 for the block that reloads c2, about 25 against our 45
for the first load). That only works if c2 IS such an interpolation. Here the fitted C2 of a parameter set is
replaced by its best piecewise-linear approximation on blocks of `blk` addresses (rounded to integers), C1 and
C0 are fitted again around it, and the whole pipeline is checked on every input.

usage: invsqrt_gen_c2exp.py blk=8 [frac=0] key=value ...     (key=value: model parameters)
  frac: extra fractional bits of the slope (the line is base + round(slope * i) with slope a multiple of 2^-frac)
"""
from __future__ import annotations

import math
import multiprocessing as mp
import sys

import numpy as np

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_gen_model as M


def structured(c2, blk, frac):
    out = list(c2)
    dev = 0.0
    for b0 in range(0, 256, blk):
        ys = np.array(c2[b0:b0 + blk], dtype=float)
        xs = np.arange(blk, dtype=float)
        sl, ic = np.polyfit(xs, ys, 1) if blk > 1 else (0.0, ys[0])
        if frac >= 0:
            sl = round(sl * (1 << frac)) / (1 << frac)
            ic = round(float(np.mean(ys - sl * xs)))
        for i in range(blk):
            v = int(math.floor(ic + sl * i + 0.5))
            dev = max(dev, abs(v - c2[b0 + i]))
            out[b0 + i] = v
    return out, dev


def main():
    opt = dict(blk=8, frac=-1)
    kw = {}
    for a in sys.argv[1:]:
        k, v = a.split("=")
        if k in opt:
            opt[k] = int(v)
        else:
            kw[k] = int(v)
    P = M.Params(**kw)
    tabs = M.fit_tables(P)
    c2s, dev = structured(tabs["C2"], opt["blk"], opt["frac"])
    print(f"blocks of {opt['blk']}, slope bits {opt['frac']}: largest change of a C2 entry {dev:.0f}", flush=True)
    M.C2_FIXED = c2s
    with mp.get_context("fork").Pool(3) as pool:
        segs = pool.map(M.fit_segment, [(P, a8) for a8 in range(256)], chunksize=8)
    t2 = dict(H=[s["H"] for s in segs], C0=[s["C0"] for s in segs], C1=[s["C1"] for s in segs],
              C2=[s["C2"] for s in segs], seg_margin=[s["margin"] for s in segs])
    assert t2["C2"] == c2s
    out, _, _ = M.check(P, t2, n_rot=41)
    out.pop("by_level")
    sm = np.array(t2["seg_margin"])
    print(out)
    print(f"segment margins: min {sm.min():.3f} at a8 = {int(sm.argmin())}, negative {int((sm < 0).sum())}, "
          f"median {np.median(sm):.3f}")


if __name__ == "__main__":
    main()
