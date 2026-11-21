#!/usr/bin/env python3
"""invsqrt_fc_inj.py -- FC: can the higher-order stage of family C (R33-INV-WIRE2) be an in-place map on U's
20 wires?  An in-place stage that does not read v computes U = f(y) of y = U1 alone; it must be injective on the
set of y values that occur (union over all inputs).  This tool measures, on all 3,133,452 inputs of the model
(re-pricer's best set L2/Q2, ho grpt2c):
  - distinct y, distinct U, U values with two or more distinct y (collisions of the model's map);
  - the same per root value v (the map restricted to one v);
  - the smallest c such that an injective integer map within +-c of the model's U exists on the occurring y
    (greedy order-preserving assignment; optimal for a monotone target).
Imports the designer's and the re-pricer's models read-only.
  invsqrt_fc_inj.py WD ['{"m":7,...}']
"""
import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


RP = _load("rp", "r33_inv_wire2_rp_tail.py")
CM = RP.CM

BEST = {"m": 7, "F1": 5, "r1": 0, "r2": 0, "ho": "grpt2c", "t": 13, "kg": 5, "thr": 25, "c0": 2, "bt": 5,
        "form": "line", "sched1": [0, 0, 1, 2, 3, 4, 5, 6], "sched2": [0, 0, 2, 3, 4, 5, 6, 7]}


def words(wd, prm):
    CM.setup(wd)
    ch = CM._STATE["ch"]
    m = int(prm["m"])
    cm = CM.chain(ch, m)
    RP.set_form(prm)
    p0 = {k: v for k, v in prm.items() if k not in ("form", "sched", "sched1", "sched2", "tag")}
    pn = dict(p0)
    pn["ho"] = "none"
    pn["c0"] = 0
    U1, _ = RP.ORIG_WORD(cm, pn)
    RP.set_form(prm)
    U, info = RP.word_cm_2c(cm, p0)
    return cm, U1.astype(np.int64), U.astype(np.int64), info


def min_c(ys, us):
    """smallest c with an order-preserving injective integer map g, |g(y) - u(y)| <= c (us sorted by ys)."""
    for c in range(0, 4097):
        prev = None
        ok = True
        lo = us - c
        hi = us + c
        g = lo.copy()
        # vectorised greedy: g_i = max(lo_i, g_{i-1} + 1) = max over j <= i of (lo_j + i - j)
        idx = np.arange(len(us))
        g = np.maximum.accumulate(lo - idx) + idx
        if np.all(g <= hi):
            return c, int((g - us).max()), int((us - g).max())
    return None, None, None


def main():
    wd = sys.argv[1]
    prm = dict(BEST)
    if len(sys.argv) > 2:
        prm.update(json.loads(sys.argv[2]))
    cm, U1, U, info = words(wd, prm)
    v = cm["v"]
    print(f"inputs {len(U1)}; y=U1 range {U1.min()}..{U1.max()}; U range {U.min()}..{U.max()}; "
          f"U err vs ideal {info['U_err']}", flush=True)
    # the model's map y -> U as a function of y: is it a function (same y -> same U)?
    pair = np.unique(np.stack([U1, U]), axis=1)
    ys, inv = np.unique(pair[0], return_inverse=True)
    print(f"distinct y {len(ys)}; distinct (y, U) pairs {pair.shape[1]} (equal: U is a function of y)", flush=True)
    order = np.argsort(pair[0], kind="stable")
    yo, uo = pair[0][order], pair[1][order]
    mono = bool(np.all(np.diff(uo) >= 0))
    du, cnt = np.unique(uo, return_counts=True)
    coll = int((cnt > 1).sum())
    print(f"distinct U {len(du)}; U values with >= 2 distinct y: {coll} (max {int(cnt.max())} y per U); "
          f"U monotone in y: {mono}", flush=True)
    # per v
    bad_v = 0
    worst = 0
    for vv in range(int(v.max()) + 1):
        sel = v == vv
        if not sel.any():
            continue
        yv = np.unique(U1[sel])
        pv = np.unique(np.stack([U1[sel], U[sel]]), axis=1)
        nu = len(np.unique(pv[1]))
        if nu < len(yv):
            bad_v += 1
            worst = max(worst, len(yv) - nu)
        sp = np.diff(yv)
    print(f"per v: roots with a collision {bad_v} of {int(v.max()) + 1} (worst {worst} lost values)", flush=True)
    # spacing of y per v (the reachable lattice)
    sps = []
    for vv in (0, 31, 63, 95, 127):
        yv = np.unique(U1[v == vv])
        sps.append((vv, int(np.diff(yv).min()), float(np.median(np.diff(yv)))))
    print(f"y spacing per v (v, min, median): {sps}", flush=True)
    # density of the union: fraction of integers in [ymin, ymax] that occur, by sub-range
    lo_, hi_ = int(yo.min()), int(yo.max())
    for a, b in ((lo_, lo_ // 2), (lo_ // 2, lo_ // 4), (lo_ // 4, 0), (0, hi_ // 4), (hi_ // 4, hi_ // 2), (hi_ // 2, hi_)):
        n = int(((yo >= a) & (yo < b)).sum())
        print(f"  y in [{a}, {b}): {n} distinct of {b - a} integers ({n / max(1, b - a):.3f})", flush=True)
    c, up, dn = min_c(yo, uo)
    print(f"smallest c for an injective map within +-c of the model on the occurring y: {c} "
          f"(max up {up}, max down {dn})", flush=True)
    print("RESULT " + json.dumps(dict(prm=prm, distinct_y=int(len(ys)), distinct_U=int(len(du)), collisions=coll,
                                      monotone=mono, bad_v=bad_v, min_c=c)), flush=True)


if __name__ == "__main__":
    main()
