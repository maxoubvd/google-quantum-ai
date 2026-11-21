#!/usr/bin/env python3
"""r63_inv_ghost_chk.py -- R63-INV-GHOST (designer): a changed circuit against its base on ALL 3,133,452 inputs, and
the live-wire count at every op boundary.  The kit's exact simulator (tools/invsqrt_gen_circ.py Sim, read-only)
on the emitted .kmx op lists.
  check BASE.kmx NEW.kmx [N]     base with outcomes 0; the changed circuit with outcomes all 0, all 1 and two random
                                 seeds: x restored, every other wire zero, the phase equal to the base's on every
                                 input (so the phase does not depend on the outcomes and the phase function is the
                                 base's); mean Toffoli-like gates a shot.  N > 0: a sample of N inputs.
  live  NEW.kmx [N] [LO HI]      live wires (nonzero on some input, random outcomes) at every op boundary: the peak,
                                 the runs at or above LO (default 39) and, inside the window of the rotations, the
                                 run table.
Safe to run: reads two files, prints."""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_lib as L          # noqa: E402
import invsqrt_gen_circ as C     # noqa: E402

MAX_X = 3133452


def load(path):
    header, ops = L.parse_kmx(os.path.expanduser(path))
    ops = [op for op in ops if op[0] != "INFO_TAG"]
    nq = 1 + max(q for op in ops for q in op[1])
    return ops, nq


def inputs(n):
    if n <= 0:
        return np.arange(1, MAX_X + 1, dtype=np.uint64)
    rng = np.random.default_rng(5)
    return np.unique(np.concatenate([rng.integers(1, MAX_X + 1, n), np.arange(1, 300),
                                     np.arange(MAX_X - 300, MAX_X + 1)])).astype(np.uint64)


def final(sim, xs, nq):
    xok = bool(np.array_equal(sim.value(range(22)), xs))
    dirty = [w for w in range(22, nq) if sim.q[w]]
    return xok, dirty


def check(base, new, n):
    xs = inputs(n)
    ops0, nq0 = load(base)
    t0 = time.time()
    s0 = C.Sim(xs, nq0, hmr="zero")
    s0.run(ops0)
    ph0 = s0.phase_units()
    xok, dirty = final(s0, xs, nq0)
    print("base %s: %d inputs, %d wires, x restored %s, dirty wires %s, Toffoli-like %.3f a shot (outcomes 0) [%.0f s]"
          % (os.path.basename(os.path.dirname(base)), len(xs), nq0, xok, dirty, s0.tof / len(xs), time.time() - t0),
          flush=True)
    ops1, nq1 = load(new)
    allok = True
    tofs = {}
    for mode, seed in (("zero", 0), ("one", 0), ("rand", 1), ("rand", 2)):
        t0 = time.time()
        s1 = C.Sim(xs, nq1, hmr=mode, seed=seed)
        s1.run(ops1)
        xok, dirty = final(s1, xs, nq1)
        bad = int(np.count_nonzero(s1.phase_units() != ph0))
        ok = xok and not dirty and bad == 0
        allok &= ok
        tofs[(mode, seed)] = s1.tof / len(xs)
        print("new, outcomes %s%s: x restored %s, dirty wires %s, phase differs from the base on %d inputs, "
              "Toffoli-like %.3f a shot, %d wires [%.0f s]" % (mode, "" if mode != "rand" else " seed %d" % seed, xok,
                                                             dirty, bad, tofs[(mode, seed)], nq1, time.time() - t0),
              flush=True)
    print("CHECK %s: Toffoli-like a shot with random outcomes = %.3f (mean of the two seeds; the exact expectation "
          "is the generator's toffolis_exp)" % ("PASS" if allok else "FAIL",
                                                (tofs[("rand", 1)] + tofs[("rand", 2)]) / 2), flush=True)


def live(new, n, lo, hi_unused):
    xs = inputs(n)
    ops, nq = load(new)
    sim = C.Sim(xs, nq, hmr="rand", seed=3)
    counts = []
    zp = [i for i, op in enumerate(ops) if op[0] == "Z_POW"]
    for i in range(len(ops)):
        sim.run(ops, i, i + 1)
        counts.append(sum(1 for w in range(nq) if sim.q[w]))
    counts = np.array(counts)
    print("%d ops, %d wires in the file; peak live %d; ops at the peak %d; rotations at ops %d..%d" % (
        len(ops), nq, counts.max(), int((counts == counts.max()).sum()), zp[0], zp[-1]), flush=True)
    for thr in range(int(counts.max()), lo - 1, -1):
        print("  op boundaries with live >= %d: %d" % (thr, int((counts >= thr).sum())))
    # the END window: from the last op before the first rotation where live was at its entry plateau, to the
    # first op after the last rotation where live is back at that plateau; printed as a run table
    a, b = max(0, zp[0] - 700), min(len(ops), zp[-1] + 700)
    print("run table of the window ops %d..%d (start, end, live):" % (a, b))
    s = a
    for i in range(a + 1, b + 1):
        if i == b or counts[i] != counts[s]:
            if i - s >= 1:
                print("   %6d..%6d  live %d%s" % (s, i - 1, counts[s], "   <-- rotations" if any(s <= z <= i - 1 for z in (zp[0], zp[-1])) else ""))
            s = i
    w = counts[zp[0]:zp[-1] + 1]
    print("between the first and the last rotation: peak live %d, at or above 39: %d boundaries, at 40: %d" % (
        w.max(), int((w >= 39).sum()), int((w >= 40).sum())), flush=True)


def main():
    if sys.argv[1] == "check":
        check(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 0)
    else:
        n = int(sys.argv[3]) if len(sys.argv) > 3 else 0
        lo = int(sys.argv[4]) if len(sys.argv) > 4 else 39
        live(sys.argv[2], n, lo, None)


if __name__ == "__main__":
    main()
