"""I39 COPY (second holder) of invsqrt_fk_relabel.py with ONE change: a busy segment ends only where the wire is zero
on every input under SEVERAL outcome settings at once (HMR outcomes all 0, all 1, random seeds 1 and 2): F1's original
used outcomes 0 only, and its relabel of the 39-wire hold file broke the outcome-1 paths (x not restored) [M I39].
--- F1's docstring ---
invsqrt_fk (F1): relabel the wires of a kmx circuit so that its qubit count equals its peak live width.

A wire is BUSY at every op that touches it and between two consecutive touches whenever it is non-zero on some
input (exact over all 3,133,452 inputs, HMR outcomes 0; conditioned blocks are value-neutral repairs and are
treated as one macro op that touches every wire used inside it). The register wires (REGISTER / APPEND_TO_
REGISTER lines: the 22 input wires) keep their label for their first and last busy segment (they hold x at the
start and the end); every other busy segment gets the lowest label that is free for its whole span (first fit
by start time). The relabelled text is the same circuit up to a time-dependent wire permutation, so its phase,
its Toffoli expectation and its exact proof are unchanged.
  python invsqrt_fk_relabel.py IN.kmx OUT.kmx
Prints: peak live width, labels before and after.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import invsqrt_lib as L          # noqa: E402
import invsqrt_gen_circ as C     # noqa: E402

MAX_X = 3133452


def segments(header, ops):
    """Busy segments [(wire, first_op, last_op)] over macro ops (a conditioned block is one macro op)."""
    nq = 1 + max(q for op in ops for q in op[1])
    # macro ops: (first index, last index, wires touched)
    macro, depth, start, wires = [], 0, None, set()
    for i, op in enumerate(ops):
        if op[0] == "PUSH_CONDITION":
            if depth == 0:
                start, wires = i, set()
            depth += 1
        elif op[0] == "POP_CONDITION":
            depth -= 1
            if depth == 0:
                macro.append((start, i, wires))
        elif depth:
            wires |= set(op[1])
        else:
            macro.append((i, i, set(op[1])))
    xs = np.arange(1, MAX_X + 1, dtype=np.uint64)
    sims = [C.Sim(xs, nq, hmr="zero"), C.Sim(xs, nq, hmr="one"), C.Sim(xs, nq, hmr="rand", seed=1),
            C.Sim(xs, nq, hmr="rand", seed=2)]           # I39: every outcome setting must see the wire at zero
    regs = set()
    for h in header:
        for t in h.split():
            if t.startswith("q") and t[1:].isdigit():
                regs.add(int(t[1:]))
    busy_from = {}                                # wire -> start op of its open segment
    segs = []
    for q in regs:
        busy_from[q] = 0                          # the input register holds x from the start
    for a, b, ws in macro:
        for q in ws:
            busy_from.setdefault(q, a)
        for sim in sims:
            sim.run(ops, a, b + 1)
        for q in ws:
            if all(not sim.q[q].any() if hasattr(sim.q[q], "any") else sim.q[q] == 0 for sim in sims):   # I39
                segs.append((q, busy_from.pop(q), b))
    end = len(ops) - 1
    for q, s in busy_from.items():
        segs.append((q, s, end))
    return segs, regs, nq


def relabel(segs, regs, nq):
    """First-fit interval colouring; register wires keep their first and last segment."""
    first = {}
    last = {}
    for q, s, e in segs:
        if q in regs:
            if q not in first or s < first[q][1]:
                first[q] = (q, s, e)
            if q not in last or e > last[q][2]:
                last[q] = (q, s, e)
    fixed = {seg: seg[0] for seg in list(first.values()) + list(last.values())}
    taken = {}                                    # label -> list of (s, e)
    for seg, lab in fixed.items():
        taken.setdefault(lab, []).append((seg[1], seg[2]))
    out = dict(fixed)
    for seg in sorted((s for s in segs if s not in fixed), key=lambda t: (t[1], t[2])):
        q, s, e = seg
        lab = 0
        while any(not (e < a or b < s) for a, b in taken.get(lab, [])):
            lab += 1
        taken.setdefault(lab, []).append((s, e))
        out[seg] = lab
    return out


def main(inp, outp):
    t0 = time.time()
    header, ops = L.parse_kmx(inp)
    ops = [op for op in ops if op[0] != "INFO_TAG"]
    segs, regs, nq = segments(header, ops)
    lab = relabel(segs, regs, nq)
    # per-op mapping: each (wire, op) -> label of the segment containing op
    by_wire = {}
    for (q, s, e), l in lab.items():
        by_wire.setdefault(q, []).append((s, e, l))
    for v in by_wire.values():
        v.sort()

    def label(q, i):
        for s, e, l in by_wire.get(q, []):
            if s <= i <= e:
                return l
        raise AssertionError(("wire used outside its segments", q, i))
    # ops inside a conditioned block use the label at the block's span (the block is inside one segment)
    new = []
    for i, (name, qs, bit, cond, ang) in enumerate(ops):
        new.append((name, tuple(label(q, i) for q in qs), bit, cond, ang))
    # peak live width over macro boundaries = max overlap of segments
    ev = sorted([(s, 1) for (q, s, e) in segs] + [(e + 1, -1) for (q, s, e) in segs])
    cur = peak = 0
    for _, d in ev:
        cur += d
        peak = max(peak, cur)
    L.write_kmx(outp, header, new)
    used = 1 + max(lab.values())
    print(f"{inp}: wires before {nq}, peak busy {peak}, labels after {used} ({time.time() - t0:.0f} s) -> {outp}")


if __name__ == "__main__":
    main(*sys.argv[1:])
