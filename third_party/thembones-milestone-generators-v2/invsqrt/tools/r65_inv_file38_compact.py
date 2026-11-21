#!/usr/bin/env python3
"""r65_inv_file38_compact.py -- R65-INV-FILE38 (designer): WIRE COMPACTION OF AN EMITTED CIRCUIT.

The generator plans every stage under a cap but ALLOCATES wires by its own book: a wire that is |0> on every input
can still be "taken" in the book, so a frame that never holds more than W nonzero wires is emitted with W + 1 or
W + 2 qubits (R64-INV-STACK's hold frame: planned at 39, 40 qubits; CH1's MAIN at 38: 39 qubits).  This pass works
on the emitted op list and needs no change of the generator:

  1. LIVENESS BY VALUE on all 3,133,452 inputs (the kit's bit-sliced simulator, four outcome modes: zero, one, two
     random seeds): before and after every op, which wires are nonzero on some input.
  2. LIFETIME SHORTENING by commutation (same gates, same counts):
       - an X-measurement (HMR, unconditioned, outside condition blocks) is moved BACK to just after the last op
         that touches its wire or its outcome bit (R63-INV-GHOST's rule, here for every HMR);
       - an unconditioned X that turns a |0> wire into the constant 1 is moved FORWARD to just before the next op
         that touches the wire; an X that returns a constant-1 wire to |0> is moved BACK likewise.
  3. Liveness again; a wire is BUSY at an op when the op touches it or it is nonzero before and after the op.
     Maximal busy runs are intervals; the intervals of the wires above the 22 input wires are packed onto qubits
     22, 23, ... by the interval rule (lowest free qubit).  The file's qubit count becomes 22 + the largest number
     of such intervals alive at once.
  4. The op list is rewritten with the new qubit numbers (every other token unchanged) and written as .kmx
     (and .kmb when kickmix is importable).

  r65_inv_file38_compact.py IN.kmx OUT.kmx [N]        N > 0: liveness on a sample of N inputs (testing only!)
Check the result with tools/r63_inv_ghost_chk.py check IN.kmx OUT.kmx 0 (all inputs, four outcome modes: x restored,
no dirty wire, the phase equal to the base's) and with the official verifier.
Safe to run: reads IN.kmx, writes OUT.kmx (and OUT's .kmb)."""
import heapq
import os
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_lib as L          # noqa: E402
import invsqrt_gen_circ as C     # noqa: E402

MAX_X = 3133452
NIO = 22
MODES = (("zero", 0), ("one", 0), ("rand", 11), ("rand", 12))


def read(path):
    """lines of the file and, per line, None (not an op) or the op tuple (name, qubits, bit, cond, angle)."""
    lines = open(os.path.expanduser(path)).read().split("\n")
    recs = []
    for ln in lines:
        s = ln.split("#")[0].strip()
        if not s:
            recs.append(None)
            continue
        tok = s.split()
        name = tok[0]
        if name in ("REGISTER", "APPEND_TO_REGISTER", "INFO_TAG"):
            recs.append(None)
            continue
        cond = None
        if "if" in tok:
            i = tok.index("if")
            cond = int(tok[i + 1][1:])
            tok = tok[:i]
        qs, bit, angle = [], None, None
        for t in tok[1:]:
            if t[0] == "q":
                qs.append(int(t[1:]))
            elif t[0] == "b":
                bit = int(t[1:])
            else:
                angle = Fraction(t)
        recs.append((name, tuple(qs), bit, cond, angle))
    return lines, recs


def inputs(n):
    if n <= 0:
        return np.arange(1, MAX_X + 1, dtype=np.uint64)
    rng = np.random.default_rng(5)
    return np.unique(np.concatenate([rng.integers(1, MAX_X + 1, n), np.arange(1, 300),
                                     np.arange(MAX_X - 300, MAX_X + 1)])).astype(np.uint64)


def liveness(ops, xs, nq):
    """(before, after): per op the bitmask of wires nonzero on some input in some outcome mode."""
    sims = [C.Sim(xs, nq, hmr=m, seed=s) for m, s in MODES]
    mask = 0
    for w in range(nq):
        if any(s.q[w] for s in sims):
            mask |= 1 << w
    before, after = [], []
    for i, op in enumerate(ops):
        before.append(mask)
        for s in sims:
            s.run(ops, i, i + 1)
        for w in op[1]:
            if any(s.q[w] for s in sims):
                mask |= 1 << w
            else:
                mask &= ~(1 << w)
        after.append(mask)
    return before, after


def shorten(ops, before, after):
    """the commuting moves of step 2; returns the new op list and the counts of moves."""
    n = len(ops)
    blk = ("PUSH_CONDITION", "POP_CONDITION")
    depth, d = [], 0
    for op in ops:
        if op[0] == "POP_CONDITION":
            d -= 1
        depth.append(d)
        if op[0] == "PUSH_CONDITION":
            d += 1
    tag = [None] * n                                   # "hmr" | "birth" | "death"
    for i, op in enumerate(ops):
        name, qs, bit, cond, _ = op
        if depth[i] != 0 or cond is not None:
            continue
        if name == "HMR":
            tag[i] = "hmr"
        elif name == "X":
            w = qs[0]
            if w >= NIO and not (before[i] >> w) & 1:
                tag[i] = "birth"
            elif w >= NIO and not (after[i] >> w) & 1:
                tag[i] = "death"
    seq = list(range(n))
    moved = {"hmr": 0, "birth": 0, "death": 0}
    dist = {"hmr": 0, "birth": 0, "death": 0}
    # backward moves (hmr, death): scan forward so that an op lands behind ops already settled
    pos = 0
    while pos < n:
        i = seq[pos]
        kind = tag[i]
        if kind in ("hmr", "death"):
            name, qs, bit, cond, _ = ops[i]
            w = qs[0]
            j = pos
            while j > 0:
                p = ops[seq[j - 1]]
                if p[0] in blk or w in p[1]:
                    break
                if kind == "hmr" and (p[3] == bit or p[2] == bit):
                    break
                seq[j - 1], seq[j] = seq[j], seq[j - 1]
                j -= 1
            if j != pos:
                moved[kind] += 1
                dist[kind] += pos - j
        pos += 1
    # forward moves (birth): scan backward
    pos = n - 1
    while pos >= 0:
        i = seq[pos]
        if tag[i] == "birth":
            w = ops[i][1][0]
            j = pos
            while j < n - 1:
                p = ops[seq[j + 1]]
                if p[0] in blk or w in p[1]:
                    break
                seq[j + 1], seq[j] = seq[j], seq[j + 1]
                j += 1
            if j != pos:
                moved["birth"] += 1
                dist["birth"] += j - pos
        pos -= 1
    out = []                                           # peephole: X q ; X q (unconditioned, adjacent) is the identity
    moved["xx"] = 0
    for i in seq:
        op = ops[i]
        if out and op[0] == "X" and op[3] is None and out[-1] == op:
            out.pop()
            moved["xx"] += 1
            continue
        out.append(op)
    return out, moved, dist


def pack(ops, before, after, nq):
    """interval packing.  A wire is BUSY from the op that first touches it while it is |0> to the op after which it
    is |0> again.  The 22 input wires start and end on their own qubits (first and last interval); every other
    interval, theirs in between included, takes the lowest qubit that is free now and is not wanted by an input
    wire's last interval before this one ends.  Returns (new ops, qubits, the largest number of busy wires)."""
    n = len(ops)
    iv = [[] for _ in range(nq)]                       # per wire: [start op, end op] (inclusive; -1 = before op 0)
    opened = [None] * nq
    init = before[0] if n else 0
    for w in range(nq):
        if (init >> w) & 1 or w < NIO:
            opened[w] = -1
    busy, peak = sum(1 for w in range(nq) if opened[w] is not None), 0
    for i, op in enumerate(ops):
        a_ = after[i]
        for w in op[1]:
            if opened[w] is None:
                opened[w] = i
                busy += 1
        peak = max(peak, busy)
        for w in op[1]:
            if not (a_ >> w) & 1 and opened[w] is not None:
                iv[w].append([opened[w], i])
                opened[w] = None
                busy -= 1
    for w in range(nq):
        if opened[w] is not None:
            iv[w].append([opened[w], n])
        elif w < NIO:
            raise AssertionError("input wire %d is |0> at the end" % w)
    last_start = [iv[w][-1][0] for w in range(NIO)]    # the input wire's last interval must sit on qubit w
    items = sorted((s_, e_, w, k) for w in range(nq) for k, (s_, e_) in enumerate(iv[w]))
    until = {}                                         # qubit -> end op of the interval on it
    color = {}
    nphys = NIO
    for s_, e_, w, k in items:
        fixed = w < NIO and (k == 0 or k == len(iv[w]) - 1)
        if fixed:
            p_ = w
            assert until.get(p_, -2) < s_, ("the input wire's own qubit is taken", w, k, s_, until.get(p_))
        else:
            p_ = None
            for c_ in range(nphys):
                if until.get(c_, -2) >= s_:
                    continue
                if c_ < NIO and last_start[c_] <= e_ and until.get(c_, -2) < last_start[c_]:
                    continue                           # the input wire c_ comes back before this interval ends
                p_ = c_
                break
            if p_ is None:
                p_ = nphys
                nphys += 1
        until[p_] = e_
        color[(w, k)] = p_
    ptr = [0] * nq
    out = []
    for i, op in enumerate(ops):
        name, qs, bit, cond, angle = op
        ph = []
        for w in qs:
            while iv[w][ptr[w]][1] < i:
                ptr[w] += 1
            assert iv[w][ptr[w]][0] <= i <= iv[w][ptr[w]][1], ("op outside its wire's interval", i, w)
            ph.append(color[(w, ptr[w])])
        assert len(set(ph)) == len(ph), ("two wires of one op on one qubit", i, op, ph)
        out.append((name, tuple(ph), bit, cond, angle))
    return out, nphys, peak


def main():
    inp, outp = sys.argv[1], sys.argv[2]
    nsamp = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    t0 = time.time()
    lines, recs = read(inp)
    ops = [r for r in recs if r is not None]
    nq = 1 + max(q for op in ops for q in op[1])
    xs = inputs(nsamp)
    print("in: %s: %d ops, %d qubits; liveness on %d inputs, modes %s" % (inp, len(ops), nq, len(xs), MODES), flush=True)
    before, after = liveness(ops, xs, nq)
    pk0 = max(bin(a).count("1") for a in after)
    print("  by value: peak %d nonzero wires; op boundaries at the peak: %d   [%.0f s]" % (
        pk0, sum(1 for a in after if bin(a).count("1") == pk0), time.time() - t0), flush=True)
    ops2, moved, dist = shorten(ops, before, after)
    print("  moved: %d X-measurements back (%d op places), %d constant births forward (%d), %d constant deaths back (%d);"
          " adjacent X X pairs removed: %d" % (moved["hmr"], dist["hmr"], moved["birth"], dist["birth"], moved["death"],
                                               dist["death"], moved["xx"]), flush=True)
    before2, after2 = liveness(ops2, xs, nq)
    pk1 = max(bin(a).count("1") for a in after2)
    hist = {}
    for a in after2:
        k = bin(a).count("1")
        hist[k] = hist.get(k, 0) + 1
    print("  after the moves: peak %d nonzero wires; boundaries per count (top 4): %s   [%.0f s]" % (
        pk1, sorted(hist.items())[-4:], time.time() - t0), flush=True)
    at_pk = [i for i, a in enumerate(after2) if bin(a).count("1") == pk1]
    if len(at_pk) <= int(os.environ.get("R65_PKSHOW", "12")):   # the last moments at the peak: the ops around them
        for i in at_pk:
            for j in range(max(0, i - 3), min(len(ops2), i + 3)):
                print("    %s op %6d  live after %2d  %s" % ("->" if j == i else "  ", j, bin(after2[j]).count("1"),
                                                         L.op_text(ops2[j])), flush=True)
    ops3, q3, pkb = pack(ops2, before2, after2, nq)
    print("  packed: %d qubits (was %d); the largest number of busy wires at one op: %d" % (q3, nq, pkb), flush=True)
    it = iter(ops3)
    with open(os.path.expanduser(outp), "w") as fh:
        for ln, r in zip(lines, recs):
            if r is None:
                fh.write(ln + "\n")
            else:
                o_ = next(it, None)                    # (fewer ops than op lines after the peephole: the last slots stay empty)
                if o_ is not None:
                    fh.write(L.op_text(o_) + "\n")
    assert next(it, None) is None
    print("out: %s" % outp, flush=True)
    try:
        sys.path.insert(0, os.path.join(os.path.expanduser("~/qcc/data"), "kickmix"))
        import kickmix as km
        kc = km.Circuit.from_file(os.path.expanduser(outp), format="kmx")
        kb = os.path.splitext(os.path.expanduser(outp))[0] + ".kmb"
        kc.to_file(kb, format="kmb")
        print("out: %s (num_qubits %s)" % (kb, getattr(kc, "num_qubits", "?")), flush=True)
    except Exception as e:                             # noqa: BLE001
        print("kmb not written (%s: %s)" % (type(e).__name__, e), flush=True)


if __name__ == "__main__":
    main()
