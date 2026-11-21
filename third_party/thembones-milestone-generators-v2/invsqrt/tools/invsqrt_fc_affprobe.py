#!/usr/bin/env python3
"""invsqrt_fc_affprobe.py -- FC (fourth holder): at the moment of every adder call in the chosen stages, are some LIVE
wires affine functions (XOR of other live wires + constant) of the others? Such a wire can be cleared and restored by
CNOTs only: a free extra clean wire for that block (and its explicit undo). Sample-based GF(2) rank test; candidates
must be re-checked on all inputs.
  invsqrt_fc_affprobe.py SRC CFG.json NSAMPLE [STAGES]      (env knobs as for the build; STAGES e.g. MAIN,H1)
"""
import json
import os
import sys

import numpy as np


def gf2_dependent(M):
    """columns of M (uint8 0/1, samples x cols) that are XOR-combinations of other columns and the all-ones column:
    returns list of (col, [cols it depends on])."""
    n, k = M.shape
    A = np.concatenate([M, np.ones((n, 1), dtype=np.uint8)], axis=1).T.copy()   # rows = columns of M (+ const)
    # reduce rows over GF(2) keeping track of combinations
    rows = [int("".join(map(str, r[::-1])), 2) if False else None for r in A]
    vecs = [int.from_bytes(np.packbits(r, bitorder="little").tobytes(), "little") for r in A]
    combos = [1 << i for i in range(len(vecs))]
    basis = []                                     # (pivot bit, vec, combo)
    deps = []
    for i, (v, cmb) in enumerate(zip(vecs, combos)):
        for pb, bv, bc in basis:
            if (v >> pb) & 1:
                v ^= bv
                cmb ^= bc
        if v == 0:
            if i < k:
                deps.append((i, [j for j in range(len(vecs)) if (cmb >> j) & 1 and j != i]))
        else:
            pb = v.bit_length() - 1
            basis.append((pb, v, cmb))
    return deps


def main():
    src = os.path.abspath(os.path.expanduser(sys.argv[1]))
    sys.path.insert(0, src)
    import invsqrt_fw_fcpolrun as R
    F8, K, w = R.load(src)
    A = sys.modules["invsqrt_fc_arith"]
    stages = (sys.argv[4] if len(sys.argv) > 4 else "MAIN,H1").split(",")
    snaps = []
    orig = A.add

    def add(c, src, tgt, *a, **kw):
        tag = c.tags[-1][1] if getattr(c, "tags", None) else "?"
        if tag in stages:
            ops_ = sorted({int(q) for q in list(src) + list(tgt) + [kw.get("cin")] + list(a[:1])
                           if isinstance(q, (int, np.integer)) and int(q) >= 0})
            snaps.append((tag, len(c.ops), sorted(set(range(c.next)) - set(c.free)), c.cap - c.live, ops_))
        return orig(c, src, tgt, *a, **kw)
    A.add = add
    NW = sys.modules["invsqrt_fc_narrow"]               # also every narrow block (explicit-undo blocks)
    oblk = NW.block

    def _block(c, s_, t_, ci, *a, **kw):
        tag = c.tags[-1][1] if getattr(c, "tags", None) else "?"
        if tag in stages:
            ops_ = sorted({int(q) for q in list(s_) + list(t_) + [ci]
                           if isinstance(q, (int, np.integer)) and int(q) >= 0})
            snaps.append((tag + "/blk", len(c.ops), sorted(set(range(c.next)) - set(c.free)), c.cap - c.live, ops_))
        return oblk(c, s_, t_, ci, *a, **kw)
    for name in list(sys.modules):
        mod = sys.modules[name]
        if mod is not None and getattr(mod, "block", None) is oblk:
            mod.block = _block
    cfg = R.cfg_for(F8, src, json.load(open(sys.argv[2])))
    c, regs = F8.build_full(dict(cfg))
    nf = regs["nfwd"]
    rng = np.random.default_rng(3)
    x = np.unique(np.concatenate([rng.integers(1, F8.MAX_X + 1, int(sys.argv[3])), np.arange(1, 64),
                                  np.arange(F8.MAX_X - 64, F8.MAX_X + 1)])).astype(np.int64)
    depth, dep = 0, []
    for op in c.ops:
        dep.append(depth)
        if op[0] == "BBEG":
            depth += 1
        elif op[0] == "BEND":
            depth -= 1
    found = 0
    for tag, oi, live, free, opw in snaps:
        if oi > nf or dep[oi] != 0 or free >= int(os.environ.get("AFF_MAXFREE", "5")):
            continue
        sim = K.Sim(len(x), c.next + 2)
        sim.fire = []
        sim.set_word(regs["X"], x)
        sim.run(c.ops, 0, oi)
        M = np.array([sim.q[wq] for wq in live], dtype=np.uint8).T     # samples x live wires
        const = [live[j] for j in range(len(live)) if M[:, j].min() == M[:, j].max()]
        deps = gf2_dependent(M)
        if deps or const:
            found += 1
            dd = [(live[i], [live[j] if j < len(live) else "1" for j in d][:8]) for i, d in deps[:4]]
            nonop = [i for i, d in deps if live[i] not in opw]
            print("%-5s op %6d live %2d free %d: constant wires %s; affine-dependent %d (not operands: %s): %s" % (
                tag, oi, len(live), free, const[:6], len(deps), [live[i] for i in nonop][:6], dd), flush=True)
    print("snapshots %d, with a free candidate %d" % (len(snaps), found))


if __name__ == "__main__":
    main()
