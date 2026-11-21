#!/usr/bin/env python3
"""chem_c4_wb.py: builder C4. The WIRE BASIS for a loader's generic outputs: basis output j := the mixed table of the
j-th non-free wire (wire q's table = XOR of X0 columns over G[q], chem_loaders_core.analyze), so G becomes the identity
on those wires. Any invertible output basis is legal (fan-outs absorb it; emit_c4 maps basis masks to wires through
an.G and checks every wire on every address); in this one every group of outputs closes exactly its own wires under
the C-4 closure, whatever the order of the groups.

  wire_basis(an)  -> modifies an.X0, an.G, an.out_level in place; returns the list of wires in basis order
  python -I chem_c4_wb.py <table.npz> --k 2 --dump d.json --out d_wb.json   (rewrites a chem_c4_dump file's X0)
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402


def wire_basis(an):
    nb = an.X0.shape[1]
    allw = [q for q, g in enumerate(an.G) if g]
    assert all(an.out_level[o] > 2 for o in range(nb)), "coded outputs: the wire basis is for wired loaders"
    # holder 2: a relation-bearing loader has more non-free wires than generic outputs (a wire = a combination of
    # others): the basis is a maximal independent set of wires (in wire order); the others are mapped onto it
    piv, tag, wires = {}, {}, []
    for q in allw:
        v, t = an.G[q], 0
        while v:
            h = v.bit_length() - 1
            if h not in piv:
                break
            v ^= piv[h]
            t ^= tag[h]
        if v:
            piv[v.bit_length() - 1] = v
            tag[v.bit_length() - 1] = t ^ (1 << len(wires))
            wires.append(q)
    assert len(wires) == nb, (len(wires), nb)
    X = np.zeros_like(an.X0)
    for j, q in enumerate(wires):
        for i in core.bits(an.G[q]):
            X[:, j] ^= an.X0[:, i]
    G = [0] * len(an.G)
    for j, q in enumerate(wires):
        G[q] = 1 << j
    for q in allw:
        if q in wires:
            continue
        v, t = an.G[q], 0                          # G[q] = XOR of the basis wires' G over the bits of t
        while v:
            h = v.bit_length() - 1
            v ^= piv[h]
            t ^= tag[h]
        G[q] = t
    # sanity: the new basis spans the same space (G was invertible on these wires)
    M = np.array([[(an.G[q] >> i) & 1 for i in range(nb)] for q in wires], np.uint8)
    rk = 0
    Mr = M.copy()
    for c in range(nb):
        piv = next((r for r in range(rk, nb) if Mr[r, c]), None)
        if piv is None:
            continue
        Mr[[rk, piv]] = Mr[[piv, rk]]
        for r in range(nb):
            if r != rk and Mr[r, c]:
                Mr[r] ^= Mr[rk]
        rk += 1
    assert rk == nb, "the wire map is not invertible on the generic outputs"
    an.X0 = X
    an.G = G
    return wires


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--dump", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    import chem_ld2_e19grp as eg
    words = np.asarray(np.load(a.table)["angle_words"], np.int64)
    V, words = eg.block_set(words, a.k)
    an = core.analyze(words, a.k, 4)
    wires = wire_basis(an)
    d = json.load(open(a.dump))
    d["X0"] = [[int(v) for v in row] for row in an.X0]
    d["wire_basis"] = wires
    json.dump(d, open(a.out, "w"))
    print(f"k {a.k}: {len(wires)} wires as the basis -> {a.out}")


if __name__ == "__main__":
    main()
