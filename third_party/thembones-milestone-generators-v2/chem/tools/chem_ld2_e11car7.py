#!/usr/bin/env python3
"""chem_ld2_e11car7.py IN.npz OUT.npz [--brot 7] [--ref REF.npz | --slots 3,5,6,9,10]: builder C2 (eleventh holder).

Puts the five J-slot carriers (top bit XOR the bit below = the slot's low-address monomial, on every row of every
column) on a finished table of ANY word width by exact gauge moves, so that a table that arrives without carriers
(R48-CHEM-Q7C's 7-bit tables: "no codes, no carriers") gets the cache the loaders' engine expects. Same move as the
loader kit's `chem_loaders_gauge.flip` (and R8's chem_own2_carrier_restore.py at 8 bits), written for width brot:
  half turn on edge g = (x, c): word g -> word g XOR 2^(brot-1)  (theta -> theta + pi),
  negation (word -> -word mod 2^brot) of the edge that created x and of every earlier edge out of x,
which leaves every decoded row vector unchanged up to its sign (checked below at width brot), so error_1, error_2,
Lambda and E_gap do not change. A negated word whose lowest set bit is bit brot-2 toggles its own carrier: the
passes repeat until every carrier is exact. Slot monomials: from --ref (the monomial each slot follows in an
8-bit table with carriers, read with R8's rule), else --slots (default 3,5,6,9,10: l0l1, l0l2, l1l2, l0l3, l1l3).
OUT.npz = IN.npz with `angle_words` replaced (all other arrays copied as they are).
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_gauge as gauge  # noqa: E402  (loader kit, read-only import: tree() and decode())

PAIRS = (3, 5, 6, 9, 10, 12)


def carrier_bits(w, brot):
    return ((w >> (brot - 1)) ^ (w >> (brot - 2))) & 1


def ref_slots(path, brot=8):
    """The low monomial each of the 5 slots follows in a table with exact carriers (R8's slot_monomials rule)."""
    w = np.array(np.load(path)["angle_words"]).astype(np.int64)
    low = np.arange(w.shape[0]) % 16
    top = carrier_bits(w, brot)
    slot = {}
    for s in range(5):
        ms = set()
        for c in range(s, w.shape[1], 5):
            for m in PAIRS:
                if np.all(top[:, c] == ((low & m) == m)):
                    ms.add(m)
                    break
        assert len(ms) == 1, f"reference slot {s}: monomials {sorted(ms)}"
        slot[s] = ms.pop()
    return [slot[s] for s in range(5)]


def flip(words, edges, g, rows, brot):
    """The exact move for word g on the given rows (boolean mask), at word width brot."""
    edges, made_by, out = gauge.tree(edges)
    w = np.array(words, dtype=np.int64, copy=True)
    x = edges[g][0]
    w[rows, g] ^= 1 << (brot - 1)
    neg = [t for t in out[x] if t < g]
    if x in made_by:
        neg.append(made_by[x])
    mask = (1 << brot) - 1
    for t in neg:
        w[rows, t] = (-w[rows, t]) & mask
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--brot", type=int, default=7)
    ap.add_argument("--ref", default=None)
    ap.add_argument("--slots", default="3,5,6,9,10")
    a = ap.parse_args()
    with np.load(a.src, allow_pickle=False) as z:
        arrays = {k: np.array(z[k]) for k in z.files}
    w0 = arrays["angle_words"].astype(np.int64)
    edges = arrays["tree_edges"].astype(np.int64)
    assert w0.min() >= 0 and w0.max() < (1 << a.brot), (w0.min(), w0.max(), a.brot)
    slots = ref_slots(a.ref) if a.ref else [int(x) for x in a.slots.split(",")]
    low = np.arange(w0.shape[0]) % 16
    want = np.stack([((low & slots[c % 5]) == slots[c % 5]).astype(np.int64) for c in range(w0.shape[1])], 1)
    print(f"slots {slots}; rows {w0.shape[0]}, brot {a.brot}; exact columns before: "
          f"{int(np.sum(np.all(carrier_bits(w0, a.brot) == want, axis=0)))}/75", flush=True)
    w = w0.copy()
    for it in range(30):
        bad = carrier_bits(w, a.brot) != want
        cols = np.flatnonzero(bad.any(axis=0))
        if cols.size == 0:
            break
        print(f"  pass {it + 1}: {cols.size} columns, {int(bad.sum())} wrong entries", flush=True)
        for c in cols:                                   # columns in edge order: later moves negate earlier edges
            r = carrier_bits(w[:, c], a.brot) != want[:, c]
            if r.any():
                w = flip(w, edges, int(c), r, a.brot)
    a0, a1 = gauge.decode(w0, edges, a.brot), gauge.decode(w, edges, a.brot)
    dev = float(np.max(np.minimum(np.abs(a0 - a1).max(axis=1), np.abs(a0 + a1).max(axis=1))))
    ok = bool(np.all(carrier_bits(w, a.brot) == want))
    print(f"decoded vectors: largest deviation {dev:.3e} (up to each row's sign); carriers exact on all 75 columns: "
          f"{ok}; bit-0 changes {int(np.sum((w ^ w0) & 1))}; words changed {int(np.sum(w != w0))} of {w.size}")
    assert dev < 1e-9 and ok, "gauge check failed: nothing written"
    arrays["angle_words"] = w.astype(arrays["angle_words"].dtype)
    np.savez(a.dst, **arrays)
    print(f"wrote {a.dst}")


if __name__ == "__main__":
    main()
