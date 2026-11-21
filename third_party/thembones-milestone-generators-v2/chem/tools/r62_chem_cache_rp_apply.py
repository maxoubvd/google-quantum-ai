"""r62_chem_cache_rp_apply: write a dead-word relief plan's TABLE CHANGES into the words and check them
(seat R62-CHEM-CACHE, re-pricer).

Input: a padded table (C2's pad.npz: 1,232 rows, the first --nreal are run) and the --chg json written by
tools/r62_chem_cache_rp_mask.py --dead (k, D, changes = [row, word slot j, bit b] of loader k's delta).
Each change flips bit b of word (k, j) of that row when that word is an exact don't-care, else of word (k-1, j)
(which must then be one); a flip of bit 5 also flips bit 6, so that bit 6 XOR bit 5 (the carrier) is unchanged.
Checks, all printed:
  1. every flipped word is dead (dead_mask of tools/r62_chem_cache_dead.py);
  2. the verifier's own evaluate_words on the run rows: error_1, error_2, Lambda, E_gap, lambda_eff before / after;
  3. loader k's analysis (chem_loaders_core.analyze): same output wires and the same output basis G;
     its delta differs exactly at the planned (row, wire) pairs;
  4. the neighbours that share a changed word (k - 1, k + 1): the plain engine's count before / after (a proxy
     for their own engines: same count = the change is generic for them).
Writes --out (the padded table with the changes) for tools/r58_chem_kept_rp_emit68.py.

  python -I tools/r62_chem_cache_rp_apply.py PAD.npz CHG.json --out PAD_mod.npz [--nreal 1219] [--noverify]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
from r62_chem_cache_dead import dead_mask  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pad")
    ap.add_argument("chg")
    ap.add_argument("--out", required=True)
    ap.add_argument("--nreal", type=int, default=1219)
    ap.add_argument("--brot", type=int, default=7)
    ap.add_argument("--noverify", action="store_true")
    ap.add_argument("--noneighbours", action="store_true")
    a = ap.parse_args()
    z = np.load(a.pad)
    d = {key: np.array(z[key]) for key in z.files}
    w = d["angle_words"].astype(np.int64)
    edges = d["tree_edges"].astype(np.int64)
    c = json.load(open(a.chg))
    k = int(c["k"])
    chg = [tuple(int(x) for x in t) for t in c["changes"]]
    dm = dead_mask(w[:a.nreal], edges, a.brot, 0)
    w2 = w.copy()
    flips = {}                                       # (row, column) -> xor mask
    for row, j, b in chg:
        own, prev = 5 * k + j, 5 * (k - 1) + j
        if dm[row, own]:
            col = own
        elif k >= 1 and dm[row, prev]:
            col = prev
        else:
            raise SystemExit(f"change at row {row} word slot {j} bit {b}: neither word ({k},{j}) nor ({k-1},{j}) is dead")
        m = 1 << b
        if b == a.brot - 2:
            m |= 1 << (a.brot - 1)                   # keep bit 6 XOR bit 5 (the carrier)
        flips[(row, col)] = flips.get((row, col), 0) ^ m
    for (row, col), m in flips.items():
        w2[row, col] ^= m
    nw = sum(1 for m in flips.values() if m)
    print(f"loader k{k}: {len(chg)} delta bits planned; {nw} words changed (all dead: checked); columns "
          f"{sorted(set(col for (_, col), m in flips.items() if m))}", flush=True)
    # 3. loader k's analysis
    core.J, core.BROT = 5, a.brot
    an0 = core.analyze(w, k, 4)
    an1 = core.analyze(w2, k, 4)
    sameG = an0.out_wires == an1.out_wires and list(an0.G) == list(an1.G)
    want = {}
    for row, j, b in chg:
        q = an0.out_wires.index((j, b))
        want[(row, q)] = want.get((row, q), 0) ^ 1
    diff = set()
    for q in range(len(an0.out_wires)):
        x = an0.delta_int[q] ^ an1.delta_int[q]
        r = 0
        while x:
            if x & 1:
                diff.add((r, q))
            x >>= 1
            r += 1
    planned = {key for key, v in want.items() if v}
    print(f"loader k{k}: same output wires and basis G: {sameG}; delta differs at {len(diff)} (row, wire) pairs; "
          f"equal to the plan's {len(planned)}: {diff == planned}; carriers before / after "
          f"{[cc['anf'] for cc in an0.carriers]} / {[cc['anf'] for cc in an1.carriers]}", flush=True)
    # 4. neighbours
    if not a.noneighbours:
        cols = sorted(set(col for (_, col), m in flips.items() if m))
        nbrs = sorted({k + 1 for col in cols if col // 5 == k and k + 1 <= 14} |
                      {k - 1 for col in cols if col // 5 == k - 1})
        for kk in nbrs:
            r0 = core.compile_circuit((w[:a.nreal], kk, 4, 0, 3, False, "auto"))
            r1 = core.compile_circuit((w2[:a.nreal], kk, 4, 0, 3, False, "auto"))
            print(f"neighbour k{kk} (plain engine, a proxy for its own): {r0['toffoli']} before, {r1['toffoli']} after",
                  flush=True)
    # 2. the verifier's own evaluation on the run rows
    if not a.noverify:
        data = os.path.expanduser("~/qcc/data")
        sys.path[:0] = [os.path.join(data, "kickmix"), os.path.join(data, "verifiers")]
        os.environ.setdefault("DFTHC_DATA", data)
        from dfthcsos import verify as V
        inputs = V.load_inputs({key: d[key] for key in ("U", "W", "H2_I", "H1_I", "Hsym_lower")})
        m0 = V.evaluate_words(inputs, w[:a.nreal], edges, a.brot)
        m1 = V.evaluate_words(inputs, w2[:a.nreal], edges, a.brot)
        for key in ("error_1", "error_2", "Lambda", "E_gap", "lambda_eff"):
            print(f"verifier {key}: {float(m0[key])!r} -> {float(m1[key])!r} (diff {float(m1[key]) - float(m0[key]):.3e})",
                  flush=True)
    d["angle_words"] = w2.astype(z["angle_words"].dtype)
    np.savez(a.out, **d)
    print(f"written {a.out}", flush=True)


if __name__ == "__main__":
    main()
