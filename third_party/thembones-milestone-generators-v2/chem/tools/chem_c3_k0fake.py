#!/usr/bin/env python3
"""chem_c3_k0fake.py IN OUT: builder C3. LOADER 0 SEEN AS A LOADER k >= 1 (for the kept-product stage with 8 helpers).

Loader 0 has no carriers (its words start at 0), so its low cache holds 11 ANDs on 11 wires and only 3 wires are
free with the 68-AND cache. Its five carrier functions g_j (bit 6 XOR bit 5 of every word of slot j, the same in
every batch) are single quadratic monomials of the low address bits on the tables of today: computed by one Toffoli
each ONTO THE CLEAN BIT-6 WIRES (tools/chem_c3_emit68.py --k0car), they leave 8 free wires and the same 79 of infra.
This tool writes a fake in-folder where 'loader 1' is loader 0 behind a fake batch that holds only g_j in bit 6:
  OUT/pad.npz   = [fake batch | the padded table] (angle_words, 1,232 x 80)
  OUT/proj.json = IN/proj.json with V['1'] = V['0'] (loader 0's own C-4 block set)
  OUT/k1*.json  = IN/k0*.json (the basis files: the same low coordinates, checked below)
Then: tools/sh/chem_c3_scan.sh OUT 1 <dir> with EMITX=--k0car.
"""
import json
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402


def main():
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    core.BROT, core.J = int(os.environ.get("CHEM_BROT", "7")), 5
    J, B = 5, core.BROT
    w = np.asarray(np.load(os.path.join(src, "pad.npz"))["angle_words"], np.int64)
    A = w.shape[0]
    l = np.arange(A) & 15
    g = []
    for j in range(J):
        x = ((w[:, j::J] >> (B - 1)) ^ (w[:, j::J] >> (B - 2))) & 1          # rows x batches
        gj = x[:, 0]
        assert (x == gj[:, None]).all(), f"slot {j}: bit6 ^ bit5 differs between batches"
        for ll in range(16):
            assert np.all(gj[l == ll] == gj[l == ll][0]), f"slot {j}: carrier is not a function of the low bits"
        g.append(gj)
    fake = np.zeros((A, J), np.int64)
    for j in range(J):
        fake[:, j] = g[j] << (B - 1)
    w2 = np.concatenate([fake, w], axis=1)
    np.savez(os.path.join(out, "pad.npz"), angle_words=w2)
    pj = json.load(open(os.path.join(src, "proj.json")))
    pj["V"] = dict(pj["V"])
    pj["V"]["1"] = pj["V"]["0"]
    json.dump(pj, open(os.path.join(out, "proj.json"), "w"), indent=1)
    # the fake loader 1 must have the same low basis (coordinates) as the real loader 0
    a0 = core.analyze(w, 0, 4)
    a1 = core.analyze(w2, 1, 4)
    b0 = [b["anf"] for b in a0.basis]
    b1 = [b["anf"] for b in a1.basis]
    assert b0 == b1, f"basis differs: {b0} vs {b1}"
    assert len(a1.carriers) == J and len(a1.out_wires) == J * (B - 1)
    assert (a0.X0 == a1.X0).all(), "the mixed tensor differs"
    n = 0
    for f in os.listdir(src):
        if f.startswith("k0") and f.endswith(".json"):
            shutil.copy(os.path.join(src, f), os.path.join(out, "k1" + f[2:]))
            n += 1
    print(f"fake in-folder {out}: carriers {[[m for m in core.bits(c['anf'])] for c in a1.carriers]}, "
          f"basis equal, X0 equal; {n} basis files; V0 = {pj['V']['0']}")


if __name__ == "__main__":
    main()
