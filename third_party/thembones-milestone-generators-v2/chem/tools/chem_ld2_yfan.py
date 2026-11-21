"""chem_ld2_yfan: for every nonzero low vector y (15 bits over the low monomials), the dimension of the fan-out
space S(y) = {s : y (x) s in F} of a dissected plan (F = span of its features).

Usage: python -I chem_ld2_yfan.py <dissect.json>
In a plan whose relations are block-diagonal in some basis V_1 + ... + V_5 of the low space, the y with a large
S(y) are the 35 nonzero vectors of the five blocks; mixed y have small S(y).
"""
import json
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from chem_ld2_rel import echelon, nullspace, NOUT  # noqa: E402


def main():
    d = json.load(open(sys.argv[1]))
    prods, fan = d["prods"], d["fan"]
    no = NOUT

    def feat(y, s):
        v = 0
        for m in range(1, 16):
            if (y >> m) & 1:
                v |= s << (no * (m - 1))
        return v
    F = []
    ys = Counter()
    for p, s in zip(prods, fan):
        if p["l2"] and p["h1"]:
            F.append(feat(p["l2"], s))
            ys[p["l2"]] += 1
        if p["l1"] and p["h2"]:
            F.append(feat(p["l1"], s))
            ys[p["l1"]] += 1
    piv = echelon(F)
    R = nullspace(list(piv.values()), 15 * no)
    # R as an array (nR, 15, 40)
    nR = len(R)
    Ra = np.zeros((nR, 15, no), np.uint8)
    for i, r in enumerate(R):
        for m in range(15):
            row = (r >> (no * m)) & ((1 << no) - 1)
            for t in range(no):
                Ra[i, m, t] = (row >> t) & 1
    dims = {}
    for y in range(1, 1 << 15):
        yb = np.array([(y >> m) & 1 for m in range(15)], np.uint8)
        M = np.tensordot(Ra, yb, axes=([1], [0])) % 2  # (nR, 40)
        vecs = [int.from_bytes(np.packbits(row, bitorder="little").tobytes(), "little") for row in M]
        rk = len(echelon(vecs))
        dims[y << 1] = no - rk  # y over monomials 1..15 -> bit m (match dissect masks: bit m = monomial m)
    hist = Counter(dims.values())
    print("dim S(y) histogram over all 32767 y:", sorted(hist.items()))
    big = sorted([(v, y) for y, v in dims.items() if v >= 8], reverse=True)
    print("y with dim S(y) >= 8:", len(big))
    print([(format(y >> 1, "015b"), v) for v, y in big[:60]])
    used = sorted(ys.items(), key=lambda kv: -kv[1])
    print("used y (count, dim S):", [(format(y >> 1, '015b'), c, dims[y]) for y, c in used[:40]])
    print("dim S over used y:", sorted(Counter(dims[y] for y in ys).items()))
    json.dump({str(y): v for y, v in dims.items()}, open(sys.argv[1] + ".yfan.json", "w"))


if __name__ == "__main__":
    main()
