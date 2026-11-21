#!/usr/bin/env python
"""r48_chem_q7c_obsigned (design seat R48-CHEM-Q7C): the one-body rows of a finished table moved into the SIGNED
RANGE (every word in [-2^(brot-2), 2^(brot-2)] = theta in [-pi/2, pi/2]) by an exact gauge move.

    python r48_chem_q7c_obsigned.py TABLE.npz OUT.npz [--brot 7] [--rows 76]

Why [M box 2, 08:27Z to 08:40Z Nov 16, logs/lpx.log, lpy.log, lpz.log, lpw.log]: the builders' loader compiler prices a
7-bit (9,127,C) table at 4,045 Toffolis per pass when the two top bits of every word agree (the rounding tool leaves
the two-body rows so: "every word in the signed range") and at 4,419 to 4,485 when the 76 one-body rows use the full
circle (r48_chem_q7c_table7.py did: its rows came from a full-circle rounding). Only the two top bits matter: my
rows' bits 0..4 under the reference rows' bits 5..6 price 4,044 to 4,045.

The move is exact on the words [P]: a row's vector is decoded from its words (the verifier's Givens convention), then
re-encoded with the sign of every rotation's parent amplitude pulled up the tree (chem_own2_round.angles_and_mass,
the signed gauge): every |cos| and |sin| is unchanged, so every new angle is again a grid point, and the decoded
vector is the old one up to ONE global sign per row (the table uses u u^T only). The tool checks both statements
numerically for every row and refuses to write otherwise. A word of +pi/2 exactly (parent amplitude 0) is written as
-pi/2 with the child's sub-vector negated by the same re-encoding where that is exact; the count of words left at
+2^(brot-2) is printed.
"""
import sys

import numpy as np


def decode_pi(th, edges):
    th = np.pi * np.asarray(th, np.float64)
    v = np.zeros((th.shape[0], edges.shape[0] + 1))
    v[:, 0] = 1.0
    for t, (p, c) in enumerate(edges):
        a, b = v[:, p].copy(), v[:, c].copy()
        v[:, p] = np.cos(th[:, t]) * a - np.sin(th[:, t]) * b
        v[:, c] = np.sin(th[:, t]) * a + np.cos(th[:, t]) * b
    return v


def signed_angles(vectors, edges):
    """theta / pi in [-1/2, 1/2] and the squared amplitude entering each rotation (the signed gauge)."""
    a = np.array(vectors, dtype=np.float64)
    theta = np.zeros((len(a), len(edges)))
    mass = np.zeros((len(a), len(edges)))
    for t in range(len(edges) - 1, -1, -1):
        p, q = edges[t]
        s = np.where(a[:, p] > 0, 1.0, np.where(a[:, p] < 0, -1.0, np.where(a[:, q] > 0, -1.0, 1.0)))
        theta[:, t] = np.arctan2(s * a[:, q], s * a[:, p])
        h = np.hypot(a[:, p], a[:, q])
        mass[:, t] = h * h
        a[:, p] = s * h
        a[:, q] = 0.0
    return theta / np.pi, mass


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    brot, nrows = 7, 76
    for i, x in enumerate(sys.argv):
        if x == "--brot":
            brot = int(sys.argv[i + 1]); args.remove(sys.argv[i + 1])
        if x == "--rows":
            nrows = int(sys.argv[i + 1]); args.remove(sys.argv[i + 1])
    src, out = args
    z = dict(np.load(src))
    words = np.asarray(z["angle_words"], np.int64)
    edges = np.asarray(z["tree_edges"], np.int64)
    half, full, q = 1 << (brot - 1), 1 << brot, 1 << (brot - 2)
    w0 = words[:nrows] % full
    v0 = decode_pi(w0 / half, edges)
    # exact zeros: cos / sin of grid multiples of pi/2 leave 1e-16 dust that would decide signs at random
    v0c = np.where(np.abs(v0) < 1e-12, 0.0, v0)
    th, mass = signed_angles(v0c, edges)
    w1 = np.rint(th * half).astype(np.int64)
    dead = mass < 1e-20                                   # no amplitude enters: the word is free; keep it small
    fold = w0 % half
    fold = np.where(fold > q, half - fold, fold)          # |angle| folded to [0, pi/2]
    w1 = np.where(dead, np.where(fold == q, -q, fold), w1)
    off = np.abs(th * half - w1)[~dead].max()
    v1 = decode_pi(w1 / half, edges)
    sgn = np.sign((v0 * v1).sum(1))
    err = np.abs(v1 * sgn[:, None] - v0).max()
    print("rows %d: max |angle - grid| %.2e words; max |decoded new - (+-)decoded old| %.2e; rows with the global "
          "sign flipped %d; words with no amplitude entering %d" % (nrows, off, err, int((sgn < 0).sum()), int(dead.sum())))
    if off > 1e-6 or err > 1e-9 or np.abs(w1).max() > q:
        print("REFUSED: the move is not exact on this table"); sys.exit(1)
    wu = w1 % full
    top = (wu >> (brot - 1)) & 1
    nxt = (wu >> (brot - 2)) & 1
    print("before: words outside the signed range %.4f; after: %.4f (words at +%d: %d, at -%d: %d); top two bits differ "
          "on %d of %d words; two-body rows: top two bits differ on %d of %d words" % (
              ((w0 > q) & (w0 < full - q)).mean(), ((wu > q) & (wu < full - q)).mean(), q, int((w1 == q).sum()), q,
              int((w1 == -q).sum()), int((top != nxt).sum()), wu.size,
              int(((((words[nrows:] % full) >> (brot - 1)) & 1) != (((words[nrows:] % full) >> (brot - 2)) & 1)).sum()),
              words[nrows:].size))
    new = words.copy()
    new[:nrows] = wu
    z["angle_words"] = new.astype(z["angle_words"].dtype)
    np.savez(out, **z)
    print("written", out)


if __name__ == "__main__":
    main()
