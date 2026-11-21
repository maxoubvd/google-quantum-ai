#!/usr/bin/env python3
"""r46_dqi_scan_heldchk.py -- R46-DQI-SCAN: can each idle single-AND wire be measured out and remade?

For every gap of tools/r46_dqi_scan_held.py whose wire t was made by one CCX a b -> t (and nothing else since its
reset), on sampled valid shots (P1's value model; half of the runs on the verifier's law, half with the error weight
uniform on 0..35 so that rare flags fire): is t == (wire a) AND (wire b) with the CURRENT values of the two birth
wires (1) just after the touch that opens the gap (then `HMR t; CZ a b if outcome` parks it) and (2) just before the
touch that closes it (then `CCX a b -> t` remakes it)? Gaps that pass both on every sample go to ok.npy.
Usage: python r46_dqi_scan_heldchk.py <dqi.kmb> <gaps.npz> <zero.npz> <ok.npy> [runs=48] [seed=5] [min gap=500]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
from dqi_pp_act import pack, v, km  # noqa: E402
from r46_dqi_scan_rank import simulate_snapop  # noqa: E402


def select(g, Z, n, nq, ming):
    zw = Z["wire"].astype(np.int64)
    s = np.clip(Z["start"].astype(np.int64), 0, n)
    e = np.clip(Z["end"].astype(np.int64), 0, n)
    ok = e > s
    zw, s, e = zw[ok], s[ok], e[ok]
    order = np.lexsort((s, zw))
    zw, s, e = zw[order], s[order], e[order]
    first = np.searchsorted(zw, np.arange(nq + 1))
    w, gs, ge, nwr, bop, bk = [g[:, k] for k in range(6)]
    keep = (ge - gs) >= ming
    mid = (gs + ge) // 2
    held = np.zeros(len(g), bool)
    for i in np.flatnonzero(keep):
        a, b = first[w[i]], first[w[i] + 1]
        k = np.searchsorted(s[a:b], mid[i], side="right") - 1
        held[i] = not (k >= 0 and e[a + k] > mid[i])
    if os.environ.get("R46_KIND", "and") == "copy":      # wires made by ONE CX a -> t: is t == wire a (a free park)?
        return keep & held & (nwr == 1) & (bk == 0x21)
    return keep & held & (nwr == 1) & (bk == 0x31)


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    g = np.load(sys.argv[2])["gaps"]
    Z = np.load(sys.argv[3])
    out = sys.argv[4]
    R = int(sys.argv[5]) if len(sys.argv) > 5 else 48
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 5
    ming = int(sys.argv[7]) if len(sys.argv) > 7 else 500
    ops = np.ascontiguousarray(np.asarray(c.ops))
    n = len(ops)
    nq = c.num_qubits
    sel = np.flatnonzero(select(g, Z, n, nq, ming))
    w, gs, ge = g[sel, 0], g[sel, 1], g[sel, 2]
    a, b = g[sel, 6], g[sel, 7]
    COPY = os.environ.get("R46_KIND", "and") == "copy"
    if COPY:
        b = a
    sites = np.unique(np.concatenate([gs + 1, ge]))
    ns = len(sites)
    print("single-AND gaps %d, sites %d, samples %d (%.0fs)" % (len(sel), ns, 64 * R, time.time() - t0), flush=True)
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    law = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    uni = np.full(spec.max_errors + 1, 1.0 / (spec.max_errors + 1))
    rng = np.random.default_rng(seed)
    i_open = np.searchsorted(sites, gs + 1)
    i_close = np.searchsorted(sites, ge)
    bad_open = np.zeros(len(sel), np.int64)
    bad_close = np.zeros(len(sel), np.int64)
    ones = np.zeros(len(sel), np.int64)
    snap = np.zeros((ns, nq), np.uint64)
    for r in range(R):
        vinit = np.zeros(nq, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        syn, mbuc, _ = v.sample_shots(spec, tables, rng, 64, law if r % 2 == 0 else uni)
        vinit[:840] = pack(v.to_bits(syn, spec.degree))
        mb = pack(v.to_bits(mbuc, spec.degree))
        binit[:len(mb)] = mb
        js = simulate_snapop(ops, *arr, vinit, binit, seed * 1000 + r, sites, snap)
        assert js == ns
        for which, idx, bad in ((0, i_open, bad_open), (1, i_close, bad_close)):
            t = snap[idx, w]
            x = snap[idx, a] if COPY else snap[idx, a] & snap[idx, b]
            diff = t ^ x
            bad += np.array([bin(int(q)).count("1") for q in diff])
            if which == 1:
                ones += np.array([bin(int(q)).count("1") for q in t])
        if r % 8 == 7 or r == R - 1:
            print("run %d: gaps failing at the open %d, at the close %d; gaps whose bit was ever 1: %d (%.0fs)" % (
                r, int((bad_open > 0).sum()), int((bad_close > 0).sum()), int((ones > 0).sum()), time.time() - t0),
                flush=True)
    good = (bad_open == 0) & (bad_close == 0)
    np.save(out, sel[good])
    np.savez(out + ".detail.npz", sel=sel, bad_open=bad_open, bad_close=bad_close, ones=ones, samples=64 * R)
    print("PASS (t == a AND b at the open and at the close on all %d samples): %d of %d gaps; of them with the bit "
          "seen at 1: %d; fail at the open %d, at the close only %d" % (
              64 * R, int(good.sum()), len(sel), int((good & (ones > 0)).sum()), int((bad_open > 0).sum()),
              int(((bad_open == 0) & (bad_close > 0)).sum())))


if __name__ == "__main__":
    main()
