#!/usr/bin/env python
"""chem_wr7_obtry (WR7): the 76 one-body rows of a finished table (its h1' held, loaded as r50_chem_lbl_obchk.py does)
rounded under a wire list N times from random initial row-address assignments (chem_wr7_table7c.one_try: the
designer's progressive_dyn + polish, wires asserted), in W processes: error_1 best / median / worst per wire list.

    python chem_wr7_obtry.py TABLE.npz --wires @a.txt [@b.txt ...] [--tries 16] [--workers 8] [--seed 0]
"""
import argparse
import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_onebody import (check_fcidump, load_fixed_shifts, load_reference,  # noqa: E402
                                  reconstruct_u_vectors, symmetric_from_lower)
from r50_chem_lbl_chain import parse  # noqa: E402
from chem_wr7_table7c import one_try  # noqa: E402


def h1_prime(table):
    z = np.load(table)
    words = np.asarray(z["angle_words"], np.int64)
    edges = np.asarray(z["tree_edges"], np.int64)
    w = np.asarray(z["W"], np.float64)
    Wo = np.asarray(z["H2_I"], np.float64) - w.sum(axis=1)
    R, B, C = w.shape
    N = edges.shape[0] + 1
    ref = load_reference(check_fcidump())
    sh = load_fixed_shifts()
    S = np.asarray(symmetric_from_lower(sh["Hsym_lower"], N), np.float64)
    beta = float(sh["H1_I"]) - 0.5 * float(np.trace(S))
    eta = int(ref["electrons"])
    tei = np.asarray(ref["tei"], np.float64)
    h1_thc = np.asarray(ref["oei"], np.float64) - 0.5 * np.einsum("prrq->pq", tei)
    del tei
    allv = reconstruct_u_vectors(words / 64.0, edges)
    allv = allv / np.linalg.norm(allv, axis=1, keepdims=True)
    u = allv[N:].reshape(R, B, N)
    k = np.einsum("rc,rbc->rb", Wo, w)
    H = h1_thc - 0.5 * eta * S + beta * np.eye(N) - np.einsum("rb,rbp,rbq->pq", k, u, u, optimize=True)
    return H, edges


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--wires", nargs="+", required=True)
    ap.add_argument("--tries", type=int, default=16)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--iters", type=int, default=300)
    ap.add_argument("--sweeps", type=int, default=40)
    ap.add_argument("--pairs", type=int, default=0, help="1: also the pair moves (chem_wr7_obpolish.polish_pairs)")
    ap.add_argument("--fine", type=int, default=1, help="each one-body stage cut into this many steps")
    a = ap.parse_args()
    t0 = time.time()
    H, edges = h1_prime(a.table)
    for wl in a.wires:
        wires = parse(wl)
        args = [(H, edges, 7, a.iters, a.sweeps, wires, t, a.seed, a.pairs, a.fine) for t in range(a.tries)]
        with mp.get_context("fork").Pool(max(1, min(a.workers, a.tries))) as pool:
            res = pool.map(one_try, args)
        es = np.array([r[0] for r in res]) * 1e3
        print(f"{os.path.basename(a.table)} {wl} ({len(wires)} wires): error_1 best {es.min():.3f} (try {int(es.argmin())}), "
              f"try 0 {es[0]:.3f}, median {np.median(es):.3f}, worst {es.max():.3f} mHa over {a.tries} placements "
              f"(pairs {a.pairs}, fine {a.fine}; {time.time() - t0:.0f} s); all " + " ".join(f"{x:.2f}" for x in es), flush=True)


if __name__ == "__main__":
    main()
