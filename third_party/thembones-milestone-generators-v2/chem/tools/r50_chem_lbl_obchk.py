"""r50_chem_lbl_obchk: THE ONE-BODY ROWS UNDER WIRE CHAINS, ALONE (R50-CHEM-LBL, construction 9b): error_1 of the 76
one-body rows of a table (its h1' held) when they are rounded by FREE's progressive rounding under a wire list of any
shape (chains allowed), the address of a row decided at freeze time (r50_chem_lbl_obchain.progressive_dyn), then
polished in steps of 2 on the wires' columns. Loading is r50_chem_lbl_obwires.py's, unchanged.

    python r50_chem_lbl_obchk.py TABLE.npz --wires "<list>"|@file [--base]
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_onebody import (angles_and_mass, check_fcidump, decode_pi, error1,  # noqa: E402
                                  load_fixed_shifts, load_reference, reconstruct_u_vectors, symmetric_from_lower)
from r49_chem_free_obround import polish_coded  # noqa: E402
from r50_chem_lbl_chain import chains_of, lin, parse  # noqa: E402
from r50_chem_lbl_obchain import progressive_dyn  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--wires", required=True)
    ap.add_argument("--iters", type=int, default=300)
    ap.add_argument("--sweeps", type=int, default=40)
    ap.add_argument("--base", action="store_true", help="also the same procedure without wires")
    a = ap.parse_args()
    t0 = time.time()
    z = np.load(a.table)
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
    e_own, perm, _, lam = error1(allv[:N], H)
    print(f"table {a.table}: its own one-body rows error_1 {1e3 * e_own:.4f} mHa ({time.time() - t0:.0f} s)", flush=True)
    lam_e, vec = np.linalg.eigh(H)
    th0, mass0 = angles_and_mass(vec[:, perm].T, edges)
    lam_at = lam_e[perm]
    wires = parse(a.wires)
    cols = set(chains_of(wires)[0])
    for name, wl in ((("no wires", []),) if a.base else ()) + ((f"{len(wires)} wires", wires),):
        wq, dyn = progressive_dyn(th0, edges, lam_at, H, 64.0, a.iters, wl)
        e_p = error1(decode_pi(wq / 64.0, edges), H)[0]
        wq3, e_d = polish_coded(wq, edges, H, 7, a.sweeps, cols if wl else set())
        addr_of = dyn.mrow.copy()
        bad = sum(int((((wq3[:, 5 * k_ + j_] ^ wq3[:, 5 * k_ - 5 + j_]) & 1) != lin(addr_of, m_)).sum()) for (k_, j_, _, m_) in wl)
        print(f"  {name}, addresses decided at freeze time: progressive error_1 {1e3 * e_p:.4f} mHa, after the polish {1e3 * e_d:.4f} mHa; "
              f"words left free {dyn.free_ok}, forced {dyn.forced_n}; relations violated {bad} ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
