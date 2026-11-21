#!/usr/bin/env python
"""r50_chem_lbl_table7w: R50-CHEM-LBL copy of r48_chem_q7c_table7 with ONE change: --wires "k:j:0:m,..." (delta wires
on address bits, construction 9): the 76 one-body rows are rounded by r50_chem_lbl_obwires.progressive_dyn (the
relation bit 0 of word 5k+j XOR bit 0 of word 5k-5+j = parity(address AND m) kept, the address of every row decided
at freeze time) and polished in steps of 2 on the wires' columns; they are written in ADDRESS order. Without --wires
the tool is the original. The two-body rows must already satisfy the wires (r50_chem_lbl_wtable.py).
Original header: r48_chem_q7c_table7 (design seat R48-CHEM-Q7C): the finishing chain of a 7-BIT rounding with no codes and no
carriers: rounded.npz of r48_chem_q7c_round.py -> a complete table.npz (the submission.npz layout: angle_words with
the 76 one-body rows first, tree_edges, U, W, H2_I, H1_I, Hsym_lower) that the hosts' verifier reads.

The steps are those of the builders' tools/chem_own2_table.py (offsets of least E_gap, provisional rows, free refit
without the error_1 limit, then `rounds` x (one-body rows for the current h1', refit with the rows fixed under both
error limits)). The differences:
  - the word width is a parameter (--brot, default 7), passed to the builders' refit problem (the verifier's
    evaluate_words at that width); r48_chem_q7c_refit.py is an unmodified copy of tools/chem_own_refit.py;
  - the 76 one-body rows come from r48_chem_q7c_onebody.py's staged rounding + exact +-1 polish on the table's own
    tree (no codes, no carriers, the full circle), not from chem_own_onebody.build_rows (8 bits, codes, carriers);
  - no code or carrier check on the two-body rows (R7g's recipe enforces none).
Every reported metric is the hosts' own evaluation (Problem.official).

    python r48_chem_q7c_table7.py ROUNDED.npz OUT_DIR [--brot 7] [--rounds 2] [--t2 19.9985] [--t1 19.99]
                                  [--perturb 1e-3] [--quick]
OUT_DIR/table.npz, OUT_DIR/table.json (the hosts' metrics after every step).
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r48_chem_q7c_onebody as OB  # noqa: E402
import r48_chem_q7c_refit as cor  # noqa: E402
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors  # noqa: E402


def log(*x):
    print(*x, flush=True)


def offsets_min_gap(prob, w, iters=600):
    """Offsets of the tightest lower bound at fixed weights (convex; copied from chem_own2_table.py)."""
    from scipy.optimize import minimize
    R, B, C = w.shape

    def f(x):
        (_, _, gap, _), (_, _, dgap, _) = prob.evaluate(w, x.reshape(R, C))
        return gap, dgap[1].ravel()

    res = minimize(f, np.zeros(R * C), jac=True, method="L-BFGS-B",
                   options={"maxiter": iters, "maxcor": 30, "ftol": 1e-15, "gtol": 1e-9})
    return res.x.reshape(R, C), float(res.fun), int(res.nit)


def h1_prime(prob, w, W):
    k = np.einsum("rc,rbc->rb", W, w)
    return prob.h1_base - np.einsum("rb,rbp,rbq->pq", k, prob.u, prob.u, optimize=True)


def set_rows(prob, rows):
    N = prob.N
    prob.words[:N] = rows
    v = reconstruct_u_vectors(prob.words[:N] / (1 << (prob.brot - 1)), prob.edges)
    prob.Vrows = v / np.linalg.norm(v, axis=1, keepdims=True)
    prob.perm = None


WIRES = []


def build_rows7(H, edges, brot, iters, sweeps):
    if WIRES and iters > 0:
        return build_rows7w(H, edges, brot, iters, sweeps)
    return build_rows7_plain(H, edges, brot, iters, sweeps)


def build_rows7w(H, edges, brot, iters, sweeps):
    """r50: the rows under the wires; returns the words in ADDRESS order."""
    import r50_chem_lbl_obwires as OW
    from r49_chem_free_obround import polish_coded
    assert brot == 7
    scale = float(1 << (brot - 1))
    lam_e, vec = np.linalg.eigh(H)
    th0, _ = OB.angles_and_mass(vec.T, edges)
    wq, dyn = OW.progressive_dyn(th0, edges, lam_e, H, scale, iters, WIRES)
    cols = set()
    for (k_, j_, _, _) in WIRES:
        cols |= {5 * k_ - 5 + j_, 5 * k_ + j_}
    wq, _ = polish_coded(wq, edges, H, brot, max(sweeps, 1), cols)
    wq = np.mod(np.array(wq, np.int64), 1 << brot)
    out = np.zeros_like(wq)
    out[dyn.mrow] = wq                                   # row -> its address
    for (k_, j_, _, m_) in WIRES:
        want = np.array([bin(a_ & m_).count("1") & 1 for a_ in range(len(out))])
        assert np.array_equal((out[:, 5 * k_ + j_] ^ out[:, 5 * k_ - 5 + j_]) & 1, want), "a wire is broken on the one-body rows"
    e = OB.error1(reconstruct_u_vectors(out / scale, edges), H)[0]
    log(f"[table7w] one-body rows under {len(WIRES)} wires: second words left free {dyn.free_ok}, forced {dyn.forced_n}")
    return out, e


def build_rows7_plain(H, edges, brot, iters, sweeps):
    """76 one-body rows at `brot` bits for h1' = H: eigenvectors -> angles on the tree -> staged rounding with a
    refit of the free angles on error_1^2 -> exact +-1 polish. Returns (words in [0, 2^brot), error_1 in Ha)."""
    scale = float(1 << (brot - 1))
    lam_e, vec = np.linalg.eigh(H)
    th0, _ = OB.angles_and_mass(vec.T, edges)
    if iters <= 0:
        wq = np.round(th0 * scale).astype(np.int64)
    else:
        wq = OB.progressive(th0, edges, lam_e, H, scale, iters)
    wq = np.mod(wq, 1 << brot)
    if sweeps > 0:
        wq, _ = OB.polish(wq, edges, H, brot, sweeps)
    e = OB.error1(reconstruct_u_vectors(np.mod(wq, 1 << brot) / scale, edges), H)[0]
    return np.mod(wq, 1 << brot), e


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rounded")
    ap.add_argument("out")
    ap.add_argument("--brot", type=int, default=7)
    ap.add_argument("--rounds", type=int, default=2)
    ap.add_argument("--t2", type=float, default=19.9985)
    ap.add_argument("--t1", type=float, default=19.99)
    ap.add_argument("--margin", type=float, default=2.5e-5)
    ap.add_argument("--outer", type=int, default=30)
    ap.add_argument("--free-outer", type=int, default=30)
    ap.add_argument("--inner", type=int, default=400)
    ap.add_argument("--rot-every", type=int, default=5)
    ap.add_argument("--ob-iters", type=int, default=300)
    ap.add_argument("--sweeps", type=int, default=40)
    ap.add_argument("--perturb", type=float, default=1e-3, help="relative jolt of the weights before each refit")
    ap.add_argument("--quick", action="store_true", help="pipeline test: short refits, 3 sweeps, 1 round")
    ap.add_argument("--wires", default="", help="r50: k:j:0:m,... (delta wires on address bits)")
    a = ap.parse_args()
    if a.wires:
        WIRES.extend(tuple(int(v) for v in s_.split(":")) for s_ in a.wires.split(","))
    if a.quick:
        a.outer, a.free_outer, a.inner, a.sweeps, a.rounds, a.ob_iters = 3, 3, 60, 3, 1, 40
    t0 = time.time()
    os.makedirs(a.out, exist_ok=True)

    def jolt(w_, seed):
        if a.perturb <= 0:
            return w_
        return w_ * (1.0 + a.perturb * np.random.default_rng(seed).standard_normal(w_.shape))

    with np.load(a.rounded) as z:
        words2, w = np.array(z["words"], np.int64), np.array(z["W"], np.float64)
        edges, U = np.array(z["tree_edges"], np.int64), np.array(z["U"], np.float64)
    R, B, C = w.shape
    N = edges.shape[0] + 1
    mod = 1 << a.brot
    full = np.zeros((N + R * B, N - 1), np.int64)
    full[N:] = np.mod(words2, mod)
    prob = cor.Problem(full, edges, w.shape, a.brot)
    dev = float(np.abs(np.abs(np.sum(prob.u * U.reshape(R, B, N), axis=2)) - 1.0).max())
    log(f"[table7] {a.rounded}: shape {(R, B, C)}, brot {a.brot}, words in [{int(words2.min())}, {int(words2.max())}] "
        f"(stored mod {mod}); decoded vectors agree with the rounding to {dev:.1e}  ({time.time() - t0:.0f} s)")
    assert dev < 1e-9
    steps = []

    def rows_step(w_, W_, iters, sweeps, tag):
        rows, e1 = build_rows7(h1_prime(prob, w_, W_), prob.edges, a.brot, iters, sweeps)
        set_rows(prob, rows)
        m = prob.official(w_, W_)
        steps.append(dict(step=tag, metrics=m, seconds=round(time.time() - t0, 1)))
        log(f"[table7] {tag}: rows' own error_1 {1e3 * e1:.4f} mHa; {json.dumps(m)}  ({time.time() - t0:.0f} s)")
        return m

    def save(w_, W_, m, done):
        out = dict(U=prob.u, W=w_, H2_I=W_ + w_.sum(axis=1), H1_I=prob.shifts["H1_I"],
                   Hsym_lower=prob.shifts["Hsym_lower"], tree_edges=edges, angle_words=prob.words)
        np.savez(os.path.join(a.out, "table.npz"), **out)
        with open(os.path.join(a.out, "table.json"), "w") as fh:
            json.dump(dict(rounded=os.path.abspath(a.rounded), shape=[R, B, C], brot=a.brot, final=m, complete=done,
                           t2=a.t2, t1=a.t1, steps=steps, seconds=round(time.time() - t0, 1)), fh, indent=1)

    set_rows(prob, np.zeros((N, N - 1), np.int64))
    m0 = prob.evaluate(w, np.zeros((R, C)), rematch=True)[0]
    W, gap, nit = offsets_min_gap(prob, w)
    m1 = prob.evaluate(w, W, rematch=True)[0]
    log(f"[table7] offsets: 0 -> E_gap {m0[2]:.4f}, Lambda {m0[3]:.2f}; least E_gap ({nit} iterations) -> E_gap "
        f"{gap:.6f} (floor {cor.GAP_FLOOR}), Lambda {m1[3]:.3f}, |W| {np.linalg.norm(W):.3f}  ({time.time() - t0:.0f} s)")
    steps.append(dict(step="offsets of least E_gap", E_gap=gap, Lambda=m1[3], seconds=round(time.time() - t0, 1)))
    m = rows_step(w, W, 0, 0, "provisional rows (nearest words)")
    log(f"[table7] THE ROUNDING'S OWN error_2 (official, its weights): {1e3 * m['error_2']:.4f} mHa = "
        f"{1e6 * m['error_2'] ** 2:.2f} mHa^2")
    log("[table7] --- free refit (no error_1 limit)")
    lam, w1, W1, m = cor.refit(prob, jolt(w, 100), W, a.t2, 1e6, a.margin, a.free_outer, a.inner, a.rot_every)
    if not np.isfinite(lam):
        log("[table7] NO FEASIBLE WEIGHTS: error_2 of this table cannot be brought under the limit")
        save(w, W, steps[-1]["metrics"], False)
        sys.exit(2)
    steps.append(dict(step="free refit", metrics=m, seconds=round(time.time() - t0, 1)))
    w, W = w1, W1
    for r in range(a.rounds):
        m = rows_step(w, W, a.ob_iters, a.sweeps, f"one-body rows, round {r + 1}")
        save(w, W, m, False)
        log(f"[table7] --- refit with the rows fixed, round {r + 1}")
        lam2, w2, W2, m2 = cor.refit(prob, jolt(w, 101 + r), W, a.t2, a.t1, a.margin, a.outer, a.inner, a.rot_every)
        if not np.isfinite(lam2):
            log("[table7] no feasible point with these rows (error_1 too high for the limit?)")
            sys.exit(3)
        steps.append(dict(step=f"refit round {r + 1}", metrics=m2, seconds=round(time.time() - t0, 1)))
        w, W, m = w2, W2, m2
        save(w, W, m, r + 1 == a.rounds)
    log(f"[table7] DONE {json.dumps(m)}  ({time.time() - t0:.0f} s) -> {os.path.join(a.out, 'table.npz')}")


if __name__ == "__main__":
    main()
