#!/usr/bin/env python
"""chem_own_table: from OUR rounded two-body table to a complete table (one-body rows, weights, offsets). CPU.

    chem_own_table.py <rounded.npz> <out_dir> [--rounds 2] [--t2 19.9985] [--t1 19.99] [--outer 30] [--inner 400]
                      [--sweeps 40] [--free-outer 30] [--quick]

<rounded.npz>: output of tools/chem_own_round.py (words of the R*B two-body rows, the weights of the rounding,
the codes and the carrier pattern). Nothing else is read except the hosts' verifier data.
<out_dir>/table.npz: a complete submission.npz (angle_words with the 76 one-body rows first, tree_edges, U, W,
H2_I, H1_I, Hsym_lower), the input of tools/chem_own_entry.py --weights keep; <out_dir>/table.json: the hosts'
metrics after every step.

Steps:
  0. checks: the codes and the carrier pattern hold on the two-body rows; the words decode to the rounding's
     vectors.
  1. the weights the table wants: provisional one-body rows for h1'(w, offsets 0) (coded, not polished), then
     the refit of weights and offsets WITHOUT the error_1 limit. (The offsets move h1' by far more than fixed
     one-body rows allow, so they are found first and the rows are built for them afterwards.)
  2. `rounds` times: one-body rows for the current h1' (tools/chem_own_onebody.py: order, carriers, codes,
     shared block, exact polish), then the refit with these rows fixed and both error limits on.
Every reported number is the hosts' own evaluation (chem_own_refit.Problem.official).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_own_onebody as ob  # noqa: E402
import chem_own_refit as cor  # noqa: E402
import chem_own_round as RD  # noqa: E402


def log(*x):
    print(*x, flush=True)


def offsets_min_gap(prob, w, W0=None, iters=600):
    """Offsets that minimise E_gap = E0 + 2 sum(negative eigenvalues of h1') + sum W^2 / 2 at fixed weights:
    the tightest lower bound these weights allow. The problem is convex (the sum of the negative parts of the
    eigenvalues of a matrix affine in W, plus a quadratic), so L-BFGS from zero reaches its minimum; the
    gradient is the one of chem_own_refit.Problem.evaluate. Returns (W, E_gap)."""
    from scipy.optimize import minimize
    R, B, C = w.shape

    def f(x):
        (_, _, gap, _), (_, _, dgap, _) = prob.evaluate(w, x.reshape(R, C))
        return gap, dgap[1].ravel()

    x0 = np.zeros(R * C) if W0 is None else np.asarray(W0, np.float64).ravel()
    res = minimize(f, x0, jac=True, method="L-BFGS-B", options={"maxiter": iters, "maxcor": 30, "ftol": 1e-15, "gtol": 1e-9})
    return res.x.reshape(R, C), float(res.fun), int(res.nit)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rounded")
    ap.add_argument("out")
    ap.add_argument("--rounds", type=int, default=2)
    ap.add_argument("--t2", type=float, default=19.9985)
    ap.add_argument("--t1", type=float, default=19.99)
    ap.add_argument("--margin", type=float, default=2.5e-5)
    ap.add_argument("--outer", type=int, default=30)
    ap.add_argument("--free-outer", type=int, default=30)
    ap.add_argument("--inner", type=int, default=400)
    ap.add_argument("--rot-every", type=int, default=5)
    ap.add_argument("--sweeps", type=int, default=40)
    ap.add_argument("--quick", action="store_true", help="pipeline test: 3 outer steps of 60 iterations, 3 sweeps, 1 round")
    a = ap.parse_args()
    if a.quick:
        a.outer, a.free_outer, a.inner, a.sweeps, a.rounds = 3, 3, 60, 3, 1
    t0 = time.time()
    os.makedirs(a.out, exist_ok=True)
    with np.load(a.rounded) as z:
        words2, w = np.array(z["words"], np.int64), np.array(z["W"], np.float64)
        deg, xdeg = np.array(z["codes"]), np.array(z["xcodes"])
        masks = [int(m) for m in z["carriers"]]
        edges, U = np.array(z["tree_edges"], np.int64), np.array(z["U"], np.float64)
    R, B, C = w.shape
    N = edges.shape[0] + 1
    nhead = (-N) % 16
    addr = N + np.arange(R * B)
    RD.check_codes(ob.signed_words(words2, RD.carrier_pattern(addr, masks)), addr, deg, xdeg)
    full = np.zeros((N + R * B, N - 1), np.int64)
    full[N:] = words2
    prob = cor.Problem(full, edges, w.shape, 8)
    dev = float(np.abs(np.abs(np.sum(prob.u * U.reshape(R, B, N), axis=2)) - 1.0).max())
    log(f"[table] {a.rounded}: shape {(R, B, C)}, carriers {masks}, {int((deg > 0).sum())} coded column-bits, "
        f"{int((xdeg > 0).sum())} XOR-coded batch-bits; codes and carriers hold on the two-body rows; decoded "
        f"vectors agree with the rounding to {dev:.1e}  ({time.time() - t0:.0f} s)")
    assert dev < 1e-9
    steps = []

    def rows_step(w_, W_, sweeps, tag):
        rows, _ = ob.build_rows(ob.h1_prime(prob, w_, W_), prob.edges, prob.words[N:N + nhead], deg, xdeg, masks,
                                sweeps, log=log)
        ob.set_rows(prob, rows)
        m = prob.official(w_, W_)
        steps.append(dict(step=tag, metrics=m, seconds=round(time.time() - t0, 1)))
        log(f"[table] {tag}: {json.dumps(m)}  ({time.time() - t0:.0f} s)")
        return m

    def save(w_, W_, m, done):
        out = dict(U=prob.u, W=w_, H2_I=W_ + w_.sum(axis=1), H1_I=prob.shifts["H1_I"],
                   Hsym_lower=prob.shifts["Hsym_lower"], tree_edges=edges, angle_words=prob.words)
        np.savez(os.path.join(a.out, "table.npz"), **out)
        with open(os.path.join(a.out, "table.json"), "w") as fh:
            json.dump(dict(rounded=os.path.abspath(a.rounded), shape=[R, B, C], carriers=masks, final=m, complete=done,
                           t2=a.t2, t1=a.t1, steps=steps, seconds=round(time.time() - t0, 1)), fh, indent=1)

    # 1. the weights the table wants. Offsets first: with offsets 0 the one-body part is the bare h1 (sum|lam| of
    #    hundreds of Ha, lambda_eff above 1,000), so the start is the offsets of the tightest lower bound.
    ob.set_rows(prob, np.zeros((N, N - 1), np.int64))
    m0 = prob.evaluate(w, np.zeros((R, C)), rematch=True)[0]
    W, gap, nit = offsets_min_gap(prob, w)
    m1 = prob.evaluate(w, W, rematch=True)[0]
    log(f"[table] offsets: 0 -> E_gap {m0[2]:.4f}, Lambda {m0[3]:.2f}; least E_gap ({nit} iterations) -> E_gap "
        f"{gap:.6f} (floor {cor.GAP_FLOOR}), Lambda {m1[3]:.3f}, |W| {np.linalg.norm(W):.3f}  ({time.time() - t0:.0f} s)")
    steps.append(dict(step="offsets of least E_gap", E_gap=gap, Lambda=m1[3], seconds=round(time.time() - t0, 1)))
    rows_step(w, W, 0, "provisional rows")
    log("[table] --- free refit (no error_1 limit)")
    lam, w1, W1, m = cor.refit(prob, w, W, a.t2, 1e6, a.margin, a.free_outer, a.inner, a.rot_every)
    if not np.isfinite(lam):
        log("[table] NO FEASIBLE WEIGHTS: error_2 of this table cannot be brought under the limit")
        save(w, W, steps[-1]["metrics"], False)
        sys.exit(2)
    steps.append(dict(step="free refit", metrics=m, seconds=round(time.time() - t0, 1)))
    w, W = w1, W1
    # 2. rows for these weights, then the refit with the rows fixed
    for r in range(a.rounds):
        m = rows_step(w, W, a.sweeps, f"one-body rows, round {r + 1}")
        save(w, W, m, False)
        log(f"[table] --- refit with the rows fixed, round {r + 1}")
        lam2, w2, W2, m2 = cor.refit(prob, w, W, a.t2, a.t1, a.margin, a.outer, a.inner, a.rot_every)
        if not np.isfinite(lam2):
            log("[table] no feasible point with these rows (error_1 too high for the limit?)")
            sys.exit(3)
        steps.append(dict(step=f"refit round {r + 1}", metrics=m2, seconds=round(time.time() - t0, 1)))
        w, W, m = w2, W2, m2
        save(w, W, m, r + 1 == a.rounds)
    log(f"[table] DONE {json.dumps(m)}  ({time.time() - t0:.0f} s) -> {os.path.join(a.out, 'table.npz')}")


if __name__ == "__main__":
    main()
