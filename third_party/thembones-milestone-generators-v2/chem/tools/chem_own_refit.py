#!/usr/bin/env python
"""chem_own_refit: the weights of a dfthcsos entry for a FIXED table of angle words (CPU, numpy/scipy only).

    chem_own_refit.py <table.npz> <out.npz> [--start table|scratch] [--t2 19.9985] [--t1 19.99]
                      [--margin 2.5e-5] [--outer 30] [--inner 400] [--rot-every 5] [--stages 2] [--quick]

<table.npz> holds `angle_words` (76 + R*B rows x 75, brot bits), `tree_edges`, and `U`, `W`, `H2_I` (the shape
(R, B, C) is read from W; with --start table its values are the starting point, with --start scratch they are
ignored and the weights are built from the decoded vectors alone). <out.npz> is a complete submission.npz: the
same words and tree, `U` = the decoded two-body vectors, the refit `W` and `H2_I`, and the hosts' `H1_I`,
`Hsym_lower`.

What it does. With every word fixed the two-body vectors u and the one-body rows V are fixed, and the metric is
a function of the weights w (R, B, C) and the offsets W (R, C) only:
    error_2^2 = || w^T K w ||_F^2 - 2 <w, Q w> + T2,   K = (u u^T)^2 elementwise, Q_r = (u u)^T T (u u)
    h1'      = h1_base - sum_rb (sum_c W_rc w_rbc) u_rb u_rb^T,   error_1 = || V^T diag(lam_perm) V - h1' ||
    E_gap    = E0 + sum|lam| - sum lam + 1/2 sum W^2,   Lambda = sum|lam| + 1/4 sum_rc (|W + s| + |w|_1)^2
It minimises log lambda_eff = 1/2 log(E_gap (2 Lambda - E_gap)) subject to error_2 <= t2, error_1 <= t1 and
E_gap >= floor + margin with an augmented Lagrangian (L-BFGS-B inner solves, exact gradients), alternated with
exact rotations of pairs of copy columns inside a group (they leave the two-body tensor, h1' and E_gap
unchanged and lower the L1 part of Lambda). Every reported number is re-evaluated with the hosts' own
`evaluate_words`.

Provenance. The formulas and the augmented-Lagrangian loop are those of tools/chem_factor_refit_fixed.py (ours),
re-written here without any import from the public pipeline: K, Q, T2 and the column rotations, which that tool
took from `~/qcc/src/marculera_pipeline/dfthc_refit` (Apache-2.0), are implemented below from the definitions.
The idea of rotating copy columns to reduce the one-norm is public (James Barlow / manwithacat "gauge
rotations", used in marculera's pipeline); the implementation (grid + golden section per column pair) is ours.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import linear_sum_assignment, minimize

DATA = os.environ.get("QCC_DATA", os.path.expanduser("~/qcc/data"))
sys.path[:0] = [os.path.join(DATA, "kickmix"), os.path.join(DATA, "verifiers")]

GAP_FLOOR = 3.914898
UPPER_ENERGY = -22140.27335


def log(*x):
    print(*x, flush=True)


def load_table(path):
    with np.load(path, allow_pickle=False) as z:
        return {k: np.array(z[k]) for k in z.files}


def host_shifts():
    path = os.path.join(DATA, "verifiers", "dfthcsos", "FeMoCo76o113e.fixed_shifts.npz")
    with np.load(path, allow_pickle=False) as z:
        return {k: np.array(z[k]) for k in ("H1_I", "Hsym_lower")}


class Problem:
    """Fixed vectors of a word table and everything the metric needs from the reference integrals."""

    def __init__(self, words, edges, shape, brot=8):
        from dfthcsos import verify as V
        from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
        self.V_mod = V
        self.words = np.asarray(words, np.int64)
        self.edges = np.asarray(edges, np.int64)
        self.brot = brot
        R, B, C = shape
        N = self.edges.shape[0] + 1
        self.R, self.B, self.C, self.N = R, B, C, N
        assert self.words.shape == (N + R * B, N - 1), "table shape does not match (R, B)"
        allv = reconstruct_u_vectors(self.words / (1 << (brot - 1)), self.edges)
        allv = allv / np.linalg.norm(allv, axis=1, keepdims=True)
        self.Vrows = allv[:N]
        self.u = allv[N:].reshape(R, B, N)
        self.shifts = host_shifts()
        dummy = {"U": self.u, "W": np.zeros((R, B, C)), "H2_I": np.zeros((R, C)), **self.shifts}
        inputs = V.load_inputs(dummy)
        ref = inputs.ref
        self.S = np.asarray(inputs.S, np.float64)
        self.beta = float(inputs.beta_1)
        self.eta = int(ref["electrons"])
        tei = np.asarray(ref["tei"], np.float64)
        oei = np.asarray(ref["oei"], np.float64)
        self.h1_thc = oei - 0.5 * np.einsum("prrq->pq", tei)
        self.h1_base = self.h1_thc - 0.5 * self.eta * self.S + self.beta * np.eye(N)
        self.energy = float(UPPER_ENERGY - ref["ecore"] + self.eta * self.beta)
        # two-body target T = tei + (S (x) I + I (x) S) / 2 as an N^2 x N^2 matrix
        Tm = tei.reshape(N * N, N * N).copy()
        del tei
        sv, iv = self.S.reshape(-1), np.eye(N).reshape(-1)
        Tm += 0.5 * (np.outer(sv, iv) + np.outer(iv, sv))
        self.T2 = float(np.sum(Tm * Tm))
        uu = (self.u[:, :, :, None] * self.u[:, :, None, :]).reshape(R, B, N * N)
        self.Q = np.stack([uu[r] @ (Tm @ uu[r].T) for r in range(R)])           # (R, B, B)
        del Tm
        flat = self.u.reshape(R * B, N)
        self.K = (flat @ flat.T) ** 2                                             # (RB, RB)
        self.perm = None

    # ---- metric and exact gradients at fixed vectors
    def evaluate(self, w, W, rematch=False):
        R, B, N, C = self.R, self.B, self.N, w.shape[2]
        KW = np.concatenate([self.K[:, r * B:(r + 1) * B] @ w[r] for r in range(R)], axis=1)       # (RB, RC)
        A = np.concatenate([w[r].T @ KW[r * B:(r + 1) * B] for r in range(R)], axis=0)             # (RC, RC)
        Qw = np.einsum("rbd,rdc->rbc", self.Q, w)
        E = np.sum(A * A) - 2 * np.sum(w * Qw) + self.T2
        e2 = float(np.sqrt(max(E, 1e-300)))
        dF = np.stack([4 * KW[r * B:(r + 1) * B] @ A[:, r * C:(r + 1) * C] for r in range(R)])
        de2 = (dF - 4 * Qw) / (2 * e2)
        k = np.einsum("rc,rbc->rb", W, w)
        H = self.h1_base - np.einsum("rb,rbp,rbq->pq", k, self.u, self.u, optimize=True)
        lam, vec = np.linalg.eigh(H)
        sgn = np.sign(lam)
        c = np.einsum("rbp,pq,rbq->rb", self.u, (vec * sgn) @ vec.T, self.u, optimize=True)
        gap = self.energy + np.abs(lam).sum() - lam.sum() + 0.5 * np.sum(W ** 2)
        s, l1 = w.sum(axis=1), np.abs(w).sum(axis=1)
        m = l1 + np.abs(W + s)
        Lam = np.abs(lam).sum() + 0.25 * np.sum(m ** 2)
        sg = np.sign(W + s)
        dgap = ((1 - c)[:, :, None] * W[:, None, :], np.einsum("rb,rbc->rc", 1 - c, w) + W)
        dLam = (-c[:, :, None] * W[:, None, :] + 0.5 * m[:, None, :] * (np.sign(w) + sg[:, None, :]),
                -np.einsum("rb,rbc->rc", c, w) + 0.5 * m * sg)
        if self.perm is None or rematch:
            rows, cols = linear_sum_assignment(-np.abs(self.Vrows @ vec))
            perm = np.empty(N, np.int64)
            perm[rows] = cols
            self.perm = perm
        D = (self.Vrows.T * lam[self.perm]) @ self.Vrows - H
        e1 = float(np.linalg.norm(D))
        g = np.einsum("jp,pq,jq->j", self.Vrows, D, self.Vrows)
        gi = np.zeros(N)
        gi[self.perm] = g
        G1 = ((vec * gi) @ vec.T - D) / max(e1, 1e-300)
        dk = -np.einsum("rbp,pq,rbq->rb", self.u, G1, self.u, optimize=True)
        de1 = (dk[:, :, None] * W[:, None, :], np.einsum("rb,rbc->rc", dk, w))
        return (e2, e1, float(gap), float(Lam)), (de2, de1, dgap, dLam)

    def official(self, w, W):
        d = {"U": self.u, "W": w, "H2_I": W + w.sum(axis=1), **self.shifts}
        m = self.V_mod.evaluate_words(self.V_mod.load_inputs(d), self.words, self.edges, self.brot)
        return {k: float(m[k]) for k in ("error_1", "error_2", "Lambda", "E_gap", "lambda_eff")}

    # ---- weights from the vectors alone
    def scratch(self, sweeps=6):
        """Least-squares weights without any starting point: block sweeps over the groups; with the other
        groups fixed, the best rank-C PSD Gram matrix G_r = w_r w_r^T of a group is the Eckart-Young
        truncation of its target in the metric K_rr (ours: tools/chem_factor_als.py)."""
        R, B, C = self.R, self.B, self.C
        G = [np.zeros((B, B)) for _ in range(R)]
        w = np.zeros((R, B, C))
        Kb = [[self.K[r * B:(r + 1) * B, s * B:(s + 1) * B] for s in range(R)] for r in range(R)]
        for sweep in range(sweeps):
            for r in range(R):
                # minimise || K^(1/2) G K^(1/2) - M ||, M = K^(-1/2) (Q_r - sum_{s != r} K_rs G_s K_sr) K^(-1/2)
                rhs = self.Q[r].copy()
                for s in range(R):
                    if s != r:
                        rhs -= Kb[r][s] @ G[s] @ Kb[s][r]
                ev, evec = np.linalg.eigh(Kb[r][r])
                ev = np.maximum(ev, 1e-12 * ev.max())
                Kih = (evec / np.sqrt(ev)) @ evec.T
                M = Kih @ rhs @ Kih
                M = 0.5 * (M + M.T)
                mv, mvec = np.linalg.eigh(M)
                top = np.argsort(mv)[::-1][:C]
                val = np.maximum(mv[top], 0.0)
                wr = Kih @ (mvec[:, top] * np.sqrt(val))
                G[r] = wr @ wr.T
                w[r] = wr
        return w


def rotate_columns(w, W, sweeps=2, grid=48):
    """Rotate pairs of copy columns of a group (w_r, W_r + s_r rotate together) to lower
    sum_c (|W_c + s_c| + |w_c|_1)^2. Exact invariances: the two-body tensor, h1', sum W^2."""
    R, B, C = w.shape
    w = w.copy()
    wB = W + w.sum(axis=1)                       # (R, C)
    th = (np.arange(grid) + 0.5) * (np.pi / 2) / grid
    cg, sg = np.cos(th), np.sin(th)

    def cost(x1, x2, a, b, c, s):
        y1 = np.abs(np.outer(c, x1) + np.outer(s, x2)).sum(axis=1) + np.abs(c * a + s * b)
        y2 = np.abs(np.outer(c, x2) - np.outer(s, x1)).sum(axis=1) + np.abs(c * b - s * a)
        return y1 ** 2 + y2 ** 2

    for _ in range(sweeps):
        for r in range(R):
            for c1 in range(C):
                for c2 in range(c1 + 1, C):
                    x1, x2, a, b = w[r, :, c1], w[r, :, c2], wB[r, c1], wB[r, c2]
                    base = (np.abs(x1).sum() + abs(a)) ** 2 + (np.abs(x2).sum() + abs(b)) ** 2
                    f = cost(x1, x2, a, b, cg, sg)
                    i = int(np.argmin(f))
                    if f[i] >= base - 1e-12:
                        continue
                    lo, hi = th[i] - th[0] * 2, th[i] + th[0] * 2
                    for _ in range(18):          # golden-section refinement on the bracket
                        m1, m2 = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
                        f1 = cost(x1, x2, a, b, np.array([np.cos(m1)]), np.array([np.sin(m1)]))[0]
                        f2 = cost(x1, x2, a, b, np.array([np.cos(m2)]), np.array([np.sin(m2)]))[0]
                        if f1 < f2:
                            hi = m2
                        else:
                            lo = m1
                    t = 0.5 * (lo + hi)
                    ct, st = np.cos(t), np.sin(t)
                    if cost(x1, x2, a, b, np.array([ct]), np.array([st]))[0] >= base - 1e-12:
                        continue
                    w[r, :, c1], w[r, :, c2] = ct * x1 + st * x2, ct * x2 - st * x1
                    wB[r, c1], wB[r, c2] = ct * a + st * b, ct * b - st * a
    return w, wB - w.sum(axis=1)


def refit(prob: Problem, w, W, t2=19.9985, t1=19.99, margin=2.5e-5, outer=30, inner=400, rot_every=5,
          rho0=1.0, rho_max=2000.0, tighten=0.003):
    """Returns (lambda_eff, w, W, official metrics) of the best point that satisfies error_2 <= t2,
    error_1 <= t1, E_gap >= floor + margin (official evaluation). The augmented Lagrangian aims `tighten` mHa
    inside the two error limits (and 2e-6 Ha inside the gap limit), because its iterates approach a
    constraint from the infeasible side."""
    lim2, lim1, limg = 1e-3 * t2, 1e-3 * t1, GAP_FLOOR + margin
    t2, t1, floor = 1e-3 * (t2 - tighten), 1e-3 * (t1 - tighten), GAP_FLOOR + margin + 2e-6
    start = prob.official(w, W)
    (e2, e1, gap, Lam), _ = prob.evaluate(w, W, rematch=True)
    log("start official", json.dumps(start))
    log(f"start formula  e2 {1e3 * e2:.5f} e1 {1e3 * e1:.5f} gap {gap:.7f} Lambda {Lam:.5f} "
        f"lambda_eff {np.sqrt(max(gap * (2 * Lam - gap), 0)):.5f}")
    assert abs(e2 - start["error_2"]) < 1e-7 and abs(e1 - start["error_1"]) < 1e-6, "formula mismatch"
    feasible0 = start["error_2"] <= lim2 and start["error_1"] <= lim1 and start["E_gap"] >= limg
    best = (start["lambda_eff"] if feasible0 else np.inf, w.copy(), W.copy(), start)
    y = np.zeros(3)
    rho = rho0
    t0 = time.time()
    for it in range(outer):
        if rot_every and it % rot_every == 0:
            w, W = rotate_columns(w, W)
            (e2, e1, gap, Lam), _ = prob.evaluate(w, W, rematch=True)
            log(f"  rotate: lambda_eff {np.sqrt(max(gap * (2 * Lam - gap), 0)):.5f} e2 {1e3 * e2:.4f} "
                f"e1 {1e3 * e1:.4f} gap {gap:.6f}  ({time.time() - t0:.0f} s)")
        scale = np.concatenate([np.full(w.size, np.std(w)), np.full(W.size, max(np.std(W), 1e-6))])

        def lagr(zv):
            x = zv * scale
            ww, WW = x[:w.size].reshape(w.shape), x[w.size:].reshape(W.shape)
            (e2, e1, gap, Lam), (de2, de1, dgap, dLam) = prob.evaluate(ww, WW)
            rad = gap * (2 * Lam - gap)
            if rad <= 0:
                return 1e12, np.zeros_like(zv)
            f = 500.0 * np.log(rad)
            gw, gW = [1e3 * ((Lam - gap) * g + gap * l) / rad for g, l in zip(dgap, dLam)]
            cons = np.array([1e3 * (e2 - t2), 1e3 * (e1 - t1), 1e3 * (floor - gap)])
            p = np.maximum(0.0, y + rho * cons)
            f += float(np.sum(p ** 2 - y ** 2) / (2 * rho))
            gw = gw + 1e3 * p[0] * de2 + 1e3 * p[1] * de1[0] - 1e3 * p[2] * dgap[0]
            gW = gW + 1e3 * p[1] * de1[1] - 1e3 * p[2] * dgap[1]
            return f, np.concatenate([gw.ravel(), gW.ravel()]) * scale

        z0 = np.concatenate([w.ravel(), W.ravel()]) / scale
        res = minimize(lagr, z0, jac=True, method="L-BFGS-B",
                       options={"maxiter": inner, "maxcor": 30, "ftol": 1e-14, "gtol": 1e-8, "maxls": 40})
        x = res.x * scale
        w, W = x[:w.size].reshape(w.shape), x[w.size:].reshape(W.shape)
        (e2, e1, gap, Lam), _ = prob.evaluate(w, W, rematch=True)
        cons = np.array([1e3 * (e2 - t2), 1e3 * (e1 - t1), 1e3 * (floor - gap)])
        lam_eff = float(np.sqrt(max(gap * (2 * Lam - gap), 0)))
        log(f"outer {it}: {res.nit} its  e2 {1e3 * e2:.4f}  e1 {1e3 * e1:.4f}  gap {gap:.6f}  Lambda {Lam:.4f}  "
            f"lambda_eff {lam_eff:.5f}  y {y.round(3)} rho {rho:.0f}  ({time.time() - t0:.0f} s)")
        if e2 <= lim2 and e1 <= lim1 and gap >= limg:
            m = prob.official(w, W)
            ok = m["error_2"] <= lim2 and m["error_1"] <= lim1 and m["E_gap"] >= limg
            if ok and m["lambda_eff"] < best[0]:
                best = (m["lambda_eff"], w.copy(), W.copy(), m)
                log("   feasible, official", json.dumps(m))
        y = np.maximum(0.0, y + rho * cons)
        if cons.max() > 0.002:
            rho = min(rho * 2.0, rho_max)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table")
    ap.add_argument("out")
    ap.add_argument("--start", default="table", choices=["table", "scratch"])
    ap.add_argument("--t2", type=float, default=19.9985)
    ap.add_argument("--t1", type=float, default=19.99)
    ap.add_argument("--margin", type=float, default=2.5e-5)
    ap.add_argument("--outer", type=int, default=30)
    ap.add_argument("--inner", type=int, default=400)
    ap.add_argument("--rot-every", type=int, default=5)
    ap.add_argument("--brot", type=int, default=8)
    ap.add_argument("--stages", type=int, default=2)
    ap.add_argument("--quick", action="store_true", help="2 outer steps of 50 iterations (pipeline test)")
    a = ap.parse_args()
    if a.quick:
        a.outer, a.inner, a.stages = 2, 50, 1
    t0 = time.time()
    data = load_table(a.table)
    shape = data["W"].shape
    prob = Problem(data["angle_words"], data["tree_edges"], shape, a.brot)
    log(f"problem built: shape {shape}, {prob.words.shape[0]} rows ({time.time() - t0:.0f} s)")
    if a.start == "table":
        w = np.asarray(data["W"], np.float64)
        W = np.asarray(data["H2_I"], np.float64) - w.sum(axis=1)
    else:
        w = prob.scratch()
        W = np.zeros((shape[0], shape[2]))
        log("scratch weights:", json.dumps(prob.official(w, W)))
    lam, m = np.inf, None
    for stage in range(a.stages):
        # every stage restarts the multipliers from the best feasible point so far (a continuation: the first
        # stage ends stiff, with a large penalty weight, and a fresh start finds more descent)
        log(f"--- stage {stage + 1} of {a.stages}")
        lam2, w2, W2, m2 = refit(prob, w, W, a.t2, a.t1, a.margin, a.outer, a.inner, a.rot_every)
        if not np.isfinite(lam2):
            break
        if lam2 < lam:
            lam, w, W, m = lam2, w2, W2, m2
        if stage and lam2 > lam - 1e-6:
            break
    if not np.isfinite(lam):
        log("NO FEASIBLE POINT: the table does not pass the gates with these weights")
        sys.exit(2)
    log("BEST", json.dumps(m))
    out = dict(U=prob.u, W=w, H2_I=W + w.sum(axis=1), H1_I=prob.shifts["H1_I"], Hsym_lower=prob.shifts["Hsym_lower"],
               tree_edges=data["tree_edges"], angle_words=data["angle_words"])
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez(a.out, **out)
    with open(os.path.splitext(a.out)[0] + ".refit.json", "w") as fh:
        json.dump(dict(metrics=m, start=a.start, t2=a.t2, t1=a.t1, margin=a.margin, outer=a.outer, inner=a.inner,
                       seconds=round(time.time() - t0, 1)), fh, indent=1)
    log("saved", a.out, f"({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
