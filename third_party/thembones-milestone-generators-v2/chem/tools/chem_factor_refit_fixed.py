#!/usr/bin/env python
"""chem_factor: refit the weights w (npz 'W') and offsets (npz 'H2_I') of a dfthcsos submission with EVERY
angle word fixed (one-body rows included), so the published ANGLE circuits stay valid byte for byte.

Minimises log lambda_eff subject to error_2 <= t2, error_1 <= t1 and E_gap >= floor + margin with an augmented
Lagrangian (L-BFGS inner solves, exact numpy gradients), alternated with the pipeline's exact copy-column
rotations (they keep the two-body tensor, h1' and E_gap, and lower the L1 term of Lambda).

Differences from the pipeline's refit_al.py: the one-body rows are NOT regenerated, and error_1 at the fixed
one-body rows is a constraint (it has 8.7 mHa of unused budget on the public table).

  python chem_factor_refit_fixed.py IN_DIR OUT_DIR [--t2 19.995] [--t1 19.95] [--margin 5e-5] [--outer 30]
"""
import argparse
import json
import os
import shutil
import sys
import time

import numpy as np
from scipy.optimize import linear_sum_assignment, minimize

DATA = os.path.expanduser("~/qcc/data")
os.environ.setdefault("DFTHC_DATA", DATA)
PIPE = os.path.expanduser("~/qcc/src/marculera_pipeline/dfthc_refit")
sys.path[:0] = [os.path.join(DATA, "kickmix"), os.path.join(DATA, "verifiers"), PIPE]
import refit_weights as rw  # noqa: E402
from refine import GAP_FLOOR, rotate_columns  # noqa: E402
from dfthcsos.verify import evaluate_words, load_inputs  # noqa: E402


class Fixed:
    """error_2, error_1, E_gap, Lambda of (w, W) at fixed two-body vectors u and fixed one-body rows V."""

    def __init__(self, u, V, K, Q, T2, h1_base, energy):
        self.u, self.V, self.T2, self.h1_base, self.energy = u, V, T2, h1_base, energy
        R, B, N = u.shape
        self.R, self.B, self.N = R, B, N
        self.K = K.reshape(R * B, R * B)
        self.Q = Q
        self.perm = None

    def evaluate(self, w, W, rematch=False):
        R, B, N = self.R, self.B, self.N
        C = w.shape[2]
        # error_2
        KW = np.concatenate([self.K[:, r * B:(r + 1) * B] @ w[r] for r in range(R)], axis=1)       # (RB, RC)
        A = np.concatenate([w[r].T @ KW[r * B:(r + 1) * B] for r in range(R)], axis=0)             # (RC, RC)
        Qw = np.einsum('rbd,rdc->rbc', self.Q, w)
        E = np.sum(A * A) - 2 * np.sum(w * Qw) + self.T2
        e2 = float(np.sqrt(max(E, 1e-300)))
        dF = np.stack([4 * KW[r * B:(r + 1) * B] @ A[:, r * C:(r + 1) * C] for r in range(R)])
        de2 = (dF - 4 * Qw) / (2 * e2)
        # h1'
        k = np.einsum('rc,rbc->rb', W, w)
        H = self.h1_base - np.einsum('rb,rbp,rbq->pq', k, self.u, self.u, optimize=True)
        lam, vec = np.linalg.eigh(H)
        sgn = np.sign(lam)
        c = np.einsum('rbp,pq,rbq->rb', self.u, (vec * sgn) @ vec.T, self.u, optimize=True)
        gap = self.energy + np.abs(lam).sum() - lam.sum() + 0.5 * np.sum(W ** 2)
        s, l1 = w.sum(axis=1), np.abs(w).sum(axis=1)
        m = l1 + np.abs(W + s)
        Lam = np.abs(lam).sum() + 0.25 * np.sum(m ** 2)
        sg = np.sign(W + s)
        dgap = ((1 - c)[:, :, None] * W[:, None, :], np.einsum('rb,rbc->rc', 1 - c, w) + W)
        dLam = (-c[:, :, None] * W[:, None, :] + 0.5 * m[:, None, :] * (np.sign(w) + sg[:, None, :]),
                -np.einsum('rb,rbc->rc', c, w) + 0.5 * m * sg)
        # error_1 at the fixed one-body rows: || V^T diag(lam[perm]) V - H ||
        if self.perm is None or rematch:
            rows, cols = linear_sum_assignment(-np.abs(self.V @ vec))
            perm = np.empty(N, np.int64)
            perm[rows] = cols
            self.perm = perm
        D = (self.V.T * lam[self.perm]) @ self.V - H
        e1 = float(np.linalg.norm(D))
        g = np.einsum('jp,pq,jq->j', self.V, D, self.V)          # per one-body row
        gi = np.zeros(N)
        gi[self.perm] = g
        G1 = ((vec * gi) @ vec.T - D) / max(e1, 1e-300)           # d e1 / d H
        dk = -np.einsum('rbp,pq,rbq->rb', self.u, G1, self.u, optimize=True)
        de1 = (dk[:, :, None] * W[:, None, :], np.einsum('rb,rbc->rc', dk, w))
        return (e2, e1, float(gap), float(Lam)), (de2, de1, dgap, dLam)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--t2", type=float, default=19.995)
    ap.add_argument("--t1", type=float, default=19.95)
    ap.add_argument("--margin", type=float, default=5e-5)
    ap.add_argument("--outer", type=int, default=30)
    ap.add_argument("--inner", type=int, default=400)
    ap.add_argument("--rot-every", type=int, default=5)
    ap.add_argument("--brot", type=int, default=8)
    ap.add_argument("--rot-sweeps", type=int, default=2)
    ap.add_argument("--rot-hops", type=int, default=4)
    ap.add_argument("--rot-hops-always", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--rho0", type=float, default=1.0)
    ap.add_argument("--rho-max", type=float, default=400.0)
    ap.add_argument("--y0", default="0,0,0", help="starting multipliers for (error_2, error_1, gap), per mille of lambda_eff per mHa")
    a = ap.parse_args()
    log = lambda *x: print(*x, flush=True)  # noqa: E731
    with np.load(os.path.join(a.inp, "submission.npz")) as z:
        data = {k: np.array(z[k]) for k in z.files}
    brot = a.brot
    edges = np.asarray(data["tree_edges"], np.int64)
    words = np.asarray(data["angle_words"], np.int64)
    w0 = np.asarray(data["W"], np.float64)
    W0 = np.asarray(data["H2_I"], np.float64) - w0.sum(axis=1)
    R, B, C = w0.shape
    N = edges.shape[0] + 1
    Vall = rw.encode.decode(words, edges, brot)
    V, u = Vall[:N], Vall[N:].reshape(R, B, N)
    ref = rw.load_reference(rw.check_fcidump())
    with np.load(rw.SHIFTS) as z:
        _, _, S, beta = rw.native_inputs({'U': u, 'W': w0, 'H2_I': data['H2_I'], 'H1_I': z['H1_I'],
                                          'Hsym_lower': z['Hsym_lower']})
    K, Q, T2 = rw.error_terms(u, np.asarray(ref['tei'], np.float64), S)
    eta = int(ref['electrons'])
    fx = Fixed(u, V, K, Q, T2, ref['h1_thc'] - 0.5 * eta * S + beta * np.eye(N),
               float(rw.UPPER_ENERGY - ref['ecore'] + eta * beta))

    def official(w, W):
        d = {'U': data['U'], 'W': w, 'H2_I': W + w.sum(axis=1), 'H1_I': data['H1_I'], 'Hsym_lower': data['Hsym_lower']}
        m = evaluate_words(load_inputs(d), words, edges, brot)
        return {k: float(m[k]) for k in ('error_1', 'error_2', 'Lambda', 'E_gap', 'lambda_eff')}

    start = official(w0, W0)
    (e2, e1, gap, Lam), _ = fx.evaluate(w0, W0, rematch=True)
    log("start official", json.dumps(start))
    log(f"start formula  e2 {1e3 * e2:.5f} e1 {1e3 * e1:.5f} gap {gap:.7f} Lambda {Lam:.5f} lambda_eff {np.sqrt(gap * (2 * Lam - gap)):.5f}")
    assert abs(e2 - start['error_2']) < 1e-7 and abs(e1 - start['error_1']) < 1e-6, "formula mismatch"

    t2, t1, floor = 1e-3 * a.t2, 1e-3 * a.t1, GAP_FLOOR + a.margin
    w, W = w0.copy(), W0.copy()
    best = (start['lambda_eff'], w0, W0, start)
    y = np.array([float(v) for v in a.y0.split(",")])
    rho = a.rho0
    t0 = time.time()
    for it in range(a.outer):
        if a.rot_every and it % a.rot_every == 0:
            w, W = rotate_columns(w, W, sweeps=a.rot_sweeps, seed=it + a.seed, hops=a.rot_hops if (it == 0 or a.rot_hops_always) else 0)
            (e2, e1, gap, Lam), _ = fx.evaluate(w, W, rematch=True)
            log(f"  rotate: lambda_eff {np.sqrt(gap * (2 * Lam - gap)):.5f} e2 {1e3 * e2:.4f} e1 {1e3 * e1:.4f} gap {gap:.6f}")
        scale = np.concatenate([np.full(w.size, np.std(w)), np.full(W.size, np.std(W))])

        def lagr(zv):
            x = zv * scale
            ww, WW = x[:w.size].reshape(w.shape), x[w.size:].reshape(W.shape)
            (e2, e1, gap, Lam), (de2, de1, dgap, dLam) = fx.evaluate(ww, WW)
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
        res = minimize(lagr, z0, jac=True, method='L-BFGS-B',
                       options={'maxiter': a.inner, 'maxcor': 30, 'ftol': 1e-14, 'gtol': 1e-8, 'maxls': 40})
        x = res.x * scale
        w, W = x[:w.size].reshape(w.shape), x[w.size:].reshape(W.shape)
        (e2, e1, gap, Lam), _ = fx.evaluate(w, W, rematch=True)
        cons = np.array([1e3 * (e2 - t2), 1e3 * (e1 - t1), 1e3 * (floor - gap)])
        lam_eff = float(np.sqrt(gap * (2 * Lam - gap)))
        log(f"outer {it}: {res.nit} its  e2 {1e3 * e2:.4f}  e1 {1e3 * e1:.4f}  gap {gap:.6f}  Lambda {Lam:.4f}  "
            f"lambda_eff {lam_eff:.5f}  y {y.round(3)} rho {rho:.0f}  ({time.time() - t0:.0f} s)")
        if cons.max() <= 0:
            m = official(w, W)
            ok = m['error_2'] <= t2 + 1e-10 and m['error_1'] <= t1 + 1e-10 and m['E_gap'] >= floor - 1e-10
            if ok and m['lambda_eff'] < best[0]:
                best = (m['lambda_eff'], w.copy(), W.copy(), m)
                log("   feasible, official", json.dumps(m))
        y = np.maximum(0.0, y + rho * cons)
        if cons.max() > 0.002:
            rho = min(rho * 2.0, a.rho_max)
    lam_eff, w, W, m = best
    log("BEST", json.dumps(m))
    os.makedirs(a.out, exist_ok=True)
    for f in os.listdir(a.inp):
        if f.endswith(".kmb"):
            shutil.copy(os.path.join(a.inp, f), os.path.join(a.out, f))
    out = dict(data)
    out['W'] = w
    out['H2_I'] = W + w.sum(axis=1)
    np.savez(os.path.join(a.out, "submission.npz"), **out)
    log("saved", a.out)


if __name__ == "__main__":
    main()
