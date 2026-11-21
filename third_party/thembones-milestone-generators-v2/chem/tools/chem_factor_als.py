#!/usr/bin/env python
"""chem_factor: exact weight solves at fixed vectors (the linear part of the DF-THC fit), float64 numpy.

In weighted pair space the represented tensor is  G = sum_r Pi_r^T Z_r Pi_r  with Pi_r (B x P) the weighted pair
vectors of group r's projectors and Z_r = w_r w_r^T (B x B, PSD, rank <= C). With K_rs = Pi_r Pi_s^T =
(u_rb . u_sd)^2 and Q0_r = Pi_r Tw Pi_r^T,

    error_2^2 = sum_rs tr(Z_r K_rs Z_s K_sr) - 2 sum_r tr(Z_r Q0_r) + ||Tw||^2 .

For one group with the others fixed, Q_r = Q0_r - sum_{s != r} K_rs Z_s K_sr and Y = K^{1/2} Z_r K^{1/2}:
error^2 = ||Y - K^{-1/2} Q_r K^{-1/2}||^2 + const, so the optimal PSD rank-C Y is the top-C positive eigen-part
(Eckart-Young): an EXACT block solve. Sweeping the groups (Gauss-Seidel) decreases error_2 monotonically.

`solve_weights` is used by chem_factor_fit.py (phase `wsolve:K`) and by the diagnostics below:
  python chem_factor_als.py TABLE.npz [--sweeps 30] [--rcond 1e-10]
reports (1) error_2 after exact weight sweeps at rank C, (2) the same with the rank constraint dropped (C = B),
(3) with the PSD constraint dropped too, i.e. what these vectors could do with unlimited copies.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_factor_lib as L  # noqa: E402

N = 76


class Weights:
    def __init__(self, refc, u):
        """u (R, B, N) unit vectors."""
        self.R, self.B, _ = u.shape
        d = np.sqrt(np.where(refc.i == refc.j, 1.0, 2.0))
        self.Tw = d[:, None] * refc.T * d[None, :]
        self.T2 = float(np.sum(self.Tw ** 2))
        R, B = self.R, self.B
        flat = u.reshape(R * B, N)
        self.K = (flat @ flat.T) ** 2                                                   # (RB, RB)
        Pi = (flat[:, refc.i] * flat[:, refc.j]) * d[None, :]                           # (RB, P)
        TP = self.Tw @ Pi.T                                                             # (P, RB)
        self.Q0 = [Pi[r * B:(r + 1) * B] @ TP[:, r * B:(r + 1) * B] for r in range(R)]
        self.Q0 = [0.5 * (q + q.T) for q in self.Q0]

    def blk(self, r, s):
        B = self.B
        return self.K[r * B:(r + 1) * B, s * B:(s + 1) * B]

    def error(self, Z):
        R, B = self.R, self.B
        KZ = [np.concatenate([self.blk(s, r) @ Z[r] for s in range(R)], axis=0) for r in range(R)]   # (RB, B) each: K[:, r] Z_r
        e = self.T2
        for r in range(R):
            e -= 2 * np.sum(Z[r] * self.Q0[r])
            for s in range(R):
                # tr(Z_r K_rs Z_s K_sr) = sum (K_sr Z_r)^T .* (Z_s K_sr) -> use KZ
                e += np.sum(KZ[r][s * B:(s + 1) * B] * (Z[s] @ self.blk(s, r)))
        return float(np.sqrt(max(e, 0.0)))

    def sweep(self, Z, C, psd=True, rcond=1e-10, order=None):
        """One Gauss-Seidel sweep of exact group solves; returns the new list Z and factors w (R, B, C)."""
        R, B = self.R, self.B
        ws = []
        for r in (range(R) if order is None else order):
            Q = self.Q0[r].copy()
            for s in range(R):
                if s != r:
                    Krs = self.blk(r, s)
                    Q -= Krs @ Z[s] @ Krs.T
            sK, vK = np.linalg.eigh(self.blk(r, r))
            keep = sK > rcond * sK[-1]
            ih = vK[:, keep] / np.sqrt(sK[keep])                    # K^{-1/2} on the kept range, (B, k)
            M = ih.T @ Q @ ih
            lam, V = np.linalg.eigh(0.5 * (M + M.T))
            idx = np.argsort(-lam if psd else -np.abs(lam))[:C]
            if psd:
                idx = idx[lam[idx] > 0]
            fac = ih @ (V[:, idx] * np.sqrt(np.abs(lam[idx])))      # (B, <=C): Z = fac diag(sign) fac^T
            Z[r] = (fac * np.sign(lam[idx])) @ fac.T
            ws.append((r, fac))
        w = np.zeros((R, B, C))
        for r, fac in ws:
            w[r, :, :fac.shape[1]] = fac
        return Z, w


def solve_weights(refc, u, w, sweeps=10, rcond=1e-10, log=None):
    """Exact weight sweeps at rank C = w.shape[2] from the current weights; returns the new w and the error trace."""
    R, B, C = w.shape
    P = Weights(refc, u)
    Z = [w[r] @ w[r].T for r in range(R)]
    trace = [P.error(Z)]
    wn = w
    for k in range(sweeps):
        Z, wn = P.sweep(Z, C, rcond=rcond)
        trace.append(P.error(Z))
        if log:
            log(f"    weight sweep {k + 1}: error_2 {1e3 * trace[-1]:.4f} mHa")
        if trace[-2] - trace[-1] < 1e-9:
            break
    return wn, trace


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--sweeps", type=int, default=30)
    ap.add_argument("--rcond", type=float, default=1e-10)
    a = ap.parse_args()
    with np.load(a.table) as z:
        U, w = np.asarray(z["U"], np.float64), np.asarray(z["W"], np.float64)
    R, B, C = w.shape
    u = U / np.linalg.norm(U, axis=2, keepdims=True)
    refc = L.Reference()
    t0 = time.time()
    P = Weights(refc, u)
    Z0 = [w[r] @ w[r].T for r in range(R)]
    print(f"{a.table} {w.shape}: error_2 {1e3 * P.error(Z0):.4f} mHa (packed {1e3 * refc.error2_np(u, w):.4f}); |w| rms {np.sqrt(np.mean(w ** 2)):.3e} "
          f"max {np.abs(w).max():.3e}; D {0.5 * np.sum(np.abs(w).sum(1) ** 2 - w.sum(1) ** 2):.2f}")
    for r in range(R):
        s = np.linalg.eigvalsh(P.blk(r, r))
        print(f"  group {r}: projector Gram eigenvalues min {s[0]:.2e} max {s[-1]:.2f}, below 1e-6: {(s < 1e-6).sum()}, below 1e-3: {(s < 1e-3).sum()}")
    for label, C_, psd in ((f"rank {C} PSD (the real problem)", C, True), (f"rank {B} PSD (unlimited copies)", B, True),
                           (f"rank {B}, signs free (best block-diagonal fit of these vectors)", B, False)):
        Z = [z.copy() for z in Z0]
        prev = P.error(Z)
        out = []
        for k in range(a.sweeps):
            Z, wn = P.sweep(Z, C_, psd=psd, rcond=a.rcond)
            e = P.error(Z)
            if k in (0, 1, 2, 4, 9, 19, 29) or prev - e < 1e-9:
                out.append(f"{k + 1}: {1e3 * e:.4f}")
            if prev - e < 1e-9:
                break
            prev = e
        print(f"{label}: error_2 after sweeps  " + "  ".join(out) + f"   |w| max {np.abs(wn).max():.2e}  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
