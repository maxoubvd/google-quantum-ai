"""r48_chem_q7c_sdiag: the weight-compensated stiffness S_kk of every angle entry of a table (chem_own2; CPU).

S_kk = H_kk - H_kw H_ww^-1 H_wk: the Gauss-Newton curvature of error_2^2 along angle k when every other angle is
held and the weights of k's group are re-solved (the group's own weights only; the other groups are held). In
the linear model of tools/r48_chem_q7c_linmodel.py (one group, every angle decided progressively) deciding the
entries in order of S_kk, largest first, costs 0.33 of the tool's impact order and 0.71 of the plain stiffness
order |w|^2 mass. Written for r48_chem_q7c_round.py --rank sdiag --rank-file OUT.npy.

    python r48_chem_q7c_sdiag.py INIT.npz OUT.npy [--groups 0,1]
OUT.npy: (R*B, 75) float64 in the row order of INIT (the rounding tool applies its own row deal to it).
"""
import argparse
import os
import sys
import time

import numpy as np
import scipy.linalg as sl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_round import N, Target, angles_and_mass, log_depth_tree  # noqa: E402
from r48_chem_q7c_absorb import tangents  # noqa: E402


def group_sdiag(r, U, W, A, I, J, s, edges):
    R, B, C = W.shape
    u, V = tangents(U[r], edges)
    Tn = V.shape[1]
    p = (u[:, I] * u[:, J]) * s
    q = ((V[:, :, I] * u[:, None, J] + u[:, None, I] * V[:, :, J]) * s).reshape(B * Tn, -1)
    g = (W[r] @ A[r]) * s
    ar = A[r] * s
    bk = np.repeat(np.arange(B), Tn)
    qq = np.sum(q * q, axis=1)
    gg = np.sum(g * g, axis=1)
    qg = np.sum(q * g[bk], axis=1)
    Haa_d = 2.0 * (gg[bk] * qq + qg * qg)
    Gqp, Gqa = q @ p.T, q @ ar.T
    Ggp, Gga = g @ p.T, g @ ar.T
    Gpp, Gaa, Gap = p @ p.T, ar @ ar.T, ar @ p.T
    wb, wc = np.repeat(np.arange(B), C), np.tile(np.arange(C), B)
    Haw = 2.0 * (Gga[np.ix_(bk, wc)] * Gqp[:, wb] + Ggp[np.ix_(bk, wb)] * Gqa[:, wc])
    Hww = 2.0 * (Gaa[np.ix_(wc, wc)] * Gpp[np.ix_(wb, wb)] + Gap[np.ix_(wc, wb)] * Gap[np.ix_(wc, wb)].T)
    nw = B * C
    lam = 1e-10 * np.trace(Hww) / nw
    X = sl.cho_solve(sl.cho_factor(Hww + lam * np.eye(nw)), Haw.T)          # (nw, na)
    sd = Haa_d - np.einsum("kj,jk->k", Haw, X)
    return np.maximum(sd, 0.0).reshape(B, Tn), Haa_d.reshape(B, Tn)


_TGT = {}


def state_table(th_, w_, K, edges, cpu=None):
    """S_kk of every entry of a rounding state (r48_chem_q7c_round.py --rank sdyn, run in a spawned worker process so
    that the rounding's main thread keeps the GIL): angles th_ (+ K half turns, units of pi) and weights w_ (R, B, C),
    rows in the caller's order; returns (R*B, 75). cpu: optional core for this worker (os.sched_setaffinity)."""
    if cpu is not None:
        os.sched_setaffinity(0, {int(cpu)})
    from r48_chem_q7c_round import decode_np
    R, B, C = w_.shape
    if "t" not in _TGT:
        t_ = Target()
        _TGT["t"] = (t_.i, t_.j)
    I, J = _TGT["t"]
    s = np.sqrt(np.where(I == J, 1.0, 2.0))
    U_ = decode_np(th_ + K, edges).reshape(R, B, N)
    U_ = U_ / np.linalg.norm(U_, axis=2, keepdims=True)
    A_ = np.einsum("rbc,rbp->rcp", w_, U_[:, :, I] * U_[:, :, J])
    out = np.zeros((R, B, N - 1))
    for r in range(R):
        out[r] = group_sdiag(r, U_, w_, A_, I, J, s, edges)[0]
    return out.reshape(R * B, N - 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("init")
    ap.add_argument("out")
    ap.add_argument("--groups", default="")
    ap.add_argument("--labels", default=None)
    a = ap.parse_args()
    t0 = time.time()
    z = np.load(a.init)
    U, W = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
    R, B, C = W.shape
    U /= np.linalg.norm(U, axis=2, keepdims=True)
    tgt = Target()
    I, J = tgt.i, tgt.j
    s = np.sqrt(np.where(I == J, 1.0, 2.0))
    A = np.einsum("rbc,rbp->rcp", W, U[:, :, I] * U[:, :, J])
    edges = log_depth_tree(N)
    if a.labels:        # orbital labels of the tree, perm[orbital] = label (as r48_chem_q7c_round.py --labels)
        perm_ = np.load(a.labels).astype(np.int64)
        inv_ = np.empty(N, np.int64)
        inv_[perm_] = np.arange(N)
        edges = inv_[edges]
    out = np.zeros((R, B, N - 1))
    raw = np.zeros((R, B, N - 1))
    groups = [int(x) for x in a.groups.split(",")] if a.groups else list(range(R))
    for r in groups:
        sd, hd = group_sdiag(r, U, W, A, I, J, s, edges)
        out[r] = sd
        raw[r] = hd
        _, mass = angles_and_mass(U[r], edges)
        c = np.sum(W[r] ** 2, axis=1)[:, None] * mass
        rho = np.corrcoef(np.log(sd.ravel() + 1e-30), np.log(c.ravel() + 1e-30))[0, 1]
        print(f"group {r}: S_kk / H_kk median {np.median(sd / np.maximum(hd, 1e-300)):.3f}; corr(log S_kk, log |w|^2 mass) "
              f"{rho:.3f} ({time.time() - t0:.0f} s)", flush=True)
    np.save(a.out, out.reshape(R * B, N - 1))
    rawf = a.out[:-4] + "_raw.npy" if a.out.endswith(".npy") else a.out + "_raw.npy"
    np.save(rawf, raw.reshape(R * B, N - 1))          # the raw Gauss-Newton diagonal H_kk, same units
    print(f"wrote {a.out} and {rawf}")


if __name__ == "__main__":
    main()
