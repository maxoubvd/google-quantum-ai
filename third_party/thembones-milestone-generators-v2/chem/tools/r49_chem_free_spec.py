"""r49_chem_free_spec: the sensitivity spectrum of the two-body WORDS of a continuous fit point, and the price of
tying one edge's word to another's in the same row (R49-CHEM-FREE kill test 1; CPU, one group at a time).

For every angle k = (row b, tree edge t) of group r, in mHa^2 per (7-bit step = pi/64)^2:
  H_kk  raw Gauss-Newton diagonal of error_2^2 (everything else held);
  S_kk  the group's weights re-solved (Schur complement on the weights; all weights free, bounds ignored);
  C_kk  = 1 / (S^-1)_kk: the group's weights AND every other angle of the group re-solved (the conditional price:
          what moving that one word by one step costs when all the rest may compensate, to second order).
The other groups are held in all three (upper readings). The model is tools/r48_chem_q7c_linmodel.py's (copied).

Ties (quadratic, per row, no cross-row compensation): the price of forcing theta_e' = f(theta_e) in a row is
  k_e k_e' / (k_e + k_e') * sin^2(theta_e' - f(theta_e))        (the stiffer one keeps its value)
for f in EQ (same word), NEG (complement: -theta), MIR (low 5 bits flipped: sgn pi/2 - theta), NMR (two top bits
flipped: theta - sgn pi/2); k = C_kk (lower) and S_kk (upper). Summed over the rows -> (75, 75, 4) per group.

    python r49_chem_free_spec.py POINT.npz OUT.npz --edges TABLE_or_EDGES.npz [--groups 0,1,..]
"""
import argparse
import os
import sys
import time

import numpy as np
import scipy.linalg as sl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_round import N, Target, angles_and_mass  # noqa: E402
from r48_chem_q7c_absorb import tangents  # noqa: E402

STEP = np.pi / 64.0


def group_spec(r, U, W, A, I, J, s, edges, want_full=True):
    R, B, C = W.shape
    u, V = tangents(U[r], edges)
    T_ = V.shape[1]
    p = (u[:, I] * u[:, J]) * s
    q = ((V[:, :, I] * u[:, None, J] + u[:, None, I] * V[:, :, J]) * s).reshape(B * T_, -1)
    g = (W[r] @ A[r]) * s
    ar = A[r] * s
    bk = np.repeat(np.arange(B), T_)
    Gqq = q @ q.T
    Gqg = q @ g.T
    Ggg = g @ g.T
    Gqp, Gqa = q @ p.T, q @ ar.T
    Ggp, Gga = g @ p.T, g @ ar.T
    Gpp, Gaa, Gap = p @ p.T, ar @ ar.T, ar @ p.T
    na, nw = B * T_, B * C
    X = Gqg[:, bk]
    Haa = 2.0 * (Ggg[np.ix_(bk, bk)] * Gqq + X.T * X)
    del Gqq, X
    wb, wc = np.repeat(np.arange(B), C), np.tile(np.arange(C), B)
    Haw = 2.0 * (Gga[np.ix_(bk, wc)] * Gqp[:, wb] + Ggp[np.ix_(bk, wb)] * Gqa[:, wc])
    Hww = 2.0 * (Gaa[np.ix_(wc, wc)] * Gpp[np.ix_(wb, wb)] + Gap[np.ix_(wc, wb)] * Gap[np.ix_(wc, wb)].T)
    lam = 1e-10 * np.trace(Hww) / nw
    cw = sl.cho_factor(Hww + lam * np.eye(nw))
    Hd = np.diag(Haa).copy()
    S = Haa - Haw @ sl.cho_solve(cw, Haw.T)
    del Haa, Haw, Hww
    S = 0.5 * (S + S.T)
    Sd = np.maximum(np.diag(S).copy(), 0.0)
    ridge = 1e-9 * np.trace(S) / na
    P = sl.cho_solve(sl.cho_factor(S + ridge * np.eye(na)), np.eye(na))
    Cd = 1.0 / np.maximum(np.diag(P), 1e-300)
    # the same with only the row's own other angles free (weights + own row): per-row block of S inverted
    Rd = np.zeros(na)
    for b in range(B):
        sl_ = slice(b * T_, (b + 1) * T_)
        Pb = np.linalg.inv(S[sl_, sl_] + ridge * np.eye(T_))
        Rd[sl_] = 1.0 / np.maximum(np.diag(Pb), 1e-300)
    sc = STEP * STEP * 1e6
    return (Hd.reshape(B, T_) * sc, Sd.reshape(B, T_) * sc, Cd.reshape(B, T_) * sc, Rd.reshape(B, T_) * sc)


def tie_prices(th, k):
    """th (B, 75) radians in the signed range, k (B, 75) stiffness per step^2 -> (75, 75, 4) summed over rows, in the
    units of k times (radians / step)^2 (i.e. mHa^2 when k is mHa^2 per step^2)."""
    B, T_ = th.shape
    a = th[:, :, None]          # theta_e
    b = th[:, None, :]          # theta_e'
    sg = np.where(a >= 0, 1.0, -1.0)
    forms = (a, -a, sg * np.pi / 2 - a, a - sg * np.pi / 2)
    ke, kf = k[:, :, None], k[:, None, :]
    h = ke * kf / np.maximum(ke + kf, 1e-300)
    out = np.zeros((T_, T_, 4))
    for f, fa in enumerate(forms):
        out[:, :, f] = np.sum(h * np.sin(b - fa) ** 2, axis=0) / (STEP * STEP)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("init")
    ap.add_argument("out")
    ap.add_argument("--edges", required=True, help="npz with tree_edges (a table) or edges")
    ap.add_argument("--groups", default="")
    a = ap.parse_args()
    t0 = time.time()
    z = np.load(a.init)
    U, W = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
    R, B, C = W.shape
    U /= np.linalg.norm(U, axis=2, keepdims=True)
    ez = np.load(a.edges)
    edges = np.array(ez["tree_edges"] if "tree_edges" in ez.files else ez["edges"], np.int64)
    tgt = Target()
    I, J = tgt.i, tgt.j
    s = np.sqrt(np.where(I == J, 1.0, 2.0))
    A = np.einsum("rbc,rbp->rcp", W, U[:, :, I] * U[:, :, J])
    groups = [int(x) for x in a.groups.split(",")] if a.groups else list(range(R))
    T_ = N - 1
    Hd = np.zeros((R, B, T_))
    Sd, Cd, Rd, TH, MS = (np.zeros((R, B, T_)) for _ in range(5))
    tieC = np.zeros((R, T_, T_, 4))
    tieS = np.zeros((R, T_, T_, 4))
    for r in groups:
        Hd[r], Sd[r], Cd[r], Rd[r] = group_spec(r, U, W, A, I, J, s, edges)
        th, mass = angles_and_mass(U[r], edges)
        TH[r], MS[r] = th * np.pi, mass
        tieC[r] = tie_prices(TH[r], Cd[r])
        tieS[r] = tie_prices(TH[r], Sd[r])
        print(f"group {r}: sum H {Hd[r].sum():.1f} S {Sd[r].sum():.1f} R {Rd[r].sum():.1f} C {Cd[r].sum():.1f} mHa^2 per step^2 "
              f"({time.time() - t0:.0f} s)", flush=True)
        np.savez(a.out, H=Hd, S=Sd, C=Cd, Rw=Rd, theta=TH, mass=MS, tieC=tieC, tieS=tieS, edges=edges,
                 wn2=np.sum(W ** 2, axis=2), groups=np.array(groups))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
