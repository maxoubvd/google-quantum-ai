"""r48_chem_q7c_absorb: idealised reading for the TREE as a lever, by what the weights cannot absorb (chem_own2).

A one-step kick d on edge t of row b (group r) moves the row's pair vector p_b = (u_i u_j) by
m_bt = d (u_i v_j + v_i u_j), v = du/dtheta_t, in every copy c with the factor w_rbc. The weights of group r can
re-express exactly the part of m_bt that lies in span{p_b' : b' in group r} (the M-weighted pair space), whatever
the copy. So
    raw  = sum_b |w_b|^2 sum_t |m_bt|^2            (the raw stiffness of all angle kicks)
    hard = sum_b |w_b|^2 sum_t |(1 - P_r) m_bt|^2  (what only other ANGLES can compensate)
are computed for a tree; hard / raw is the share of the kick energy that needs angle compensation. Trees: the
log-depth tree of chem_own_round*.py on the natural orbital labels, and the same tree shape with the labels
1..75 permuted (orbital 0 stays the root), random or chosen by --perm.

    python r48_chem_q7c_absorb.py INIT.npz [--random K] [--groups 0,1] [--perm PERM.npy]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_round import N, angles_and_mass, decode_np, log_depth_tree  # noqa: E402

I, J = np.tril_indices(N)
SM = np.sqrt(np.where(I == J, 1.0, 2.0))


def tangents(U, edges):
    """v[b, t] = d u_b / d theta_t (exact central difference of the decode at a small step)."""
    th, _ = angles_and_mass(U, edges)
    h = 1e-6
    B, T = th.shape
    V = np.zeros((B, T, N))
    for t in range(T):
        tp, tm = th.copy(), th.copy()
        tp[:, t] += h / np.pi
        tm[:, t] -= h / np.pi
        V[:, t] = (decode_np(tp, edges) - decode_np(tm, edges)) / (2 * h)
    u = decode_np(th, edges)
    sgn = np.sign(np.sum(u * U, axis=1))          # the decode reproduces each row up to a global sign
    return u * sgn[:, None], V * sgn[:, None, None]


def group_reading(Ug, wn2, edges):
    u, V = tangents(Ug, edges)
    P = (u[:, I] * u[:, J]) * SM[None, :]                      # (B, pairs)
    Q, _ = np.linalg.qr(P.T)                                   # orthonormal basis of the group's span
    M = (u[:, None, I] * V[:, :, J] + V[:, :, I] * u[:, None, J]) * SM[None, None, :]     # (B, T, pairs)
    tot = np.sum(M * M, axis=2)
    Bn, T, Pn = M.shape
    inside = np.sum((M.reshape(Bn * T, Pn) @ Q) ** 2, axis=1).reshape(Bn, T)
    hard = np.maximum(tot - inside, 0.0)
    return float(np.sum(wn2[:, None] * tot)), float(np.sum(wn2[:, None] * hard)), hard, tot


def relabel(edges, perm):
    """perm[old label] = new label (perm[0] = 0)."""
    return perm[edges]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("init")
    ap.add_argument("--random", type=int, default=2)
    ap.add_argument("--groups", default="")
    ap.add_argument("--perm", default=None)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    z = np.load(a.init)
    U, W = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
    R, B, C = W.shape
    U /= np.linalg.norm(U, axis=2, keepdims=True)
    wn2 = np.sum(W ** 2, axis=2)
    groups = [int(x) for x in a.groups.split(",")] if a.groups else list(range(R))
    base = log_depth_tree(N)
    trees = [("log-depth, natural labels", base)]
    rng = np.random.default_rng(a.seed)
    for k in range(a.random):
        perm = np.concatenate([[0], 1 + rng.permutation(N - 1)])
        trees.append((f"log-depth, random labels {k}", relabel(base, perm)))
    if a.perm:
        trees.append(("log-depth, labels from " + os.path.basename(a.perm), relabel(base, np.load(a.perm))))
    for name, E in trees:
        raw = hard = 0.0
        for r in groups:
            r_raw, r_hard, _, _ = group_reading(U[r], wn2[r], E)
            raw += r_raw
            hard += r_hard
        print(f"{name}: groups {groups}: raw {raw:.4f}, needs angle compensation {hard:.4f} ({hard / raw:.3f})", flush=True)


if __name__ == "__main__":
    main()
