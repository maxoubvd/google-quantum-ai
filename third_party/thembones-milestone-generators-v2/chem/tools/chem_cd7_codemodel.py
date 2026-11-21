"""chem_cd7_codemodel: CD7 kill test (a). The coset form of MY rounding tool (tools/chem_cd7_round.py: shared code
d == 11 handled by its residues(), the tool's exact weighted decoder: per (block, slot) the cheapest common syndrome
of RM(1,4) given the frozen words) inside the designer's linear staged model (tools/r49_chem_free_codemodel.py:
S = the group's Gauss-Newton matrix with its weights eliminated, staged freezing in S_kk order, exact Gaussian
conditioning of the free words). For each coded range b..14 it prints the staged damage of
  (1) no code, (2) the designer's CosetCoder (first-come), (3) my tool's residues with weights S_kk,
  (4) my tool's residues with weights 1 / P_jj (the present conditional stiffness).

    python chem_cd7_codemodel.py POINT.npz --edges TABLE.npz --perm cost_perm.npy --group 2 --batches 7,5,3
"""
import argparse
import os
import sys
import time

import numpy as np
import scipy.linalg as sl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_cd7_round as CR  # noqa: E402
from r49_chem_free_codemodel import CosetCoder, build_S, staged_coded  # noqa: E402
from r48_chem_q7c_round import LONG, N, Target, angles_and_mass  # noqa: E402


def staged_tool(S, x, fracs, sk, deg, xdeg, blk_row, pos_row, nrows, ncols, wmode):
    """The staged procedure of staged_coded with the picked words set by the tool's residues() (all words at their
    present means, frozen words with their values)."""
    n = len(x)
    P = sl.cho_solve(sl.cho_factor(S + 1e-9 * np.trace(S) / n * np.eye(n)), np.eye(n))
    mu = x.copy()
    free = np.arange(n)
    frozen = np.zeros(n, bool)
    val = np.zeros(n, np.int64)
    E, done, nflip = 0.0, 0, 0
    cP = np.zeros(n)
    for f in fracs:
        need = int(round(f * n)) - done
        if need <= 0:
            continue
        pick = np.argsort(-sk[free], kind="stable")[:need]
        keep = np.setdiff1d(np.arange(len(free)), pick)
        t = np.where(frozen, val.astype(np.float64), mu)
        if wmode == "S":
            c = sk
        else:
            cP[:] = 0.0
            cP[free] = 1.0 / np.maximum(np.diag(P), 1e-300)
            c = np.where(frozen, 1.0, cP)
        r_, m_ = CR.residues(t.reshape(nrows, ncols), c.reshape(nrows, ncols) + 1e-30, blk_row, pos_row, deg, xdeg,
                             frozen.reshape(nrows, ncols), val.reshape(nrows, ncols))
        allowed = CR.into_range(CR.residue_round(t.reshape(nrows, ncols), r_, m_), m_).reshape(-1)
        k_ = free[pick]
        zt = allowed[k_]
        nflip += int(np.sum(zt != np.round(mu[k_])))
        e = zt - mu[k_]
        Pdd = P[np.ix_(pick, pick)]
        cf = sl.cho_factor(Pdd + 1e-14 * np.trace(Pdd) / len(pick) * np.eye(len(pick)))
        sol = sl.cho_solve(cf, e)
        E += float(e @ sol)
        Pud = P[np.ix_(keep, pick)]
        mu[free[keep]] += Pud @ sol
        P = P[np.ix_(keep, keep)] - Pud @ sl.cho_solve(cf, Pud.T)
        frozen[k_] = True
        val[k_] = zt.astype(np.int64)
        mu[k_] = zt
        free = free[keep]
        done += len(pick)
    # the code must hold on the finished words (only this group's rows; other rows of a block are absent = free)
    return E, nflip, val


def coset_ok(val, nrows, ncols, blk_row, pos_row, cols):
    """True if on every (block, slot) the bit-0 patterns of the coded columns share one syndrome of RM(1,4) on the
    rows present (checked as: XOR of two columns of a slot is affine on the present positions)."""
    v = val.reshape(nrows, ncols)
    CW = np.array([[(a0 ^ (bin(a & l).count("1") & 1)) for l in range(16)] for a0 in (0, 1) for a in range(16)])
    bad = 0
    for b in np.unique(blk_row):
        rr = np.flatnonzero(blk_row == b)
        ps = pos_row[rr]
        for j in range(5):
            cs = [c for c in cols if c % 5 == j]
            for c1, c2 in zip(cs[:-1], cs[1:]):
                bits = (v[rr, c1] ^ v[rr, c2]) & 1
                if not (CW[:, ps] == bits[None, :]).all(axis=1).any():
                    bad += 1
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("init")
    ap.add_argument("--edges", required=True)
    ap.add_argument("--perm", required=True)
    ap.add_argument("--group", type=int, default=2)
    ap.add_argument("--batches", default="7,5,3")
    ap.add_argument("--modes", default="first,S,P")
    a = ap.parse_args()
    t0 = time.time()
    CR.HALF = 1 << 20                                  # the model's angles are not gauged into the signed range: no
    #                                                    range limit here (every gauge move keeps bit 0)
    z = np.load(a.init)
    U, W = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
    R, B, C = W.shape
    U /= np.linalg.norm(U, axis=2, keepdims=True)
    edges = np.array(np.load(a.edges)["tree_edges"], np.int64)
    perm = np.load(a.perm)
    tgt = Target()
    I, J = tgt.i, tgt.j
    s = np.sqrt(np.where(I == J, 1.0, 2.0))
    A = np.einsum("rbc,rbp->rcp", W, U[:, :, I] * U[:, :, J])
    r = a.group
    S = build_S(r, U, W, A, I, J, s, edges)
    T_ = N - 1
    th, _ = angles_and_mass(U[r], edges)
    x = (th * 64.0).reshape(-1)
    addr = np.zeros(B, np.int64)
    addr[perm[r]] = 76 + r * B + np.arange(B)
    block_of = np.repeat(addr // 16, T_)
    low_of = np.repeat(addr % 16, T_)
    col_of = np.tile(np.arange(T_), B)
    sk = np.diag(S).copy()
    blk_row = addr // 16
    blk_row = blk_row - blk_row.min()                  # the group's own blocks, 0..k-1 (speed)
    pos_row = addr % 16
    print(f"group {r}: S built ({time.time() - t0:.0f} s); {B * T_} words; blocks {blk_row.max() + 1}", flush=True)
    base, _, _ = staged_tool(S, x, LONG, sk, CR.parse_codes(""), CR.parse_xcodes(""), blk_row, pos_row, B, T_, "S")
    E0, _, _ = staged_coded(S, x, LONG, sk, None, block_of, low_of, col_of)
    print(f"  none: tool path {base:.2f}, designer path {E0:.2f} mHa^2 ({time.time() - t0:.0f} s)", flush=True)
    for bb in (int(v) for v in a.batches.split(",")):
        cols = list(range(5 * bb, 75))
        for mode in a.modes.split(","):
            if mode == "first":
                E, nf, nfl = staged_coded(S, x, LONG, sk, CosetCoder(set(cols)), block_of, low_of, col_of)
                tag = "designer CosetCoder, first-come"
                bad = -1
            else:
                E, nfl, val = staged_tool(S, x, LONG, sk, CR.parse_codes(f"0:11:{bb}-14"), CR.parse_xcodes(""),
                                          blk_row, pos_row, B, T_, mode)
                bad = coset_ok(val, B, T_, blk_row, pos_row, cols)
                tag = f"MY TOOL residues d=11, weights {'S_kk' if mode == 'S' else '1/P_jj'}"
            print(f"  batches {bb}-14 [{tag}]: staged damage {E:.2f} mHa^2 ({E - E0:+.2f}, x{E / E0:.3f}); moved off nearest "
                  f"{nfl}; coset violations {bad} ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
