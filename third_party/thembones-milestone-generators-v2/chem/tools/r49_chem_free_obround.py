"""r49_chem_free_obround: the ONE-BODY rows rounded to 7 bits UNDER THE BIT-0 CODE, with refits (R49-CHEM-FREE).

A copy of the progressive rounding of tools/r48_chem_q7c_onebody.py (26 stages, impact order lam^2 x mass x offset^2,
an L-BFGS refit of the free angles on error_1^2 after every stage, then exact coordinate descent), with the loaders'
bit-0 code on the columns --cols: per 16-row block (rows 16p..16p+15 of the 76) and column, bit 0 of the words is an
affine function of the four low address bits. In the staged order a word's bit 0 is free until its low address lies in
the affine span of the lows already frozen in its (block, column); afterwards it is forced, and the word goes to the
nearest value of that parity. --weighted: a free word takes the parity of the best codeword for its whole set (weights
lam^2 x mass, penalties from the present refitted angles). The polish moves coded words by +-2 (parity kept).
The two-body part of the table (h1') is held. Rows 76-79 (two-body, block 4) are not constrained here.

    python r49_chem_free_obround.py TABLE.npz [--cols 15-74] [--weighted] [--out OUT.npz]
"""
import argparse
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r48_chem_q7c_onebody as OB  # noqa: E402
from r48_chem_q7c_onebody import (STAGES, angles_and_mass, check_fcidump, decode_jac, decode_pi, error1,  # noqa: E402
                                  load_fixed_shifts, load_reference, reconstruct_u_vectors, symmetric_from_lower)

CW = np.array([[(a0 ^ (bin(a & l).count("1") & 1)) for l in range(16)] for a0 in (0, 1) for a in range(16)], np.int64)


def near_par(x, p):
    return 2.0 * np.round((x - p) / 2.0) + p


def coset_vec(l, cidx):
    """Equation of the word at low address l, column number cidx of its slot, in the unknowns (v_0..v_15, then five
    affine coefficients per column): parity = v_l + a0 + a.l (C2's weaker requirement: the columns of one slot in ONE
    coset of RM(1,4) per block; the XOR of two neighbouring columns is then affine)."""
    v = 1 << l
    for i, m in enumerate((0, 1, 2, 4, 8)):
        if (l & m) == m:
            v |= 1 << (16 + 5 * cidx + i)
    return v


def progressive_coded(th0, edges, lam_at, H, scale, iters, coded, weighted, coset=False):
    theta = th0.copy()
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros(theta.shape)
    lam2 = lam_at[:, None] ** 2
    nrow = theta.shape[0]
    fb = {}                                    # (block, col) -> {low: bit}
    nforced = nmoved = 0

    def fun(x, mask, base):
        t = base.copy()
        t[mask] = x
        V, J = decode_jac(t, edges)
        D = (V.T * lam_at) @ V - H
        gV = 4.0 * lam_at[:, None] * (V @ D)
        return float(np.sum(D * D)), np.einsum("jn,jkn->jk", gV, J)[mask]

    def consistent(d):
        ok = np.ones(len(CW), bool)
        for l_, b_ in d.items():
            ok &= CW[:, l_] == b_
        return ok

    for frac in STAGES:
        need = int(round(frac * theta.size)) - int(frozen.sum())
        if need > 0:
            now = np.where(frozen, qv, theta)
            mass = angles_and_mass(decode_pi(now, edges), edges)[1]
            q = np.round(now * scale) / scale
            impact = lam2 * mass * (q - now) ** 2
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            wgt = lam2 * mass
            for fl in pick:
                j, t = divmod(int(fl), theta.shape[1])
                xs = now[j, t] * scale
                z = np.round(xs)
                if t in coded and coset:
                    key, l = (j // 16, t % 5), j % 16
                    bas = fb.setdefault(key, [])
                    v, val = coset_vec(l, (t - min(coded)) // 5), 0
                    for bv, bval in bas:
                        if v & (bv & -bv):
                            v ^= bv
                            val ^= bval
                    if v == 0:
                        z = near_par(xs, val)
                        nforced += 1
                        nmoved += int(z != np.round(xs))
                    else:
                        bas.append((v, val ^ (int(z) & 1)))
                elif t in coded:
                    key, l = (j // 16, t), j % 16
                    d = fb.setdefault(key, {})
                    ok = consistent(d)
                    vals = np.unique(CW[ok, l])
                    if len(vals) == 1:                       # forced
                        z = near_par(xs, int(vals[0]))
                        nforced += 1
                        nmoved += int(z != np.round(xs))
                    elif weighted:
                        cost = np.zeros(len(CW))
                        for l2 in range(16):
                            j2 = 16 * key[0] + l2
                            if j2 >= nrow or l2 in d:
                                continue
                            x2 = now[j2, t] * scale
                            pen = [(near_par(x2, par) - x2) ** 2 for par in (0, 1)]
                            cost += wgt[j2, t] * np.where(CW[:, l2] == 0, pen[0], pen[1])
                        cost[~ok] = np.inf
                        z = near_par(xs, int(CW[int(np.argmin(cost)), l]))
                    d[l] = int(z) & 1
                frozen[j, t] = True
                qv[j, t] = z / scale
        free = ~frozen
        if free.any():
            base = np.where(frozen, qv, theta)
            r = minimize(fun, theta[free], args=(free, base), jac=True, method="L-BFGS-B",
                         options={"maxiter": iters, "maxcor": 20, "ftol": 1e-14, "gtol": 1e-10})
            theta = base.copy()
            theta[free] = r.x
        else:
            theta = qv.copy()
    return np.round(qv * scale).astype(np.int64), nforced, nmoved


def progressive_parity(th0, edges, lam_at, H, scale, iters, par):
    """The progressive rounding with the parity of every word prescribed (par >= 0) or free (par = -1)."""
    theta = th0.copy()
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros(theta.shape)
    lam2 = lam_at[:, None] ** 2

    def fun(x, mask, base):
        t = base.copy()
        t[mask] = x
        V, J = decode_jac(t, edges)
        D = (V.T * lam_at) @ V - H
        gV = 4.0 * lam_at[:, None] * (V @ D)
        return float(np.sum(D * D)), np.einsum("jn,jkn->jk", gV, J)[mask]

    for frac in STAGES:
        need = int(round(frac * theta.size)) - int(frozen.sum())
        if need > 0:
            now = np.where(frozen, qv, theta)
            mass = angles_and_mass(decode_pi(now, edges), edges)[1]
            xs = now * scale
            q = np.where(par >= 0, near_par(xs, np.maximum(par, 0)), np.round(xs)) / scale
            impact = lam2 * mass * (q - now) ** 2
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            frozen.flat[pick] = True
            qv.flat[pick] = q.flat[pick]
        free = ~frozen
        if free.any():
            base = np.where(frozen, qv, theta)
            r = minimize(fun, theta[free], args=(free, base), jac=True, method="L-BFGS-B",
                         options={"maxiter": iters, "maxcor": 20, "ftol": 1e-14, "gtol": 1e-10})
            theta = base.copy()
            theta[free] = r.x
        else:
            theta = qv.copy()
    return np.round(qv * scale).astype(np.int64)


def two_pass(th0, edges, lam_at, H, iters, coded, use_v):
    """Pass A: uncoded progressive rounding (every word at its own nearest value, with refits). The parities of the
    five rows at the low addresses 0, 1, 2, 4, 8 of each block then DEFINE each column's affine function; every other
    row takes the forced parity, with (use_v) one free bit per row and slot (the coset form): the weighted majority of
    its own preferences over the columns of the slot (weights lam^2 x mass x (1 - 2 |offset|) at the exact angles).
    Pass B: the progressive rounding with all those parities prescribed."""
    nrow, ncol = th0.shape
    wa = progressive_parity(th0, edges, lam_at, H, 64.0, iters, np.full(th0.shape, -1, np.int64))
    pa = wa & 1
    xs = th0 * 64.0
    mass = angles_and_mass(decode_pi(th0, edges), edges)[1]
    pref = np.round(xs).astype(np.int64) & 1
    wgt = lam_at[:, None] ** 2 * mass * (1.0 - 2.0 * np.abs(xs - np.round(xs)))
    par = np.full(th0.shape, -1, np.int64)
    cols = sorted(coded)
    wrong = tot = 0.0
    for blk in range(5):
        rows = [j for j in range(16 * blk, min(16 * blk + 16, nrow))]
        aff = {}
        for t in cols:
            q0 = pa[16 * blk, t]
            aff[t] = {j: int(q0 ^ np.bitwise_xor.reduce([pa[16 * blk + (1 << i), t] ^ q0 for i in range(4) if (j % 16) >> i & 1] or [0]))
                      for j in rows}
        for j in rows:
            heavy = (j % 16) in (0, 1, 2, 4, 8)
            for slot in range(5):
                ts = [t for t in cols if t % 5 == slot]
                if heavy:
                    for t in ts:
                        par[j, t] = pa[j, t]
                    continue
                f = np.array([aff[t][j] for t in ts])
                mis = (pref[j, ts] != f).astype(np.float64) * wgt[j, ts]
                v = 0
                if use_v and mis.sum() > 0.5 * wgt[j, ts].sum():
                    v = 1
                for t, ft in zip(ts, f):
                    par[j, t] = ft ^ v
                wrong += float(((pref[j, ts] != (f ^ v)) * wgt[j, ts]).sum())
                tot += float(wgt[j, ts].sum())
    wb = progressive_parity(th0, edges, lam_at, H, 64.0, iters, par)
    return wb, wrong / max(tot, 1e-300)


def polish_coded(wq, edges, H, brot, sweeps, coded, inner=25):
    """r48's exact coordinate descent on error_1 with the step of a coded column = 2 (bit 0 kept)."""
    scale, mod = float(1 << (brot - 1)), 1 << brot
    N = edges.shape[0] + 1
    ncol = wq.shape[1]
    wq = np.mod(np.array(wq, np.int64), mod)
    V = decode_pi(wq / scale, edges)
    e1, perm, D, lam = error1(V, H)
    ar = np.arange(ncol)
    step = np.array([2 if t in coded else 1 for t in range(ncol)], np.int64)
    for sweep in range(sweeps):
        moved = 0
        for j in range(N):
            lj = lam[perm[j]]
            base = D - lj * np.outer(V[j], V[j])
            for _ in range(inner):
                v0 = V[j]
                cand = np.repeat(wq[j][None, :], 2 * ncol, axis=0)
                cand[ar, ar] += step
                cand[ncol + ar, ar] -= step
                cand = np.mod(cand, mod)
                Vc = decode_pi(cand / scale, edges)
                quad = np.einsum("ip,pq,iq->i", Vc, base, Vc)
                n4 = np.sum(Vc * Vc, axis=1) ** 2
                val = 2 * lj * quad + lj * lj * n4
                cur = 2 * lj * float(v0 @ base @ v0) + lj * lj * float(v0 @ v0) ** 2
                i = int(np.argmin(val))
                if val[i] >= cur - 1e-13:
                    break
                wq[j] = cand[i]
                V[j] = Vc[i]
                moved += 1
            D = base + lj * np.outer(V[j], V[j])
        e_prev = e1
        e1, perm, D, lam = error1(V, H)
        if not moved or e_prev - e1 < 1e-7 * e1:
            break
    return wq, e1


def check_coset(wq, coded):
    bad = 0
    for blk in range(5):
        rows = np.arange(16 * blk, min(16 * blk + 16, wq.shape[0]))
        for t in coded:
            if t + 5 in coded:
                bits = (wq[rows, t] ^ wq[rows, t + 5]) & 1
                if not (CW[:, :len(rows)] == bits[None, :]).all(axis=1).any():
                    bad += 1
    return bad


def check_code(wq, coded):
    bad = 0
    for blk in range(5):
        rows = np.arange(16 * blk, min(16 * blk + 16, wq.shape[0]))
        for t in coded:
            bits = wq[rows, t] & 1
            if not (CW[:, :len(rows)] == bits[None, :]).all(axis=1).any():
                bad += 1
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--cols", default="15-74")
    ap.add_argument("--weighted", action="store_true")
    ap.add_argument("--iters", type=int, default=300)
    ap.add_argument("--sweeps", type=int, default=40)
    ap.add_argument("--out", default="")
    ap.add_argument("--order", default="table", help="table | heavy5: the 25 rows of largest |lam| at the low addresses "
                    "0, 1, 2, 4, 8 of the five blocks (affinely independent: the code can give all of them their own "
                    "parity in every column), the others by |lam| over the remaining addresses")
    ap.add_argument("--only-weighted", action="store_true")
    ap.add_argument("--coset", action="store_true", help="the weaker code: the columns of a slot in one coset per block")
    ap.add_argument("--two-pass", action="store_true", help="heavy rows first (their parities define the code), then the rest forced")
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
    th0, mass0 = angles_and_mass(vec[:, perm].T, edges)          # rows in the table's one-body order
    lam_at = lam_e[perm]
    l2 = np.sort(lam_at ** 2)[::-1]
    print("  lam^2 shares: top 5 / 10 / 25 / 40 rows carry " + " / ".join(f"{100 * l2[:k_].sum() / l2.sum():.1f}%" for k_ in (5, 10, 25, 40))
          + f"; |lam| max {np.abs(lam_at).max():.3f} median {np.median(np.abs(lam_at)):.3f} min {np.abs(lam_at).min():.4f} Ha", flush=True)
    lo, hi = (int(v) for v in a.cols.split("-"))
    if a.order == "assign":
        # heavy rows as in heavy5; every light row goes to the address whose FORCED parities (the affine extension of
        # its block's five heavy rows, taken from their unconstrained nearest words) cost it least:
        # sum over the coded columns of lam^2 x mass x (1 - 2 |offset|) where its own nearest parity differs.
        from scipy.optimize import linear_sum_assignment as lsa
        rank = np.argsort(-np.abs(lam_at), kind="stable")
        ind = [16 * b_ + l_ for l_ in (0, 1, 2, 4, 8) for b_ in range(5)]
        heavy = rank[:25]
        row_at = np.full(N, -1, np.int64)
        row_at[np.array(ind)] = heavy
        xs = th0 * 64.0
        par = (np.round(xs).astype(np.int64) & 1)[:, lo:hi + 1]
        pen = (lam_at[:, None] ** 2 * mass0 * (1.0 - 2.0 * np.abs(xs - np.round(xs))))[:, lo:hi + 1]
        slots = [a_ for a_ in range(N) if a_ not in ind]
        light = rank[25:]
        cost = np.zeros((len(light), len(slots)))
        for si, a_ in enumerate(slots):
            b_, l_ = divmod(a_, 16)
            q0 = par[row_at[16 * b_]]
            f = q0.copy()
            for i_ in range(4):
                if l_ >> i_ & 1:
                    f = f ^ par[row_at[16 * b_ + (1 << i_)]] ^ q0
            cost[:, si] = ((par[light] != f[None, :]) * pen[light]).sum(axis=1)
        ri, ci = lsa(cost)
        base_cost = float(np.mean(cost) * len(light))
        for r_, c_ in zip(ri, ci):
            row_at[slots[c_]] = light[r_]
        print(f"  assign: light rows placed by forced-parity cost: {cost[ri, ci].sum():.4e} against {base_cost:.4e} for a random placement", flush=True)
        th0, lam_at = th0[row_at], lam_at[row_at]
    if a.order == "heavy5":
        rank = np.argsort(-np.abs(lam_at), kind="stable")
        ind = [16 * b_ + l_ for l_ in (0, 1, 2, 4, 8) for b_ in range(5)]
        rest = [a_ for a_ in range(N) if a_ not in ind]
        # light rows: spread by |lam| over the blocks too (round robin over the remaining addresses, block by block)
        rest.sort(key=lambda a_: (a_ % 16, a_ // 16))
        slots = ind + rest
        row_at = np.empty(N, np.int64)                           # row_at[address] = old row index
        row_at[np.array(slots)] = rank
        th0, lam_at = th0[row_at], lam_at[row_at]
    variants = [("uncoded", set(), False, False), (f"bit-0 code on columns {a.cols}, first-come", set(range(lo, hi + 1)), False, False),
                (f"bit-0 code on columns {a.cols}, weighted", set(range(lo, hi + 1)), True, False)]
    if a.coset:
        variants = [(f"bit-0 COSET code (one coset per slot and block) on columns {a.cols}, first-come", set(range(lo, hi + 1)), False, True)]
    if a.two_pass:
        variants = []
        for use_v in (False, True):
            coded = set(range(lo, hi + 1))
            wq, share = two_pass(th0, edges, lam_at, H, a.iters, coded, use_v)
            e_p = error1(decode_pi(wq / 64.0, edges), H)[0]
            wq2, e_f = polish_coded(wq, edges, H, 7, a.sweeps, coded)
            bad = check_coset(wq2, coded) if use_v else check_code(wq2, coded)
            print(f"  two-pass, {'coset form (one free bit per row and slot)' if use_v else 'every column a codeword'}, columns {a.cols}: weighted share of "
                  f"forced-wrong parities among the light rows {100 * share:.1f}%; progressive error_1 {1e3 * e_p:.4f} mHa, after the polish "
                  f"{1e3 * e_f:.4f} mHa; code violations {bad} ({time.time() - t0:.0f} s)", flush=True)
            if a.out:
                out = words.copy()
                out[:N] = np.mod(wq2, 128)
                keep = {k_: z[k_] for k_ in z.files if k_ != "angle_words"}
                np.savez(a.out[:-4] + ("_v.npz" if use_v else "_c.npz"), angle_words=out, **keep)
    for name, coded, wtd, cst in variants:
        if name.endswith("weighted") and not a.weighted:
            continue
        if a.only_weighted and not name.endswith("weighted") and not cst:
            continue
        wq, nf, nm = progressive_coded(th0, edges, lam_at, H, 64.0, a.iters, coded, wtd, coset=cst)
        V = decode_pi(wq / 64.0, edges)
        e_p = error1(V, H)[0]
        wq2, e_f = polish_coded(wq, edges, H, 7, a.sweeps, coded)
        bad = check_coset(wq2, coded) if cst else check_code(wq2, coded)
        print(f"  {name}: progressive error_1 {1e3 * e_p:.4f} mHa, after the polish {1e3 * e_f:.4f} mHa; forced words {nf}, moved off nearest {nm}; "
              f"code violations {bad} ({time.time() - t0:.0f} s)", flush=True)
        if a.out and coded:
            out = words.copy()
            out[:N] = np.mod(wq2, 128)
            keep = {k_: z[k_] for k_ in z.files if k_ != "angle_words"}
            np.savez(a.out if not wtd else a.out[:-4] + "_w.npz", angle_words=out, **keep)


if __name__ == "__main__":
    main()
