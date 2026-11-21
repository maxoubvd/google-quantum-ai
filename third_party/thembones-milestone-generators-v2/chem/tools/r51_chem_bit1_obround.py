"""r51_chem_bit1_obround: the ONE-BODY rows rounded to 7 bits under a code LAYOUT (any bits, chains of columns in
one coset of RM(1,4) per block), with refits (R51-CHEM-BIT1). A copy of tools/r49_chem_free_obround.py (R49-CHEM-FREE:
the progressive rounding of tools/r48_chem_q7c_onebody.py, 26 stages, impact order, an L-BFGS refit of the free angles
on error_1^2 after every stage, then exact coordinate descent) with the code generalised: the layout is the spec of
tools/r51_chem_bit1_codemodel.py (first-come: a word's coded bit is free while its equation is independent of those
frozen in its (block, chain)); a word with forced bits goes to the nearest integer carrying them; the polish moves a
word by 2^(coded bits of its column), which keeps every coded bit. h1' (the two-body part of the table) is held.
Rows 76-79 (two-body, block 4) are not constrained here.

    python r51_chem_bit1_obround.py TABLE.npz --spec "0:1:B5-14,1:1:G2@7-14" [--spec2 ...] [--order heavy5]
"""
import argparse
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_onebody import (STAGES, angles_and_mass, check_fcidump, decode_jac, decode_pi, error1,  # noqa: E402
                                  load_fixed_shifts, load_reference, reconstruct_u_vectors, symmetric_from_lower)
from r51_chem_bit1_codemodel import GenCoder, parse_spec  # noqa: E402
import r51_chem_bit1_round as RD  # noqa: E402

RD.HALF = 64                                   # no wall here: the rows' words use the full circle (-64..63); the signed
                                               # range is reached afterwards by the exact gauge move (obsigned)

CW = np.array([[(a0 ^ (bin(a & l).count("1") & 1)) for l in range(16)] for a0 in (0, 1) for a in range(16)], np.int64)


def progressive_layout(th0, edges, lam_at, H, scale, iters, coder, tool=None):
    """tool = (deg, xdeg): the builder's own decoder (tools/r51_chem_bit1_round.py residues with the layout's chains:
    per stage the cheapest residues of every unfrozen word given the frozen ones, weights lam^2 x mass) instead of
    first-come."""
    theta = th0.copy()
    val = np.zeros(theta.shape, np.int64)
    rows_ = np.arange(theta.shape[0])
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros(theta.shape)
    lam2 = lam_at[:, None] ** 2
    nforced = nmoved = 0

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
            q = np.round(now * scale) / scale
            impact = lam2 * mass * (q - now) ** 2
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            if tool is not None:
                r_, m_ = RD.residues(now * scale, lam2 * mass + 1e-30, rows_ // 16, rows_ % 16, tool[0], tool[1], frozen, val)
            for fl in pick:
                j, t = divmod(int(fl), theta.shape[1])
                xs = now[j, t] * scale
                z = int(np.round(xs))
                if tool is not None:
                    z2 = int(RD.residue_round(xs, r_[j, t], m_[t]))
                    nforced += int(m_[t] > 1)
                    nmoved += int(z2 != z)
                    z = z2
                    val[j, t] = z
                elif coder is not None and t in coder.where:
                    info = coder.forced(j // 16, t, j % 16)
                    cons = {bit: fv for _, bit, fv, _, _ in info if fv is not None}
                    if cons:
                        span = 1 << (max(cons) + 1)
                        best = None
                        for dz in range(-span, span + 1):
                            zc = z + dz
                            if all(((zc >> b) & 1) == v for b, v in cons.items()):
                                if best is None or abs(zc - xs) < abs(best - xs):
                                    best = zc
                        nforced += 1
                        nmoved += int(best != z)
                        z = best
                    coder.commit(info, z)
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


def polish(wq, edges, H, brot, sweeps, step, inner=25):
    """r48's exact coordinate descent on error_1; the step of a column = 2^(its coded bits)."""
    scale, mod = float(1 << (brot - 1)), 1 << brot
    N = edges.shape[0] + 1
    ncol = wq.shape[1]
    wq = np.array(wq, np.int64)
    V = decode_pi(wq / scale, edges)
    e1, perm, D, lam = error1(V, H)
    ar = np.arange(ncol)
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


def violations(wq, codes):
    bad = 0
    for bit, deg, chains in codes:
        for ch in chains:
            for c0, c1 in zip(ch[:-1], ch[1:]):
                for blk in range(5):
                    rows = np.arange(16 * blk, min(16 * blk + 16, wq.shape[0]))
                    d = ((wq[rows, c0] ^ wq[rows, c1]) >> bit) & 1
                    bad += int(not (CW[:, :len(rows)] == d[None, :]).all(axis=1).any())
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--spec", action="append", default=[])
    ap.add_argument("--specfile", action="append", default=[])
    ap.add_argument("--iters", type=int, default=300)
    ap.add_argument("--sweeps", type=int, default=40)
    ap.add_argument("--order", default="heavy5")
    ap.add_argument("--out", default="")
    ap.add_argument("--dumpw", default="")
    ap.add_argument("--decoder", default="tool", help="tool (the rounding tool's exact residues) | first (first-come)")
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
    th0, mass0 = angles_and_mass(vec[:, perm].T, edges)
    lam_at = lam_e[perm]
    if a.order == "heavy5":
        rank = np.argsort(-np.abs(lam_at), kind="stable")
        ind = [16 * b_ + l_ for l_ in (0, 1, 2, 4, 8) for b_ in range(5)]
        rest = [a_ for a_ in range(N) if a_ not in ind]
        rest.sort(key=lambda a_: (a_ % 16, a_ // 16))
        row_at = np.empty(N, np.int64)
        row_at[np.array(ind + rest)] = rank
        th0, lam_at = th0[row_at], lam_at[row_at]
    if a.dumpw:      # the rows' static weights lam^2 x mass in the placed order (for the layout proxy)
        np.savez(a.dumpw, dyn=lam_at[:, None] ** 2 * mass0[row_at] if a.order == "heavy5" else lam_at[:, None] ** 2 * mass0,
                 sk=lam_at[:, None] ** 2 * (mass0[row_at] if a.order == "heavy5" else mass0), addr=np.arange(N), w2=lam_at ** 2,
                 group=-1, damage=1e6 * e_own ** 2)
        return
    specs = list(a.spec) + [open(f).read().strip() for f in a.specfile]
    for si, spec in enumerate(specs):
        _, codes = parse_spec(spec, np.arange(75))
        coder = GenCoder(codes) if codes else None
        tool = None
        if a.decoder == "tool" and codes:
            deg, chain = RD.layout_arrays([[b_, 1, ch] for b_, d_, chains in codes for ch in chains], (75, 8))
            RD.CHAIN = chain
            tool = (deg, np.zeros((15, 8), np.int64))
        wq, nf, nm = progressive_layout(th0, edges, lam_at, H, 64.0, a.iters, coder, tool)
        e_p = error1(decode_pi(wq / 64.0, edges), H)[0]
        step = np.ones(75, np.int64)
        for bit, deg, chains in codes:
            for ch in chains:
                for c in ch:
                    step[c] = max(step[c], 1 << (bit + 1))
        wq2, e_f = polish(wq, edges, H, 7, a.sweeps, step)
        bad = violations(wq2, codes)
        links = sum(len(ch) - 1 for _, _, chains in codes for ch in chains)
        print(f"  [{spec[:90]}{'...' if len(spec) > 90 else ''}] links {links}: progressive error_1 {1e3 * e_p:.4f} mHa, after the polish "
              f"{1e3 * e_f:.4f} mHa; forced words {nf}, moved off nearest {nm}; code violations {bad} ({time.time() - t0:.0f} s)", flush=True)
        if a.out:
            out = words.copy()
            out[:N] = np.mod(wq2, 128)
            keep = {k_: z[k_] for k_ in z.files if k_ != "angle_words"}
            np.savez(a.out[:-4] + f"_{si}.npz", angle_words=out, **keep)


if __name__ == "__main__":
    main()
