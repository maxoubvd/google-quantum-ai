#!/usr/bin/env python
"""chem_own2_onebody (chem_own2, 2026-11-14): a copy of tools/chem_own_onebody.py whose codes come from
tools/chem_own2_round.py (its residues and check_codes know code 3 = even 16-row parity and code 4 = even parity
on both 8-row halves, L2S's low-degree deltas); nothing else changed.  Original header follows.

chem_own_onebody: the 76 one-body rows of a dfthcsos table, from our own two-body table (Level B piece, CPU).

    chem_own_onebody.py <table.npz> <out.npz> [--codes "0:1:3-14"] [--xcodes ""] [--carriers 3,5,6,9,10]
                        [--sweeps 40] [--order deal|eig]

<table.npz>: `angle_words` (76 + R*B rows; rows 0..75 are ignored and rebuilt), `tree_edges`, `W`, `H2_I`.
<out.npz>: the same table with rows 0..75 replaced.

Method. h1' = h1_base - sum_rb (sum_c W_rc w_rbc) u_rb u_rb^T depends only on the two-body rows and the
weights. Its eigenvectors are the ideal one-body rows; the verifier pairs every submitted row with an exact
eigenvector by |overlap| and measures error_1 = || V^T diag(lam) V - h1' ||_F. So the order of the rows is free
and every row may be written in any exact form. Here:
  1. order: the 64 eigenvectors with the largest |eigenvalue| are dealt, strongest first, over the four full
     address blocks 0..3; the 12 weakest go to addresses 64..75, the block they share with the two-body rows
     76..79 (whose words already fix most of that block's codewords);
  2. form: the carrier pattern of the table (bit7 xor bit6 of word position j = j-th low-address monomial),
     put on exactly as a gauge of the continuous angles (chem_own_round.angles_and_mass);
  3. codes: the same code set as the two-body rows (chem_own_round.residues: the cheapest codeword per block,
     column and bit; on the shared block only codewords that agree with the frozen two-body rows 76..79);
  4. staged rounding with a refit after every stage (`progressive`): the two-body method on error_1. The entries
     with the largest impact are frozen first, the free angles are refitted on error_1^2 with the exact gradient
     (`decode_jac`), every row is re-gauged exactly at the end of a stage;
  5. exact coordinate descent on error_1: single-word moves that keep every code (step = 2^(coded low bits of
     the column, per-column and XOR-of-five alike)) and keep the word in the signed range, priced exactly from
     D = V^T diag(lam) V - h1'; the best move of a row is applied, up to 25 per row and sweep.
The eigenvalues are re-matched to the rows after every sweep, as the verifier does.

Measured (2026-11-11): first version (round at once, then polish; public two-body rows, bit-0 code, no carriers,
shared block not enforced) 64.1 -> 16.63 mHa; with carriers, code and shared block on our rehearsal table 58.1
-> 18.54; **with the staged rounding 9.85, after the polish 9.61 mHa** (the public pipeline's own rows: 9.3-9.7
on its (8,143,45) table, 11.26 on the public table), 163 s on one core.
All code here is ours; the staged idea (round, then exact +-1 polish) is the public one (mrmeows "one-body
refresh", marculera `onebody.polish`), re-implemented from its description. Codes, carriers and the gauge come
from tools/chem_own_round.py (credits there).
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_own2_round as RD  # noqa: E402  (chem_own2: codes 3 and 4)

HALF = RD.HALF


def decode(words, edges, brot=8):
    th = np.pi * np.asarray(words, np.float64) / (1 << (brot - 1))
    v = np.zeros((th.shape[0], edges.shape[0] + 1))
    v[:, 0] = 1.0
    for t, (p, c) in enumerate(edges):
        a, b = v[:, p].copy(), v[:, c].copy()
        v[:, p] = np.cos(th[:, t]) * a - np.sin(th[:, t]) * b
        v[:, c] = np.sin(th[:, t]) * a + np.cos(th[:, t]) * b
    return v


def signed_words(stored, K):
    """Stored 8-bit words -> signed words in -64..63 after removing the half turns K (raises if a word does
    not carry its pattern)."""
    s = np.mod(np.asarray(stored, np.int64) - 128 * np.asarray(K, np.int64) + HALF, 256) - HALF
    if s.max() > HALF - 1:
        raise ValueError("a word does not carry the carrier pattern (bit7 xor bit6 differs)")
    return s


def error1(V, lam, H):
    rows, cols = linear_sum_assignment(-np.abs(V @ np.linalg.eigh(H)[1]))
    perm = np.empty(len(lam), np.int64)
    perm[rows] = cols
    D = (V.T * lam[perm]) @ V - H
    return float(np.linalg.norm(D)), perm, D


def polish(wq, K, edges, lam, H, step, sweeps=40, log=print, inner=25):
    """Exact coordinate descent on error_1 over single-word moves of +-step[column] inside the signed range.
    wq: signed words (76, 75); K: half turns. Returns (wq, error_1)."""
    N = edges.shape[0] + 1
    ncol = wq.shape[1]
    wq = np.array(wq, np.int64)
    K128 = 128 * np.asarray(K, np.int64)
    V = decode(np.mod(wq + K128, 256), edges)
    e1, perm, D = error1(V, lam, H)
    log(f"  polish start: error_1 {1e3 * e1:.3f} mHa")
    ar = np.arange(ncol)
    for sweep in range(sweeps):
        moved = 0
        for j in range(N):
            lj = lam[perm[j]]
            base = D - lj * np.outer(V[j], V[j])                  # D without row j
            for _ in range(inner):                                # best single-word move of the row, repeated
                v0 = V[j]
                cand = np.repeat(wq[j][None, :], 2 * ncol, axis=0)
                cand[ar, ar] += step
                cand[ncol + ar, ar] -= step
                ok = (cand.max(axis=1) <= HALF - 1) & (cand.min(axis=1) >= -HALF)
                Vc = decode(np.mod(cand + K128[j], 256), edges)    # (2 ncol, N)
                # || base + lj v v^T ||^2 = ||base||^2 + 2 lj v^T base v + lj^2 |v|^4
                quad = np.einsum("ip,pq,iq->i", Vc, base, Vc)
                n4 = np.sum(Vc * Vc, axis=1) ** 2
                val = np.where(ok, 2 * lj * quad + lj * lj * n4, np.inf)
                cur = 2 * lj * float(v0 @ base @ v0) + lj * lj * float(v0 @ v0) ** 2
                i = int(np.argmin(val))
                if val[i] >= cur - 1e-13:                         # (a smaller threshold accepts rounding noise and cycles)
                    break
                wq[j] = cand[i]
                V[j] = Vc[i]
                moved += 1
            D = base + lj * np.outer(V[j], V[j])
        e_prev = e1
        e1, perm, D = error1(V, lam, H)
        log(f"  sweep {sweep}: {moved} moves, error_1 {1e3 * e1:.3f} mHa")
        if not moved or e_prev - e1 < 1e-7 * e1:
            break
    return wq, e1


STAGES = (0.005, 0.01, 0.02, 0.035, 0.05, 0.075, 0.1, 0.14, 0.18, 0.23, 0.28, 0.34, 0.4, 0.46, 0.52, 0.58, 0.64, 0.7,
          0.76, 0.82, 0.87, 0.91, 0.94, 0.97, 0.99, 1.0)


def decode_jac(th, K, edges):
    """V (rows, N) and J (rows, ncol, N) = dV / d(theta / pi) for the angles theta + K (units of pi)."""
    rows, ncol = th.shape
    ang = np.pi * (th + K)
    c, s = np.cos(ang), np.sin(ang)
    # every row is decoded ncol + 1 times; in copy k the rotation k is replaced by its derivative, which is zero
    # outside the plane of that rotation: the state of the copy is cleared except for the two rotated entries
    v = np.zeros((rows * (ncol + 1), ncol + 1))
    v[:, 0] = 1.0
    first = np.arange(rows) * (ncol + 1)
    for t, (p, q) in enumerate(edges):
        a = v[:, p].copy()
        v[:, p] = np.repeat(c[:, t], ncol + 1) * a   # the child is new: its old value is zero
        v[:, q] = np.repeat(s[:, t], ncol + 1) * a
        i = first + t
        v[i] = 0.0
        v[i, p] = -s[:, t] * a[i]
        v[i, q] = c[:, t] * a[i]
    v = v.reshape(rows, ncol + 1, ncol + 1)
    return v[:, ncol], np.pi * v[:, :ncol]


def progressive(th, K, edges, lam_at, H, rest_T, rest_value, deg, xdeg, blk, pos, iters=300, log=print):
    """Staged rounding of the one-body rows with a refit of the free angles after every stage (the two-body
    method, on error_1): the entries with the largest impact lam^2 x mass x offset^2 are frozen first on the
    cheapest codewords that agree with what is frozen (rows of the shared block included), the free angles are
    refitted on error_1^2 = || V^T diag(lam) V - h1' ||^2 with the exact gradient (L-BFGS-B, box only where a
    crossing is not allowed), and every row is re-gauged exactly at the end of the stage. Returns signed words."""
    from scipy.optimize import minimize
    N, ncol = th.shape[0], th.shape[1]
    lo, hi = -0.5 - 0.5 / 128, 0.5 - 0.5 / 128
    NB = RD.neighbour_matrix(edges)
    pcbits = (deg > 0).sum(axis=1)
    xcol = np.repeat(xdeg[:, 1:].any(axis=1), 5)
    theta = th.copy()
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros(theta.shape)
    nrest = len(rest_T)
    lam2 = lam_at[:, None] ** 2

    def fun(x, mask, base):
        t = base.copy()
        t[mask] = x
        V, J = decode_jac(t, K, edges)
        D = (V.T * lam_at) @ V - H
        gV = 4.0 * lam_at[:, None] * (V @ D)
        return float(np.sum(D * D)), np.einsum("jn,jkn->jk", gV, J)[mask]

    for frac in STAGES:
        need = int(round(frac * theta.size)) - int(frozen.sum())
        if need > 0:
            now = np.where(frozen, qv, theta)
            mass = RD.angles_and_mass(RD.decode_np(now + K, edges), edges)[1]
            c = lam2 * mass + 1e-30
            T = np.vstack([now * 128.0, rest_T])
            fr = np.vstack([frozen, np.ones((nrest, ncol), bool)])
            val = np.vstack([np.round(qv * 128).astype(np.int64), rest_value])
            cc = np.vstack([c, np.full((nrest, ncol), float(c.max()))])
            res = None
            if xdeg.any():
                res = RD.xor_slack(cc, blk, pos, xdeg)
                res[N:] = False
            r_, m_ = RD.residues(T, cc, blk, pos, deg, xdeg, fr, val, res)
            q = RD.into_range(RD.residue_round(T, r_, m_), m_)[:N] / 128.0
            impact = c * (q - now) ** 2
            if res is not None:
                impact = np.where(res[:N], -1.0 / (1.0 + impact), impact)
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            frozen.flat[pick] = True
            qv.flat[pick] = q.flat[pick]
        wv = np.round(qv * 128).astype(np.int64)
        modh = np.where(pcbits >= 2, 1 << np.maximum(pcbits - 1, 0), 1)[None, :]
        protected = frozen & (((pcbits[None, :] >= 2) & (np.mod(wv, modh) != 0)) | (xcol[None, :] & (wv != 0)))
        negat, flip_ok = RD.negatable(frozen, wv, protected, NB)
        free = ~frozen
        if free.any():
            base = np.where(frozen, qv, theta)
            wall = (free & ~flip_ok)[free]
            x0 = np.where(wall, np.clip(theta[free], lo, hi - 1e-9), theta[free])
            bounds = [(lo, hi - 1e-9) if w_ else (None, None) for w_ in wall]
            f0 = fun(x0, free, base)[0]
            r = minimize(fun, x0, args=(free, base), jac=True, method="L-BFGS-B", bounds=bounds,
                         options={"maxiter": iters, "maxcor": 20, "ftol": 1e-14, "gtol": 1e-10})
            theta = base.copy()
            theta[free] = r.x
            if os.environ.get("ONEBODY_TRACE"):
                log(f"    stage {frac:.3f}: error_1 {1e3 * np.sqrt(f0):.3f} -> {1e3 * np.sqrt(r.fun):.3f} mHa, {r.nit} its, "
                    f"{int(wall.sum())} walled, {str(r.message)[:40]}")
        else:
            theta = qv.copy()
        theta, nf, nn, ncl = RD.regauge(theta, frozen, flip_ok, negat, edges, lo, hi)
        qv = np.where(frozen, theta, qv)
    V = RD.decode_np(qv + K, edges)
    D = (V.T * lam_at) @ V - H
    log(f"  progressive rounding ({len(STAGES)} stages x {iters}): error_1 {1e3 * np.linalg.norm(D):.3f} mHa "
        f"(row order as built, before the exact polish)")
    return np.round(qv * 128).astype(np.int64)


def build_rows(H, edges, two_body_head, deg, xdeg, masks, sweeps=40, order="deal", log=print, staged=True):
    """The 76 one-body rows for h1' = H as stored 8-bit words, and their error_1.
    two_body_head: the stored words of the two-body rows that share the last one-body block (addresses 76..79)."""
    N = edges.shape[0] + 1
    nhead = len(two_body_head)
    assert nhead == (-N) % 16, "the two-body rows that complete the last one-body block are needed"
    lam, vec = np.linalg.eigh(H)
    idx = np.argsort(-lam ** 2, kind="stable")
    at = np.arange(N)
    if order == "deal":
        nb = N // 16
        for i, e in enumerate(idx[:16 * nb]):
            at[16 * (i % nb) + i // nb] = e
        at[16 * nb:] = idx[16 * nb:]
    elif order == "d1":     # R8, Nov 15: R41-CHEM-CODE Part 3's dependent block D1 = the 16 WEAKEST rows in block 0
        nb = N // 16        # (C2's fill then rewrites their bits 1-2 for little error_1); the next 12 weakest go to
        at[:16] = idx[-16:]  # the shared block (addresses 64..75); the 48 strongest are dealt over blocks 1..3
        rest = idx[:-16]
        for i, e in enumerate(rest[:16 * (nb - 1)]):
            at[16 * (1 + i % (nb - 1)) + i // (nb - 1)] = e
        at[16 * nb:] = rest[16 * (nb - 1):]
    V0, lam_at = vec.T[at], lam[at]
    addr = np.arange(N + nhead)
    blk, pos = addr // 16, addr % 16
    Kall = RD.carrier_pattern(addr, masks)
    K = Kall[:N]
    lo, hi = -0.5 - 0.5 / 128, 0.5 - 0.5 / 128
    NB = RD.neighbour_matrix(edges)
    th, mass = RD.angles_and_mass(V0, edges, K)
    fr0 = np.zeros(th.shape, bool)
    ng, fo = RD.negatable(fr0, np.zeros(th.shape, np.int64), fr0, NB)
    th = RD.regauge(th, fr0, fo, ng, edges, lo, hi)[0]
    c = lam_at[:, None] ** 2 * mass + 1e-30
    s2 = signed_words(two_body_head, Kall[N:])
    T = np.vstack([th * 128.0, s2.astype(np.float64)])
    frozen = np.zeros(T.shape, bool)
    frozen[N:] = True
    value = np.zeros(T.shape, np.int64)
    value[N:] = s2
    cc = np.vstack([c, np.full(s2.shape, float(c.max()))])
    res = None
    if xdeg.any():
        res = RD.xor_slack(cc, blk, pos, xdeg)
        res[N:] = False
    r_, m_ = RD.residues(T, cc, blk, pos, deg, xdeg, frozen, value, res)
    q = RD.into_range(RD.residue_round(T, r_, m_), m_).astype(np.int64)
    if not np.array_equal(q[N:], s2):
        raise ValueError("the two-body rows of the shared block do not lie on a codeword")
    RD.check_codes(q, addr, deg, xdeg)
    wq = q[:N]
    e0 = error1(decode(np.mod(wq + 128 * K.astype(np.int64), 256), edges), lam, H)[0]
    log(f"  one-body rows: sum|lam| {np.abs(lam).sum():.4f} Ha, order {order}; coded rounding: error_1 "
        f"{1e3 * e0:.3f} mHa ({int((m_ > 1).sum())} coded columns, largest step {int(m_.max())})")
    if staged and sweeps:
        wq = progressive(th, K, edges, lam_at, H, s2.astype(np.float64), s2, deg, xdeg, blk, pos, log=log)
        RD.check_codes(np.vstack([wq, s2]), addr, deg, xdeg)
    if sweeps:
        wq, e0 = polish(wq, K, edges, lam, H, m_.astype(np.int64), sweeps, log)
    RD.check_codes(np.vstack([wq, s2]), addr, deg, xdeg)
    assert wq.min() >= -HALF and wq.max() <= HALF - 1
    return np.mod(wq + 128 * K.astype(np.int64), 256), e0


def h1_prime(prob, w, W):
    k = np.einsum("rc,rbc->rb", W, w)
    return prob.h1_base - np.einsum("rb,rbp,rbq->pq", k, prob.u, prob.u, optimize=True)


def set_rows(prob, rows):
    """Put new one-body rows into a chem_own_refit.Problem (the verifier's own decoding)."""
    from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
    N = prob.N
    prob.words[:N] = rows
    v = reconstruct_u_vectors(prob.words[:N] / (1 << (prob.brot - 1)), prob.edges)
    prob.Vrows = v / np.linalg.norm(v, axis=1, keepdims=True)
    prob.perm = None


def main():
    import chem_own_refit as cor
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table")
    ap.add_argument("out")
    ap.add_argument("--sweeps", type=int, default=40)
    ap.add_argument("--codes", default="0:1:3-14")
    ap.add_argument("--xcodes", default="")
    ap.add_argument("--carriers", default="3,5,6,9,10")
    ap.add_argument("--order", default="deal", choices=["deal", "eig", "d1"])
    a = ap.parse_args()
    t0 = time.time()
    data = cor.load_table(a.table)
    prob = cor.Problem(data["angle_words"], data["tree_edges"], data["W"].shape, 8)
    w = np.asarray(data["W"], np.float64)
    W = np.asarray(data["H2_I"], np.float64) - w.sum(axis=1)
    N = prob.N
    masks = [] if a.carriers in ("", "none") else [int(x) for x in a.carriers.split(",")]
    deg, xdeg = RD.parse_codes(a.codes), RD.parse_xcodes(a.xcodes)
    nhead = (-N) % 16
    rows, e1 = build_rows(h1_prime(prob, w, W), prob.edges, prob.words[N:N + nhead], deg, xdeg, masks, a.sweeps, a.order)
    set_rows(prob, rows)
    out = dict(data)
    out["angle_words"] = prob.words.astype(data["angle_words"].dtype)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez(a.out, **out)
    print("official with the new one-body rows:", prob.official(w, W), f"({time.time() - t0:.0f} s)", flush=True)
    print("saved", a.out, flush=True)


if __name__ == "__main__":
    main()
