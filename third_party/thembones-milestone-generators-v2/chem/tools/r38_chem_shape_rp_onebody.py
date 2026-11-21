#!/usr/bin/env python
"""R38-CHEM-SHAPE re-pricer: error_1 of the 76 one-body rows at a given word width, on the h1' of a real table.

Why: the verifier gates error_1 <= 20 mHa with the SAME brot as the two-body rows (one table, one brot). The
designer's 7-bit readings price only the two-body damage; the builders record error_1 at 7 bits as "a second gate
that no run has passed" (projections 15-22 mHa). This tool measures it on the entry's own h1'.

Method (a copy of tools/chem_own2_onebody.py's staged rounding + exact polish, with the bit width as a parameter
and with NO codes and NO carriers, so the word may take any of the 2^brot values on the full circle):
  h1' = H1_THC - (eta/2) S + beta_1 I - sum_rb (sum_c W_rc w_rbc) u_rb u_rb^T   (verify.evaluate_words; S, beta_1
        are the hosts' fixed shifts; u = the table's two-body rows decoded at the table's own brot);
  rows = eigenvectors of h1' -> angles on the table's tree (signed gauge) -> staged freezing (largest impact
  lam^2 x mass x offset^2 first, 26 stages) with an L-BFGS refit of the free angles on error_1^2 after every stage
  (exact gradient) -> exact coordinate descent of +-1 word moves (rows re-paired to eigenvalues by |overlap| after
  every sweep, as the verifier does).
Check first: error_1 of the table's OWN one-body rows at its own brot must equal the verifier's record.

  python r38_chem_shape_rp_onebody.py TABLE.npz --brots 8,7,6 [--table-brot 8] [--iters 300] [--sweeps 40]
         [--json out.json]
TABLE.npz: angle_words (76 + R*B, 75), tree_edges, W (R,B,C) = lowercase w, H2_I (R,C) = w_B.
Runs on a card without kickmix (only dfthcsos.scoreinfra.reference / dfthc_math are imported).
"""
import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import linear_sum_assignment, minimize

DATA = os.environ.get("QCC_DATA", os.path.expanduser("~/qcc/data"))
for p in (os.path.join(DATA, "verifiers"),):
    if p not in sys.path:
        sys.path.insert(0, p)
os.environ.setdefault("DFTHC_DATA", DATA)
from dfthcsos.scoreinfra.reference import check_fcidump, load_fixed_shifts, load_reference  # noqa: E402
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors, symmetric_from_lower  # noqa: E402

STAGES = (0.005, 0.01, 0.02, 0.035, 0.05, 0.075, 0.1, 0.14, 0.18, 0.23, 0.28, 0.34, 0.4, 0.46, 0.52, 0.58, 0.64, 0.7,
          0.76, 0.82, 0.87, 0.91, 0.94, 0.97, 0.99, 1.0)


def log(*a):
    print(time.strftime("%H:%M:%SZ", time.gmtime()), *a, flush=True)


def decode_pi(th, edges):
    """Rows of unit vectors from angles theta / pi on the tree (the verifier's Givens convention)."""
    th = np.pi * np.asarray(th, np.float64)
    v = np.zeros((th.shape[0], edges.shape[0] + 1))
    v[:, 0] = 1.0
    for t, (p, c) in enumerate(edges):
        a, b = v[:, p].copy(), v[:, c].copy()
        v[:, p] = np.cos(th[:, t]) * a - np.sin(th[:, t]) * b
        v[:, c] = np.sin(th[:, t]) * a + np.cos(th[:, t]) * b
    return v


def angles_and_mass(vectors, edges):
    """tools/chem_own2_round.angles_and_mass with K = None (signed gauge): theta / pi in [-1/2, 1/2] and the
    squared amplitude entering each rotation."""
    a = np.array(vectors, dtype=np.float64)
    a /= np.linalg.norm(a, axis=1, keepdims=True)
    theta = np.zeros((len(a), len(edges)))
    mass = np.zeros((len(a), len(edges)))
    for t in range(len(edges) - 1, -1, -1):
        p, q = edges[t]
        s = np.where(a[:, p] > 0, 1.0, np.where(a[:, p] < 0, -1.0, np.where(a[:, q] > 0, -1.0, 1.0)))
        theta[:, t] = np.arctan2(s * a[:, q], s * a[:, p])
        h = np.hypot(a[:, p], a[:, q])
        mass[:, t] = h * h
        a[:, p] = s * h
        a[:, q] = 0.0
    return theta / np.pi, mass


def decode_jac(th, edges):
    """V (rows, N) and J (rows, ncol, N) = dV / d(theta / pi) (tools/chem_own2_onebody.decode_jac, K = 0)."""
    rows, ncol = th.shape
    ang = np.pi * th
    c, s = np.cos(ang), np.sin(ang)
    v = np.zeros((rows * (ncol + 1), ncol + 1))
    v[:, 0] = 1.0
    first = np.arange(rows) * (ncol + 1)
    for t, (p, q) in enumerate(edges):
        a = v[:, p].copy()
        v[:, p] = np.repeat(c[:, t], ncol + 1) * a
        v[:, q] = np.repeat(s[:, t], ncol + 1) * a
        i = first + t
        v[i] = 0.0
        v[i, p] = -s[:, t] * a[i]
        v[i, q] = c[:, t] * a[i]
    v = v.reshape(rows, ncol + 1, ncol + 1)
    return v[:, ncol], np.pi * v[:, :ncol]


def error1(V, H):
    """The verifier's error_1: rows paired one-to-one with exact eigenvectors of H by |overlap|."""
    lam, vec = np.linalg.eigh(H)
    rows, cols = linear_sum_assignment(-np.abs(V @ vec))
    perm = np.empty(len(lam), np.int64)
    perm[rows] = cols
    D = (V.T * lam[perm]) @ V - H
    return float(np.linalg.norm(D)), perm, D, lam


def progressive(th0, edges, lam_at, H, scale, iters):
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
            q = np.round(now * scale) / scale
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


def polish(wq, edges, H, brot, sweeps, inner=25):
    """Exact coordinate descent on error_1 over single-word moves +-1 (full circle, words mod 2^brot)."""
    scale, mod = float(1 << (brot - 1)), 1 << brot
    N = edges.shape[0] + 1
    ncol = wq.shape[1]
    wq = np.mod(np.array(wq, np.int64), mod)
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
                cand[ar, ar] += 1
                cand[ncol + ar, ar] -= 1
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
        if sweep % 5 == 0 or not moved:
            log(f"    brot {brot} sweep {sweep}: {moved} moves, error_1 {1e3 * e1:.3f} mHa")
        if not moved or e_prev - e1 < 1e-7 * e1:
            break
    return wq, e1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--brots", default="8,7")
    ap.add_argument("--table-brot", type=int, default=8)
    ap.add_argument("--iters", type=int, default=300)
    ap.add_argument("--sweeps", type=int, default=40)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    t0 = time.time()
    z = np.load(a.table)
    words = np.asarray(z["angle_words"], np.int64)
    edges = np.asarray(z["tree_edges"], np.int64)
    w = np.asarray(z["W"], np.float64)
    W = np.asarray(z["H2_I"], np.float64) - w.sum(axis=1)
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
    allv = reconstruct_u_vectors(words / float(1 << (a.table_brot - 1)), edges)
    allv = allv / np.linalg.norm(allv, axis=1, keepdims=True)
    u = allv[N:].reshape(R, B, N)
    k = np.einsum("rc,rbc->rb", W, w)
    H = h1_thc - 0.5 * eta * S + beta * np.eye(N) - np.einsum("rb,rbp,rbq->pq", k, u, u, optimize=True)
    e_own, _, _, lam = error1(allv[:N], H)
    out = dict(table=a.table, shape=[R, B, C], table_brot=a.table_brot, own_rows_error_1_mHa=1e3 * e_own,
               sum_abs_lam=float(np.abs(lam).sum()), fro_h1p=float(np.linalg.norm(H)), runs={})
    log(f"table {a.table} {R}x{B}x{C}: own one-body rows error_1 {1e3 * e_own:.4f} mHa at brot {a.table_brot} "
        f"(check against the verifier's record); sum|lam| {np.abs(lam).sum():.4f} Ha, ||h1'||_F {np.linalg.norm(H):.4f}")
    lam_e, vec = np.linalg.eigh(H)
    th0, _ = angles_and_mass(vec.T, edges)
    for brot in [int(x) for x in a.brots.split(",")]:
        scale = float(1 << (brot - 1))
        tb = time.time()
        wn = np.mod(np.round(th0 * scale).astype(np.int64), 1 << brot)
        e_near = error1(decode_pi(wn / scale, edges), H)[0]
        wq = progressive(th0, edges, lam_e, H, scale, a.iters)
        e_prog = error1(decode_pi(np.mod(wq, 1 << brot) / scale, edges), H)[0]
        log(f"  brot {brot}: nearest {1e3 * e_near:.3f} mHa; staged {1e3 * e_prog:.3f} mHa; polishing")
        wq, e_pol = polish(wq, edges, H, brot, a.sweeps)
        # the verifier's own decoding of the final words (same convention) as a check
        e_chk = error1(reconstruct_u_vectors(np.mod(wq, 1 << brot) / scale, edges), H)[0]
        out["runs"][str(brot)] = dict(nearest=1e3 * e_near, staged=1e3 * e_prog, polished=1e3 * e_pol,
                                      verifier_decode=1e3 * e_chk, seconds=round(time.time() - tb))
        log(f"  brot {brot}: nearest {1e3 * e_near:.3f} | staged {1e3 * e_prog:.3f} | polished {1e3 * e_pol:.3f} mHa "
            f"(verifier decode {1e3 * e_chk:.3f}; {time.time() - tb:.0f} s)")
    out["seconds"] = round(time.time() - t0)
    print(json.dumps(out), flush=True)
    if a.json:
        json.dump(out, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
