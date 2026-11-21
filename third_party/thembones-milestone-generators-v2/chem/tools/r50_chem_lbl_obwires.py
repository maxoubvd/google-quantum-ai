"""r50_chem_lbl_obwires: THE ONE-BODY ROWS UNDER DELTA WIRES ON ADDRESS BITS (R50-CHEM-LBL, construction 9).

A wire k:j:b:m asks in every row: bit 0 of word 5k+j XOR bit 0 of word 5k-5+j = parity(address AND m) (m = an address
mask: 1, 2, 4, .. = one address bit; any mask is as free for the loader). The 76 one-body rows sit at addresses
0..75 in any order (the verifier matches them by assignment), so: (1) the rows are ASSIGNED to addresses on the cost
lam^2 x mass x (1 - 2 |offset|) of the cheaper word of every relation that their natural parities break; (2) they
are rounded by FREE's copy of the progressive rounding (tools/r49_chem_free_obround.py: 26 stages, an L-BFGS refit
of the free angles on error_1^2 after every stage), the relation imposed AT FREEZE TIME: the second word of a pair
takes the parity the relation asks for; (3) the exact coordinate descent with steps of 2 on the wires' columns.
h1' (the two-body part of the table) is held. Loading and all numerics are FREE's / r48's functions, imported.

    python r50_chem_lbl_obwires.py TABLE.npz --wires "14:0:0:1,14:1:0:2,14:2:0:4,14:3:0:8,14:4:0:16,13:0:0:32" [--out OUT.npz]
Prints error_1 uncoded (the same procedure, no wires), with the wires in the table's row order, and with the rows
assigned. --out: the table with its one-body rows replaced (words in address order).
"""
import argparse
import os
import sys
import time

import numpy as np
from scipy.optimize import linear_sum_assignment, minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_onebody import (STAGES, angles_and_mass, check_fcidump, decode_jac, decode_pi, error1,  # noqa: E402
                                  load_fixed_shifts, load_reference, reconstruct_u_vectors, symmetric_from_lower)
from r49_chem_free_obround import near_par, polish_coded  # noqa: E402


def par_of(q):
    return np.round(q).astype(np.int64) & 1


def progressive_rel(th0, edges, lam_at, H, scale, iters, rels):
    """FREE's progressive_parity with RELATIONS: rels = [(e, f, want (rows,))]: parity(word e) ^ parity(word f) = want."""
    theta = th0.copy()
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
            xs = now * scale
            q = np.round(xs)
            for (e, f, want) in rels:
                for (a_, b_) in ((e, f), (f, e)):
                    m = frozen[:, a_] & ~frozen[:, b_]
                    if m.any():
                        q[m, b_] = near_par(xs[m, b_], par_of(qv[m, a_] * scale) ^ want[m])
            q = q / scale
            impact = lam2 * mass * (q - now) ** 2
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            newly = np.zeros(theta.shape, bool)
            newly.flat[pick] = True
            frozen |= newly
            qv[newly] = q[newly]
            for (e, f, want) in rels:                     # both words of a pair frozen in this stage: the softer one follows
                both = newly[:, e] & newly[:, f]
                viol = both & ((par_of(qv[:, e] * scale) ^ par_of(qv[:, f] * scale)) != want)
                for j in np.where(viol)[0]:
                    c, o = (f, e) if impact[j, f] <= impact[j, e] else (e, f)
                    qv[j, c] = near_par(xs[j, c], int(par_of(qv[j, o] * scale)) ^ int(want[j])) / scale
            for (e, f, want) in rels:
                done = frozen[:, e] & frozen[:, f] & (newly[:, e] | newly[:, f])
                nforced += int(done.sum())
        free = ~frozen
        if free.any():
            base = np.where(frozen, qv, theta)
            r = minimize(fun, theta[free], args=(free, base), jac=True, method="L-BFGS-B",
                         options={"maxiter": iters, "maxcor": 20, "ftol": 1e-14, "gtol": 1e-10})
            theta = base.copy()
            theta[free] = r.x
        else:
            theta = qv.copy()
    return np.round(qv * scale).astype(np.int64), nforced


def progressive_dyn(th0, edges, lam_at, H, scale, iters, wires):
    """The same progressive rounding with THE ADDRESS DECIDED AT FREEZE TIME (DynWires of r50_chem_lbl_wround.py): the
    second word of a pair takes its own nearest word whenever the rows can still be given distinct addresses."""
    from r50_chem_lbl_wround import DynWires
    nrow = th0.shape[0]
    dyn = DynWires(wires, np.arange(nrow), nrow)
    theta = th0.copy()
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros(theta.shape)
    lam2 = lam_at[:, None] ** 2
    ncol = theta.shape[1]

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
            q = np.round(xs)
            impact = lam2 * mass * (q / scale - now) ** 2
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            for flat in pick:
                j, c = divmod(int(flat), ncol)
                zc = int(q[j, c])
                if c in dyn.pair:
                    p_, w_ = dyn.pair[c]
                    if frozen[j, p_]:
                        zp = int(round(qv[j, p_] * scale))
                        if not dyn._restrict(j, w_, (zc & 1) ^ (zp & 1), True):
                            dyn.forced_n += 1
                            vcur = int(dyn.L[w_][dyn.mrow[j]])
                            zc = int(near_par(xs[j, c], (zp & 1) ^ vcur))
                            assert dyn._restrict(j, w_, vcur, True)
                        else:
                            dyn.free_ok += 1
                frozen[j, c] = True
                qv[j, c] = zc / scale
        free = ~frozen
        if free.any():
            base = np.where(frozen, qv, theta)
            r = minimize(fun, theta[free], args=(free, base), jac=True, method="L-BFGS-B",
                         options={"maxiter": iters, "maxcor": 20, "ftol": 1e-14, "gtol": 1e-10})
            theta = base.copy()
            theta[free] = r.x
        else:
            theta = qv.copy()
    return np.round(qv * scale).astype(np.int64), dyn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--wires", required=True)
    ap.add_argument("--iters", type=int, default=300)
    ap.add_argument("--sweeps", type=int, default=40)
    ap.add_argument("--out", default="")
    ap.add_argument("--skip-base", action="store_true")
    ap.add_argument("--dynamic", action="store_true", help="also: the rows' addresses decided at freeze time")
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
    wires = [tuple(int(v) for v in s.split(":")) for s in a.wires.split(",")]
    assert all(b == 0 for (_, _, b, _) in wires), "bit 0 only (the sign gauge keeps a word's parity, not its bit 1)"
    A = np.arange(N)
    want_at = np.stack([np.array([bin(int(x) & m).count("1") & 1 for x in A]) for (_, _, _, m) in wires], axis=1)   # (address, wire)
    cols = sorted({5 * k_ - 5 + j for (k_, j, _, _) in wires} | {5 * k_ + j for (k_, j, _, _) in wires})
    assert len(cols) == 2 * len(wires), "a word in two relations: use a chain-free layout"

    def run(name, th, la, rels):
        wq, nf = progressive_rel(th, edges, la, H, 64.0, a.iters, rels)
        e_p = error1(decode_pi(wq / 64.0, edges), H)[0]
        wq2, e_f = polish_coded(wq, edges, H, 7, a.sweeps, set(cols) if rels else set())
        bad = sum(int((((wq2[:, e] ^ wq2[:, f]) & 1) != wt).sum()) for (e, f, wt) in rels)
        print(f"  {name}: progressive error_1 {1e3 * e_p:.4f} mHa, after the polish {1e3 * e_f:.4f} mHa; relations {len(rels) * N}, violated {bad} "
              f"({time.time() - t0:.0f} s)", flush=True)
        return wq2, e_f

    if not a.skip_base:
        run("no wires (the same procedure)", th0, lam_at, [])
    rels0 = [(5 * k_ - 5 + j, 5 * k_ + j, want_at[:, wi]) for wi, (k_, j, _, _) in enumerate(wires)]
    run(f"{len(wires)} wires, rows in the table's order", th0, lam_at, rels0)
    # the assignment
    xs = th0 * 64.0
    par = np.round(xs).astype(np.int64) & 1
    pen = lam_at[:, None] ** 2 * mass0 * (1.0 - 2.0 * np.abs(xs - np.round(xs)))
    nat = np.stack([par[:, e] ^ par[:, f] for (e, f, _) in rels0], axis=1)                 # (row, wire)
    fp = np.stack([np.minimum(pen[:, e], pen[:, f]) for (e, f, _) in rels0], axis=1)
    cost = np.einsum("jw,jaw->ja", fp, (nat[:, None, :] != want_at[None, :, :]).astype(np.float64))
    ri, ci = linear_sum_assignment(cost)
    row_at = np.empty(N, np.int64)
    row_at[ci] = ri
    mis = nat[ri] != want_at[ci]
    print(f"  assignment: cost {cost[ri, ci].sum():.3e} against {cost[np.arange(N), np.arange(N)].sum():.3e} in the table's order; relations broken by the natural "
          f"parities {int(mis.sum())} of {mis.size} in {int(mis.any(axis=1).sum())} rows; of them in the 25 heaviest rows "
          f"{int(mis[np.isin(ri, np.argsort(-np.abs(lam_at))[:25])].sum())}", flush=True)
    wq2, e_f = run(f"{len(wires)} wires, rows ASSIGNED", th0[row_at], lam_at[row_at], rels0)
    if a.dynamic:
        wq, dyn = progressive_dyn(th0, edges, lam_at, H, 64.0, a.iters, wires)
        e_p = error1(decode_pi(wq / 64.0, edges), H)[0]
        wq3, e_d = polish_coded(wq, edges, H, 7, a.sweeps, set(cols))
        addr_of = dyn.mrow.copy()                                  # row (table order) -> address
        bad = sum(int((((wq3[:, e] ^ wq3[:, f]) & 1) != want_at[addr_of, wi]).sum()) for wi, (e, f, _) in enumerate(rels0))
        print(f"  {len(wires)} wires, ADDRESSES DECIDED AT FREEZE TIME: progressive error_1 {1e3 * e_p:.4f} mHa, after the polish {1e3 * e_d:.4f} mHa; "
              f"second words left free {dyn.free_ok}, forced {dyn.forced_n}; relations violated {bad} ({time.time() - t0:.0f} s)", flush=True)
        if e_d < e_f:
            row_at = np.empty(N, np.int64)
            row_at[addr_of] = np.arange(N)
            wq2 = wq3[row_at]
    if a.out:
        out = words.copy()
        out[:N] = np.mod(wq2, 128)
        keep = {k_: z[k_] for k_ in z.files if k_ != "angle_words"}
        np.savez(a.out, angle_words=out, ob_row_at=row_at, **keep)
        print("wrote", a.out, flush=True)


if __name__ == "__main__":
    main()
