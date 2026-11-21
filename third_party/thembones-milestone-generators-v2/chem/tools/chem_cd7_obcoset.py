"""chem_cd7_obcoset (CD7): the ONE-BODY rows rounded to 7 bits under the bit-0 code in its COSET form, decided by an
EXACT WEIGHTED decoder: a copy of the staged rounding of tools/r49_chem_free_obround.py (26 stages, impact order
lam^2 x mass x offset^2, an L-BFGS refit of the free angles on error_1^2 after every stage), where at every stage the
parities of the words frozen in that stage come from tools/chem_cd7_round.py's shared_coset_bits (code d == 11): per
(block, slot) the cheapest common RM(1,4) syndrome over the slot's coded columns given the words already frozen
(weights lam^2 x mass, penalties from the present angles). The column form ("every coded column a codeword") is the
special case syndrome 0, so this decoder has strictly more freedom. Polish: +-2 on coded columns (parity kept).

The loaders need only the coset form: delta_k = column 5(k-1)+j XOR column 5k+j affine on every block, true when both
columns share the block's syndrome. Block 4 holds the two-body rows 76-79 at low addresses 12-15: in the test mode
they are left out (as obround does); in the chain (ROWS76 given) their words enter the decoder with their own
parities as prices (cost of a parity change = w76 x a one-step move), so the block's syndrome can be met by them.

    python chem_cd7_obcoset.py TABLE.npz --cols 15-74 [--iters 300] [--sweeps 40]      (fixed-h1' test, as obround)
"""
import argparse
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_cd7_round as CR  # noqa: E402
import r49_chem_free_obround as OR  # noqa: E402
from r48_chem_q7c_onebody import (STAGES, angles_and_mass, check_fcidump, decode_jac, decode_pi, error1,  # noqa: E402
                                  load_fixed_shifts, load_reference, reconstruct_u_vectors, symmetric_from_lower)

CW = OR.CW


def near_par(x, p):
    return 2.0 * np.round((x - p) / 2.0) + p


def coset_targets(xs, wgt, frozen, val, coded, nrow, extra=None):
    """Bits (nrow, ncol) chosen by the exact shared-coset decoder on the coded columns (others 0). xs: present
    angles in steps (nrow, ncol); wgt: weights; frozen / val: frozen words and their integer values. extra: optional
    (pos 12..15 of block 4) = (bits (4, ncol), price (4, ncol)) of the two-body rows 76-79."""
    cols = np.array(sorted(coded), np.int64)
    nb = 5
    t0 = np.zeros((nb, 16, len(cols)))
    t1 = np.zeros((nb, 16, len(cols)))
    x = xs[:, cols]
    c = wgt[:, cols]
    p0 = c * (near_par(x, 0) - x) ** 2
    p1 = c * (near_par(x, 1) - x) ** 2
    big = max(1e4 * max(float(p0.max()), float(p1.max()), 1e-30) + 1.0, 1e7)    # above any price of rows 76-79
    fb = (val[:, cols] & 1)
    fz = frozen[:, cols]
    p0 = np.where(fz, np.where(fb == 0, 0.0, big), p0)
    p1 = np.where(fz, np.where(fb == 1, 0.0, big), p1)
    rows = np.arange(nrow)
    t0[rows // 16, rows % 16] = p0
    t1[rows // 16, rows % 16] = p1
    if extra is not None:
        eb, ep = extra
        for i in range(4):
            t0[4, 12 + i] = np.where(eb[i, cols] == 0, 0.0, ep[i, cols])
            t1[4, 12 + i] = np.where(eb[i, cols] == 1, 0.0, ep[i, cols])
    bits = CR.shared_coset_bits(t0, t1, cols, 11)                       # (5, 16, ncols)
    out = np.zeros(xs.shape, np.int64)
    out[:, cols] = bits[rows // 16, rows % 16]
    b76 = bits[4, 12:16] if extra is not None else None
    return out, b76, cols


def progressive_cosetdp(th0, edges, lam_at, H, scale, iters, coded, extra=None, wires=None):
    """CD7 (third holder, 20:55Z): wires = [(c1, c0, m), ..] bit-0 delta wires at the rows' FIXED addresses (row index =
    address): the second word of a relation to freeze takes bit 0 = partner's bit 0 XOR parity(address AND m)
    (R50's DynChain without the address freedom, which a coded table does not have). Disjoint from the coded columns."""
    theta = th0.copy()
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros(theta.shape)
    val = np.zeros(theta.shape, np.int64)
    lam2 = lam_at[:, None] ** 2
    nrow = theta.shape[0]
    cmask = np.zeros(theta.shape[1], bool)
    cmask[sorted(coded)] = True
    partner = {}                                   # UNCODED column -> list of (partner column, mask)
    for c1, c0, m in (wires or []):
        assert not (cmask[c1] and cmask[c0]), "a wire with both columns coded"
        if not cmask[c1]:                          # 21:25Z: a coded partner's bit comes from the decoder (frozen
            partner.setdefault(c1, []).append((c0, m))     # value, else its present target); only the uncoded side
        if not cmask[c0]:                          # of a relation is ever forced by the wire
            partner.setdefault(c0, []).append((c1, m))
    nforced = nmoved = 0
    b76 = None

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
            xs = now * scale
            tb, b76, _ = coset_targets(xs, lam2 * mass + 1e-30, frozen, val, coded, nrow, extra)
            for fl in pick:
                j, t = divmod(int(fl), theta.shape[1])
                z = np.round(xs[j, t])
                if cmask[t]:
                    z2 = near_par(xs[j, t], int(tb[j, t]))
                    nforced += 1
                    nmoved += int(z2 != z)
                    z = z2
                elif t in partner:
                    for (t2, m) in partner[t]:
                        if frozen[j, t2] or cmask[t2]:
                            b2 = (int(val[j, t2]) & 1) if frozen[j, t2] else int(tb[j, t2])
                            want = b2 ^ (bin(j & m).count("1") & 1)
                            z2 = near_par(xs[j, t], want)
                            nforced += 1
                            nmoved += int(z2 != z)
                            z = z2
                            break
                frozen[j, t] = True
                qv[j, t] = z / scale
                val[j, t] = int(z)
        free = ~frozen
        if free.any():
            base = np.where(frozen, qv, theta)
            r = minimize(fun, theta[free], args=(free, base), jac=True, method="L-BFGS-B",
                         options={"maxiter": iters, "maxcor": 20, "ftol": 1e-14, "gtol": 1e-10})
            theta = base.copy()
            theta[free] = r.x
        else:
            theta = qv.copy()
    if extra is not None:                      # the final syndrome choice with everything frozen
        _, b76, _ = coset_targets(qv * scale, np.ones(qv.shape), frozen, val, coded, nrow, extra)
    wq = np.round(qv * scale).astype(np.int64)
    for c1, c0, m in (wires or []):            # 21:25Z: a relation whose coded side changed after its uncoded side
        cu = c0 if cmask[c1] else c1           # froze: the uncoded word moves one step (the polish follows)
        for j in range(nrow):
            if ((int(wq[j, c1]) ^ int(wq[j, c0])) & 1) != (bin(j & m).count("1") & 1):
                z = int(wq[j, cu])
                wq[j, cu] = z - 1 if abs(z - 1) < abs(z + 1) else z + 1
                nmoved += 1
    return wq, nforced, nmoved, b76


def progressive_fp(th0, edges, lam_at, H, scale, iters, chains, wires, w76):
    """CD7 (21:40Z): the 76 one-body rows rounded with their PLACES decided at freeze time by R50-CHEM-LBL's
    FreezePlacer (tools/r50_chem_lbl_hook.py, unchanged) for the bit-0 code chains `chains` (one coset per 16-address
    block) and the bit-0 wires `wires` [(k, j, 0, m)], on addresses 0..79: places 0..75 are the one-body rows'
    (free; the verifier pairs one-body rows with eigenvectors by overlap, so their order costs nothing), places 76-79
    hold the two-body rows 76-79 (w76: their signed words, all frozen, PINNED: each its own group). Staged rounding as
    progressive_cosetdp (impact order, L-BFGS refit of the free angles on error_1^2). Returns the words in place order
    (76, ncol), the place of every input row, and the counts."""
    from r50_chem_lbl_hook import FreezePlacer
    import r51_chem_bit1_round as RD
    RD.HALF = 1 << 20
    n, nc = th0.shape
    fp = FreezePlacer(1, n + 4, chains, wires, n_one=0)
    fp.R = 5
    fp.group = np.concatenate([np.zeros(n, np.int64), np.arange(1, 5)])         # rows 76-79 pinned
    fpcol = np.zeros(nc, bool)
    fpcol[list(fp.pl.cc)] = True
    theta = th0.copy()
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros(theta.shape)
    val = np.zeros(theta.shape, np.int64)
    lam2 = lam_at[:, None] ** 2
    w76s = np.asarray(w76, np.int64)
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
            t_e = np.vstack([now * scale, w76s.astype(np.float64)])
            c_e = np.vstack([lam2 * mass + 1e-30, np.full((4, nc), 1e6)])
            f_e = np.vstack([frozen, np.ones((4, nc), bool)])
            v_e = np.vstack([val, w76s])
            r_, m_ = fp.residues(t_e, c_e, f_e, v_e)
            for fl in pick:
                j, t = divmod(int(fl), nc)
                xs = now[j, t] * scale
                z = np.round(xs)
                if fpcol[t]:
                    z2 = near_par(xs, int(r_[j, t]) & 1)
                    nforced += 1
                    nmoved += int(z2 != z)
                    z = z2
                frozen[j, t] = True
                qv[j, t] = z / scale
                val[j, t] = int(z)
        free = ~frozen
        if free.any():
            base = np.where(frozen, qv, theta)
            r = minimize(fun, theta[free], args=(free, base), jac=True, method="L-BFGS-B",
                         options={"maxiter": iters, "maxcor": 20, "ftol": 1e-14, "gtol": 1e-10})
            theta = base.copy()
            theta[free] = r.x
        else:
            theta = qv.copy()
    place = fp.place[:n]
    assert sorted(place.tolist()) == list(range(n)) and np.all(fp.place[n:] == np.arange(n, n + 4)), "pinned rows moved"
    wq = np.round(qv * scale).astype(np.int64)
    out = np.empty_like(wq)
    out[place] = wq
    return out, place, nforced, nmoved


def check_coset_full(wq, coded, b76=None):
    """Violations of the coset form on blocks 0-4 (with rows 76-79's bits b76 on block 4 when given)."""
    bad = 0
    cols = sorted(coded)
    for blk in range(5):
        rows = np.arange(16 * blk, min(16 * blk + 16, wq.shape[0]))
        for t in cols:
            if t + 5 in coded:
                bits = (wq[rows, t] ^ wq[rows, t + 5]) & 1
                if blk == 4 and b76 is not None:
                    ci, ci5 = cols.index(t), cols.index(t + 5)
                    bits = np.concatenate([bits, (b76[:, ci] ^ b76[:, ci5]) & 1])
                n = len(bits)
                if not (CW[:, :n] == bits[None, :]).all(axis=1).any():
                    bad += 1
    return bad


def heavy5(lam_at, N):
    rank = np.argsort(-np.abs(lam_at), kind="stable")
    ind = [16 * b_ + l_ for l_ in (0, 1, 2, 4, 8) for b_ in range(5)]
    rest = [a_ for a_ in range(N) if a_ not in ind]
    rest.sort(key=lambda a_: (a_ % 16, a_ // 16))
    row_at = np.empty(N, np.int64)
    row_at[np.array(ind + rest)] = rank
    return row_at


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--cols", default="15-74")
    ap.add_argument("--iters", type=int, default=300)
    ap.add_argument("--sweeps", type=int, default=40)
    ap.add_argument("--also-column", action="store_true", help="also run obround's weighted column form")
    ap.add_argument("--eig-order", action="store_true", help="start from eigh's order instead of the table's rows")
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
    e_own, perm, _, _ = error1(allv[:N], H)
    print(f"table {a.table}: its own one-body rows error_1 {1e3 * e_own:.4f} mHa ({time.time() - t0:.0f} s)", flush=True)
    lam_e, vec = np.linalg.eigh(H)
    if a.eig_order:
        th0, _ = angles_and_mass(vec.T, edges)
        lam_at = lam_e
    else:                                   # obround's order: eigenvectors matched to the table's own rows
        th0, _ = angles_and_mass(vec[:, perm].T, edges)
        lam_at = lam_e[perm]
    row_at = heavy5(lam_at, N)
    th0, lam_at = th0[row_at], lam_at[row_at]
    lo, hi = (int(v) for v in a.cols.split("-"))
    coded = set(range(lo, hi + 1))
    if a.also_column:
        wq, nf, nm = OR.progressive_coded(th0, edges, lam_at, H, 64.0, a.iters, coded, True)
        e_p = error1(decode_pi(wq / 64.0, edges), H)[0]
        wq2, e_f = OR.polish_coded(wq, edges, H, 7, a.sweeps, coded)
        print(f"  column form, weighted (obround): progressive {1e3 * e_p:.4f} mHa, polished {1e3 * e_f:.4f} mHa; "
              f"code violations {OR.check_code(wq2, coded)} ({time.time() - t0:.0f} s)", flush=True)
    wq, nf, nm, _ = progressive_cosetdp(th0, edges, lam_at, H, 64.0, a.iters, coded)
    e_p = error1(decode_pi(wq / 64.0, edges), H)[0]
    wq2, e_f = OR.polish_coded(wq, edges, H, 7, a.sweeps, coded)
    print(f"  COSET form, exact weighted decoder (shared_coset_bits d=11): progressive {1e3 * e_p:.4f} mHa, polished "
          f"{1e3 * e_f:.4f} mHa; coded words {nf}, moved off nearest {nm}; coset violations {check_coset_full(wq2, coded)}"
          f" ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
