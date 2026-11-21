#!/usr/bin/env python
"""chem_lens_project: put further XOR-of-five relations on an EXISTING table of angle words (CPU only, no GPU).

Why this works (research/ws_chem_lens.md 6a): a relation "the XOR of the five words of a batch at bit b is an
affine (or degree-2) function of the low address inside every 16-address block" can be met by changing ONE word
in each of at most 11 (5) rows of the block, and the weakest word of a row hardly matters: in the two-body rows
those entries hold 4e-5 of the table's rounding sensitivity.  So the two-body rows need no new rounding; the 76
one-body rows (blocks 0-4) must carry the same relations and are handled with the exact error_1.

    chem_lens_project.py BASE.npz OUT.npz --xcodes 1:1:3-14,2:1:3-14 [--codes 0:1:3-14]
                         [--base-xcodes 1:2:3-7,1:1:8-14,2:2:8-14] [--polish 6] [--w2 30] [--room]

BASE.npz / OUT.npz: complete submission tables (angle_words, tree_edges, U, W, H2_I, H1_I, Hsym_lower).  Words keep
bit 0 where a per-column code holds, keep bits 6 and 7 always (the loaders' carriers), and only move inside
their 64-word window.  Steps:
  1. two-body rows of blocks 5..: per block and batch the codewords are chosen by coordinate descent on
     cost = sensitivity x (word change)^2, the change going to the cheapest of the five words of a row;
  2. one-body rows (blocks 0-3) and block 4 (12 one-body + 4 two-body rows): the same choice with the EXACT change of
     error_1^2 = || V^T diag(lam) V - h1' ||^2 for every candidate (weights fixed at the stored ones);
  3. polish of the one-body rows: single-word moves that keep every coded bit and pair moves inside a batch whose
     flips cancel in every XOR relation, accepted when the exact error_1 falls;
  4. report: both errors at the stored weights from chem_own_refit.Problem (the hosts' formulas), optionally the
     minimum of error_2 over the weights (--room).  The weights are NOT refitted here: tools/sh/chem_own_build.sh does.
Credits: Reed-Muller codewords and checks from the public pipeline's rmcode.py (marculera, Apache-2.0, read-only
import); metric formulas through tools/chem_own_refit.py (ours).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

os.environ.setdefault("JAX_PLATFORMS", "cpu")
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem_factor_lib as L  # noqa: E402
from chem_lens_orbit import angles_and_mass  # noqa: E402
import chem_lens_codes as cc  # noqa: E402
import chem_own_refit as cor  # noqa: E402  (read-only import: Problem)

rmcode = cc.rmcode
N = 76
_NEAR: dict = {}


def nearest_table(codedmask: int, bits: tuple):
    """d1, d2[f, v + 128]: the smallest word changes (two on a tie) that flip exactly the coded bits of pattern f
    (bit j of f <-> bits[j]), keep every other bit of `codedmask`, and stay inside v's 64-word window."""
    key = (codedmask, bits)
    if key in _NEAR:
        return _NEAR[key]
    npat = 1 << len(bits)
    d1 = np.zeros((npat, 256), np.int64)
    d2 = np.zeros((npat, 256), np.int64)
    for f in range(npat):
        fm = sum(((f >> j) & 1) << b for j, b in enumerate(bits))
        for v in range(-128, 128):
            lo = v - (v % 64)
            best = []
            for dist in range(64):
                for dd in ((0,) if dist == 0 else (dist, -dist)):
                    u = v + dd
                    if lo <= u <= lo + 63 and ((u ^ v) & codedmask) == fm:
                        best.append(dd)
                if best:
                    break
            assert best, (codedmask, bits, f, v)
            d1[f, v + 128], d2[f, v + 128] = best[0], best[-1]
    _NEAR[key] = (d1, d2)
    return d1, d2


def choose_codewords(x, posr, degs, cost):
    """x (nbits, nrows): XOR-of-five bits now; posr: low addresses of the rows; degs: code degree per bit;
    cost (nrows, npat): cost of flip pattern f in each row (cost[:, 0] = 0).  Returns the pattern per row that makes
    every bit a codeword at the least summed cost (coordinate descent over the bits, exhaustive in the codewords)."""
    nb, nr = x.shape
    ar = np.arange(nr)
    cands = [np.unique(rmcode.code(int(d))[:, posr].astype(np.int64) ^ x[j][None, :], axis=0) for j, d in enumerate(degs)]
    m = np.zeros((nb, nr), np.int64)
    for j in range(nb):
        m[j] = cands[j][int(np.argmin(cands[j] @ cost[:, 1 << j]))]
    for _ in range(10):
        changed = False
        for j in range(nb):
            others = np.zeros(nr, np.int64)
            for i in range(nb):
                if i != j:
                    others |= m[i] << i
            delta = cost[ar, others | (1 << j)] - cost[ar, others]
            tot = cands[j] @ delta
            b = int(np.argmin(tot))
            if tot[b] < float(m[j] @ delta) - 1e-18:
                m[j], changed = cands[j][b], True
        if not changed:
            break
    f = np.zeros(nr, np.int64)
    for j in range(nb):
        f |= m[j] << j
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base")
    ap.add_argument("out")
    ap.add_argument("--codes", default="0:1:3-14")
    ap.add_argument("--base-xcodes", default="1:2:3-7,1:1:8-14,2:2:8-14")
    ap.add_argument("--xcodes", required=True)
    ap.add_argument("--polish", type=int, default=6)
    ap.add_argument("--w2", type=float, default=30.0, help="weight of a two-body mHa^2 against a one-body mHa^2 in block 4")
    ap.add_argument("--room", action="store_true", help="also minimise error_2 over the weights (about a minute)")
    ap.add_argument("--proxy2", action="store_true", help="two-body costs from sensitivity x change^2 (the first version)")
    ap.add_argument("--polish2", type=int, default=0, help="sweeps of the two-body polish (exact, sequential; 0 = none)")
    ap.add_argument("--per-row", type=int, default=6, help="accepted moves per row and sweep at most")
    ap.add_argument("--tol2", type=float, default=1e-4, help="smallest gain of a two-body move, in mHa^2")
    a = ap.parse_args()
    t0 = time.time()

    def log(*x):
        print(*x, f"({time.time() - t0:.0f} s)", flush=True)

    z = cor.load_table(os.path.expanduser(a.base))
    words0 = np.array(z["angle_words"]).astype(np.int64)
    edges = np.array(z["tree_edges"]).astype(np.int64)
    w0, H2 = np.array(z["W"], np.float64), np.array(z["H2_I"], np.float64)
    W0 = H2 - w0.sum(axis=1)
    R, B, C = w0.shape
    s = ((words0 + 128) % 256) - 128                                       # signed words, (1220, 75)
    addr = np.arange(s.shape[0])
    blk, pos = addr // 16, addr % 16
    deg = rmcode.parse(a.codes, N - 1, 8)
    xbase, xnew = cc.parse_x(a.base_xcodes), cc.parse_x(a.xcodes)
    rmcode.check(words0, addr, deg)
    cc.check_x(words0, addr, xbase)
    prob = cor.Problem(words0, edges, (R, B, C))
    (e2, e1, gap, Lam), _ = prob.evaluate(w0, W0, rematch=True)
    log(f"base {a.base}: {(R, B, C)}, codes [{a.codes}] + [{a.base_xcodes}] hold; stored weights: error_2 {1e3 * e2:.4f}, "
        f"error_1 {1e3 * e1:.4f} mHa, lambda_eff {np.sqrt(gap * (2 * Lam - gap)):.5f}")

    # ---------------------------------------------------------------- two-body rows, blocks 5..
    mass = angles_and_mass(prob.u.reshape(-1, N), edges)[1]
    c2 = np.sum(w0 ** 2, axis=2).reshape(-1)[:, None] * mass + 1e-300     # (R*B, 75)
    kappa = 200e-6 / (c2.sum() / 12.0)                                    # raw units -> Ha^2 (plain rounding ~ 200 mHa^2)
    c2 = c2 * kappa
    batches = [k for k in range(15) if xnew[k].any()]
    info = {}
    for k in batches:
        bits = tuple(int(b) for b in np.flatnonzero(xnew[k]))
        cols = np.arange(5 * k, 5 * k + 5)
        codedmask = sum(1 << b for b in bits) | (1 if deg[cols[0], 0] > 0 else 0)
        info[k] = (bits, cols, codedmask, nearest_table(codedmask, bits), [int(xnew[k, b]) for b in bits])
    def decode(rows_words):
        v = L.decode_np(np.asarray(rows_words, np.float64) / 128.0, edges)
        return v / np.linalg.norm(v, axis=1, keepdims=True)

    # Exact model of error_2^2 at fixed weights, kept current after every accepted change:
    #   M = sum_rc a_rc a_rc^T,  a_rc = sum_b w_rbc u_rb (x) u_rb,  Res = M - T,  g_i = sum_c w_ic a_rc,  y_i = Res g_i
    # Row i takes u': d = u'(x)u' - u(x)u,  dM = d g_i^T + g_i d^T + |w_i|^2 d d^T,  dE = 2 <Res, dM> + |dM|^2  (exact).
    # The first version took every change's first-order gain separately and lost: 4,730 changes predicted -7.9
    # mHa^2 and gave +34.6, because they all used the same residual.
    E_now = 0.0
    if not a.proxy2:
        dummy = {"U": prob.u, "W": np.zeros((R, B, C)), "H2_I": np.zeros((R, C)), **prob.shifts}
        ref = prob.V_mod.load_inputs(dummy).ref
        Tm = np.asarray(ref["tei"], np.float64).reshape(N * N, N * N).copy()
        sv, iv = prob.S.reshape(-1), np.eye(N).reshape(-1)
        Tm += 0.5 * (np.outer(sv, iv) + np.outer(iv, sv))
        Amat = np.einsum("rbc,rbp,rbq->rcpq", w0, prob.u, prob.u, optimize=True).reshape(R * C, N * N)
        Gmat = np.einsum("rbc,rcx->rbx", w0, Amat.reshape(R, C, N * N), optimize=True).reshape(R * B, N * N)
        Ymat = (Gmat @ Amat.T) @ Amat - (Tm @ Gmat.T).T
        E_now = float(np.sum((Amat.T @ Amat - Tm) ** 2))
        g2 = np.sum(Gmat ** 2, axis=1)
        wflat = w0.reshape(R * B, C)
        log(f"two-body model ready: error_2 from the explicit residual {1e3 * np.sqrt(E_now):.4f} mHa")

    def delta_e2(i2, cand_words):
        """Change of error_2^2 (Ha^2) when two-body row i2 takes each row of cand_words (third-order terms left out)."""
        vc = decode(cand_words)
        u = decode(s[N + i2][None, :])[0]
        Y, G = Ymat[i2].reshape(N, N), Gmat[i2].reshape(N, N)
        dy = np.einsum("np,np->n", vc @ Y, vc) - u @ Y @ u
        dg = np.einsum("np,np->n", vc @ G, vc) - u @ G @ u
        dd = 2.0 - 2.0 * (vc @ u) ** 2
        return 4.0 * dy + 2.0 * dd * g2[i2] + 2.0 * dg ** 2

    def apply2(i2, new_row):
        """Commit a change of two-body row i2 and bring the model up to date; returns the exact change of error_2^2."""
        nonlocal E_now
        u = decode(s[N + i2][None, :])[0]
        un = decode(new_row[None, :])[0]
        s[N + i2] = new_row
        if a.proxy2:
            return 0.0
        d = (np.outer(un, un) - np.outer(u, u)).reshape(-1)
        r = i2 // B
        grp = slice(r * B, (r + 1) * B)
        wi = wflat[i2]
        gi = Gmat[i2].copy()
        Rd = Amat.T @ (Amat @ d) - Tm @ d
        dgj, gig = Gmat @ d, Gmat @ gi
        w2, dd, gid = float(wi @ wi), float(d @ d), float(gi @ d)
        dE = (4.0 * float(d @ Ymat[i2]) + 2.0 * w2 * float(d @ Rd) + 2.0 * dd * float(gi @ gi) + 2.0 * gid ** 2
              + 4.0 * w2 * dd * gid + w2 * w2 * dd * dd)
        Ymat[:] += np.outer(gig + w2 * dgj, d) + np.outer(dgj, gi)
        delta = wflat[grp] @ wi
        Ymat[grp] += np.outer(delta, Rd + gid * d + dd * gi + w2 * dd * d)
        Gmat[grp] += np.outer(delta, d)
        Amat[r * C:(r + 1) * C] += np.outer(wi, d)
        g2[grp] = np.sum(Gmat[grp] ** 2, axis=1)
        E_now += dE
        return dE

    def unit_costs2(r, cols, npat, d1, d2):
        """Cheapest way for two-body row r (address) to carry each flip pattern of this batch."""
        vals = s[r, cols]
        cost, pick = np.zeros(npat), np.zeros((npat, 2), np.int64)
        if a.proxy2:
            for f in range(1, npat):
                d = d1[f, vals + 128]
                cst = c2[r - N, cols] * d ** 2
                j = int(np.argmin(cst))
                cost[f], pick[f] = cst[j], (j, vals[j] + d[j])
            return cost, pick
        cand, tag = [], []
        for f in range(1, npat):
            for j in range(5):
                for d in {int(d1[f, vals[j] + 128]), int(d2[f, vals[j] + 128])}:
                    row = s[r].copy()
                    row[cols[j]] += d
                    cand.append(row)
                    tag.append((f, j, vals[j] + d))
        de = delta_e2(r - N, np.array(cand))
        cost[1:] = np.inf
        for (f, j, v), dv in zip(tag, de):
            if dv < cost[f]:
                cost[f], pick[f] = dv, (j, v)
        return cost, pick

    raw2, changed2 = 0.0, 0
    for k in batches:
        bits, cols, codedmask, (d1, d2), degs = info[k]
        npat = 1 << len(bits)
        for b in np.unique(blk):
            rows = np.flatnonzero(blk == b)
            if b <= 4:
                continue                                                   # blocks 0-4 hold one-body rows: step 2
            nr = len(rows)
            cost = np.zeros((nr, npat))
            pick = np.zeros((nr, npat, 2), np.int64)
            for i, r in enumerate(rows):
                cost[i], pick[i] = unit_costs2(r, cols, npat, d1, d2)
            x = np.stack([np.bitwise_xor.reduce((s[rows][:, cols] >> bit) & 1, axis=1) for bit in bits])
            # a flip is never counted as a gain when the codeword is chosen (gains are taken by the polish, one by one)
            safe = np.where(np.arange(npat)[None, :] > 0, np.maximum(cost, 0.0) + 1e-13, 0.0)
            f = choose_codewords(x, pos[rows], degs, safe)
            for i in np.flatnonzero(f):
                row = s[rows[i]].copy()
                row[cols[pick[i, f[i], 0]]] = pick[i, f[i], 1]
                raw2 += apply2(rows[i] - N, row)
                changed2 += 1
    log(f"two-body rows (blocks 5..): {changed2} words changed for the relations, change of error_2^2 {1e6 * raw2:+.3f} mHa^2 "
        f"({'not tracked (proxy costs)' if a.proxy2 else 'exact, sequential'})")

    # polish of the two-body rows at fixed weights: the compensation a rounding would have given (exact, sequential)
    cmask_c = np.array([(1 if deg[c, 0] > 0 else 0) | sum(1 << int(b_) for b_ in np.flatnonzero(xnew[c // 5]))
                        for c in range(N - 1)], np.int64)
    xmask_c = np.array([sum(1 << int(b_) for b_ in np.flatnonzero(xnew[c // 5])) for c in range(N - 1)], np.int64)

    def moves(v, span, single_free, single_coded):
        """Candidate rows around word row v that keep every code: single moves (uncoded columns: |change| <= single_free;
        coded columns: changes that flip no coded bit, up to single_coded) and pair moves inside a batch whose flips
        cancel in every XOR relation (|change| <= span)."""
        top = max(span, single_free, single_coded)
        dl = np.array([d for d in range(-top, top + 1) if d != 0], np.int64)
        u = v[:, None] + dl[None, :]
        ok = (u >= (v - (v % 64))[:, None]) & (u <= (v - (v % 64) + 63)[:, None])
        fm = (u ^ v[:, None]) & cmask_c[:, None]
        cand = []
        lim = np.where(cmask_c > 0, single_coded, single_free)[:, None]
        for c, d in zip(*np.nonzero(ok & (fm == 0) & (np.abs(dl)[None, :] <= lim))):
            row = v.copy()
            row[c] = u[c, d]
            cand.append(row)
        small = np.abs(dl) <= span
        for k in batches:
            cols = np.arange(5 * k, 5 * k + 5)
            good = [ok[c] & small & ((fm[c] & ~xmask_c[c]) == 0) & (fm[c] != 0) for c in cols]
            for i1 in range(5):
                for i2_ in range(i1 + 1, 5):
                    c1, c2_ = cols[i1], cols[i2_]
                    eq = (fm[c1][:, None] == fm[c2_][None, :]) & good[i1][:, None] & good[i2_][None, :]
                    for p_, q_ in zip(*np.nonzero(eq)):
                        row = v.copy()
                        row[c1], row[c2_] = u[c1, p_], u[c2_, q_]
                        cand.append(row)
        return np.array(cand) if cand else np.zeros((0, N - 1), np.int64)

    if not a.proxy2:
        for sweep in range(a.polish2):
            nmove, gain, t_s = 0, 0.0, time.time()
            for i2 in range(R * B):
                for _ in range(a.per_row):
                    cand = moves(s[N + i2], 4, 2, 8)
                    if not len(cand):
                        break
                    de = delta_e2(i2, cand)
                    best = int(np.argmin(de))
                    if de[best] >= -a.tol2 * 1e-6:
                        break
                    gain += apply2(i2, cand[best])
                    nmove += 1
            log(f"two-body polish sweep {sweep}: {nmove} moves, change of error_2^2 {1e6 * gain:+.3f} mHa^2, "
                f"error_2 now {1e3 * np.sqrt(E_now):.4f} mHa at the stored weights ({time.time() - t_s:.0f} s for the sweep)")
            if nmove == 0:
                break
        del Tm

    # ---------------------------------------------------------------- one-body rows and block 4
    prob = cor.Problem(np.mod(s, 256), edges, (R, B, C))
    (e2, e1, gap, Lam), _ = prob.evaluate(w0, W0, rematch=True)
    log(f"after the two-body step, stored weights: error_2 {1e3 * e2:.4f}, error_1 {1e3 * e1:.4f} mHa")
    kk = np.einsum("rc,rbc->rb", W0, w0)
    H = prob.h1_base - np.einsum("rb,rbp,rbq->pq", kk, prob.u, prob.u, optimize=True)
    lam_all = np.linalg.eigvalsh(H)
    lam = lam_all[prob.perm]                                               # eigenvalue of every one-body row
    V = prob.Vrows.copy()
    D = (V.T * lam) @ V - H

    assert np.allclose(np.abs(np.sum(decode(s[:N].astype(np.float64)) * V, axis=1)), 1.0, atol=1e-9), "decoder differs from the hosts'"
    V = decode(s[:N].astype(np.float64))
    D = (V.T * lam) @ V - H
    log(f"one-body rows: error_1 {1e3 * np.linalg.norm(D):.4f} mHa at the stored weights (own decoder)")

    def delta_e1(r, cand_words):
        """Exact change of error_1^2 when one-body row r is replaced by each row of cand_words."""
        vc = decode(cand_words.astype(np.float64))
        q = np.einsum("np,np->n", vc @ D, vc) - V[r] @ D @ V[r]
        return 2.0 * lam[r] * q + lam[r] ** 2 * (2.0 - 2.0 * (vc @ V[r]) ** 2), vc

    def commit(r, new_row):
        nonlocal D
        vn = decode(new_row[None, :].astype(np.float64))[0]
        D = D + lam[r] * (np.outer(vn, vn) - np.outer(V[r], V[r]))
        V[r] = vn
        s[r] = new_row

    changed1 = 0
    for k in batches:
        bits, cols, codedmask, (d1, d2), degs = info[k]
        npat = 1 << len(bits)
        for b in range(5):
            rows = np.flatnonzero(blk == b)
            nr = len(rows)
            cost = np.zeros((nr, npat))
            pick = np.zeros((nr, npat, 2), np.int64)
            for i, r in enumerate(rows):
                vals = s[r, cols]
                if r >= N:                                                 # the four two-body rows of block 4
                    cost[i], pick[i] = unit_costs2(r, cols, npat, d1, d2)
                    cost[i] *= a.w2
                    continue
                cand, tag = [], []
                for f in range(1, npat):
                    for j in range(5):
                        for d in {int(d1[f, vals[j] + 128]), int(d2[f, vals[j] + 128])}:
                            row = s[r].copy()
                            row[cols[j]] += d
                            cand.append(row)
                            tag.append((f, j, vals[j] + d))
                de, _ = delta_e1(r, np.array(cand))
                cost[i, 1:] = np.inf
                for (f, j, v), dv in zip(tag, de):
                    if dv < cost[i, f]:
                        cost[i, f], pick[i, f] = dv, (j, v)
            x = np.stack([np.bitwise_xor.reduce((s[rows][:, cols] >> bit) & 1, axis=1) for bit in bits])
            f = choose_codewords(x, pos[rows], degs, cost)
            for i in np.flatnonzero(f):
                r = rows[i]
                if r >= N:
                    s[r, cols[pick[i, f[i], 0]]] = pick[i, f[i], 1]
                else:
                    row = s[r].copy()
                    row[cols[pick[i, f[i], 0]]] = pick[i, f[i], 1]
                    commit(r, row)
                changed1 += 1
    log(f"blocks 0-4: {changed1} words changed; error_1 {1e3 * np.linalg.norm(D):.4f} mHa before the polish")

    # ---------------------------------------------------------------- polish of the one-body rows
    cmask = np.zeros(N - 1, np.int64)
    for c in range(N - 1):
        cmask[c] = (1 if deg[c, 0] > 0 else 0) | sum(1 << int(b) for b in np.flatnonzero(xnew[c // 5]))
    xmask = np.array([sum(1 << int(b) for b in np.flatnonzero(xnew[c // 5])) for c in range(N - 1)], np.int64)
    dl = np.array([d for d in range(-16, 17) if d != 0], np.int64)
    for sweep in range(a.polish):
        moves = 0
        for r in range(N):
            for _ in range(200):
                v = s[r]
                u = v[:, None] + dl[None, :]                                # (75, 32) candidate words
                ok = (u >= (v - (v % 64))[:, None]) & (u <= (v - (v % 64) + 63)[:, None])
                fm = (u ^ v[:, None]) & cmask[:, None]                      # flipped coded bits
                cand = []
                ci, di = np.nonzero(ok & (fm == 0) & (np.abs(dl)[None, :] <= np.where(cmask > 0, 16, 3)[:, None]))
                for c, d in zip(ci, di):                                    # single moves that keep every coded bit
                    row = v.copy()
                    row[c] = u[c, d]
                    cand.append(row)
                small = np.abs(dl) <= 8
                for k in batches:                                           # pair moves: flips cancel in the XOR relations
                    cols = np.arange(5 * k, 5 * k + 5)
                    for i1 in range(5):
                        for i2 in range(i1 + 1, 5):
                            c1, c2_ = cols[i1], cols[i2]
                            g1 = ok[c1] & small & ((fm[c1] & ~xmask[c1]) == 0) & (fm[c1] != 0)
                            g2 = ok[c2_] & small & ((fm[c2_] & ~xmask[c2_]) == 0) & (fm[c2_] != 0)
                            eq = (fm[c1][:, None] == fm[c2_][None, :]) & g1[:, None] & g2[None, :]
                            for p, q in zip(*np.nonzero(eq)):
                                row = v.copy()
                                row[c1], row[c2_] = u[c1, p], u[c2_, q]
                                cand.append(row)
                if not cand:
                    break
                cand = np.array(cand)
                de, _ = delta_e1(r, cand)
                best = int(np.argmin(de))
                if de[best] >= -1e-13:
                    break
                commit(r, cand[best])
                moves += 1
        log(f"polish sweep {sweep}: {moves} moves, error_1 {1e3 * np.linalg.norm(D):.4f} mHa")
        if moves == 0:
            break

    # ---------------------------------------------------------------- checks, report, save
    table = np.mod(s, 256).astype(np.int64)
    rmcode.check(table, addr, deg)
    cc.check_x(table, addr, xnew)
    assert np.array_equal(table >> 6, words0 >> 6), "a word left its 64-word window"
    prob = cor.Problem(table, edges, (R, B, C))
    (e2, e1, gap, Lam), _ = prob.evaluate(w0, W0, rematch=True)
    off = prob.official(w0, W0)
    rep = {"base": a.base, "codes": a.codes, "xcodes": a.xcodes, "words_changed": int(np.sum(table != words0)),
           "words_changed_two_body": int(np.sum(table[N:] != words0[N:])), "raw_two_body_mHa2": 1e6 * raw2,
           "stored_weights": {"error_2_mHa": 1e3 * e2, "error_1_mHa": 1e3 * e1, "lambda_eff": float(np.sqrt(gap * (2 * Lam - gap)))},
           "official_stored_weights": {k_: float(v_) for k_, v_ in off.items()}}
    log(f"new table: codes [{a.codes}] + [{a.xcodes}] hold on every block; {rep['words_changed']} words differ "
        f"({rep['words_changed_two_body']} two-body); stored weights: error_2 {1e3 * e2:.4f}, error_1 {1e3 * e1:.4f} mHa; "
        f"hosts' evaluate_words: error_2 {1e3 * off['error_2']:.4f}, error_1 {1e3 * off['error_1']:.4f}")
    if a.room:
        from scipy.optimize import minimize

        def fobj(xv):
            (e2_, _, _, _), (de2, _, _, _) = prob.evaluate(xv.reshape(w0.shape), W0)
            return e2_ * 1e3, de2.reshape(-1) * 1e3
        res = minimize(fobj, w0.reshape(-1), jac=True, method="L-BFGS-B", options={"maxiter": 600, "maxfun": 1200, "maxcor": 30})
        rep["min_error_2_mHa"] = float(res.fun)
        log(f"minimum of error_2 over the weights: {res.fun:.4f} mHa (room to 19.9985: {19.9985 ** 2 - res.fun ** 2:.2f} mHa^2)")
    out = {k_: np.array(v_) for k_, v_ in z.items()}
    out["angle_words"] = table
    out["U"] = prob.u
    np.savez(os.path.expanduser(a.out), **out)
    json.dump(rep, open(os.path.splitext(os.path.expanduser(a.out))[0] + ".json", "w"), indent=1)
    log(f"saved {a.out}")


if __name__ == "__main__":
    main()
