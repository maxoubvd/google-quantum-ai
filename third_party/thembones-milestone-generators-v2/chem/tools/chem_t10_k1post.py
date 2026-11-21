"""chem_t10_k1post.py (T10 seat 3, Nov 18): KILL TEST K1 of ledger row O14.

Bit-1 relations put on a FINISHED, GAUGED table (all_c.npz: C2's carrier gauge applied, carriers exact) by the
cheapest repairs, then the weights re-solved with the words fixed.

Per relation (loader k, outputs (j, 1); output (j, b) of loader k = bit b of w[5k+j] XOR w[5k-5+j]): the rows'
present parity on the gauged words, the best target f(addr & 15) XOR g(addr >> 4) (exhaustive f, majority g, as
C3's chem_c3_relfit.py). Per row the broken relations are repaired JOINTLY by the cheapest set (1 to 3) of +-2 moves
of support words (a +-2 move flips bit 1 and keeps bit 0; a move that changes the carrier bit6^bit5 is refused; a
move of a batch-k word also toggles loader k+1's relations on the same slot: every toggle is accounted). Prices are
EXACT at fixed weights: two-body rows by Delta error_2^2 (incremental K / Q update), one-body rows by Delta error_1^2.
Then the weights w minimise error_2 (L-BFGS, exact gradient of r48_chem_q7c_refit.Problem, words fixed): the
END-like reading to set against the wall 399.9 mHa^2; the same re-solve on the unmodified table gives the base.

    python chem_t10_k1post.py all_c.npz --layout rel.txt [--out OUT.npz] [--brot 7] [--iters 400]
"""
import argparse
import itertools
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r48_chem_q7c_refit as RF  # noqa: E402  (read-only import)


def log(*x):
    print(*x, flush=True)


def read_layout(path):
    """bit-1 (and higher) relation lines of an RK7-format layout file -> [(k, [(j, b), ...])]."""
    rels = []
    for line in open(path):
        line = line.split("#")[0].strip()
        if not line or "|" not in line:
            continue
        k, outs, _t = [x.strip() for x in line.split("|")]
        o = sorted(tuple(int(v) for v in t.split(":")) for t in outs.split())
        if all(b == 0 for (_, b) in o):
            continue
        rels.append((int(k), o))
    return rels


def parity(w, k, outs):
    p = np.zeros(len(w), np.int64)
    for (j, b) in outs:
        p ^= ((w[:, 5 * k + j] >> b) & 1) ^ ((w[:, 5 * k - 5 + j] >> b) & 1)
    return p


def fit(p, addr):
    """best f (16 bits, exhaustive) and g (majority per block) for parity p (C3's chem_c3_relfit.fit)."""
    lo, blk = addr & 15, addr >> 4
    nb = int(blk.max()) + 1
    fs = ((np.arange(1 << 16)[:, None] >> np.arange(16)[None, :]) & 1).astype(np.int8)
    q = (p[None, :].astype(np.int8) ^ fs[:, lo]).astype(np.int16)
    ones = np.zeros((fs.shape[0], nb), np.int32)
    for b in range(nb):
        ones[:, b] = q[:, blk == b].sum(axis=1)
    nrow = np.bincount(blk, minlength=nb)
    mis = np.minimum(ones, nrow[None, :] - ones).sum(axis=1)
    i = int(np.argmin(mis))
    g = (ones[i] * 2 > nrow).astype(np.int64)
    return int(mis[i]), fs[i][lo].astype(np.int64) ^ g[blk]


class Pricer:
    """Exact Delta error_2^2 of one two-body row's new vector at fixed weights; Delta error_1^2 of a one-body row."""

    def __init__(self, words, edges, shape, brot, w, W):
        from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
        self.rec = reconstruct_u_vectors
        self.brot = brot
        self.edges = edges
        P = RF.Problem(words, edges, shape, brot)       # K, Q, T2, u, Vrows, h1 parts
        self.P = P
        R, B, C = shape
        N = P.N
        self.R, self.B, self.C, self.N = R, B, C, N
        # the two-body target as an N^2 x N^2 matrix, rebuilt (Problem drops it)
        from dfthcsos import verify as V
        dummy = {"U": P.u, "W": np.zeros((R, B, C)), "H2_I": np.zeros((R, C)), **P.shifts}
        ref = V.load_inputs(dummy).ref
        Tm = np.asarray(ref["tei"], np.float64).reshape(N * N, N * N).copy()
        sv, iv = P.S.reshape(-1), np.eye(N).reshape(-1)
        Tm += 0.5 * (np.outer(sv, iv) + np.outer(iv, sv))
        self.Tm = Tm
        self.flat = P.u.reshape(R * B, N).copy()
        self.uu = (self.flat[:, :, None] * self.flat[:, None, :]).reshape(R * B, N * N)
        self.w = w
        Wf = np.zeros((R * B, R * C))
        for r in range(R):
            Wf[r * B:(r + 1) * B, r * C:(r + 1) * C] = w[r]
        self.Wf = Wf
        self.A = Wf.T @ P.K @ Wf
        self.wn = np.einsum("xc,yc->xy", w.reshape(R * B, C), w.reshape(R * B, C))   # (RB, RB) w_x . w_y
        # one-body side at fixed weights
        k = np.einsum("rc,rbc->rb", W, w)
        H = P.h1_base - np.einsum("rb,rbp,rbq->pq", k, P.u, P.u, optimize=True)
        lam, vec = np.linalg.eigh(H)
        from scipy.optimize import linear_sum_assignment
        rows, cols = linear_sum_assignment(-np.abs(P.Vrows @ vec))
        perm = np.empty(N, np.int64)
        perm[rows] = cols
        self.l1 = lam[perm]
        self.D = (P.Vrows.T * self.l1) @ P.Vrows - H
        self.V = P.Vrows.copy()

    def vec(self, wordrows):
        v = self.rec(np.asarray(wordrows, np.float64) / (1 << (self.brot - 1)), self.edges)
        return v / np.linalg.norm(v, axis=1, keepdims=True)

    def dE2(self, xs, U2):
        """xs: two-body row indices (0..RB-1), U2: their new unit vectors (n, N) -> exact Delta (error_2^2)."""
        R, B, C, N = self.R, self.B, self.C, self.N
        K = self.P.K
        n = len(xs)
        out = np.zeros(n)
        for s in range(0, n, 256):
            xi, Ui = xs[s:s + 256], U2[s:s + 256]
            dk = (Ui @ self.flat.T) ** 2 - K[xi]
            dk[np.arange(len(xi)), xi] = 0.0
            bv = dk @ self.Wf
            av = self.Wf[xi]
            t1 = 4 * np.einsum("ni,ij,nj->n", av, self.A, bv)
            t2 = 2 * ((av * av).sum(1) * (bv * bv).sum(1) + (av * bv).sum(1) ** 2)
            uu2 = (Ui[:, :, None] * Ui[:, None, :]).reshape(len(xi), N * N)
            tp = uu2 @ self.Tm                                           # (n, N^2)
            dq = np.zeros(len(xi))
            for i, x in enumerate(xi):
                r, b = divmod(int(x), B)
                newrow = self.uu[r * B:(r + 1) * B] @ tp[i]
                newrow[b] = uu2[i] @ tp[i]
                d = newrow - self.P.Q[r, b]
                wx = self.wn[x, r * B:(r + 1) * B]
                dq[i] = 2 * (d * wx).sum() - d[b] * wx[b]
            out[s:s + 256] = t1 + t2 - 2 * dq
        return out

    def dE1(self, i, v2):
        """one-body row i -> new vector v2: Delta (error_1^2) at fixed weights and matching."""
        v = self.V[i]
        d = self.l1[i] * (np.outer(v2, v2) - np.outer(v, v))
        return float(np.sum((self.D + d) ** 2) - np.sum(self.D ** 2))


def min_e2(prob, w, W, iters):
    """minimise error_2 over the weights w (words fixed): L-BFGS with the exact gradient."""
    sc = max(float(np.std(w)), 1e-9)

    def f(z):
        ww = z.reshape(w.shape) * sc
        (e2, e1, gap, Lam), (de2, de1, dgap, dLam) = prob.evaluate(ww, W)
        return 1e6 * e2 * e2, (1e6 * 2 * e2 * de2 * sc).ravel()

    res = minimize(f, (w / sc).ravel(), jac=True, method="L-BFGS-B",
                   options={"maxiter": iters, "maxcor": 30, "ftol": 1e-15, "gtol": 1e-10, "maxls": 40})
    w2 = res.x.reshape(w.shape) * sc
    return w2, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--layout", required=True)
    ap.add_argument("--out", default="")
    ap.add_argument("--brot", type=int, default=7)
    ap.add_argument("--iters", type=int, default=400)
    ap.add_argument("--base-only", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    z = np.load(a.table)
    d = {k: np.array(z[k]) for k in z.files}
    w0 = np.mod(np.asarray(d["angle_words"], np.int64), 1 << a.brot)
    edges = np.asarray(d["tree_edges"], np.int64)
    wt = np.asarray(d["W"], np.float64)
    shape = wt.shape
    Woff = np.asarray(d["H2_I"], np.float64) - wt.sum(axis=1)
    R, B, C = shape
    nrow = w0.shape[0]
    N1 = edges.shape[0] + 1
    addr = np.arange(nrow)
    car0 = ((w0 >> 6) ^ (w0 >> 5)) & 1
    rels = read_layout(a.layout)
    log(f"table {a.table}: {nrow} rows, shape {shape}; {len(rels)} bit-1 relations {[(k, len(o)) for k, o in rels]}")

    # ---- base: the official reading and the min-error_2 re-solve of the unmodified table
    prob0 = RF.Problem(w0, edges, shape, a.brot)
    m0 = prob0.official(wt, Woff)
    log(f"BASE official: e2 {1e3 * m0['error_2']:.5f} mHa ({1e6 * m0['error_2'] ** 2:.2f} mHa^2) e1 "
        f"{1e3 * m0['error_1']:.4f} lambda_eff {m0['lambda_eff']:.6f} ({time.time() - t0:.0f} s)")
    wb, res = min_e2(prob0, wt, Woff, a.iters)
    mb = prob0.official(wb, Woff)
    log(f"BASE min-e2 re-solve: e2^2 {1e6 * mb['error_2'] ** 2:.3f} mHa^2 e1 {1e3 * mb['error_1']:.4f} "
        f"({res.nit} its, {time.time() - t0:.0f} s)")
    if a.base_only:
        return

    # ---- targets and broken rows on the gauged words
    nrel = len(rels)
    targets, broken = [], np.zeros((nrow, nrel), np.int64)
    for i, (k, o) in enumerate(rels):
        p = parity(w0, k, o)
        mis, t = fit(p, addr)
        targets.append(t)
        broken[:, i] = p ^ t
        log(f"  rel {i} k{k} {o}: mismatched after the best target {mis} ({100 * mis / nrow:.1f}%)")

    # toggle mask of a move of column c
    def toggles(c):
        kb, j = divmod(c, 5)
        m = 0
        for i, (k, o) in enumerate(rels):
            if (j, 1) in o and (k == kb or k == kb + 1):
                m |= 1 << i
        return m

    tog = np.array([toggles(c) for c in range(w0.shape[1])], np.int64)
    supp = {}
    for i, (k, o) in enumerate(rels):
        supp[i] = sorted({5 * k + j for (j, b) in o} | {5 * k - 5 + j for (j, b) in o})

    pr = Pricer(w0, edges, shape, a.brot, wt, Woff)
    log(f"pricer built ({time.time() - t0:.0f} s)")
    bad_rows = np.nonzero(broken.any(axis=1))[0]
    # candidates per bad row: support columns of its broken relations, +-2, carriers kept
    cand = []        # (row, col, newword, togglemask)
    for r in bad_rows:
        cols = sorted(set(itertools.chain.from_iterable(supp[i] for i in np.nonzero(broken[r])[0])))
        for c in cols:
            for dlt in (2, -2):
                nw = (w0[r, c] + dlt) % (1 << a.brot)
                if (((nw >> 6) ^ (nw >> 5)) & 1) != car0[r, c]:
                    continue
                cand.append((int(r), int(c), int(nw), int(tog[c])))
    log(f"bad rows {len(bad_rows)} (two-body {int((bad_rows >= N1).sum())}, one-body {int((bad_rows < N1).sum())}); "
        f"candidates {len(cand)}")
    # exact prices
    price = np.zeros(len(cand))
    tb = [i for i, (r, c, nw, m) in enumerate(cand) if r >= N1]
    ob = [i for i, (r, c, nw, m) in enumerate(cand) if r < N1]
    rows2 = np.array([w0[cand[i][0]] for i in tb])
    for jj, i in enumerate(tb):
        rows2[jj, cand[i][1]] = cand[i][2]
    U2 = pr.vec(rows2) if len(tb) else np.zeros((0, pr.N))
    xs = np.array([cand[i][0] - N1 for i in tb], np.int64)
    price[tb] = pr.dE2(xs, U2) * 1e6                     # mHa^2 of error_2^2
    for i in ob:
        r, c, nw, m = cand[i]
        row = w0[r].copy()
        row[c] = nw
        price[i] = pr.dE1(r, pr.vec(row[None, :])[0]) * 1e6
    log(f"priced ({time.time() - t0:.0f} s): two-body moves median {np.median(price[tb]) if tb else 0:.4f} mHa^2")

    # joint cheapest repair per row (1 to 3 moves, prices added)
    by_row = {}
    for i, (r, c, nw, m) in enumerate(cand):
        by_row.setdefault(r, []).append(i)
    w1 = w0.copy()
    tot2, tot1, nmov, unfixed = 0.0, 0.0, 0, 0
    for r in bad_rows:
        need = 0
        for i in np.nonzero(broken[r])[0]:
            need |= 1 << int(i)
        lst = sorted(by_row.get(int(r), []), key=lambda i: price[i])
        best = None
        for i in lst:                                        # one move
            if cand[i][3] == need and (best is None or price[i] < best[0]):
                best = (price[i], (i,))
        if best is None:                                     # two moves on different columns
            for x in range(len(lst)):
                for y in range(x + 1, len(lst)):
                    i, j = lst[x], lst[y]
                    if cand[i][1] != cand[j][1] and (cand[i][3] ^ cand[j][3]) == need:
                        cst = price[i] + price[j]
                        if best is None or cst < best[0]:
                            best = (cst, (i, j))
        if best is None:                                     # three moves among the 14 cheapest
            for comb in itertools.combinations(lst[:14], 3):
                if len({cand[i][1] for i in comb}) < 3:
                    continue
                if (cand[comb[0]][3] ^ cand[comb[1]][3] ^ cand[comb[2]][3]) == need:
                    cst = sum(price[i] for i in comb)
                    if best is None or cst < best[0]:
                        best = (cst, comb)
        if best is None:
            unfixed += 1
            continue
        for i in best[1]:
            w1[r, cand[i][1]] = cand[i][2]
            nmov += 1
        if r >= N1:
            tot2 += best[0]
        else:
            tot1 += best[0]
    # checks on the repaired words
    car1 = ((w1 >> 6) ^ (w1 >> 5)) & 1
    assert np.array_equal(car0, car1), "carriers changed"
    assert np.array_equal(w0 & 1, w1 & 1), "bit 0 changed"
    nb = 0
    for i, (k, o) in enumerate(rels):
        nb += int((parity(w1, k, o) ^ targets[i]).sum())
    log(f"REPAIRED: {nmov} moves, rows unfixed {unfixed}, relations broken after {nb}; first-order sums: "
        f"two-body {tot2:.2f} mHa^2 (error_2^2), one-body {tot1:.2f} mHa^2 (error_1^2) ({time.time() - t0:.0f} s)")

    # ---- the repaired table: official at the old weights, then the min-error_2 re-solve
    prob1 = RF.Problem(w1, edges, shape, a.brot)
    m1 = prob1.official(wt, Woff)
    log(f"REPAIRED official at the table's weights: e2^2 {1e6 * m1['error_2'] ** 2:.3f} mHa^2 e1 "
        f"{1e3 * m1['error_1']:.4f}")
    w2, res = min_e2(prob1, wt, Woff, a.iters)
    m2 = prob1.official(w2, Woff)
    ex = 1e6 * (m2["error_2"] ** 2 - mb["error_2"] ** 2)
    log(f"REPAIRED min-e2 re-solve: e2^2 {1e6 * m2['error_2'] ** 2:.3f} mHa^2 (base {1e6 * mb['error_2'] ** 2:.3f}: "
        f"excess {ex:.3f} = {ex / max(nrel, 1):.3f} a relation) e1 {1e3 * m2['error_1']:.4f} ({res.nit} its, "
        f"{time.time() - t0:.0f} s)")
    verdict = "UNDER THE WALL" if 1e6 * m2["error_2"] ** 2 <= 399.9 and m2["error_1"] <= 0.02 else "OVER THE WALL"
    log(f"K1 VERDICT: {verdict} (wall 399.9 mHa^2)")
    if a.out:
        out = dict(d)
        out["angle_words"] = w1
        out["W"] = w2
        out["H2_I"] = Woff + w2.sum(axis=1)
        np.savez(a.out, **out)
        with open(os.path.splitext(a.out)[0] + ".k1.json", "w") as fh:
            json.dump(dict(base=mb, repaired_at_old=m1, repaired=m2, moves=nmov, excess=ex, relations=nrel), fh,
                      indent=1)
        log("saved", a.out)


if __name__ == "__main__":
    main()
