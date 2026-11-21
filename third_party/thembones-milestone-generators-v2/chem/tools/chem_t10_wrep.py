#!/usr/bin/env python
"""chem_t10_wrep (T10 seat 5, Nov 19): THE WIRE-CAPABLE END REPAIR of a relation table.

THE DEFECT IT REPAIRS (seat 5, kill test K0 [M]): a late eps change (a word entering or leaving the pole z = -32 moves
the carrier recursion's n, so eps, of the columns that read it) on a relation whose words are all frozen leaves a row with
no relation-consistent pattern at its place; the placer's backtrack then hands that row's still-free words arbitrary
residues, so the row ends with several relations broken on the STORED words, wires (bit 0) among them.

THE REPAIR, per row with a relation broken on the stored words (rows in address order 76 + i, as the END check reads
them): the placer's own programme (tools/chem_t10_relge.RelPlacerGE._dp, diag, at the row's place) over EVERY relation of
the row (wires and bit 1, on y = eps x z, eps from the row's present words) finds the cheapest residue pattern; a word's
price for each residue is the EXACT change of error_2^2 at fixed weights of the best value of that residue within +-4 of
the present word (R8's one-row screen, as tools/chem_t10_ppol.py). The pattern is applied, the row re-checked on its
stored words and the step repeated (a moved word can change a pole and so eps) until the row breaks nothing. Then the
weights are re-solved (StageLM7, weights only). Rows that do not converge are put back and reported.

    python chem_t10_wrep.py SRC.npz REL_USED_G.txt OUT.npz [--maxit 8] [--span 4] [--wsteps 12]
SRC / OUT have the format of the rounding's rounded.npz / the hold's held.npz (words mod 2^bits, W, U, tree_edges, ...)."""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("R8_LM_WONLY", "1")
from chem_own2_lmround import Target  # noqa: E402  (read-only imports)
import chem_own2_lmstage_w as LW  # noqa: E402
from r65_chem_all15_relg import parse_g  # noqa: E402
from chem_rk7_relplace import BIG  # noqa: E402
from chem_t10_relge import RelPlacerGE, PERM, stored_words  # noqa: E402
from chem_t10_k0rows import per_row  # noqa: E402

N = 76


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("rels")
    ap.add_argument("out")
    ap.add_argument("--maxit", type=int, default=8, help="programme + re-check rounds per row")
    ap.add_argument("--span", type=int, default=4, help="candidate values within +-span of the present word")
    ap.add_argument("--wsteps", type=int, default=12, help="weights-only LM steps of the final re-solve")
    ap.add_argument("--no-keep-pole", dest="keep_pole", action="store_false", help="v1: allow moves onto / off z = -32")
    a = ap.parse_args()
    import jax
    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp

    t0 = time.time()
    z = dict(np.load(a.src))
    bits = int(z["bits"])
    top = 1 << (bits - 1)
    HALF = top // 2
    words = np.asarray(z["words"], np.int64)
    wv = np.mod(words + top, 2 * top) - top                     # signed words (no carriers in these runs)
    assert wv.min() >= -HALF and wv.max() <= HALF - 1, "a word outside the signed range"
    W = np.asarray(z["W"], np.float64)
    R_, B_, C_ = W.shape
    nrows = R_ * B_
    ncol = N - 1
    edges = np.asarray(z["tree_edges"], np.int64)
    RELS = parse_g(open(a.rels).read())
    addr = N + np.arange(nrows)

    def bad_of(wv_, rows):
        rows = np.asarray(rows, np.int64)
        st = stored_words(np.mod(wv_[rows], 2 * top), edges, addr[rows], bits)
        return per_row(st, addr[rows], RELS)

    bad0 = bad_of(wv, np.arange(nrows))
    rows_bad = np.flatnonzero(bad0.sum(1))
    print(f"[wrep] {a.src}: {nrows} rows, {len(RELS)} relations a row; broken on the stored words {int(bad0.sum())} in "
          f"{len(rows_bad)} rows {list(rows_bad)}", flush=True)

    # the placer of the run, with the fitted targets of rel_used_g.txt, every row at its own place (address order)
    fp = RelPlacerGE(R_, B_, RELS, edges, n_one=N, ncols=ncol, mode="exact", brot=bits)
    assert np.array_equal(fp.paddr[fp.place], addr)
    relcols = np.flatnonzero(fp.mod > 1)

    # ---- the exact one-row screen (tools/chem_t10_ppol.py, same formulas)
    refc = Target()
    I, J = np.asarray(refc.i), np.asarray(refc.j)
    m = np.where(I == J, 1.0, 2.0)
    sm = np.sqrt(m)
    Tt = np.outer(sm, sm) * refc.T
    lv = LW._levels(edges)

    def decode(th):
        c, s = jnp.cos(th), jnp.sin(th)
        u = jnp.zeros((th.shape[0], N), th.dtype).at[:, 0].set(1.0)
        for ts in lv:
            p, q = edges[ts, 0], edges[ts, 1]
            up, uq = u[:, p], u[:, q]
            u = u.at[:, p].set(c[:, ts] * up - s[:, ts] * uq).at[:, q].set(s[:, ts] * up + c[:, ts] * uq)
        return u

    smj, Ij, Jj, Ttj = jnp.asarray(sm), jnp.asarray(I), jnp.asarray(J), jnp.asarray(Tt)

    @jax.jit
    def pvec(th):
        u = decode(th)
        return smj * u[:, Ij] * u[:, Jj]

    @jax.jit
    def state(Pr, Wf):
        A = jnp.einsum("rbc,rbp->rcp", Wf, Pr.reshape(R_, B_, -1)).reshape(R_ * C_, -1)
        Rr = A.T @ A - Ttj
        return A, Rr, 0.5 * jnp.sum(Rr * Rr)

    @jax.jit
    def screen(th_rows, Prows, arow, wnorm2, Rr):
        d = pvec(th_rows) - Prows
        Rd = d @ Rr
        ad = jnp.sum(arow * d, axis=1)
        dd = jnp.sum(d * d, axis=1)
        return (2.0 * jnp.sum(arow * Rd, axis=1) + wnorm2 * jnp.sum(d * Rd, axis=1) + jnp.sum(arow * arow, axis=1) * dd
                + ad * ad + 2.0 * wnorm2 * ad * dd + 0.5 * wnorm2 * wnorm2 * dd * dd)

    @jax.jit
    def rank2(Rr, a_, d_, c_):
        return Rr + jnp.outer(a_, d_) + jnp.outer(d_, a_) + c_ * jnp.outer(d_, d_)

    th = jnp.asarray(np.pi * wv / top)
    Pr = pvec(th)
    Wf = jnp.asarray(W)
    A, Rr, f = state(Pr, Wf)
    e2_0 = 2e6 * float(f)
    print(f"[wrep] start error_2^2 {e2_0:.4f} mHa^2 [{time.time() - t0:.0f} s]", flush=True)

    def row_a(r):
        g_, b_ = divmod(int(r), B_)
        w_r = Wf[g_, b_]
        return g_, w_r, A[g_ * C_:(g_ + 1) * C_, :].T @ w_r

    D = np.array([d for d in range(-a.span, a.span + 1) if d != 0], np.int64)
    wvc = wv.copy()
    fixed, failed, cost_tot, moved = 0, 0, 0.0, 0
    for r in rows_bad:
        orig = wvc[r].copy()
        g_, w_r, a_r = row_a(r)
        wn2r = jnp.sum(w_r * w_r)
        ok = False
        for it in range(a.maxit):
            br = np.flatnonzero(bad_of(wvc, [r])[0])
            if len(br) == 0:
                ok = True
                break
            cur = wvc[r].copy()
            # exact deltas of every single-word change d in D on the relation columns (from the row's present words);
            # v2 (seat 5, 05:03Z): no move onto or off the pole z = -32 (with --keep-pole, the default): such a move changes
            # no carrier bit, so the eps of the row stays what the programme assumed [P: R53's recursion]; v1 could flip it
            # and loop (t10_y46a row 1027: 'NOT repaired after 8 rounds' [M])
            cand = [(c, d) for c in relcols for d in D if -HALF <= cur[c] + d <= HALF - 1
                    and (not a.keep_pole or (cur[c] != -HALF and cur[c] + d != -HALF))]
            thr = np.repeat((np.pi * cur / top)[None, :], len(cand), axis=0)
            for i_, (c, d) in enumerate(cand):
                thr[i_, c] += np.pi * d / top
            th_cur = jnp.asarray(np.pi * cur / top)[None, :]
            base = float(screen(th_cur, Pr[r][None, :], a_r[None, :], wn2r[None], Rr)[0])      # the row's present value
            dl = np.asarray(screen(jnp.asarray(thr), jnp.repeat(Pr[r][None, :], len(cand), axis=0),
                                   jnp.repeat(a_r[None, :], len(cand), axis=0), jnp.repeat(wn2r[None], len(cand)), Rr)) - base
            # the placer's word costs by z-residue: the best value of each residue (0 for the present word)
            wc, bestd = {}, {}
            for c in relcols:
                mc = int(fp.mod[c])
                cst = np.full((nrows, 4), BIG)
                bd = np.zeros(4, np.int64)
                cst[r, :] = BIG
                q0 = int(np.mod(cur[c], mc))
                cst[r, q0] = 0.0
                wc[c], bestd[c] = cst, bd
            for i_, (c, d) in enumerate(cand):
                mc = int(fp.mod[c])
                q = int(np.mod(cur[c] + d, mc))
                if q == int(np.mod(cur[c], mc)):
                    continue                     # the present residue keeps its present word (moves of +-4: the polish's)
                if dl[i_] < wc[c][r, q]:
                    wc[c][r, q] = dl[i_]
                    bestd[c][q] = d
            # eps of every (row, low address, column) from the present words; the row's programme at its place
            fp.update_eps(wvc.astype(np.float64), np.ones(wvc.shape, bool), wvc)
            rows1 = np.array([r], np.int64)
            cost, arg = fp._dp(wc, rows1, rows1, diag=True)
            if cost[0] >= 1e17:
                break
            codes = fp._backtrack(arg, rows1)
            L = fp.lowp[rows1]
            new = cur.copy()
            for b in fp.batches:
                x = codes[b]
                for j in range(5):
                    col = 5 * b + j
                    mc = int(fp.mod[col])
                    if mc <= 1:
                        continue
                    c0, c1 = fp.need[b][j]
                    b0 = c0[x] if c0 is not None else None
                    b1 = c1[x] if c1 is not None else None
                    if mc == 2:
                        new[col] = cur[col] + bestd[col][int(b0[0])]
                        continue
                    cst = wc[col][rows1]
                    neg = fp.E[rows1, L, col] < 0
                    cy = np.where(neg[:, None], cst[:, PERM], cst)
                    best_q, best_v = 0, np.inf
                    for q in range(4):
                        if b0 is not None and (q & 1) != int(b0[0]):
                            continue
                        if b1 is not None and ((q >> 1) & 1) != int(b1[0]):
                            continue
                        if cy[0, q] < best_v:
                            best_q, best_v = q, cy[0, q]
                    qz = int(PERM[best_q]) if neg[0] else best_q
                    new[col] = cur[col] + bestd[col][qz]
            wvc[r] = new
        if not ok and len(np.flatnonzero(bad_of(wvc, [r])[0])) == 0:
            ok = True
        if not ok:
            wvc[r] = orig
            failed += 1
            print(f"[wrep] row {r}: NOT repaired after {a.maxit} rounds (put back)", flush=True)
            continue
        thr = np.pi * wvc[r] / top
        exact = float(screen(jnp.asarray(thr[None, :]), Pr[r][None, :], a_r[None, :], wn2r[None], Rr)[0])
        Pnew = pvec(jnp.asarray(thr[None, :]))[0]
        dP = Pnew - Pr[r]
        A = A.at[g_ * C_:(g_ + 1) * C_, :].add(jnp.outer(w_r, dP))
        Pr = Pr.at[r].set(Pnew)
        th = th.at[r].set(jnp.asarray(thr))
        Rr = rank2(Rr, a_r, dP, wn2r)
        f = f + exact
        fixed += 1
        cost_tot += 2e6 * exact
        nmv = int(np.sum(wvc[r] != orig))
        moved += nmv
        print(f"[wrep] row {r}: repaired in {it + 1} round(s), {nmv} words moved "
              f"({', '.join(f'c{c}:{orig[c]}->{wvc[r][c]}' for c in np.flatnonzero(wvc[r] != orig))}), "
              f"+{2e6 * exact:.4f} mHa^2 at fixed weights", flush=True)
    print(f"[wrep] {fixed} rows repaired (+{cost_tot:.4f} mHa^2 at fixed weights, {moved} words), {failed} left; "
          f"error_2^2 {2e6 * float(f):.4f} mHa^2 [{time.time() - t0:.0f} s]", flush=True)
    from chem_own2_lmstage7 import StageLM7
    lm = StageLM7(refc.T, refc.i, refc.j, edges, (R_, B_, C_), -0.5 - 0.5 / top, 0.5 - 0.5 / top, cache_dir="jax_cache")
    if "allowed" in z:      # t10 seat 6 (the lead's 16:20Z): the table's zero mask (True = a free weight)
        lm.wmask = np.asarray(z["allowed"], bool)
    frozen_all = np.ones((nrows, ncol), bool)
    thp = np.asarray(th) / np.pi
    _, w2, info = lm.refit(thp, np.asarray(Wf), frozen_all, np.zeros_like(frozen_all), steps=a.wsteps, seconds=1e9,
                           stop=0.001, lam_start=1e-3, stop_lam=0.1)
    Wf = jnp.asarray(w2)
    A, Rr, f = state(Pr, Wf)
    e2 = 2e6 * float(f)
    bf = int(bad_of(wvc, np.arange(nrows)).sum())
    print(f"[wrep] weights re-solved: {info['e2_in']:.4f} -> {e2:.4f} mHa^2", flush=True)
    assert wvc.min() >= -HALF and wvc.max() <= HALF - 1
    U = np.asarray(decode(th)).reshape(R_, B_, N)
    z["words"] = np.mod(wvc, 2 * top)
    z["W"] = np.asarray(Wf)
    z["U"] = U
    if lm.wmask is not None:    # t10 seat 6: the mask holds exactly
        assert np.all(np.asarray(z["W"])[~lm.wmask] == 0.0), "a masked weight moved"
    np.savez(a.out, **z)
    print(f"[wrep] END: error_2^2 {e2_0:.4f} -> {e2:.4f} mHa^2 ({e2 - e2_0:+.4f}); relations broken on the stored words "
          f"{bf}; rows repaired {fixed}, left {failed}; written {a.out} [{time.time() - t0:.0f} s]", flush=True)
    if bf:
        sys.exit(3)


if __name__ == "__main__":
    main()
