#!/usr/bin/env python
"""chem_t10_wpol (T10 seat 5, Nov 19): THE WIRE-CAPABLE ROW POLISH of a relation table (after its hold / polish).

WHY A NEW MOVE CLASS: tools/chem_t10_ppol.py keeps bit 0 of every word (only even moves; its row-DP toggles bit 1 in the
kernel of the bit-1 relations). In the P57x tables of card 30 bit 0 is constrained ONLY by the 45 wires (the --codes
flag is inactive under --no-relations [P, seat 5 01:08Z]), so odd moves are legal when every wire chain they cross keeps
its parity. Nothing has priced them since the rounding chose each row's bit-0 pattern at freeze time.

ONE PASS: (1) the exact change of error_2^2 at fixed weights of every single-word move d in +-1..+-span on every column
of every row (R8's one-row screen, batched as ppol); (2) per row and column the best value of each residue (2 or 4
residues; the present residue keeps 0 or its best +-4); (3) THE PLACER'S OWN PROGRAMME (tools/chem_t10_relge.RelPlacerGE
._dp, eps from the present words, the run's fitted targets) over all rows at once: the cheapest pattern that keeps all
57 relations on y = eps x z; the columns no relation reads add their best single move; (4) every row whose predicted
change is below -tol: the pattern is applied, checked on the row's STORED words (a moved pole changes eps) and by the
exact row screen; kept only if both pass (rank-2 update of the residual); (5) the weights re-solved (StageLM7).

    python chem_t10_wpol.py SRC.npz REL_USED_G.txt OUT.npz [--passes 4] [--minutes 50] [--span 4] [--tol 1e-4]
SRC / OUT: the format of held.npz / ppol's output (words mod 2^bits, W, U, tree_edges, ...)."""
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
    ap.add_argument("--passes", type=int, default=4)
    ap.add_argument("--minutes", type=float, default=50.0)
    ap.add_argument("--span", type=int, default=4)
    ap.add_argument("--tol", type=float, default=1e-4, help="mHa^2: a row is changed only if the exact screen gains this")
    ap.add_argument("--wsteps", type=int, default=12)
    ap.add_argument("--chunk", type=int, default=16384)
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
    wv = np.mod(words + top, 2 * top) - top
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

    b0 = int(bad_of(wv, np.arange(nrows)).sum())
    print(f"[wpol] {a.src}: {nrows} rows, {len(RELS)} relations a row; broken on the stored words at the start {b0}",
          flush=True)
    assert b0 == 0, "the input breaks a relation on its stored words (repair it first: tools/chem_t10_wrep.py)"
    fp = RelPlacerGE(R_, B_, RELS, edges, n_one=N, ncols=ncol, mode="exact", brot=bits)
    assert np.array_equal(fp.paddr[fp.place], addr)
    relcols = np.flatnonzero(fp.mod > 1)
    freecols = np.flatnonzero(fp.mod <= 1)
    allr = np.arange(nrows)

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

    def arows(A_, Wf_):
        Ag = A_.reshape(R_, C_, -1)
        return jnp.einsum("rcp,rbc->rbp", Ag, Wf_).reshape(nrows, -1)

    from chem_own2_lmstage7 import StageLM7
    lm = StageLM7(refc.T, refc.i, refc.j, edges, (R_, B_, C_), -0.5 - 0.5 / top, 0.5 - 0.5 / top, cache_dir="jax_cache")
    if "allowed" in z:      # t10 seat 6 (the lead's 16:20Z): the table's zero mask (True = a free weight)
        lm.wmask = np.asarray(z["allowed"], bool)
    frozen_all = np.ones((nrows, ncol), bool)

    def resolve(th_, Wf_):
        _, w2, info = lm.refit(np.asarray(th_) / np.pi, np.asarray(Wf_), frozen_all, np.zeros_like(frozen_all),
                               steps=a.wsteps, seconds=1e9, stop=0.001, lam_start=1e-3, stop_lam=0.1)
        return jnp.asarray(w2), info

    th = jnp.asarray(np.pi * wv / top)
    Pr = pvec(th)
    Wf = jnp.asarray(W)
    A, Rr, f = state(Pr, Wf)
    e2_0 = 2e6 * float(f)
    Wf, info = resolve(th, Wf)
    A, Rr, f = state(Pr, Wf)
    print(f"[wpol] start {e2_0:.4f} mHa^2; weights re-solved {2e6 * float(f):.4f} [{time.time() - t0:.0f} s]", flush=True)
    e2_w = 2e6 * float(f)
    D = np.array([d for d in range(-a.span, a.span + 1) if d != 0], np.int64)
    nd = len(D)
    wvc = wv.copy()
    total_rows = 0
    for p in range(a.passes):
        if time.time() - t0 > 60.0 * a.minutes:
            break
        tp = time.time()
        ar = arows(A, Wf)
        wn2 = jnp.sum(Wf.reshape(nrows, -1) ** 2, axis=1)
        thn = np.asarray(th)
        # (1) exact single-word deltas of every (row, column, d)
        R1, C1, K1 = np.meshgrid(allr, np.arange(ncol), np.arange(nd), indexing="ij")
        R1, C1, K1 = R1.ravel(), C1.ravel(), K1.ravel()
        nv = wvc[R1, C1] + D[K1]
        dd = np.full(len(R1), np.inf)
        idx = np.flatnonzero((nv >= -HALF) & (nv <= HALF - 1))
        for s0 in range(0, len(idx), a.chunk):
            sl = idx[s0:s0 + a.chunk]
            thr = thn[R1[sl]].copy()
            thr[np.arange(len(sl)), C1[sl]] += np.pi * D[K1[sl]] / top
            dd[sl] = np.asarray(screen(jnp.asarray(thr), Pr[R1[sl]], ar[R1[sl]], wn2[R1[sl]], Rr))
        dd = dd.reshape(nrows, ncol, nd)
        t_scr = time.time() - tp
        # (2) the best value of each residue, per (row, relation column)
        wc, bestd = {}, {}
        for c in relcols:
            mc = int(fp.mod[c])
            cst = np.full((nrows, 4), BIG)
            bd = np.zeros((nrows, 4), np.int64)
            q0 = np.mod(wvc[:, c], mc)
            cst[allr, q0] = 0.0
            for k, d in enumerate(D):
                q = np.mod(wvc[:, c] + d, mc)
                v = dd[:, c, k]
                better = v < cst[allr, q]
                cst[allr[better], q[better]] = v[better]
                bd[allr[better], q[better]] = d
            wc[c], bestd[c] = cst, bd
        # (3) the programme over all rows at their places, eps from the present words
        fp.update_eps(wvc.astype(np.float64), np.ones(wvc.shape, bool), wvc)
        places = fp.place[allr]
        cost, arg = fp._dp(wc, allr, places, diag=True)
        codes = fp._backtrack(arg, places)
        L = fp.lowp[places]
        newv = wvc.copy()
        for b in fp.batches:
            x = codes[b]
            for j in range(5):
                col = 5 * b + j
                mc = int(fp.mod[col])
                if mc <= 1:
                    continue
                c0, c1 = fp.need[b][j]
                b0_ = c0[x] if c0 is not None else None
                b1_ = c1[x] if c1 is not None else None
                if mc == 2:
                    qz = b0_
                else:
                    cst = wc[col]
                    neg = fp.E[allr, L, col] < 0
                    cy = np.where(neg[:, None], cst[:, PERM], cst)
                    best_q = np.zeros(nrows, np.int64)
                    best_v = np.full(nrows, np.inf)
                    for q in range(4):
                        ok = np.ones(nrows, bool)
                        if b0_ is not None:
                            ok &= (q & 1) == b0_
                        if b1_ is not None:
                            ok &= ((q >> 1) & 1) == b1_
                        v = np.where(ok, cy[:, q], np.inf)
                        sel = v < best_v
                        best_q[sel], best_v[sel] = q, v[sel]
                    qz = np.where(neg, PERM[best_q], best_q)
                newv[:, col] = wvc[:, col] + bestd[col][allr, qz]
        pred = np.where(cost < 1e17, cost, np.inf)
        for c in freecols:                      # the columns no relation reads: their best single move, if it gains
            kbest = np.argmin(dd[:, c, :], axis=1)
            vbest = dd[allr, c, kbest]
            take = vbest < 0
            newv[take, c] = wvc[take, c] + D[kbest[take]]
            pred = pred + np.where(take, vbest, 0.0)
        cand = np.flatnonzero((pred < -a.tol / 2e6) & np.any(newv != wvc, axis=1))
        cand = cand[np.argsort(pred[cand])]
        took, undone_rel, undone_ex, gain = 0, 0, 0, 0.0
        for r in cand:
            old = wvc[r].copy()
            wvc[r] = newv[r]
            if per_row(stored_words(np.mod(wvc[[r]], 2 * top), edges, addr[[r]], bits), addr[[r]], RELS).sum():
                wvc[r] = old
                undone_rel += 1
                continue
            thr = np.pi * wvc[r] / top
            g_, b_ = divmod(int(r), B_)
            w_r = Wf[g_, b_]
            a_r = A[g_ * C_:(g_ + 1) * C_, :].T @ w_r
            exact = float(screen(jnp.asarray(thr[None, :]), Pr[r][None, :], a_r[None, :], jnp.sum(w_r * w_r)[None], Rr)[0])
            if 2e6 * exact > -a.tol:
                wvc[r] = old
                undone_ex += 1
                continue
            Pnew = pvec(jnp.asarray(thr[None, :]))[0]
            dP = Pnew - Pr[r]
            A = A.at[g_ * C_:(g_ + 1) * C_, :].add(jnp.outer(w_r, dP))
            Pr = Pr.at[r].set(Pnew)
            th = th.at[r].set(jnp.asarray(thr))
            Rr = rank2(Rr, a_r, dP, jnp.sum(w_r * w_r))
            f = f + exact
            took += 1
            gain -= 2e6 * exact
        total_rows += took
        nodd = int(np.sum(np.mod(wvc - wv, 2) != 0))
        Wf, info = resolve(th, Wf)
        A, Rr, f = state(Pr, Wf)
        print(f"[wpol] pass {p + 1}: screens {len(idx)} in {t_scr:.0f} s; rows with a predicted gain {len(cand)}; taken {took} "
              f"(-{gain:.4f} mHa^2 exact at fixed weights), undone by the stored-word check {undone_rel}, by the exact "
              f"screen {undone_ex}; words with bit 0 changed so far {nodd}; weights re-solved: {2e6 * float(f):.4f} mHa^2 "
              f"({2e6 * float(f) - e2_w:+.4f}) [{time.time() - t0:.0f} s]", flush=True)
        if took == 0:
            break
    e2 = 2e6 * float(f)
    bf = int(bad_of(wvc, allr).sum())
    assert bf == 0, f"{bf} relations broken on the stored words"
    assert wvc.min() >= -HALF and wvc.max() <= HALF - 1
    U = np.asarray(decode(th)).reshape(R_, B_, N)
    z["words"] = np.mod(wvc, 2 * top)
    z["W"] = np.asarray(Wf)
    z["U"] = U
    if lm.wmask is not None:    # t10 seat 6: the mask holds exactly
        assert np.all(np.asarray(z["W"])[~lm.wmask] == 0.0), "a masked weight moved"
    np.savez(a.out, **z)
    print(f"[wpol] END: error_2^2 {e2_0:.4f} -> {e2:.4f} mHa^2 ({e2 - e2_0:+.4f}; {e2 - e2_w:+.4f} after the first weights "
          f"re-solve); rows changed {total_rows}, words changed {int(np.sum(wvc != wv))} (bit 0 in "
          f"{int(np.sum(np.mod(wvc - wv, 2) != 0))}); relations broken on the stored words 0; written {a.out} "
          f"[{time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
