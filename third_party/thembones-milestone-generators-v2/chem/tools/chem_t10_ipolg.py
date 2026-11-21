#!/usr/bin/env python
"""T10 copy (tools/chem_t10_mkipolg.py): wire relations checked against the input's own bit-0 XORs.
chem_own2_ipolish (R8 holder 7, Nov 17 04:50Z): THE INTEGER POLISH AFTER THE END (LEAD_NOTES 03:49Z, R8's second
item). On a finished 7-bit wired table (rounded.npz of the staged tool), greedy descent over discrete word moves that
keep every wire relation, accepted on the EXACT objective:
  * single-word moves of one step (+-1) on columns that carry no wire, +-2 on wire columns (bit 0 kept);
  * pair moves (+-1, +-1) on the two words of every wire (their bit-0 XOR kept);
  * every word stays in the signed range [-HALF, HALF - 1].
Exact screening: a move changes ONE row's vector u -> u', so with d = P(u') - P(u) (P = sqrt(m) u_i u_j), w the row's
weights, a = A_g^T w (A_g = the C rows of the row's group) and R = A^T A - T (all scaled as in chem_own2_lmstage):
  delta f = 2 a.Rd + |w|^2 d.Rd + |a|^2 |d|^2 + (a.d)^2 + 2 |w|^2 (a.d) |d|^2 + |w|^4 |d|^4 / 2      [P: expansion]
at FIXED weights; a move with delta f < 0 lowers the error, re-solving the weights afterwards only lowers it more.
Every --resolve accepted moves the weights are re-solved (weights-only LM of tools/chem_own2_lmstage_w.py, hard sign).
  python chem_own2_ipolish.py ROUNDED.npz WIRES.txt OUT.npz [--minutes 60] [--resolve 25] [--tol 1e-4]
OUT.npz has rounded.npz's format (words mod 2^bits, W, U, order, ...), the relations re-checked."""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("R8_LM_WONLY", "1")
from chem_own2_lmround import Target  # noqa: E402
import chem_own2_lmstage_w as LW  # noqa: E402
from r50_chem_lbl_chain import parse as wparse, lin as wlin  # noqa: E402

N = 76


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rounded")
    ap.add_argument("wires")
    ap.add_argument("out")
    ap.add_argument("--minutes", type=float, default=60.0)
    ap.add_argument("--resolve", type=int, default=100, help="weights re-solved every this many accepted moves")
    ap.add_argument("--wsteps", type=int, default=12, help="weights-only LM steps per re-solve")
    ap.add_argument("--tol", type=float, default=1e-4, help="accept a move only if it lowers error_2^2 by this (mHa^2)")
    ap.add_argument("--burst", type=int, default=256, help="moves accepted from one full screen (each re-checked exactly)")
    ap.add_argument("--chunk", type=int, default=16384)
    ap.add_argument("--weights-only", type=int, default=0, help="no moves: this many rounds of 20 weights-only LM steps")
    a = ap.parse_args()
    import jax
    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp

    t0 = time.time()
    z = dict(np.load(a.rounded))
    bits = int(z["bits"])
    top = 1 << (bits - 1)
    HALF = top // 2
    words = np.asarray(z["words"], np.int64)
    wv = np.mod(words + top, 2 * top) - top                     # signed words (no carriers in these runs)
    assert wv.min() >= -HALF and wv.max() <= HALF - 1, "a word outside the signed range"
    W = np.asarray(z["W"], np.float64)
    R_, B_, C_ = W.shape
    edges = np.asarray(z["tree_edges"], np.int64)
    WIRES = wparse(open(a.wires).read().strip())
    addr = N + np.arange(R_ * B_)
    wcol = np.zeros(N - 1, bool)
    pairs = []
    for (k_, j_, _, m_) in WIRES:
        c1, c2 = 5 * k_ - 5 + j_, 5 * k_ + j_
        wcol[c1] = wcol[c2] = True
        pairs.append((c1, c2, m_))

    REF_ = {(c1, c2): ((wv[:, c2] ^ wv[:, c1]) & 1) for c1, c2, m_ in pairs}    # t10: the input's own bit-0 XORs (G table)
    def broken(wv_):
        return sum(int(np.sum(((wv_[:, c2] ^ wv_[:, c1]) & 1) != REF_[(c1, c2)])) for c1, c2, m_ in pairs)
    assert broken(wv) == 0, "the input breaks a wire relation"

    refc = Target()
    I, J = np.asarray(refc.i), np.asarray(refc.j)
    m = np.where(I == J, 1.0, 2.0)
    sm = np.sqrt(m)
    Tt = np.outer(sm, sm) * refc.T
    lv = LW._levels(edges)

    def decode(th):                     # radians (n, 75) -> (n, 76)
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
    def state(Pr, Wf):                  # Pr (R*B, P), Wf (R, B, C) -> A (R*C, P), Rres (P, P), f
        A = jnp.einsum("rbc,rbp->rcp", Wf, Pr.reshape(R_, B_, -1)).reshape(R_ * C_, -1)
        Rr = A.T @ A - Ttj
        return A, Rr, 0.5 * jnp.sum(Rr * Rr)

    @jax.jit
    def screen(th_rows, Prows, arow, wnorm2, Rr):
        """th_rows (n, 75) the moved rows (radians), Prows (n, P) their present P, arow (n, P) = A_g^T w, wnorm2 (n,)."""
        d = pvec(th_rows) - Prows
        Rd = d @ Rr
        ad = jnp.sum(arow * d, axis=1)
        dd = jnp.sum(d * d, axis=1)
        return (2.0 * jnp.sum(arow * Rd, axis=1) + wnorm2 * jnp.sum(d * Rd, axis=1) + jnp.sum(arow * arow, axis=1) * dd
                + ad * ad + 2.0 * wnorm2 * ad * dd + 0.5 * wnorm2 * wnorm2 * dd * dd)

    @jax.jit
    def rank2(Rr, a_, d_, c_):
        return Rr + jnp.outer(a_, d_) + jnp.outer(d_, a_) + c_ * jnp.outer(d_, d_)

    # candidate list: (row, col1, delta1, col2, delta2); col2 = -1 for a single move
    cand = []
    rows_all = np.arange(R_ * B_)
    for c in range(N - 1):
        for dlt in ((1, -1) if not wcol[c] else (2, -2)):
            cand.append(np.stack([rows_all, np.full_like(rows_all, c), np.full_like(rows_all, dlt),
                                  np.full_like(rows_all, -1), np.zeros_like(rows_all)], 1))
    from collections import Counter
    use_ = Counter([c for c1_, c2_, _ in pairs for c in (c1_, c2_)])
    for c1, c2, _ in pairs:
        if use_[c1] > 1 or use_[c2] > 1:       # a chained wire (a column in two wires): no pair move here
            continue
        for d1 in (1, -1):
            for d2 in (1, -1):
                cand.append(np.stack([rows_all, np.full_like(rows_all, c1), np.full_like(rows_all, d1),
                                      np.full_like(rows_all, c2), np.full_like(rows_all, d2)], 1))
    cand = np.concatenate(cand)
    print(f"[ipolish] {a.rounded}: {R_ * B_} rows, {len(WIRES)} wires on {int(wcol.sum())} columns, {len(cand)} moves",
          flush=True)

    th = jnp.asarray(np.pi * wv / top)          # radians: word * pi / 2^(bits-1)
    Pr = pvec(th)
    Wf = jnp.asarray(W)
    A, Rr, f = state(Pr, Wf)
    e2_0 = 2e6 * float(f)
    print(f"[ipolish] start error_2^2 {e2_0:.4f} mHa^2 (error_2 {np.sqrt(e2_0):.4f} mHa) [{time.time() - t0:.0f} s]", flush=True)

    # R8 holder 7 (06:18Z): the weight re-solves use the backtracking step (the old step stalls on a finished table)
    from chem_own2_lmstage7 import StageLM7
    lm = StageLM7(refc.T, refc.i, refc.j, edges, (R_, B_, C_), -0.5 - 0.5 / top, 0.5 - 0.5 / top, cache_dir="jax_cache")
    frozen_all = np.ones((R_ * B_, N - 1), bool)

    def resolve(th_, Wf_):
        thp = np.asarray(th_) / np.pi
        _, w2, info = lm.refit(thp, np.asarray(Wf_), frozen_all, np.zeros_like(frozen_all), steps=a.wsteps, seconds=1e9,
                               stop=0.001, lam_start=1e-3, stop_lam=0.1)
        return jnp.asarray(w2), info

    Wf, info = resolve(th, Wf)
    A, Rr, f = state(Pr, Wf)
    print(f"[ipolish] weights re-solved: {info['e2_in']:.4f} -> {2e6 * float(f):.4f} mHa^2 [{time.time() - t0:.0f} s]", flush=True)

    def arows(A_, Wf_):                 # a_r = A_g^T w_r for every row
        Ag = A_.reshape(R_, C_, -1)
        return jnp.einsum("rcp,rbc->rbp", Ag, Wf_).reshape(R_ * B_, -1)

    wvc = wv.copy()
    accepted, gain_acc, it = 0, 0.0, 0
    log_every = time.time()
    if a.weights_only > 0:      # R8 holder 7 (06:15Z): NO moves: the END's weights alone, a.weights_only rounds of 20
        # with the backtracking step (chem_own2_lmstage7): the old step stalls on the finished table [M 06:12Z: 0.0000
        # in 20 steps, where lmhold7 took -0.24]; every word is frozen, so the step moves the weights only
        from chem_own2_lmstage7 import StageLM7
        lm7 = StageLM7(refc.T, refc.i, refc.j, edges, (R_, B_, C_), -0.5 - 0.5 / top, 0.5 - 0.5 / top, cache_dir="jax_cache")
        for i_ in range(a.weights_only):
            _, w2, info = lm7.refit(np.asarray(th) / np.pi, np.asarray(Wf), frozen_all, np.zeros_like(frozen_all),
                                    steps=20, seconds=1e9, stop=-1.0)
            Wf = jnp.asarray(w2)
            A, Rr, f = state(Pr, Wf)
            print(f"[ipolish] weights only, round {i_ + 1}: {2e6 * float(f):.4f} mHa^2 [{time.time() - t0:.0f} s]", flush=True)
        a.minutes = 0.0
    while time.time() - t0 < 60.0 * a.minutes:
        it += 1
        ar = arows(A, Wf)
        wn2 = jnp.sum(Wf.reshape(R_ * B_, -1) ** 2, axis=1)
        # the moves that keep the signed range
        r_, c1_, d1_, c2_, d2_ = cand.T
        ok = (wvc[r_, c1_] + d1_ >= -HALF) & (wvc[r_, c1_] + d1_ <= HALF - 1)
        two = c2_ >= 0
        ok &= ~two | ((wvc[r_, np.maximum(c2_, 0)] + d2_ >= -HALF) & (wvc[r_, np.maximum(c2_, 0)] + d2_ <= HALF - 1))
        idx = np.flatnonzero(ok)
        df = np.empty(len(idx))
        thn = np.asarray(th)
        for s0 in range(0, len(idx), a.chunk):
            sl = idx[s0:s0 + a.chunk]
            rr, cc1, dd1, cc2, dd2 = cand[sl].T
            thr = thn[rr].copy()
            thr[np.arange(len(sl)), cc1] += np.pi * dd1 / top
            m2 = cc2 >= 0
            thr[np.flatnonzero(m2), cc2[m2]] += np.pi * dd2[m2] / top
            df[s0:s0 + len(sl)] = np.asarray(screen(jnp.asarray(thr), Pr[rr], ar[rr], wn2[rr], Rr))
        order_ = np.argsort(df)
        best = 2e6 * df[order_[0]]
        if best > -a.tol:
            Wn, info = resolve(th, Wf)
            An, Rn, fn = state(Pr, Wn)
            if 2e6 * (float(f) - float(fn)) > a.tol:
                Wf, A, Rr, f = Wn, An, Rn, fn
                print(f"[ipolish] no improving move; weights re-solved to {2e6 * float(f):.4f} mHa^2", flush=True)
                continue
            print(f"[ipolish] converged: no move lowers the error by {a.tol} mHa^2", flush=True)
            break
        # accept up to --burst moves from this screen, each re-checked exactly against the updated state
        took, rows_taken = 0, set()
        for o in order_[:4 * a.burst]:
            if took >= a.burst or 2e6 * df[o] > -a.tol:
                break
            rr, cc1, dd1, cc2, dd2 = cand[idx[o]]
            if rr in rows_taken:
                continue
            thr = np.asarray(th)[rr].copy()
            thr[cc1] += np.pi * dd1 / top
            if cc2 >= 0:
                thr[cc2] += np.pi * dd2 / top
            g_, b_ = divmod(int(rr), B_)
            w_r = Wf[g_, b_]
            a_r = A[g_ * C_:(g_ + 1) * C_, :].T @ w_r                  # a = A_g^T w of this row, present state
            exact = float(screen(jnp.asarray(thr[None, :]), Pr[rr][None, :], a_r[None, :], wn2[rr][None], Rr)[0])
            if 2e6 * exact > -a.tol:
                continue
            Pnew = pvec(jnp.asarray(thr[None, :]))[0]
            dP = Pnew - Pr[rr]
            A = A.at[g_ * C_:(g_ + 1) * C_, :].add(jnp.outer(w_r, dP))
            Pr = Pr.at[rr].set(Pnew)
            th = th.at[rr].set(jnp.asarray(thr))
            wvc[rr, cc1] += dd1
            if cc2 >= 0:
                wvc[rr, cc2] += dd2
            Rr = rank2(Rr, a_r, dP, jnp.sum(w_r * w_r))             # R += a d^T + d a^T + |w|^2 d d^T (exact)
            gain_acc -= 2e6 * exact
            f = f + exact
            accepted += 1
            took += 1
            rows_taken.add(int(rr))
            if accepted % a.resolve == 0:
                Wf, info = resolve(th, Wf)
                A, Rr, f = state(Pr, Wf)            # also clears any drift of the rank-2 updates
                wn2 = jnp.sum(Wf.reshape(R_ * B_, -1) ** 2, axis=1)
        if time.time() - log_every > 60 or it <= 3:
            log_every = time.time()
            print(f"[ipolish] screen {it}: best {best:+.4f}, {int(np.sum(df < -a.tol / 2e6))} improving of {len(idx)}; "
                  f"accepted {accepted}; error_2^2 {2e6 * float(f):.4f} mHa^2 ({2e6 * float(f) - e2_0:+.3f}) "
                  f"[{time.time() - t0:.0f} s]", flush=True)
    Wf, info = resolve(th, Wf)
    A, Rr, f = state(Pr, Wf)
    e2 = 2e6 * float(f)
    assert broken(wvc) == 0, "a wire relation broke"
    assert wvc.min() >= -HALF and wvc.max() <= HALF - 1
    U = np.asarray(decode(th)).reshape(R_, B_, N)
    z["words"] = np.mod(wvc + top * 0, 2 * top)
    z["W"] = np.asarray(Wf)
    z["U"] = U
    np.savez(a.out, **z)
    nch = int(np.sum(wvc != wv))
    print(f"[ipolish] END: error_2^2 {e2_0:.4f} -> {e2:.4f} mHa^2 ({e2 - e2_0:+.3f}); error_2 {np.sqrt(e2):.4f} mHa; "
          f"{accepted} moves accepted, {nch} words changed, relations broken 0; written {a.out} [{time.time() - t0:.0f} s]",
          flush=True)


if __name__ == "__main__":
    main()
