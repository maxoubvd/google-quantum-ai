#!/usr/bin/env python
"""chem_t10_ppol (T10 seat 4): THE INTEGER POLISH OF A RELATION (P) TABLE AFTER ITS HOLD, every relation kept on the
STORED words. A copy of R8's polish (tools/chem_t10_ipolg.py = chem_own2_ipolish with G checks: the exact one-row screen,
bursts, rank-2 updates, weights re-solved by StageLM7) with the move set of a table under bit-1 relations:
  * bit 0 of every word is kept (wires, the bit-0 codes, the carriers' parity): only even moves;
  * +-2 on one word of a column that no bit-1 relation reads;
  * BIT-1 KERNEL MOVES (seat 3's one-body polish, tools/chem_t10_obcar.py): +-2 on two words of one row whose bit-1
    syndromes are equal (and, with --triples, three words whose syndromes cancel): every bit-1 relation of the row keeps
    its parity on the signed words;
  * every accepted move is re-checked EXACTLY on the row's STORED words (R53's carrier form, tools/chem_t10_relge.py
    stored_words) against the run's fitted targets (rel_used_g.txt): a move that breaks one relation there is undone.
Rows of the held / rounded table are in address order (address = 76 + row), as the END check of the rounding reads them.
  python chem_t10_ppol.py HELD.npz REL_USED_G.txt OUT.npz [--minutes 25] [--resolve 100] [--triples]
OUT.npz has the input's format (words mod 2^bits, W, U, fp_G, ...)."""
import argparse
import itertools
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("R8_LM_WONLY", "1")
from chem_own2_lmround import Target  # noqa: E402
import chem_own2_lmstage_w as LW  # noqa: E402
from r65_chem_all15_relg import parse_g, check_g  # noqa: E402
from chem_t10_relge import stored_words  # noqa: E402

N = 76


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rounded")
    ap.add_argument("rels")
    ap.add_argument("out")
    ap.add_argument("--minutes", type=float, default=25.0)
    ap.add_argument("--resolve", type=int, default=100, help="weights re-solved every this many accepted moves")
    ap.add_argument("--wsteps", type=int, default=12, help="weights-only LM steps per re-solve")
    ap.add_argument("--tol", type=float, default=1e-4, help="accept a move only if it lowers error_2^2 by this (mHa^2)")
    ap.add_argument("--burst", type=int, default=256, help="moves accepted from one full screen (each re-checked exactly)")
    ap.add_argument("--chunk", type=int, default=16384)
    ap.add_argument("--triples", action="store_true", help="also the bit-1 kernel triples")
    ap.add_argument("--four", action="store_true", help="also +-4 on every word (bits 0 and 1 kept: no relation can break)")
    ap.add_argument("--dp", type=int, default=0, help="seat 4 fix round 3: this many ROW-DP rounds after the greedy phase: "
                    "per row the cheapest bit-1 toggle pattern in the kernel of the relations (Viterbi over the syndromes) "
                    "from exact single-word deltas (+-2 toggles, 0 / +-4 keeps), applied when the exact row screen improves")
    ap.add_argument("--dp-minutes", type=float, default=30.0)
    ap.add_argument("--repair", action="store_true", help="seat 4: rows whose relations are broken on the stored words at "
                    "the start (late eps changes) are repaired first: the cheapest bit-1 toggle set whose syndrome is the "
                    "row's violation pattern (Viterbi from exact single-word deltas), checked on the stored words")
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
    RELS = parse_g(open(a.rels).read())
    addr = N + np.arange(R_ * B_)
    ncol = N - 1

    def broken_rows(wv_, rows):         # relations broken on the STORED words of these rows (address order)
        rows = np.asarray(rows, np.int64)
        st = stored_words(np.mod(wv_[rows], 2 * top), edges, addr[rows], bits)
        return int(check_g(st, addr[rows], RELS))

    b0 = broken_rows(wv, np.arange(R_ * B_))
    print(f"[ppol] {a.rounded}: {R_ * B_} rows, {len(RELS)} relations; broken on the stored words at the start: {b0}",
          flush=True)
    assert b0 == 0 or a.repair, "the input breaks a relation on its stored words (use --repair)"

    def viol(wv_, r):                   # the row's violation pattern on the stored words: bit n = relation n broken
        st = stored_words(np.mod(wv_[[r]], 2 * top), edges, addr[[r]], bits)
        per = check_g(st, addr[[r]], RELS, per=True)
        return [i_ for i_, x in enumerate(per) if x]

    # the bit-1 syndromes of the columns (every relation that reads bit 1 of a column's word, as loader k or k-1)
    syn1 = [0] * ncol
    for i_, (k_, o_, _t) in enumerate(RELS):
        for (j_, b_) in o_:
            if b_ == 1:
                syn1[5 * k_ + j_] ^= 1 << i_
                syn1[5 * k_ - 5 + j_] ^= 1 << i_
            assert b_ <= 1, "a relation on bit 2 or higher: not covered by this polish"
    moves = []                          # (cols, deltas)
    for c in range(ncol):
        if syn1[c] == 0:
            moves += [((c,), (2,)), ((c,), (-2,))]
        if a.four:
            moves += [((c,), (4,)), ((c,), (-4,))]
    by1 = {}
    for c in range(ncol):
        if syn1[c]:
            by1.setdefault(syn1[c], []).append(c)
    npair = 0
    for grp in by1.values():
        for c1, c2 in itertools.combinations(grp, 2):
            for sg in itertools.product((2, -2), repeat=2):
                moves.append(((c1, c2), sg))
                npair += 1
    ntri = 0
    if a.triples:
        cols1 = [c for c in range(ncol) if syn1[c]]
        seen = set()
        for c1, c2 in itertools.combinations(cols1, 2):
            for c3 in by1.get(syn1[c1] ^ syn1[c2], []):
                t3 = tuple(sorted((c1, c2, c3)))
                if c3 in (c1, c2) or t3 in seen:
                    continue
                seen.add(t3)
                for sg in itertools.product((2, -2), repeat=3):
                    moves.append((t3, sg))
                    ntri += 1
    nm = len(moves)
    Dm = np.zeros((nm, ncol), np.int64)
    for i, (S, sg) in enumerate(moves):
        for c, d in zip(S, sg):
            Dm[i, c] += d
    print(f"[ppol] moves a row: {nm} (single {nm - npair - ntri}, kernel pairs {npair}, triples {ntri}); "
          f"columns without a bit-1 relation {sum(1 for s in syn1 if s == 0)}", flush=True)

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
    print(f"[ppol] start error_2^2 {e2_0:.4f} mHa^2 [{time.time() - t0:.0f} s]", flush=True)
    from chem_own2_lmstage7 import StageLM7
    lm = StageLM7(refc.T, refc.i, refc.j, edges, (R_, B_, C_), -0.5 - 0.5 / top, 0.5 - 0.5 / top, cache_dir="jax_cache")
    if "allowed" in z:      # t10 seat 6 (the lead's 16:20Z): the table's zero mask (True = a free weight)
        lm.wmask = np.asarray(z["allowed"], bool)
    frozen_all = np.ones((R_ * B_, ncol), bool)

    def resolve(th_, Wf_):
        thp = np.asarray(th_) / np.pi
        _, w2, info = lm.refit(thp, np.asarray(Wf_), frozen_all, np.zeros_like(frozen_all), steps=a.wsteps, seconds=1e9,
                               stop=0.001, lam_start=1e-3, stop_lam=0.1)
        return jnp.asarray(w2), info

    Wf, info = resolve(th, Wf)
    A, Rr, f = state(Pr, Wf)
    print(f"[ppol] weights re-solved: {info['e2_in']:.4f} -> {2e6 * float(f):.4f} mHa^2 [{time.time() - t0:.0f} s]",
          flush=True)

    def arows(A_, Wf_):
        Ag = A_.reshape(R_, C_, -1)
        return jnp.einsum("rcp,rbc->rbp", Ag, Wf_).reshape(R_ * B_, -1)

    wvc = wv.copy()
    accepted, undone, it = 0, 0, 0
    nrows = R_ * B_
    if a.repair and b0:     # ---- seat 4: the END repair on the signed words (stored-word relations as the target)
        from r65_chem_all15_relg import ttable_g
        rel1_ = [i_ for i_, (k_, o_, _t) in enumerate(RELS) if any(b_ == 1 for (_, b_) in o_)]
        bit_ = {i_: n_ for n_, i_ in enumerate(rel1_)}
        syc = np.zeros(ncol, np.int64)
        for c in range(ncol):
            for i_ in rel1_:
                if (syn1[c] >> i_) & 1:
                    syc[c] |= 1 << bit_[i_]
        NS_ = 1 << len(rel1_)
        st_all = stored_words(np.mod(wvc, 2 * top), edges, addr, bits)
        bad = np.zeros((nrows, len(RELS)), bool)
        for i_, (k_, o_, t_) in enumerate(RELS):
            acc = np.zeros(nrows, np.int64)
            for (j_, b_) in o_:
                acc ^= ((st_all[:, 5 * k_ + j_] ^ st_all[:, 5 * k_ - 5 + j_]) >> b_) & 1
            bad[:, i_] = acc != ttable_g(t_, addr)
        ar = arows(A, Wf)
        wn2 = jnp.sum(Wf.reshape(nrows, -1) ** 2, axis=1)
        D1r = np.array([2, -2, 4, -4], np.int64)
        fixed, failed, cost_r = 0, 0, 0.0
        for r in np.flatnonzero(bad.any(axis=1)):
            brk = np.flatnonzero(bad[r])
            if any(i_ not in bit_ for i_ in brk):
                failed += 1
                print(f"[ppol] repair: row {r} breaks a relation that is not bit-1 ({list(brk)}): left", flush=True)
                continue
            T_ = 0
            for i_ in brk:
                T_ |= 1 << bit_[i_]
            thn = np.asarray(th)
            nv = wvc[r][None, :] + D1r[:, None]                       # (4, cols)
            dd = np.full((4, ncol), np.inf)
            for q in range(4):
                okq = (nv[q] >= -HALF) & (nv[q] <= HALF - 1)
                cs = np.flatnonzero(okq)
                thr = np.repeat(thn[r][None, :], len(cs), axis=0)
                thr[np.arange(len(cs)), cs] += np.pi * D1r[q] / top
                dd[q, cs] = np.asarray(screen(jnp.asarray(thr), jnp.repeat(Pr[r][None, :], len(cs), axis=0),
                                              jnp.repeat(ar[r][None, :], len(cs), axis=0), jnp.repeat(wn2[r][None], len(cs)), Rr))
            togc, togd = np.minimum(dd[0], dd[1]), np.where(dd[0] <= dd[1], 2, -2)
            kp = np.stack([np.zeros(ncol), dd[2], dd[3]])
            keepc, keepd = kp.min(axis=0), np.array([0, 4, -4])[kp.argmin(axis=0)]
            cost = np.full(NS_, np.inf)
            cost[0] = 0.0
            back = np.zeros((ncol, NS_), bool)
            ars = np.arange(NS_)
            for c in range(ncol):
                stay = cost + keepc[c]
                flip = cost[ars ^ syc[c]] + togc[c]
                back[c] = flip < stay
                cost = np.where(back[c], flip, stay)
            if not np.isfinite(cost[T_]):
                failed += 1
                continue
            s, dv = T_, np.zeros(ncol, np.int64)
            for c in range(ncol - 1, -1, -1):
                if back[c, s]:
                    dv[c] = togd[c]
                    s ^= syc[c]
                else:
                    dv[c] = keepd[c]
            old = wvc[r].copy()
            wvc[r] = old + dv
            if viol(wvc, r):
                wvc[r] = old
                failed += 1
                continue
            thr = np.pi * wvc[r] / top
            g_, b_ = divmod(int(r), B_)
            w_r = Wf[g_, b_]
            a_r = A[g_ * C_:(g_ + 1) * C_, :].T @ w_r
            exact = float(screen(jnp.asarray(thr[None, :]), Pr[r][None, :], a_r[None, :], wn2[r][None], Rr)[0])
            Pnew = pvec(jnp.asarray(thr[None, :]))[0]
            dP = Pnew - Pr[r]
            A = A.at[g_ * C_:(g_ + 1) * C_, :].add(jnp.outer(w_r, dP))
            Pr = Pr.at[r].set(Pnew)
            th = th.at[r].set(jnp.asarray(thr))
            Rr = rank2(Rr, a_r, dP, jnp.sum(w_r * w_r))
            f = f + exact
            fixed += 1
            cost_r += 2e6 * exact
        Wf, info = resolve(th, Wf)
        A, Rr, f = state(Pr, Wf)
        print(f"[ppol] repair: {fixed} rows repaired (+{cost_r:.4f} mHa^2 at fixed weights), {failed} left; broken on the "
              f"stored words now {broken_rows(wvc, np.arange(nrows))}; error_2^2 {2e6 * float(f):.4f} mHa^2", flush=True)
    log_every = time.time()
    while time.time() - t0 < 60.0 * a.minutes:
        it += 1
        ar = arows(A, Wf)
        wn2 = jnp.sum(Wf.reshape(nrows, -1) ** 2, axis=1)
        # every (row, move) inside the signed range
        new = wvc[:, None, :] + Dm[None, :, :]                       # (rows, moves, cols)
        okm = np.all((new >= -HALF) & (new <= HALF - 1), axis=2)
        rr_all, mm_all = np.nonzero(okm)
        df = np.empty(len(rr_all))
        thn = np.asarray(th)
        for s0 in range(0, len(rr_all), a.chunk):
            rr, mm = rr_all[s0:s0 + a.chunk], mm_all[s0:s0 + a.chunk]
            thr = thn[rr] + np.pi * Dm[mm] / top
            df[s0:s0 + len(rr)] = np.asarray(screen(jnp.asarray(thr), Pr[rr], ar[rr], wn2[rr], Rr))
        order_ = np.argsort(df)
        best = 2e6 * df[order_[0]] if len(df) else 0.0
        if best > -a.tol:
            Wn, info = resolve(th, Wf)
            An, Rn, fn = state(Pr, Wn)
            if 2e6 * (float(f) - float(fn)) > a.tol:
                Wf, A, Rr, f = Wn, An, Rn, fn
                print(f"[ppol] no improving move; weights re-solved to {2e6 * float(f):.4f} mHa^2", flush=True)
                continue
            print(f"[ppol] converged: no move lowers the error by {a.tol} mHa^2", flush=True)
            break
        took, rows_taken = 0, set()
        for o in order_[:4 * a.burst]:
            if took >= a.burst or 2e6 * df[o] > -a.tol:
                break
            rr, mm = int(rr_all[o]), int(mm_all[o])
            if rr in rows_taken:
                continue
            old = wvc[rr].copy()
            wvc[rr] = old + Dm[mm]
            if broken_rows(wvc, [rr]) != 0:         # the stored words of this row must keep every relation
                wvc[rr] = old
                undone += 1
                continue
            thr = np.pi * wvc[rr] / top
            g_, b_ = divmod(rr, B_)
            w_r = Wf[g_, b_]
            a_r = A[g_ * C_:(g_ + 1) * C_, :].T @ w_r
            exact = float(screen(jnp.asarray(thr[None, :]), Pr[rr][None, :], a_r[None, :], wn2[rr][None], Rr)[0])
            if 2e6 * exact > -a.tol:
                wvc[rr] = old
                continue
            Pnew = pvec(jnp.asarray(thr[None, :]))[0]
            dP = Pnew - Pr[rr]
            A = A.at[g_ * C_:(g_ + 1) * C_, :].add(jnp.outer(w_r, dP))
            Pr = Pr.at[rr].set(Pnew)
            th = th.at[rr].set(jnp.asarray(thr))
            Rr = rank2(Rr, a_r, dP, jnp.sum(w_r * w_r))
            f = f + exact
            accepted += 1
            took += 1
            rows_taken.add(rr)
            if accepted % a.resolve == 0:
                Wf, info = resolve(th, Wf)
                A, Rr, f = state(Pr, Wf)
                wn2 = jnp.sum(Wf.reshape(nrows, -1) ** 2, axis=1)
        if time.time() - log_every > 60 or it <= 3:
            log_every = time.time()
            print(f"[ppol] screen {it}: best {best:+.4f}, {int(np.sum(df < -a.tol / 2e6))} improving of {len(df)}; "
                  f"accepted {accepted}, undone by the stored-word check {undone}; error_2^2 {2e6 * float(f):.4f} mHa^2 "
                  f"({2e6 * float(f) - e2_0:+.3f}) [{time.time() - t0:.0f} s]", flush=True)
    # ---- seat 4 fix round 3: ROW-DP rounds (every relation-keeping toggle pattern of a row, not only pairs / triples)
    rel1 = [i_ for i_, (k_, o_, _t) in enumerate(RELS) if any(b_ == 1 for (_, b_) in o_)]
    bitof = {i_: n_ for n_, i_ in enumerate(rel1)}
    sync = np.zeros(ncol, np.int64)                  # compact syndrome of each column over the bit-1 relations
    for c in range(ncol):
        for i_ in rel1:
            if (syn1[c] >> i_) & 1:
                sync[c] |= 1 << bitof[i_]
    NS = 1 << len(rel1)
    D1 = np.array([2, -2, 4, -4], np.int64)          # single-word deltas: 2 / -2 toggle bit 1, 4 / -4 keep bits 0-1
    t_dp = time.time()
    for rnd in range(a.dp):
        if time.time() - t_dp > 60.0 * a.dp_minutes:
            break
        ar = arows(A, Wf)
        wn2 = jnp.sum(Wf.reshape(nrows, -1) ** 2, axis=1)
        thn = np.asarray(th)
        # exact single-word deltas of every (row, column, d)
        R1, C1, K1 = np.meshgrid(np.arange(nrows), np.arange(ncol), np.arange(4), indexing="ij")
        R1, C1, K1 = R1.ravel(), C1.ravel(), K1.ravel()
        nv = wvc[R1, C1] + D1[K1]
        okv = (nv >= -HALF) & (nv <= HALF - 1)
        dd = np.full(len(R1), np.inf)
        idx = np.flatnonzero(okv)
        for s0 in range(0, len(idx), a.chunk):
            sl = idx[s0:s0 + a.chunk]
            thr = thn[R1[sl]].copy()
            thr[np.arange(len(sl)), C1[sl]] += np.pi * D1[K1[sl]] / top
            dd[sl] = np.asarray(screen(jnp.asarray(thr), Pr[R1[sl]], ar[R1[sl]], wn2[R1[sl]], Rr))
        dd = dd.reshape(nrows, ncol, 4)
        tog_c = np.minimum(dd[:, :, 0], dd[:, :, 1])
        tog_d = np.where(dd[:, :, 0] <= dd[:, :, 1], 2, -2)
        kp = np.stack([np.zeros((nrows, ncol)), dd[:, :, 2], dd[:, :, 3]], axis=2)
        keep_c = kp.min(axis=2)
        keep_d = np.array([0, 4, -4])[kp.argmin(axis=2)]
        took_r, tried, gain_r = 0, 0, 0.0
        for r in range(nrows):
            # Viterbi over the columns: state = accumulated syndrome; choice 1 = toggle bit 1
            cost = np.full(NS, np.inf)
            cost[0] = 0.0
            back = np.zeros((ncol, NS), bool)
            ar_s = np.arange(NS)
            for c in range(ncol):
                stay = cost + keep_c[r, c]
                flip = cost[ar_s ^ sync[c]] + tog_c[r, c]
                back[c] = flip < stay
                cost = np.where(back[c], flip, stay)
            if not np.isfinite(cost[0]) or cost[0] > -a.tol / 2e6:
                continue
            s, dv = 0, np.zeros(ncol, np.int64)
            for c in range(ncol - 1, -1, -1):
                if back[c, s]:
                    dv[c] = tog_d[r, c]
                    s ^= sync[c]
                else:
                    dv[c] = keep_d[r, c]
            if not np.any(dv) or np.count_nonzero(dv) < 2:
                continue                                  # single-word moves are the greedy phase's
            tried += 1
            old = wvc[r].copy()
            wvc[r] = old + dv
            if broken_rows(wvc, [r]) != 0:
                wvc[r] = old
                undone += 1
                continue
            thr = np.pi * wvc[r] / top
            g_, b_ = divmod(r, B_)
            w_r = Wf[g_, b_]
            a_r = A[g_ * C_:(g_ + 1) * C_, :].T @ w_r
            exact = float(screen(jnp.asarray(thr[None, :]), Pr[r][None, :], a_r[None, :], wn2[r][None], Rr)[0])
            if 2e6 * exact > -a.tol:
                wvc[r] = old
                continue
            Pnew = pvec(jnp.asarray(thr[None, :]))[0]
            dP = Pnew - Pr[r]
            A = A.at[g_ * C_:(g_ + 1) * C_, :].add(jnp.outer(w_r, dP))
            Pr = Pr.at[r].set(Pnew)
            th = th.at[r].set(jnp.asarray(thr))
            Rr = rank2(Rr, a_r, dP, jnp.sum(w_r * w_r))
            f = f + exact
            accepted += 1
            took_r += 1
            gain_r -= 2e6 * exact
        Wf, info = resolve(th, Wf)
        A, Rr, f = state(Pr, Wf)
        print(f"[ppol] row-DP round {rnd + 1}: {tried} rows with a predicted gain, {took_r} taken (-{gain_r:.4f} mHa^2 "
              f"exact); error_2^2 {2e6 * float(f):.4f} mHa^2 ({2e6 * float(f) - e2_0:+.3f}) [{time.time() - t0:.0f} s]",
              flush=True)
        if took_r == 0:
            break
    Wf, info = resolve(th, Wf)
    A, Rr, f = state(Pr, Wf)
    e2 = 2e6 * float(f)
    bf = broken_rows(wvc, np.arange(nrows))
    assert bf == 0, f"{bf} relations broken on the stored words"
    assert wvc.min() >= -HALF and wvc.max() <= HALF - 1
    assert np.all(np.mod(wvc, 2) == np.mod(wv, 2)), "a bit 0 changed"
    U = np.asarray(decode(th)).reshape(R_, B_, N)
    z["words"] = np.mod(wvc, 2 * top)
    z["W"] = np.asarray(Wf)
    z["U"] = U
    if lm.wmask is not None:    # t10 seat 6: the mask holds exactly
        assert np.all(np.asarray(z["W"])[~lm.wmask] == 0.0), "a masked weight moved"
    np.savez(a.out, **z)
    nch = int(np.sum(wvc != wv))
    print(f"[ppol] END: error_2^2 {e2_0:.4f} -> {e2:.4f} mHa^2 ({e2 - e2_0:+.3f}); error_2 {np.sqrt(e2):.4f} mHa; "
          f"{accepted} moves accepted ({undone} undone by the stored-word check), {nch} words changed, relations broken "
          f"on the stored words 0, bit 0 kept; written {a.out} [{time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
