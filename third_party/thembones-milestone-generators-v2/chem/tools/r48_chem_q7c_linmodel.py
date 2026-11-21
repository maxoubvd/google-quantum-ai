"""r48_chem_q7c_linmodel: idealised reading for JOINT rounding decisions in a linear model of one group (chem_own2).

Model. At a continuous optimum the error changes by d' H d for small changes d of the parameters, H = J'J the
Gauss-Newton matrix of the residual r = sqrt(M) * (A'A - T). For one group r every derivative of A'A is a symmetric
rank-two matrix x y' + y x' (s-scaled pair space, s = sqrt(m)): a weight (b, c) has x = a_rc (row of A), y = p_b
(pair vector of row b); an angle (b, t) has x = g_b = sum_c w_rbc a_rc, y = q_bt = d p_b / d theta_t. Two of them
have the inner product 2 (x1.x2)(y1.y2) + 2 (x1.y2)(y1.x2), so H comes from small Gram matrices. The other groups
are held fixed (no cross-group compensation: an upper reading of every damage below).

What is read (7-bit steps, uncoded, the group's weights always free):
  S      = Schur complement of H onto the group's angles (weights eliminated);
  stiff  = the first FRAC of the angles by impact |w_b|^2 mass_bt offset^2 (the entries the long schedule freezes
           first), soft = the rest; S_hard = S on the stiff set with the soft angles eliminated (they compensate);
  damage of the stiff decisions d' S_hard d (mHa^2) for
    (a) independent nearest rounding,
    (b) progressive rounding with full compensation in the impact order (Babai on S_hard, first frozen first),
    (c) (b) followed by a joint local search (single +-1 moves on the stiff words, exact in the model) and a
        simulated-annealing pass of SA_SWEEPS sweeps.
  If (c) is far below (b), the stiff words are worth deciding jointly; if not, decoding is not the lever.

    python r48_chem_q7c_linmodel.py INIT.npz [--group G] [--frac 0.46] [--bits 7] [--sa 20]
"""
import argparse
import os
import sys
import time

import numpy as np
import scipy.linalg as sl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_round import N, Target, angles_and_mass, log_depth_tree  # noqa: E402
from r48_chem_q7c_absorb import tangents  # noqa: E402


def babai(S, x, order):
    """Progressive rounding with full compensation: entries decided in `order` (first = decided first), each
    rounded to the nearest integer of its position corrected by the decisions made before it, the undecided ones
    free. Returns the integer vector z (same indexing as x)."""
    rev = order[::-1]                       # Cholesky order: the first decided entry is the LAST one
    Sp = S[np.ix_(rev, rev)]
    Rm = sl.cholesky(Sp, lower=False)       # Sp = R'R, R upper triangular
    xp = x[rev]
    n = len(xp)
    z = np.zeros(n)
    e = np.zeros(n)                         # z - x in the permuted order, filled from the end
    for k in range(n - 1, -1, -1):
        corr = Rm[k, k + 1:] @ e[k + 1:] / Rm[k, k] if k + 1 < n else 0.0
        z[k] = np.round(xp[k] - corr)
        e[k] = z[k] - xp[k]
    out = np.zeros(n)
    out[rev] = z
    return out


def staged(S, x, c, fracs, rank, sk=None, info=None):
    """The real tool's procedure in the model: at every stage the entries still free sit at their conditional
    means (converged refits); the next ones are chosen by `rank` ('impact' = c * (offset of the present mean)^2,
    re-ranked every stage; 'stiff' = c; 'sdiag' = sk) and frozen at their nearest integers. Exact Gaussian
    conditioning: P = covariance of the free entries (S^-1 restricted), mu = their means.
    Returns (total damage, number of stages)."""
    n = len(x)
    # r47: a ridge (1e-9 of the mean diagonal), so that a table with twin rows (S singular on the pair's difference) factors
    P = sl.cho_solve(sl.cho_factor(S + 1e-9 * np.trace(S) / n * np.eye(n)), np.eye(n))
    mu = x.copy()
    free = np.arange(n)
    E = 0.0
    done = 0
    for f in fracs:
        need = int(round(f * n)) - done
        if need <= 0:
            continue
        m = mu[free]
        off = m - np.round(m)
        if rank == "impact":
            score = c[free] * off * off
        elif rank == "stiff":
            score = c[free].copy()
        elif rank == "cond":           # conditional stiffness given every other free entry free: 1 / P_kk
            score = 1.0 / np.maximum(np.diag(P), 1e-300)
        elif rank == "sdiag_off":      # S_kk times the present offset squared
            score = sk[free] * off * off
        else:
            score = sk[free].copy()
        pick = np.argsort(-score, kind="stable")[:need]
        keep = np.setdiff1d(np.arange(len(free)), pick)
        e = np.round(m[pick]) - m[pick]                   # the kick of the frozen batch
        Pdd = P[np.ix_(pick, pick)]
        cf = sl.cho_factor(Pdd + 1e-14 * np.trace(Pdd) / len(pick) * np.eye(len(pick)))
        sol = sl.cho_solve(cf, e)
        E += float(e @ sol)
        if info is not None:        # r47: the damage word by word (e_k sol_k sums to the stage's damage)
            info["contrib"][free[pick]] = e * sol
            info["share"][free[pick]] = f
            info["word"][free[pick]] = np.round(m[pick])
        Pud = P[np.ix_(keep, pick)]
        mu[free[keep]] += Pud @ sol
        P = P[np.ix_(keep, keep)] - Pud @ sl.cho_solve(cf, Pud.T)
        free = free[keep]
        done += len(pick)
    return E, len(fracs)


def local_search(S, x, z, sweeps, rng, temp0):
    """Single +-1 moves on z, exact energy d'Sd with d = z - x: greedy descent, then simulated annealing (best
    state kept), then greedy again."""
    d = z - x
    g = S @ d
    E = float(d @ g)
    diag = np.diag(S).copy()

    def greedy(d, g, E):
        while True:
            dm = -2 * g + diag                    # move -1
            dp = 2 * g + diag                     # move +1
            k = int(np.argmin(np.minimum(dm, dp)))
            best = min(dm[k], dp[k])
            if best >= -1e-18:
                return d, g, E
            s = 1.0 if dp[k] < dm[k] else -1.0
            d[k] += s
            g += s * S[:, k]
            E += best
    d, g, E = greedy(d, g, E)
    bestd, bestE = d.copy(), E
    n = len(d)
    for sw in range(sweeps):
        T = temp0 * (1.0 - sw / max(sweeps, 1))
        for k in rng.permutation(n):
            s = 1.0 if rng.random() < 0.5 else -1.0
            dE = 2 * s * g[k] + diag[k]
            if dE < 0 or (T > 0 and rng.random() < np.exp(-dE / T)):
                d[k] += s
                g += s * S[:, k]
                E += dE
                if E < bestE:
                    bestE, bestd = E, d.copy()
    d = bestd.copy()
    g = S @ d
    d, g, E = greedy(d, g, float(d @ g))
    return d, E


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("init")
    ap.add_argument("--group", type=int, default=-1)
    ap.add_argument("--frac", type=float, default=0.46)
    ap.add_argument("--bits", type=int, default=7)
    ap.add_argument("--sa", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tree", default="log", help="log | random:<seed> | <edges.npy>")
    ap.add_argument("--staged", default="", help="long | full | fine2: the tool's staged procedure (dynamic impact etc.)")
    ap.add_argument("--ranks", default="", help="comma list for --staged: impact, stiff, sdiag, cond, sdiag_off")
    ap.add_argument("--rows", type=int, default=0, help="r47: sub-model with the ROWS heaviest rows of the group as variables")
    ap.add_argument("--classes", default="", help="r47: with --staged: the damage by word class (1 = print; a path = also save npz)")
    ap.add_argument("--order", default="impact", choices=("impact", "stiff", "reverse", "random", "sdiag"),
                    help="order of the progressive decisions on the stiff set")
    a = ap.parse_args()
    t0 = time.time()
    z = np.load(a.init)
    U, W = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
    R, B, C = W.shape
    U /= np.linalg.norm(U, axis=2, keepdims=True)
    tgt = Target()
    I, J = tgt.i, tgt.j
    s = np.sqrt(np.where(I == J, 1.0, 2.0))
    A = np.einsum("rbc,rbp->rcp", W, U[:, :, I] * U[:, :, J])                    # (R, C, P)
    wn2 = np.sum(W ** 2, axis=2)
    r = int(np.argmax(wn2.sum(axis=1))) if a.group < 0 else a.group
    edges = log_depth_tree(N)
    if a.tree.startswith("random:"):
        perm = np.concatenate([[0], 1 + np.random.default_rng(int(a.tree[7:])).permutation(N - 1)])
        edges = perm[edges]
    elif a.tree != "log":
        edges = np.load(a.tree).astype(np.int64)
    Ur, Wr, wn2r = U[r], W[r], wn2[r]
    if a.rows and a.rows < B:      # r47: a sub-model: only the heaviest rows of the group are variables (the others
        keep_ = np.sort(np.argsort(-wn2r, kind="stable")[:a.rows])     # stay in the table A, fixed and unrounded)
        Ur, Wr, wn2r = Ur[keep_], Wr[keep_], wn2r[keep_]
        B = a.rows
    u, V = tangents(Ur, edges)                                                    # (B, 76), (B, 75, 76) per radian
    T_ = V.shape[1]
    p = (u[:, I] * u[:, J]) * s
    q = ((V[:, :, I] * u[:, None, J] + u[:, None, I] * V[:, :, J]) * s).reshape(B * T_, -1)
    g = (Wr @ A[r]) * s
    ar = A[r] * s
    bk = np.repeat(np.arange(B), T_)
    # Gram matrices
    Gqq = q @ q.T
    Gqg = q @ g.T                      # (na, B)
    Ggg = g @ g.T
    Gqp, Gqa = q @ p.T, q @ ar.T       # (na, B), (na, C)
    Ggp, Gga = g @ p.T, g @ ar.T       # (B, B), (B, C)
    Gpp, Gaa, Gap = p @ p.T, ar @ ar.T, ar @ p.T
    na, nw = B * T_, B * C
    # angle-angle: 2[(g_bk.g_bl)(q_k.q_l) + (g_bk.q_l)(q_k.g_bl)]
    X = Gqg[:, bk]                     # X[k, l] = q_k . g_{b(l)}
    Haa = 2.0 * (Ggg[np.ix_(bk, bk)] * Gqq + X.T * X)
    del Gqq, X
    # weights indexed (b, c) -> b * C + c ; angle-weight: x1 = g_bk, y1 = q_k ; x2 = a_c, y2 = p_b'
    wb, wc = np.repeat(np.arange(B), C), np.tile(np.arange(C), B)
    Haw = 2.0 * (Gga[np.ix_(bk, wc)] * Gqp[:, wb] + Ggp[np.ix_(bk, wb)] * Gqa[:, wc])
    # weight-weight: x = a_c, y = p_b
    Hww = 2.0 * (Gaa[np.ix_(wc, wc)] * Gpp[np.ix_(wb, wb)] + Gap[np.ix_(wc, wb)] * Gap[np.ix_(wc, wb)].T)
    lam = 1e-10 * np.trace(Hww) / nw
    cw = sl.cho_factor(Hww + lam * np.eye(nw))
    S = Haa - Haw @ sl.cho_solve(cw, Haw.T)
    S = 0.5 * (S + S.T)
    del Haa, Haw, Hww
    step = np.pi / (1 << (a.bits - 1))           # radians per step
    S *= step * step * 1e6                       # energy in mHa^2 per (steps)^2
    th, mass = angles_and_mass(Ur, edges)
    x = (th * (1 << (a.bits - 1))).reshape(-1)   # positions in steps
    off = np.abs(x - np.round(x))
    impact = (wn2r[:, None] * mass).reshape(-1) * off ** 2
    if a.staged:                                  # the tool's staged procedure, every angle decided
        from r48_chem_q7c_round import FULL, LONG
        fr = {"long": LONG, "full": FULL}.get(a.staged)
        if fr is None:                            # fine2: every long stage halved
            fr, prev = [], 0.0
            for f_ in LONG:
                fr += [prev + 0.5 * (f_ - prev), f_]
                prev = f_
        cvec = (wn2r[:, None] * mass).reshape(-1)
        for rk in (a.ranks.split(",") if a.ranks else ("impact", "stiff", "sdiag")):
            info = dict(contrib=np.zeros(na), share=np.zeros(na), word=np.zeros(na)) if a.classes else None
            Es, ns = staged(S, x, cvec, fr, rk, np.diag(S).copy(), info)
            print(f"  staged {a.staged} ({ns} stages), rank {rk}: {Es:.2f} mHa^2 ({time.time() - t0:.0f} s)", flush=True)
            if info is not None:    # r47: where the model's damage sits
                d = info["contrib"]
                quarter = 1 << (a.bits - 2)
                triv = (np.round(info["word"]) % quarter) == 0
                ms = mass.reshape(-1)
                rrank = np.empty(B, np.int64)
                rrank[np.argsort(-wn2r)] = np.arange(B)
                rr = np.repeat(rrank, T_)
                lev = np.tile(np.searchsorted(np.array([1, 3, 7, 15, 31, 63]), np.arange(T_), side="right"), B)
                sh = info["share"]

                def line(name, parts):
                    print("    " + name + ": " + "  ".join(f"{k} {d[m_].sum():.1f} ({int(m_.sum())})" for k, m_ in parts), flush=True)
                line("kind", [("trivial", triv), ("generic", ~triv)])
                wq = np.round(info["word"]) % (2 * quarter)                      # distance to 0 or a quarter turn, steps
                dq = np.minimum(wq % quarter, quarter - wq % quarter)
                line("generic by distance to a trivial word (steps)", [("1", dq == 1), ("2", dq == 2), ("3", dq == 3), ("4-6", (dq >= 4) & (dq <= 6)), ("7+", dq >= 7)])
                line("the same, mass > 0.05", [("1", (dq == 1) & (ms > 0.05)), ("2", (dq == 2) & (ms > 0.05)), ("3", (dq == 3) & (ms > 0.05)), ("4-6", (dq >= 4) & (dq <= 6) & (ms > 0.05)), ("7+", (dq >= 7) & (ms > 0.05))])
                cv = cvec * (step * step * 1e6) / 12.0
                print("    expected kick c/12 (same units x1e0) of those classes, mass > 0.05: " + "  ".join(
                    f"{k} {cv[m_ & (ms > 0.05)].sum():.0f}" for k, m_ in [("1", dq == 1), ("2", dq == 2), ("3", dq == 3), ("4-6", (dq >= 4) & (dq <= 6)), ("7+", dq >= 7)]), flush=True)
                line("mass", [(">0.5", ms > 0.5), ("0.05-0.5", (ms > 0.05) & (ms <= 0.5)), ("0.005-0.05", (ms > 0.005) & (ms <= 0.05)), ("<0.005", ms <= 0.005)])
                line("generic by mass", [(">0.5", ~triv & (ms > 0.5)), ("0.05-0.5", ~triv & (ms > 0.05) & (ms <= 0.5)), ("<0.05", ~triv & (ms <= 0.05))])
                line("row rank", [("0-5", rr < 6), ("6-24", (rr >= 6) & (rr < 25)), ("25-49", (rr >= 25) & (rr < 50)), ("50-99", (rr >= 50) & (rr < 100)), ("100+", rr >= 100)])
                line("level", [(str(l_), lev == l_) for l_ in range(7)])
                line("frozen share", [("<1%", sh <= 0.01), ("1-4%", (sh > 0.01) & (sh <= 0.04)), ("4-10%", (sh > 0.04) & (sh <= 0.10)), ("10-25%", (sh > 0.10) & (sh <= 0.25)), ("25%+", sh > 0.25)])
                if a.classes != "1":
                    np.savez(a.classes, contrib=d, share=sh, word=info["word"], mass=ms, rowrank=rr, level=lev, x=x)
        return
    order = np.argsort(-impact, kind="stable")
    ns = int(round(a.frac * na))
    stiff, soft = np.sort(order[:ns]), np.sort(order[ns:])
    Sss = S[np.ix_(stiff, stiff)]
    if len(soft):
        Ssf, Sff = S[np.ix_(stiff, soft)], S[np.ix_(soft, soft)]
        lf = 1e-12 * np.trace(Sff) / len(soft)
        Sh = Sss - Ssf @ sl.cho_solve(sl.cho_factor(Sff + lf * np.eye(len(soft))), Ssf.T)
    else:                                        # --frac 1: every angle decided, only the weights compensate at the end
        Sh = Sss.copy()
    Sh = 0.5 * (Sh + Sh.T)
    xs = x[stiff]
    # (a) nearest
    dn = np.round(xs) - xs
    Ea = float(dn @ Sh @ dn)
    # (b) progressive in the impact order (restricted to the stiff set)
    rank = np.empty(na, np.int64)
    rank[order] = np.arange(na)
    ord_s = np.argsort(rank[stiff], kind="stable")
    if a.order == "stiff":
        ord_s = np.argsort(-(wn2r[:, None] * mass).reshape(-1)[stiff], kind="stable")
    elif a.order == "reverse":
        ord_s = ord_s[::-1]
    elif a.order == "sdiag":          # weights-compensated diagonal S_kk, largest first
        ord_s = np.argsort(-np.diag(S)[stiff], kind="stable")
    elif a.order == "random":
        ord_s = np.random.default_rng(a.seed).permutation(len(stiff))
    zb = babai(Sh, xs, ord_s)
    db = zb - xs
    Eb = float(db @ Sh @ db)
    # (c) joint local search + annealing from (b)
    rng = np.random.default_rng(a.seed)
    dc, Ec = local_search(Sh, xs, zb.copy(), a.sa, rng, temp0=0.3 * float(np.median(np.diag(Sh))))
    raw = float(dn @ (Sss @ dn))
    moved = int(np.sum(np.abs((dc + xs) - zb) > 0.5))
    print(f"group {r} (sum|w|^2 {wn2[r].sum():.3f}), {na} angles, stiff {ns} ({a.frac:.2f}), {a.bits} bits; "
          f"{time.time() - t0:.0f} s")
    print(f"  stiff words at nearest, nothing compensating but the weights: {raw:.2f} mHa^2")
    print(f"  (a) nearest, soft angles compensating:           {Ea:.2f} mHa^2")
    print(f"  (b) progressive with compensation (impact order): {Eb:.2f} mHa^2  ({Eb / Ea:.3f} of a)")
    print(f"  (c) (b) + joint local search and annealing:       {Ec:.2f} mHa^2  ({Ec / Eb:.3f} of b), {moved} words moved")
    w_ = np.linalg.eigvalsh(Sh)
    print(f"  S_hard eigenvalues: min {w_[0]:.3e} median {np.median(w_):.3e} max {w_[-1]:.3e}; "
          f"diag median {np.median(np.diag(Sh)):.3e}; mean |offdiag| / diag {np.mean(np.abs(Sh - np.diag(np.diag(Sh)))) * len(Sh) / np.trace(Sh):.3f}")


if __name__ == "__main__":
    main()
