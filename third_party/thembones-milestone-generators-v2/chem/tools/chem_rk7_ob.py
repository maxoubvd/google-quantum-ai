"""chem_rk7_ob (RK7): the 76 one-body rows under GENERAL RELATIONS (a relation file, chem_rk7_rel).

progressive_rel: the designer's staged one-body rounding (as chem_wr7_obb.progressive_dynb) with chem_rk7_rel.RelRank:
  the stiffest entries freeze first; a word is forced only when the words still free cannot pay (rank rule).
polish_rel: exact coordinate descent on error_1 (polishb's objective) over candidate moves of one row: every word +-1,
  +-2, +-4, and every small KERNEL SET of the relation system (2 to 4 words whose bit-0 flips keep every relation;
  enumerated from the bit-0 syndromes) with all +-1 sign patterns; EVERY candidate is filtered by an explicit check
  of all relations on its row (so the moves are exact for relations on any bit).
one_try_rel: one rounding (fine schedule, optional random placement), polish, relations checked on the result.

    python chem_rk7_ob.py TABLE.npz --rel @rel.txt [--fine 3] [--tries 1] [--workers 1] [--range full|signed]
"""
import argparse
import itertools
import multiprocessing as mp
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_onebody import STAGES, angles_and_mass, decode_jac, decode_pi, error1  # noqa: E402
from chem_rk7_rel import RelRank, check, columns, fmt, parse, ttable  # noqa: E402


def fine_stages(stages, k):
    out, prev = [], 0.0
    for f in stages:
        out += [round(prev + (f - prev) * j / k, 7) for j in range(1, k + 1)]
        prev = f
    return tuple(out)


class RelFilter:
    """Explicit relation check of candidate rows (any bits): parity per relation = sum_b bits_b @ M_b^T mod 2."""

    def __init__(self, rels, ncol=75):
        self.rels = rels
        bits = sorted({b for (_, o, _) in rels for (_, b) in o})
        self.M = {}
        for b in bits:
            M = np.zeros((len(rels), ncol), np.int64)
            for i, (k, o, _) in enumerate(rels):
                for (j, bb) in o:
                    if bb == b:
                        M[i, 5 * k + j] ^= 1
                        M[i, 5 * k - 5 + j] ^= 1
            self.M[b] = M

    def parities(self, rows):
        rows = np.mod(np.asarray(rows, np.int64), 128)
        acc = np.zeros((rows.shape[0], len(self.rels)), np.int64)
        for b, M in self.M.items():
            acc += ((rows >> b) & 1) @ M.T
        return acc & 1


def kernel_sets(rels, ncol=75, maxsize=4):
    """Small sets of words whose bit-0 flips keep every relation (from the bit-0 syndromes): (free columns, sets)."""
    syn = [0] * ncol
    for i, (k, o, _) in enumerate(rels):
        for (j, b) in o:
            if b == 0:
                syn[5 * k + j] ^= 1 << i
                syn[5 * k - 5 + j] ^= 1 << i
    rel_cols = set(columns(rels))
    free = [c for c in range(ncol) if c not in rel_cols]
    cols = [c for c in range(ncol) if c in rel_cols and syn[c]]
    sets = set()
    by = {}
    for c in cols:
        by.setdefault(syn[c], []).append(c)
    for s, cs in by.items():           # pairs with equal syndrome
        for a_, b_ in itertools.combinations(cs, 2):
            sets.add((a_, b_))
    pairsyn = {}
    for a_, b_ in itertools.combinations(cols, 2):
        pairsyn.setdefault(syn[a_] ^ syn[b_], []).append((a_, b_))
    if maxsize >= 3:
        for (a_, b_) in itertools.combinations(cols, 2):
            for c in by.get(syn[a_] ^ syn[b_], []):
                if c not in (a_, b_):
                    sets.add(tuple(sorted((a_, b_, c))))
    if maxsize >= 4:
        for s, ps in pairsyn.items():
            if s == 0 or len(ps) > 400:
                continue
            for (p, q) in itertools.combinations(ps, 2):
                u = set(p) | set(q)
                if len(u) == 4:
                    t = tuple(sorted(u))
                    if not any(tuple(sorted(x)) in sets for x in itertools.combinations(t, 2)) and \
                            not any(tuple(sorted(x)) in sets for x in itertools.combinations(t, 3)):
                        sets.add(t)
    for b in range(ncol // 5):          # whole batches (five words): the kernel moves of all-equal / wire layouts
        t = tuple(c for c in range(5 * b, 5 * b + 5) if c in rel_cols)
        if len(t) > 4:
            s_ = 0
            for c in t:
                s_ ^= syn[c]
            if s_ == 0:
                sets.add(t)
    return free, sorted(sets, key=lambda t: (len(t), t))


def progressive_rel(th0, edges, lam_at, H, scale, iters, rels, stages=STAGES, lo=-64, hi=63, zone=0):
    nrow = th0.shape[0]
    dyn = RelRank(rels, np.arange(nrow), nrow, zone)
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

    for frac in stages:
        need = int(round(frac * theta.size)) - int(frozen.sum())
        if need > 0:
            now = np.where(frozen, qv, theta)
            mass = angles_and_mass(decode_pi(now, edges), edges)[1]
            xs = now * scale
            if hi - lo + 1 >= 128:
                q = np.mod(np.round(xs) - lo, hi - lo + 1) + lo
            else:
                q = np.clip(np.round(xs), lo, hi)
            impact = lam2 * mass * (q / scale - now) ** 2
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            for flat in pick:
                j, c = divmod(int(flat), ncol)
                zc = int(q[j, c])
                v = dyn.choose(j, c, float(xs[j, c]), lo, hi)
                if v is not None:
                    zc = int(v)
                dyn.commit(j, c, zc)
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


PAR32 = np.array([bin(x).count("1") & 1 for x in range(32)], np.int64)
BITS32 = (np.arange(32)[:, None] >> np.arange(5)[None, :]) & 1          # (pattern, slot)


def dp_reparity(cost, by_k, tgt):
    """Viterbi over the 15 batches: state = bit-0 pattern of a batch (32); a transition b-1 -> b must keep every bit-0
    relation of loader b (parity(slotmask AND (P_{b-1} XOR P_b)) = target). cost (75, 2): the cost of word c taking
    bit 0 = p. -> (pattern per batch (15,), total cost)."""
    nb = cost.shape[0] // 5
    sl = np.arange(5)[None, :]
    best = cost[0:5][sl, BITS32].sum(axis=1)
    xs = np.arange(32)[:, None] ^ np.arange(32)[None, :]
    back = []
    for b in range(1, nb):
        ok = np.ones((32, 32), bool)
        for (sm, ri) in by_k.get(b, ()):
            ok &= PAR32[xs & sm] == tgt[ri]
        tot = np.where(ok, best[:, None], np.inf)
        arg = np.argmin(tot, axis=0)
        best = tot[arg, np.arange(32)] + cost[5 * b:5 * b + 5][sl, BITS32].sum(axis=1)
        back.append(arg)
    p = int(np.argmin(best))
    tc = float(best[p])
    pats = [p]
    for arg in reversed(back):
        p = int(arg[p])
        pats.append(p)
    return np.array(pats[::-1], np.int64), tc


def polish_rel(wq, edges, H, brot, sweeps, rels, addrs, inner=25, lo=-64, hi=63, maxsize=4, dp=True, zone=0):
    """Coordinate descent on error_1 with relation-keeping moves; wq rows at addresses `addrs` (targets per row).
    Each row first tries the Viterbi re-parity move (bit-0 relations), then the local moves; the best state by the
    verifier's error_1 (rows re-paired by overlap) over the sweeps is returned."""
    scale = float(1 << (brot - 1))
    N = edges.shape[0] + 1
    ncol = wq.shape[1]
    wq = np.array(wq, np.int64)
    filt = RelFilter(rels, ncol)
    target = np.stack([ttable(m, addrs) for (_, _, m) in rels], axis=1) if rels else np.zeros((len(wq), 0), np.int64)
    free, ksets = kernel_sets(rels, ncol, maxsize)
    full = hi - lo + 1 >= 128
    dp = dp and any(all(b == 0 for (_, b) in o) for (_, o, _) in rels)
    by_k = {}                       # the Viterbi move speaks of the pure bit-0 relations only
    for ri, (k, o, _) in enumerate(rels):
        if all(b == 0 for (_, b) in o):
            by_k.setdefault(k, []).append((sum(1 << j for (j, _) in o), ri))
    hib = np.zeros(ncol, np.int64)  # per column: the highest relation bit above 0 (its moves must keep bits 1..hib)
    for (k, o, _) in rels:
        for (j, b) in o:
            for c_ in (5 * k + j, 5 * k - 5 + j):
                hib[c_] = max(hib[c_], b)
    keepm = ((1 << hib) - 1)        # bits 1..hib of a word, shifted down by one
    sdel = np.array([d for d in range(-3, 4) if d], np.int64)
    dpn = dpa = 0
    moves = []                                   # (columns, deltas)
    for t in range(ncol):
        for d in (1, -1, 2, -2, 4, -4):
            moves.append(((t,), (d,)))
    for S in ksets:
        for sg in itertools.product((1, -1), repeat=len(S)):
            moves.append((S, sg))
    Dm = np.zeros((len(moves), ncol), np.int64)
    for i, (S, sg) in enumerate(moves):
        for c, d in zip(S, sg):
            Dm[i, c] += d
    V = decode_pi(np.mod(wq, 1 << brot) / scale, edges)
    e1, perm, D, lam = error1(V, H)
    best_e, best_wq = e1, wq.copy()
    stats = dict(kernel_sets=len(ksets), free_cols=len(free), moves=len(moves), moved=0, sweeps=0, e_in=float(e1))

    def rowval(Vc, base, lj):
        quad = np.einsum("ip,pq,iq->i", Vc, base, Vc)
        n4 = np.sum(Vc * Vc, axis=1) ** 2
        return 2 * lj * quad + lj * lj * n4

    for sweep in range(sweeps):
        moved = 0
        for j in range(N):
            lj = lam[perm[j]]
            base = D - lj * np.outer(V[j], V[j])
            if dp and addrs[j] >= zone:      # the Viterbi re-parity move: every word's best value of each parity, the row's best pattern
                v0 = V[j]
                cur = 2 * lj * float(v0 @ base @ v0) + lj * lj * float(v0 @ v0) ** 2
                Cs = np.repeat(wq[j][None, :], ncol * len(sdel), axis=0)
                ci = np.repeat(np.arange(ncol), len(sdel))
                Cs[np.arange(len(Cs)), ci] += np.tile(sdel, ncol)
                if full:
                    Cs = np.mod(Cs - lo, 128) + lo
                    okc = np.ones(len(Cs), bool)
                else:
                    okc = (Cs[np.arange(len(Cs)), ci] >= lo) & (Cs[np.arange(len(Cs)), ci] <= hi)
                vals = np.full(len(Cs), np.inf)
                vals[okc] = rowval(decode_pi(np.mod(Cs[okc], 1 << brot) / scale, edges), base, lj) - cur
                vals = vals.reshape(ncol, len(sdel))
                newv = Cs[np.arange(len(Cs)), ci].reshape(ncol, len(sdel))
                if np.any(hib):         # columns with relation bits above 0: only moves that keep those bits
                    sameh = (((np.mod(newv, 128) >> 1) ^ (np.mod(wq[j], 128)[:, None] >> 1)) & keepm[:, None]) == 0
                    vals = np.where(sameh, vals, np.inf)
                par = np.mod(newv, 2)
                cost = np.zeros((ncol, 2))
                pick = np.zeros((ncol, 2), np.int64)
                p0 = np.mod(wq[j], 2)
                for p in (0, 1):
                    vv = np.where(par == p, vals, np.inf)
                    a_ = np.argmin(vv, axis=1)
                    cost[:, p] = vv[np.arange(ncol), a_]
                    pick[:, p] = newv[np.arange(ncol), a_]
                    same = p0 == p              # the present word itself: cost 0 at its own parity
                    better = same & (cost[:, p] > 0)
                    cost[better, p] = 0.0
                    pick[better, p] = wq[j][better]
                pats, tc = dp_reparity(cost, by_k, target[j])
                newbits = BITS32[pats].reshape(-1)[:ncol]
                if np.any(newbits != p0) and tc < -1e-13:
                    dpn += 1
                    c = wq[j].copy()
                    ch = newbits != p0
                    c[ch] = pick[ch, newbits[ch]]
                    if not rels or np.all(filt.parities(c[None, :])[0] == target[j]):
                        Vn = decode_pi(np.mod(c[None, :], 1 << brot) / scale, edges)
                        if rowval(Vn, base, lj)[0] < cur - 1e-13:
                            wq[j] = c
                            V[j] = Vn[0]
                            dpa += 1
                            moved += 1
            for _ in range(inner):
                C = wq[j][None, :] + Dm
                if full:
                    C = np.mod(C - lo, 128) + lo
                    ok = np.ones(len(C), bool)
                else:
                    ok = np.all((C >= lo) & (C <= hi), axis=1)
                if rels and addrs[j] >= zone:      # exempt rows (R53-CHEM-EXEMPT's zone) need no relation
                    ok &= np.all(filt.parities(C) == target[j][None, :], axis=1)
                C = C[ok]
                if not len(C):
                    break
                Vc = decode_pi(np.mod(C, 1 << brot) / scale, edges)
                v0 = V[j]
                quad = np.einsum("ip,pq,iq->i", Vc, base, Vc)
                n4 = np.sum(Vc * Vc, axis=1) ** 2
                val = 2 * lj * quad + lj * lj * n4
                cur = 2 * lj * float(v0 @ base @ v0) + lj * lj * float(v0 @ v0) ** 2
                i = int(np.argmin(val))
                if val[i] >= cur - 1e-13:
                    break
                wq[j] = C[i]
                V[j] = Vc[i]
                moved += 1
            D = base + lj * np.outer(V[j], V[j])
        e_prev = e1
        e1, perm, D, lam = error1(V, H)
        stats["moved"] += moved
        stats["sweeps"] = sweep + 1
        if e1 < best_e:
            best_e, best_wq = e1, wq.copy()
        if not moved or e_prev - e1 < 1e-7 * e1:
            break
    stats.update(dp_tried=dpn, dp_taken=dpa, e_out=float(best_e))
    return best_wq, best_e, stats


def one_try_rel(args):
    """(H, edges, brot, iters, sweeps, rels, t, seed, fine, lo, hi[, btab]) -> (error_1, words in address order, info).
    btab: on tries t > 0 every target-0 relation takes a random 'B:x' table (free per block on blocks 0..3 = addresses
    0..63, 0 elsewhere: in the family F and 0 on every two-body address), so the one-body rows get an address choice;
    info['rels'] = the relation file actually kept (the chain checks the table against it)."""
    H, edges, brot, iters, sweeps, rels, t, seed, fine, lo, hi = args[:11]
    btab = bool(args[11]) if len(args) > 11 else False
    zone = int(args[12]) if len(args) > 12 else 0
    if btab and t > 0:
        rb = np.random.default_rng([int(seed), int(t), 77])
        rels = [(k, o, (f"B:{int(rb.integers(0, 16))}" if (not isinstance(m, str) and m == 0) else m)) for (k, o, m) in rels]
    from r48_chem_q7c_onebody import reconstruct_u_vectors
    scale = float(1 << (brot - 1))
    lam_e, vec = np.linalg.eigh(H)
    if hi < (1 << (brot - 1)) - 1:
        from r48_chem_q7c_obsigned import signed_angles
        th0, _ = signed_angles(vec.T, edges)
    else:
        th0, _ = angles_and_mass(vec.T, edges)
    n = th0.shape[0]
    perm = np.arange(n) if t == 0 else np.random.default_rng(1000 * seed + t).permutation(n)
    st = fine_stages(STAGES, fine) if fine > 1 else STAGES
    t0 = time.time()
    wq, dyn = progressive_rel(th0[perm], edges, lam_e[perm], H, scale, iters, rels, st, lo, hi, zone)
    out = np.zeros_like(wq)
    out[dyn.mrow] = wq                       # row -> its address
    bad0 = check(np.mod(out, 128), np.arange(n), rels, zone=zone)
    assert bad0 == 0, f"{bad0} relations broken after the staged rounding"
    e_pre = error1(reconstruct_u_vectors(np.mod(out, 1 << brot) / scale, edges), H)[0]
    out, _, st_ = polish_rel(out, edges, H, brot, max(sweeps, 1), rels, np.arange(n), lo=lo, hi=hi, zone=zone)
    bad = check(np.mod(out, 128), np.arange(n), rels, zone=zone)
    assert bad == 0, f"{bad} relations broken on the one-body rows"
    out = np.mod(out, 1 << brot)
    e = error1(reconstruct_u_vectors(out / scale, edges), H)[0]
    info = dict(free=int(dyn.free_ok), forced=int(dyn.forced_n), completions=int(dyn.ncons), e_staged=float(e_pre),
                seconds=round(time.time() - t0, 1), **st_)
    info["rels"] = fmt(rels)
    return float(e), out, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--rel", nargs="+", required=True, help="relation files (@file); '-' = no relations")
    ap.add_argument("--tries", type=int, default=1)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--iters", type=int, default=300)
    ap.add_argument("--sweeps", type=int, default=40)
    ap.add_argument("--fine", type=int, default=3)
    ap.add_argument("--range", default="full", choices=["full", "signed"])
    ap.add_argument("--save", default="", help="npz with the best rows of the last relation file")
    ap.add_argument("--btab", type=int, default=0, help="1: tries > 0 draw random B:x targets for target-0 relations")
    ap.add_argument("--zone", type=int, default=0, help="R53-CHEM-EXEMPT's zone: addresses below it carry no relation")
    a = ap.parse_args()
    from chem_wr7_obtry import h1_prime
    t0 = time.time()
    H, edges = h1_prime(a.table)
    lo, hi = (-64, 63) if a.range == "full" else (-32, 31)
    for rf in a.rel:
        rels = [] if rf == "-" else parse(rf)
        args = [(H, edges, 7, a.iters, a.sweeps, rels, t, a.seed, a.fine, lo, hi, a.btab, a.zone) for t in range(a.tries)]
        if a.workers > 1 and a.tries > 1:
            with mp.get_context("fork").Pool(max(1, min(a.workers, a.tries))) as pool:
                res = pool.map(one_try_rel, args)
        else:
            res = [one_try_rel(x) for x in args]
        es = np.array([r[0] for r in res]) * 1e3
        b = int(es.argmin())
        print(f"{os.path.basename(a.table)} {rf} ({len(rels)} relations): error_1 best {es.min():.3f} (try {b}), "
              f"median {np.median(es):.3f} mHa over {a.tries} (fine {a.fine}, range {a.range}, btab {a.btab}, zone {a.zone}); all "
              f"{' '.join(f'{x:.2f}' for x in es)}; best: { {k_: v_ for k_, v_ in res[b][2].items() if k_ != 'rels'} }; "
              f"{time.time() - t0:.0f} s", flush=True)
        if a.save:
            np.savez(a.save, rows=res[b][1], e1=res[b][0], rels=res[b][2]["rels"])


if __name__ == "__main__":
    main()
