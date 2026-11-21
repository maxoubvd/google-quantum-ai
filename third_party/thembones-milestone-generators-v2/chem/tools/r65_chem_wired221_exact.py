#!/usr/bin/env python
"""r65_chem_wired221_exact.py: designer R65-CHEM-WIRED221. EXACT SINGLE-FLIP PRICES of a rounded table.

For one changed word of a two-body row b of group r (decoded vector u -> u', dP = u'u'^T - uu^T) the verifier's
error_2^2 = || sum_rc X_rc (x) X_rc - T' ||^2 moves by

    4 <dP, G_b> + 2 w2 <dP (x) dP, Res> + 2 |dP|^2 |Y_b|^2 + 2 <dP, Y_b>^2 + 4 w2 |dP|^2 <dP, Y_b> + w2^2 |dP|^4

with Y_b = sum_c w_brc X_rc, G_b = Res . Y_b, w2 = sum_c w_brc^2 (Res = the verifier's h2_delta). Everything but the
second term is O(N^2) a flip once G is built (one 5,776 x 5,776 product for all rows); the second term is added for
chosen flips only (--full). error_1^2 moves to first order by 2 kappa_b (u'^T F u' - u^T F u), kappa_b = sum_c W_rc
w_brc, F = E1 - sum_j (v_j^T E1 v_j) e_j e_j^T (eigenvalues re-matched).

    python tools/r65_chem_wired221_exact.py SUBMISSION.npz --k 14 [--bits 1,2] [--lmin 3] [--check 12] [--out DIR]

Prints: the check of the formula against the verifier on random flips; the price statistics of every bit-b flip of
the loader's ten columns; the kill test (a) with the EXACT-cheapest live word of every uncovered (block, word).
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.expanduser("~/qcc/data")
sys.path[:0] = [os.path.join(DATA, "kickmix"), os.path.join(DATA, "verifiers")]
os.environ.setdefault("DFTHC_DATA", DATA)
os.environ.setdefault("CHEM_BROT", "7")
os.environ.setdefault("CHEM_HB", "77")
from dfthcsos import verify as V  # noqa: E402
from dfthcsos.scoreinfra import dfthc_math as dm  # noqa: E402
from scipy.optimize import linear_sum_assignment  # noqa: E402

BROT = 7


def decode(words, edges):
    return dm.reconstruct_u_vectors(words / float(1 << (BROT - 1)), edges)


class Pricer:
    def __init__(self, d):
        self.inputs = V.load_inputs({k: d[k] for k in ("U", "W", "H2_I", "H1_I", "Hsym_lower")})
        inp = self.inputs
        self.words = d["angle_words"].astype(np.int64)
        self.edges = d["tree_edges"].astype(np.int64)
        self.w = np.asarray(inp.w, float)                    # (R, B, C)
        self.W = np.asarray(inp.W, float)                    # (R, C)
        self.S = np.asarray(inp.S, float)
        ref = inp.ref
        oei = np.asarray(ref["oei"], float)
        tei = np.asarray(ref["tei"], float)
        N = oei.shape[0]
        self.N = N
        R, B, C = self.w.shape
        self.R, self.B, self.C = R, B, C
        self.n1 = N
        u_all = decode(self.words, self.edges)
        u_all = u_all / np.linalg.norm(u_all, axis=1, keepdims=True)
        self.u_all = u_all
        ob, u2 = u_all[:N], u_all[N:].reshape(R, B, N)
        X = np.einsum("tjp,tjq,tjl->tlpq", u2, u2, self.w, optimize=True)      # (R, C, N, N)
        res = np.einsum("rcpq,rcst->pqst", X, X, optimize=True)
        res -= tei
        dg = np.arange(N)
        res[:, :, dg, dg] -= 0.5 * self.S[:, :, None]
        res[dg, dg, :, :] -= 0.5 * self.S[None, :, :]
        self.e2sq = float(np.sum(res * res))
        self.res2 = res.reshape(N * N, N * N)
        Y = np.einsum("rbc,rcpq->rbpq", self.w, X, optimize=True).reshape(R * B, N * N)   # rows (r, b)
        self.Y = Y
        t0 = time.time()
        self.G = Y @ self.res2                               # (R B, N^2)
        self.tG = time.time() - t0
        self.w2 = np.sum(self.w ** 2, axis=2).reshape(R * B)
        self.Yn = np.sum(Y * Y, axis=1)
        self.kap = np.einsum("rc,rbc->rb", self.W, self.w).reshape(R * B)
        # one-body side
        h1_thc = oei - 0.5 * np.einsum("prrq->pq", tei)
        h1p = dm.build_h1_prime(h1_thc, self.S, inp.beta_1, ref["electrons"], self.W, X)
        lam, E = np.linalg.eigh(h1p)
        ov = np.abs(ob @ E)
        qi, ii = linear_sum_assignment(-ov)
        matched = np.empty_like(lam)
        matched[qi] = lam[ii]
        M = np.einsum("j,jp,jq->pq", matched, ob, ob, optimize=True)
        E1 = M - h1p
        self.e1sq = float(np.sum(E1 * E1))
        q = np.einsum("jp,pq,jq->j", ob, E1, ob, optimize=True)            # per one-body row
        qe = np.empty_like(lam)
        qe[ii] = q[qi]                                                   # per exact eigenvector
        self.F = E1 - (E * qe[None, :]) @ E.T

    def flip_prices(self, col, xormask=None, newwords=None):
        """prices of changing column col of EVERY two-body row (XOR by xormask, or to newwords[rows]):
        returns (dE2 in Ha^2 without the w2 <dPdP,Res> term, dE1sq in Ha^2 first order, uprime)."""
        N, n1 = self.N, self.n1
        w2 = self.words.copy()
        if newwords is not None:
            w2[:, col] = newwords
        else:
            w2[:, col] ^= xormask
        up = decode(w2, self.edges)[n1:]
        up = up / np.linalg.norm(up, axis=1, keepdims=True)
        u = self.u_all[n1:]
        nb = u.shape[0]
        Gm = self.G.reshape(nb, N, N)
        Ym = self.Y.reshape(nb, N, N)
        dG = np.einsum("bp,bpq,bq->b", up, Gm, up) - np.einsum("bp,bpq,bq->b", u, Gm, u)
        dY = np.einsum("bp,bpq,bq->b", up, Ym, up) - np.einsum("bp,bpq,bq->b", u, Ym, u)
        dp2 = 2.0 - 2.0 * np.sum(u * up, axis=1) ** 2
        dE2 = 4 * dG + 2 * dp2 * self.Yn + 2 * dY ** 2 + 4 * self.w2 * dp2 * dY + self.w2 ** 2 * dp2 ** 2
        dF = np.einsum("bp,pq,bq->b", up, self.F, up) - np.einsum("bp,pq,bq->b", u, self.F, u)
        dE1sq = 2 * self.kap * dF
        return dE2, dE1sq, up

    def res_term(self, row2, up):
        """the missing term 2 w2 <dP (x) dP, Res> of one flip (two-body row index row2, new vector up)."""
        u = self.u_all[self.n1 + row2]
        dP = (np.outer(up, up) - np.outer(u, u)).ravel()
        return 2 * self.w2[row2] * float(dP @ (self.res2 @ dP))

    def exact(self, words2):
        m = V.evaluate_words(self.inputs, words2, self.edges, BROT)
        return float(m["error_2"]) ** 2, float(m["error_1"]), float(m["lambda_eff"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, default=14)
    ap.add_argument("--bits", default="1")
    ap.add_argument("--lmin", type=int, default=3)
    ap.add_argument("--check", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    ap.add_argument("--pricesonly", default=None, help="loaders (e.g. 2,4,6): write OUT/prices_kK.npz for --bits and stop")
    ap.add_argument("--evalonly", nargs="*", default=None, help="tables (npz with angle_words): the verifier's evaluation of each")
    ap.add_argument("--scan", default=None, help="moves (e.g. 2,-2,4,-4): exact price of word+move for every live two-body word of the table; the negative ones are listed (a free polish)")
    a = ap.parse_args()
    from r62_chem_cache_dead import dead_mask
    z = np.load(a.table)
    d = {k: np.array(z[k]) for k in z.files}
    t0 = time.time()
    P = Pricer(d)
    words, edges = P.words, P.edges
    A, ncol = words.shape
    N, n1 = P.N, P.n1
    print(f"pricer built in {time.time() - t0:.1f} s (G in {P.tG:.1f} s): error_2^2 {P.e2sq * 1e6:.4f} mHa^2, "
          f"error_1 {np.sqrt(P.e1sq) * 1e3:.5f} mHa", flush=True)
    e2b, e1b, lb = P.exact(words)
    print(f"verifier: error_2^2 {e2b * 1e6:.4f} mHa^2, error_1 {e1b * 1e3:.5f} mHa, lambda_eff {lb:.6f}", flush=True)
    if a.evalonly is not None:
        for f in a.evalonly:
            w2 = np.load(f)["angle_words"].astype(np.int64)
            e2, e1, lam = P.exact(w2)
            print(f"EVAL {f}: words changed {int((w2 != words).sum())}; d(error_2^2) {(e2 - e2b) * 1e6:+.5f} mHa^2 "
                  f"(error_2 {np.sqrt(e2) * 1e3:.5f} mHa); error_1 {e1 * 1e3:.5f} mHa; lambda_eff {lam:.6f}", flush=True)
        return
    if a.scan:
        dead = dead_mask(words, edges, BROT)
        moves = [int(x) for x in a.scan.split(",")]
        best = np.full((A - n1, ncol), np.inf)
        bmv = np.zeros((A - n1, ncol), np.int64)
        for c in range(ncol):
            w = words[:, c]
            low6 = w & 63
            car = ((w >> 6) ^ (w >> 5)) & 1
            for mv in moves:
                nl = (low6 + mv) % 64
                nw = nl | ((((nl >> 5) & 1) ^ car) << 6)
                dE2 = P.flip_prices(c, newwords=nw)[0] * 1e6
                better = dE2 < best[:, c]
                best[better, c] = dE2[better]
                bmv[better, c] = mv
        best[dead[n1:]] = np.inf
        neg = best < 0
        print(f"SCAN moves {moves}: live two-body words {int((~dead[n1:]).sum())}; words with a NEGATIVE-price move "
              f"{int(neg.sum())}, sum {best[neg].sum():+.4f} mHa^2; per batch "
              f"{[round(float(best[:, 5 * k:5 * k + 5][neg[:, 5 * k:5 * k + 5]].sum()), 4) for k in range(15)]}", flush=True)
        rb = best.min(axis=1)
        print(f"   one move a row (the best): rows with a negative move {int((rb < 0).sum())}, sum {rb[rb < 0].sum():+.4f} mHa^2",
              flush=True)
        # apply one move a row and evaluate exactly
        w2 = words.copy()
        for r in np.nonzero(rb < 0)[0]:
            c = int(np.argmin(best[r]))
            w = int(words[n1 + r, c])
            car = ((w >> 6) ^ (w >> 5)) & 1
            nl = ((w & 63) + int(bmv[r, c])) % 64
            w2[n1 + r, c] = nl | ((((nl >> 5) & 1) ^ car) << 6)
        e2, e1, lam = P.exact(w2)
        print(f"   applied (one a row, {int((w2 != words).sum())} words): EXACT d(error_2^2) {(e2 - e2b) * 1e6:+.4f} mHa^2; "
              f"error_1 {e1 * 1e3:.5f} mHa; lambda_eff {lam:.6f}", flush=True)
        if a.out:
            os.makedirs(a.out, exist_ok=True)
            np.savez(os.path.join(a.out, "scan_polish.npz"), angle_words=w2, tree_edges=edges, best=best, bmv=bmv)
        return
    if a.pricesonly:
        os.makedirs(a.out, exist_ok=True)
        bl = [int(x) for x in a.bits.split(",")]
        for kk in [int(x) for x in a.pricesonly.split(",")]:
            cl = [5 * kk + j for j in range(5)] + ([5 * (kk - 1) + j for j in range(5)] if kk > 0 else [])
            arr = np.stack([np.stack([P.flip_prices(c, xormask=(1 << b) | ((1 << 6) if b == 5 else 0))[0] * 1e6
                                      for c in cl]) for b in bl])
            np.savez(os.path.join(a.out, f"prices_k{kk}.npz"), cols=np.array(cl), bits=np.array(bl), dE2=arr, n1=P.n1)
            print(f"prices_k{kk}.npz written: bits {bl}, {arr.shape}", flush=True)
        return
    dead = dead_mask(words, edges, BROT)
    J, k = 5, a.k
    rng = np.random.default_rng(a.seed)
    cols = [J * k + j for j in range(J)] + ([J * (k - 1) + j for j in range(J)] if k > 0 else [])
    bits = [int(x) for x in a.bits.split(",")]
    price = {}
    for b in bits:
        for c in cols:
            dE2, dE1sq, up = P.flip_prices(c, xormask=(1 << b))
            price[(c, b)] = (dE2 * 1e6, dE1sq, up)
    # ---- check against the verifier on random live flips
    for i in range(a.check):
        b = bits[i % len(bits)]
        c = cols[int(rng.integers(len(cols)))]
        r2 = int(rng.integers(A - n1))
        w2 = words.copy()
        w2[n1 + r2, c] ^= (1 << b)
        e2, e1, _ = P.exact(w2)
        dE2, dE1sq, up = price[(c, b)]
        extra = P.res_term(r2, up[r2]) * 1e6
        pe1 = (np.sqrt(max(P.e1sq + dE1sq[r2], 0.0)) - np.sqrt(P.e1sq)) * 1e3
        print(f"check row {n1 + r2} col {c} bit {b} dead {bool(dead[n1 + r2, c])}: formula {dE2[r2]:+.6f} + res term "
              f"{extra:+.6f} = {dE2[r2] + extra:+.6f}; verifier {(e2 - e2b) * 1e6:+.6f} mHa^2 | error_1 formula "
              f"{pe1:+.6f} verifier {(e1 - e1b) * 1e3:+.6f} mHa", flush=True)
    # ---- statistics and the kill test (a)
    nblk = (A + 15) // 16
    last = nblk - 1
    for b in bits:
        allp = []
        for c in cols:
            dE2 = price[(c, b)][0]
            live = ~dead[n1:, c]
            allp.append(dE2[live])
        allp = np.concatenate(allp)
        qs = np.quantile(allp, [0.001, 0.01, 0.05, 0.25, 0.5, 0.75])
        print(f"bit {b}: {allp.size} live two-body words of loader {k}: price (mHa^2) quantiles 0.1% {qs[0]:+.2e} 1% "
              f"{qs[1]:+.2e} 5% {qs[2]:+.2e} 25% {qs[3]:+.2e} 50% {qs[4]:+.2e} 75% {qs[5]:+.2e}; negative "
              f"{int((allp < 0).sum())} (sum {allp[allp < 0].sum():+.4f}); under 1e-4: {int((allp < 1e-4).sum())}, "
              f"under 1e-3: {int((allp < 1e-3).sum())}", flush=True)
        picks = []
        nunc = 0
        for j in range(J):
            cj = [J * k + j] + ([J * (k - 1) + j] if k > 0 else [])
            for h in range(last):
                rows = [16 * h + l for l in range(a.lmin, 16) if 16 * h + l < A]
                if any(dead[r, c] for r in rows for c in cj):
                    continue
                nunc += 1
                best = None
                for r in rows:
                    if r < n1:
                        continue
                    for c in cj:
                        p2 = price[(c, b)][0][r - n1]
                        if best is None or p2 < best[0]:
                            best = (p2, r, c)
                if best is not None:
                    picks.append(best)
        pr = np.array([p[0] for p in picks])
        print(f"bit {b}: uncovered (block, word) {nunc}, picks {len(picks)}; exact-formula price a pick: min "
              f"{pr.min():+.2e} median {np.median(pr):+.2e} mean {pr.mean():+.2e} max {pr.max():+.2e}; sum "
              f"{pr.sum():+.4f} mHa^2; negative picks {int((pr < 0).sum())}", flush=True)
        for name, sel in (("all", picks), ("half", [picks[i] for i in rng.permutation(len(picks))[:len(picks) // 2]])):
            w2 = words.copy()
            pe1 = 0.0
            for _, r, c in sel:
                w2[r, c] ^= (1 << b)
                pe1 += price[(c, b)][1][r - n1]
            e2, e1, lam = P.exact(w2)
            print(f"   bit {b} {name}: flips {len(sel)}; formula sum {sum(p[0] for p in sel):+.4f}; EXACT d(error_2^2) "
                  f"{(e2 - e2b) * 1e6:+.4f} mHa^2; error_1 {e1 * 1e3:.5f} mHa (base {e1b * 1e3:.5f}; first-order "
                  f"{(np.sqrt(max(P.e1sq + pe1, 0)) - np.sqrt(P.e1sq)) * 1e3:+.5f}); lambda_eff {lam:.6f}", flush=True)
            if a.out and name == "all":
                os.makedirs(a.out, exist_ok=True)
                np.savez(os.path.join(a.out, f"exact_k{k}_b{b}.npz"), angle_words=w2, tree_edges=edges)
    if a.out:
        os.makedirs(a.out, exist_ok=True)
        np.savez(os.path.join(a.out, f"prices_k{k}.npz"), cols=np.array(cols), bits=np.array(bits),
                 dE2=np.stack([np.stack([price[(c, b)][0] for c in cols]) for b in bits]),
                 dE1sq=np.stack([np.stack([price[(c, b)][1] for c in cols]) for b in bits]))


if __name__ == "__main__":
    main()
