"""chem_t10_relge (T10 seat 3, Nov 18; ledger row O14): R53-CHEM-GAUGE's CARRIER RULE inside R65's RelPlacerG.

THE RULE (R53-CHEM-GAUGE 'For the builders', exact through C2's own tool on fixed tables): C2's carrier gauge
(chem_ld2_e11car7.py) turns a signed-range word z into the stored word (-1)^n z + 64 h, with h = the slot monomials of the
row's low address (+ a correction at poles) and n fixed by h and the tree; so the bits 0..5 the loaders read are those of
y = eps x z, eps = (-1)^n. A relation above bit 0 asked of z breaks in the gauge (the P tables of Nov 17 lost their
carriers to keep them); asked of y it survives the gauge and every carrier stays exact.

HERE: every relation of the relation file is asked of the bits of y. eps is taken from the EXACT recursion
(tools/r53_chem_gauge_car.carrier_gauge) on the rows' words as they stand at each freeze (frozen words at their values,
free words at their nearest integers in [-32, 31]), for all 16 low addresses, because the address of a row is decided by
the placer's assignment. In the placer's programme a word with eps = -1 costs, for y-residue q, what its z-residue -q mod 4
costs (bit 0 is the same; bit 1 of an odd word flips); the residues handed back to the rounding are z-residues, so the
rounding itself is unchanged (nearest z of the residue). Bit-0 wires are untouched by the rule.
    fp = RelPlacerGE(R, B, rels, edges)            # rels as RelPlacerG's (chem_rk7_rel.parse / parse_g)
    r, m = fp.residues(t, c, frozen, value)        # z-residues, as RelPlacer's
    fp.broken_stored(words)                        # relations broken on C2's stored words (rows in INIT order)
    python chem_t10_relge.py --selftest
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chem_rk7_relplace import BIG, HI_W, LO_W  # noqa: E402  (read-only imports)
from r53_chem_gauge_car import MONO, carrier_gauge  # noqa: E402
from r65_chem_all15_relg import RelPlacerG, check_g, parse_g, fmt_g, P51M, _wrap  # noqa: E402

PERM = np.array([0, 3, 2, 1], np.int64)      # residue q mod 4 -> -q mod 4 (an involution)


def stored_words(z, edges, addrs, brot=7):
    """C2's carrier form of signed-range words z (rows at addresses addrs), by R53's exact recursion."""
    z = np.asarray(z, np.int64) % (1 << brot)
    low = np.asarray(addrs, np.int64) % 16
    T = z.shape[1]
    want = np.stack([((low & MONO[t % 5]) == MONO[t % 5]).astype(np.int64) for t in range(T)], 1)
    st, _n, _h = carrier_gauge(z, edges, want, brot)
    return st


class RelPlacerGE(RelPlacerG):
    def __init__(self, R, B, rels, edges, n_one=76, ncols=75, fixed=False, gfit=True, rounds=1, pin_blocks=(4,),
                 mode="exact", brot=7):
        """mode: 'exact' (eps from the rows' words at every freeze), 'blind' (eps = the pole-blind table N c(l)),
        'off' (eps = +1: RelPlacerG exactly, the twin)."""
        super().__init__(R, B, rels, n_one=n_one, ncols=ncols, fixed=fixed, gfit=gfit, rounds=rounds,
                         pin_blocks=pin_blocks)
        self.edges = np.asarray(edges, np.int64)
        self.mode, self.brot = mode, brot
        self.lowp = self.paddr & 15                                      # low address of every place
        self.want16 = np.array([[int((l & MONO[t % 5]) == MONO[t % 5]) for t in range(ncols)] for l in range(16)],
                               np.int64)
        self.E = np.ones((R * B, 16, ncols), np.int8)
        self.eps_changes = 0          # frozen relation words whose eps changed after their freeze (poles)
        self._Eprev = None
        self.assign_every = int(os.environ.get("T10_ASSIGN_EVERY", "1"))   # the rows' reassignment every k-th freeze
        self.ncalls = 0

    def step(self, wc):
        """RelPlacerG.step; with assign_every k > 1 the reassignment of rows to places (the programme over every place
        of a group) runs on every k-th call only; the block constants are fitted on every call."""
        self.ncalls += 1
        if self.assign_every > 1 and self.ncalls > 1 and self.ncalls % self.assign_every \
                and not getattr(self, "infeasible", 0):     # never skipped on the first call or while a row is inconsistent
            if self.gfit:
                self.fit_g(wc)
            return
        return super().step(wc)

    # ---- eps of every (row, low address, column)
    def update_eps(self, t, frozen, value):
        if self.mode == "off":
            return
        n_rows = t.shape[0]
        if self.mode == "blind":
            z = np.zeros((n_rows, self.ncols), np.int64)
        else:
            z = np.where(frozen, value, np.clip(np.round(t), LO_W, HI_W)).astype(np.int64)
        zm = z % (1 << self.brot)
        for l in range(16):
            want = np.broadcast_to(self.want16[l], zm.shape)
            _st, n, _h = carrier_gauge(zm, self.edges, want, self.brot)
            self.E[:, l, :] = (1 - 2 * n).astype(np.int8)
        if self._Eprev is not None:
            rel = self.mod > 2
            ar = np.arange(n_rows)
            lr = self.lowp[self.place]
            ch = (self.E[ar, lr] != self._Eprev[ar, lr]) & frozen & rel[None, :] & (np.mod(value, 2) == 1)
            self.eps_changes += int(ch.sum())
        self._Eprev = self.E.copy()

    # ---- the batch cost in y (per row at one low address each, or at all 16)
    def _batchcost_e(self, wc, rows, b, L=None):
        nv = len(self.V[b])
        nc = 1 << nv
        n = len(rows)
        tot = np.zeros((n, nc) if L is not None else (n, 16, nc))
        for j in range(5):
            col = 5 * b + j
            m = int(self.mod[col])
            if m <= 1:
                continue
            cst = np.minimum(wc[col][rows], BIG)                              # (n, 4): cost by z-residue
            c0, c1 = self.need[b][j]
            if m == 2:
                add = cst[:, c0]
                tot += add if L is not None else add[:, None, :]
                continue
            if L is not None:
                neg = self.E[rows, L, col] < 0
                cy = np.where(neg[:, None], cst[:, PERM], cst)                # (n, 4): cost by y-residue
            else:
                neg = self.E[rows, :, col] < 0
                cy = np.where(neg[:, :, None], cst[:, None, PERM], cst[:, None, :])     # (n, 16, 4)
            k0 = np.zeros(nc, np.int64) if c0 is None else c0 + 1
            k1 = np.zeros(nc, np.int64) if c1 is None else c1 + 1
            key = k0 * 3 + k1
            uk, inv = np.unique(key, return_inverse=True)
            mins = []
            for kv in uk:
                a0, a1 = kv // 3 - 1, kv % 3 - 1
                qs = [q for q in range(4) if (a0 < 0 or (q & 1) == a0) and (a1 < 0 or ((q >> 1) & 1) == a1)]
                mins.append(cy[..., qs].min(axis=-1))
            M = np.stack(mins, axis=-1)
            tot += M[..., inv]
        return np.minimum(tot, BIG)

    def _batchcost(self, wc, rows, b):         # (not used by the overridden paths; kept exact for rows at their places)
        return self._batchcost_e(wc, rows, b, self.lowp[self.place[rows]])

    def _messages(self, wc):
        rows = np.arange(self.R * self.B)
        L = self.lowp[self.place]
        return {b: self._pair_min(self._batchcost_e(wc, rows, b, L), b) for b in self.batches}

    def _dp(self, wc, rows, places, diag=False):
        """the programme, vectorised over the codes of one out-syndrome class (a loop of at most 2^r classes, not 2^nv
        codes): the same minima and the same argmins (the smallest code among ties) as _dp_loop."""
        I = len(rows)
        P = 1 if diag else len(places)
        lo, hi = self.lo, self.hi
        lP = self.lowp[places]
        so = self.syn_out[lo]
        ns = self._ns(lo + 1)
        arg = None
        if diag:
            cb = self._batchcost_e(wc, rows, lo, lP)
            f = np.full((I, ns), BIG)
            a0 = np.zeros((I, ns), np.int64)
            for s_ in range(ns):
                xs = np.flatnonzero(so == s_)
                if len(xs):
                    k_ = np.argmin(cb[:, xs], axis=1)
                    f[:, s_] = cb[np.arange(I), xs[k_]]
                    a0[:, s_] = xs[k_]
            arg = {lo: a0}
            f = f[:, None, :]
        else:
            cb16 = self._batchcost_e(wc, rows, lo, None)
            f16 = np.full((I, 16, ns), BIG)
            for s_ in range(ns):
                xs = np.flatnonzero(so == s_)
                if len(xs):
                    f16[:, :, s_] = cb16[:, :, xs].min(axis=2)
            f = f16[:, lP, :]
        for b in range(lo + 1, hi + 1):
            si = self.syn_in[b]
            so = self.syn_out[b] if b < hi else np.zeros(len(si), np.int64)
            ns2 = self._ns(b + 1)
            od = np.argsort(so, kind="stable")                                     # codes by out-syndrome, ascending inside
            so_s = so[od]
            st_ = np.flatnonzero(np.r_[True, so_s[1:] != so_s[:-1]])
            cls = so_s[st_]
            if b in self.T:
                Tb = self.T[b][places]                                             # (P,) or, diag, (I,)
            else:
                Tb = np.zeros(I if diag else P, np.int64)
            g = np.full((I, P, ns2), BIG)
            if diag:
                cb = self._batchcost_e(wc, rows, b, lP)[:, od]                     # (I, nc) sorted codes
                idx = si[od][None, :] ^ Tb[:, None]                                # (I, nc)
                F = np.take_along_axis(f[:, 0, :], idx, axis=1) + cb
                ag = np.zeros((I, P, ns2), np.int64)
                ends = np.r_[st_[1:], len(od)]
                for s2, a_, e_ in zip(cls, st_, ends):
                    k_ = np.argmin(F[:, a_:e_], axis=1)
                    g[:, 0, s2] = F[np.arange(I), a_ + k_]
                    ag[:, 0, s2] = od[a_ + k_]
            else:
                cb16 = self._batchcost_e(wc, rows, b, None)[:, :, od]              # (I, 16, nc)
                idx = si[od][None, :] ^ Tb[:, None]                                # (P, nc)
                F = np.take_along_axis(f, idx[None, :, :], axis=2)                 # (I, P, nc)
                F += cb16[:, lP, :]
                g[:, :, cls] = np.minimum.reduceat(F, st_, axis=2)
                del F, cb16
            f = np.minimum(g, BIG)
            if diag:
                arg[b] = ag[:, 0, :]
        if diag:
            return f[:, 0, 0], arg
        return f[:, :, 0], None

    def _dp_loop(self, wc, rows, places, diag=False):
        I = len(rows)
        P = 1 if diag else len(places)
        lo, hi = self.lo, self.hi
        lP = self.lowp[places]
        so = self.syn_out[lo]
        ns = self._ns(lo + 1)
        arg = None
        if diag:
            cb = self._batchcost_e(wc, rows, lo, lP)
            f = np.full((I, ns), BIG)
            a0 = np.zeros((I, ns), np.int64)
            for s_ in range(ns):
                xs = np.flatnonzero(so == s_)
                if len(xs):
                    k_ = np.argmin(cb[:, xs], axis=1)
                    f[:, s_] = cb[np.arange(I), xs[k_]]
                    a0[:, s_] = xs[k_]
            arg = {lo: a0}
            f = f[:, None, :]
        else:
            cb16 = self._batchcost_e(wc, rows, lo, None)
            f16 = np.full((I, 16, ns), BIG)
            for s_ in range(ns):
                xs = np.flatnonzero(so == s_)
                if len(xs):
                    f16[:, :, s_] = cb16[:, :, xs].min(axis=2)
            f = f16[:, lP, :]
        ar = np.arange(P)
        ii = np.arange(I)[:, None]
        for b in range(lo + 1, hi + 1):
            if diag:
                cb = self._batchcost_e(wc, rows, b, lP)                      # (I, nc)
                nc = cb.shape[1]
            else:
                cbP = np.ascontiguousarray(self._batchcost_e(wc, rows, b, None).transpose(2, 0, 1))[:, :, lP]   # (nc, I, P)
                nc = cbP.shape[0]
            si = self.syn_in[b]
            so = self.syn_out[b] if b < hi else np.zeros(len(si), np.int64)
            ns2 = self._ns(b + 1)
            if b in self.T:
                Tb = (self.T[b][places] if diag else self.T[b][places][None, :]).reshape(I if diag else 1, P if not diag else 1)
            else:
                Tb = np.zeros((1, 1), np.int64)
            g = np.full((I, P, ns2), BIG)
            ag = np.zeros((I, P, ns2), np.int64) if diag else None
            for x in range(nc):
                idx = np.broadcast_to(si[x] ^ Tb, (I, P))
                cand = f[ii, ar[None, :], idx] + (cb[:, x][:, None] if diag else cbP[x])
                s2 = so[x]
                upd = cand < g[:, :, s2]
                g[:, :, s2] = np.where(upd, cand, g[:, :, s2])
                if diag:
                    ag[:, :, s2] = np.where(upd, x, ag[:, :, s2])
            f = np.minimum(g, BIG)
            if diag:
                arg[b] = ag[:, 0, :]
        if diag:
            return f[:, 0, 0], arg
        return f[:, :, 0], None

    def residues(self, t, c, frozen, value):
        t = np.asarray(t, np.float64)
        c = np.asarray(c, np.float64)
        value = np.asarray(value, np.int64)
        self.update_eps(t, frozen, value)
        wc = self._wordcost(t, c, frozen, value)
        self.step(wc)
        r = np.zeros(t.shape, np.int64)
        self.infeasible = 0
        for g in range(self.R):
            rows = np.flatnonzero(self.group == g)
            pls = self.place[rows]
            cost, arg = self._dp(wc, rows, pls, diag=True)
            self.infeasible += int(np.sum(cost >= 1e17))
            codes = self._backtrack(arg, pls)
            L = self.lowp[pls]
            for b in self.batches:
                x = codes[b]
                for j in range(5):
                    col = 5 * b + j
                    m = int(self.mod[col])
                    if m <= 1:
                        continue
                    c0, c1 = self.need[b][j]
                    b0 = c0[x] if c0 is not None else None
                    b1 = c1[x] if c1 is not None else None
                    if m == 2:
                        r[rows, col] = b0
                        continue
                    cst = wc[col][rows]
                    neg = self.E[rows, L, col] < 0
                    cy = np.where(neg[:, None], cst[:, PERM], cst)
                    best_r = np.zeros(len(rows), np.int64)
                    best_v = np.full(len(rows), np.inf)
                    for q in range(4):
                        ok = np.ones(len(rows), bool)
                        if b0 is not None:
                            ok &= (q & 1) == b0
                        if b1 is not None:
                            ok &= ((q >> 1) & 1) == b1
                        v = np.where(ok, cy[:, q], np.inf)
                        sel = v < best_v
                        best_r[sel], best_v[sel] = q, v[sel]
                    r[rows, col] = np.where(neg, PERM[best_r], best_r)       # back to the z-residue
        return r, self.mod.copy()

    def broken_stored(self, words):
        """words[row, column] (signed range, rows in INIT order): relations broken on C2's stored words at the rows'
        addresses (with the fitted constants)."""
        od = self.order(np.arange(self.R * self.B))
        st = stored_words(np.asarray(words, np.int64)[od], self.edges, self.paddr, self.brot)
        return check_g(st, self.paddr, self.rels_used())


def _tree(rng, n=76):
    return np.array([(int(rng.integers(0, t + 1)), t + 1) for t in range(n - 1)], np.int64)


def selftest(seed=0, big=True):
    rng = np.random.default_rng(seed)
    edges = _tree(rng)
    # (1) mode 'off' is RelPlacerG exactly (same residues and places on the same inputs)
    R, B = 2, 30
    small = parse_g("13 | 0:0 1:0 | F:7;14 | 1:0 2:0 | F:9;14 | 0:1 1:1 2:1 | F:37;14 | 3:0 4:0 | 3;13 | 2:1 3:1 | F:5")
    t = rng.uniform(-30, 30, (R * B, 75))
    c = rng.uniform(0.5, 2.0, (R * B, 75))
    fz = (rng.random((R * B, 75)) < 0.3) & np.isin(np.arange(75) % 5, (0, 2, 4))[None, :]
    z = np.round(t).astype(np.int64)
    a = RelPlacerG(R, B, small, pin_blocks=())
    e = RelPlacerGE(R, B, small, edges, pin_blocks=(), mode="off")
    ra, _ = a.residues(t, c, fz, z)
    re_, _ = e.residues(t, c, fz, z)
    assert np.array_equal(ra, re_) and np.array_equal(a.place, e.place) and np.array_equal(a.G, e.G)
    print(f"selftest {seed} (1): mode off = RelPlacerG (residues, places and constants equal)")
    # (2) the programme's least cost with eps equals brute force over y-patterns for a row at a place
    small = parse_g("13 | 0:0 1:0 | F:7;14 | 1:0 2:0 | F:9;14 | 0:1 1:1 | F:10;14 | 3:0 4:0 | 3")
    e = RelPlacerGE(R, B, small, edges, pin_blocks=(), mode="exact")
    e.update_eps(t, np.zeros(t.shape, bool), z)
    assert (e.E < 0).any() and (e.E > 0).any()
    wc = e._wordcost(t, c, np.zeros(t.shape, bool), z)
    M, _ = e._dp(wc, np.arange(B), np.arange(B))
    import itertools
    cols = list(e.pl.cc)
    mods = [int(e.mod[cl]) for cl in cols]
    for (row, pl) in ((0, 0), (3, 7), (5, 17)):
        best = BIG
        l = int(e.lowp[pl])
        for combo in itertools.product(*[range(m) for m in mods]):           # combo = y-residues
            w = np.zeros((1, 75), np.int64)
            cost = 0.0
            for cl, q in zip(cols, combo):
                w[0, cl] = q
                qz = q if (int(e.mod[cl]) == 2 or e.E[row, l, cl] > 0) else int(PERM[q])
                cost += wc[cl][row, qz]
            if check_g(w, np.array([e.paddr[pl]]), e.rels_used()) == 0:
                best = min(best, cost)
        assert abs(M[row, pl] - best) < 1e-9 * max(1.0, best), (row, pl, M[row, pl], best)
    print(f"selftest {seed} (2): the programme with eps = brute force over y-patterns (3 rows x places, {len(cols)} columns)")
    if not big:
        return
    # (3) a mock staged freeze, P51m's relations (45 wires + 6 bit-1), 3 groups of 127: relations on C2's STORED words
    R, B, T_ = 3, 127, 75
    rels = parse_g(P51M)
    for mode in ("off", "exact"):
        r_ = np.random.default_rng(seed + 1)
        x = r_.uniform(-31, 31, (R * B, T_))
        c = np.exp(r_.normal(0, 2.5, (R * B, T_))) * np.exp(r_.normal(0, 1.5, (R * B, 1)))
        fp = RelPlacerGE(R, B, rels, edges, mode=mode)
        frozen = np.zeros(x.shape, bool)
        val = np.zeros(x.shape, np.int64)
        rank = np.argsort(-c, axis=None)
        cost = 0.0
        for st in np.array_split(rank, 12):
            t = np.where(frozen, val, x + (0.0 if not frozen.any() else 0.05 * r_.normal(size=x.shape)))
            r, m = fp.residues(t, c, frozen, val)
            assert fp.infeasible == 0
            allowed = _wrap(r + m * np.round((t - r) / m), m).astype(np.int64)
            i, j = np.unravel_index(st, x.shape)
            val[i, j] = allowed[i, j]
            cost += float((c[i, j] * ((allowed[i, j] - t[i, j]) ** 2 - (np.round(t[i, j]) - t[i, j]) ** 2)).sum())
            frozen[i, j] = True
        bz = fp.broken(val)
        bs = fp.broken_stored(val)
        print(f"selftest {seed} (3) [{mode}]: relations broken on z {bz}, on C2's stored words {bs}; excess cost {cost:.1f}; "
              f"late eps changes on frozen odd relation words {fp.eps_changes}", flush=True)
        if mode == "off":
            assert bz == 0
        else:
            assert bs <= fp.eps_changes, (bs, fp.eps_changes)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--small", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest(big=not a.small)
        print("selftest passed")
