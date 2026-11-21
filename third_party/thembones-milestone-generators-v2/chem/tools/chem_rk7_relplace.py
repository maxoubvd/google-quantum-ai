"""chem_rk7_relplace (RK7 holder 2, Nov 17): R50's FreezePlacer GENERALISED TO RELATION FILES (any GF(2) relation
among the delta outputs of a loader, any bit, F-table / mask targets: the format of tools/chem_rk7_rel.py).

WR7 measured R50's FreezePlacer (rows re-placed inside their groups by an assignment before every freeze) 4.1 mHa^2
under DynChain on the same eight wires [M, LEAD_NOTES 06:08Z]. The Placer holds single bit-0 wires only; this class
holds every relation of a relation file EXACTLY:
  * per word of a relation column a cost per residue r mod m (m = 2 for bit-0 columns, 4 when a relation reads bit 1):
    c (z_r - t)^2, z_r the nearest word = r mod m in the signed range [-32, 31]; a frozen word: 0 at its own residue;
  * per row and per place of its group the EXACT least cost of a relation-consistent pattern at that place's targets:
    a dynamic programme over the batches (batch b = columns 5b..5b+4; loader k reads batches k-1, k; the state after
    batch b is loader b+1's syndrome of batch b's bits), vectorised over rows and places;
  * one assignment per group (scipy linear_sum_assignment) on that cost matrix, kept if it lowers the total;
  * the residues of every relation word at the rows' places (the programme's argmin, backtracked).
Same interface as r50_chem_lbl_hook.FreezePlacer: .place, residues(t, c, frozen, value) -> (r, m), order(order_in),
.pl.cc (the decided columns), addr().
    python chem_rk7_relplace.py --selftest
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

BIG = 1e30
LO_W, HI_W = -32, 31          # the signed range of a 7-bit word


class _PL:
    def __init__(self, cc):
        self.cc = cc


class RelPlacer:
    def __init__(self, R, B, rels, n_one=76, ncols=75, fixed=False):
        from chem_rk7_rel import ttable
        self.R, self.B, self.n_one, self.ncols, self.fixed = R, B, n_one, ncols, fixed
        self.paddr = n_one + np.arange(R * B)
        self.place = np.arange(R * B)
        self.group = np.arange(R * B) // B
        self.moves = 0
        self.rels = rels
        # variables of batch b: the (slot, bit) pairs read by loader b (as its second batch) or loader b + 1 (first)
        vb = {}
        for (k, outs, _t) in rels:
            for (j, bit) in outs:
                vb.setdefault(k - 1, set()).add((j, bit))
                vb.setdefault(k, set()).add((j, bit))
        self.batches = sorted(vb)
        self.V = {b: sorted(vb[b]) for b in self.batches}
        mod = np.ones(ncols, np.int64)
        for b in self.batches:
            for (j, bit) in self.V[b]:
                mod[5 * b + j] = max(mod[5 * b + j], 1 << (bit + 1))
        self.mod = mod
        self.pl = _PL(np.flatnonzero(mod > 1))
        lo, hi = self.batches[0], self.batches[-1]
        assert self.batches == list(range(lo, hi + 1)), "relation batches must be contiguous"
        self.lo, self.hi = lo, hi
        # loader k's relations and their syndrome tables on batch k-1 (out) and batch k (in), per code x of the batch
        self.rk = {k: [r for r in rels if r[0] == k] for k in range(lo + 1, hi + 1)}
        assert not any(r[0] == lo for r in rels), "loader lo has relations: its first batch would carry variables"
        self.syn_out, self.syn_in = {}, {}
        for b in self.batches:
            nv = len(self.V[b])
            codes = np.arange(1 << nv)
            bits = (codes[:, None] >> np.arange(nv)[None, :]) & 1          # (ncode, nv)
            pos = {v: i for i, v in enumerate(self.V[b])}

            def syn(rel_list):
                s = np.zeros(1 << nv, np.int64)
                for ri, (_k, outs, _t) in enumerate(rel_list):
                    par = np.zeros(1 << nv, np.int64)
                    for (j, bit) in outs:
                        par ^= bits[:, pos[(j, bit)]]
                    s |= par << ri
                return s
            self.syn_in[b] = syn(self.rk.get(b, []))           # loader b reads batch b as its second batch
            self.syn_out[b] = syn(self.rk.get(b + 1, []))      # loader b + 1 reads batch b as its first batch
        # per slot j of batch b: for every code x the residue bits (bit0, bit1) it fixes and which bits it fixes
        self.need = {}
        for b in self.batches:
            pos = {v: i for i, v in enumerate(self.V[b])}
            nv = len(self.V[b])
            codes = np.arange(1 << nv)
            per = []
            for j in range(5):
                c0 = ((codes >> pos[(j, 0)]) & 1) if (j, 0) in pos else None
                c1 = ((codes >> pos[(j, 1)]) & 1) if (j, 1) in pos else None
                per.append((c0, c1))
            self.need[b] = per
        # targets: T[k][place] = the bit vector of loader k's relation targets at that place's address
        self.T = {}
        for k, rl in self.rk.items():
            t = np.zeros(R * B, np.int64)
            for ri, (_k, _o, tg) in enumerate(rl):
                t |= ttable(tg, self.paddr) << ri
            self.T[k] = t

    def addr(self):
        return self.paddr[self.place]

    # ---- costs
    def _wordcost(self, t, c, frozen, value):
        """cost[row, col, r] for r < 4 (r < m used) of every relation column."""
        n = t.shape[0]
        out = {}
        for col in self.pl.cc:
            m = int(self.mod[col])
            tc, cc_ = t[:, col], c[:, col]
            cst = np.full((n, 4), BIG)
            for r in range(m):
                z = r + m * np.round((tc - r) / m)
                z = np.where(z > HI_W, z - m, np.where(z < LO_W, z + m, z))
                cst[:, r] = cc_ * (z - tc) ** 2
            fz = frozen[:, col]
            if fz.any():
                vr = np.mod(value[fz, col], m)
                cst[fz] = BIG
                cst[np.flatnonzero(fz), vr] = 0.0
            out[col] = cst
        return out

    def _batchcost(self, wc, rows, b):
        """cost_b[row, code] of batch b's variables for the given rows (the cheapest residue a code allows)."""
        nv = len(self.V[b])
        tot = np.zeros((len(rows), 1 << nv))
        for j in range(5):
            col = 5 * b + j
            if self.mod[col] <= 1:
                continue
            cst = wc[col][rows]                                  # (n, 4)
            m = int(self.mod[col])
            c0, c1 = self.need[b][j]
            if m == 2:
                assert c0 is not None
                tot += cst[:, c0] if True else 0
            else:                                                # m == 4: residue = bit0 + 2 bit1
                opts = []
                for b0 in (0, 1):
                    for b1 in (0, 1):
                        opts.append((b0, b1, cst[:, b0 + 2 * b1]))
                best = np.full((len(rows), 1 << nv), BIG)
                for (b0, b1, v) in opts:
                    ok = np.ones(1 << nv, bool)
                    if c0 is not None:
                        ok &= c0 == b0
                    if c1 is not None:
                        ok &= c1 == b1
                    best = np.where(ok[None, :], np.minimum(best, v[:, None]), best)
                tot += best
        return np.minimum(tot, BIG)

    # ---- the programme
    def _ns(self, k):
        """states after the batch that loader k reads first: 2^(number of loader k's relations)."""
        return 1 << len(self.rk.get(k, [])) if k <= self.hi else 1

    def _dp(self, wc, rows, places, diag=False):
        """diag=False: least cost [row, place] of a relation-consistent pattern for every row of `rows` at every place of
        `places` (one group). diag=True: row i only at places[i]; returns (cost [row], argmins per batch)."""
        I = len(rows)
        P = 1 if diag else len(places)
        lo, hi = self.lo, self.hi
        cb = self._batchcost(wc, rows, lo)
        so = self.syn_out[lo]
        ns = self._ns(lo + 1)
        f = np.full((I, ns), BIG)
        a0 = np.zeros((I, ns), np.int64)
        for s_ in range(ns):
            xs = np.flatnonzero(so == s_)
            if len(xs):
                k_ = np.argmin(cb[:, xs], axis=1)
                f[:, s_] = cb[np.arange(I), xs[k_]]
                a0[:, s_] = xs[k_]
        arg = {lo: a0}
        f = np.repeat(f[:, None, :], P, axis=1)                  # (I, P, ns)
        ar = np.arange(P)
        for b in range(lo + 1, hi + 1):
            cb = self._batchcost(wc, rows, b)
            si = self.syn_in[b]
            so = self.syn_out[b] if b < hi else np.zeros(len(si), np.int64)
            ns2 = self._ns(b + 1)
            if b in self.T:
                Tb = (self.T[b][places] if diag else self.T[b][places][None, :]).reshape(I if diag else 1, P if not diag else 1)
            else:
                Tb = np.zeros((1, 1), np.int64)
            g = np.full((I, P, ns2), BIG)
            ag = np.zeros((I, P, ns2), np.int64) if diag else None
            ii = np.arange(I)[:, None]
            for x in range(cb.shape[1]):
                idx = si[x] ^ Tb                                   # the state batch b-1 must have left
                idx = np.broadcast_to(idx, (I, P))
                cand = f[ii, ar[None, :], idx] + cb[:, x][:, None]  # (I, P)
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

    def _backtrack(self, arg, places):
        """codes x_b per row, each row at places[row] (the diag programme's argmins)."""
        lo, hi = self.lo, self.hi
        I = len(places)
        ri = np.arange(I)
        out = {}
        s = np.zeros(I, np.int64)
        for b in range(hi, lo, -1):
            x = arg[b][ri, s]
            out[b] = x
            Tb = self.T[b][places] if b in self.T else 0
            s = self.syn_in[b][x] ^ Tb
        out[lo] = arg[lo][ri, s]
        return out

    def step(self, wc):
        from scipy.optimize import linear_sum_assignment
        for g in range(self.R):
            rows = np.flatnonzero(self.group == g)
            places = g * self.B + np.arange(self.B)
            M, _ = self._dp(wc, rows, places)
            cur = self.place[rows] - g * self.B
            cur_tot = float(M[np.arange(len(rows)), cur].sum())
            if self.fixed:
                continue
            ri, ci = linear_sum_assignment(np.minimum(M, 1e18))
            tot = float(M[ri, ci].sum())
            if tot < cur_tot - 1e-12 * abs(cur_tot) and tot < 1e17:
                new = self.place.copy()
                new[rows[ri]] = places[ci]
                self.moves += int((new != self.place).sum())
                self.place = new

    def residues(self, t, c, frozen, value):
        t = np.asarray(t, np.float64)
        c = np.asarray(c, np.float64)
        value = np.asarray(value, np.int64)
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
                    else:
                        cst = wc[col][rows]
                        best_r = np.zeros(len(rows), np.int64)
                        best_v = np.full(len(rows), np.inf)
                        for q in range(4):
                            ok = np.ones(len(rows), bool)
                            if b0 is not None:
                                ok &= (q & 1) == b0
                            if b1 is not None:
                                ok &= ((q >> 1) & 1) == b1
                            v = np.where(ok, cst[:, q], np.inf)
                            sel = v < best_v
                            best_r[sel], best_v[sel] = q, v[sel]
                        r[rows, col] = best_r
        return r, self.mod.copy()

    def order(self, order_in):
        out = np.empty_like(order_in)
        out[self.place] = order_in
        return out


def selftest(seed=0):
    """(1) the programme's least cost equals brute force over every pattern (a small relation set); (2) residues over
    three freeze phases keep every relation at the rows' final places (a richer set, F tables and masks)."""
    import itertools
    from chem_rk7_rel import check, parse
    rng = np.random.default_rng(seed)
    R, B = 2, 6
    n = R * B
    small = parse("13 | 0:0 1:0 | F:7;14 | 1:0 2:0 | F:9;14 | 0:1 1:1 | F:10;14 | 3:0 4:0 | 3")
    fp = RelPlacer(R, B, small)
    t = rng.uniform(-20, 20, (n, 75))
    c = rng.uniform(0.5, 2.0, (n, 75))
    wc = fp._wordcost(t, c, np.zeros((n, 75), bool), np.round(t).astype(np.int64))
    M, _ = fp._dp(wc, np.arange(B), np.arange(B))
    cols = list(fp.pl.cc)
    mods = [int(fp.mod[cl]) for cl in cols]
    for (row, pl) in ((0, 0), (2, 5)):
        best = BIG
        for combo in itertools.product(*[range(m) for m in mods]):
            w = np.zeros((1, 75), np.int64)
            for cl, q in zip(cols, combo):
                w[0, cl] = q
            if check(w, np.array([fp.paddr[pl]]), small) == 0:
                best = min(best, sum(wc[cl][row, q] for cl, q in zip(cols, combo)))
        print(f"selftest {seed}: row {row} place {pl}: programme {M[row, pl]:.6f}  brute force {best:.6f}  ({len(cols)} columns)")
        assert abs(M[row, pl] - best) < 1e-9 * max(1.0, best)
    rich = parse("12 | 0:0 1:0 | F:5;12 | 2:0 3:0 | F:6;13 | 0:0 1:0 | F:7;13 | 0:1 1:1 2:1 3:1 4:1 | F:8;"
                 "14 | 1:0 2:0 | F:9;14 | 0:1 1:1 | F:10;14 | 3:0 4:0 | 3")
    fp = RelPlacer(R, B, rich)
    frozen = np.zeros((n, 75), bool)
    z = np.round(t).astype(np.int64)
    for phase in range(4):
        r, m = fp.residues(t, c, frozen, z)
        assert fp.infeasible == 0
        for col in fp.pl.cc:
            mm = int(m[col])
            zz = r[:, col] + mm * np.round((t[:, col] - r[:, col]) / mm)
            zz = np.where(zz > HI_W, zz - mm, np.where(zz < LO_W, zz + mm, zz))
            z[:, col] = np.where(frozen[:, col], z[:, col], zz)
        frozen = frozen | (rng.random((n, 75)) < (0.4 if phase < 3 else 1.0))
        t = np.where(frozen, z, t + rng.normal(0, 0.5, t.shape))
    perm = np.argsort(fp.place)
    bad = check(np.mod(z[perm], 128), fp.paddr, rich)
    print(f"selftest {seed}: four freeze phases: relations broken {bad}; rows off their first place "
          f"{int((fp.place != np.arange(n)).sum())}; words outside [-32, 31] {int(np.sum((z < LO_W) | (z > HI_W)))}")
    assert bad == 0


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        for s in range(3):
            selftest(s)
        print("selftest passed")
