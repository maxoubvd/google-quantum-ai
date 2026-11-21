"""r65_chem_all15_relg (R65-CHEM-ALL15, second hold): THE BLOCK CONSTANTS OF EVERY RELATION FITTED AT FREEZE TIME, inside
RK7's RelPlacer (tools/chem_rk7_relplace.py: the placer of the relation runs `chem_rk7_frel.py` / `chem_t10_frel.py
--fp --relp @rel/<file>`). A subclass: same calls, same state, plus one array.

  A relation (k | outs | target) asks, in every row,  XOR over its outputs = target(address).  The loader takes any
  target of the family F = f(address & 15) XOR g(address >> 4) for nothing (R52 / R53: the plain compiler reads the
  same count with masks, with random F tables and with f, g of any kind; the relation files already carry `F:seed`
  targets). HERE g IS FITTED: the target becomes  target(address) XOR g_rel(address >> 4),  one free bit per
  (relation, block of 16 addresses), decided together with the rows' places before every freeze:
    (1) one backward and one forward pass of the placer's own programme over all rows at their places give, for every
        loader k and every row, the row's exact least cost as a function of loader k's whole target vector (an XOR
        min-plus product, 2^r x 2^r for r relations on the loader);
    (2) per block the constant vector of loader k that minimises the sum over the block's rows is taken (the exact
        joint optimum over the 2^r vectors of that loader; loaders in turn, forward messages refreshed: Gauss-Seidel);
    (3) RK7's assignment step (rows to places inside their group) under the new targets.
  A frozen word of the wrong residue costs BIG in the placer, so a constant that rows have frozen under never flips.

    fp = RelPlacerG(R, B, rels)                       # rels: chem_rk7_rel.parse(...) or parse_g(...) (with G/ targets)
    r, m = fp.residues(t, c, frozen, value)           # as RelPlacer's
    fp.rels_used()                                    # the relations with their fitted targets: `G/<base>/<hex>`
    fp.get_state() / fp.set_state(place, G)           # a resumed run keeps fp.place AND fp.G
    fp.broken(words)                                  # words[row, column], rows in INIT order

TARGET KIND `G/<base>/<hex>` (for chem_rk7_rel.ttable, four lines; here: ttable_g): base = a mask or F:seed / B:x as
before; hex = the block constants, bit i = the constant of block i (i = address >> 4, 0..76).

    python r65_chem_all15_relg.py --selftest
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chem_rk7_rel import parse, ttable  # noqa: E402
from chem_rk7_relplace import BIG, HI_W, LO_W, RelPlacer  # noqa: E402

NBLK = 128


def split_target(t, low=False):
    """target -> (base target for chem_rk7_rel.ttable, the block constants as a python int[, the low table as an int]).
    G/<base>/<hex>/<hex16>: the optional fourth field is a table over address & 15 (RelPlacerGF: the model's f part)."""
    base, gv, fv = t, 0, 0
    if isinstance(t, str) and t.startswith("G/"):
        parts = t.split("/")
        base = parts[1] if parts[1][:2] in ("F:", "B:") else int(parts[1], 0)
        gv = int(parts[2], 16)
        fv = int(parts[3], 16) if len(parts) > 3 else 0
    return (base, gv, fv) if low else (base, gv)


def ttable_g(t, addrs):
    """chem_rk7_rel.ttable with the target kind G/<base>/<hex>."""
    addrs = np.asarray(addrs, np.int64)
    base, gv, fv = split_target(t, low=True)
    gb = np.array([(gv >> i) & 1 for i in range(NBLK)], np.int64)
    fb = np.array([(fv >> i) & 1 for i in range(16)], np.int64)
    return ttable(base, addrs) ^ gb[addrs >> 4] ^ fb[addrs & 15]


def parse_g(s):
    """A relation file with or without G/ targets -> [(k, outs, target)] (chem_rk7_rel.parse refuses G/ targets)."""
    if s.startswith("@"):
        s = open(s[1:]).read()
    held, lines = {}, []
    for ln in s.splitlines():
        ln = ln.split("#")[0]
        for part in ln.split(";"):
            part = part.strip()
            if not part:
                continue
            f = [v.strip() for v in part.split("|")]
            if len(f) == 3 and f[2].startswith("G/"):
                base = f[2].split("/")[1]
                held[len(lines)] = f[2]
                part = f"{f[0]} | {f[1]} | {base}"
            lines.append(part)
    rels = parse("\n".join(lines))
    assert len(rels) == len(lines)
    return [(k, o, held.get(i, t)) for i, (k, o, t) in enumerate(rels)]


def fmt_g(rels):
    return "\n".join(f"{k} | {' '.join(f'{j}:{b}' for (j, b) in o)} | {t}" for (k, o, t) in rels) + "\n"


def check_g(words, addrs, rels, per=False):
    """chem_rk7_rel.check with G/ targets: (relation, row) pairs broken on a table (rows in address order)."""
    w = np.asarray(words, np.int64)
    addrs = np.asarray(addrs, np.int64)
    bad = []
    for (k, o, t) in rels:
        acc = np.zeros(len(w), np.int64)
        for (j, b) in o:
            acc ^= ((w[:, 5 * k + j] ^ w[:, 5 * k - 5 + j]) >> b) & 1
        bad.append(int(np.sum(acc != ttable_g(t, addrs))))
    return np.array(bad) if per else int(sum(bad))


class RelPlacerG(RelPlacer):
    def __init__(self, R, B, rels, n_one=76, ncols=75, fixed=False, gfit=True, rounds=1, pin_blocks=(4,)):
        """rounds: alternations (constants, assignment) per freeze; the constants are fitted once more after the last
        assignment. pin_blocks: blocks (address >> 4) whose constants stay as given (default block 4 = addresses 64..79,
        shared with the one-body rows 64..75: the one-body tools then read the relation files unchanged)."""
        sp = [split_target(t) for (_k, _o, t) in rels]
        super().__init__(R, B, [(k, o, b) for (k, o, _t), (b, _g) in zip(rels, sp)], n_one=n_one, ncols=ncols, fixed=fixed)
        self.gfit, self.rounds = gfit, int(rounds)
        self.blk = self.paddr >> 4                                  # of a place
        self.nblk = int(self.blk.max()) + 1
        self.pinned = np.zeros(self.nblk, bool)
        for b in pin_blocks:
            if 0 <= int(b) < self.nblk:
                self.pinned[int(b)] = True
        self.T0 = {k: v.copy() for k, v in self.T.items()}
        self.G = np.zeros((15, self.nblk), np.int64)               # G[k, block] = bit vector over loader k's relations
        cnt = {}
        self.rel_pos = []                                           # of every relation of the file: (k, index among loader k's)
        for (k, _o, _t), (_b, gv) in zip(rels, sp):
            ri = cnt.get(k, 0)
            cnt[k] = ri + 1
            self.rel_pos.append((k, ri))
            for b in range(self.nblk):
                self.G[k, b] |= ((gv >> b) & 1) << ri
        self.gflips = 0
        self._sync()

    # ---- the constants
    def _sync(self):
        for k in self.T0:
            self.T[k] = self.T0[k] ^ self.G[k][self.blk]

    def _pair_min(self, cb, b):
        """A[row, sigma, s2] = the least batch cost over the codes of batch b with syn_in = sigma and syn_out = s2."""
        si = self.syn_in[b]
        so = self.syn_out[b] if b < self.hi else np.zeros(len(si), np.int64)
        ns, ns2 = self._ns(b), self._ns(b + 1)
        key = si * ns2 + so
        od = np.argsort(key, kind="stable")
        ks = key[od]
        st = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
        A = np.full((cb.shape[0], ns * ns2), BIG)
        A[:, ks[st]] = np.minimum.reduceat(cb[:, od], st, axis=1)
        return A.reshape(cb.shape[0], ns, ns2)

    def _messages(self, wc):
        """the programme's pieces for all rows at their places: A[b], the backward tables Bq[b][row, sigma] (least cost
        of batches b..hi when loader b's syndrome of batch b is sigma)."""
        rows = np.arange(self.R * self.B)
        A = {b: self._pair_min(self._batchcost(wc, rows, b), b) for b in self.batches}
        return A

    def fit_g(self, wc):
        """one Gauss-Seidel sweep over the loaders; returns the number of (loader, block) constant vectors changed."""
        lo, hi = self.lo, self.hi
        I = self.R * self.B
        ar = np.arange(I)
        places = self.place
        A = self._messages(wc)
        Hn = np.zeros((I, 1))
        Bq = {}
        for b in range(hi, lo, -1):
            Bq[b] = np.minimum((A[b] + Hn[:, None, :]).min(axis=2), BIG)
            sig = np.arange(Bq[b].shape[1])
            Hn = Bq[b][ar[:, None], sig[None, :] ^ self.T[b][places][:, None]]
        F = A[lo][:, 0, :]
        changed = 0
        for k in range(lo + 1, hi + 1):
            ns = self._ns(k)
            sig = np.arange(ns)
            if ns > 1 and self.gfit:
                cost = np.empty((I, ns))
                for tau in range(ns):
                    cost[:, tau] = (F[:, sig ^ tau] + Bq[k]).min(axis=1)
                cost = np.minimum(cost, 1e25)
                ch = self._choose(k, cost, places, ar)
                if ch:
                    changed += ch
                    self._sync()
            Tb = self.T[k][places]
            F = np.minimum((F[ar[:, None], sig[None, :] ^ Tb[:, None]][:, :, None] + A[k]).min(axis=1), BIG)
        self.gflips += changed
        self.last_total = float(np.minimum(F[:, 0], 1e25).sum())
        return changed

    def _choose(self, k, cost, places, ar):
        """per block the constant vector of loader k with the least sum of the block's rows; returns the number changed."""
        ns = cost.shape[1]
        blk = self.blk[places]
        t0 = self.T0[k][places]
        tot = np.zeros((self.nblk, ns))
        for gam in range(ns):
            tot[:, gam] = np.bincount(blk, weights=cost[ar, t0 ^ gam], minlength=self.nblk)
        bi = np.arange(self.nblk)
        curg = self.G[k]
        curtot = tot[bi, curg]
        best = tot.argmin(axis=1)
        besttot = tot[bi, best]
        take = (besttot < curtot - 1e-12 * np.abs(curtot)) & (besttot < 1e17) & (np.bincount(blk, minlength=self.nblk) > 0) & ~self.pinned
        if take.any():
            self.G[k] = np.where(take, best, curg)
        return int(take.sum())

    def step(self, wc):
        if not self.gfit:
            return super().step(wc)
        for _ in range(max(1, self.rounds)):
            ch = self.fit_g(wc)
            before = self.place.copy()
            super().step(wc)
            if ch == 0 and np.array_equal(before, self.place):
                return
        self.fit_g(wc)

    # ---- what the tools around the placer need
    def rels_used(self):
        """the relations of the file, in file order, with their fitted targets."""
        out = []
        for (k, o, t), (k2, ri) in zip(self.rels, self.rel_pos):
            assert k == k2
            gv = 0
            for b in range(self.nblk):
                gv |= int((self.G[k, b] >> ri) & 1) << b
            out.append((k, o, t if gv == 0 else f"G/{t}/{gv:x}"))
        return out

    def get_state(self):
        return dict(fp_place=self.place.copy(), fp_G=self.G.copy())

    def set_state(self, place, G=None):
        self.place = np.array(place, np.int64).copy()
        if G is not None:
            self.G = np.array(G, np.int64).copy()
        self._sync()

    def broken(self, words):
        """words[row, column] of the rows in INIT order: relations broken at the rows' addresses with the constants."""
        od = self.order(np.arange(self.R * self.B))
        return check_g(np.asarray(words, np.int64)[od], self.paddr, self.rels_used())


class RelPlacerGF(RelPlacerG):
    """THE MODEL'S NEXT STEP (not wired into the rounding tool): the low table too. The loader's free family is
    f(address & 15) XOR g(address >> 4); RelPlacerG fits g; this class also fits f, one bit per (relation, low value),
    shared by the about 71 rows of the table with that low value (all groups: in a one-group model 8 rows). The
    one-body rows see f at every address: the one-body tools would have to follow it (not built)."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.pos = self.paddr & 15
        self.Fq = np.zeros((15, 16), np.int64)
        self.fflips = 0

    def _sync(self):
        fq = getattr(self, "Fq", None)
        for k in self.T0:
            self.T[k] = self.T0[k] ^ self.G[k][self.blk] ^ (fq[k][self.paddr & 15] if fq is not None else 0)

    def _choose(self, k, cost, places, ar):
        ns = cost.shape[1]
        ch = 0
        for kind in ("g", "f"):
            cls = (self.blk if kind == "g" else self.pos)[places]
            ncls = self.nblk if kind == "g" else 16
            other = self.Fq[k][self.pos[places]] if kind == "g" else self.G[k][self.blk[places]]
            t0 = self.T0[k][places] ^ other
            tot = np.zeros((ncls, ns))
            for gam in range(ns):
                tot[:, gam] = np.bincount(cls, weights=cost[ar, t0 ^ gam], minlength=ncls)
            ci = np.arange(ncls)
            cur = self.G[k] if kind == "g" else self.Fq[k]
            curtot = tot[ci, cur]
            best = tot.argmin(axis=1)
            take = (tot[ci, best] < curtot - 1e-12 * np.abs(curtot)) & (tot[ci, best] < 1e17) & (np.bincount(cls, minlength=ncls) > 0)
            if kind == "g":
                take &= ~self.pinned
            if take.any():
                if kind == "g":
                    self.G[k] = np.where(take, best, cur)
                else:
                    self.Fq[k] = np.where(take, best, cur)
                    self.fflips += int(take.sum())
                ch += int(take.sum())
        return ch

    def rels_used(self):
        out = []
        for (k, o, t), (k2, ri) in zip(self.rels, self.rel_pos):
            gv = sum(int((self.G[k, b] >> ri) & 1) << b for b in range(self.nblk))
            fv = sum(int((self.Fq[k, q] >> ri) & 1) << q for q in range(16))
            out.append((k, o, t if (gv == 0 and fv == 0) else f"G/{t}/{gv:x}/{fv:x}"))
        return out


P51M = ("2|0:0|93;2|1:0|114;2|2:0|21;2|3:0|109;2|4:0|3;4|0:0|27;4|1:0|111;4|2:0|91;4|3:0|81;4|4:0|14;6|0:0|66;6|1:0|42;"
        "6|2:0|38;6|3:0|22;6|4:0|79;8|0:0|101;8|1:0|110;8|2:0|5;8|3:0|99;8|4:0|49;9|0:0|44;9|1:0|20;9|2:0|4;9|3:0|125;"
        "9|4:0|62;10|0:0|16;10|1:0|125;10|2:0|119;10|3:0|63;10|4:0|52;12|0:0|106;12|1:0|53;12|2:0|121;12|3:0|113;"
        "12|4:0|67;13|0:0|93;13|1:0|93;13|2:0|93;13|3:0|21;13|4:0|21;14|0:0|74;14|1:0|61;14|2:0|70;14|3:0|41;14|4:0|54;"
        "6|1:1 4:1 2:1 0:1 3:1|F:6150;9|3:1 0:1 4:1 2:1 1:1|F:9150;12|0:1 1:1 4:1 2:1 3:1|F:12150;"
        "13|4:1 1:1 2:1 0:1 3:1|F:13150;14|4:1 2:1 3:1|F:14150;14|0:1 2:1 1:1|F:14151")


def _wrap(z, m):
    return np.where(z > HI_W, z - m, np.where(z < LO_W, z + m, z))


def selftest(seed=0, big=True):
    import itertools
    rng = np.random.default_rng(seed)
    # (1) the cost tables of fit_g against the programme itself, and the block optimum against brute force
    R, B = 2, 14
    n = R * B
    small = parse("13 | 0:0 1:0 | F:7;14 | 1:0 2:0 | F:9;14 | 0:1 1:1 | 37;14 | 3:0 4:0 | 3;13 | 2:1 3:1 | 5")
    t = rng.uniform(-20, 20, (n, 75))
    c = rng.uniform(0.5, 2.0, (n, 75))
    fz = (rng.random((n, 75)) < 0.5) & np.isin(np.arange(75) % 5, (0, 2, 4))[None, :]     # no relation is ever complete
    z = np.round(t).astype(np.int64)
    fp = RelPlacerG(R, B, small, pin_blocks=())
    wc = fp._wordcost(t, c, fz, z)
    base, _ = fp._dp(wc, np.arange(n), fp.place, diag=True)
    assert float(base.max()) < 1e17
    fp.gfit = False
    fp.fit_g(wc)
    assert abs(fp.last_total - float(base.sum())) < 1e-9 * max(1.0, float(base.sum())), (fp.last_total, float(base.sum()))
    fp.gfit = True
    for _ in range(4):
        if fp.fit_g(wc) == 0:
            break
    got, _ = fp._dp(wc, np.arange(n), fp.place, diag=True)
    assert abs(fp.last_total - float(got.sum())) < 1e-9 * max(1.0, float(got.sum()))
    # brute force: no single (loader, block) vector lowers the total any further
    for k in (13, 14):
        nr = len(fp.rk[k])
        for b in np.unique(fp.blk):
            keep = fp.G[k, b]
            best = float(got.sum())
            for gam in range(1 << nr):
                fp.G[k, b] = gam
                fp._sync()
                v, _ = fp._dp(wc, np.arange(n), fp.place, diag=True)
                best = min(best, float(v.sum()))
            fp.G[k, b] = keep
            fp._sync()
            assert best >= float(got.sum()) - 1e-9 * max(1.0, float(got.sum())), (k, b, best, float(got.sum()))
    print(f"selftest {seed} (1): fit_g's total = the programme's total ({float(base.sum()):.4f} -> {float(got.sum()):.4f} with the constants; "
          f"{int(fp.gflips)} vectors changed); no single (loader, block) vector lowers it (brute force)")
    # (2) the target kind G/: text round trip
    ru = fp.rels_used()
    assert parse_g(fmt_g(ru)) == ru and any(isinstance(t_, str) and t_.startswith("G/") for (_k, _o, t_) in ru)
    fp2 = RelPlacerG(R, B, parse_g(fmt_g(ru)), pin_blocks=())
    assert all(np.array_equal(fp2.T[k], fp.T[k]) for k in fp.T)
    print(f"selftest {seed} (2): rels_used -> text -> parse_g -> the same targets")
    if not big:
        return
    # (3) a mock staged freeze on P51m's relations (45 wires + 6 bit-1 relations), 3 groups of 127
    R, B, T_ = 3, 127, 75
    rels = parse(P51M)
    out = {}
    for name in ("RelPlacer", "RelPlacerG rounds 1", "RelPlacerG rounds 3"):
        r_ = np.random.default_rng(seed + 1)
        x = r_.uniform(-30, 30, (R * B, T_))
        c = np.exp(r_.normal(0, 2.5, (R * B, T_))) * np.exp(r_.normal(0, 1.5, (R * B, 1)))
        fp = RelPlacerG(R, B, rels, gfit=(name != "RelPlacer"), rounds=(3 if name.endswith("3") else 1))
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
        bad = fp.broken(val)
        od = fp.order(np.arange(R * B))
        assert bad == 0, (name, bad)
        assert sorted(od.tolist()) == list(range(R * B)) and np.all(od // B == np.arange(R * B) // B)
        st_ = fp.get_state()
        fp3 = RelPlacerG(R, B, rels)
        fp3.set_state(st_["fp_place"], st_["fp_G"])
        assert fp3.broken(val) == 0
        fp4 = RelPlacerG(R, B, parse_g(fmt_g(fp.rels_used())), gfit=False)
        fp4.set_state(st_["fp_place"])
        assert fp4.broken(val) == 0 and check_g(val[od], fp.paddr, fp.rels_used()) == 0
        assert not fp.G[:, 4].any()
        out[name] = cost
        print(f"selftest {seed} (3) [{name}]: relations broken 0 on {len(rels)} relations x {R * B} rows; rows kept in their groups; excess cost "
              f"{cost:.1f}; constant vectors changed {fp.gflips}; constants set to 1: {int(sum(bin(int(v)).count('1') for v in fp.G.reshape(-1)))}; "
              f"state and relation text restored: relations broken 0; block 4 untouched", flush=True)
    print(f"selftest {seed} (3): the constants x{out['RelPlacerG rounds 1'] / out['RelPlacer']:.3f} of RelPlacer (rounds 1), "
          f"x{out['RelPlacerG rounds 3'] / out['RelPlacer']:.3f} (rounds 3)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--small", action="store_true", help="parts (1) and (2) only (seconds)")
    a = ap.parse_args()
    if a.selftest:
        selftest(big=not a.small)
