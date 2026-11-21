"""r65_chem_all15_hook: THE FREEZE-TIME PLACER WITH FITTED BLOCK CONSTANTS, as a drop-in for the builders' rounding
tools (R65-CHEM-ALL15). A subclass of R50's FreezePlacer (tools/r50_chem_lbl_hook.py; the hook of chem_wr7_froundl.py
`--fp`): same calls, same state, plus one array.

  A bit-0 wire k:j:0:m asks  bit0(word 5k+j) XOR bit0(word 5k-5+j) = parity(address AND m) XOR g_w(block),
  block = address // 16. The loader takes any target f(low 4 bits) XOR g(block) for nothing (R52's family F: the
  plain compiler reads 240 with masks and 240 with random f, g [M R52-CHEM-WIRE2 C2]); g_w is a CNOT fan-in of the
  cached one-hot block wires. HERE g IS FITTED: before every freeze the alternation is (1) the block constants of
  every wired column in turn, by the decoder's own cost at the present places; (2) the linear assignment of every
  group's rows to its places (R50's step). A frozen word of the wrong parity is prohibitive, so a constant rows have
  frozen under never flips.

    fp = FreezePlacerG(R, B, wires)                 # wires only (no block code); G lives in fp.pl.G[block, column]
    r, m = fp.residues(t, c, frozen, value)         # as R50's: the stage's words take r + m * round((t - r) / m)
    fp.block_constants()                            # {(k, j): 0/1 array over the blocks 0..(76 + R B - 1) // 16}
    a resumed run keeps fp.place AND fp.pl.G (set_state / get_state).

--selftest: a mock staged freeze on random angles (3 groups): relations broken 0 with the fitted constants; the
excess cost of fixed places / R50's placer / this placer.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r51_chem_bit1_round as RD  # noqa: E402
from r50_chem_lbl_hook import FreezePlacer  # noqa: E402


class FreezePlacerG(FreezePlacer):
    def __init__(self, R, B, wires, n_one=76, ncols=75, fixed=False, gfit=True, tau=None):
        super().__init__(R, B, [], wires, n_one=n_one, ncols=ncols, fixed=fixed)
        pl = self.pl
        self.gfit = gfit
        self.tau = tau          # None, or the discount length (share of all words) of a word's weight by its distance to its freeze
        pl.G = np.zeros((pl.nb, len(pl.cc)), np.int64)
        first = {}
        for i, e in enumerate(pl.e_of):
            first.setdefault(int(e), i)
        pl.isroot = np.zeros(len(pl.cc), bool)
        pl.isroot[list(first.values())] = True
        pl.Wp0 = pl.Wp.copy()
        self.gflips = 0
        self.wires = [tuple(int(v) for v in w) for w in wires]

    def _sync(self):
        self.pl.Wp = self.pl.Wp0 ^ self.pl.G[self.pl.pblk]

    def fit_g(self, cost, sweeps=2):
        pl, place = self.pl, self.place
        blk_of_row = pl.pblk[place]
        changed = 0
        for g in pl.groups:
            if all(pl.isroot[c] for c in g):
                continue
            c0g, c1g = cost[0][:, g], cost[1][:, g]
            for b in np.unique(blk_of_row):
                rows = np.flatnonzero(blk_of_row == b)
                q = (pl.Wp0[place[rows]][:, g] ^ pl.G[b, g][None, :]).astype(bool)
                a0, a1 = c0g[rows], c1g[rows]

                def tot(qq):
                    return float(np.minimum(np.where(qq, a1, a0).sum(axis=1), np.where(qq, a0, a1).sum(axis=1)).sum())
                cur = tot(q)
                for _ in range(sweeps):
                    better = False
                    for ci, c in enumerate(g):
                        if pl.isroot[c]:
                            continue
                        q[:, ci] ^= True
                        t = tot(q)
                        if t < cur - 1e-12 * (abs(cur) + 1e-300):
                            cur = t
                            pl.G[b, c] ^= 1
                            changed += 1
                            better = True
                        else:
                            q[:, ci] ^= True
                    if not better:
                        break
        self.gflips += changed
        self._sync()
        return changed

    def step(self, cost, maxit=40):
        if not self.gfit:
            return super().step(cost, maxit)
        for _ in range(6):
            ch = self.fit_g(cost)
            before = self.place.copy()
            super().step(cost, maxit)
            if ch == 0 and np.array_equal(before, self.place):
                break
        self.fit_g(cost)

    def residues(self, t, c, frozen, value):
        if self.tau:            # a word far from its freeze has an unknown parity: its weight fades with its rank among the free words
            flat, fr = np.asarray(c, np.float64).reshape(-1), np.asarray(frozen).reshape(-1)
            idx = np.flatnonzero(~fr)
            o = idx[np.argsort(-flat[idx], kind="stable")]
            q = np.zeros(flat.shape)
            q[o] = np.arange(len(o)) / float(flat.size)
            c = (flat * np.exp(-q / self.tau)).reshape(np.shape(c))
        return super().residues(t, c, frozen, value)

    def block_constants(self):
        """the constant of every wire on every block: g_w = G[column 5k+j] XOR G[column 5k-5+j]."""
        pl = self.pl
        idx = {int(c): i for i, c in enumerate(pl.cc)}
        return {(k, j): (pl.G[:, idx[5 * k + j]] ^ pl.G[:, idx[5 * k - 5 + j]]).astype(np.uint8) for (k, j, _, _) in self.wires}

    def get_state(self):
        return dict(fp_place=self.place.copy(), fp_G=self.pl.G.copy())

    def set_state(self, place, G):
        self.place = np.array(place, np.int64).copy()
        self.pl.G = np.array(G, np.int64).copy()
        self._sync()

    def broken(self, words):
        """words[row, column] of the rows (in INIT order): relations broken at the rows' addresses."""
        a = self.addr()
        blk = a // 16 - self.n_one // 16
        gc = self.block_constants()
        bad = 0
        for (k, j, _, m) in self.wires:
            par = np.array([bin(int(x) & int(m)).count("1") & 1 for x in a]) ^ gc[(k, j)][blk]
            bad += int(np.sum(((words[:, 5 * k + j] ^ words[:, 5 * k - 5 + j]) & 1) != par))
        return bad


def selftest(seed=0):
    rng = np.random.default_rng(seed)
    R, B, T_ = 3, 127, 75
    RD.HALF = 1 << 20
    wires = [(k, j, 0, int(rng.integers(1, 128))) for k in (2, 4, 6, 8, 9, 10, 12, 13, 14) for j in range(5)]
    out = {}
    for name in ("fixed", "placer", "placer + g"):
        r_ = np.random.default_rng(seed + 1)
        x = r_.uniform(-30, 30, (R * B, T_))
        c = np.exp(r_.normal(0, 2.5, (R * B, T_))) * np.exp(r_.normal(0, 1.5, (R * B, 1)))
        fp = FreezePlacerG(R, B, wires, fixed=(name == "fixed"), gfit=(name == "placer + g"))
        frozen = np.zeros(x.shape, bool)
        val = np.zeros(x.shape, np.int64)
        rank = np.argsort(-c, axis=None)
        cost = 0.0
        for st in np.array_split(rank, 12):
            t = np.where(frozen, val, x + (0.0 if not frozen.any() else 0.05 * r_.normal(size=x.shape)))
            r, m = fp.residues(t, c, frozen, val)
            allowed = (r + m * np.round((t - r) / m)).astype(np.int64)
            i, j = np.unravel_index(st, x.shape)
            val[i, j] = allowed[i, j]
            cost += float((c[i, j] * ((allowed[i, j] - t[i, j]) ** 2 - (np.round(t[i, j]) - t[i, j]) ** 2)).sum())
            frozen[i, j] = True
        bad = fp.broken(val)
        od = fp.order(np.arange(R * B))
        assert bad == 0, (name, bad)
        assert sorted(od.tolist()) == list(range(R * B)) and np.all(od // B == np.arange(R * B) // B)
        st = fp.get_state()
        fp2 = FreezePlacerG(R, B, wires)
        fp2.set_state(st["fp_place"], st["fp_G"])
        assert fp2.broken(val) == 0
        out[name] = cost
        print(f"selftest [{name}]: relations broken 0 on {len(wires)} wires x {R * B} rows; rows kept in their groups; excess cost {cost:.1f}; "
              f"constants set to 1: {int(fp.pl.G.sum())}; state restored: relations broken 0")
    print(f"selftest: placer x{out['placer'] / out['fixed']:.3f} of fixed places; placer + fitted block constants x{out['placer + g'] / out['fixed']:.3f} "
          f"(x{out['placer + g'] / out['placer']:.3f} of the placer)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
