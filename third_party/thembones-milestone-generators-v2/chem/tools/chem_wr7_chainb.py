"""chem_wr7_chainb (WR7): delta wires on BIT 0 AND BIT 1, the address of a row decided at freeze time.

A wire k:j:b:m asks in every row: bit b of word 5k+j XOR bit b of word 5k-5+j = parity(address AND m), b in {0, 1}
(the loader XORs the two batches bitwise, so each bit is its own relation). Wires of one bit in one slot form chains
exactly as in r50_chem_lbl_chain.py (its chains_of is reused, bit by bit); a word belongs to at most one chain per bit.

DynChainB(wires, addrs, nrows): the generalisation of the designer's DynChain to several bits per word:
  choose(row, col, x, lo, hi) -> the integer word nearest to x inside [lo, hi] whose bits keep the row's relations
      feasible (the group's rows still matched one-to-one to its addresses); None if the column carries no wire;
  commit(row, col, z) records the word (its relations restrict the row's addresses; asserts feasibility).
With bit-0 wires only it takes the same decisions as DynChain (nearest word if feasible, else the nearest word of
the other parity). Python ints: (z >> b) & 1 is the two's-complement bit, the same as the stored word mod 128.

    python chem_wr7_chainb.py --loaders 4,6,8,10,12,14 --bits 0,1 [--seed0 11 --seed1 3]   prints a wire list
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r50_chem_lbl_chain import chains_of, lin, mk  # noqa: E402


def parse(s):
    if s.startswith("@"):
        s = open(s[1:]).read().strip()
    return [tuple(int(v) for v in t.split(":")) for t in s.split(",") if t]


def fmt(wires):
    return ",".join(":".join(str(v) for v in w) for w in wires)


def wbit(z, b):
    return (int(z) >> b) & 1


class DynChainB:
    def __init__(self, wires, addrs, nrows):
        self.col = {}                       # col -> [(chain id, cumulative mask, bit)]
        self.chains = []                    # every chain, all bits ([(col, M), ..] as chains_of gives them)
        for b in sorted({w[2] for w in wires}):
            wb = [(k, j, 0, m) for (k, j, bb, m) in wires if bb == b]
            colb, chainsb = chains_of(wb)
            for c, (ci, M) in colb.items():
                self.col.setdefault(c, []).append((len(self.chains) + ci, M, b))
            self.chains += chainsb
        self.nchains = len(self.chains)
        self.addrs = np.asarray(addrs, np.int64)
        self.na = len(self.addrs)
        self.allowed = np.ones((nrows, self.na), bool)
        self.mrow = np.arange(nrows)
        self.maddr = np.full(self.na, -1, np.int64)
        self.maddr[:nrows] = np.arange(nrows)
        self.ref = {}
        self._L = {}
        self.free_ok = self.forced_n = 0

    def L(self, mask):
        if mask not in self._L:
            self._L[mask] = lin(self.addrs, mask)
        return self._L[mask]

    def _augment(self, row, allowed_row, apply):
        mrow, maddr = self.mrow.copy(), self.maddr.copy()
        maddr[mrow[row]] = -1
        mrow[row] = -1
        allowed = self.allowed
        seen = np.zeros(self.na, bool)

        def dfs(r, al):
            for a_ in np.where(al & ~seen)[0]:
                if seen[a_]:
                    continue
                seen[a_] = True
                o = maddr[a_]
                if o < 0 or dfs(o, allowed[o]):
                    maddr[a_] = r
                    mrow[r] = a_
                    return True
            return False

        ok = dfs(row, allowed_row)
        if ok and apply:
            self.mrow, self.maddr = mrow, maddr
        return ok

    def _conds(self, row, col, z):
        out = []
        for (ci, M, b) in self.col.get(col, ()):
            rf = self.ref.get((row, ci))
            if rf is not None:
                Mr, pr = rf
                out.append((M ^ Mr, wbit(z, b) ^ pr))
        return out

    def _restrict(self, row, conds, apply):
        if not conds:
            return True
        al = self.allowed[row].copy()
        for (mask, v) in conds:
            al &= (self.L(mask) == v)
        ok = True if al[self.mrow[row]] else self._augment(row, al, apply)
        if ok and apply:
            self.allowed[row] = al
        return ok

    def choose(self, row, col, x, lo, hi):
        row, col = int(row), int(col)
        if col not in self.col:
            return None
        span = hi - lo + 1
        if span >= 128:     # the full circle: the window wraps (bits 0 and 1 are kept mod 128)
            f = int(np.floor(x))
            cands = [((v - lo) % span) + lo for v in sorted(range(f - 3, f + 5), key=lambda v: (abs(v - x), v))]
        else:               # a range (the signed range): the window around x clipped into it (at least 4 words)
            f = int(np.floor(min(max(x, lo), hi)))
            cands = sorted((v for v in range(f - 3, f + 5) if lo <= v <= hi), key=lambda v: (abs(v - x), v))
        for i, v in enumerate(cands):
            if self._restrict(row, self._conds(row, col, v), False):
                if i == 0:
                    self.free_ok += 1
                else:
                    self.forced_n += 1
                return v
        raise AssertionError("no feasible word (cannot happen: the row's present address satisfies one residue)")

    def commit(self, row, col, z, *_):
        row, col = int(row), int(col)
        if col not in self.col:
            return
        conds = self._conds(row, col, z)
        assert self._restrict(row, conds, True), "commit of an infeasible word"
        for (ci, M, b) in self.col[col]:
            if (row, ci) not in self.ref:
                self.ref[(row, ci)] = (M, wbit(z, b))


def check(words, addrs, wires):
    """Relations broken on a word table (rows in address order): count per wire, summed."""
    w = np.asarray(words, np.int64)
    bad = 0
    for (k, j, b, m) in wires:
        bad += int(np.sum((((w[:, 5 * k + j] ^ w[:, 5 * k - 5 + j]) >> b) & 1) != lin(addrs, m)))
    return bad


def make(loaders, bits=(0, 1), seeds=(11, 3)):
    out = []
    for b, s in zip(bits, seeds):
        out += [(k, j, b, m) for (k, j, _, m) in mk(loaders, seed=s)]
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--loaders", required=True)
    ap.add_argument("--bits", default="0,1")
    ap.add_argument("--seeds", default="11,3", help="mask seed per bit (r50_chem_lbl_chain.mk)")
    a = ap.parse_args()
    print(fmt(make([int(v) for v in a.loaders.split(",")], tuple(int(v) for v in a.bits.split(",")),
                   tuple(int(v) for v in a.seeds.split(",")))))
