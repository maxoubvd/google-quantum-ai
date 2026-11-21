"""r50_chem_lbl_chain: DELTA WIRES IN CHAINS (R50-CHEM-LBL, construction 9b): more bit-0 wires than address bits.

A wire k:j:0:m asks in every row: bit 0 of word 5k+j XOR bit 0 of word 5k-5+j = parity(address AND m). Wires of
neighbouring loaders in one slot form a CHAIN of words c_0, c_1, .. with cumulative masks M_0 = 0, M_i = M_{i-1} XOR m_i:
    parity(word c_i) = t XOR parity(address AND M_i),   t = one free bit per (row, chain).
A row has 7 free address bits (its place in the group) and one free bit per chain; every further relation forces a
word's parity. This module: the wire lists (mk), the chain structure (chains_of), and DynChain = the address decided
at freeze time for chains (the generalisation of DynWires of r50_chem_lbl_wround.py; same matching).

    python r50_chem_lbl_chain.py --loaders 12,13,14 [--slots 0,1,2,3,4] [--seed 0]     prints a wire list
Cumulative masks: distinct non-zero masks on address bits 0..6, the most balanced on the one-body addresses 0..75 first.
"""
import argparse

import numpy as np


def lin(addr, mask):
    x = np.asarray(addr, np.int64) & int(mask)
    out = np.zeros_like(x)
    while np.any(x):
        out ^= x & 1
        x = x >> 1
    return out


def parse(s):
    if s.startswith("@"):
        s = open(s[1:]).read().strip()
    return [tuple(int(v) for v in t.split(":")) for t in s.split(",") if t]


def chains_of(wires):
    """-> (col -> (chain id, cumulative mask)), list of chains [(col, M), ..] in word order."""
    nxt, tails = {}, set()
    for (k, j, b, m) in wires:
        assert b == 0 and k >= 1
        e, f = 5 * k - 5 + j, 5 * k + j
        assert e not in nxt
        nxt[e] = (f, m)
        tails.add(f)
    col, chains = {}, []
    for s in sorted(e for e in nxt if e not in tails):
        c, M, ch = s, 0, [(s, 0)]
        while c in nxt:
            c, m = nxt[c]
            M ^= m
            ch.append((c, M))
        for (c, M) in ch:
            col[c] = (len(chains), M)
        chains.append(ch)
    return col, chains


def mk(loaders, slots=(0, 1, 2, 3, 4), seed=0, nmask=0):
    """A wire list for whole loaders: cumulative masks distinct, balanced on addresses 0..75 first."""
    rng = np.random.default_rng(seed)
    a = np.arange(76)
    cand = sorted(range(1, 128), key=lambda m: (abs(int(lin(a, m).sum()) - 38), rng.random()))
    loaders = sorted(loaders)
    need = len(loaders) * len(slots)
    assert need <= 127, "more wires than masks"
    basis, sel, rest = [], [], []

    def reduce(v):
        for b_ in basis:
            v = min(v, v ^ b_)
        return v

    for m in cand:                                       # the first seven: independent
        r_ = reduce(m)
        if len(basis) < 7 and r_:
            basis.append(r_)
            basis.sort(reverse=True)
            sel.append(m)
        else:
            rest.append(m)
    sel = (sel + rest)[:need]
    cum, out, it = {}, {}, iter(sel)
    for k in reversed(loaders):                          # late loaders take the first masks
        for j in slots:
            cum[(k, j)] = next(it)
    for k in loaders:
        for j in slots:
            out[(k, j)] = (k, j, 0, cum[(k, j)] ^ cum.get((k - 1, j), 0))
    return [out[key] for key in sorted(out)]


def fmt(wires):
    return ",".join(":".join(str(v) for v in w) for w in wires)


class DynChain:
    """The address decided at freeze time, for chains. ask(row, col): None if the word may take its nearest value
    (no chain / first word of its chain in the row / both parities still leave the group a perfect matching), else
    the parity it must take. commit(row, col, z) records the word."""

    def __init__(self, wires, addrs, nrows, static=False):
        self.static = static                             # True: rows stay at their addresses (no row-order freedom)
        self.col, self.chains = chains_of(wires)
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
        stack_guard = [0]

        def dfs(r, al):
            stack_guard[0] += 1
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

    def _restrict(self, row, mask, v, apply):
        al = self.allowed[row] & (self.L(mask) == v)
        ok = True if al[self.mrow[row]] else self._augment(row, al, apply)
        if ok and apply:
            self.allowed[row] = al
        return ok

    def ask(self, row, col):
        row, col = int(row), int(col)
        if col not in self.col:
            return None
        ci, M = self.col[col]
        rf = self.ref.get((row, ci))
        if rf is None:
            return None
        Mr, pr = rf
        vcur = int(self.L(M ^ Mr)[self.mrow[row]])
        if not self.static and self._restrict(row, M ^ Mr, 1 - vcur, False):
            self.free_ok += 1
            return None
        self.forced_n += 1
        return pr ^ vcur

    def commit(self, row, col, z, *_):
        row, col = int(row), int(col)
        if col not in self.col:
            return
        ci, M = self.col[col]
        rf = self.ref.get((row, ci))
        if rf is None:
            self.ref[(row, ci)] = (M, int(z) & 1)
            return
        Mr, pr = rf
        if self.static:
            assert ((int(z) & 1) ^ pr) == int(self.L(M ^ Mr)[self.mrow[row]])
            return
        assert self._restrict(row, M ^ Mr, (int(z) & 1) ^ pr, True)


class StagedChain:
    """DynChain behind the interface that FREE's staged_coded expects (forced / commit)."""

    def __init__(self, wires, addrs, nrows, static=False):
        self.d = DynChain(wires, addrs, nrows, static)

    def forced(self, row, col, low):
        if int(col) not in self.d.col:
            return {}
        p = self.d.ask(row, col)
        return {0: (p, None, None)}

    def commit(self, row, col, low, z, info):
        self.d.commit(row, col, int(z))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--loaders", required=True)
    ap.add_argument("--slots", default="0,1,2,3,4")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    print(fmt(mk([int(v) for v in a.loaders.split(",")], tuple(int(v) for v in a.slots.split(",")), a.seed)))
