"""chem_rk7_rel (RK7): GENERAL GF(2) RELATIONS among the delta outputs of a loader, kept by the rank rule.

A relation of loader k (k = 1..14) is a set of outputs (j, b) (slot j = 0..4, bit b) and an address mask m: in every
row, XOR over its outputs of [bit b of word 5k+j XOR bit b of word 5k-5+j] = parity(address AND m). (m = 0: target 0.)

RELATION FILE: one relation per line, `k | j:b j:b ... | m` (slots 0-BASED: the lead's d1..d5 = slots 0..4); '#' starts
a comment. Example (the lead's first layout on loader 7):  `7 | 0:0 3:0 4:0 | 0`.

RelRank(rels, addrs, nrows): the freeze-time chooser of one group of rows (generalises WR7's DynChainB / R50's DynChain).
  Per row the OPEN relations are an echelon basis (pivot = highest open bit) over the bits of the words not yet
  frozen; each basis row carries (open bits, combined address mask, parity of the frozen bits so far). When a word
  freezes, every combination of open relations whose open bits all lie in that word is COMPLETED: it becomes the
  constraint parity(word AND inner bits) XOR acc = parity(address AND mask) on the row's address (the rank rule: a
  word is forced only when the words still free cannot pay; the stiff words freeze first, so the soft ones pay).
  choose(row, col, x, lo, hi) -> the integer word nearest to x in [lo, hi] whose constraints leave the group a perfect
      matching of rows to addresses (None if the column carries no relation);
  commit(row, col, z) records the word (asserts feasibility). Any commit order of a feasible final table is feasible
      (a resume replays frozen words in table order).
check(words, addrs, rels) -> (relation, row) pairs broken on a table (rows in address order, words mod 2^bits).

    python chem_rk7_rel.py --layout first --loaders 2-14 --out FILE     writes the lead's first layout
    python chem_rk7_rel.py --selftest                                  random freezes: 0 broken, forced = rank theory
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SH = 8          # bit offset per word position in a row's bit vector: bit b of word c is bit SH * c + b


def lin(addr, mask):
    """parity(addr AND mask) per address (same as r50_chem_lbl_chain.lin)."""
    x = np.asarray(addr, np.int64) & int(mask)
    out = np.zeros_like(x)
    while np.any(x):
        out ^= x & 1
        x = x >> 1
    return out


def ttable(t, addrs):
    """The target of a relation per address: an int = parity(address AND mask); 'F:<seed>' = R53-CHEM-RANK's random
    table of the family F, f(address & 15) XOR g(address >> 4) (r53_chem_rank_rel.table, same generator)."""
    addrs = np.asarray(addrs, np.int64)
    if isinstance(t, str) and t.startswith("B:"):      # g(block) only on blocks 0..3 (one-body addresses 0..63), else 0:
        x = int(t[2:], 0)                               # in F, and 0 on every two-body address (a target-0 run keeps it)
        blk = addrs >> 4
        return np.where(blk < 4, (x >> np.minimum(blk, 3)) & 1, 0).astype(np.int64)
    if isinstance(t, str):
        assert t.startswith("F:"), t
        rng = np.random.default_rng([int(t[2:]), 5])
        f, g = rng.integers(0, 2, 16), rng.integers(0, 2, 128)
        return (f[addrs & 15] ^ g[addrs >> 4]).astype(np.int64)
    return lin(addrs, t)


def parse(s):
    """A relation file (or '@file', or the inline form 'k|j:b j:b|m;k|...') -> [(k, ((j, b), ..), target)].
    Lines `k | j:b j:b | mask-or-F:seed` (mine) or `k j:b,j:b F:seed` (R53-CHEM-RANK's text form)."""
    if s.startswith("@"):
        s = open(s[1:]).read()
    lines = []
    for ln in s.splitlines():                    # comments first, then the inline separator ';'
        lines += ln.split("#")[0].split(";")
    out = []
    for ln in lines:
        ln = ln.strip()
        if not ln:
            continue
        if "|" in ln:
            f = [t.strip() for t in ln.split("|")]
            assert len(f) == 3, ("bad relation line", ln)
            k, so, ts = int(f[0]), f[1].split(), f[2]
        else:
            f = ln.split()
            assert len(f) == 3, ("bad relation line", ln)
            k, so, ts = int(f[0]), f[1].split(","), f[2]
        outs = tuple(sorted(tuple(int(v) for v in t.split(":")) for t in so))
        t = ts if ts[:2] in ("F:", "B:") else int(ts, 0)
        assert 1 <= k <= 14 and all(0 <= j <= 4 and 0 <= b <= 6 for (j, b) in outs) and len(set(outs)) == len(outs)
        out.append((k, outs, t))
    return out


def fmt(rels):
    return "\n".join(f"{k} | {' '.join(f'{j}:{b}' for (j, b) in o)} | {m}" for (k, o, m) in rels) + "\n"


def first_layout(loaders):
    """The lead's first layout: on every loader, d1^d4^d5, d2^d4^d5, d3^d4 on bit 0, target 0 (slots 0-based)."""
    out = []
    for k in loaders:
        out += [(k, ((0, 0), (3, 0), (4, 0)), 0), (k, ((1, 0), (3, 0), (4, 0)), 0), (k, ((2, 0), (3, 0)), 0)]
    return out


CODES0 = {      # designed bit-0 codes of one loader (slots 0-based): r relations, any r of the five pairs pay where possible
    1: [((0, 1, 2, 3, 4),)],
    2: [((0, 3, 4),), ((1, 2, 3, 4),)],
    3: [((0, 3, 4),), ((1, 3, 4),), ((2, 3),)],               # the lead's first layout (8 of 10 triples pay)
    4: [((0, 1),), ((1, 2),), ((2, 3),), ((3, 4),)],          # all five deltas equal: any 4 of 5 pairs pay
    5: [((0,),), ((1,),), ((2,),), ((3,),), ((4,),)],         # the five wires
}


def profile_layout(prof):
    """prof 'a-b:r,c-d:r' (r relations of the designed bit-0 code on each loader of the range) -> relations, target 0."""
    out = []
    for part in prof.split(","):
        ks, r = part.split(":")
        lo_, hi_ = (int(v) for v in ks.split("-")) if "-" in ks else (int(ks), int(ks))
        for k in range(lo_, hi_ + 1):
            for (sl,) in CODES0[int(r)]:
                out.append((k, tuple((j, 0) for j in sl), 0))
    return out


def with_ftables(rels, seed=0):
    """Every relation gets its own random table of the family F (R53's generator, seed seed*1000 + index)."""
    return [(k, o, f"F:{seed * 1000 + i}") for i, (k, o, _) in enumerate(rels)]


def with_masks(rels, seed=0, nbits=7):
    """Every relation gets an address mask (the address of a row decided at freeze time absorbs completions): masks on
    the low `nbits` address bits, balanced on the one-body addresses 0..75 (38 +- 2 ones), the masks of one loader
    linearly independent, drawn at random (seed)."""
    rng = np.random.default_rng(seed)
    a = np.arange(76)
    cand = [m for m in range(1, 1 << nbits) if abs(int(lin(a, m).sum()) - 38) <= 2]
    out, by = [], {}
    for (k, o, _) in rels:
        used = by.setdefault(k, [])
        while True:
            m = int(cand[int(rng.integers(len(cand)))])
            basis = []
            for v in used + [m]:
                for b_ in basis:
                    v = min(v, v ^ b_)
                if v:
                    basis.append(v)
            if len(basis) == len(used) + 1:
                break
        used.append(m)
        out.append((k, o, m))
    return out


def relmask(k, outs):
    m = 0
    for (j, b) in outs:
        m ^= (1 << (SH * (5 * k + j) + b)) ^ (1 << (SH * (5 * k - 5 + j) + b))
    return m


def columns(rels):
    """The word columns that carry at least one relation bit."""
    cs = set()
    for (k, o, _) in rels:
        for (j, _b) in o:
            cs |= {5 * k + j, 5 * k - 5 + j}
    return sorted(cs)


def echelon(rels):
    """Echelon basis {pivot: (bits, combination of relation indices, acc)} of the relations (asserts independence).
    The target of a combination is the XOR of its relations' tables."""
    ech = {}
    for i, (k, o, _) in enumerate(rels):
        v, cb = relmask(k, o), 1 << i
        while v:
            h = v.bit_length() - 1
            e = ech.get(h)
            if e is None:
                ech[h] = (v, cb, 0)
                break
            v ^= e[0]
            cb ^= e[1]
        assert v, "dependent relations in the file"
    return ech


class RelRank:
    def __init__(self, rels, addrs, nrows, zone=0):
        self.rels = list(rels)
        self.addrs = np.asarray(addrs, np.int64)
        self.ex = self.addrs < int(zone)        # R53-CHEM-EXEMPT's zone: no relation need hold at these addresses
        self.na = len(self.addrs)
        base = echelon(self.rels)
        self.act = [dict(base) for _ in range(nrows)]
        self.col = {c: True for c in columns(self.rels)}       # the hook marks these columns
        self.chains = self.rels                                 # for the hook's log line
        self.allowed = np.ones((nrows, self.na), bool)
        self.mrow = np.arange(nrows)
        self.maddr = np.full(self.na, -1, np.int64)
        self.maddr[:nrows] = np.arange(nrows)
        self.tabs = [ttable(t, self.addrs).astype(bool) for (_, _, t) in self.rels]
        self.zero = not any(np.any(t) for t in self.tabs)       # every target 0 on these addresses
        self._L = {}
        self.free_ok = self.forced_n = self.ncons = 0
        self.maxbit = max([b for (_, o, _) in self.rels for (_, b) in o] + [0])

    def L(self, combo):
        """The target table of a combination of relations (bit i of combo = relation i)."""
        t = self._L.get(combo)
        if t is None:
            t = np.zeros(self.na, bool)
            c, i = combo, 0
            while c:
                if c & 1:
                    t ^= self.tabs[i]
                c >>= 1
                i += 1
            if len(self._L) < 300000:
                self._L[combo] = t
        return t

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

    def _close(self, row, col):
        """-> (comp, tmp): the combinations COMPLETED by freezing word `col` [(inner bits, address mask, acc)] and the
        touching basis rows reduced to bits outside the word {pivot: (outer bits, inner bits, address mask, acc)}."""
        act = self.act[row]
        sh = SH * col
        pm = 0x7F << sh
        tmp, comp = {}, []
        for h0 in [h for h, r in act.items() if r[0] & pm]:
            v, am, acc = act[h0]
            out, inn = v & ~pm, (v >> sh) & 0x7F
            while out:
                h = out.bit_length() - 1
                e = act.get(h)
                if e is not None and not (e[0] & pm):
                    out ^= e[0]
                    am ^= e[1]
                    acc ^= e[2]
                    continue
                e = tmp.get(h)
                if e is None:
                    break
                out ^= e[0]
                inn ^= e[1]
                am ^= e[2]
                acc ^= e[3]
            if out:
                tmp[out.bit_length() - 1] = (out, inn, am, acc)
            else:
                comp.append((inn, am, acc))
        return comp, tmp

    def _allowed_for(self, row, comp, z):
        al = self.allowed[row].copy()
        for (inn, am, acc) in comp:
            par = (bin(z & inn).count("1") & 1) ^ acc
            al &= (self.L(am) if par else ~self.L(am)) | self.ex
        return al

    def choose(self, row, col, x, lo, hi):
        row, col = int(row), int(col)
        if col not in self.col:
            return None
        comp, _ = self._close(row, col)
        span = hi - lo + 1
        W = max(4, 1 << (self.maxbit + 1))      # every residue mod 2^(maxbit+1) (R53-CHEM-EXEMPT's 22:35Z find)
        if span >= 128:
            f = int(np.floor(x))
            cands = [((v - lo) % span) + lo for v in sorted(range(f - W + 1, f + W + 1), key=lambda v: (abs(v - x), v))]
        else:
            f = int(np.floor(min(max(x, lo), hi)))
            cands = sorted((v for v in range(f - W, f + W + 1) if lo <= v <= hi), key=lambda v: (abs(v - x), v))
        if not comp:
            self.free_ok += 1
            return cands[0]
        for i, v in enumerate(cands):
            al = self._allowed_for(row, comp, v)
            if al.any() and (al[self.mrow[row]] or self._augment(row, al, False)):
                if i == 0:
                    self.free_ok += 1
                else:
                    self.forced_n += 1
                return v
        raise AssertionError("no feasible word")

    def commit(self, row, col, z, *_):
        row, col = int(row), int(col)
        if col not in self.col:
            return
        z = int(z)
        comp, tmp = self._close(row, col)
        act = self.act[row]
        if comp:
            al = self._allowed_for(row, comp, z)
            ok = al.any() and (al[self.mrow[row]] or self._augment(row, al, True))
            assert ok, "commit of an infeasible word"
            self.allowed[row] = al
            self.ncons += len(comp)
        sh = SH * col
        pm = 0x7F << sh
        for h in [h for h, r in act.items() if r[0] & pm]:
            del act[h]
        for h, (out, inn, am, acc) in tmp.items():
            assert h not in act
            act[h] = (out, am, acc ^ (bin(z & inn).count("1") & 1))

    def open_count(self, row):
        return len(self.act[int(row)])


def check(words, addrs, rels, per=False, zone=0):
    """Relations broken on a word table (rows in address order): (relation, row) pairs; per=True -> array per relation.
    zone: addresses below it are exempt (R53-CHEM-EXEMPT's zone) and not checked."""
    w = np.asarray(words, np.int64)
    addrs = np.asarray(addrs, np.int64)
    keep = addrs >= int(zone)
    w, addrs = w[keep], addrs[keep]
    bad = []
    for (k, o, m) in rels:
        acc = np.zeros(len(w), np.int64)
        for (j, b) in o:
            acc ^= ((w[:, 5 * k + j] ^ w[:, 5 * k - 5 + j]) >> b) & 1
        bad.append(int(np.sum(acc != ttable(m, addrs))))
    return np.array(bad) if per else int(sum(bad))


def selftest(seed=0, nrows=40, trials=3):
    """Random freeze orders on random words, with the first layout and with masked relations: 0 broken; with target-0
    relations every row ends with exactly as many forced words as ... (counted) and the system is never infeasible."""
    rng = np.random.default_rng(seed)
    for t in range(trials):
        for kind in ("first", "masked", "bit01", "F"):
            if kind == "first":
                rels = first_layout(range(2, 15))
            elif kind == "F":
                rels = [(k, o, f"F:{int(rng.integers(0, 10 ** 6))}") for (k, o, _) in first_layout(range(2, 15))]
            elif kind == "masked":
                rels = [(k, o, int(rng.integers(0, 128))) for (k, o, _) in first_layout(range(2, 15))]
            else:
                rels = first_layout(range(2, 15)) + [(k, ((0, 1), (1, 1), (2, 1)), int(rng.integers(0, 128)))
                                                     for k in range(2, 15)]
            addrs = 76 + rng.permutation(4 * nrows)[:nrows]
            dyn = RelRank(rels, addrs, nrows)
            x = rng.uniform(-30, 30, (nrows, 75))
            order = rng.permutation(nrows * 75)
            words = np.zeros((nrows, 75), np.int64)
            for flat in order:
                r, c = divmod(int(flat), 75)
                v = dyn.choose(r, c, float(x[r, c]), -32, 31)
                z = int(np.clip(np.round(x[r, c]), -32, 31)) if v is None else int(v)
                assert -32 <= z <= 31
                dyn.commit(r, c, z)
                words[r, c] = z
            perm = np.argsort(dyn.mrow[:nrows])            # new place -> old row
            bad = check(np.mod(words[perm], 128), addrs, rels)
            # replay in table order on a fresh object: no assert may fire
            d2 = RelRank(rels, addrs, nrows)
            for r in range(nrows):
                for c in range(75):
                    d2.commit(r, c, int(words[r, c]))
            print(f"selftest {kind} trial {t}: relations {len(rels)} x rows {nrows}: broken {bad}; free {dyn.free_ok}, "
                  f"forced {dyn.forced_n}, completions {dyn.ncons}; rows moved {int(np.sum(dyn.mrow[:nrows] != np.arange(nrows)))}"
                  f"; replay ok")
            assert bad == 0
            if kind == "first":
                assert dyn.ncons == len(rels) * nrows, (dyn.ncons, len(rels) * nrows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--layout", default="")
    ap.add_argument("--loaders", default="2-14")
    ap.add_argument("--out", default="")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--masks", type=int, default=-1, help="seed: every relation gets a random balanced address mask")
    ap.add_argument("--ftab", type=int, default=-1, help="seed: every relation gets a random table of the family F")
    a = ap.parse_args()
    if a.selftest:
        selftest()
        sys.exit(0)
    if a.layout.startswith("prof:"):
        rels = profile_layout(a.layout[5:])
        txt = f"# profile {a.layout[5:]}: designed bit-0 codes per loader (CODES0), slots 0-based\n"
        if a.ftab >= 0:
            rels = with_ftables(rels, a.ftab)
            txt += f"# targets: random tables of the family F, seed {a.ftab}\n"
        elif a.masks >= 0:
            rels = with_masks(rels, a.masks)
            txt += f"# address masks: random balanced, independent per loader, seed {a.masks}\n"
        txt += fmt(rels)
        if a.out:
            open(a.out, "w").write(txt)
        echelon(rels)
        print(f"# {len(rels)} relations, independent; columns {len(columns(rels))} -> {a.out}")
        sys.exit(0)
    lo_, hi_ = (int(v) for v in a.loaders.split("-"))
    if a.layout == "first":
        rels = first_layout(range(lo_, hi_ + 1))
        txt = "# the lead's first layout (LEAD_NOTES 21:41Z): d1^d4^d5, d2^d4^d5, d3^d4 on bit 0, target 0, slots 0-based\n"
        if a.masks >= 0:
            rels = with_masks(rels, a.masks)
            txt += f"# address masks: random balanced, independent per loader, seed {a.masks}\n"
        txt += fmt(rels)
        if a.out:
            open(a.out, "w").write(txt)
        print(txt, end="")
        echelon(rels)
        print(f"# {len(rels)} relations, independent; columns {len(columns(rels))}")
