"""chem_loaders compiler library: ANGLE_k lookup circuits for dfthcsos angle-word tables (word width 8, J = 5).

Scheme (the one the shipped public loaders use, re-derived by dissection; see research/ws_chem_loaders.md):
  * address a = 2**lb * h + l.  Every nonconstant function of l and of h is made available on qubits:
    the address bits, cached monomials (one AND each, erased by HMR + conditioned CZ), and "carriers":
    word wires that hold the same pure function of l before and after the circuit in both traversal
    directions (bit7 xor bit6 of a word, or a constant bit column).  A carrier whose function is already
    available is cleaned and reused as a workspace qubit.
  * every other Toffoli is a product (A1(h) + y2(l)) (A2(h) + y1(l)) written straight into the word wires
    with CNOT fan-out; its mixed part is A1 (x) y1 + A2 (x) y2.
  * output bits are first put in a basis adapted to the table's codes (combinations that are affine or
    quadratic in l), found by linear algebra, so coded bits are cheap without being told where they are.
  * two engines find the products and the cheaper result is kept:
      blocks  (the shipped scheme) the low coordinates are cut into blocks of about three; a feature is
              (y, s) with y a nonzero vector of one block and s a fan-out over the output basis, taken inside
              a span that holds the block's row space; two features with the same s share one Toffoli.
      groups  output bits in groups of about six, paired features (y1 (x) s, y2 (x) s) inside a span of the
              group's row space, after "kills" that reuse earlier products. Better for small tables.
  * if the caches do not fit under the 128-qubit cap, the low cache is split into phases.

Everything is exact GF(2) algebra on Python integers; a classical check of every wire on every address is
part of emission.  Public entry points: analyze(), decompose(), emit(), compile_circuit().
"""
from __future__ import annotations

import random
from functools import lru_cache

import numpy as np

BROT, J = 8, 5
QUBIT_CAP = 128
MARGIN = 2          # spare pool unknowns required beyond the kill equations
MAX_GROUP = 8       # largest group of generic output bits (2**b fan-out patterns are enumerated)


# ----------------------------------------------------------------------------------------------- GF(2)
def bits(x: int) -> list[int]:
    out = []
    while x:
        low = x & -x
        out.append(low.bit_length() - 1)
        x ^= low
    return out


class Span:
    """Echelon basis of GF(2) vectors (Python ints) with a tag (combination) carried along."""
    __slots__ = ("piv",)

    def __init__(self, piv=None):
        self.piv = {} if piv is None else piv

    def copy(self):
        return Span(dict(self.piv))

    def reduce(self, v, c=0):
        piv = self.piv
        while v:
            e = piv.get(v.bit_length() - 1)
            if e is None:
                return v, c
            v ^= e[0]
            c ^= e[1]
        return 0, c

    def add(self, v, tag=0):
        r, c = self.reduce(v, tag)
        if r:
            self.piv[r.bit_length() - 1] = (r, c)
            return True
        return False

    def __len__(self):
        return len(self.piv)


def nullspace(rows, n):
    """A basis of {x in GF(2)^n : parity(x & r) = 0 for every r in rows}."""
    piv = {}
    for v in rows:
        while v:
            p = v.bit_length() - 1
            e = piv.get(p)
            if e is None:
                piv[p] = v
                break
            v ^= e
    order = sorted(piv)
    for p in order:                       # clear the other pivot columns from every pivot row
        v = piv[p]
        for q in order:
            if q >= p:
                break
            if v >> q & 1:
                v ^= piv[q]
        piv[p] = v
    out = []
    for f in range(n):
        if f in piv:
            continue
        x = 1 << f
        for p, r in piv.items():
            if r >> f & 1:
                x |= 1 << p
        out.append(x)
    return out


def mobius(arr, nbits):
    """In-place GF(2) Moebius transform along the last axis (length 2**nbits)."""
    sh = arr.shape[:-1]
    for i in range(nbits):
        v = arr.reshape(sh + (-1, 2, 1 << i))
        v[..., 1, :] ^= v[..., 0, :]
    return arr


def to_int(bits_arr) -> int:
    """0/1 numpy vector -> Python int with bit i = element i."""
    return int.from_bytes(np.packbits(np.asarray(bits_arr, np.uint8), bitorder="little").tobytes(), "little")


def rows_to_ints(mat) -> list[int]:
    packed = np.packbits(np.asarray(mat, np.uint8), axis=1, bitorder="little")
    return [int.from_bytes(r.tobytes(), "little") for r in packed]


def int_to_bits(x: int, n: int):
    return np.unpackbits(np.frombuffer(x.to_bytes((n + 7) // 8, "little"), np.uint8), bitorder="little")[:n]


# ------------------------------------------------------------------------------------------- analysis
class Analysis:
    """Everything about circuit k that does not depend on the decomposition."""


def analyze(words: np.ndarray, k: int, lb: int = 4) -> Analysis:
    words = np.asarray(words, np.int64)
    A, ncol = words.shape
    K = -(-ncol // J)
    pad = np.zeros((A, K * J), np.int64)
    pad[:, :ncol] = words
    cur = pad[:, J * k:J * k + J]
    prev = pad[:, J * k - J:J * k] if k else np.zeros_like(cur)
    an = Analysis()
    an.k, an.A, an.lb = k, A, lb
    an.n_sel = max(1, (A - 1).bit_length())
    L = 1 << lb
    if A < 2 * L or an.n_sel <= lb:
        raise ValueError("table too small for this address split")
    an.L, an.nl = L, L - 1
    an.H = -(-A // L)
    an.nr = an.H - 1
    an.hbits = an.n_sel - lb
    a = np.arange(A)
    lidx = a & (L - 1)

    def bit(arr, j, b):
        return ((arr[:, j] >> b) & 1).astype(np.uint8)

    def low_anf(f):
        """ANF mask (bit m = monomial m, bit 0 = constant) if f depends on l only, else None."""
        tt = f[:L].copy()
        if not np.array_equal(f, tt[lidx]):
            return None
        return to_int(mobius(tt, lb))

    # --- carriers: wires that hold the same pure function of l on entry in both directions
    carriers = []          # dicts: wire (j, b), conj, anf
    taken = set()
    for j in range(J):
        for b in range(BROT):
            p, c = bit(prev, j, b), bit(cur, j, b)
            if np.array_equal(p, c):
                f = low_anf(p)
                if f is not None:
                    carriers.append(dict(wire=(j, b), conj=False, anf=f))
                    taken.add((j, b))
    for j in range(J):
        if (j, BROT - 1) in taken:
            continue
        p = bit(prev, j, BROT - 1) ^ bit(prev, j, BROT - 2)
        c = bit(cur, j, BROT - 1) ^ bit(cur, j, BROT - 2)
        if np.array_equal(p, c):
            f = low_anf(p)
            if f is not None:
                carriers.append(dict(wire=(j, BROT - 1), conj=True, anf=f))
                taken.add((j, BROT - 1))
    an.carriers = carriers
    an.out_wires = [(j, b) for j in range(J) for b in range(BROT) if (j, b) not in taken]
    delta = np.stack([bit(prev, j, b) ^ bit(cur, j, b) for j, b in an.out_wires]) if an.out_wires else \
        np.zeros((0, A), np.uint8)
    an.delta_int = [to_int(d) for d in delta]

    # --- low basis: address bits, independent carriers, computed monomials
    basis = []             # dicts: anf (nonconstant mask), kind, ...
    sp = Span()
    for i in range(lb):
        basis.append(dict(anf=1 << (1 << i), kind="lin", bit=i))
        sp.add(1 << (1 << i), 1 << i)
    for c in carriers:
        nc = c["anf"] & ~1
        r, combo = sp.reduce(nc)
        c["const"] = c["anf"] & 1
        if r == 0:
            c["clean"], c["combo"] = True, bits(combo)
        else:
            c["clean"] = False
            c["basis"] = len(basis)
            basis.append(dict(anf=nc, kind="car", wire=c["wire"]))
            sp.add(nc, 1 << (len(basis) - 1))
    comp = []
    for m in sorted(range(1, L), key=lambda m: (m.bit_count(), m)):
        if m.bit_count() >= 2 and sp.reduce(1 << m)[0]:
            basis.append(dict(anf=1 << m, kind="comp", mono=m))
            sp.add(1 << m, 1 << (len(basis) - 1))
            comp.append(len(basis) - 1)
    assert len(basis) == an.nl
    an.basis, an.comp = basis, comp
    an.level = [max(m.bit_count() for m in bits(bs["anf"])) for bs in basis]
    mono_coords = [0] * L
    for m in range(1, L):
        r, combo = sp.reduce(1 << m)
        assert r == 0
        mono_coords[m] = combo
    an.mono_coords = mono_coords

    # --- caches and qubit budget
    an.hi_list = [M for M in range(1, an.H) if M & (M - 1)]
    an.n_clean = sum(1 for c in carriers if c["clean"])
    an.slots = QUBIT_CAP - an.n_sel - J * BROT + an.n_clean
    an.slots_low = an.slots - len(an.hi_list)
    if an.slots_low < 0 or (comp and an.slots_low < 1):
        raise ValueError(f"caches do not fit: {len(an.hi_list)} high monomials, {an.slots} slots (lb={lb})")
    plan_phases(an)

    # --- mixed ANF tensor of the output wires in low-basis coordinates: T[q, M-1] (M = 1..H-1)
    nw = len(an.out_wires)
    full = np.zeros((nw, 1 << an.n_sel), np.uint8)
    full[:, :A] = delta
    mobius(full, an.n_sel)
    full[:, A:] = 0                      # rows past the table are don't-cares: drop their monomials
    C = full.reshape(nw, 1 << an.hbits, L)[:, 1:an.H, :]
    T = np.zeros((nw, an.nr), np.uint32)
    for m in range(1, L):
        T ^= (C[:, :, m].astype(np.uint32)) * np.uint32(mono_coords[m])
    # --- output basis adapted to the low-degree filtration (coded bits)
    newpos = sorted(range(an.nl), key=lambda t: (an.level[t], t))
    pos_of = {t: i for i, t in enumerate(newpos)}
    nr = an.nr
    spo = Span()
    G = []                 # per wire: mask over basis outputs
    Xcols, levels = [], []
    for q in range(nw):
        vec = 0
        for t in range(an.nl):
            col = to_int((T[q] >> np.uint32(t)) & 1)
            if col:
                vec |= col << (pos_of[t] * nr)
        r, c = spo.reduce(vec)
        if r:
            b = len(Xcols)
            spo.piv[r.bit_length() - 1] = (r, 1 << b)
            G.append(c | (1 << b))
            col = np.zeros(nr, np.uint32)
            for t in range(an.nl):
                cbits = (r >> (pos_of[t] * nr)) & ((1 << nr) - 1)
                if cbits:
                    col |= int_to_bits(cbits, nr).astype(np.uint32) << np.uint32(t)
            Xcols.append(col)
            levels.append(an.level[newpos[(r.bit_length() - 1) // nr]])
        else:
            G.append(c)
    an.G = G
    an.X0 = np.stack(Xcols, axis=1) if Xcols else np.zeros((nr, 0), np.uint32)
    an.out_level = levels
    an.infra = len(an.hi_list) + len(an.comp)
    return an


def plan_phases(an: Analysis) -> None:
    """Split the computed low monomials into a core (always cached) and per-phase exclusive sets so that at
    most an.slots_low are cached at once and every AND has an available parent. Sets an.core, an.excl,
    an.recipe (basis index -> (parent as list of basis indices, address bit)) and an.ymasks."""
    comp, basis = an.comp, an.basis
    base = Span()
    for t, bs in enumerate(basis):
        if bs["kind"] != "comp":
            base.add(bs["anf"], 1 << t)

    def recipe_for(t, sp):
        m = basis[t]["mono"]
        for i in bits(m):
            p = m ^ (1 << i)
            r, combo = sp.reduce(1 << p)
            if r == 0:
                return bits(combo), i
        return None

    def try_plan(core, excl):
        rec = {}
        sp = base.copy()
        for t in core:
            rc = recipe_for(t, sp)
            if rc is None:
                return None
            rec[t] = rc
            sp.add(basis[t]["anf"], 1 << t)
        for ex in excl:
            s2 = sp.copy()
            for t in ex:
                rc = recipe_for(t, s2)
                if rc is None:
                    return None
                rec[t] = rc
                s2.add(basis[t]["anf"], 1 << t)
        return rec

    def by_level(ts):
        return sorted(ts, key=lambda t: (an.level[t], t))

    n_low, cap = len(comp), an.slots_low
    if n_low <= cap:
        core, excl = list(comp), [[]]
        rec = try_plan(core, excl)
        assert rec is not None
    else:
        P = 2
        while P * cap < n_low:
            P += 1
        c = max(0, min(cap - 1, (P * cap - n_low) // (P - 1)))
        rng = random.Random(1234 + an.k)
        rec = None
        for trial in range(4000):
            if trial == 0:
                order = list(comp)
            else:
                order = sorted(comp, key=lambda t: (an.level[t] + 3 * rng.random() * (trial > 200), rng.random()))
            core, rest = by_level(order[:c]), order[c:]
            if trial:
                rng.shuffle(rest)
            sizes, left = [], len(rest)
            for p in range(P):
                s = min(cap - c, -(-left // (P - p)))
                sizes.append(s)
                left -= s
            if left:
                continue
            excl, i = [], 0
            for s in sizes:
                excl.append(by_level(rest[i:i + s]))
                i += s
            rec = try_plan(core, excl)
            if rec is not None:
                break
        if rec is None:
            raise ValueError("no feasible low-cache phase plan")
    an.core, an.excl, an.recipe = core, excl, rec
    always = 0
    for t, bs in enumerate(basis):
        if bs["kind"] != "comp" or t in core:
            always |= 1 << t
    an.ymasks = [always | sum(1 << t for t in ex) for ex in excl]


# ---------------------------------------------------------------------------------------- group model
def _bonus(live: int, ymasks) -> int:
    dims = [(live & ym).bit_count() for ym in ymasks]
    if len(ymasks) == 1:
        return max(0, dims[0] - 2)
    union_other = [0] * len(ymasks)
    for p in range(len(ymasks)):
        for q, ym in enumerate(ymasks):
            if q != p:
                union_other[p] |= ym
    excl = [p for p, ym in enumerate(ymasks) if live & ym & ~union_other[p]]
    if excl:
        return max(0, min(dims[p] for p in excl) - 2)
    return max(0, max(dims) - 2)


def kill_order(coords: int, an: Analysis) -> list[int]:
    """Low coordinates in the order they are killed: phase-exclusive ones first, then highest degree."""
    everywhere = an.ymasks[0]
    for y in an.ymasks:
        everywhere &= y
    return sorted(bits(coords), key=lambda t: (bool(everywhere >> t & 1), -an.level[t], -t))


def group_model(b: int, n: int, coords: int, an: Analysis):
    """(predicted products, column kills, row kills) for a group of b outputs with n pool products.
    A product has a zero column in a killed coordinate one time in four, hence the 3/4."""
    order = kill_order(coords, an)
    nr = an.nr
    best = None
    live = coords
    cap = (3 * n) // 4
    for dl in range(0, min(4, len(order) - 2) + 1):
        if dl:
            live &= ~(1 << order[dl - 1])
            if nr * dl > (n * (4 ** dl - 1)) // 4 ** dl - MARGIN:
                break
        a = live.bit_count()
        bon = _bonus(live, an.ymasks)
        for d in range(0, 13):
            if dl + d and nr * dl + a * d > cap - MARGIN:
                break
            amb = a * b
            dim_w = min(amb, nr - d)
            dim_s = max(dim_w, amb - bon)
            cost = (dim_s + 1) // 2
            if best is None or cost < best[0]:
                best = (cost, dl, d)
    return best


def make_planner(an: Analysis, coords: int):
    @lru_cache(None)
    def cost(b, n):
        return group_model(b, n, coords, an)

    @lru_cache(None)
    def solve(left, n):
        if left == 0:
            return 0, 0
        best = None
        for b in range(1, min(MAX_GROUP, left) + 1):
            c = cost(b, min(n, 400))[0]
            tot = c + solve(left - b, min(n + c, 400))[0]
            if best is None or tot < best[0]:
                best = (tot, b)
        return best

    def plan(left, n):
        b = solve(left, min(n, 400))[1]
        _, dl, d = cost(b, min(n, 400))
        return b, dl, d

    plan.total = lambda left, n: solve(left, min(n, 400))[0]
    return plan


# -------------------------------------------------------------------------------------- group engine
class Decomp:
    def __init__(self, an: Analysis, seed: int = 0, effort: int = 4):
        self.an = an
        self.nr = an.nr
        self.X = an.X0.copy()
        self.prods = []                       # [a1, b1, a2, b2, s, phase]
        self.PM = np.zeros((64, an.nr), np.uint32)
        self.rng = random.Random(1000003 * seed + an.k)
        self.effort = effort
        self.log = []

    def add_product(self, a1, b1, a2, b2, s, phase):
        nr = self.nr
        n = len(self.prods)
        if n == len(self.PM):
            self.PM = np.concatenate([self.PM, np.zeros_like(self.PM)])
        row = int_to_bits(b1, nr).astype(np.uint32) * np.uint32(a1)
        if a2:
            row ^= int_to_bits(b2, nr).astype(np.uint32) * np.uint32(a2)
        self.PM[n] = row
        for o in bits(s):
            self.X[:, o] ^= row
        self.prods.append([a1, b1, a2, b2, s, phase])

    def apply_kills(self, outs, coords, dl, d, shift=0):
        an, nr = self.an, self.nr
        order = kill_order(coords, an)
        killed = order[shift:shift + dl] if shift + dl <= len(order) else order[:dl]
        live = coords
        for t in killed:
            live &= ~(1 << t)
        usable = [i for i, p in enumerate(self.prods) if not ((p[0] | p[2]) & ~coords)]
        if not usable:
            return None
        P = self.PM[usable]
        cols = [0] * len(usable)
        tgts = [0] * len(outs)
        Xg = self.X[:, outs].T                     # (b, nr)
        width = 0
        for t in killed:
            for i, v in enumerate(rows_to_ints((P >> np.uint32(t)) & 1)):
                cols[i] |= v << width
            for i, v in enumerate(rows_to_ints((Xg >> np.uint32(t)) & 1)):
                tgts[i] |= v << width
            width += nr
        for _ in range(d):
            z = np.array([self.rng.random() < 0.5 for _ in range(nr)])
            if not z.any():
                z[self.rng.randrange(nr)] = True
            red = np.bitwise_xor.reduce(P[:, z], axis=1) & np.uint32(live)
            for i, v in enumerate(red.tolist()):
                cols[i] |= v << width
            redx = np.bitwise_xor.reduce(Xg[:, z], axis=1) & np.uint32(live)
            for i, v in enumerate(redx.tolist()):
                tgts[i] |= v << width
            width += an.nl
        sp = Span()
        for i, c in enumerate(cols):
            sp.add(c, 1 << i)
        sols = []
        for t in tgts:
            r, c = sp.reduce(t)
            if r:
                return None
            sols.append(c)
        for o, c in zip(outs, sols):
            idx = bits(c)
            if idx:
                self.X[:, o] ^= np.bitwise_xor.reduce(P[idx], axis=0)
                for i in idx:
                    self.prods[usable[i]][4] ^= 1 << o
        return live

    def solve_group(self, outs, coords, dl=0, d=0, tag=""):
        an, nr, rng = self.an, self.nr, self.rng
        b = len(outs)
        live = None
        while dl or d:
            for shift in range(6):       # a rank-deficient system: try other coordinates / directions
                live = self.apply_kills(outs, coords, dl, d, shift)
                if live is not None:
                    break
            if live is not None:
                break
            if d:
                d -= 1
            else:
                dl -= 1
        if live is None:
            live = coords
        lc = bits(live)
        a = len(lc)
        Xg = self.X[:, outs]
        assert not (Xg & ~np.uint32(live)).any(), "group has content outside its live coordinates"
        rows = []
        for M in range(nr):
            v = 0
            for gi in range(b):
                x = int(Xg[M, gi])
                if x:
                    for ci, t in enumerate(lc):
                        if x >> t & 1:
                            v |= 1 << (gi * a + ci)
            rows.append(v)
        W = Span()
        wb = [v for v in rows if W.add(v)]
        dim_w, amb = len(wb), a * b
        if dim_w == 0:
            return 0
        rel = nullspace(wb, amb)
        ylocal, yphase = [], []
        for q, ym in enumerate(an.ymasks):
            m = sum(1 << ci for ci, t in enumerate(lc) if ym >> t & 1)
            if m and m not in ylocal:
                ylocal.append(m)
                yphase.append(q)
        full_a = (1 << a) - 1
        bon = _bonus(live, an.ymasks)
        e_top = min(len(rel), bon)
        lower = (dim_w + 1) // 2
        best = None
        for e in range(e_top, -1, -1):
            if best is not None and e < e_top - 1:
                break                         # codimensions below the top two are only a fallback
            trials = self.effort if e < len(rel) or len(ylocal) > 1 else max(2, self.effort // 2)
            for _ in range(trials):
                if e == len(rel):
                    chosen_rel = rel
                else:
                    chosen_rel, ind = [], Span()
                    while len(chosen_rel) < e:
                        v = 0
                        for r in rel:
                            if rng.random() < 0.5:
                                v ^= r
                        if ind.add(v):
                            chosen_rel.append(v)
                keys = {}
                for s in range(1, 1 << b):
                    gis = bits(s)
                    cons = []
                    for r in chosen_rel:
                        v = 0
                        for gi in gis:
                            v ^= (r >> (gi * a)) & full_a
                        if v:
                            cons.append(v)
                    for p, ym in enumerate(ylocal):
                        if ym == full_a:
                            ys = nullspace(cons, a)
                        else:
                            ys = nullspace(cons + [1 << ci for ci in range(a) if not ym >> ci & 1], a)
                        if ys:
                            keys[(s, p)] = [(y, sum(y << (gi * a) for gi in gis)) for y in ys]
                res = self._select(keys, wb, dim_w)
                if res is not None and (best is None or res[0] < best[0]):
                    best = res
                if best is not None and best[0] <= lower:
                    break
            if best is not None and best[0] <= max(lower, (amb - e_top + 1) // 2):
                break
        if best is None:
            raise RuntimeError("no feature set spans the group")
        _, chosen = best
        fs = Span()
        for idx, (_, _, flat) in enumerate(chosen):
            ok = fs.add(flat, 1 << idx)
            assert ok
        Acoef = [0] * len(chosen)
        for M, v in enumerate(rows):
            r, c = fs.reduce(v)
            assert r == 0
            for idx in bits(c):
                Acoef[idx] |= 1 << M
        bykey = {}
        for idx, (key, y, _) in enumerate(chosen):
            if Acoef[idx]:
                bykey.setdefault(key, []).append(idx)
        n0 = len(self.prods)
        for (s, p), idxs in bykey.items():
            phase = yphase[p]
            sfull = sum(1 << outs[gi] for gi in bits(s))
            for i in range(0, len(idxs), 2):
                i1 = idxs[i]
                y1 = sum(1 << lc[ci] for ci in bits(chosen[i1][1]))
                if i + 1 < len(idxs):
                    i2 = idxs[i + 1]
                    y2 = sum(1 << lc[ci] for ci in bits(chosen[i2][1]))
                    self.add_product(y1, Acoef[i1], y2, Acoef[i2], sfull, phase)
                else:
                    self.add_product(y1, Acoef[i1], 0, 0, sfull, phase)
        assert not self.X[:, outs].any(), "group not cleared"
        made = len(self.prods) - n0
        self.log.append(f"{tag} b={b} a={a} kills=({dl},{d}) dimW={dim_w} amb={amb} bonus={bon} -> {made}")
        return made

    def _select(self, keys, wb, dim_w):
        """Greedy: independent feature pairs of one (fan-out, phase) first, then singles, until the span
        holds every row. Returns (products, chosen features)."""
        rng = self.rng
        order = list(keys)
        rng.shuffle(order)
        F = Span()
        chosen = []
        used = {k: 0 for k in order}
        piv = F.piv

        def covered():
            if len(piv) < dim_w:
                return False
            red = F.reduce
            for v in wb:
                if red(v)[0]:
                    return False
            return True

        done = False
        progress = True
        while progress and not done:
            progress = False
            for key in order:
                first = None
                for y, flat in keys[key]:
                    r, _ = F.reduce(flat)
                    if not r:
                        continue
                    if first is None:
                        first = (y, flat, r.bit_length() - 1)
                        piv[first[2]] = (r, 0)
                    else:
                        piv[r.bit_length() - 1] = (r, 0)
                        chosen.append((key, first[0], first[1]))
                        chosen.append((key, y, flat))
                        used[key] += 2
                        first = "paired"
                        progress = True
                        break
                if first is not None and first != "paired":
                    del piv[first[2]]
                if first == "paired" and covered():
                    done = True
                    break
        if not done:
            for key in order:
                for y, flat in keys[key]:
                    r, _ = F.reduce(flat)
                    if r:
                        piv[r.bit_length() - 1] = (r, 0)
                        chosen.append((key, y, flat))
                        used[key] += 1
                        if covered():
                            done = True
                            break
                if done:
                    break
        if not done:
            return None
        return sum((c + 1) // 2 for c in used.values()), chosen

    def run(self):
        an = self.an
        nb = self.X.shape[1]
        allc = (1 << an.nl) - 1
        classes = {}
        for o in range(nb):
            lv = an.out_level[o]
            classes.setdefault(lv if lv <= 2 else 3, []).append(o)
        for lv in (1, 2):
            outs = classes.get(lv, [])
            if not outs:
                continue
            coords = sum(1 << t for t in range(an.nl) if an.level[t] <= lv)
            size = 10 if lv == 1 else 2
            for i in range(0, len(outs), size):
                self.solve_group(outs[i:i + size], coords, tag=f"level{lv}")
        gen = classes.get(3, [])
        if gen:
            plan = make_planner(an, allc)
            while gen:
                b, dl, d = plan(len(gen), len(self.prods))
                self.solve_group(gen[:b], allc, dl, d, tag="generic")
                gen = gen[b:]
        assert not self.X.any()
        return self.prods


# ------------------------------------------------------------------------- block engine (shipped scheme)
def make_blocks(an: Analysis, supp, nb: int):
    """Partition the low coordinates into blocks of at most 4 inside each phase class, balancing the number
    of (coordinate, output) cells per block; the number of blocks is the one the feature model likes best."""
    count = [sum(1 for o in range(nb) if supp[o] >> t & 1) for t in range(an.nl)]
    classes = {}
    for t in range(an.nl):
        if count[t]:
            sig = tuple(bool(ym >> t & 1) for ym in an.ymasks)
            classes.setdefault(sig, []).append(t)
    nr = an.nr

    def feat(ts):
        amb = sum(count[t] for t in ts)
        lim = min(count[t] for t in ts) // 2
        return amb - min(amb - min(nr, amb), lim)

    blocks = []
    for sig, ts in classes.items():
        ts.sort(key=lambda t: (-count[t], an.level[t], t))
        n = len(ts)
        best = None
        for m in range(-(-n // 4), n + 1):
            bins = [[] for _ in range(m)]
            load = [0] * m
            ok = True
            for t in ts:
                order = sorted(range(m), key=lambda i: (load[i], i))
                i = next((i for i in order if len(bins[i]) < 4), None)
                if i is None:
                    ok = False
                    break
                bins[i].append(t)
                load[i] += count[t]
            if not ok or any(not b for b in bins):
                continue
            c = sum(feat(b) for b in bins)
            if best is None or c < best[0]:
                best = (c, bins)
        blocks += best[1]
    return blocks


class BlockDecomp:
    """Products from low-coordinate blocks: features (y, s) with y a nonzero vector of one block and s a
    fan-out over the output basis, inside a span S_b that holds the block's row space; two features with the
    same s (from any two blocks, same phase) share one Toffoli."""

    def __init__(self, an: Analysis, seed: int = 0, effort: int = 4):
        self.an, self.effort = an, effort
        self.rng = random.Random(7919 * seed + an.k + 17)
        self.log = []

    def prepare(self):
        an = self.an
        X = an.X0
        nr, nb = X.shape
        self.nb = nb
        supp = [int(np.bitwise_or.reduce(X[:, o])) if nr else 0 for o in range(nb)]
        self.supp = supp
        self.blocks = make_blocks(an, supp, nb)
        self.binfo = []
        offset = 0
        for ts in self.blocks:
            cell = {}
            for t in ts:
                for o in range(nb):
                    if supp[o] >> t & 1:
                        cell[(t, o)] = len(cell)
            amb = len(cell)
            cols = {key: to_int((X[:, key[1]] >> np.uint32(key[0])) & 1) for key in cell}
            rows = []
            for M in range(nr):
                v = 0
                for key, idx in cell.items():
                    if cols[key] >> M & 1:
                        v |= 1 << idx
                rows.append(v)
            W = Span()
            wb = [v for v in rows if W.add(v)]
            rel = nullspace(wb, amb)
            phases = [p for p, ym in enumerate(an.ymasks) if all(ym >> t & 1 for t in ts)]
            assert phases, "block mixes exclusive coordinates of different phases"
            self.binfo.append(dict(ts=ts, cell=cell, amb=amb, rows=rows, wb=wb, rel=rel, offset=offset,
                                   phases=phases,
                                   nout=min(sum(1 for o in range(nb) if supp[o] >> t & 1) for t in ts)))
            offset += amb

    def default_e(self):
        return [min(len(bi["rel"]), bi["nout"] // 2) for bi in self.binfo]

    def attempt(self, evec):
        """One randomized plan for the given span codimensions. Returns (products, pairs, singles, points)."""
        rng, nb = self.rng, self.nb
        points = []        # dicts: b, yt, y (coordinate mask), outs (mask), cons (masks over outputs)
        for b, bi in enumerate(self.binfo):
            ts, cell, rel = bi["ts"], bi["cell"], bi["rel"]
            e = min(len(rel), max(0, evec[b]))
            if e == len(rel):
                chosen = rel
            else:
                chosen, ind = [], Span()
                while len(chosen) < e:
                    v = 0
                    for r in rel:
                        if rng.random() < 0.5:
                            v ^= r
                    if ind.add(v):
                        chosen.append(v)
            for ymask in range(1, 1 << len(ts)):
                yt = [ts[i] for i in bits(ymask)]
                outs = [o for o in range(nb) if all(self.supp[o] >> t & 1 for t in yt)]
                if not outs:
                    continue
                cons = []
                for r in chosen:
                    v = 0
                    for o in outs:
                        par = 0
                        for t in yt:
                            par ^= r >> cell[(t, o)] & 1
                        if par:
                            v |= 1 << o
                    if v:
                        cons.append(v)
                points.append(dict(b=b, yt=yt, y=sum(1 << t for t in yt), outs=sum(1 << o for o in outs), cons=cons))

        def flat(p, s):
            bi = self.binfo[p["b"]]
            cell = bi["cell"]
            v = 0
            for t in p["yt"]:
                for o in bits(s):
                    v |= 1 << cell[(t, o)]
            return v << bi["offset"]

        for p in points:
            dead = [1 << o for o in range(nb) if not p["outs"] >> o & 1]
            p["sigma"] = nullspace(p["cons"] + dead, nb)
        cands = []
        npnt = len(points)
        phases = [set(self.binfo[p["b"]]["phases"]) for p in points]
        for i in range(npnt):
            pi = points[i]
            sig_i = pi["sigma"]
            if not sig_i:
                continue
            for j in range(i + 1, npnt):
                pj = points[j]
                if not phases[i] & phases[j] or not pi["outs"] & pj["outs"]:
                    continue
                # vectors of sigma_i that also satisfy point j's constraints and stay inside its outputs
                cj, nj, outside = pj["cons"], len(pj["cons"]), ~pj["outs"]
                syn = Span()
                inter = []
                for kk, u in enumerate(sig_i):
                    sv = (u & outside) << nj
                    for ci, c in enumerate(cj):
                        if (u & c).bit_count() & 1:
                            sv |= 1 << ci
                    r, combo = syn.reduce(sv, 1 << kk)
                    if r:
                        syn.piv[r.bit_length() - 1] = (r, combo)
                    else:
                        x = 0
                        for q in bits(combo):
                            x ^= sig_i[q]
                        inter.append(x)
                if not inter:
                    continue
                if len(inter) <= 3:
                    vecs = []
                    for m in range(1, 1 << len(inter)):
                        x = 0
                        for q in bits(m):
                            x ^= inter[q]
                        vecs.append(x)
                else:
                    vecs = list(inter)
                    for _ in range(len(inter)):
                        x = 0
                        for u in inter:
                            if rng.random() < 0.5:
                                x ^= u
                        if x:
                            vecs.append(x)
                for x in vecs:
                    cands.append((i, j, x))
        return self._pair(points, cands, flat) + (points,)

    def _pair(self, points, cands, flat):
        rng = self.rng
        F = Span()
        piv = F.piv
        nblk = len(self.binfo)
        need = [True] * nblk
        cnt = [0] * nblk

        def check(b):
            bi = self.binfo[b]
            if cnt[b] < len(bi["wb"]):
                return
            off = bi["offset"]
            for v in bi["wb"]:
                if F.reduce(v << off)[0]:
                    return
            need[b] = False

        for b in range(nblk):
            if not self.binfo[b]["wb"]:
                need[b] = False
        deg = [0] * len(points)
        for i, j, _ in cands:
            deg[i] += 1
            deg[j] += 1
        keyed = [(min(deg[i], deg[j]) * (0.5 + rng.random()), ci) for ci, (i, j, _) in enumerate(cands)]
        keyed.sort()
        pairs, singles = [], []
        for _, ci in keyed:
            i, j, s = cands[ci]
            bi_, bj_ = points[i]["b"], points[j]["b"]
            if not (need[bi_] and need[bj_]):
                continue
            v1 = flat(points[i], s)
            r1, _ = F.reduce(v1)
            if not r1:
                continue
            k1 = r1.bit_length() - 1
            piv[k1] = (r1, 0)
            v2 = flat(points[j], s)
            r2, _ = F.reduce(v2)
            if not r2:
                del piv[k1]
                continue
            piv[r2.bit_length() - 1] = (r2, 0)
            pairs.append((i, j, s, v1, v2))
            cnt[bi_] += 1
            cnt[bj_] += 1
            check(bi_)
            if bj_ != bi_:
                check(bj_)
        pidx = list(range(len(points)))
        rng.shuffle(pidx)
        for i in pidx:                      # singles for the blocks whose row space is not yet covered
            b = points[i]["b"]
            if not need[b]:
                continue
            for s in points[i]["sigma"]:
                v = flat(points[i], s)
                r, _ = F.reduce(v)
                if r:
                    piv[r.bit_length() - 1] = (r, 0)
                    singles.append((i, s, v))
                    cnt[b] += 1
                    check(b)
                    if not need[b]:
                        break
        if any(need):
            return (1 << 30, pairs, singles)
        return (len(pairs) + len(singles), pairs, singles)

    def run(self):
        self.prepare()
        evec = self.default_e()
        best = self.attempt(evec)
        if self.effort > 1:
            # coordinate search on the span codimension of every block
            improved, evals = True, 0
            while improved and evals < 12 * len(evec):
                improved = False
                for b in range(len(evec)):
                    for step in (-1, 1, -2, 2, -4):
                        e2 = list(evec)
                        e2[b] += step
                        if e2[b] < 0 or e2[b] > len(self.binfo[b]["rel"]):
                            continue
                        res = self.attempt(e2)
                        evals += 1
                        if res[0] < best[0]:
                            best, evec, improved = res, e2, True
                            break
        if self.effort > 2:
            # polish: fresh random spans at the codimensions found, and one block pushed one further
            for it in range(3 * self.effort):
                e2 = list(evec)
                if it % 2:
                    b = self.rng.randrange(len(evec))
                    e2[b] = min(len(self.binfo[b]["rel"]), e2[b] + 1)
                res = self.attempt(e2)
                if res[0] < best[0]:
                    best, evec = res, e2
        while best[0] >= 1 << 30 and any(evec):
            evec = [max(0, e - 2) for e in evec]      # spans too small to be spanned by features: widen
            best = self.attempt(evec)
        self.evec = evec
        count, pairs, singles, points = best
        if count >= 1 << 30:
            raise RuntimeError("block plan does not cover the table")
        feats = []          # (vector, key): key = ("p", pair index, 0/1) or ("s", single index)
        for n, (i, j, s, v1, v2) in enumerate(pairs):
            feats.append((v1, ("p", n, 0)))
            feats.append((v2, ("p", n, 1)))
        for n, (i, s, v) in enumerate(singles):
            feats.append((v, ("s", n)))
        fs = Span()
        for idx, (v, _) in enumerate(feats):
            ok = fs.add(v, 1 << idx)
            assert ok
        coef = [0] * len(feats)
        for bi in self.binfo:
            for M, v in enumerate(bi["rows"]):
                r, c = fs.reduce(v << bi["offset"])
                assert r == 0
                for idx in bits(c):
                    coef[idx] |= 1 << M
        A = {key: coef[idx] for idx, (_, key) in enumerate(feats)}
        prods = []

        def phase_of(*pts):
            common = None
            for p in pts:
                ph = set(self.binfo[p["b"]]["phases"])
                common = ph if common is None else common & ph
            return min(common)

        for n, (i, j, s, _, _) in enumerate(pairs):
            a1, a2 = A[("p", n, 0)], A[("p", n, 1)]
            pi, pj = points[i], points[j]
            if a1 and a2:
                prods.append([pi["y"], a1, pj["y"], a2, s, phase_of(pi, pj)])
            elif a1:
                prods.append([pi["y"], a1, 0, 0, s, phase_of(pi)])
            elif a2:
                prods.append([pj["y"], a2, 0, 0, s, phase_of(pj)])
        for n, (i, s, _) in enumerate(singles):
            a1 = A[("s", n)]
            if a1:
                prods.append([points[i]["y"], a1, 0, 0, s, phase_of(points[i])])
        self.log.append(f"blocks {[len(b['ts']) for b in self.binfo]} cells {[b['amb'] for b in self.binfo]} "
                        f"rowspace {[len(b['wb']) for b in self.binfo]} codim {self.evec} "
                        f"pairs {len(pairs)} singles {len(singles)} -> {len(prods)}")
        return prods


def decompose(an: Analysis, seed: int = 0, effort: int = 4, engine: str = "auto"):
    """Products [y1, A1, y2, A2, fan-out over the output basis, phase] and a log. engine: blocks | groups | auto
    (both, keep the cheaper)."""
    best = None
    log = []
    if an.X0.shape[1] == 0:
        return [], ["no output depends on both address parts"]
    if engine in ("auto", "blocks"):
        try:
            bd = BlockDecomp(an, seed, effort)
            best = bd.run()
            log += bd.log
        except RuntimeError as exc:
            log.append(f"block engine failed: {exc}")
    if engine in ("auto", "groups") or best is None:
        gd = Decomp(an, seed, max(2, effort))
        try:
            prods = gd.run()
        except RuntimeError as exc:
            if best is None:
                raise
            prods = None
            log.append(f"group engine failed: {exc}")
        log += gd.log
        if prods is None:
            log.append("engine: blocks")
        elif best is None or len(prods) < len(best):
            best = prods
            log.append("engine: groups")
        else:
            log.append("engine: blocks")
    return best, log


# -------------------------------------------------------------------------------------------- emission
def emit(an: Analysis, prods) -> tuple[str, dict]:
    """kmx text of the circuit and a stats dict. Checks every wire on every address classically."""
    A, lb, L, nl, n_sel = an.A, an.lb, an.L, an.nl, an.n_sel
    a = np.arange(A)
    lidx, hidx = a & (L - 1), a >> lb
    ALL = (1 << A) - 1
    LBf = []                           # truth tables of the low basis functions
    for bs in an.basis:
        tt = np.zeros(L, np.uint8)
        for m in bits(bs["anf"]):
            tt[m] = 1
        tt = mobius(tt, lb)            # the transform is an involution: ANF -> truth table
        LBf.append(to_int(tt[lidx]))
    HMf = [None] + [to_int((hidx & M) == M) for M in range(1, an.H)]

    def wire_q(w):
        return n_sel + BROT * w[0] + w[1]

    ops = []
    free = sorted(wire_q(c["wire"]) for c in an.carriers if c["clean"]) + \
        list(range(n_sel + J * BROT, QUBIT_CAP))
    used_max = n_sel + J * BROT - 1

    def take():
        nonlocal used_max
        if not free:
            raise ValueError("out of qubits")
        q = free.pop(0)
        used_max = max(used_max, q)
        return q

    def give(q):
        free.append(q)
        free.sort()

    lowq = [None] * nl
    for t, bs in enumerate(an.basis):
        if bs["kind"] == "lin":
            lowq[t] = bs["bit"]
        elif bs["kind"] == "car":
            lowq[t] = wire_q(bs["wire"])
    # carrier preparation (self-inverse gates; undone in reverse at the end)
    prep = []
    for c in an.carriers:
        q = wire_q(c["wire"])
        if c["conj"]:
            prep.append(f"CX q{q - 1} q{q}")
        if c["const"]:
            prep.append(f"X q{q}")
        if c["clean"]:
            for t in c["combo"]:
                prep.append(f"CX q{lowq[t]} q{q}")
    ops += prep
    # high cache
    highq = {}
    made_high = []
    for i in range(an.hbits):
        if (1 << i) < an.H:
            highq[1 << i] = lb + i
    for M in an.hi_list:
        top = M.bit_length() - 1
        t = take()
        c1, c2 = highq[M ^ (1 << top)], lb + top
        ops.append(f"CCX q{c1} q{c2} q{t}")
        highq[M] = t
        made_high.append((c1, c2, t))

    def low_and(t):
        par, i = an.recipe[t]
        qs = [lowq[x] for x in par]
        piv = next(q for q in qs if q != i)
        pre = [f"CX q{q} q{piv}" for q in qs if q != piv]
        tq = take()
        ops.extend(pre)
        ops.append(f"CCX q{piv} q{i} q{tq}")
        ops.extend(reversed(pre))
        lowq[t] = tq
        return (pre, piv, i, tq, t)

    def low_del(rec):
        pre, piv, i, tq, t = rec
        ops.extend(pre)
        ops.append(f"HMR q{tq} b0")
        ops.append(f"CZ q{piv} q{i} if b0")
        ops.extend(reversed(pre))
        lowq[t] = None
        give(tq)

    nw = len(an.out_wires)
    wq = [wire_q(w) for w in an.out_wires]

    def func(low_mask, high_mask):
        f = 0
        for t in bits(low_mask):
            f ^= LBf[t]
        for M in bits(high_mask):
            f ^= HMf[M + 1]
        return f

    # pre-pass: what every product adds to every wire, then the pure correction each wire still needs
    acc = [0] * nw
    targets_of = []
    for a1, b1, a2, b2, s, _ in prods:
        P = func(a1, b2) & func(a2, b1)
        tg = [q for q in range(nw) if (an.G[q] & s).bit_count() & 1]
        assert tg, "product without a target"
        targets_of.append(tg)
        for q in tg:
            acc[q] ^= P
    pure = Span()
    pure.add(ALL, 1)
    for t in range(nl):
        pure.add(LBf[t], 1 << (1 + t))
    for M in range(1, an.H):
        pure.add(HMf[M], 1 << (nl + M))
    corr = []
    for q in range(nw):
        r, c = pure.reduce(an.delta_int[q] ^ acc[q])
        assert r == 0, f"wire {an.out_wires[q]} is not loaded correctly (circuit {an.k})"
        corr.append(c)

    core_made = [low_and(t) for t in an.core]
    by_phase = [[] for _ in an.excl]
    for i, p in enumerate(prods):
        by_phase[p[5]].append(i)
    n_emitted = 0

    def do_product(i):
        nonlocal n_emitted
        a1, b1, a2, b2, s, _ = prods[i]
        s1 = {lowq[t] for t in bits(a1)} | {highq[M + 1] for M in bits(b2)}
        s2 = {lowq[t] for t in bits(a2)} | {highq[M + 1] for M in bits(b1)}
        assert None not in s1 and None not in s2, "product uses a low function that is not cached in its phase"
        assert s1 and s2 and s1 != s2, "degenerate product"
        p1 = min(s1 - s2) if s1 - s2 else min(s1)
        pre = [f"CX q{q} q{p1}" for q in sorted(s1) if q != p1]
        if p1 in s2:
            s2 = (s1 ^ s2) | {p1}
        p2 = min(q for q in s2 if q != p1)
        pre += [f"CX q{q} q{p2}" for q in sorted(s2) if q != p2]
        tg = targets_of[i]
        first = wq[tg[0]]
        fan = [f"CX q{first} q{wq[q]}" for q in tg[1:]]
        ops.extend(pre)
        ops.extend(fan)
        ops.append(f"CCX q{p1} q{p2} q{first}")
        ops.extend(reversed(fan))
        ops.extend(reversed(pre))
        n_emitted += 1

    def low_corrections(ts):
        for q in range(nw):
            for t in ts:
                if corr[q] >> (1 + t) & 1:
                    ops.append(f"CX q{lowq[t]} q{wq[q]}")

    excl_all = {t for ex in an.excl for t in ex}
    for ph, ex in enumerate(an.excl):
        made = [low_and(t) for t in ex]
        for i in by_phase[ph]:
            do_product(i)
        low_corrections(ex)
        if ph == len(an.excl) - 1:
            low_corrections([t for t in range(nl) if t not in excl_all])
            for q in range(nw):
                c = corr[q]
                if c & 1:
                    ops.append(f"X q{wq[q]}")
                for M in range(1, an.H):
                    if c >> (nl + M) & 1:
                        ops.append(f"CX q{highq[M]} q{wq[q]}")
        for rec in reversed(made):
            low_del(rec)
    for rec in reversed(core_made):
        low_del(rec)
    for c1, c2, t in reversed(made_high):
        ops.append(f"HMR q{t} b0")
        ops.append(f"CZ q{c1} q{c2} if b0")
    ops.extend(reversed(prep))

    head = []
    for i in range(n_sel):
        head.append(f"APPEND_TO_REGISTER q{i} r0")
    head.append('REGISTER r0 "selection"')
    for j in range(J):
        for b in range(BROT):
            head.append(f"APPEND_TO_REGISTER q{n_sel + BROT * j + b} r{j + 1}")
        head.append(f'REGISTER r{j + 1} "word_{j}"')
    stats = dict(k=an.k, infra=an.infra, products=n_emitted, toffoli=an.infra + n_emitted,
                 qubits=used_max + 1, phases=len(an.excl), carriers=sum(1 for c in an.carriers if not c["clean"]),
                 clean_wires=an.n_clean, outputs=len(an.out_level), ops=len(ops))
    return "\n".join(head + ops) + "\n", stats


# --------------------------------------------------------------------------------------- entry points
def compile_circuit(args):
    """(words, k, lb, seed, effort, want_text[, engine]) -> stats dict (with 'kmx' when want_text)."""
    words, k, lb, seed, effort, want_text = args[:6]
    engine = args[6] if len(args) > 6 else "auto"
    an = analyze(words, k, lb)
    prods, log = decompose(an, seed, effort, engine)
    stats = dict(k=k, lb=lb, infra=an.infra, products=len(prods), toffoli=an.infra + len(prods), log=log,
                 phases=len(an.excl), carriers=sum(1 for c in an.carriers if not c["clean"]),
                 clean_wires=an.n_clean,
                 levels=[an.out_level.count(1), an.out_level.count(2), sum(1 for x in an.out_level if x > 2)])
    if want_text:
        text, st = emit(an, prods)
        assert st["toffoli"] == stats["toffoli"]
        stats.update(st)
        stats["kmx"] = text
    return stats
