"""invsqrt_ph: GF(2) helpers and the special-bit synthesiser of the small-x path (own code, workstream invsqrt_ph).

Small-x path (idea credited to the best public circuit, manwithacat "HW-25" = gengsr + maoxin123, read in
research/ws_invsqrt_decode.md 3.6): for x < 64 the lo word reaches above position 10 where the lo stream has no
rotation; instead of three more rotations the hi stream delivers hi + (lo >> (17 + j)) at positions 0..9. The
difference to the plain hi bits is a set of "special bits", functions of the 63 small inputs only. Here they are
synthesised from scratch for any tables: each special bit is  L xor sum_i A_i B_i  on the 63 inputs, with L, A_i, B_i
XOR-forms over wires that are alive during the hi stream (the caller supplies their truth tables), one Toffoli
per product and one for the gate with [x < 64]; everything is erased by measurement (Gidney 1709.06648).
"""
from __future__ import annotations

import random


def solve(cols, target):
    """Indices of columns (ints as bit vectors) whose XOR is target, or None."""
    piv = {}
    for i, v in enumerate(cols):
        comb = 1 << i
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h][0]
                comb ^= piv[h][1]
            else:
                piv[h] = (v, comb)
                break
    v, comb = target, 0
    while v:
        h = v.bit_length() - 1
        if h not in piv:
            return None
        v ^= piv[h][0]
        comb ^= piv[h][1]
    return [i for i in range(len(cols)) if (comb >> i) & 1]


class Form:
    """XOR of pool wires (indices into the pool list) plus an optional constant 1."""

    def __init__(self, idx, const=0):
        self.idx, self.const = frozenset(idx), const & 1

    def tt(self, tts, dom):
        v = dom if self.const else 0
        for i in self.idx:
            v ^= tts[i]
        return v


def synth(f, tts, dom, seed, max_r=6, tries=60):
    """f, tts[i]: truth tables as ints over the domain mask dom. Returns (L, [(A, B), ...]) with the fewest
    products found: f = L xor sum A_i B_i on the domain."""
    n = len(tts)
    rng = random.Random(seed)
    base = list(tts) + [dom]
    for r in range(0, max_r + 1):
        for _ in range(1 if r == 0 else tries):
            As = []
            for _k in range(r):
                idx = [i for i in range(n) if rng.random() < 0.5]
                if len(idx) < 1:
                    idx = [rng.randrange(n)]
                As.append(Form(idx, rng.randrange(2)))
            cols = list(base)
            for A in As:
                a = A.tt(tts, dom)
                cols += [a & c for c in base]
            sol = solve(cols, f & dom)
            if sol is None:
                continue
            L = Form([i for i in sol if i < n], int(n in sol))
            prods = []
            ok = True
            for k, A in enumerate(As):
                o = (n + 1) * (k + 1)
                B = Form([i - o for i in sol if o <= i < o + n], int(o + n in sol))
                if not B.idx and not B.const:
                    continue                                   # product not used
                if not B.idx:                                  # A AND 1 = A: linear
                    L = Form(L.idx ^ A.idx, L.const ^ A.const)
                elif B.idx == A.idx:                           # A AND (A xor c): A if c differs from A's constant
                    if B.const == A.const:
                        L = Form(L.idx ^ A.idx, L.const ^ A.const)
                elif len(A.idx) == 0:
                    ok = False
                    break
                else:
                    prods.append((A, B))
            if ok:
                return L, prods
    raise RuntimeError("special bit not synthesised; give more pool wires")


class Collected:
    """In-place collection of two forms onto two distinct pool wires (CNOTs only), undone by undo()."""

    def __init__(self, c, wires, A, B):
        self.c, self.ops = c, []
        ra = min(A.idx)
        for i in sorted(A.idx - {ra}):
            self._cx(wires[i], wires[ra])
        if A.const:
            self._x(wires[ra])
        bi, bc = set(B.idx), B.const
        if ra in bi:
            bi = (bi - {ra}) ^ (set(A.idx) - {ra})
            bi.add(ra)
            bc ^= A.const
        cand = sorted(bi - {ra})
        assert cand, "degenerate product"
        rb = cand[0]
        for i in sorted(bi - {rb}):
            self._cx(wires[i], wires[rb])
        if bc:
            self._x(wires[rb])
        self.a, self.b = wires[ra], wires[rb]

    def _cx(self, a, b):
        self.c.cx(a, b)
        self.ops.append(("cx", a, b))

    def _x(self, a):
        self.c.x(a)
        self.ops.append(("x", a))

    def undo(self):
        for op in reversed(self.ops):
            if op[0] == "cx":
                self.c.cx(op[1], op[2])
            else:
                self.c.x(op[1])


def emit_special(c, wires, L, prods, gate, r):
    """XOR  gate AND (L xor sum A_i B_i)  onto r. Returns fix(m): repair for the outcome bit m of r and erasure of
    the two temporaries. Toffolis: len(prods) + 1."""
    w = c.alloc()
    for i in sorted(L.idx):
        c.cx(wires[i], w)
    if L.const:
        c.x(w)
    for A, B in prods:
        col = Collected(c, wires, A, B)
        c.ccx(col.a, col.b, w)
        col.undo()
    v = c.and_new(gate, w)
    c.cx(v, r)

    def fix(m):
        c.z(v, cond=m)
        c.and_del(v, gate, w)
        m2 = c.hmr(w)
        for i in sorted(L.idx):
            c.z(wires[i], cond=m2)
        if L.const:
            c.neg(cond=m2)
        for A, B in prods:
            col = Collected(c, wires, A, B)
            c.cz(col.a, col.b, cond=m2)
            col.undo()
        c.bit_free(m2)
        c.free(w)

    return fix
