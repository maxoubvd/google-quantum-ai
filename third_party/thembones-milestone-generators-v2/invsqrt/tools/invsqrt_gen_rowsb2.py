"""invsqrt_rows_b: add-subtract rows of stage B (lo -+= u >> k), of stage A (u window -+= c2 >> k) and their exact
inverses, for our own generator.

Own code of workstream invsqrt_rows, written against tools/invsqrt_gen_circ.py (builder + exact simulator).
Method credits (public design, implemented here from its description in research/ws_invsqrt_decode.md; no gate
list of any public circuit was read into this file):
  - truncated row windows in groups with one pending wire per group, pending wire measured before the reverse rows
    and repaired by one CZ per row, forward carries measured without fix-up and repaired by the mirrored row,
    boundary carries of the reverse rows repaired under their outcome bit, a row emitted twice under
    complementary outcomes: manwithacat "HW-25" (join of gengsr and maoxin123), best public invsqrt circuit.
  - temporary-AND adder, measurement-based uncomputation: Gidney, arXiv:1709.06648.
  - add-subtract by conjugating the accumulator with CNOTs from the digit: Litinski, arXiv:2410.00899 (variant).
Own devices (research/ws_invsqrt_rows.md):
  - the pending wire starts as a copy of the first digit of its group (not constant 1), so its merge is a
    one-bit increment in the frame of that digit: no overflow AND, no constant left in lo;
  - merges are nested (pending wire of group g+1 goes into the window of group g): 10 ANDs instead of 14;
  - forward rows keep two carry wires (every forward carry is measured at once, repaired by the reverse row);
  - tree=True: every forward row works on its own window (operand bits + 1) with a temporary pending wire on
    top; rows of a group run narrowest first, the node of the row below is merged one bit up into the next
    node and measured at once; the reverse rows, which ripple over the group's full window, repair these
    outcomes with one Z each. A forward row then costs operand bits + 1 (carry-in) + 2, whatever its distance
    to the top of its group, with two temporary nodes alive at one time.

Arithmetic of B (exact for every u, t, lo; see rows_model): row with shift k and digit wire t_k does, mod 2^len(lo),
    t_k = 0:  lo += (u >> k) + u[k-2]          t_k = 1:  lo -= (u >> k) + u[k-1]          (u[-1] = u[-2] = 0)
which is the generator's row_virt arithmetic (invsqrt_gen_parts.py) on the same wires.
Arithmetic of A (see rows_model_A): the generator's parts.row_win on the window wires.
"""
from __future__ import annotations

import numpy as np

PUBLIC_GROUPS = ((0, 1, 2), (3, 4, 5), (6, 7, 8), (9, 10, 11, 12, 13))     # partition of the public file
BEST_GROUPS = ((0, 1), (2, 3, 4), (5, 6, 7), (8, 9, 10, 11, 12, 13))      # measured best of 10 partitions
DEFAULT_GROUPS = PUBLIC_GROUPS


# ------------------------------------------------------------------ reference models

def rows_model(u, t, lo, nu=17, nt=14, wlo=20):
    """Bit-exact content of the lo register after emit_B, for integer arrays u (nu bits), t (nt digits, digit i is
    bit i; its row has shift k = nt - 1 - i) and lo (wlo bits). Windows, groups and tree change nothing: the
    result equals full-width rows."""
    u = np.asarray(u, dtype=np.int64)
    t = np.asarray(t, dtype=np.int64)
    acc = np.asarray(lo, dtype=np.int64).copy()
    for k in range(nt):
        i = nt - 1 - k
        add = (u >> k) + (((u >> (k - 2)) & 1) if k >= 2 else 0)
        sub = (u >> k) + (((u >> (k - 1)) & 1) if k >= 1 else 0)
        acc = np.where((t >> i) & 1, acc - sub, acc + add)
    return acc & ((1 << wlo) - 1)


def rows_model_A(c2, t, acc, nt=14, nrows=8, kg=7, wacc=11):
    """Bit-exact content of (window wires + pending wire as top bit), an integer of wacc bits, after emit_A.
    c2: operand (its wires LSB first), t: digits (row with shift k uses digit nt - 1 - k), acc: initial content.
    Row k: acc -+= (c2 >> (k + 1)) + cin, cin = c2[k-1] xor (c2[kg] xor t)(c2[k] xor c2[k-1]), minus when t = 1."""
    c2 = np.asarray(c2, dtype=np.int64)
    t = np.asarray(t, dtype=np.int64)
    acc = np.asarray(acc, dtype=np.int64).copy()
    for k in range(nrows):
        tk = (t >> (nt - 1 - k)) & 1
        lowb = ((c2 >> (k - 1)) & 1) if k >= 1 else 0
        cin = lowb ^ ((((c2 >> kg) & 1) ^ tk) & (((c2 >> k) & 1) ^ lowb))
        v = (c2 >> (k + 1)) + cin
        acc = np.where(tk == 1, acc - v, acc + v)
    return acc & ((1 << wacc) - 1)


# ------------------------------------------------------------------ one adder stage (Gidney temporary AND)
# A stage has an accumulator wire a, an operand wire b (or None) and a carry-in wire ci (or None).
# "entered" state: a holds a ^ ci, b holds b ^ ci, so that carry-out = a * b ^ ci (or a * b, or a * ci).
# `written` says whether a already holds the sum bit a ^ b ^ ci.

def _enter(c, a, b, ci, written):
    if b is not None and ci is not None:
        c.cx(b if written else ci, a)
        c.cx(ci, b)
    elif b is not None:
        if written:
            c.cx(b, a)
    elif ci is not None:
        if written:
            c.cx(ci, a)


def _exit(c, a, b, ci, written):
    if b is not None and ci is not None:
        c.cx(ci, b)
        c.cx(b if written else ci, a)
    elif b is not None:
        if written:
            c.cx(b, a)
    elif ci is not None:
        if written:
            c.cx(ci, a)


def _ctl(a, b, ci):
    return (a, b) if b is not None else (a, ci)


def _carry_into(c, a, b, ci, tgt):
    x, y = _ctl(a, b, ci)
    c.ccx(x, y, tgt)
    if b is not None and ci is not None:
        c.cx(ci, tgt)


def _carry_phase(c, a, b, ci, cond):
    x, y = _ctl(a, b, ci)
    c.cz(x, y, cond=cond)
    if b is not None and ci is not None:
        c.z(ci, cond=cond)


def _carry_erase(c, a, b, ci, w):
    m = c.hmr(w)
    _carry_phase(c, a, b, ci, m)
    c.bit_free(m)
    c.free(w)


def _write(c, a, b, ci):
    if b is not None:
        c.cx(b, a)
    if ci is not None:
        c.cx(ci, a)


# ------------------------------------------------------------------ row description

class RowSpec:
    """Operand and carry-in of one row. Window bit j receives src[base + j]; the carry-in is
    cin = cin_add xor s (cin_add xor cin_sub) with s = t (xor sel): one AND. t is the digit wire (the row
    subtracts when t = 1); cin_add / cin_sub / sel are wires or None (constant 0)."""

    def __init__(self, key, src, base, t, cin_add, cin_sub, sel=None):
        self.key, self.src, self.base, self.t = key, src, base, t
        self.ca, self.cs, self.sel = cin_add, cin_sub, sel
        self.m = len(src) - base              # operand bits
        self.has_cin = cin_add is not None or cin_sub is not None

    def opb(self, j):
        return self.src[self.base + j] if self.base + j < len(self.src) else None

    def _carrier(self):
        if self.ca is None or self.cs is None:
            return (self.ca if self.cs is None else self.cs), None
        return (self.ca, self.cs) if self.cs == self.sel else (self.cs, self.ca)

    def _forms(self, c, undo=False):
        """CNOTs that put s on the t wire and (cin_add xor cin_sub) on one carrier wire; returns the carrier."""
        car, oth = self._carrier()
        ops = []
        if self.sel is not None:
            ops.append((self.sel, self.t))
        if oth is not None:
            ops.append((oth, car))
        for a, b in (reversed(ops) if undo else ops):
            c.cx(a, b)
        return car

    def cin_make(self, c):
        if not self.has_cin:
            return None
        w = c.alloc()
        car = self._forms(c)
        c.ccx(self.t, car, w)
        self._forms(c, undo=True)
        if self.ca is not None:
            c.cx(self.ca, w)
        return w

    def cin_erase(self, c, w):
        if self.ca is not None:
            c.cx(self.ca, w)
        m = c.hmr(w)
        car = self._forms(c)
        c.cz(self.t, car, cond=m)
        self._forms(c, undo=True)
        c.bit_free(m)
        c.free(w)


def _spec_B(u, k, t):
    return RowSpec(k, u, k, t, u[k - 2] if k >= 2 else None, u[k - 1] if k >= 1 else None)


def _spec_A(c2, k, t, kg):
    return RowSpec(("A", k), c2, k + 1, t, c2[k - 1] if k >= 1 else None, c2[k], c2[kg])


# ------------------------------------------------------------------ forward row

def _row_fwd(c, win, top, sp, ghosts):
    """Row `sp` in the frame of its digit on the window `win` (LSB first); the carry out of win[-1] is XORed into
    `top`. Every carry is measured as soon as the next one exists; the outcome bit of the carry into bit j is
    stored in ghosts[(sp.key, j)] and repaired by the reverse row. Two carry wires at one time."""
    t = sp.t
    for a in win:
        c.cx(t, a)
    car = sp.cin_make(c)
    is_cin = True
    n = len(win)
    for j in range(n):
        a, b = win[j], sp.opb(j)
        assert b is not None or car is not None
        _enter(c, a, b, car, False)
        if j == n - 1:
            new = None
            _carry_into(c, a, b, car, top)
        else:
            new = c.alloc()
            _carry_into(c, a, b, car, new)
        _exit(c, a, b, car, False)
        _write(c, a, b, car)
        if car is not None:
            if is_cin:
                sp.cin_erase(c, car)
            else:
                ghosts[(sp.key, j)] = c.hmr(car)
                c.free(car)
        car = new
        is_cin = False
    for a in win:
        c.cx(t, a)


# ------------------------------------------------------------------ one-bit increment (merge of a pending wire)

def _incr(c, bits, x, top):
    """bits += x at bits[0]; the carry out of bits[-1] is XORed into top. len(bits) ANDs, all temporaries erased
    for free."""
    car = [x]
    n = len(bits)
    for j in range(n):
        if j == n - 1:
            c.ccx(bits[j], car[j], top)
        else:
            w = c.alloc()
            c.ccx(bits[j], car[j], w)
            car.append(w)
    for j in range(n - 1, -1, -1):
        if j < n - 1:
            m = c.hmr(car[j + 1])
            c.cz(bits[j], car[j], cond=m)
            c.bit_free(m)
            c.free(car[j + 1])
        c.cx(car[j], bits[j])


def _merge(c, P, tf, bits, top, inverse):
    """Forward: bits +- (pending digit) in the frame of the group's first digit tf; inverse: the same adder on the
    complemented frame, which undoes it with the same carries."""
    c.cx(tf, P)
    for a in bits:
        c.cx(tf, a)
        if inverse:
            c.x(a)
    _incr(c, bits, P, top)
    for a in bits:
        if inverse:
            c.x(a)
        c.cx(tf, a)
    c.cx(tf, P)


# ------------------------------------------------------------------ forward group

def _group_fwd(c, reg, specs, top_wire, top_bit, ghosts, qghost, tree):
    """Rows `specs` (widest first) of one group whose pending wire / top target `top_wire` sits at bit top_bit of
    reg. tree=False: rows in the given order, each on reg[:top_bit]. tree=True: narrowest row first; row i works
    on reg[:m_i + 1] with its own node on top (top_wire itself if m_i + 1 == top_bit, else a temporary wire that
    starts as a copy of the row's digit); before a row runs, the node of the row below is merged into it through
    one bit and measured (outcome in qghost[key] = (bit, node position))."""
    if not tree:
        for sp in specs:
            _row_fwd(c, reg[:top_bit], top_wire, sp, ghosts)
        return
    prev = None
    for sp in reversed(specs):
        a = sp.m + 1
        assert a <= top_bit
        if a == top_bit:
            node, temp = top_wire, False
        else:
            node, temp = c.alloc(), True
            c.cx(sp.t, node)
        if prev is not None:
            pw, pa, psp = prev
            _merge(c, pw, psp.t, reg[pa:a], node, False)
            qghost[psp.key] = (c.hmr(pw), pa)
            c.free(pw)
        _row_fwd(c, reg[:a], node, sp, ghosts)
        prev = (node, a, sp) if temp else None
    if prev is not None:
        pw, pa, psp = prev
        _merge(c, pw, psp.t, reg[pa:top_bit], top_wire, False)
        qghost[psp.key] = (c.hmr(pw), pa)
        c.free(pw)


# ------------------------------------------------------------------ reverse row with a bounded number of carry wires

def plan_blocks(n, p):
    """Sizes of consecutive blocks of carry nodes for n nodes and p wires. One block if n <= p. Otherwise block i
    (0-based) may hold p - i nodes (i boundary carries of earlier blocks are alive next to it); the last blocks
    are filled first. The mean price is (n - p) / 2 Toffolis for any such plan."""
    if n <= p:
        return [n]
    q, cap = 1, p
    while cap < n:
        q += 1
        assert p - q + 1 >= 1, "pool too small"
        cap += p - q + 1
    sizes = [0] * q
    rest = n
    for i in range(q - 1, 0, -1):
        sizes[i] = p - i
        rest -= sizes[i]
    sizes[0] = rest
    assert 1 <= sizes[0] <= p
    return sizes


def _row_rev_core(c, win, top, sp, ghosts, mP, pool, zfix, direct, consume):
    """One emission of the reverse row. direct=True: the carry into the top window bit gets no wire (its AND is
    written straight into win[-1]); only valid when nothing needs the carry out of the window."""
    t = sp.t
    n = len(win)
    for a in win:
        c.x(a)
        c.cx(t, a)
    car = [None] * n                     # car[j]: wire with the carry into bit j (car[0] = carry-in of the row)
    written = [False] * n
    opb = sp.opb
    j_min = 0 if sp.has_cin else 1
    nodes = list(range(j_min, n - 1 if direct else n))

    def compute(j):
        if j == 0:
            car[0] = sp.cin_make(c)
            return
        a, b, ci = win[j - 1], opb(j - 1), car[j - 1]
        w = c.alloc()
        _enter(c, a, b, ci, written[j - 1])
        _carry_into(c, a, b, ci, w)
        _exit(c, a, b, ci, written[j - 1])
        car[j] = w

    def erase(j):
        if j == 0:
            sp.cin_erase(c, car[0])
        else:
            a, b, ci = win[j - 1], opb(j - 1), car[j - 1]
            _enter(c, a, b, ci, written[j - 1])
            _carry_erase(c, a, b, ci, car[j])
            _exit(c, a, b, ci, written[j - 1])
        car[j] = None

    def write(j):
        _write(c, win[j], opb(j), car[j])
        written[j] = True

    def fix_bits(j):
        out = list((zfix or {}).get(j, ()))
        if (sp.key, j) in ghosts:
            out.append(ghosts[(sp.key, j)])
            if consume:
                c.bit_free(ghosts.pop((sp.key, j)))
        return out

    sizes = plan_blocks(len(nodes), pool)
    blocks, pos = [], 0
    for s in sizes:
        blocks.append(nodes[pos:pos + s])
        pos += s
    for bi, blk in enumerate(blocks):
        last = bi == len(blocks) - 1
        for j in blk:
            compute(j)
            if j >= 1:
                for bq in fix_bits(j):
                    c.z(car[j], cond=bq)
        if last:
            if direct:
                assert opb(n - 1) is None and top is None
                a, b, ci = win[n - 2], opb(n - 2), car[n - 2]
                _enter(c, a, b, ci, False)
                _carry_into(c, a, b, ci, win[n - 1])
                for bq in fix_bits(n - 1):
                    _carry_phase(c, a, b, ci, bq)
                _exit(c, a, b, ci, False)
                written[n - 1] = True
            else:
                a, b, ci = win[n - 1], opb(n - 1), car[n - 1]
                if top is not None:
                    _enter(c, a, b, ci, False)
                    _carry_into(c, a, b, ci, top)
                    _exit(c, a, b, ci, False)
                elif mP is not None:
                    _enter(c, a, b, ci, False)
                    _carry_phase(c, a, b, ci, mP)
                    _exit(c, a, b, ci, False)
        for j in reversed(blk):
            if not written[j]:
                write(j)
            if last or j != blk[-1]:
                erase(j)
    if j_min == 1:
        write(0)
    for bi in range(len(blocks) - 2, -1, -1):
        blk = blocks[bi]
        jb, inner = blk[-1], blk[:-1]
        if not inner:
            erase(jb)
            continue
        m = c.hmr(car[jb])
        c.free(car[jb])
        car[jb] = None
        c.push(m)
        for j in inner:
            compute(j)
        a, b, ci = win[jb - 1], opb(jb - 1), car[jb - 1]
        _enter(c, a, b, ci, written[jb - 1])
        _carry_phase(c, a, b, ci, None)
        _exit(c, a, b, ci, written[jb - 1])
        for j in reversed(inner):
            erase(j)
        c.pop()
        c.bit_free(m)
    for a in win:
        c.cx(t, a)
        c.x(a)


def _row_rev_spec(c, win, top, sp, ghosts, mP, pool, zfix=None, tail=False):
    """Inverse of the row `sp` on the window `win`. top: wire that receives the carry out of win[-1], or None when
    that carry only went into a pending wire whose measured outcome is the classical bit mP (repaired here by the
    phase of the same carry). Repairs the forward outcomes ghosts[(sp.key, j)] and, for zfix = {j: [bits]}, the
    outcomes of forward temporary nodes that sat at bit j (Z on the carry into bit j). At most `pool` carry wires.
    tail=True: when the pool is too small for all carries and the top carry is needed only if mP = 1, the row is
    emitted twice under complementary conditions (with / without a wire for the top carry)."""
    n_full = len(win) - (0 if sp.has_cin else 1)
    if tail and top is None and mP is not None and n_full > pool and sp.opb(len(win) - 1) is None \
            and len(win) >= 3:
        c.push(mP)
        _row_rev_core(c, win, None, sp, ghosts, mP, pool, zfix, False, False)
        c.pop()
        c._op("BIT_INVERT", (), bit=mP)
        c.push(mP)
        _row_rev_core(c, win, None, sp, ghosts, None, pool, zfix, True, True)
        c.pop()
        c._op("BIT_INVERT", (), bit=mP)
    else:
        _row_rev_core(c, win, top, sp, ghosts, mP, pool, zfix, False, True)


# ------------------------------------------------------------------ rows joined in the reverse pass

class _Trk:
    """Proxy of the builder that tracks the content of every wire it touches as an affine form over symbols
    (one symbol per wire it first sees, one per AND), so that every temporary can be erased by measurement
    with the phase of its content rebuilt from the wires alive at that moment (Clifford gates only)."""

    def __init__(self, c):
        self.c, self.f, self.nv, self.ands, self.temps = c, {}, 1, {}, []

    def _m(self, w):
        if w not in self.f:
            self.f[w] = 1 << self.nv
            self.nv += 1
        return self.f[w]

    def alloc(self):
        w = self.c.alloc()
        self.f[w] = 0
        self.temps.append(w)
        return w

    def x(self, w):
        self._m(w)
        self.f[w] ^= 1
        self.c.x(w)

    def cx(self, a, b):
        ma = self._m(a)
        self._m(b)
        self.f[b] ^= ma
        self.c.cx(a, b)

    def ccx(self, a, b, tgt):
        v = 1 << self.nv
        self.nv += 1
        self.ands[v] = (self._m(a), self._m(b))
        self._m(tgt)
        self.f[tgt] ^= v
        self.c.ccx(a, b, tgt)

    def _solve(self, target, skip):
        """target (symbol mask) as XOR of wire contents (+ constant): returns (wires, constant bit)."""
        basis = {}                                   # leading symbol -> (mask, set of wires)
        for w, m in self.f.items():
            if w == skip:
                continue
            comb = {w}
            while m >> 1:
                h = m.bit_length() - 1
                if h not in basis:
                    basis[h] = (m, comb)
                    break
                bm, bc = basis[h]
                m ^= bm
                comb = comb ^ bc
        comb = set()
        while target >> 1:
            h = target.bit_length() - 1
            assert h in basis, "content is not a function of the live wires"
            bm, bc = basis[h]
            target ^= bm
            comb = comb ^ bc
        return sorted(comb), target & 1

    def _phase_prod(self, A, ca, B, cb, cond):
        c = self.c
        for a in A:
            for b in B:
                if a == b:
                    c.z(a, cond=cond)
                else:
                    c.cz(a, b, cond=cond)
        if ca:
            for b in B:
                c.z(b, cond=cond)
        if cb:
            for a in A:
                c.z(a, cond=cond)
        if ca and cb:
            c.neg(cond=cond)

    def erase(self, w):
        mask = self.f.pop(w)
        c = self.c
        m = c.hmr(w)
        lin = mask
        for v, (X, Y) in self.ands.items():
            if mask & v:
                lin ^= v
                A, ca = self._solve(X, w)
                B, cb = self._solve(Y, w)
                self._phase_prod(A, ca, B, cb, m)
        A, ca = self._solve(lin, w)
        for a in A:
            c.z(a, cond=m)
        if ca:
            c.neg(cond=m)
        c.bit_free(m)
        c.free(w)


def _rows_rev_joint(c, reg, specs, top_bit, mP, ghosts, qghost, zfix_up):
    """Inverse of the consecutive rows `specs` (widest first) of a tree group, done together: every row is undone
    on its own window reg[:m + 1] with its top carry on a digit wire; the digits are merged one bit up, lowest
    first, and ONE chain takes the result to the top of the group window (the carry out of it is repaired under
    mP). Every carry stays alive until the end and is then erased for free in the one order in which each
    content is still a function of the live wires: chain (top first), then per row, widest first, its digit
    and its carries (top first). zfix_up =
    {bit position: [bits]}: Z on the chain carry into that position (outcomes of forward nodes of rows above).
    Price: sum of (operand bits + 2) + rows - 1 + (top_bit - 1 - widest window) Toffolis."""
    T = _Trk(c)
    D, cars, chain = [], [], []
    for sp in specs:
        a = sp.m + 1
        win = reg[:a]
        for w in win:
            T.x(w)
            T.cx(sp.t, w)
        car = [sp.cin_make(T)]
        for j in range(a):
            _enter(T, win[j], sp.opb(j), car[j], False)
            w = T.alloc()
            _carry_into(T, win[j], sp.opb(j), car[j], w)
            _exit(T, win[j], sp.opb(j), car[j], False)
            if j + 1 < a and (sp.key, j + 1) in ghosts:
                c.z(w, cond=ghosts[(sp.key, j + 1)])
                c.bit_free(ghosts.pop((sp.key, j + 1)))
            car.append(w)
        for j in range(a):
            _write(T, win[j], sp.opb(j), car[j])
        for w in win:
            T.cx(sp.t, w)
            T.x(w)
        D.append(car[a])
        cars.append(car[:a])
    for i in range(len(specs) - 1, -1, -1):
        sp = specs[i]
        if i + 1 < len(specs):
            low = specs[i + 1]
            assert low.m + 2 == sp.m + 1, "joined rows must have consecutive shifts"
            bit = reg[low.m + 1]
            T.x(bit)
            T.cx(low.t, bit)
            T.ccx(bit, D[i + 1], D[i])
            T.cx(D[i + 1], bit)
            T.cx(low.t, bit)
            T.x(bit)
        if sp.key in qghost:
            c.z(D[i], cond=qghost[sp.key][0])
    sp = specs[0]
    a = sp.m + 1
    bits = reg[a:top_bit]
    if not bits:
        if mP is not None:
            c.z(D[0], cond=mP)
    else:
        for w in bits:
            T.x(w)
            T.cx(sp.t, w)
        cw = [D[0]]
        for idx, w in enumerate(bits):
            pos = a + idx + 1                        # position the carry out of this bit goes into
            if pos < top_bit:
                nw = T.alloc()
                T.ccx(w, cw[idx], nw)
                for bq in zfix_up.get(pos, ()):
                    c.z(nw, cond=bq)
                cw.append(nw)
                chain.append(nw)
            elif mP is not None:
                c.cz(w, cw[idx], cond=mP)
        for idx, w in enumerate(bits):
            T.cx(cw[idx], w)
        for w in bits:
            T.cx(sp.t, w)
            T.x(w)
    for w in reversed(chain):
        T.erase(w)
    for i in range(len(specs)):
        T.erase(D[i])
        for w in reversed(cars[i]):
            if w is not None:
                T.erase(w)


def plan_join(specs, top_bit, pool):
    """Blocks of consecutive rows (index ranges into specs, widest first) to undo jointly so that the Toffoli
    count is smallest with at most `pool` clean wires per block. Returns a list of (first, last); single rows
    (first == last) are undone flat."""
    r = len(specs)
    a = [sp.m + 1 for sp in specs]
    flat = [(1 if sp.has_cin else 0) + top_bit - 1 for sp in specs]

    def price(i, j):
        if any(a[k] != a[k + 1] + 1 for k in range(i, j)) or not all(specs[k].has_cin for k in range(i, j + 1)):
            return None
        wires = sum(a[i:j + 1]) + (j - i + 1) + (top_bit - 1 - a[i])
        if wires > pool:
            return None
        return sum(x + 1 for x in a[i:j + 1]) + (j - i) + (top_bit - 1 - a[i])

    best = [(0, [])] * (r + 1)
    for i in range(r - 1, -1, -1):
        cand = (flat[i] + best[i + 1][0], [(i, i)] + best[i + 1][1])
        for j in range(i + 1, r):
            pr = price(i, j)
            if pr is not None and pr + best[j + 1][0] < cand[0]:
                cand = (pr + best[j + 1][0], [(i, j)] + best[j + 1][1])
        best[i] = cand
    return best[0][1]


def _row_rev(c, win, top, u, k, t, ghosts, mP, pool, zfix=None, tail=False):
    """Stable entry for a B row (signature used by the generator's own reverse loop)."""
    _row_rev_spec(c, win, top, _spec_B(u, k, t), ghosts, mP, pool, zfix, tail)


def _group_rev(c, reg, specs, top_bit, top, mP, ghosts, qghost, pool, tree, tail, join=False):
    """Reverse rows of one group in mirror order of _group_fwd, each on reg[:top_bit]. join=True (tree groups
    without a top wire only): consecutive rows are undone jointly where that is cheaper within `pool` wires."""
    keys = [sp.key for sp in specs]

    def zfix_for(idx):
        zfix = {}
        for k2 in keys[:idx + 1]:
            if k2 in qghost:
                zfix.setdefault(qghost[k2][1], []).append(qghost[k2][0])
        return zfix

    if join and tree and top is None:
        for i, j in plan_join(specs, top_bit, pool):
            if i == j:
                _row_rev_spec(c, reg[:top_bit], top, specs[i], ghosts, mP, pool, zfix_for(i), tail)
            else:
                _rows_rev_joint(c, reg, specs[i:j + 1], top_bit, mP, ghosts, qghost, zfix_for(i - 1))
    else:
        order = list(range(len(specs))) if tree else list(range(len(specs) - 1, -1, -1))
        for i in order:
            _row_rev_spec(c, reg[:top_bit], top, specs[i], ghosts, mP, pool, zfix_for(i), tail)
    for sp in specs:
        if sp.key in qghost:
            bq, _ = qghost.pop(sp.key)
            c.z(sp.t, cond=bq)
            c.bit_free(bq)


# ------------------------------------------------------------------ stage B

def _layout(groups, nu, wlo):
    """Per group: (shifts, window width in lo wires, has pending wire)."""
    out = []
    for gi, grp in enumerate(groups):
        if gi == 0:
            out.append((tuple(grp), wlo - 1, False))
        else:
            out.append((tuple(grp), nu + 1 - grp[0], True))
    return out


def _pool_of(pool, gi):
    if callable(pool):
        return pool(gi)
    if isinstance(pool, dict):
        return pool.get(gi, pool.get("default", 10))
    return pool


def emit_B(c, t, u, lo, groups=DEFAULT_GROUPS, tree=False, hook=None):
    """Forward rows and merge. t: digit wires (t[i] = digit i, row shift k = len(t) - 1 - i), u: operand wires LSB
    first, lo: accumulator wires LSB first. Returns the context for emit_B_inv; ctx['P'] are the pending wires that
    stay allocated (and non-zero) until emit_B_inv. hook(event, group index) is called at 'group_start',
    'group_end' and ('merge', None); nothing but the pending wires is allocated at these points."""
    nt, nu, wlo = len(t), len(u), len(lo)
    assert sorted(k for g in groups for k in g) == list(range(nt))
    lay = _layout(groups, nu, wlo)
    ghosts, pend, qghost, specs = {}, {}, {}, {}
    for gi, (grp, w, has_p) in enumerate(lay):
        if hook:
            hook("group_start", gi)
        specs[gi] = [_spec_B(u, k, t[nt - 1 - k]) for k in grp]
        if has_p:
            tf = t[nt - 1 - grp[0]]
            P = c.alloc()
            c.cx(tf, P)
            pend[gi] = (P, tf, w)
            _group_fwd(c, lo, specs[gi], P, w, ghosts, qghost, tree)
        else:
            _group_fwd(c, lo, specs[gi], lo[w], w, ghosts, qghost, tree)
        if hook:
            hook("group_end", gi)
    if hook:
        hook("merge", None)
    for gi in range(len(lay) - 1, 0, -1):
        P, tf, w = pend[gi]
        if gi == 1:
            _merge(c, P, tf, lo[w:wlo - 1], lo[wlo - 1], False)
        else:
            Pup, _, wup = pend[gi - 1]
            _merge(c, P, tf, lo[w:wup], Pup, False)
    return dict(ghosts=ghosts, pend=pend, lay=lay, qghost=qghost, specs=specs, tree=tree,
                P=[pend[g][0] for g in sorted(pend)])


def emit_B_inv(c, t, u, lo, ctx, pool=10, hook=None, tail=True, join=False):
    """Exact inverse of emit_B: un-merge, pending wires measured (repaired inside the rows), rows in mirror order.
    pool: number of clean wires the rows of a group may hold at one time: an int, a dict {group index: n,
    'default': n} or a function of the group index. hook(event, group index) is called at ('unmerged', None)
    (pending wires measured), 'group_start' and 'group_end' (groups run from the last one to group 0); nothing
    of this module is allocated at these points."""
    nt, nu, wlo = len(t), len(u), len(lo)
    lay, pend, ghosts, qghost = ctx["lay"], ctx["pend"], ctx["ghosts"], ctx["qghost"]
    for gi in range(1, len(lay)):
        P, tf, w = pend[gi]
        if gi == 1:
            _merge(c, P, tf, lo[w:wlo - 1], lo[wlo - 1], True)
        else:
            Pup, _, wup = pend[gi - 1]
            _merge(c, P, tf, lo[w:wup], Pup, True)
    mP = {}
    for gi in range(1, len(lay)):
        P, tf, w = pend[gi]
        mP[gi] = c.hmr(P)
        c.z(tf, cond=mP[gi])
        c.free(P)
    if hook:
        hook("unmerged", None)
    for gi in range(len(lay) - 1, -1, -1):
        grp, w, has_p = lay[gi]
        if hook:
            hook("group_start", gi)
        _group_rev(c, lo, ctx["specs"][gi], w, None if has_p else lo[w], mP.get(gi), ghosts, qghost,
                   _pool_of(pool, gi), ctx["tree"], tail, join)
        if has_p:
            c.bit_free(mP[gi])
        if hook:
            hook("group_end", gi)
    assert not ghosts and not qghost, "unrepaired forward outcomes"


# ------------------------------------------------------------------ stage A (rows on the low window of u)

def emit_A(c, t, win, pend, c2, nrows=8, kg=7, tree=True):
    """Stage A rows: (win + pend as top bit) -+= (c2 >> (k + 1)) + cin for k = 0..nrows-1, digit of row k =
    t[len(t) - 1 - k]; same arithmetic as the generator's parts.row_win loop (see rows_model_A). win: window
    wires LSB first (bits 1.. of u18), pend: the pending wire on top of the window (any content), c2: operand
    wires LSB first. Returns the context for emit_A_inv."""
    nt = len(t)
    specs = [_spec_A(c2, k, t[nt - 1 - k], kg) for k in range(nrows)]
    ghosts, qghost = {}, {}
    _group_fwd(c, win, specs, pend, len(win), ghosts, qghost, tree)
    return dict(ghosts=ghosts, qghost=qghost, specs=specs, tree=tree)


def emit_A_inv(c, t, win, c2, ctx, mP, pool=10, tail=True, join=True):
    """Exact inverse of the rows of emit_A on the window wires alone. mP: classical bit = X-basis outcome of the
    pending wire (measured by the caller once everything added to it after emit_A is taken out again); the part
    of its phase that comes from the rows is repaired here, the caller repairs the initial content.
    join=True: narrow rows are undone jointly (one propagation chain for several rows) as far as `pool` clean
    wires allow."""
    _group_rev(c, win, ctx["specs"], len(win), None, mP, ctx["ghosts"], ctx["qghost"], pool, ctx["tree"], tail,
               join)
    assert not ctx["ghosts"] and not ctx["qghost"], "unrepaired forward outcomes"
