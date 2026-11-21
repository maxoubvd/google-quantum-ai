"""invsqrt_gen: stage emitters (own code, workstream invsqrt_gen). Every gate of the generated circuit comes from
the functions below. Techniques and their sources:
  - temporary AND + measurement erasure, ripple adder with temporary carries: Gidney, arXiv:1709.06648
  - add-subtract by conjugating the accumulator with CNOTs from the control: Litinski, arXiv:2410.00899 (variant)
  - unary-iteration table lookup: Babbush et al., arXiv:1805.03662
  - thermometer flags and conditional shifts for normalisation: companion paper of the competition, App. B
  - phase stream through delay cells that are zero low bits of the normalised input: idea read in the dissection
    of the best public circuit (research/ws_invsqrt.md); the copy / measure / one-hot fix-up form is our own.
"""
from __future__ import annotations

from fractions import Fraction


# ------------------------------------------------------------------ normalisation

def thermo(c, xq):
    """f[k] = [x < 2^k] for k = 1..21 on fresh wires (19 ANDs; x21 and x20 are never both set on the domain)."""
    f = {}
    w = c.alloc()
    c.cx(xq[21], w)
    c.x(w)
    f[21] = w
    w = c.alloc()
    c.cx(xq[21], w)
    c.cx(xq[20], w)
    c.x(w)
    f[20] = w
    for k in range(19, 0, -1):
        w = c.alloc()
        c.x(xq[k])
        c.ccx(f[k + 1], xq[k], w)
        c.x(xq[k])
        f[k] = w
    return f


def thermo_erase(c, xq, f):
    for k in range(1, 20):                      # f_k = f_{k+1} AND NOT x_k
        m = c.hmr(f[k])
        c.cz(f[k + 1], xq[k], cond=m)
        c.z(f[k + 1], cond=m)
        c.bit_free(m)
        c.free(f[k])
    c.x(f[20])
    c.cx(xq[20], f[20])
    c.cx(xq[21], f[20])
    c.free(f[20])
    c.x(f[21])
    c.cx(xq[21], f[21])
    c.free(f[21])


def linear_of_e(c, f, g, w=None):
    """XOR g(e) into wire w (fresh if None), g any 0/1 function of the exponent e in 0..21: CNOTs only, because
    g(e) = g(21) xor sum_k f_k (g(k-1) xor g(k)). Applying it twice restores the wire."""
    if w is None:
        w = c.alloc()
    if g(21):
        c.x(w)
    for k in range(1, 22):
        if g(k - 1) != g(k):
            c.cx(f[k], w)
    return w


def shifter(c, xq, stages, inverse=False):
    """stages: list of (D, control wire, negate). Forward: for each stage whose control is set, x <<= D (the top
    D bits are zero then). 22 - D CSWAPs per stage."""
    n = len(xq)
    for D, ctl, neg in (stages[::-1] if inverse else stages):
        if neg:
            c.x(ctl)
        idx = range(0, n - D) if inverse else range(n - D - 1, -1, -1)
        for i in idx:
            c.cswap(ctl, xq[i], xq[i + D])
        if neg:
            c.x(ctl)


def shifter_refl(c, xq, p, db, inverse=False):
    """x <<= s for s = 2 d + 1 - p (d = db in binary), the top s bits being zero, with 64 CSWAPs instead of 79.
    Idea credited to the best public circuit (read there: its shifter is built from reflections): two
    reflections q -> c - q and q -> c' - q shift by c' - c, and a reflection of an interval costs half as many
    swaps as a shift of it. Our own network keeps binary controls: three reflections with sums c, c + a,
    c + a + b and controls s_a, s_a xor s_b, s_b realise the shifts 0, a, b, a + b [P: the taken reflections are
    the two ends of the run of set bits]. Stages: shift by 16 (6 swaps); sums 17, 21, 29 for the shifts 4 and 8
    (9 + 11 + 7 swaps); sums 20, 21, 23 for the shifts 1 and 2 (10 + 11 + 10 swaps). The data stays inside
    q0..q21 at every step (checked on every input by the simulator)."""
    def refl(csum, lo, hi, ctl):
        for q in range(lo, hi + 1):
            c.cswap(ctl, xq[q], xq[csum - q])

    def st16():
        for i in (range(0, 6) if inverse else range(5, -1, -1)):
            c.cswap(db[3], xq[i], xq[i + 16])

    def r17():
        refl(17, 0, 8, db[1])

    def r21a():
        c.cx(db[1], db[2])
        refl(21, 0, 10, db[2])
        c.cx(db[1], db[2])

    def r29():
        refl(29, 8, 14, db[2])

    def r20():
        c.x(p)
        refl(20, 0, 9, p)
        c.x(p)

    def r21b():
        c.cx(p, db[0])
        c.x(db[0])
        refl(21, 0, 10, db[0])
        c.x(db[0])
        c.cx(p, db[0])

    def r23():
        refl(23, 2, 11, db[0])

    seq = [st16, r17, r21a, r29, r20, r21b, r23]
    for fn in (reversed(seq) if inverse else seq):
        fn()


# ------------------------------------------------------------------ lookup by unary iteration

def qrom(c, addr, data, outs):
    """outs ^= data[address]; addr MSB first; 2^n - 2 ANDs, n - 1 temporary wires."""
    n = len(addr)

    def leaf(ctrl, v):
        for i, w in enumerate(outs):
            if (v >> i) & 1:
                c.cx(ctrl, w)

    def rec(ctrl, level, base):
        if level == n:
            leaf(ctrl, data[base])
            return
        a = addr[level]
        half = 1 << (n - level - 1)
        t = c.alloc()
        c.x(a)
        c.ccx(ctrl, a, t)
        c.x(a)
        rec(t, level + 1, base)
        c.cx(ctrl, t)
        rec(t, level + 1, base + half)
        m = c.hmr(t)
        c.cz(ctrl, a, cond=m)
        c.bit_free(m)
        c.free(t)

    a = addr[0]
    c.x(a)
    rec(a, 1, 0)
    c.x(a)
    rec(a, 1, 1 << (n - 1))


# ------------------------------------------------------------------ adders and rows

def add_into(c, acc, opnd, cin):
    """acc += value(opnd) + cin (mod 2^len(acc)); opnd: operand wires LSB first (may be shorter than acc), cin:
    carry-in wire or None. One AND per stage that has a carry to produce, none for the top stage."""
    n = len(acc)
    car = [None] * n
    car[0] = cin
    for j in range(n - 1):
        b = opnd[j] if j < len(opnd) else None
        ci, a = car[j], acc[j]
        if b is not None and ci is not None:
            c.cx(ci, a)
            c.cx(ci, b)
            w = c.and_new(a, b)
            c.cx(ci, w)
        elif b is not None:
            w = c.and_new(a, b)
        elif ci is not None:
            w = c.and_new(a, ci)
        else:
            w = None
        car[j + 1] = w
    j = n - 1
    if j < len(opnd):
        c.cx(opnd[j], acc[j])
    if car[j] is not None:
        c.cx(car[j], acc[j])
    for j in range(n - 2, -1, -1):
        b = opnd[j] if j < len(opnd) else None
        ci, a, w = car[j], acc[j], car[j + 1]
        if b is not None and ci is not None:
            c.cx(ci, w)
            c.and_del(w, a, b)
            c.cx(ci, b)
            c.cx(b, a)
        elif b is not None:
            c.and_del(w, a, b)
            c.cx(b, a)
        elif ci is not None:
            c.and_del(w, a, ci)
            c.cx(ci, a)


def row(c, acc, src, k, tbit, sub_when):
    """acc -= rr(src, k) if tbit == sub_when else acc += rr(src, k); rr(v, k) = floor(v / 2^k) + bit k-1 of v."""
    if sub_when == 0:
        c.x(tbit)
    for a in acc:
        c.cx(tbit, a)
    add_into(c, acc, src[k:], src[k - 1] if k >= 1 else None)
    for a in acc:
        c.cx(tbit, a)
    if sub_when == 0:
        c.x(tbit)


# ------------------------------------------------------------------ blocked adders (few carry wires)

def _stage_up(c, a, b, ci):
    """Carry out of one adder stage on a fresh wire (None if the stage cannot carry). Leaves a and b XORed with
    ci when both an operand bit and a carry-in are present (undone by _stage_down)."""
    if b is not None and ci is not None:
        c.cx(ci, a)
        c.cx(ci, b)
        w = c.and_new(a, b)
        c.cx(ci, w)
        return w
    if b is not None:
        return c.and_new(a, b)
    if ci is not None:
        return c.and_new(a, ci)
    return None


def _stage_down(c, a, b, ci, w, erase, write):
    """Second half of a stage: erase the carry-out w (free) unless it is kept, restore b, and either leave the
    sum bit in a (write) or restore a."""
    if b is not None and ci is not None:
        if erase:
            c.cx(ci, w)
            c.and_del(w, a, b)
        c.cx(ci, b)
        c.cx(b if write else ci, a)
    elif b is not None:
        if erase:
            c.and_del(w, a, b)
        if write:
            c.cx(b, a)
    elif ci is not None:
        if erase:
            c.and_del(w, a, ci)
        if write:
            c.cx(ci, a)


def and_forms(c, t, A, ca, B, cb):
    """t ^= (XOR of the wires in A, xor ca) AND (XOR of the wires in B, xor cb): one Toffoli between CNOT fans.
    A and B are non-empty sets of wires, A != B, t in neither; the two forms may share wires."""
    A, B = set(A), set(B)
    assert A and B and A != B and t not in A and t not in B
    if A - B:
        X, cx_, Y, cy = A, ca, B, cb
    else:
        X, cx_, Y, cy = B, cb, A, ca
    p1, p2 = min(X - Y), min(Y)
    ops = [("cx", w, p2) for w in sorted(Y - {p2})]
    if cy:
        ops.append(("x", p2))
    rest, cst = set(X - {p1}), cx_
    if p2 in rest:                               # wire p2 now holds the whole form Y
        rest.remove(p2)
        rest ^= (Y - {p2})
        ops.append(("cx", p2, p1))
        cst ^= cy
    ops += [("cx", w, p1) for w in sorted(rest)]
    if cst:
        ops.append(("x", p1))
    for o in ops:
        c.cx(o[1], o[2]) if o[0] == "cx" else c.x(o[1])
    c.ccx(p1, p2, t)
    for o in reversed(ops):
        c.cx(o[1], o[2]) if o[0] == "cx" else c.x(o[1])


def phase_forms(c, A, ca, B, cb, cond):
    """Phase (-1)^((XOR of A, xor ca) AND (XOR of B, xor cb)) under the classical bit cond: Z, CZ and NEG only."""
    for a in sorted(A):
        for b in sorted(B):
            if a == b:
                c.z(a, cond=cond)
            else:
                c.cz(a, b, cond=cond)
    if ca:
        for b in sorted(B):
            c.z(b, cond=cond)
    if cb:
        for a in sorted(A):
            c.z(a, cond=cond)
    if ca and cb:
        c.neg(cond=cond)


def add_blocked(c, acc, opnd, cin, blk, mode, pend=None, top_fix=None):
    """acc += value(opnd) + cin (mod 2^len(acc)) with at most blk + (number of blocks - 1) carry wires.
    The stages are cut into blocks of blk; inside a block the carries are erased for free as in add_into. The
    carry between two blocks cannot be erased for free once the lower block is finished, so:
      mode 'fwd': it is measured as soon as the upper block is done and the outcome bit is returned in a dict
                  {carry index: bit}. The phase is repaired for free in the mirrored row (the row that later
                  undoes this one has the same carry chain [P]), which must be given the dict as `pend`.
      mode 'rev': boundary carries stay alive to the end of the row, then each is measured and, only if the
                  outcome is 1, the block below is recomputed on the complemented sum bits to apply the phase of
                  the carry (blk - 1 Toffolis under the condition, mean (blk - 1) / 2).
    blk may be a list of block sizes (see below).
    top_fix: list of classical bits. For each, the phase (-1)^(bit AND carry out of the top stage) is applied by
    a conditioned CZ between the carry into the top stage and the top accumulator wire (the top stage must have
    no operand bit). This is how a pending wire that took the top carries of a forward row is repaired after it
    was measured: the mirrored row has the same carries [P], and for a stage without operand bit the carry out
    is (accumulator bit before the sum) AND (carry in) in both rows."""
    n = len(acc)
    ob = lambda j: opnd[j] if j < len(opnd) else None
    car = [None] * (n + 1)
    car[0] = cin
    out = {}
    bounds = []
    if isinstance(blk, int):
        starts = list(range(0, n, blk))
    else:
        # list of block sizes from the bottom, the last one repeated. In mode 'rev' the wires in use are
        # (boundaries kept so far) + (carries of the block being built), and the mean repair is
        # (stages below the top block - number of blocks below it) / 2: shrinking blocks trade wires for repairs.
        starts, j = [], 0
        while j < n:
            starts.append(j)
            j += blk[min(len(starts) - 1, len(blk) - 1)]
    ends = starts[1:] + [n]
    for bi, j0 in enumerate(starts):
        j1 = ends[bi] - 1
        for j in range(j0, j1 + 1):
            if j == n - 1:
                break
            car[j + 1] = _stage_up(c, acc[j], ob(j), car[j])
            if pend and (j + 1) in pend:
                c.z(car[j + 1], cond=pend[j + 1])
                c.bit_free(pend[j + 1])
        for j in range(j1, j0 - 1, -1):
            if j == n - 1:
                if top_fix:
                    assert ob(j) is None and car[j] is not None
                    for bfix in top_fix:
                        c.cz(car[j], acc[j], cond=bfix)
                if ob(j) is not None:
                    c.cx(ob(j), acc[j])
                if car[j] is not None:
                    c.cx(car[j], acc[j])
                continue
            _stage_down(c, acc[j], ob(j), car[j], car[j + 1], erase=(j != j1), write=True)
        if bi > 0 and car[j0] is not None:
            if mode == "fwd":
                out[j0] = c.hmr(car[j0])          # phase repaired later by the mirrored row
                c.free(car[j0])
            else:
                bounds.append((j0, starts[bi - 1]))
    for J, jb in reversed(bounds):
        m = c.hmr(car[J])
        c.free(car[J])
        c.push(m)
        for a in acc[jb:J]:
            c.x(a)
        tmp = {jb: car[jb]}
        for j in range(jb, J - 1):
            tmp[j + 1] = _stage_up(c, acc[j], ob(j), tmp[j])
        a, b, ci = acc[J - 1], ob(J - 1), tmp[J - 1]
        if b is not None and ci is not None:
            c.cz(a, b)
            c.cz(a, ci)
            c.cz(b, ci)
        elif b is not None:
            c.cz(a, b)
        elif ci is not None:
            c.cz(a, ci)
        for j in range(J - 2, jb - 1, -1):
            _stage_down(c, acc[j], ob(j), tmp[j], tmp[j + 1], erase=True, write=False)
        for a in acc[jb:J]:
            c.x(a)
        c.pop()
        c.bit_free(m)
    return out


def row_blocked(c, acc, src, k, tbit, sub_when, blk, mode, pend=None):
    """Same arithmetic as row(), with add_blocked."""
    if sub_when == 0:
        c.x(tbit)
    for a in acc:
        c.cx(tbit, a)
    out = add_blocked(c, acc, src[k:], src[k - 1] if k >= 1 else None, blk, mode, pend)
    for a in acc:
        c.cx(tbit, a)
    if sub_when == 0:
        c.x(tbit)
    return out


def row_virt(c, acc, u, k, tbit, sub_when, blk, mode, pend=None):
    """Row on an accumulator whose lowest bit is virtual (own device, see research/ws_invsqrt_gen.md).
    Arithmetic: acc21 -+= rr(u18, k), where acc21 = (acc, g0) and u18 = (u, 0): u holds bits 1..17 of the slope
    register (its bit 0 is |0> after A) and acc holds bits 1.. of lo. With C0 even, g0 equals u18 bit k after the
    row with shift k (rows in increasing k), so it needs no wire; the carry out of it into acc[0] is
    u18_(k-1) xor t (u18_(k-1) xor u18_k) for forward and reverse rows alike (one Toffoli)."""
    c1w = None
    if k >= 1:
        b = u[k - 1]                       # u18 bit k
        c1w = c.alloc()
        if k >= 2:
            a = u[k - 2]                   # u18 bit k-1
            c.cx(a, b)
            c.ccx(tbit, b, c1w)
            c.cx(a, b)
            c.cx(a, c1w)
        else:
            c.ccx(tbit, b, c1w)
    if sub_when == 0:
        c.x(tbit)
    for w in acc:
        c.cx(tbit, w)
    out = add_blocked(c, acc, u[k:], c1w, blk, mode, pend)
    for w in acc:
        c.cx(tbit, w)
    if sub_when == 0:
        c.x(tbit)
    if c1w is not None:
        m = c.hmr(c1w)
        if k >= 2:
            c.z(a, cond=m)
            c.cz(tbit, a, cond=m)
        c.cz(tbit, b, cond=m)
        c.bit_free(m)
        c.free(c1w)
    return out


def row_win(c, acc, c2, k, kg, tbit, sub_when, blk, mode, pend=None, top_fix=None):
    """Stage A row on the low window of u (own construction; the public file's stage A has this shape, read in
    research/ws_invsqrt_decode.md 2.1): acc -+= rr(c2, k) on the window, whose bit 0 is virtual.
    acc holds bits 1.. of the 18-bit u (LSB first; the last wire may be a pending wire that takes the top carry).
    Bit 0 of u is c2[kg] before the first row (the tie of the model) and changes by c2[k] xor c2[k-1] in the
    row with shift k, so before that row it equals c2[kg] xor c2[k-1] and needs no wire [P]. The carry out of it
    into acc[0] is the majority of (that bit in the row's frame, c2[k], c2[k-1]):
        c2[k-1] xor (c2[kg] xor t) AND (c2[k] xor c2[k-1])            (one Toffoli)
    with t the digit wire: the same expression in the forward row (subtract when t = 1) and in its mirror,
    because the mirrored row has the same carries [P]. Rows must run with k increasing forward, decreasing in
    the mirror. sub_when: value of t for which the row subtracts."""
    A = {c2[kg], tbit}
    B = {c2[k], c2[k - 1]} if k >= 1 else {c2[0]}
    c1w = c.alloc()
    and_forms(c, c1w, A, 0, B, 0)
    if k >= 1:
        c.cx(c2[k - 1], c1w)
    if sub_when == 0:
        c.x(tbit)
    for w in acc:
        c.cx(tbit, w)
    out = add_blocked(c, acc, c2[k + 1:], c1w, blk, mode, pend, top_fix)
    for w in acc:
        c.cx(tbit, w)
    if sub_when == 0:
        c.x(tbit)
    if k >= 1:
        c.cx(c2[k - 1], c1w)
    m = c.hmr(c1w)
    phase_forms(c, A, 0, B, 0, m)
    c.bit_free(m)
    c.free(c1w)
    return out


# ------------------------------------------------------------------ monomials of the delay bits

def monomial_basis(nbits, valid):
    """Monomials (bit masks over the delay bits, low degree first) that are independent as functions on the
    valid delay values, and for every valid d' the mask of chosen monomials whose XOR is the indicator [d = d']."""
    valid = list(valid)
    cands = sorted(range(1 << nbits), key=lambda m: (bin(m).count("1"), m))
    lead, chosen, vecs = {}, [], []
    for mk in cands:
        vec = 0
        for idx, d in enumerate(valid):
            if d & mk == mk:
                vec |= 1 << idx
        v = vec
        while v:
            h = v.bit_length() - 1
            if h in lead:
                v ^= lead[h]
            else:
                lead[h] = v
                chosen.append(mk)
                vecs.append(vec)
                break
        if len(chosen) == len(valid):
            break
    assert len(chosen) == len(valid)
    rows = [(v, 1 << i) for i, v in enumerate(vecs)]
    used, piv = set(), {}
    for bit in range(len(valid)):
        r = next(k for k in range(len(rows)) if k not in used and (rows[k][0] >> bit) & 1)
        used.add(r)
        piv[bit] = r
        for k in range(len(rows)):
            if k != r and (rows[k][0] >> bit) & 1:
                rows[k] = (rows[k][0] ^ rows[r][0], rows[k][1] ^ rows[r][1])
    ind = {d: rows[piv[idx]][1] for idx, d in enumerate(valid)}
    return chosen, ind


class Monos:
    """Wires holding the products of delay bits named by `masks` (one AND per product of two or more bits)."""

    def __init__(self, c, ctrl, masks):
        self.c, self.ctrl, self.w, self.made = c, ctrl, {}, []
        for mk in sorted(masks, key=lambda m: (bin(m).count("1"), m)):
            self.get(mk)

    def get(self, mk):
        if mk in self.w:
            return self.w[mk]
        bits = [b for b in range(len(self.ctrl)) if (mk >> b) & 1]
        if len(bits) == 1:
            w = self.ctrl[bits[0]]
        else:
            sub = self.get(mk & ~(1 << bits[-1]))
            w = self.c.and_new(sub, self.ctrl[bits[-1]])
            self.made.append((w, sub, self.ctrl[bits[-1]]))
        self.w[mk] = w
        return w

    def erase(self):
        for w, a, b in reversed(self.made):
            self.c.and_del(w, a, b)
        self.made = []


# ------------------------------------------------------------------ phase stream

def stream_plan(n_src, tiers, delays):
    """tiers: list of (D, on) with on(d) -> bool, in cascade order. Symbolic run for every delay d.
    Returns (steps, need[tau][tier], out[d][tau]) and checks out[d][tau] == tau - d and that all cells end empty."""
    steps = n_src + max(delays)
    need = [[False] * len(tiers) for _ in range(steps)]
    out = {}
    for d in delays:
        cells = [[None] * D for D, _ in tiers]
        o = []
        for tau in range(steps):
            r = tau if tau < n_src else None
            for ti, (D, on) in enumerate(tiers):
                if not on(d):
                    continue
                cell = cells[ti][tau % D]
                if r is not None or cell is not None:
                    need[tau][ti] = True
                cells[ti][tau % D], r = r, cell
            o.append(r)
            want = tau - d if 0 <= tau - d < n_src else None
            assert r == want, (d, tau, r, want)
        assert all(x is None for cs in cells for x in cs), (d, cells)
        out[d] = o
    return steps, need, out


def rotate(c, q, pos, angles=None):
    """Z_POW of weight 2^-pos half turns (pos 0, 1, 2 are Z, S, T)."""
    if angles is not None and pos in angles:
        c.zpow(q, angles[pos])
    else:
        c.zpow(q, Fraction(1, 1 << pos))


def stream(c, src, tiers, ctrl, cells, delays, onehot, pos0, top_pos, angles=None):
    """Deliver bit k of src at position pos0 - (k + d) for the delay d encoded by the tier controls.
    tiers: list of (D, on(d)); ctrl[i]: control wire of tier i; cells[i]: D wires that are |0> whenever tier i is on;
    onehot[d]: wire that is 1 iff the delay is d (used only as CZ partner in the fix-ups)."""
    steps, need, out = stream_plan(len(src), tiers, delays)
    r = c.alloc()
    for tau in range(steps):
        if tau < len(src):
            c.cx(src[tau], r)
        for ti, (D, _) in enumerate(tiers):
            if need[tau][ti]:
                c.cswap(ctrl[ti], r, cells[ti][tau % D])
        pos = pos0 - tau
        if 0 <= pos <= top_pos:
            rotate(c, r, pos, angles)
        m = c.hmr(r)
        for d in delays:
            k = out[d][tau]
            if k is not None:
                c.cz(onehot[d], src[k], cond=m)
        c.bit_free(m)
    c.free(r)


def stream_m(c, src, tier_bits, ctrl, cells, delays, pos0, top_pos, angles=None):
    """Same stream as stream(), with binary delay wires only: tier i has delay 2^tier_bits[i] and is on iff that
    bit of d is set; ctrl[b] is the wire of delay bit b. The measured copy at step tau equals
    sum over d' of [d = d'] src[tau - d'], and every indicator is a XOR of monomials of the delay bits, so the
    fix-up is a set of conditioned CZs between monomial wires and source wires (and conditioned Zs)."""
    tiers = [(1 << b, (lambda d, b=b: (d >> b) & 1)) for b in tier_bits]
    steps, need, out = stream_plan(len(src), tiers, delays)
    monos, ind = monomial_basis(len(ctrl), delays)
    mw = Monos(c, ctrl, [mk for mk in monos if mk])
    r = c.alloc()
    for tau in range(steps):
        if tau < len(src):
            c.cx(src[tau], r)
        for ti, (D, _) in enumerate(tiers):
            if need[tau][ti]:
                c.cswap(ctrl[tier_bits[ti]], r, cells[ti][tau % D])
        pos = pos0 - tau
        if 0 <= pos <= top_pos:
            rotate(c, r, pos, angles)
        m = c.hmr(r)
        for i, mk in enumerate(monos):
            for k in range(len(src)):
                par = 0
                for d in delays:
                    if out[d][tau] == k and (ind[d] >> i) & 1:
                        par ^= 1
                if par:
                    if mk:
                        c.cz(mw.get(mk), src[k], cond=m)
                    else:
                        c.z(src[k], cond=m)
        c.bit_free(m)
    c.free(r)
    mw.erase()


def stream_d(c, src, tier_bits, ctrl, cells, delays, pos0, top_pos, angles=None):
    """stream_m with the product fix-ups deferred to the end of the stream (the public file repairs its streams
    after the fact as well, research/ws_invsqrt_decode.md 3.5; own construction): the outcome bit of every step is
    kept, fix-ups that need no product wire (Z, CZ with a control wire) are applied at once, and each product
    of control bits is then put on a wire once, used for all its conditioned CZs and erased again. The source
    wires and the control wires do not change in between, so the deferred phase is the same [P].
    Wires in use: the rotation wire during the stream, at most two product wires after it (stream_m: six product
    wires next to the rotation wire). Same Toffoli count as stream_m.
    An entry None in src is a source bit that is 1 on every input (no wire): an X instead of the CNOT, and its
    fix-ups lose one factor (CZ -> Z, Z -> NEG)."""
    tiers = [(1 << b, (lambda d, b=b: (d >> b) & 1)) for b in tier_bits]
    steps, need, out = stream_plan(len(src), tiers, delays)
    monos, ind = monomial_basis(len(ctrl), delays)
    deg = [bin(mk).count("1") for mk in monos]
    fix = {i: [] for i in range(len(monos))}          # product index -> [(outcome bit, source wire)]
    kept = []

    def cz_src(w, s, m):
        if s is None:
            c.z(w, cond=m)
        else:
            c.cz(w, s, cond=m)

    r = c.alloc()
    for tau in range(steps):
        if tau < len(src):
            if src[tau] is None:
                c.x(r)
            else:
                c.cx(src[tau], r)
        for ti, (D, _) in enumerate(tiers):
            if need[tau][ti]:
                c.cswap(ctrl[tier_bits[ti]], r, cells[ti][tau % D])
        pos = pos0 - tau
        if 0 <= pos <= top_pos:
            rotate(c, r, pos, angles)
        m = c.hmr(r)
        kept.append(m)
        for i, mk in enumerate(monos):
            for k in range(len(src)):
                par = 0
                for d in delays:
                    if out[d][tau] == k and (ind[d] >> i) & 1:
                        par ^= 1
                if not par:
                    continue
                if deg[i] == 0:
                    if src[k] is None:
                        c.neg(cond=m)
                    else:
                        c.z(src[k], cond=m)
                elif deg[i] == 1:
                    cz_src(ctrl[mk.bit_length() - 1], src[k], m)
                else:
                    fix[i].append((m, src[k]))
    c.free(r)
    todo = sorted((i for i in range(len(monos)) if deg[i] >= 2), key=lambda i: (deg[i], monos[i]))
    assert all(deg[i] <= 3 for i in todo)
    done = set()
    for i in todo:
        if i in done or deg[i] != 2:
            continue
        b0, b1 = [b for b in range(len(ctrl)) if (monos[i] >> b) & 1]
        w = c.and_new(ctrl[b0], ctrl[b1])
        for m, s in fix[i]:
            cz_src(w, s, m)
        done.add(i)
        for i2 in todo:                               # products of three bits that contain this pair
            if i2 in done or deg[i2] != 3 or monos[i2] & monos[i] != monos[i]:
                continue
            b2 = (monos[i2] ^ monos[i]).bit_length() - 1
            w2 = c.and_new(w, ctrl[b2])
            for m, s in fix[i2]:
                cz_src(w2, s, m)
            c.and_del(w2, w, ctrl[b2])
            done.add(i2)
        c.and_del(w, ctrl[b0], ctrl[b1])
    for i in todo:                                    # a product of three bits none of whose pairs is in the basis
        if i in done:
            continue
        b0, b1, b2 = [b for b in range(len(ctrl)) if (monos[i] >> b) & 1]
        w = c.and_new(ctrl[b0], ctrl[b1])
        w2 = c.and_new(w, ctrl[b2])
        for m, s in fix[i]:
            cz_src(w2, s, m)
        c.and_del(w2, w, ctrl[b2])
        c.and_del(w, ctrl[b0], ctrl[b1])
    for m in kept:
        c.bit_free(m)
