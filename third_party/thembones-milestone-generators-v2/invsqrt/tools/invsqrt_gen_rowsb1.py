"""invsqrt_rows_b: add-subtract rows of stage B (lo -+= u >> k) and their exact inverse, for our own generator.

Own code of workstream invsqrt_rows, written against tools/invsqrt_gen_circ.py (builder + exact simulator).
Method credits (public design, implemented here from its description in research/ws_invsqrt_decode.md; no gate
list of any public circuit was read into this file):
  - truncated row windows in groups with one pending wire per group, pending wire measured before the reverse rows
    and repaired by one CZ per row, forward carries measured without fix-up and repaired by the mirrored row:
    manwithacat "HW-25" (join of gengsr and maoxin123), best public invsqrt circuit of the competition.
  - temporary-AND adder, measurement-based uncomputation: Gidney, arXiv:1709.06648.
  - add-subtract by conjugating the accumulator with CNOTs from the digit: Litinski, arXiv:2410.00899 (variant).
Own changes against that design (research/ws_invsqrt_rows.md):
  - the pending wire starts as a copy of the first digit of its group (not constant 1), so its merge is a
    one-bit increment in the frame of that digit: no overflow AND, no constant left in lo;
  - merges are nested (pending wire of group g+1 is merged into the window of group g, its carry goes into the
    pending wire of group g): 10 ANDs instead of 14, every merge temporary is erased for free in both passes;
  - forward rows keep two carry wires only (every forward carry is measured at once and repaired by the reverse row).

Arithmetic (exact, for every u, t, lo; see rows_model): row with shift k and digit wire t_k does, mod 2^len(lo),
    t_k = 0:  lo += (u >> k) + u[k-2]          t_k = 1:  lo -= (u >> k) + u[k-1]          (u[-1] = u[-2] = 0)
which is the generator's row_virt arithmetic (invsqrt_gen_parts.py) on the same wires.
"""
from __future__ import annotations

import numpy as np

DEFAULT_GROUPS = ((0, 1), (2, 3, 4), (5, 6, 7), (8, 9, 10, 11, 12, 13))   # measured best of 10 partitions
PUBLIC_GROUPS = ((0, 1, 2), (3, 4, 5), (6, 7, 8), (9, 10, 11, 12, 13))     # partition of the public file


# ------------------------------------------------------------------ reference model

def rows_model(u, t, lo, nu=17, nt=14, wlo=20):
    """Bit-exact content of the lo register after emit_B, for integer arrays u (nu bits), t (nt digits, digit i is
    bit i; its row has shift k = nt - 1 - i) and lo (wlo bits). The truncation into windows changes nothing:
    the result equals full-width rows."""
    u = np.asarray(u, dtype=np.int64)
    t = np.asarray(t, dtype=np.int64)
    acc = np.asarray(lo, dtype=np.int64).copy()
    for k in range(nt):
        i = nt - 1 - k
        add = (u >> k) + (((u >> (k - 2)) & 1) if k >= 2 else 0)
        sub = (u >> k) + (((u >> (k - 1)) & 1) if k >= 1 else 0)
        acc = np.where((t >> i) & 1, acc - sub, acc + add)
    return acc & ((1 << wlo) - 1)


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


# ------------------------------------------------------------------ carry-in of a row (one AND)

def _cin_make(c, u, k, t):
    """cin = u[k-2] xor t (u[k-1] xor u[k-2]) on a fresh wire (None for k = 0)."""
    if k == 0:
        return None
    w = c.alloc()
    if k >= 2:
        c.cx(u[k - 2], u[k - 1])
    c.ccx(t, u[k - 1], w)
    if k >= 2:
        c.cx(u[k - 2], u[k - 1])
        c.cx(u[k - 2], w)
    return w


def _cin_erase(c, u, k, t, w):
    if k >= 2:
        c.cx(u[k - 2], w)
        c.cx(u[k - 2], u[k - 1])
    m = c.hmr(w)
    c.cz(t, u[k - 1], cond=m)
    c.bit_free(m)
    c.free(w)
    if k >= 2:
        c.cx(u[k - 2], u[k - 1])


# ------------------------------------------------------------------ forward row

def _row_fwd(c, win, top, u, k, t, ghosts):
    """Row with shift k in the frame of digit wire t on the window `win` (lo wires, LSB first); the carry out of
    win[-1] is XORed into `top`. Every carry is measured as soon as the next one exists; the outcome bit of the
    carry into bit j is stored in ghosts[(k, j)] and repaired by _row_rev."""
    nu = len(u)
    for a in win:
        c.cx(t, a)
    car = _cin_make(c, u, k, t)
    is_cin = True
    n = len(win)
    for j in range(n):
        a = win[j]
        b = u[k + j] if k + j < nu else None
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
                _cin_erase(c, u, k, t, car)
            else:
                ghosts[(k, j)] = c.hmr(car)
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


# ------------------------------------------------------------------ forward stage

def _layout(groups, nu, wlo):
    """Per group: (shifts, window width in lo wires, has pending wire)."""
    out = []
    for gi, grp in enumerate(groups):
        if gi == 0:
            out.append((tuple(grp), wlo - 1, False))
        else:
            out.append((tuple(grp), nu + 1 - grp[0], True))
    return out


def emit_B(c, t, u, lo, groups=DEFAULT_GROUPS, tree=True):
    """Forward rows and merge. t: digit wires (t[i] = digit i, row shift k = len(t) - 1 - i), u: operand wires LSB
    first, lo: accumulator wires LSB first. Returns the context for emit_B_inv; ctx['P'] are the pending wires that
    stay allocated (and non-zero) until emit_B_inv.
    tree=True (own device): inside a group every row after the first gets its own, narrower window (its operand
    bits plus one) with a temporary pending wire on top; after the group's rows these wires are merged upwards
    one bit each (lowest first) and measured at once. Their outcomes are repaired by the reverse rows, which
    recompute the same carries on the group's full window. A row then costs its operand bits + 2 instead of
    operand bits + 1 + its distance to the group top."""
    nt, nu, wlo = len(t), len(u), len(lo)
    assert sorted(k for g in groups for k in g) == list(range(nt))
    lay = _layout(groups, nu, wlo)
    ghosts = {}
    pend = {}
    qghost = {}                              # shift k -> (outcome bit of its temporary wire, bit position)
    for gi, (grp, w, has_p) in enumerate(lay):
        nodes = []                           # (k, digit wire, node bit position, node wire, is temporary)
        for idx, k in enumerate(grp):
            tk = t[nt - 1 - k]
            if has_p and idx == 0:
                P = c.alloc()
                c.cx(tk, P)
                pend[gi] = (P, tk, w)
                nodes.append((k, tk, w, P, False))
            elif tree and nu + 1 - k < w + (0 if has_p else 1):
                q = c.alloc()
                c.cx(tk, q)
                nodes.append((k, tk, nu + 1 - k, q, True))
            else:
                nodes.append((k, tk, w, pend[gi][0] if has_p else lo[w], False))
            _, _, a, wire, _ = nodes[-1]
            _row_fwd(c, lo[:a], wire, u, k, tk, ghosts)
        for idx in range(len(nodes) - 1, -1, -1):
            k, tk, a, wire, temp = nodes[idx]
            if not temp:
                continue
            up = [nd for nd in nodes[:idx] if nd[4] or has_p]
            if up:
                _, _, a_up, w_up, _ = up[-1]
            else:
                a_up, w_up = w, lo[w]
            _merge(c, wire, tk, lo[a:a_up], w_up, False)
            qghost[k] = (c.hmr(wire), a)
            c.free(wire)
    for gi in range(len(lay) - 1, 0, -1):
        P, tf, w = pend[gi]
        if gi == 1:
            _merge(c, P, tf, lo[w:wlo - 1], lo[wlo - 1], False)
        else:
            Pup, _, wup = pend[gi - 1]
            _merge(c, P, tf, lo[w:wup], Pup, False)
    return dict(ghosts=ghosts, pend=pend, lay=lay, qghost=qghost, P=[pend[g][0] for g in sorted(pend)])


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


def _row_rev(c, win, top, u, k, t, ghosts, mP, pool, zfix=None):
    """Inverse of _row_fwd on the same window. top: wire that receives the carry out of win[-1] (first group) or
    None (windowed groups: that carry only went into the pending wire, whose measured outcome mP is repaired here
    by the phase of the same carry). Repairs the forward outcomes ghosts[(k, j)] and, for zfix = {j: [bits]},
    the outcomes of the forward temporary pending wires that sat at bit j (Z on the carry into bit j).
    At most `pool` carry wires."""
    nu = len(u)
    n = len(win)
    for a in win:
        c.x(a)
        c.cx(t, a)
    car = [None] * n                     # car[j]: wire with the carry into bit j (car[0] = carry-in of the row)
    written = [False] * n
    opb = lambda j: u[k + j] if k + j < nu else None
    j_min = 0 if k >= 1 else 1
    nodes = list(range(j_min, n))

    def compute(j):
        if j == 0:
            car[0] = _cin_make(c, u, k, t)
            return
        a, b, ci = win[j - 1], opb(j - 1), car[j - 1]
        w = c.alloc()
        _enter(c, a, b, ci, written[j - 1])
        _carry_into(c, a, b, ci, w)
        _exit(c, a, b, ci, written[j - 1])
        car[j] = w

    def erase(j):
        if j == 0:
            _cin_erase(c, u, k, t, car[0])
        else:
            a, b, ci = win[j - 1], opb(j - 1), car[j - 1]
            _enter(c, a, b, ci, written[j - 1])
            _carry_erase(c, a, b, ci, car[j])
            _exit(c, a, b, ci, written[j - 1])
        car[j] = None

    def write(j):
        _write(c, win[j], opb(j), car[j])
        written[j] = True

    sizes = plan_blocks(len(nodes), pool)
    blocks, pos = [], 0
    for s in sizes:
        blocks.append(nodes[pos:pos + s])
        pos += s
    for bi, blk in enumerate(blocks):
        last = bi == len(blocks) - 1
        for j in blk:
            compute(j)
            if (k, j) in ghosts:
                c.z(car[j], cond=ghosts[(k, j)])
                c.bit_free(ghosts.pop((k, j)))
            for bq in (zfix or {}).get(j, ()):
                c.z(car[j], cond=bq)
        if last:
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


def emit_B_inv(c, t, u, lo, ctx, pool=10):
    """Exact inverse of emit_B: un-merge, pending wires measured (repaired inside the rows), rows in mirror order.
    `pool`: number of clean wires this stage may hold at one time (all taken with c.alloc(), all released)."""
    nt, nu, wlo = len(t), len(u), len(lo)
    lay, pend, ghosts = ctx["lay"], ctx["pend"], ctx["ghosts"]
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
    qghost = ctx.get("qghost", {})
    for k, (bq, a) in qghost.items():
        c.z(t[nt - 1 - k], cond=bq)
    for gi in range(len(lay) - 1, -1, -1):
        grp, w, has_p = lay[gi]
        for k in reversed(grp):
            zfix = {}
            for k2 in grp:
                if k2 <= k and k2 in qghost:
                    zfix.setdefault(qghost[k2][1], []).append(qghost[k2][0])
            if has_p:
                _row_rev(c, lo[:w], None, u, k, t[nt - 1 - k], ghosts, mP[gi], pool, zfix)
            else:
                _row_rev(c, lo[:w], lo[w], u, k, t[nt - 1 - k], ghosts, None, pool, zfix)
        if has_p:
            c.bit_free(mP[gi])
        for k in grp:
            if k in qghost:
                c.bit_free(qghost.pop(k)[0])
    assert not ghosts, "unrepaired forward outcomes"
