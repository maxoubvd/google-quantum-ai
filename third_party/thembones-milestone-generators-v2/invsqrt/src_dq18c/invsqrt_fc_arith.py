#!/usr/bin/env python3
"""invsqrt_fc_arith.py -- FC arithmetic on the kit's Circ (own code).

Operands are lists, LSB first. A source / carry-in entry is a wire id (int >= 0), None (constant 0) or ONE
(constant 1). Targets are wire ids.
  add(c, src, tgt, cin)      tgt += src + cin  (mod 2^len(tgt)); Gidney temporary-AND carries when the wire cap
                             allows (one AND per carry that is not a constant / alias, erased for free), else the
                             Cuccaro ripple (2 CCX per position, one ancilla, constants materialised)
  fadd(c, ctl, src, tgt, cin) tgt += src + cin if ctl == 1, tgt -= src + cin if ctl == 0 (frame by ctl)
Every routine is self-contained (no measurement debt left), so Circ.mirror() inverts it.
"""
ONE = "1"


def _is_w(x):
    return isinstance(x, int)


def _maj_compute(c, a, b, cc):
    """carry MAJ(a, b, cc) for target bit b; returns (carry, kind) with kind 'w' (fresh wire, AND used),
    'alias' (an existing wire), 'c0' / 'c1' (constant)."""
    if not _is_w(a) and not _is_w(cc):
        av = 1 if a == ONE else 0
        cv = 1 if cc == ONE else 0
        if av == cv:
            return (ONE if av else None), ("c1" if av else "c0")
        return b, "alias"                      # MAJ(1, b, 0) = MAJ(0, b, 1) = b
    if not _is_w(a):                           # a constant, cc wire
        t = c.alloc(1)[0]
        c.and_compute(b, cc, t)
        if a == ONE:
            c.cx(b, t)
            c.cx(cc, t)
        return t, "w"
    if not _is_w(cc):                          # a wire, cc constant
        t = c.alloc(1)[0]
        c.and_compute(a, b, t)
        if cc == ONE:
            c.cx(a, t)
            c.cx(b, t)
        return t, "w"
    t = c.alloc(1)[0]
    c.cx(cc, a)
    c.cx(cc, b)
    c.and_compute(a, b, t)
    c.cx(cc, t)
    c.cx(cc, a)
    c.cx(cc, b)
    return t, "w"


def _maj_erase(c, a, b, cc, t, kind):
    if kind != "w":
        return
    if not _is_w(a):
        if a == ONE:
            c.cx(cc, t)
            c.cx(b, t)
        c.and_erase(b, cc, t)
        return
    if not _is_w(cc):
        if cc == ONE:
            c.cx(b, t)
            c.cx(a, t)
        c.and_erase(a, b, t)
        return
    c.cx(cc, b)
    c.cx(cc, a)
    c.cx(cc, t)
    c.and_erase(a, b, t)
    c.cx(cc, b)
    c.cx(cc, a)


def _sum_into(c, a, cc, b):
    if _is_w(a):
        c.cx(a, b)
    elif a == ONE:
        c.x(b)
    if _is_w(cc):
        c.cx(cc, b)
    elif cc == ONE:
        c.x(b)


def gidney_need(src, tgt, cin):
    """upper bound of fresh carry wires (n - 1)."""
    return max(0, len(tgt) - 1)


def gidney_fresh(src, tgt, cin):
    """exact number of fresh carry wires add_gidney holds at its peak (constant / alias carries need none)."""
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    cur_const = not _is_w(cin)
    curv = 1 if cin == ONE else 0
    k = 0
    for i in range(n - 1):
        a = src[i]
        if cur_const and not _is_w(a):
            av = 1 if a == ONE else 0
            if av == curv:
                continue
            cur_const = False
            continue
        cur_const = False
        k += 1
    return k


def add_gidney(c, src, tgt, cin=None):
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    carries = [(cin, "in")]
    for i in range(n - 1):
        cc = carries[-1][0]
        t, kind = _maj_compute(c, src[i], tgt[i], cc)
        carries.append((t, kind))
    for i in range(n - 1, -1, -1):
        _sum_into(c, src[i], carries[i][0], tgt[i])
        if i >= 1:
            t, kind = carries[i]
            _maj_erase(c, src[i - 1], tgt[i - 1], carries[i - 1][0], t, kind)
    return


def add_cuccaro(c, src, tgt, cin=None):
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    if n == 0:
        return
    mat = []
    A = []
    used = set(tgt) | ({cin} if _is_w(cin) else set())
    for a in src:
        if _is_w(a) and a not in used:
            A.append(a)
            used.add(a)
        elif _is_w(a):
            w = c.alloc(1)[0]
            c.cx(a, w)
            mat.append((w, ("copy", a)))
            A.append(w)
        else:
            w = c.alloc(1)[0]
            if a == ONE:
                c.x(w)
            mat.append((w, a))
            A.append(w)
    anc = None
    if _is_w(cin):
        c0 = cin
    else:
        anc = c.alloc(1)[0]
        if cin == ONE:
            c.x(anc)
        c0 = anc
    if n == 1:
        c.cx(A[0], tgt[0])
        c.cx(c0, tgt[0])
    else:
        def maj(x, y, z):
            c.cx(z, y)
            c.cx(z, x)
            c.ccx(x, y, z)

        def uma(x, y, z):
            c.ccx(x, y, z)
            c.cx(z, x)
            c.cx(x, y)

        chain = [c0] + A
        for i in range(n - 1):
            maj(chain[i], tgt[i], A[i])
        c.cx(A[n - 1], tgt[n - 1])
        c.cx(A[n - 2], tgt[n - 1])
        for i in range(n - 2, -1, -1):
            uma(chain[i], tgt[i], A[i])
    if anc is not None:
        if cin == ONE:
            c.x(anc)
        c.free_([anc])
    for w, a in reversed(mat):
        if a == ONE:
            c.x(w)
        elif isinstance(a, tuple):
            c.cx(a[1], w)
        c.free_([w])


def _gidney_seg(c, src, tgt, cin, keep_out):
    """Gidney add on one segment; if keep_out, also compute the carry-out into a fresh wire and keep it."""
    n = len(tgt)
    carries = [(cin, "in")]
    upto = n if keep_out else n - 1
    for i in range(upto):
        cc = carries[-1][0]
        t, kind = _maj_compute(c, src[i], tgt[i], cc)
        carries.append((t, kind))
    out = carries[n] if keep_out else None
    if keep_out and out[1] == "alias":
        t = c.alloc(1)[0]
        c.cx(out[0], t)
        out = (t, "w")
    elif keep_out and out[1] == "c1":
        t = c.alloc(1)[0]
        c.x(t)
        out = (t, "w")
    for i in range(n - 1, -1, -1):
        _sum_into(c, src[i], carries[i][0], tgt[i])
        if i >= 1:
            t, kind = carries[i]
            _maj_erase(c, src[i - 1], tgt[i - 1], carries[i - 1][0], t, kind)
    return out


def _erase_boundary(c, src, sums, cin, out):
    """out = carry-out of the segment = [sums < src + cin] = carry-out of (NOT sums + src + cin): recompute
    that carry chain with temporary ANDs, XOR it into out, uncompute."""
    t_out, kind_out = out
    if kind_out != "w":
        return
    n = len(sums)
    for w in sums:
        c.x(w)
    carries = [(cin, "in")]
    for i in range(n):
        t, kind = _maj_compute(c, src[i], sums[i], carries[-1][0])
        carries.append((t, kind))
    fin = carries[n][0]
    if isinstance(fin, int):
        c.cx(fin, t_out)
    elif fin == ONE:
        c.x(t_out)
    for i in range(n, 0, -1):
        t, kind = carries[i]
        _maj_erase(c, src[i - 1], sums[i - 1], carries[i - 1][0], t, kind)
    for w in sums:
        c.x(w)
    c.free_([t_out])


def add_segmented(c, src, tgt, cin, M):
    """ripple add in segments of M positions (Gidney inside a segment, boundary carries kept, then erased from
    the top down by a recomputed comparison): about n - 1 + (M + 1) per boundary ANDs, M + boundaries wires."""
    n = len(tgt)
    segs = [(s, min(n, s + M)) for s in range(0, n, M)]
    outs = []
    cur = cin
    for k, (s, e) in enumerate(segs):
        last = k == len(segs) - 1
        out = _gidney_seg(c, src[s:e], tgt[s:e], cur, keep_out=not last)
        if not last:
            outs.append((s, e, cur, out))
            cur = out[0]
    for s, e, cin_k, out in reversed(outs):
        _erase_boundary(c, src[s:e], tgt[s:e], cin_k, out)


def add(c, src, tgt, cin=None, force=None):
    n = len(tgt)
    if n == 0:
        return
    src = list(src) + [None] * (n - len(src))
    need = gidney_need(src, tgt, cin)
    cap = getattr(c, "cap", 10 ** 9)
    if force == "cuccaro":
        return add_cuccaro(c, src, tgt, cin)
    if force == "gidney" or c.live + need <= cap:
        return add_gidney(c, src, tgt, cin)
    if getattr(c, "narrow", False):
        if c.live + gidney_fresh(src, tgt, cin) <= cap:
            return add_gidney(c, src, tgt, cin)
        import invsqrt_fc_narrow as NW
        if c.live + 2 <= cap and NW.feasible(src, tgt, cin, cap - c.live + len(getattr(c, "_fk_lend", None) or [])):   # CH1 P4
            return NW.block(c, src, tgt, cin)
        if getattr(c, "cuccaro_fwd", False) and NW.cuccaro_need(src, tgt, cin) <= cap - c.live:
            NW.cuccaro_add(c, src, tgt, cin, {})       # in place, unitary: its mirror is its undo
            return
    free = cap - c.live
    best = None
    for M in range(2, max(3, free)):
        nb = -(-n // M) - 1
        if M + 1 + nb <= free:
            best = M                          # largest feasible segment
    if best is not None:
        return add_segmented(c, src, tgt, cin, best)
    return add_cuccaro(c, src, tgt, cin)


def frame(c, ctl, tgt):
    """flip every target bit iff ctl == 0."""
    for t in tgt:
        c.cx(ctl, t)
        c.x(t)


def fadd(c, ctl, src, tgt, cin=None, force=None):
    frame(c, ctl, tgt)
    add(c, src, tgt, cin, force)
    frame(c, ctl, tgt)
