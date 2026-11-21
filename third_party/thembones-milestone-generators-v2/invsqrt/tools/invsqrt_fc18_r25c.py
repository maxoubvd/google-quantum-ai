#!/usr/bin/env python3
"""invsqrt_fc18_r25c.py -- FC eighteenth holder: R25's bits 7 and 8 re-made ONLY UNDER THE OUTCOME that needs them.
Copied into a source snapshot; active only with cfg "_fc18r25c": 1 (on top of "_r63chain", "_ch1r24", "_fc17skipe").

TODAY (Egq6n25): row 24's undo ends by re-making R25's bits 7 (w7', 1 AND, also paying skip E's m8 part by Z under
m8) and 8 (w8', 1 AND); their only readers are the chain block's mirror (ROW9's block 37), which erases bit 9 (digit
24) by its one-AND form (factors read bit 8), then bit 8 (factors read bit 7), then bit 7, each by measurement (free).
So the two re-makes (2.0 a shot) exist only to pay the PHASE of three measured erases.
HERE: bits 7 and 8 are never re-made unconditionally.  The m8 part (-1)^(m8 w7') is paid by Cliffords (w7' is one AND
of two affine forms of live wires: quadratic).  The chain mirror X-measures bit 9 FIRST (outcome g); only where g = 1
it re-makes bits 7 and 8 by the chain's own forward ops (2 ANDs gated by g: CANDC), applies bit 9's phase (CZ of its two
folded factors + Z on its linear terms, all under g), and erases bits 8 and 7 again (measured, gated by g: XANDE).
Count: 2 x 0.5 = 1.0 a shot against 2.0: - 1.0 a shot [P].  Bit 9 has degree >= 4 in the wires live without bits 7
and 8 [M ALL, tools/invsqrt_fc_s18r25.py], so two gated non-Clifford gates is this route's floor [P].
The frame segment between row 24 and the chain (X, CX(digit 24) on row 24's ten window wires) still runs on bits 7
and 8's wires, which are |0> here: the mirror's head clears that residue (CX(digit 24, w), X(w)) on both.
"""
import sys


def phase_and(c, left, right, cond):
    """(-1)^(cond * AND(XOR left, XOR right)) by folds + CZ (left / right: wires and "1"), no Toffoli."""
    lw = [w for w in left if w != "1"]
    rw = [w for w in right if w != "1"]
    hl = [w for w in lw if w not in rw]
    hr = [w for w in rw if w not in lw]
    assert hl and hr, ("no host wires", left, right)
    hl, hr = hl[0], hr[0]
    ops = []
    for w in lw:
        if w != hl:
            ops.append(("CX", w, hl, None))
    if "1" in left:
        ops.append(("X", hl, None))
    for w in rw:
        if w != hr:
            ops.append(("CX", w, hr, None))
    if "1" in right:
        ops.append(("X", hr, None))
    c.ops.extend(ops)
    c.cz(hl, hr, cond)
    c.ops.extend(ops[::-1])


def pay_w7_phase(c, tgt, d, v, m8):
    """(-1)^(m8 w7'), w7' = v3^v4^v5^v6^w6' ^ (v3^d^w6')(v3^v5) (CH1's form A, as _ch1_w7 makes it), by Cliffords."""
    w6 = tgt[7]
    for w in (v[3], v[4], v[5], v[6], w6):
        c.z(w, m8)
    phase_and(c, [v[3], d, w6], [v[3], v[5]], m8)


def _gate(c, o, g):
    """the op o executed only where classical bit g is 1."""
    k = o[0]
    if k == "CX":
        assert o[3] is None, o
        return ("CX", o[1], o[2], g)
    if k == "X":
        assert o[2] is None, o
        return ("X", o[1], g)
    if k == "ANDC":
        return ("CANDC", o[1], o[2], o[3], g)
    if k == "ANDE":
        return ("XANDE", o[1], o[2], o[3], g, c.bit())
    if k in ("CZ", "Z", "NEG"):
        return o                                   # already gated by g (built below)
    raise ValueError("op %s cannot be gated" % (o,))


def chain_mirror(c, K_, prefix, ops7, ops8, ops9, wire):
    """the chain block's undo with bits 7 and 8 absent (see the module text)."""
    t7, t8, t9 = wire[7], wire[8], wire[9]
    save = c.ops
    c.ops = []
    for w in (t7, t8):                             # the frame residue on bits 7 / 8's wires: 1 ^ d -> 0
        c.cx(t9, w)
        c.x(w)
    g = c.hmr(t9)                                  # bit 9 (digit 24) X-measured
    head = c.ops
    c.ops = save
    ph9 = []
    for o in ops9:                                 # bit 9's phase: its AND -> CZ under g, its linear XORs -> Z under g
        if o[0] == "ANDC":
            assert o[3] == t9
            ph9.append(("CZ", o[1], o[2], g))
        elif o[0] == "CX" and o[2] == t9:
            ph9.append(("Z", o[1], g))
        elif o[0] == "X" and o[1] == t9:
            ph9.append(("NEG", g))
        else:
            assert t9 not in o[1:3], o
            ph9.append(o)
    body = list(ops7) + list(ops8) + ph9 + K_.Circ.mirror(ops8) + K_.Circ.mirror(ops7)
    gated = [_gate(c, o, g) for o in body]
    return head + gated + K_.Circ.mirror(prefix)


if __name__ == "__main__":
    print(__doc__)
    sys.exit(0)
