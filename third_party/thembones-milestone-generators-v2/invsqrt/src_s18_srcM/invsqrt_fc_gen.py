#!/usr/bin/env python3
"""invsqrt_fc_gen.py -- FC: generator of the family-C invsqrt circuit (own code; kit + arith modules).

Stages (forward; the inverse half is Circ.mirror of the forward list):
  NORM    x -> n = x 4^k, flags F0..F3 = binary k
  ROOT    non-restoring root of n >> 6, 7 steps after the leading 1, in place; root bits stored as remainder signs
  CENTRE  delta = n - (8Q + 4)^2 in one framed add (includes the non-restoring correction)
  DIV1    theta = (delta << 5) / d   exact rows t = 10..0, line tail rows L2, final sign   (20 digits)
  DIV2    P = theta / d              exact rows t = 10..0, line tail rows Q2, final sign   (20 digits)
  HO      per-digit non-restoring higher-order stage (fq fraction bits), then the merge to binary U8 (+ held g)
  MAIN    N = 2^35 + 2^5 U8, NR rows b = 26..5 (feed), tail rows b = 4..0 (ship form, b_t = 4), final sign
  END     27-slot word routed by k (cyclic), slot rotations + 10 wrap features
Each stage is checked against an integer reference on the given inputs.
  invsqrt_fc_gen.py test N_SAMPLE [stage ...]     (local sample)
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_fc_kit as K   # noqa: E402
import invsqrt_fc_arith as A  # noqa: E402
import invsqrt_fc_gnr as GNR  # noqa: E402

ONE = A.ONE
MAX_X = 3133452


# ------------------------------------------------------------------ references (integers per input)
def ref_norm(x):
    n = x.astype(np.int64).copy()
    k = np.zeros_like(n)
    for _ in range(11):
        sh = n < (1 << 20)
        n = np.where(sh, n << 2, n)
        k = k + sh
    return n, k


def isqrt(a):
    r = np.floor(np.sqrt(a.astype(np.float64))).astype(np.int64)
    r = np.where(r * r > a, r - 1, r)
    r = np.where((r + 1) * (r + 1) <= a, r + 1, r)
    return r


def ref_root(n):
    """non-restoring root of N' = n >> 6, 8 digits: returns Q, the signed final remainder r_nr."""
    Np = n >> 6
    rem = ((Np >> 14) & 3) - 1          # step 7: q7 = 1
    Q = np.ones_like(n)
    for j in range(6, -1, -1):
        pair = (Np >> (2 * j)) & 3
        neg = rem < 0
        rem = np.where(neg, 4 * rem + pair + 4 * Q + 3, 4 * rem + pair - 4 * Q - 1)
        Q = 2 * Q + (rem >= 0)
    return Q, rem


# ------------------------------------------------------------------ helpers
def multi_and(c, ins, out):
    """out ^= AND(ins) (out a zero wire), temps erased; len(ins) >= 2."""
    if len(ins) == 2:
        c.and_compute(ins[0], ins[1], out)
        return
    temps = []
    cur = ins[0]
    for k in range(1, len(ins) - 1):
        t = c.and_compute(cur, ins[k])
        temps.append((cur, ins[k], t))
        cur = t
    c.and_compute(cur, ins[-1], out)
    for a, b, t in reversed(temps):
        c.and_erase(a, b, t)


# ------------------------------------------------------------------ stages
def st_norm(c, X):
    F = c.alloc(4)
    for f, lo, s in ((3, 6, 16), (2, 14, 8), (1, 18, 4), (0, 20, 2)):
        bits = X[lo:22]
        for b in bits:
            c.x(b)
        multi_and(c, bits, F[f])
        for b in bits:
            c.x(b)
        for i in range(21 - s, -1, -1):
            c.cswap(F[f], X[i], X[i + s])
    return F


def st_root(c, X):
    """W(p) = X[p + 6] (N' positions); position 16 is a constant 0. After: root bits q_t stored as NOT q at
    W(t + 9) (t = 0..6), remainder r_nr (signed 10 bits) at W(0..9) (its sign = W(9) = NOT q0)."""
    W = lambda p: X[p + 6]  # noqa: E731
    # step 7: rem = pair - 1 on W(14), W(15)
    c.cx(W(14), W(15))
    c.x(W(15))
    c.x(W(14))
    for j in range(6, -1, -1):
        tgt = [W(p) for p in range(2 * j, j + 10)]          # window [2j, j+9]
        if j == 6:
            # sign of rem_7 is the constant 0: subtract 4*1 + 1 = 5
            for t in tgt:
                c.x(t)
            A.add(c, [ONE, None, ONE], tgt)
            for t in tgt:
                c.x(t)
            continue
        sig = W(j + 10)
        # operand 4Q + 1 + 2 sigma; Q_0 = q_{j+1} = NOT sigma, Q_k = q_{j+1+k} = NOT W(j+10+k), leading Q = 1
        cp = c.alloc(1)[0]
        c.cx(sig, cp)
        qsrc = [W(j + 10 + k) for k in range(0, 6 - j)]      # NOT of these are Q_0 .. Q_{5-j}
        A.frame(c, sig, tgt)
        for w in qsrc:
            c.x(w)
        src = [ONE, cp] + qsrc + [ONE]
        A.add(c, src, tgt)
        for w in qsrc:
            c.x(w)
        A.frame(c, sig, tgt)
        c.cx(sig, cp)
        c.free_([cp])
    # root bits: v_t = q_t = NOT W(t + 9)


def st_centre(c, X):
    """delta (15-bit two's complement) on X[0..14]; c = W(9) = X[15] = NOT q0 (sign of r_nr).
    delta = R' + (c ? 64Q + 48 : -(64Q + 16)); operand 64Q + 16 + 32c on X[4..14]."""
    cw = X[15]
    tgt = X[4:15]                                          # bits 4 .. 14
    cp = c.alloc(1)[0]
    c.cx(cw, cp)
    # Q_0 = q0 = NOT cw (use cw itself, X'ed), Q_t = NOT X[15 + t] for t = 1..6, Q_7 = 1
    A.frame(c, cw, tgt)
    qs = [X[15 + t] for t in range(0, 7)]
    for w in qs:
        c.x(w)
    src = [ONE, cp] + qs + [ONE]                           # bits 4, 5, 6..12, 13
    A.add(c, src, tgt)
    for w in qs:
        c.x(w)
    A.frame(c, cw, tgt)
    c.cx(cw, cp)
    c.free_([cp])
    # v in true form on X[15..21]
    for w in X[15:22]:
        c.x(w)


def ref_div(reg, d, a_top, sched):
    """the re-pricer's nr_line arithmetic: returns the list of sign arrays (1 = negative) in row order
    (exact rows, tail rows, final sign) and the quotient value sum s 2^(n-1-j)."""
    sg = []
    reg = reg.astype(np.int64).copy()
    for t in range(a_top, -1, -1):
        s = np.where(reg < 0, -1, 1)
        sg.append((reg < 0).astype(np.int64))
        reg = reg - s * (d << t)
    rr = reg.copy()
    V = d >> 1
    for kk in sched:
        s = np.where(rr < 0, -1, 1)
        sg.append((rr < 0).astype(np.int64))
        if kk:
            Av = (((V >> (kk + 1)) + ((V >> kk) & 1)) << (kk + 1)) - (V & 1)
            dt = 2 * Av + 1
        else:
            dt = d
        rr = 2 * rr - s * dt
    sg.append((rr < 0).astype(np.int64))
    n = len(sg)
    q = np.zeros(len(d), dtype=np.int64)
    for j, g in enumerate(sg):
        q += (1 - 2 * g) << (n - 1 - j)
    return sg, q


def st_div(c, R, v, a_top, sched, u7v=()):
    """in-place division of the register R (positions 0 .. a_top + 10, two's complement, sign at the top) by
    d = 2^8 + 2v + 1. Returns the digit wires in row order (weight 2^(n-1-j)) and the 9 final register wires."""
    dsrc = [ONE] + list(v) + [ONE]
    Vb = list(v) + [ONE]
    digits = []
    for t in range(a_top, -1, -1):
        sig = R[t + 10]
        A.fadd(c, sig, dsrc + [None], R[t:t + 10])
        digits.append(sig)
    # first tail row (kk = 0): u' = r - s V - [s = +1] on R[0..8], digit R[9]
    assert sched[0] == 0
    sig = R[9]
    tgt = R[0:9]
    A.frame(c, sig, tgt)
    c.x(sig)
    A.add(c, Vb + [None], tgt, cin=sig)
    c.x(sig)
    A.frame(c, sig, tgt)
    digits.append(sig)
    u = list(R[0:9])
    for kk in sched[1:]:
        u = [u[8]] + u[0:8]
        su = u[0]
        src = Vb[kk + 1:]
        cin = Vb[kk]
        A.fadd(c, su, src + [None] * (8 - kk - len(src)), u[kk + 1:9], cin=cin)
        c.cx(v[0], u[0])
        digits.append(u[0])
    digits.append(u[8])
    # clean-up of the window: rotation-row digits hold sigma xor v0; u[7] is NOT (first tail digit)
    for w in digits[a_top + 2:-1]:
        c.cx(v[0], w)
    c.cx(digits[a_top + 1], u[7])
    c.x(u[7])
    for t in u7v:
        c.cx(v[t], u[7])
    return digits, u[7]


def div_input(c, dig, free_w):
    """register for the next division from the previous division's 20 digit wires (dig[idx] = sigma of weight
    2^(19 - idx)): theta = 2B + 1 - 2^20 as 21-bit two's complement: [1, NOT s_0 .. NOT s_18, s_19]."""
    c.x(free_w)                                   # constant 1 at position 0
    R = [free_w]
    for pos in range(1, 20):
        w = dig[20 - pos]
        c.x(w)
        R.append(w)
    R.append(dig[0])
    return R


def apply_terms(c, reg, dq, i, sw):
    """reg (LSB first, units 2^-fq) += sum dq[key] * XOR(sigma over key) + constant; sw(j) = wire of sigma_j.
    Pairs {i, j}: sigma_i is XORed into sigma_j in place; other XOR sets on temps. Positive and negative parts
    in greedy layers of distinct positions; constant bits merged into the positive layers."""
    W = len(reg)
    pairs = sorted({j for (S, E) in dq if len(S) == 2 and i in S for j in S if j != i})
    for j in pairs:
        c.cx(sw(i), sw(j))
    temps = []

    def bitwire(S):
        S = set(S)
        if len(S) == 2 and i in S:
            return sw(next(j for j in S if j != i))
        if S == {i}:
            return sw(i)
        if len(S) == 1 and next(iter(S)) not in pairs:
            return sw(next(iter(S)))
        t = c.alloc(1)[0]
        cn = []
        par = 0
        for j in sorted(S):
            cn.append(sw(j))
            if j in pairs:
                par ^= 1
        if par:
            cn.append(sw(i))
        for w in cn:
            c.cx(w, t)
        temps.append((t, cn))
        return t

    const = dq.get((frozenset(), frozenset()), 0)
    pos_bits, neg_bits = [], []
    for (S, E), v in sorted(dq.items(), key=lambda kv: (len(kv[0][0]), sorted(kv[0][0]))):
        if not S:
            continue
        assert not E
        w = bitwire(S)
        a, b = abs(v), 0
        while a:
            if a & 1:
                (pos_bits if v > 0 else neg_bits).append((b, w))
            a >>= 1
            b += 1
    if neg_bits and not pos_bits:
        cm = (-const) % (1 << W)                 # r + C - N = r - (N - C): constant joins the negative layers
        for b in range(W):
            if (cm >> b) & 1:
                neg_bits.append((b, ONE))
    else:
        cm = const % (1 << W)
        for b in range(W):
            if (cm >> b) & 1:
                pos_bits.append((b, ONE))
    for bits, neg in ((pos_bits, False), (neg_bits, True)):
        layers = []
        for b, w in bits:
            if b >= W:
                continue
            for L in layers:
                if b not in L:
                    L[b] = w
                    break
            else:
                layers.append({b: w})
        for L in layers:
            lo = min(L)
            src = [L.get(b, None) for b in range(lo, W)]
            tgt = reg[lo:W]
            if neg:
                for t in tgt:
                    c.x(t)
            A.add(c, src, tgt)
            if neg:
                for t in tgt:
                    c.x(t)
    for t, cn in reversed(temps):
        for w in reversed(cn):
            c.cx(w, t)
        c.free_([t])
    for j in reversed(pairs):
        c.cx(sw(i), sw(j))


def st_ho(c, Y, fq=3, keep=-10.0, K=0):
    """per-digit non-restoring HO on Y (20 wires, U1) + guard G + fq fraction wires. Returns (G, Afr, pos)."""
    G = c.alloc(1)[0]
    c.cx(Y[19], G)
    Afr = c.alloc(fq)

    def pos(p):
        return Afr[p + fq] if p < 0 else (Y[p] if p <= 19 else G)

    q_prev = {}
    digits = []
    for i in range(18, -1, -1):
        S, Pb, Qb = pos(i + 2), pos(i + 1), pos(i)
        c.cx(S, Pb)
        c.cx(Qb, Pb)
        c.x(Pb)
        c.x(Qb)
        digits.append(i)
        q_now = GNR.quant(GNR.h_basis(digits, [], 0, fq), fq, keep)
        keys = set(q_now) | set(q_prev)
        dq = {k: q_now.get(k, 0) - q_prev.get(k, 0) for k in keys}
        if i == 0 and K:
            ck = (frozenset(), frozenset())
            dq[ck] = dq.get(ck, 0) + K * (1 << fq)
        dq = {k: v for k, v in dq.items() if v}
        reg = [pos(p) for p in range(-fq, i + 2)]
        apply_terms(c, reg, dq, i, lambda j: pos(j + 2))
        q_prev = q_now
    return G, Afr, pos


def st_merge(c, pos, fq=3):
    """HO register (digits sigma_i at pos(i+2), remainder r at pos(-fq..1)) -> U8 + 2^22 = 16 B' + L' in binary
    on 23 wires (fraction, NOT l, B'_0..B'_18) and the held bit g. Returns (U8 wires LSB first, g, freed S)."""
    # with K = 1: U8 + 2^22 = 16 B + 8 (r_int + 2) + frac, r_int = -2 S + l in [-2, 1]: the carry into B is
    # NOT S; the held bit is S itself (B' = B + NOT S, low nibble = 8 l + frac).
    S, l = pos(1), pos(0)
    B = [pos(i + 2) for i in range(19)]
    for w in B:
        c.x(w)                                  # b_i = NOT sigma_i
    c.x(S)
    A.add(c, [], B, cin=S)                      # B + NOT S
    c.x(S)
    U8 = [pos(p) for p in range(-fq, 0)] + [l] + B
    c.x(U8[-1])                                 # U8 = (16 B' + 8 l + frac) - 2^22 as 23-bit two's complement
    return U8, S, None


def ref_main(U8, d, bt=4):
    """N = 2^35 + 2^5 U8; NR rows b = 26..0, divisor rounded (ship form) for b < bt; returns the sign list
    (row order), Y (integer word) and the final remainder."""
    N = (np.int64(1) << 35) + (U8.astype(np.int64) << 5)
    R = N.copy()
    V = d >> 1
    sg = []
    Qw = np.zeros(len(d), dtype=np.int64)
    for b in range(26, -1, -1):
        t = max(0, bt - b)
        if t:
            Av = ((((V >> (t + 1)) + ((V >> t) & 1)) << (t + 1)) - (V & 1))
            dt = 2 * Av + 1
        else:
            dt = d
        s = np.where(R >= 0, 1, -1)
        sg.append((R < 0).astype(np.int64))
        R = R - s * (dt << b)
        Qw = Qw + s * (1 << b)
    fin = (R < 0).astype(np.int64)
    Y = Qw - fin
    return sg, fin, Y


def st_main(c, U8, v):
    """main divider on N = 2^35 + 2^5 U8 (U8 at N positions 5..27, 8 fresh wires at 28..35). Returns the word
    wires w_0 .. w_26 (Y = sum w_j 2^j) after the clean-up, and the remainder wires left."""
    E = c.alloc(8)
    us = U8[22]
    for t in range(7):
        c.cx(us, E[t])
    c.cx(us, E[7])
    c.x(E[7])

    def pos(p):
        return U8[p - 5] if p <= 27 else E[p - 28]

    dsrc = [ONE] + list(v) + [ONE]
    Vb = list(v) + [ONE]
    win = [pos(p) for p in range(26, 36)]
    for t in win:
        c.x(t)
    A.add(c, dsrc + [None], win)                    # row 26: sign is the constant 0 -> subtract d 2^26
    for t in win:
        c.x(t)
    digit = {}
    for b in range(25, 4, -1):
        sig = pos(b + 10)
        A.fadd(c, sig, dsrc + [None], [pos(p) for p in range(b, b + 10)])
        digit[b] = sig
    sig = pos(14)                                    # b = 4: first tail row (kk = 0)
    tgt = [pos(p) for p in range(5, 14)]
    A.frame(c, sig, tgt)
    c.x(sig)
    A.add(c, Vb + [None], tgt, cin=sig)
    c.x(sig)
    A.frame(c, sig, tgt)
    digit[4] = sig
    u = list(tgt)
    for b, kk in ((3, 1), (2, 2), (1, 3), (0, 4)):
        u = [u[8]] + u[0:8]
        src = Vb[kk + 1:]
        A.fadd(c, u[0], src + [None] * (8 - kk - len(src)), u[kk + 1:9], cin=Vb[kk])
        c.cx(v[0], u[0])
        digit[b] = u[0]
    fin = u[8]
    for b in range(0, 4):
        c.cx(v[0], digit[b])
    # word: w_0 = NOT fin, w_{b+1} = NOT sigma_b
    word = [fin] + [digit[b] for b in range(0, 26)]
    for w in word:
        c.x(w)
    rest = [w for w in u if w not in word]
    return word, rest, digit, fin


THR = {1: "f0|f1|f2|f3", 2: "f1|f2|f3", 3: "(f0&f1)|f2|f3", 4: "f2|f3", 5: "(f2&(f0|f1))|f3", 6: "(f2&f1)|f3",
       7: "(f2&f1&f0)|f3", 8: "f3", 9: "f3&(f0|f1)", 10: "f3&f1"}


def _or2(c, a, b):
    """t = a OR b on a fresh wire (1 AND); returns (t, undo ops recipe)."""
    c.x(a)
    c.x(b)
    t = c.and_compute(a, b)
    c.x(t)
    c.x(a)
    c.x(b)
    return t, ("or", a, b, t)


def _and2(c, a, b):
    t = c.and_compute(a, b)
    return t, ("and", a, b, t)


def _undo(c, rec):
    kind, a, b, t = rec
    if kind == "or":
        c.x(a)
        c.x(b)
        c.x(t)
        c.and_erase(a, b, t)
        c.x(a)
        c.x(b)
    else:
        c.and_erase(a, b, t)


def threshold(c, F, cst):
    """wire holding [k >= cst] for binary k on F (k <= 10), with the undo records."""
    f0, f1, f2, f3 = F
    recs = []

    def o(a, b):
        t, r = _or2(c, a, b)
        recs.append(r)
        return t

    def n(a, b):
        t, r = _and2(c, a, b)
        recs.append(r)
        return t

    if cst == 1:
        w = o(o(o(f0, f1), f2), f3)
    elif cst == 2:
        w = o(o(f1, f2), f3)
    elif cst == 3:
        w = o(o(n(f0, f1), f2), f3)
    elif cst == 4:
        w = o(f2, f3)
    elif cst == 5:
        w = o(n(f2, o(f0, f1)), f3)
    elif cst == 6:
        w = o(n(f2, f1), f3)
    elif cst == 7:
        w = o(n(n(f2, f1), f0), f3)
    elif cst == 8:
        w = f3
    elif cst == 9:
        w = n(f3, o(f0, f1))
    else:
        w = n(f3, f1)
    return w, recs


def st_end(c, word, F, ang=None, feat=None):
    """route the 27-bit word cyclically by k, slot rotations, wrap features, unroute. ang / feat: angle
    numerators (units 2^-64) per slot / feature (default ideal). Returns the list of ZP op indices."""
    route = []
    c0 = len(c.ops)
    for f in range(4):
        s = 1 << f
        cyc = [(t * s) % 27 for t in range(27)]
        for t in range(1, 27):
            c.cswap(F[f], word[cyc[0]], word[cyc[t]])
    route = c.ops[c0:]
    zps = []
    for p in range(27):
        num = ang[p] if ang is not None else (1 << (p + 27))
        zps.append(c.zp(word[p], num))
    for p in range(10):
        thr, recs = threshold(c, F, p + 1)
        t2 = c.and_compute(thr, word[p])
        num = feat[p] if feat is not None else (((1 << (p + 54)) - (1 << (p + 27))) % (1 << 64))
        zps.append(c.zp(t2, num))
        c.and_erase(thr, word[p], t2)
        for r in reversed(recs):
            _undo(c, r)
    c.ops.extend(K.Circ.mirror(route))
    return zps


def build_full(cap=52, cap_ho=52, ang=None, feat=None):
    stages = ["root", "centre", "div1", "div2", "ho", "merge"]
    K.Circ.cap_ho = cap_ho
    c, regs = build(stages, cap=cap)
    c.tag("MAIN")
    regs["word"], regs["rest"], regs["digit"], regs["fin"] = st_main(c, regs["U8"], regs["v"])
    nfwd = len(c.ops)
    fwd = list(c.ops)
    c.tag("END")
    regs["zps"] = st_end(c, regs["word"], regs["F"], ang, feat)
    c.tag("UNDO")
    c.ops.extend(K.Circ.mirror(fwd))
    regs["nfwd"] = nfwd
    return c, regs


# ------------------------------------------------------------------ test harness
def sample_x(nsamp, seed=7):
    rng = np.random.default_rng(seed)
    xs = np.concatenate([np.arange(1, 200), rng.integers(1, MAX_X + 1, nsamp),
                         np.array([MAX_X, 1 << 20, (1 << 20) - 1, (1 << 18), 4177936 // 4 + 1])])
    return np.unique(xs).astype(np.int64)


def build(stages, cap=45):
    c = K.Circ(n_fixed=22)
    c.cap = cap
    X = list(range(22))
    regs = {"X": X}
    c.tag("NORM")
    regs["F"] = st_norm(c, X)
    if "root" in stages:
        c.tag("ROOT")
        st_root(c, X)
    if "centre" in stages:
        c.tag("CENTRE")
        st_centre(c, X)
    if "div1" in stages:
        c.tag("DIV1")
        Z = c.alloc(5)
        G = c.alloc(1)[0]
        c.cx(X[14], G)
        R1 = Z + X[0:15] + [G]
        v = X[15:22]
        regs["v"] = v
        regs["R1"] = R1
        regs["dig1"], regs["free1"] = st_div(c, R1, v, 10, [0, 0, 1, 2, 3, 4, 5, 6])
    if "div2" in stages:
        c.tag("DIV2")
        R2 = div_input(c, regs["dig1"], regs["free1"])
        regs["dig2"], regs["free2"] = st_div(c, R2, regs["v"], 10, [0, 0, 2, 3, 4, 5, 6, 7], u7v=(1,))
        # U1 = -(P + 1) / 2: bit j = sigma_j (j <= 18), bit 19 = NOT sigma_19
        Y = [regs["dig2"][19 - j] for j in range(20)]
        c.x(Y[19])
        regs["Y"] = Y
    if "ho" in stages:
        c.tag("HO")
        cap0 = c.cap
        c.cap = getattr(c, "cap_ho", cap0)
        regs["G"], regs["Afr"], regs["pos"] = st_ho(c, regs["Y"])
        c.cap = cap0
    if "merge" in stages:
        c.tag("MERGE")
        regs["U8"], regs["g"], regs["Sfree"] = st_merge(c, regs["pos"])
    return c, regs


def run_test(nsamp, stages):
    x = sample_x(nsamp)
    c, regs = build(stages)
    sim = K.Sim(len(x), c.next + 2)
    sim.set_word(regs["X"], x)
    sim.run(c.ops)
    X, F = regs["X"], regs["F"]
    n_ref, k_ref = ref_norm(x)
    k_c = sim.get_word(F)
    print(f"stages {stages}: ops {len(c.ops)} toffolis {c.toffolis()} peak {c.peak} wires used {c.next}")
    print("  per stage:", c.stage_counts())
    if "root" not in stages:
        n_c = sim.get_word(X)
        print(f"  NORM: n ok {np.array_equal(n_c, n_ref)}  k ok {np.array_equal(k_c, k_ref)}")
        return
    Q, rnr = ref_root(n_ref)
    if "centre" not in stages:
        rem = sim.get_word(X[6:16], signed=True)
        qb = sim.get_word(X[15:22])
        v_c = (~qb) & 127
        print(f"  ROOT: k ok {np.array_equal(k_c, k_ref)}  v ok {np.array_equal(v_c, Q - 128)}  r_nr ok "
              f"{np.array_equal(rem, rnr)}  low6 ok {np.array_equal(sim.get_word(X[0:6]), n_ref & 63)}")
        return
    delta_ref = n_ref - (8 * Q + 4) ** 2
    if "div1" not in stages:
        delta_c = sim.get_word(X[0:15], signed=True)
        v_c = sim.get_word(X[15:22])
        print(f"  CENTRE: v ok {np.array_equal(v_c, Q - 128)}  delta ok {np.array_equal(delta_c, delta_ref)} "
              f"range {delta_ref.min()}..{delta_ref.max()}")
        return
    d = 2 * Q + 1
    sg, q = ref_div(delta_ref << 5, d, 10, [0, 0, 1, 2, 3, 4, 5, 6])
    if "div2" in stages:
        sg2, P = ref_div(q, d, 10, [0, 0, 2, 3, 4, 5, 6, 7])
        U1 = -((P + 1) >> 1)
        U1_c = sim.get_word(regs["Y"], signed=True)
        ok2 = all(np.array_equal(sim.q[w].astype(np.int64) ^ (1 if j == 0 else 0), g)
                  for j, (w, g) in enumerate(zip(regs["dig2"], sg2)))
        if "ho" in stages:
            U0, r, ratio, terms = GNR.run_pd(U1, {"fq": 3, "lo": 0})
            if "merge" not in stages:
                pos = regs["pos"]
                rc = sim.get_word([pos(p) for p in range(-3, 2)], signed=True)
                print(f"  HO: remainder ok {np.array_equal(rc, r)}  (values {np.unique(r)[:12]})")
                return
            Phi = U0 - np.floor_divide(r, 8)
            U8_ref = 8 * (Phi + 1) + r                 # K = 1
            U8_c = sim.get_word(regs["U8"], signed=True)
            print(f"  MERGE: U8 ok {np.array_equal(U8_c, U8_ref)}  held S ones {int(sim.q[regs['g']].sum())}"
                  f"  U8 range {U8_ref.min()}..{U8_ref.max()}")
            return
        print(f"  DIV2: digits ok {ok2}  U1 ok {np.array_equal(U1_c, U1)}  free wire zero "
              f"{not sim.q[regs['free2']].any()}  U1 range {U1.min()}..{U1.max()}")
        return
    vv = Q - 128
    vbits = [((vv >> t) & 1) for t in range(7)]
    print("  DIV1 digit wires vs reference sign (match, match xor v_t):")
    for j, (w, g) in enumerate(zip(regs["dig1"], sg)):
        got = sim.q[w].astype(np.int64)
        m = [t for t in range(7) if np.array_equal(got ^ vbits[t], g)]
        print(f"    digit {j:2d} wire q{w}: exact {np.array_equal(got, g)} not {np.array_equal(got ^ 1, g)} xor_v {m}")
    print(f"    free wire q{regs['free1']} zero: {not sim.q[regs['free1']].any()}")


def target_tol(x):
    """verifier units 2^-64 turn: target = round(2^64 / sqrt(x)) mod 2^64, tol = 2^40 / sqrt(x)."""
    from decimal import Decimal, getcontext
    getcontext().prec = 50
    tg = np.array([int((Decimal(2) ** 64 / Decimal(int(xx)).sqrt()).to_integral_value()) % (1 << 64) for xx in x],
                  dtype=np.uint64)
    tol = np.ldexp(1.0, 40) / np.sqrt(x.astype(np.float64))
    return tg, tol


def full_check(x, cap=52, cap_ho=52, ang=None, feat=None, check_main=True, tgt_tol=None):
    c, regs = build_full(cap, cap_ho, ang, feat)
    sim = K.Sim(len(x), c.next + 2)
    sim.fire = []
    sim.set_word(regs["X"], x)
    nf = regs["nfwd"]
    sim.run(c.ops, 0, nf)
    out = {}
    if check_main:
        n_ref, k_ref = ref_norm(x)
        Q, _ = ref_root(n_ref)
        d = 2 * Q + 1
        delta = n_ref - (8 * Q + 4) ** 2
        sg1, q1 = ref_div(delta << 5, d, 10, [0, 0, 1, 2, 3, 4, 5, 6])
        sg2, P = ref_div(q1, d, 10, [0, 0, 2, 3, 4, 5, 6, 7])
        U1 = -((P + 1) >> 1)
        U0, r, _, _ = GNR.run_pd(U1, {"fq": 3, "lo": 0})
        Phi = U0 - np.floor_divide(r, 8)
        U8 = 8 * (Phi + 1) + r
        sgm, fin, Y = ref_main(U8, d)
        Yc = sim.get_word(regs["word"])
        out["Y_ok"] = bool(np.array_equal(Yc, Y))
        out["Y"] = Y
        out["k"] = k_ref
    sim.run(c.ops, nf, None)
    others = [w for w in range(c.next) if w not in regs["X"]]
    out["x_ok"] = bool(np.array_equal(sim.get_word(regs["X"]), x))
    out["anc_ok"] = bool(not sim.q[others].any())
    out["sign_ok"] = bool(not sim.sign.any())
    tg, tol = target_tol(x) if tgt_tol is None else tgt_tol
    err = (sim.ang - tg).view(np.int64).astype(np.float64)
    out["err"] = err
    out["tol"] = tol
    out["margin"] = float((tol - np.abs(err)).min())
    out["fire"] = np.array(sim.fire, dtype=np.uint8)
    out["c"] = c
    out["regs"] = regs
    return out


if __name__ == "__main__":
    if sys.argv[1] == "test":
        run_test(int(sys.argv[2]), sys.argv[3:])
    elif sys.argv[1] == "full":
        x = sample_x(int(sys.argv[2]))
        o = full_check(x)
        c = o["c"]
        print(f"full: ops {len(c.ops)} toffolis {c.toffolis()} peak {c.peak} wires {c.next}")
        print("  per stage:", c.stage_counts())
        print(f"  Y ok {o['Y_ok']}  x restored {o['x_ok']}  ancillas zero {o['anc_ok']}  sign ok {o['sign_ok']}  "
              f"margin (ideal angles) {o['margin']:.4g}  err range {o['err'].min():.4g}..{o['err'].max():.4g}")
