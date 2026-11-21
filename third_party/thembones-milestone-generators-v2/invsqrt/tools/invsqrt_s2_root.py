#!/usr/bin/env python3
"""invsqrt_s2_root.py -- S2 family R: the root split. The in-place root keeps m bits (v) and a residual of
22 - m bits; d = 2 Q + 1 with Q = isqrt(n 4^(m-10)) in [2^m, 2^(m+1)); delta = n 4^(m-9) - d^2; the word
Y = 2^(28+m) / sqrt(d^2 + delta) = (2^(28+m) + 2^(m+1) U) / d with U = 2^27 (d / D - 1) (m = 12 is shipped).

Generic exact model (own code), used for m = 12 and m = 11 alike so that their DIFFERENCE is measured in one model:
  chain_m      v, d, delta, h, c (delta = 2^(2m-20) h + c; c a function of v's low bits) for every input
  U models     ideal (rounded), first-order only (no w^2), sign-kinked first order, table fits
  divide_m     non-restoring rows for quotient weights 2^26..2^0; in rows with weight 2^b < 2^b_t the divisor is
               rounded to drop t_b = b_t - b low bits (the shipped tail drops 11 - b bits for b <= 9)
Commands:
  invsqrt_s2_root.py WD ideal M          ideal U and full divider: LP margin (the arithmetic's own limit)
  invsqrt_s2_root.py WD terms M          U without the w^2 term; with a sign-kinked first-order coefficient
  invsqrt_s2_root.py WD tail M B_T       ideal U, tail truncation from weight 2^B_T down
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_s2_core as CORE  # noqa: E402


def isqrt_vec(a):
    r = np.floor(np.sqrt(a.astype(np.float64))).astype(np.int64)
    r = np.where(r * r > a, r - 1, r)
    r = np.where((r + 1) * (r + 1) <= a, r + 1, r)
    return r


def chain_m(ch, m):
    """root split m: v (m bits), d (m+2 bits, odd), delta, h, c; n recovered from the shipped chain."""
    N64 = ch["d"] * ch["d"] + ch["delta"]           # 64 n
    n = N64 // 64
    sc = n << (2 * (m - 10)) if m >= 10 else n >> (2 * (10 - m))
    Q = isqrt_vec(sc)
    assert Q.min() >= (1 << m) and Q.max() < (1 << (m + 1))
    v = Q - (1 << m)
    d = 2 * Q + 1
    nm = n << (2 * (m - 9))
    delta = nm - d * d
    g = 1 << (2 * m - 20)                            # 16 for m = 12 ... wait: see below
    # delta = 4^(m-9) n - d^2: d^2 mod 4^(m-9) is fixed by v's low bits, so delta = 4^(m-9) h + c
    g = 4 ** (m - 9)
    c = (-(d * d)) % g
    h = (delta - c) // g
    assert np.array_equal(g * h + c, delta)
    Dr = np.sqrt(nm.astype(np.float64))
    return dict(m=m, v=v, d=d, delta=delta, h=h, c=c, g=g, D=Dr, n=n)


def U_ideal(cm):
    """U* = 2^27 (d / D - 1) (float)."""
    return np.ldexp(cm["d"].astype(np.float64) / cm["D"] - 1.0, 27)


def kk_drop(kk):
    """S1's tail schedule kk (rows j = 25..35) as {word weight exponent b: bits dropped}. In the shipped word
    p = 2^28 + Qw/2 + 1 - 2 s36 (phase formula of the divider), row j's digit has word weight 2^(36-j) [P]."""
    return {36 - j: int(kk[j - 25]) for j in range(25, 36)}


TAILFORM = {"form": "ship"}      # "ship": the shipped 2A + 1; "round": d rounded to a multiple of 2^t


def divide_m(cm, U, b_t=10, rows_drop=None, b_min=0):
    """Y = NR quotient of N = 2^(28+m) + 2^(m+1) U by d; rows b = 26..b_min (weight 2^b); the divisor is rounded
    to drop t_b low bits in row b (rows_drop per b, else t_b = b_t - b). With b_min = 2 the last digit is the
    remainder's sign at weight 2 (the shipped word's resolution). Returns the word p (int64)."""
    m = cm["m"]
    d = cm["d"].astype(np.int64)
    N = (np.int64(1) << (28 + m)) + (U.astype(np.int64) << (m + 1))
    R = N.copy()
    Qw = np.zeros(len(d), dtype=np.int64)
    V = d >> 1
    for b in range(26, b_min - 1, -1):
        t = rows_drop.get(b, 0) if rows_drop is not None else max(0, b_t - b)
        if t and TAILFORM["form"] == "round":
            dt = ((d + (1 << (t - 1))) >> t) << t
        elif t:
            # the shipped tail row's divisor [P, the divider recurrence]: 2 A + 1 with
            # A = round(V to a multiple of 2^(t+1)) - (V & 1), V = d // 2
            A = ((((V >> (t + 1)) + ((V >> t) & 1)) << (t + 1)) - (V & 1))
            dt = 2 * A + 1
        else:
            dt = d
        s = np.where(R >= 0, 1, -1)
        R = R - s * (dt << b)
        Qw = Qw + s * (1 << b)
    if b_min == 0:
        # NR quotient -> integer: Qw is odd; the remainder's sign picks Qw or Qw - 1
        return np.where(R < 0, Qw - 1, Qw)
    # the final sign digit at weight 2^(b_min - 1) (shipped: rows end at weight 2, then 1 - 2 s36 = +-1)
    s = np.where(R >= 0, 1, -1)
    return Qw + s * (1 << (b_min - 1))


def evaluate(ph, ch, p):
    r = CORE.evaluate_word(ph, ch, p)
    return r


def report(name, ph, ch, p):
    r = evaluate(ph, ch, p)
    e = p.astype(np.float64) - ch["Wf"]
    print(f"{name:52s} word err {e.min():7.2f}..{e.max():6.2f}  m(pub angles) {r['m_pa']:.4g} fail "
          f"{r['nfail']}  LP m {r['m_lp']:.4g}  T/rot {CORE.tper(r['m']) if r['m'] > 0 else float('nan'):.3f}",
          flush=True)
    return r


def main():
    wd, cmd, m = sys.argv[1], sys.argv[2], int(sys.argv[3])
    ch = CORE.load(wd)
    ph = CORE.phase(ch, wd)
    cm = chain_m(ch, m)
    if m == 12:
        assert np.array_equal(cm["v"], ch["v"]) and np.array_equal(cm["d"], ch["d"])
    print(f"m = {m}: v {cm['v'].min()}..{cm['v'].max()}, d {cm['d'].min()}..{cm['d'].max()}, h "
          f"{cm['h'].min()}..{cm['h'].max()} ({int(math.log2(cm['h'].max() - cm['h'].min() + 1) + 0.999)} bits), "
          f"c values {sorted(set(np.unique(cm['c']).tolist()))[:8]}", flush=True)
    Us = U_ideal(cm)
    print(f"U* range {Us.min():.1f}..{Us.max():.1f}", flush=True)
    if cmd == "ideal":
        U = np.round(Us).astype(np.int64)
        report(f"m {m}: ideal U, exact divider (no tail truncation)", ph, ch, divide_m(cm, U, b_t=0))
        report(f"m {m}: ideal U, tail drops 10 - b bits (b < 10)", ph, ch, divide_m(cm, U, b_t=10))
        report(f"m {m}: ideal U, tail drops 11 - b bits (b < 11)", ph, ch, divide_m(cm, U, b_t=11))
    elif cmd == "terms":
        d = cm["d"].astype(np.float64)
        w = cm["delta"].astype(np.float64) / (d * d)
        U1 = -np.ldexp(w / 2, 27)                           # first order only
        U2 = np.ldexp(-w / 2 + 3 * w * w / 8, 27)           # second order
        for name, Uf in (("first order only", U1), ("second order", U2), ("exact", Us)):
            report(f"m {m}: U {name}, exact divider", ph, ch, divide_m(cm, np.round(Uf).astype(np.int64), 0))
        # sign kink: the w^2 term replaced by its best fit linear in |w| on the root's range: x^2 ~ |x| - 1/8
        # (x = w / wm, wm = 2 / d per root; error 1/8 of the term) [P]; the form used by word_R
        q2 = np.ldexp(3 * w * w / 8, 27)
        wm = 2.0 / d                                          # |w| bound per root
        kink = np.ldexp(3 * wm / 8, 27) * np.abs(w)
        off = np.ldexp(3 * wm * wm / 8, 27) / 8
        Uk = U1 + kink - off
        report(f"m {m}: U first order + sign kink (no product of w^2)", ph, ch,
               divide_m(cm, np.round(Uk).astype(np.int64), 0))
        print(f"   w^2 term max {q2.max():.2f} units; kink residual max {np.abs(Uk - Us).max():.2f} units", flush=True)
    elif cmd == "tail":
        bt = int(sys.argv[4])
        U = np.round(Us).astype(np.int64)
        report(f"m {m}: ideal U, tail from 2^{bt}", ph, ch, divide_m(cm, U, b_t=bt))


if __name__ == "__main__":
    main()
