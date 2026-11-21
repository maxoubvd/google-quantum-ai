#!/usr/bin/env python3
"""r33_inv-wire2_cm.py -- R33 seat INV-WIRE2: family C (in-place root, two in-place correction divisions by the
root, main divider; NO table, coefficient or product) swept over the root split m (own code; S1's and S2's
tools imported read-only).

  n (22 bits, normalised)  ->  Q = floor(sqrt(n) 2^(m-10)) (m bits + leading one), d = 2 Q + 1,
  delta = n 4^(m-9) - d^2 (m >= 9)   or   n - 4^(9-m) d^2 (m <= 9: an integer with no dependent low bits),
  w = delta / (d^2 4^S), S = max(0, 9 - m);  word Y = (2^(28+m) + 2^(m+1) U) / d,  U = 2^27 ((1+w)^(-1/2) - 1).
  division 1: theta = delta / d      non-restoring, in place on delta's register (F1 zero bits appended), the
                                     divisor tapered (rounded right shifts) once it has passed the register's end
  division 2: P = theta / d          the same on theta's digits;  U1 = -2^26 w = -P (scaled)
  higher order: U = c0 + U1 + q(top t bits of U1),  q(w) = 2^27 ((1+w)^(-1/2) - 1 + w/2)
  main divider: invsqrt_s2_root.divide_m (generic, tail drops m - 1 - b divisor bits in the row of weight 2^b)

Exact integer arithmetic on all 3,133,452 inputs; feasibility = LP margin of S1's endpoint (invsqrt_s2_core).
Wires: ledger per stage (see price()). Toffolis: the totals printed by this file use my FIRST row law and
higher-order estimate; tools/r33_inv-wire2_report.py re-prices every result line with the checked laws
(tools/r33_inv-wire2_rowlaw.py) and the pair counts of mode "grp": use the report.
ho modes: "none"; "ideal" (the exact term); "top" (the term from the top t bits of U1: accuracy proxy);
"grp" (the stage as it has to be built: an in-place recurrence, digit groups of kg bits, t digits fixed).

Commands:
  r33_inv-wire2_cm.py WD one  '{"m":8,"xd":1,"ho":"top","t":10}'       one set (LP)
  r33_inv-wire2_cm.py WD grid OUT.jsonl NPROC [m list]                  the sweep (one JSON line per set)
  r33_inv-wire2_cm.py WD chain                                          ranges per m (no LP)
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_s2_core as CORE  # noqa: E402
import invsqrt_s2_root as RT    # noqa: E402

BASE = dict(cost=75436.96, q=46, tof=1391.4375, n_t=1, n_rot=34)


# ------------------------------------------------------------------------------------------------ helpers
def bits_signed(a):
    lo, hi = int(a.min()), int(a.max())
    b = 1
    while not (-(1 << (b - 1)) <= lo and hi < (1 << (b - 1))):
        b += 1
    return b


def nbits_vals(a):
    n = len(np.unique(a))
    return n, (0 if n <= 1 else int(math.ceil(math.log2(n))))


def fwd_row(span):
    """one non-restoring row, forward. span = signed bits of the register before the row (the doubled remainder)
    above the divisor's lowest bit. Calibrated on the line's main divider [M F1]: remainder 14 wires, doubled 16
    signed bits, odd divisor: 12 Toffolis = span - 4 (the lowest bit is constant, the top two carry no carry)."""
    return max(0, span - 4)


def inv_row(span, c):
    """its undo with c clean wires: (n - 1) + 0.5 per carry not held, n = span - 3 adder bits [P R31 / R32]."""
    n = span - 3
    if n <= 1:
        return 0.0
    return (n - 1) + 0.5 * max(0, n - 2 - c)


def chain(ch, m):
    """root split m (6..12): v, d, delta (integer), S, the ideal U (float)."""
    N64 = ch["d"] * ch["d"] + ch["delta"]            # 64 n
    n = N64 // 64
    sc = n << (2 * (m - 10)) if m >= 10 else n >> (2 * (10 - m))
    Q = RT.isqrt_vec(sc)
    assert Q.min() >= (1 << m) and Q.max() < (1 << (m + 1)), (m, int(Q.min()), int(Q.max()))
    d = 2 * Q + 1
    S = max(0, 9 - m)
    if m >= 9:
        delta = (n << (2 * (m - 9))) - d * d
        Dr = np.sqrt(n.astype(np.float64)) * 2.0 ** (m - 9)
    else:
        delta = n - ((d * d) << (2 * S))
        Dr = np.sqrt(n.astype(np.float64)) / 2.0 ** (9 - m)
    Us = np.ldexp(d.astype(np.float64) / Dr - 1.0, 27)
    w = delta.astype(np.float64) / (d.astype(np.float64) ** 2 * 4.0 ** S)
    return dict(m=m, v=Q - (1 << m), d=d, delta=delta, S=S, Us=Us, w=w, n=n)


def nr_tapered(reg, d, a_top, r, F, mbits, adj0=0):
    """in-place non-restoring division. Row i subtracts s_i T_i, s_i = sign of the register; T_i = d << t for
    t = a_top - i >= 0 (exact rows), then K = mbits + 2 + r tapered rows with T = max(1, rhu(d / 2^k)),
    k = 1..K (k = mbits + 2 is the first row whose T is 1 on every input; r > 0 adds rows with T = 1, r < 0
    stops earlier and leaves a wider final register). F = zero bits appended to the register (they need a wire
    only from the first row that reaches them). adj0 = 1 when the register is a string of +-1 digits (its
    integer needs one bit more than it has wires). Returns (quotient integer sum s_i 2^(n-1-i), final register,
    rows: list of dict(t, span, live), last digit)."""
    rows = []
    digs = []
    K = mbits + 2 + r
    ts = list(range(a_top, -1, -1)) + [-k for k in range(1, K + 1)]
    for i, t in enumerate(ts):
        if t >= 0:
            T = d << t
            tz = t
        else:
            k = -t
            T = np.maximum((d + (1 << (k - 1))) >> k, 1)
            tz = 0
        b_reg = bits_signed(reg) - (adj0 if i == 0 else 0)
        b_T = int(T.max()).bit_length() + 1
        virgin = min(max(t, 0), F)
        span = max(b_reg, b_T) - tz
        live = max(b_reg, b_T) - virgin + i            # register wires in use + digits stored before this row
        s = np.where(reg < 0, -1, 1).astype(np.int64)
        reg = reg - s * T
        digs.append(s)
        rows.append(dict(t=t, span=int(span), live=int(live)))
    n = len(digs)
    q = np.zeros(len(reg), dtype=np.int64)
    for j, s in enumerate(digs):
        q += s << (n - 1 - j)
    return q, reg, rows, digs[-1]


def a_top_of(reg, d):
    """smallest a with |reg| <= 2 d 2^a on every input (the non-restoring invariant at the first row)."""
    r = (np.abs(reg).astype(np.float64) / (2.0 * d.astype(np.float64))).max()
    return max(0, int(math.ceil(math.log2(r))) if r > 0 else 0)


def q_ho(w):
    """2^27 ((1+w)^(-1/2) - 1 + w/2): everything beyond the first order."""
    return np.ldexp(1.0 / np.sqrt(1.0 + w) - 1.0 + w / 2.0, 27)


# ------------------------------------------------------------------------------------------------ the word
def word_cm(cm, prm):
    m, d, delta, S = cm["m"], cm["d"], cm["delta"], cm["S"]
    xd = int(prm.get("xd", 0))                 # extra low digits of P beyond U (precision bought with wires)
    r1, r2 = int(prm.get("r1", 0)), int(prm.get("r2", 0))
    ho = prm.get("ho", "top")
    t_ho = int(prm.get("t", 10))
    c0 = int(prm.get("c0", 0))
    info = {}
    b_delta = bits_signed(delta)
    # --- division 1: theta_int = delta 2^F1 / (d 2^(a1 - n1 + 1))
    n_target = 27 - m + xd                                           # wires of U (+ xd extra)
    F1 = int(prm.get("F1", max(0, n_target - b_delta)))
    reg1 = delta << F1
    a1 = a_top_of(reg1, d)
    th, l1, rows1, last1 = nr_tapered(reg1, d, a1, r1, F1, m)
    n1 = len(rows1)
    # theta_int ~ (delta / d) 2^(F1 - a1 + n1 - 1)   [T_last = d 2^(a1 - n1 + 1)]
    e1 = F1 - a1 + n1 - 1
    # --- division 2: P_int = theta_int / (d 2^(a2 - n2 + 1))
    a2 = a_top_of(th, d)
    P, l2, rows2, last2 = nr_tapered(th.copy(), d, a2, r2, 0, m, adj0=1)
    n2 = len(rows2)
    e2 = -a2 + n2 - 1
    # w = delta / (d^2 4^S) = theta_int 2^-e1 / (d 4^S) = P_int 2^-(e1 + e2 + 2 S);  U1 = -2^26 w
    g = 26 - e1 - e2 - 2 * S
    if g >= 0:
        U1 = -(P << g)
        dropped = 0
    else:
        U1 = -((P + (1 << (-g - 1))) >> (-g))
        dropped = max(0, -g - 1)       # n digits +-1 are n wires and one bit more as an integer
    grp_info = None
    if ho == "grp":
        # the higher-order stage as the in-place recurrence it has to be (integer model of the circuit):
        # U solves y = G(U) = 2^26 (1 - (1 + U / 2^27)^-2) = U - 3 U^2 / 2^28 + U^3 / 2^53 - ...   (y = U1).
        # Digits of U are fixed in groups of kg bits from the top: a group is the remainder's own top bits
        # (the "divisor" is within 2^-m of 1), F collects the fixed digits, the remainder loses
        # G(F + group) - G(F) rounded to a unit; after t_ho digits the rest of the remainder passes through.
        kg = int(prm.get("kg", max(2, m - 2)))
        Lb = 27 - m
        y = U1.astype(np.int64)
        F = np.zeros_like(y)
        R = y.copy()
        Gprev = np.zeros(len(y))
        hi = Lb
        lo_end = max(0, Lb - t_ho)
        groups = []
        while hi > lo_end:
            lo = max(hi - kg, lo_end)
            D = R >> lo                               # floor: the remainder's bits from lo up (carry included)
            Fn = F + (D << lo)
            e = Fn.astype(np.float64) / 2.0 ** 27
            Gn = np.ldexp(1.0 - 1.0 / (1.0 + e) ** 2, 26)
            R = R - np.round(Gn - Gprev).astype(np.int64)
            groups.append([int(lo), int(hi), int(np.abs(D).max()).bit_length(), bits_signed(R)])
            F, Gprev, hi = Fn, Gn, lo
        Ug = F + R
        top = Lb - 1
        pairs = {c: sum(1 for ee in range(lo_end, Lb) for ff in range(ee + 1, Lb) if ee + ff >= c) for c in (24, 25, 26)}
        grp_info = dict(kg=kg, groups=groups, pairs=pairs, lo_end=int(lo_end))
        qv = (Ug - U1).astype(np.float64)
    elif ho == "grpt":
        # mode "grp" with the square built from digit pairs, the small pairs DROPPED as a circuit would:
        # S = sum over e < f with e + f >= thr of 2 a_e a_f 2^(e+f) + sum over e with 2e >= thr - 1 of a_e 4^e
        # (a = bits of |F|); G_t(F) = F - 3 S / 2^28 + F^3 / 2^53 - 5 F^4 / 2^82 (the cubic and quartic exact).
        kg = int(prm.get("kg", max(2, m - 2)))
        thr = int(prm.get("thr", 26))
        Lb = 27 - m
        y = U1.astype(np.int64)
        F = np.zeros_like(y)
        R = y.copy()
        hi = Lb
        lo_end = max(0, Lb - t_ho)
        groups = []
        npairs = sum(1 for ee in range(lo_end, Lb) for ff in range(ee + 1, Lb) if ee + ff >= thr)
        while hi > lo_end:
            lo = max(hi - kg, lo_end)
            D = R >> lo
            F = F + (D << lo)
            a = np.abs(F)
            bits = [(a >> ee) & 1 for ee in range(Lb + 2)]
            S = np.zeros_like(a)
            for ee in range(Lb + 2):
                if 2 * ee >= thr - 1:
                    S += bits[ee] << (2 * ee)
                for ff in range(ee + 1, Lb + 2):
                    if ee + ff >= thr:
                        S += (bits[ee] & bits[ff]) << (ee + ff + 1)
            Ff = F.astype(np.float64)
            Gt = Ff - 3.0 * S.astype(np.float64) / 2.0 ** 28 + Ff ** 3 / 2.0 ** 53 - 5.0 * Ff ** 4 / 2.0 ** 82
            R = y - np.round(Gt).astype(np.int64)
            groups.append([int(lo), int(hi), int(np.abs(D).max()).bit_length(), bits_signed(R)])
            hi = lo
        Ug = F + R
        grp_info = dict(kg=kg, groups=groups, pairs={thr: npairs}, lo_end=int(lo_end), thr=thr)
        qv = (Ug - U1).astype(np.float64)
    elif ho == "none":
        qv = np.zeros(len(d))
    elif ho == "ideal":
        qv = q_ho(-U1.astype(np.float64) / 2.0 ** 26)
    else:                                       # from the top t bits of U1 (two's complement, midpoint)
        bU = bits_signed(U1)
        sdrop = max(0, bU - t_ho)
        U1q = ((U1 >> sdrop) << sdrop) + ((1 << sdrop) >> 1)
        qv = q_ho(-U1q.astype(np.float64) / 2.0 ** 26)
    q_int = np.round(qv).astype(np.int64)
    U = c0 + U1 + q_int
    lim = 1 << (26 - m)
    n_clip = int(((U < -lim) | (U >= lim)).sum())
    U = np.clip(U, -lim, lim - 1)
    ue = U - cm["Us"]
    nl1 = nbits_vals(l1)
    nl2 = nbits_vals(l2)
    # is the leftover a function of the last digit (then it costs no wire)?
    def fn_of_last(l, last):
        pairs = np.unique(np.stack([last, l]), axis=1)
        return bool(len(np.unique(pairs[0])) == pairs.shape[1])
    info.update(m=m, xd=xd, F1=F1, a1=a1, n1=n1, e1=e1, a2=a2, n2=n2, e2=e2, g=g, dropped=dropped,
                b_delta=b_delta, b_theta=bits_signed(th), b_P=bits_signed(P), b_U=bits_signed(U),
                left1=list(nl1), left2=list(nl2), l1_fn_last=fn_of_last(l1, last1), l2_fn_last=fn_of_last(l2, last2),
                rows1=rows1, rows2=rows2, U_err=[float(ue.min()), float(ue.max())],
                q_max=float(qv.max()), q_bits=int(max(1, math.ceil(math.log2(max(1.0, float(np.abs(qv).max())) + 1)))),
                ho=ho, t=t_ho, c0=c0, n_clip=n_clip, grp=grp_info,
                cubic_max=float(np.ldexp(5 * np.abs(cm["w"]) ** 3 / 16, 27).max()))
    return U, info


# ------------------------------------------------------------------------------------------------ price
W12 = [13] * 15 + [11, 10, 9, 8, 7, 6, 5, 4, 3]     # main divider rows 12..35 of the line (bits) [S R32 table]
MAIN_FWD12, MAIN_CLEAN12 = 260.0, 309.0              # sweep2's reference m = 12 stages [M S1 / F1]
SHARED = dict(NORM=51.0, NORM_INV=43.25, ENDPOINT=114.0, PASSES=-8.0)


def price(info, q_hold=None):
    """stages (expected Toffolis) and the wire ledger. [P] rows from measured spans; [E] where marked."""
    m = info["m"]
    st = dict(SHARED)
    root = 96.0 * sum(i + 1 for i in range(1, m + 1)) / 90.0       # [E] in-place root, step i costs i + 1
    st["ROOT"] = root
    st["ROOT_INV"] = root
    if m <= 9:
        st["CENTRE"] = 2.0 * (m + 2)          # [E] delta centred on d = 2Q + 1: one (m+3)-bit row each way
    else:
        st["DELTA_LO"] = 22.0                 # [E sweep2] the dependent low bits of delta, both ways
    red = sum(w - min(w, m + 1) for w in W12)
    st["DIV_FWD"] = MAIN_FWD12 - red
    st["DIV_CLEAN"] = MAIN_CLEAN12 - (MAIN_CLEAN12 / MAIN_FWD12) * red
    # wires. An in-place division keeps its wire count: digits + what is left in the register.
    lb1 = 0 if info["l1_fn_last"] else info["left1"][1]
    lb2 = 0 if info["l2_fn_last"] else info["left2"][1]
    xd_w = info["dropped"]
    extras = lb1 + lb2 + xd_w
    w_div1 = m + 4 + max(r["live"] for r in info["rows1"]) + 2
    w_div2 = m + 4 + lb1 + max(r["live"] for r in info["rows2"]) + 2
    w_ho = m + 4 + info["n2"] + lb1 + lb2 + 3           # [E] the term accumulated in place (one guard wire)
    w_main = m + 34 + extras                             # [P] v m + flags 4 + dividend / quotient / remainder 29 + 1
    w_end = m + 32 + extras                              # [P] endpoint: 3 root wires retired as today (44 at m 12)
    q = max(w_div1, w_div2, w_main, w_ho, w_end)
    if q_hold is not None:
        q = max(q, q_hold)
    c1 = max(0, q - (w_div1 - 2))
    c2 = max(0, q - (w_div2 - 2))
    st["DIV1"] = float(sum(fwd_row(r["span"]) for r in info["rows1"]))
    st["DIV1_INV"] = float(sum(inv_row(r["span"], c1) for r in info["rows1"]))
    st["DIV2"] = float(sum(fwd_row(r["span"]) for r in info["rows2"]))
    st["DIV2_INV"] = float(sum(inv_row(r["span"], c2) for r in info["rows2"]))
    if info["ho"] != "none":
        t = info["t"] if info["ho"] == "top" else info["q_bits"] + 1
        sq = t * t / 2.0 + t                  # [E] truncated square of t bits, times 3/8, accumulated in place
        add = 27 - m                          # [P] its carry into U
        cb = 0.0
        if info["cubic_max"] > 0.5:           # [E] the cubic term: bits^2
            cb = float(int(math.ceil(math.log2(info["cubic_max"] + 1))) ** 2)
        st["HO"] = 2.0 * (sq + add + cb)      # [E] forward and back
    tof = sum(st.values())
    led = dict(div1=w_div1, div2=w_div2, ho=w_ho, main=w_main, end=w_end, left_bits=[lb1, lb2], xd=xd_w)
    return st, tof, q, led


def cost_of(q, tof, margin, n_rot=34, n_t=1):
    if margin is None or margin <= 0:
        return float("inf"), float("nan")
    tr = 0.56 * math.log2(n_rot * 2.0 ** 64 / (math.pi * margin)) + 5.3
    return q * (tof + (n_t + n_rot * tr) / 4.0), tr


# ------------------------------------------------------------------------------------------------ commands
_STATE = {}


def setup(wd):
    ch = CORE.load(wd)
    ph = CORE.phase(ch, wd)
    _STATE.update(ch=ch, ph=ph, cm={})
    return ch, ph


def run_one(prm, lp=True):
    ch, ph = _STATE["ch"], _STATE["ph"]
    m = int(prm["m"])
    if m not in _STATE["cm"]:
        _STATE["cm"][m] = chain(ch, m)
    cm = _STATE["cm"][m]
    U, info = word_cm(cm, prm)
    bt = int(prm.get("bt", m - 1))
    p = RT.divide_m(cm, U, b_t=bt)
    r = CORE.evaluate_word(ph, ch, p, lp=lp)
    st, tof, q, led = price(info)
    marg = r["m"] if r["m"] is not None else None
    c, tr = cost_of(q, tof, marg, n_rot=33)
    out = dict(prm=prm, m=m, q=q, tof=round(tof, 2), margin=float(marg), m_pa=float(r["m_pa"]), nfail=int(r["nfail"]),
               m_lp=(None if r["m_lp"] is None else float(r["m_lp"])), tper=tr, cost=c,
               gain_pct=(None if not math.isfinite(c) else round(100.0 * (1 - c / BASE["cost"]), 2)),
               ledger=led, stages={k: round(v, 2) for k, v in st.items()},
               n1=info["n1"], n2=info["n2"], F1=info["F1"], b_delta=info["b_delta"], b_U=info["b_U"],
               left1=info["left1"], left2=info["left2"], l1_fn_last=info["l1_fn_last"], l2_fn_last=info["l2_fn_last"],
               dropped=info["dropped"], g=info["g"], n_clip=info["n_clip"], U_err=info["U_err"], grp=info["grp"],
               q_max=info["q_max"],
               spans1=[r_["span"] for r_ in info["rows1"]], spans2=[r_["span"] for r_ in info["rows2"]])
    return out


def _worker(args):
    wd, prm = args
    if "ch" not in _STATE:
        setup(wd)
    try:
        return run_one(prm)
    except Exception as e:      # a set that does not build is a result too
        return dict(prm=prm, error=repr(e))


def grid(ms, mode="a"):
    """mode a: ideal higher-order term; F1 (appended zero bits) x where each division stops (r1, r2).
    mode b: sets read from the environment variable R33_SETS (JSON list), swept over the higher-order term."""
    sets = []
    if mode == "a":
        for m in ms:
            b_delta = {6: 16, 7: 15, 8: 14, 9: 13, 10: 14, 11: 15, 12: 16}[m]
            f0 = 27 - m - b_delta
            for F1 in (f0 - 1, f0, f0 + 1, f0 + 2):
                if F1 < 0:
                    continue
                for r1 in (-3, -2, -1, 0, 1):
                    for r2 in (-3, -2, -1, 0, 1):
                        sets.append(dict(m=m, F1=F1, r1=r1, r2=r2, ho="ideal"))
    else:
        base = json.loads(os.environ["R33_SETS"])
        for b in base:
            for ho, ts in (("none", (0,)), ("top", (4, 5, 6, 7, 8, 9, 10, 11, 12, 14))):
                for t in ts:
                    for c0 in ((0,) if ho != "none" else (0, 1, 2, 4, 8, 16, 32, 64, 128, 256)):
                        p = dict(b)
                        p.update(ho=ho, t=t, c0=c0)
                        sets.append(p)
    return sets


def main():
    wd, cmd = sys.argv[1], sys.argv[2]
    if cmd == "chain":
        ch, ph = setup(wd)
        for m in range(6, 13):
            cm = chain(ch, m)
            print(f"m {m}: d {int(cm['d'].min())}..{int(cm['d'].max())}  delta {int(cm['delta'].min())}.."
                  f"{int(cm['delta'].max())} ({bits_signed(cm['delta'])} bits)  |w| max {np.abs(cm['w']).max():.3e} "
                  f"(2^{math.log2(np.abs(cm['w']).max()):.2f})  U* {cm['Us'].min():.1f}..{cm['Us'].max():.1f}  "
                  f"q_ho max {q_ho(cm['w']).max():.1f}  cubic max {np.ldexp(5 * np.abs(cm['w']) ** 3 / 16, 27).max():.2f}",
                  flush=True)
    elif cmd == "one":
        setup(wd)
        prm = json.loads(sys.argv[3])
        out = run_one(prm, lp=(len(sys.argv) < 5 or sys.argv[4] != "nolp"))
        print("RESULT " + json.dumps(out), flush=True)
    elif cmd == "grid":
        out_path, nproc = sys.argv[3], int(sys.argv[4])
        ms = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else list(range(6, 13))
        mode = sys.argv[6] if len(sys.argv) > 6 else "a"
        if mode == "list":                       # an explicit list of sets (JSON file)
            sets = json.load(open(sys.argv[7]))
        else:
            sets = grid(ms, mode)
        import multiprocessing as mp
        with mp.Pool(nproc) as pool, open(out_path, "a") as f:
            for res in pool.imap_unordered(_worker, [(wd, s) for s in sets]):
                f.write(json.dumps(res) + "\n")
                f.flush()
        print(f"done {len(sets)} sets", flush=True)


if __name__ == "__main__":
    main()
