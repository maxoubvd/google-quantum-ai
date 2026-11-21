#!/usr/bin/env python3
"""r31_inv-hyb_model.py -- seat INV-HYB (R31 panel): exact models for research/panel/R31-INV-HYB.md.

Our own code. The public model JSON (tables C and U0 of the public 46-wire circuit, Syed Asad Ali, Apache-2.0) is
READ AS DATA from the path given on the command line (only by the commands that say so); no public code is
imported or run. The chain x -> (k, v, d, h) and the divider recurrence follow the public oracle's formulas as
read (native/candidate01/whole_model.py, parse-only reading).

Units: phases in 2^-64 turn as in the official verifier: r = isqrt(2^130 // x); target = (r + 1) >> 1;
tol = r >> 25. The phase word p of the public circuit is in 2^-37 turn at k = 0 (its slots step by 4 units).

  r31_inv-hyb_model.py CACHE chain                 build the reference and the root chain (all 3,133,452 inputs)
  r31_inv-hyb_model.py CACHE base MODEL.json       public word from its own tables: worst |err| / tol
  r31_inv-hyb_model.py CACHE ideal                 public divider fed with the ideal (rounded) product
  r31_inv-hyb_model.py CACHE hpoly [opts]          piecewise form on (v, h): least table bits that pass
  r31_inv-hyb_model.py CACHE dd                    table-free h path: two divisions by d, then the public divider
  r31_inv-hyb_model.py CACHE theta F               divisor form: quotient of 2^40 by d + theta (theta with F bits)
  r31_inv-hyb_model.py CACHE clite MODEL.json R    coefficient from a table of 2^R rows + interpolation (ideal rounding)
  r31_inv-hyb_model.py CACHE cq MODEL.json R       the same with quantised A and S, and the coefficient coarsened
  r31_inv-hyb_model.py CACHE tail MODEL.json       slack of the divider's truncated tail
  r31_inv-hyb_model.py CACHE tailfit MODEL.json 1 [seed]  the tail one row earlier, table words re-picked per segment
  r31_inv-hyb_model.py CACHE tailrows MODEL.json   one more bit dropped only in the last tail rows
  r31_inv-hyb_model.py CACHE tailfit2 MODEL.json 1 [span]  the same repair in the public table format exactly
  r31_inv-hyb_model.py CACHE tailfit2 MODEL.json 0 4 J     the public tail with one more bit dropped from row J on
  r31_inv-hyb_model.py CACHE variants MODEL.json [span]    four coefficient rows; product rows one bit shorter
  r31_inv-hyb_model.py CACHE prodrows MODEL.json MASK ..   the product rows of each mask one bit shorter
  r31_inv-hyb_model.py CACHE seedrows MODEL.json           two seed rows instead of three
"""
import json
import math
import os
import sys
import time

import numpy as np

MAX_X = 3133452
M64 = (1 << 64) - 1


def load_chain(cache):
    path = os.path.join(cache, "chain.npz")
    if os.path.exists(path):
        z = np.load(path)
        return {key: z[key] for key in z.files}
    t0 = time.time()
    tg = np.empty(MAX_X, dtype=np.uint64)
    tl = np.empty(MAX_X, dtype=np.uint64)
    for x in range(1, MAX_X + 1):
        r = math.isqrt((1 << 130) // x)
        tg[x - 1] = ((r + 1) >> 1) & M64
        tl[x - 1] = min(r >> 25, 1 << 63)
    x = np.arange(1, MAX_X + 1, dtype=np.int64)
    lead = np.zeros_like(x)
    for b in range(22):
        lead = np.where((x >> b) > 0, b, lead)
    k = (21 - lead) // 2
    n = x << (2 * k)
    Q = np.floor(np.sqrt((n << 4).astype(np.float64))).astype(np.int64)
    Q = np.where(Q * Q > (n << 4), Q - 1, Q)
    Q = np.where((Q + 1) * (Q + 1) <= (n << 4), Q + 1, Q)
    assert np.all(Q * Q <= n << 4) and np.all((Q + 1) * (Q + 1) > n << 4)
    v = Q - 4096
    d = 2 * Q + 1
    h = ((n << 6) - d * d) // 64
    rc = (v * (v + 1) // 2) & 7
    delta = 64 * h + 63 - 8 * rc
    assert np.array_equal(d * d + delta, n << 6), "delta identity"
    assert v.min() >= 0 and v.max() < 4096 and h.min() >= -512 and h.max() < 512
    out = dict(x=x, k=k, n=n, v=v, d=d, h=h, rc=rc, delta=delta, target=tg, tol=tl)
    np.savez(path, **out)
    print(f"chain built in {time.time() - t0:.0f} s", flush=True)
    return out


def public_word(d, U, k, tail_shift=0, from_row=0):
    """Divider + endpoint word (2^-37 turn), the public oracle's recurrence as read. U: signed 15-bit product.
    tail_shift = 1 starts the truncated tail one row earlier (every tail row adds one bit less of the root);
    from_row = J > 0 instead keeps the public tail and drops one more bit only in rows j >= J."""
    assert U.min() >= -16384 and U.max() < 16384
    product = U & 32767
    r = 4096 - d
    phase = np.full(len(d), -((1 << 37) - (1 << 29)), dtype=np.int64)
    for j in range(9, 36):
        sign = (r < 0).astype(np.int64)
        phase -= sign << (37 - j)
        if j <= 23:
            bit = (product >> (23 - j)) & 1
            if j == 9:
                bit = -bit
            raw = 2 * r + (2 * sign - 1) * d + bit
        elif j == 24:
            raw = 2 * r + (2 * sign - 1) * d
        else:
            if from_row:
                kk = (j - 26 if j >= 28 else 0) + (1 if j >= from_row else 0)
                if j < 28 and kk:
                    kk = 2
            else:
                kk = j - 26 + tail_shift if j >= 28 - tail_shift else 0
            V = d // 2
            A = ((((V >> (kk + 1)) + ((V >> kk) & 1)) << (kk + 1)) - (V & 1)) if kk else V
            u = (r - 1) // 2
            raw = 2 * (2 * u + sign + (2 * sign - 1) * A) + 1
        nextu = raw // 2
        assert np.all((nextu >= -8192) & (nextu < 8192)), j
        r = 2 * nextu + (raw & 1)
    phase -= (r < 0).astype(np.int64) << 1
    p = (phase - 1) & ((1 << 37) - 1)
    word = ((p << k) + (p >> (27 - k))) & ((1 << 37) - 1)
    return word


def ratio_of_word(word, c):
    e = ((word.astype(np.uint64) << np.uint64(27)) - c["target"]).view(np.int64)
    return np.abs(e.astype(np.float64)) / c["tol"].astype(np.float64)


def report(name, ratio):
    i = int(np.argmax(ratio))
    print(f"{name:58s} worst |err|/tol {ratio[i]:.5f} at x={i + 1}; failing {int((ratio >= 1).sum())}",
          flush=True)
    return float(ratio[i])


def ideal_U(c):
    """U* = 2^27 (d / sqrt(N) - 1), N = d^2 + delta: the product the divider would need (float)."""
    d = c["d"].astype(np.float64)
    N = (c["n"] << 6).astype(np.float64)
    return np.ldexp(d / np.sqrt(N) - 1.0, 27)


def public_product(c, model_path):
    m = json.load(open(model_path))
    C = np.array(m["C"], dtype=np.int64)[c["v"]]
    U0 = np.array(m["U0"], dtype=np.int64)[c["v"]]
    return C, U0, product_from_C(c, C, U0)


def product_from_C(c, C, U0, coarse=0):
    """coarse: bit i set = product row i (of rows 0..8) one bit shorter (its lowest bit cleared after
    rounding); 1 alone is kept as "all nine rows" for the variants command."""
    if coarse == 1:
        coarse = 0x1FF
    q = c["h"] + 512
    F = np.zeros_like(q)
    for i in range(10):
        shift = 8 - i
        if (coarse >> i) & 1 and shift > -1:
            sh = shift + 1
            term = ((C + (1 << (sh - 1))) >> sh) << 1
        elif shift > 0:
            term = (C + (1 << (shift - 1))) >> shift
        elif shift == 0:
            term = C
        else:
            term = 2 * C + 1
        F += (2 * ((q >> i) & 1) - 1) * term
    return U0 - F // 2


# ----------------------------------------------------------------------------------------------------------
def cmd_base(c, model_path):
    C, U0, U = public_product(c, model_path)
    report("public tables + public divider (nominal angles)", ratio_of_word(public_word(c["d"], U, c["k"]), c))
    Us = ideal_U(c)
    e = U - Us
    print(f"public product minus ideal product: min {e.min():.2f} max {e.max():.2f} (product units)")
    print(f"C range {C.min()}..{C.max()}; U range {U.min()}..{U.max()}")


def cmd_ideal(c):
    Us = ideal_U(c)
    for off in (0.0, 0.5, 1.0, 1.5, 2.0):
        U = np.clip(np.floor(Us + off).astype(np.int64), -16384, 16383)
        report(f"ideal product floor(U* + {off}) + public divider", ratio_of_word(public_word(c["d"], U, c["k"]), c))


# ----------------------------------------------------------------------------------------------------------
def cmd_dd(c):
    """Table-free h path. delta = N - d^2 (odd, 16 bits, on the h wires + 3 root-class bits + constants).
    Division 1: delta -> digits of theta = delta / d (weights 2, 1, 1/2, ...).  Division 2: theta -> P digits,
    P = 2^13 theta / d = 2^26 delta / d^2 scaled by 1/2 (the public product is U = -P up to a constant).
    Each division is the non-restoring recurrence reg -= sign(reg) * T_i with T_i = round(d / 2^i) (truncated
    divisor in the tail), which is reversible given the digits and the final register (the leftover)."""
    d, delta, k = c["d"], c["delta"], c["k"]
    # division 1: digits weights 2, 1, 1/2, ..., 2^-12 (15 rows)
    reg = delta.copy()
    theta_int = np.zeros_like(reg)          # theta * 2^12
    rowbits1 = 0
    exps = list(range(1, -13, -1))
    for e in exps:
        T = 2 * d if e == 1 else (d if e == 0 else (d + (1 << (-e - 1))) >> (-e))
        s = np.where(reg < 0, -1, 1)
        reg = reg - s * T
        theta_int += s << (e + 12)
        rowbits1 += int(T.max()).bit_length() + 1
    l1 = reg
    print(f"division 1: {len(exps)} rows, sum of row widths {rowbits1}; leftover range {l1.min()}..{l1.max()} "
          f"({len(np.unique(l1))} values)")
    th_err = theta_int.astype(np.float64) / 4096.0 - delta.astype(np.float64) / d.astype(np.float64)
    print(f"  theta error (units of 1): min {th_err.min():.5f} max {th_err.max():.5f}; |theta_int| max "
          f"{np.abs(theta_int).max()}")
    # division 2: reg2 = theta * 2^12 (odd integer), P = sum rho_i 2^(13-i), reg2 2^14 / d ... rows T = round(d / 2^i)
    for extra in (0, 1, 2):
        reg2 = theta_int << extra
        P = np.zeros_like(reg2)
        rowbits2 = 0
        nrows = 14 + extra
        for i in range(nrows):
            sh = i - extra
            if sh < 0:
                T = d << (-sh)
            elif sh == 0:
                T = d
            else:
                T = (d + (1 << (sh - 1))) >> sh
            s = np.where(reg2 < 0, -1, 1)
            reg2 = reg2 - s * T
            P += s << (nrows - 1 - i)
            rowbits2 += int(T.max()).bit_length() + 1
        l2 = reg2
        # P ~ theta_int 2^extra 2^(nrows-1) / (d 2^extra) = theta 2^12 2^13 / d = 2^25 delta / d^2  (extra = 0)
        scale = 2.0 ** (25 + extra)
        Pid = scale * delta.astype(np.float64) / (d.astype(np.float64) ** 2)
        pe = P - Pid
        print(f"division 2 (+{extra} low rows): {nrows} rows, sum of row widths {rowbits2}; leftover "
              f"{l2.min()}..{l2.max()} ({len(np.unique(l2))} values); P - ideal: {pe.min():.2f}..{pe.max():.2f}; "
              f"|P| max {np.abs(P).max()}")
        # feed the public divider: U = -(2 P / 2^extra) + c0, 15-bit; the public product unit is 2^26 delta/d^2
        Us = ideal_U(c)
        for c0 in (0, 1, 2):
            U = -((2 * P) >> extra) + c0
            U = np.clip(U, -16384, 16383)
            report(f"  DD product (+{extra} rows, offset {c0}) + public divider",
                   ratio_of_word(public_word(d, U, k), c))
        ue = (-(2.0 * P) / 2 ** extra) - Us
        print(f"  DD product minus ideal product: {ue.min():.2f}..{ue.max():.2f} (product units)")


# ----------------------------------------------------------------------------------------------------------
def cmd_theta(c, F):
    """Divisor form: word = 2^40 / (d + theta), theta = delta / (2 d) with F fraction bits (first-order root).
    Exact quotient (float, error < 2^-10 unit) against the tolerance: what the arithmetic form can reach before
    row truncation."""
    d = c["d"].astype(np.float64)
    delta = c["delta"].astype(np.float64)
    th = np.round(delta / (2 * d) * 2.0 ** F) / 2.0 ** F
    word = 2.0 ** 40 / (d + th)                       # units 2^-37 turn at k = 0
    ph = word * np.ldexp(1.0, 27) * (2.0 ** c["k"])   # 2^-64 turn
    tf = c["target"].astype(np.float64)
    tl = c["tol"].astype(np.float64)

    def rat(phase):
        e = np.mod(phase - tf + 2.0 ** 63, 2.0 ** 64) - 2.0 ** 63
        return np.abs(e) / tl

    report(f"divisor form, theta with {F} fraction bits, exact quotient", rat(ph))
    for cc in (0.25, 0.5, 0.75):
        # second-order term of the root: sqrt(d^2 + delta) = d + delta/(2d) - delta^2/(8 d^3) + ...
        word2 = 2.0 ** 40 / (d + th - cc * (2.0 / d))
        report(f"  same with divisor offset -{cc} * 2/d", rat(word2 * np.ldexp(1.0, 27) * (2.0 ** c["k"])))
    tolw = c["tol"].astype(np.float64) / (np.ldexp(1.0, 27) * 2.0 ** c["k"])
    print(f"tolerance in word units (2^-37 turn at k = 0): {tolw.min():.3f}..{tolw.max():.3f}; one slot = 4 units")


# ----------------------------------------------------------------------------------------------------------
def seg_fit(a, cols, y, nseg, wts=None):
    """Least squares per segment a: y ~ cols @ coef[a]. Returns coef (nseg, ncol)."""
    nc = len(cols)
    A = np.zeros((nseg, nc, nc))
    b = np.zeros((nseg, nc))
    for i in range(nc):
        ci = cols[i] if wts is None else cols[i] * wts
        b[:, i] = np.bincount(a, weights=ci * y, minlength=nseg)
        for j in range(i, nc):
            s = np.bincount(a, weights=ci * cols[j], minlength=nseg)
            A[:, i, j] = s
            A[:, j, i] = s
    coef = np.zeros((nseg, nc))
    for s in range(nseg):
        if A[s, 0, 0] > 0:
            coef[s] = np.linalg.lstsq(A[s], b[s], rcond=None)[0]
    return coef


def cmd_hpoly(c, opts):
    """Piecewise form on (v, h): word(a, t, delta) = A0 + A1 t + A2 t^2 - delta (B0 + B1 t), a = v >> (12 - R),
    t = centred low root bits, delta = N - d^2. Word in slot units (2^-35 turn at k = 0; the public circuit's
    slots). Tables fitted by weighted least squares on the full domain, then rounded to g fraction bits each;
    products exact (no row truncation), so the table widths found are lower bounds for this form."""
    R = int(opts.get("R", 7))
    lowbits = 12 - R
    nseg = 1 << R
    v = c["v"]
    a = (v >> lowbits).astype(np.int64)
    t = ((v & ((1 << lowbits) - 1)) - ((1 << lowbits) - 1) / 2.0).astype(np.float64)
    delta = c["delta"].astype(np.float64)
    N = (c["n"] << 6).astype(np.float64)
    y = 2.0 ** 38 / np.sqrt(N)                     # slot units: 2^-35 turn at k = 0
    tol = y * 2.0 ** -24
    w = 1.0 / tol ** 2
    forms = {
        "q": [np.ones_like(t), t, t * t, -delta, -delta * t],
        "q+": [np.ones_like(t), t, t * t, -delta, -delta * t, -delta * t * t],
        "c": [np.ones_like(t), t, t * t, t ** 3, -delta, -delta * t],
        "lin": [np.ones_like(t), t, -delta, -delta * t],
    }
    names = {"q": ["A0", "A1", "A2", "B0", "B1"], "q+": ["A0", "A1", "A2", "B0", "B1", "B2"],
             "c": ["A0", "A1", "A2", "A3", "B0", "B1"], "lin": ["A0", "A1", "B0", "B1"]}
    form = opts.get("form", "q")
    cols = forms[form]
    coef = seg_fit(a, cols, y, nseg, w)
    fit = sum(coef[a, i] * cols[i] for i in range(len(cols)))
    ratio = np.abs(fit - y) / tol
    print(f"R = {R} ({nseg} rows), form {form}: unrounded fit worst |err|/tol {ratio.max():.4f}")
    # greedy: for each coefficient the fewest fraction bits g such that the total worst ratio stays < limit
    limit = float(opts.get("limit", 0.90))
    scale_of = {}
    for i, nm in enumerate(names[form]):
        scale_of[nm] = float(np.abs(cols[i]).max())
    g = {nm: 40 for nm in names[form]}

    def worst(gg):
        acc = np.zeros_like(y)
        for i, nm in enumerate(names[form]):
            q = np.round(coef[:, i] * 2.0 ** gg[nm]) / 2.0 ** gg[nm]
            acc += q[a] * cols[i]
        # the word itself is rounded to slots
        return float((np.abs(np.round(acc) - y) / tol).max())

    base = worst(g)
    print(f"  with the word rounded to slots, tables unrounded: {base:.4f}")
    if base >= limit:
        print("  form does not pass at this limit; no widths reported")
        return
    for nm in names[form]:
        lo = -8
        while True:
            trial = dict(g)
            trial[nm] = lo
            if worst(trial) < limit:
                break
            lo += 1
        g[nm] = lo
    total = 0
    print(f"  limit {limit}: fraction bits and widths (bits = ceil(log2(range)) + fraction bits):")
    for i, nm in enumerate(names[form]):
        used = np.bincount(a, minlength=nseg) > 0
        cv = coef[used, i]
        rng = float(cv.max() - cv.min())
        mag = float(np.abs(cv).max())
        bits_rng = max(0, math.ceil(math.log2(rng * 2.0 ** g[nm] + 1)))
        bits_mag = max(0, math.ceil(math.log2(mag * 2.0 ** g[nm] + 1)))
        total += bits_mag
        print(f"    {nm}: fraction bits {g[nm]:3d}; |value| max {mag:.6g}; width {bits_mag} bits "
              f"(span over the table {bits_rng} bits); multiplicand range of its column {scale_of[nm]:.6g}")
    print(f"  final worst ratio {worst(g):.4f}; table bits in all {total} per row, {nseg} rows")


# ----------------------------------------------------------------------------------------------------------
def cmd_clite(c, model_path, R):
    """Coefficient C(v) of the public product from a table of 2^R rows + linear interpolation (least squares on
    the public C table, read as data), rounded to integers; everything else as the public circuit."""
    m = json.load(open(model_path))
    Ct = np.array(m["C"], dtype=np.float64)
    U0t = np.array(m["U0"], dtype=np.int64)
    vv = np.arange(4096)
    low = 12 - R
    a = vv >> low
    t = (vv & ((1 << low) - 1)) - ((1 << low) - 1) / 2.0
    for deg in (1, 2):
        cols = [np.ones(4096), t] + ([t * t] if deg == 2 else [])
        coef = seg_fit(a, cols, Ct, 1 << R)
        fit = sum(coef[a, i] * cols[i] for i in range(len(cols)))
        Cn = np.round(fit).astype(np.int64)
        dC = Cn - Ct
        C = Cn[c["v"]]
        U = product_from_C(c, C, U0t[c["v"]])
        U = np.clip(U, -16384, 16383)
        tag = f"C from {1 << R} rows, degree {deg} in the low {low} root bits"
        print(f"{tag}: C_new - C_public in {dC.min():.0f}..{dC.max():.0f}")
        report("  " + tag + " + public divider", ratio_of_word(public_word(c["d"], U, c["k"]), c))


def cmd_cq(c, model_path, R):
    """As clite, but with the table values quantised as the public format does: D = A + S t / 8 with integer
    A and a slope S of g fraction bits per root step, C = 2 D; and with the coefficient coarsened (low bits
    dropped) to see what the product tolerates. All checked on every input."""
    m = json.load(open(model_path))
    Ct = np.array(m["C"], dtype=np.float64)
    U0t = np.array(m["U0"], dtype=np.int64)
    vv = np.arange(4096)
    low = 12 - R
    a = vv >> low
    t = (vv & ((1 << low) - 1)) - ((1 << low) - 1) / 2.0
    cols = [np.ones(4096), t]
    coef = seg_fit(a, cols, Ct / 2.0, 1 << R)       # D = C / 2
    for g in (3, 4, 5):
        for half in (0, 1):
            A = np.round(coef[:, 0] * (1 << half)) / (1 << half)
            S = np.round(coef[:, 1] * (1 << g)) / (1 << g)
            Dn = np.floor(A[a] + S[a] * t + 0.5).astype(np.int64)
            dD = Dn - Ct / 2.0
            U = np.clip(product_from_C(c, 2 * Dn[c["v"]], U0t[c["v"]]), -16384, 16383)
            sb = int(np.ceil(np.log2(np.abs(S).max() * (1 << g) + 1)))
            tag = (f"{1 << R} rows, D = A + S t, A with {half} fraction bit(s), S with {g} fraction bits "
                   f"({sb} bits)")
            print(f"{tag}: D_new - D_public in {dD.min():.1f}..{dD.max():.1f}")
            report("  " + tag, ratio_of_word(public_word(c["d"], U, c["k"]), c))
    if R == 7:
        C0 = np.array(m["C"], dtype=np.int64)
        for drop in (2, 3, 4):
            Cc = ((C0 + (1 << (drop - 1))) >> drop) << drop
            U = np.clip(product_from_C(c, Cc[c["v"]], U0t[c["v"]]), -16384, 16383)
            report(f"public C with its low {drop} bits rounded away (D with {14 - drop} bits)",
                   ratio_of_word(public_word(c["d"], U, c["k"]), c))


def cmd_tail(c, model_path):
    """Is there slack in the divider's truncated tail? Public tables and the ideal product, with the tail
    starting one and two rows earlier (each tail row then adds one or two bits less)."""
    C, U0, U = public_product(c, model_path)
    Ui = np.clip(np.floor(ideal_U(c)).astype(np.int64), -16384, 16383)
    for ts in (0, 1, 2):
        try:
            report(f"public product, tail shifted by {ts}", ratio_of_word(public_word(c["d"], U, c["k"], ts), c))
            report(f"ideal product,  tail shifted by {ts}", ratio_of_word(public_word(c["d"], Ui, c["k"], ts), c))
        except AssertionError as err:
            print(f"tail shifted by {ts}: remainder leaves its 14 bits at row {err}")


def seed_from_C(C, rc, offset=-1):
    """The public seed word as a function of the coefficient and the root class (formula as read)."""
    seed = np.full(len(C), offset, dtype=np.int64)
    for i in range(3):
        seed += (2 * ((rc >> i) & 1) - 1) * ((C + (1 << (11 - i))) >> (12 - i))
    return seed


def cmd_tailfit(c, model_path, ts, so=-1):
    """Tail started ts rows earlier: can a change of the table words A and S (the segment's coefficient line
    moved by an integer and tilted by one or two units of its slope) repair the failing inputs? Checked on
    every input of every touched segment; the seed follows its formula."""
    m = json.load(open(model_path))
    Ct = np.array(m["C"], dtype=np.int64)
    U0t = np.array(m["U0"], dtype=np.int64)
    v, rc = c["v"], c["rc"]
    assert np.array_equal(seed_from_C(Ct[v], rc), U0t[v]), "seed formula"
    U = np.clip(product_from_C(c, Ct[v], seed_from_C(Ct[v], rc, so)), -16384, 16383)
    ratio = ratio_of_word(public_word(c["d"], U, c["k"], ts), c)
    report(f"public tables, seed constant {so}, tail shifted by {ts}", ratio)
    a = v >> 5
    bad_segs = sorted(set(a[ratio >= 1].tolist()))
    print(f"failing inputs {int((ratio >= 1).sum())} in segments {bad_segs}")
    worst_seg = np.zeros(128)
    np.maximum.at(worst_seg, a, ratio)
    fixed = 0
    newC = Ct.copy()
    for sgm in bad_segs:
        idx = np.nonzero(a == sgm)[0]
        sub = {key: c[key][idx] for key in ("d", "k", "h", "target", "tol")}
        best = None
        tt = (v[idx] & 31) - 15.5
        for dS in (0, 1, -1, 2, -2, 3, -3):
            for dA in (0, 1, -1, 2, -2, 3, -3):
                Cs = Ct[v[idx]] + 2 * (dA + np.floor(dS * tt / 16.0 + 0.5).astype(np.int64))
                Us = np.clip(product_from_C(sub, Cs, seed_from_C(Cs, rc[idx], so)), -16384, 16383)
                r = float(ratio_of_word(public_word(sub["d"], Us, sub["k"], ts), sub).max())
                if best is None or r < best[0]:
                    best = (r, dA, dS)
        print(f"  segment {sgm}: {len(idx)} inputs, worst before {worst_seg[sgm]:.4f}; best change A {best[1]:+d}, "
              f"S {best[2]:+d} gives {best[0]:.4f}")
        if best[0] < 1:
            fixed += 1
            lo = sgm << 5
            w32 = np.arange(32) - 15.5
            newC[lo:lo + 32] += 2 * (best[1] + np.floor(best[2] * w32 / 16.0 + 0.5).astype(np.int64))
    print(f"segments repaired by A and S: {fixed} of {len(bad_segs)}")
    U2 = np.clip(product_from_C(c, newC[v], seed_from_C(newC[v], rc, so)), -16384, 16383)
    report(f"all inputs, repaired table, tail shifted by {ts}", ratio_of_word(public_word(c["d"], U2, c["k"], ts), c))


def cmd_tailrows(c, model_path):
    """Public tables; one more root bit dropped only in the tail rows j >= J (J = 28 .. 35)."""
    C, U0, U = public_product(c, model_path)
    for J in range(28, 36):
        saved = 36 - J
        report(f"public tables, rows {J}..35 one bit shorter ({saved} row bits each way)",
               ratio_of_word(public_word(c["d"], U, c["k"], 0, J), c))


def D_from_AS(A, S, w, first_row=0):
    """The public coefficient rows as read (lifecycle.coefficient_forward): five signed, separately rounded
    additions of the 6-bit slope S into the 13-bit table word A, one per low root bit. first_row = 1 leaves
    out the row of the lowest root bit."""
    D = np.array(A, dtype=np.int64) + 0 * w
    for i in range(first_row, 5):
        kk = 4 - i
        bit = (w >> i) & 1
        up = (S >> (kk + 1)) + ((S >> kk) & 1)
        dn = (S >> (kk + 1)) + (((S >> (kk - 1)) & 1) if kk else 0)
        D = D + np.where(bit == 0, up, -dn)
    return D


def cmd_tailfit2(c, model_path, ts, span=4, J=0):
    """As tailfit, in the public table format exactly: A (13 bits) and the stored slope (5 bits; the slope's
    low bit stays the implicit one) re-picked for the segments that fail when the tail starts ts rows earlier."""
    m = json.load(open(model_path))
    A0 = np.array(m["A"], dtype=np.int64)
    S0 = np.array(m["S"], dtype=np.int64)
    Dt = np.array(m["D"], dtype=np.int64)
    vv = np.arange(4096)
    assert np.array_equal(D_from_AS(A0[vv >> 5], S0[vv >> 5], vv & 31), Dt), "coefficient rows as read"
    # the slope's low bit is a wire the circuit already has (the model's implicit slope bit): keep its value
    print("coefficient formula reproduces all 4,096 public values of D")
    v, rc = c["v"], c["rc"]
    a = v >> 5
    C = 2 * Dt[v]
    ratio = ratio_of_word(public_word(c["d"], np.clip(product_from_C(c, C, seed_from_C(C, rc)), -16384, 16383),
                                      c["k"], ts, J), c)
    report(f"public tables, tail shifted by {ts}, extra bit from row {J}", ratio)
    bad = sorted(set(a[ratio >= 1].tolist()))
    newD = Dt.copy()
    fixed = 0
    for sgm in bad:
        idx = np.nonzero(a == sgm)[0]
        sub = {key: c[key][idx] for key in ("d", "k", "h", "target", "tol")}
        w = v[idx] & 31
        best = None
        for dH in range(-span, span + 1):
            S1 = int(S0[sgm]) + 2 * dH
            if not 0 <= S1 < 64:
                continue
            for dA in range(-span, span + 1):
                A1 = int(A0[sgm]) + dA
                if not 0 <= A1 < 8192:
                    continue
                Cs = 2 * D_from_AS(A1, S1, w)
                Us = np.clip(product_from_C(sub, Cs, seed_from_C(Cs, rc[idx])), -16384, 16383)
                r = float(ratio_of_word(public_word(sub["d"], Us, sub["k"], ts, J), sub).max())
                if best is None or r < best[0]:
                    best = (r, A1, S1)
        print(f"  segment {sgm}: A {int(A0[sgm])} -> {best[1]}, S {int(S0[sgm])} -> {best[2]}: worst {best[0]:.4f}")
        if best[0] < 1:
            fixed += 1
            newD[sgm << 5:(sgm + 1) << 5] = D_from_AS(best[1], best[2], np.arange(32))
    print(f"segments repaired in the public format: {fixed} of {len(bad)}")
    C2 = 2 * newD[v]
    report(f"all inputs, re-picked table, tail shifted by {ts}, extra bit from row {J}",
           ratio_of_word(public_word(c["d"], np.clip(product_from_C(c, C2, seed_from_C(C2, rc)), -16384, 16383),
                                     c["k"], ts, J), c))


def cmd_variants(c, model_path, span=4):
    """Two cheaper variants of the public coefficient and product, each with A and the stored slope re-picked
    per segment in the public format (all 128 segments): (1) four coefficient rows instead of five (the lowest
    root bit not interpolated); (2) every product row one bit shorter."""
    m = json.load(open(model_path))
    A0 = np.array(m["A"], dtype=np.int64)
    S0 = np.array(m["S"], dtype=np.int64)
    v, rc = c["v"], c["rc"]
    a = v >> 5
    order = np.argsort(a, kind="stable")
    bounds = np.searchsorted(a[order], np.arange(129))
    for name, first_row, coarse in (("four coefficient rows", 1, 0), ("product rows one bit shorter", 0, 1)):
        worst = 0.0
        failing_segments = 0
        for sgm in range(128):
            idx = order[bounds[sgm]:bounds[sgm + 1]]
            if len(idx) == 0:
                continue
            sub = {key: c[key][idx] for key in ("d", "k", "h", "target", "tol")}
            w = v[idx] & 31
            best = None
            for dH in range(-span, span + 1):
                S1 = int(S0[sgm]) + 2 * dH
                if not 0 <= S1 < 64:
                    continue
                for dA in range(-span, span + 1):
                    A1 = int(A0[sgm]) + dA
                    Cs = 2 * D_from_AS(A1, S1, w, first_row)
                    Us = np.clip(product_from_C(sub, Cs, seed_from_C(Cs, rc[idx]), coarse), -16384, 16383)
                    r = float(ratio_of_word(public_word(sub["d"], Us, sub["k"]), sub).max())
                    if best is None or r < best:
                        best = r
            worst = max(worst, best)
            failing_segments += best >= 1
        print(f"{name}: best re-picked table, worst |err|/tol {worst:.4f}; segments that cannot pass: "
              f"{failing_segments} of 128", flush=True)


def cmd_prodrows(c, model_path, masks, span=4):
    """Product rows of the mask one bit shorter; public tables first, then A and the stored slope re-picked
    (public format) for the segments that fail."""
    m = json.load(open(model_path))
    A0 = np.array(m["A"], dtype=np.int64)
    S0 = np.array(m["S"], dtype=np.int64)
    Dt = np.array(m["D"], dtype=np.int64)
    v, rc = c["v"], c["rc"]
    a = v >> 5
    C = 2 * Dt[v]
    for mask in masks:
        rows = [i for i in range(9) if (mask >> i) & 1]
        U = np.clip(product_from_C(c, C, seed_from_C(C, rc), mask), -16384, 16383)
        ratio = ratio_of_word(public_word(c["d"], U, c["k"]), c)
        bad = sorted(set(a[ratio >= 1].tolist()))
        left = []
        worst = float(ratio[np.isin(a, bad, invert=True)].max()) if len(bad) < 128 else 0.0
        for sgm in bad:
            idx = np.nonzero(a == sgm)[0]
            sub = {key: c[key][idx] for key in ("d", "k", "h", "target", "tol")}
            w = v[idx] & 31
            best = 9.0
            for dH in range(-span, span + 1):
                S1 = int(S0[sgm]) + 2 * dH
                if not 0 <= S1 < 64:
                    continue
                for dA in range(-span, span + 1):
                    Cs = 2 * D_from_AS(int(A0[sgm]) + dA, S1, w)
                    Us = np.clip(product_from_C(sub, Cs, seed_from_C(Cs, rc[idx]), mask), -16384, 16383)
                    best = min(best, float(ratio_of_word(public_word(sub["d"], Us, sub["k"]), sub).max()))
            worst = max(worst, best)
            if best >= 1:
                left.append((sgm, round(best, 4)))
        print(f"product rows {rows} one bit shorter: public tables fail on {int((ratio >= 1).sum())} inputs in "
              f"{len(bad)} segments; after re-picking A, S: worst {worst:.4f}, segments left {left}", flush=True)


def cmd_seedrows(c, model_path):
    """Two seed rows instead of three (public tables, the seed's constant re-picked among -2, -1, 0)."""
    m = json.load(open(model_path))
    Dt = np.array(m["D"], dtype=np.int64)
    v, rc = c["v"], c["rc"]
    C = 2 * Dt[v]
    for rows in ((1, 2), (0, 2), (0, 1)):
        for off in (-2, -1, 0):
            sd = np.full(len(C), off, dtype=np.int64)
            for i in rows:
                sd += (2 * ((rc >> i) & 1) - 1) * ((C + (1 << (11 - i))) >> (12 - i))
            U = np.clip(product_from_C(c, C, sd), -16384, 16383)
            report(f"seed rows {rows}, constant {off}", ratio_of_word(public_word(c["d"], U, c["k"]), c))


def main():
    cache, cmd = sys.argv[1], sys.argv[2]
    os.makedirs(cache, exist_ok=True)
    c = load_chain(cache)
    if cmd == "chain":
        print("inputs", len(c["x"]), "v max", int(c["v"].max()), "k max", int(c["k"].max()))
    elif cmd == "base":
        cmd_base(c, sys.argv[3])
    elif cmd == "ideal":
        cmd_ideal(c)
    elif cmd == "dd":
        cmd_dd(c)
    elif cmd == "theta":
        cmd_theta(c, int(sys.argv[3]))
    elif cmd == "hpoly":
        opts = dict(s.split("=", 1) for s in sys.argv[3:])
        cmd_hpoly(c, opts)
    elif cmd == "tailfit2":
        cmd_tailfit2(c, sys.argv[3], int(sys.argv[4]), int(sys.argv[5]) if len(sys.argv) > 5 else 4,
                     int(sys.argv[6]) if len(sys.argv) > 6 else 0)
    elif cmd == "seedrows":
        cmd_seedrows(c, sys.argv[3])
    elif cmd == "prodrows":
        cmd_prodrows(c, sys.argv[3], [int(x, 0) for x in sys.argv[4:]])
    elif cmd == "variants":
        cmd_variants(c, sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 4)
    elif cmd == "tailrows":
        cmd_tailrows(c, sys.argv[3])
    elif cmd == "tailfit":
        cmd_tailfit(c, sys.argv[3], int(sys.argv[4]), int(sys.argv[5]) if len(sys.argv) > 5 else -1)
    elif cmd == "tail":
        cmd_tail(c, sys.argv[3])
    elif cmd == "cq":
        cmd_cq(c, sys.argv[3], int(sys.argv[4]))
    elif cmd == "clite":
        cmd_clite(c, sys.argv[3], int(sys.argv[4]))
    else:
        raise SystemExit("unknown command")


if __name__ == "__main__":
    main()
