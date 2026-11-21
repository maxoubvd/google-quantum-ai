"""invsqrt_gen: bit-exact numeric model of the pipeline our generator builds, and the table fit.

Own design (workstream invsqrt_gen). Public ideas used: piecewise-quadratic evaluation with a hi/lo split of the
phase word (structure of the best public invsqrt circuit, read in research/ws_invsqrt.md), add-subtract rows
(Litinski, arXiv:2410.00899), coefficient fitting against the exact integer pipeline (Gu et al.,
arXiv:2610.22367). No table, angle or gate list of any public circuit is used: every number below is computed
from 1/sqrt(x) and the tolerance 2^-24/sqrt(x).

Definitions (x in [1, MAX_X]):
  e = floor(log2 x), p = e & 1, j = e >> 1, s = 21 - e, xn = x << s, a7 = xn[14..20], t = xn[0..13],
  a8 = 128 p + a7, F = 2^(-p/2) / sqrt(xn / 2^21), exact phase = 2^-j F turns.
  "reg unit" at level j = 2^-(26 + glo + j) turn; exact phase = 2^(26+glo) F reg units; tolerance 2^(2+glo) F.
Pipeline (integers; P = Params, sigma_i = +1 if t_i else -1, rr(v, k) = (v + 2^(k-1)) >> k):
  u  = C1[a8];  for i in P.a_rows:  u  -= sigma_i * rr(C2[a8], P.a_shift(i))
  lo = C0[a8];  for i in P.b_rows:  lo -= sigma_i * rr(u, P.b_shift(i))
  Wreg = H[a8] * 2^(wlo-1) + lo   (lo in [0, 2^wlo), u in [0, 2^wu), wlo = 19 + glo, wu = 17 + gu)
  bit k of lo sits at position 25 + glo + j - k (weight 2^-position half turns), bit h of H at position
  7 + j - h; bits at positions > top_pos are not delivered.
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import os

import numpy as np

MAX_X = 3 * (2 ** 10 - 2) ** 2  # 3,133,452
WS = os.path.expanduser("~/qcc/ws/invsqrt_gen")
C2_FIXED = None                 # experiments only: list of 256 curvature values that the fit must use as they are


@dataclasses.dataclass(frozen=True)
class Params:
    nt: int = 14          # bits of t
    glo: int = 1          # guard bits of lo (lo has 19 + glo wires)
    gu: int = 0           # guard bits of u (u has 17 + gu wires)
    fc2: int = 0          # extra fractional bits of the curvature table (c2 has about 8 + fc2 bits)
    wh: int = 8           # hi word bits
    a_lo: int = 6         # A uses t bits nt-1 .. a_lo
    top_pos: int = 35     # deepest rotation position (half-turn exponent)
    c1_span: int = 2      # fit search half-widths, in units of the 17-bit / 8-bit LSB
    c2_span: int = 2
    w_span: int = 3
    want: int = 32        # per-segment margin target of the fit, in 1/64 reg units
    virt: int = 0         # 1: bit 0 of lo is virtual (C0 even, never delivered, no wire; it equals a bit of u)
    tie_c1: int = 0       # 1: bit 0 of C1 equals bit a_shift(a_lo) of C2, so the guard bit of u is 0 after A
    fit_rev: int = 3
    tie_t0: int = 0       # 1: window form of stage A (invsqrt_gen_g6.py): bit 1 of C0 equals bit 0 of the top
                          #    table T = (C1 >> win) - P0, P0 = 1 - bit (win-1) of C1, and the A rows cannot
                          #    leave the window: sum_k rr(C2, k) <= 2^(win-1)
    win: int = 11         # low bits of the 18-bit u (bit 0 included) that the A rows work on
    c2_even: int = 0      # 1: C2 is even, so its bit 0 needs no wire (one table output less in two lookups)
    tie_mask: int = 2     # tie_t0: the XOR of the C0 bits in this mask equals bit 0 of the top table (2 = bit 1)
    b_free: int = 0       # K > 0: B rows with shift k >= K use the "free carry-in" of the rows module
                          #        (invsqrt_rows_b.rows_model_R21): lo -= sigma * 2 * rr(u >> 1, k), i.e. the row
                          #        is rounded to the 20 lo wires instead of to 21 bits (no AND for its carry-in)
    neg: int = 0          # 1: "negated word" form (idea and module side: workstream invsqrt_ph, 40 rotations).
                          #    The arithmetic computes Hc 2^(wlo-1) - (phase word): lo = C0 + sum sigma_i rr(u, k)
                          #    (the B rows run on complemented digit wires), C0 = (Hc 2^(wlo-1) - W) mod 2^wlo,
                          #    Hc = ceil(max over t of the phase word / 2^(wlo-1)); the table "H" then holds
                          #    N = 2^wh - Hc (wh - 1 wires). Delivered: -(N bit h at position 7 + j - h), -(1 at
                          #    positions 0..j-1), -(lo bit k at position 25 + glo + j - k); undelivered low bits
                          #    of lo round the phase up. The top lo wire is 0 at t = 0, so x in 64..255 needs no
                          #    rotation at position 10.
    zone: int = 0         # 1: zone format "Z9c" (format, fit and tables: workstream invsqrt_tab; hand-over note in
                          #    research/ws_invsqrt_tab.md). Where 1/sqrt is flat (zone 1: a7 >= 76 for even
                          #    exponents, a7 >= 64 for odd ones) two neighbouring segments share one table entry:
                          #    the address is a8 with its lowest bit cleared and the digit word is
                          #    (a7 bit 0) 2^13 + (t >> 1); everything else is the pipeline above (used with
                          #    a_lo = 5, nine A rows). 198 of the 256 addresses occur. The tables come from a
                          #    file (tabs=<name>); our own fit below does not know the zones.

    @property
    def wlo(self):
        return 19 + self.glo

    @property
    def wu(self):
        return 17 + self.gu

    @property
    def wc2(self):
        return 8 + self.fc2

    @property
    def a_rows(self):
        return list(range(self.nt - 1, self.a_lo - 1, -1))

    @property
    def b_rows(self):
        return list(range(self.nt - 1, -1, -1))

    def a_shift(self, i):
        return self.nt - 1 - i + self.fc2 - self.gu

    def b_shift(self, i):
        return self.nt - 1 - i + self.gu - (self.glo - 1)

    def key(self):
        d = dataclasses.asdict(self)
        if not d["tie_t0"]:                  # fields added later are left out at their default, so that the keys
            d.pop("tie_t0")                  # (and cached tables) of earlier builds do not change
            d.pop("win")
        if not d["c2_even"]:
            d.pop("c2_even")
        if d["tie_mask"] == 2:
            d.pop("tie_mask")
        if not d["b_free"]:
            d.pop("b_free")
        if not d["neg"]:
            d.pop("neg")
        if not d["zone"]:
            d.pop("zone")
        return hashlib.sha256(json.dumps(d, sort_keys=True).encode()).hexdigest()[:12]

    def top_table(self, c1):
        """Window form of u: (P0, T) with C1 = (C1 mod 2^win) + 2^win (P0 + T)."""
        p0 = 1 - ((c1 >> (self.win - 1)) & 1)
        return p0, (c1 >> self.win) - p0


def rr(v, k):
    """Rounded right shift: floor(v / 2^k) + bit (k-1) of v."""
    if k == 0:
        return v
    return (v + (1 << (k - 1))) >> k


def rows(P: Params, c1, c2, t):
    """Returns (u, R): u after A and R = lo - C0 after B, for arrays (broadcast) of c1, c2, t (int64)."""
    u = c1 + np.zeros_like(t)
    for i in P.a_rows:
        v = rr(c2, P.a_shift(i))
        u = np.where((t >> i) & 1, u - v, u + v)
    R = np.zeros_like(u)
    for i in P.b_rows:
        k = P.b_shift(i)
        if P.b_free and k >= max(P.b_free, 1):
            v = 2 * ((u >> (k + 1)) + ((u >> k) & 1))
        else:
            v = rr(u, k)
        R = np.where((t >> i) & 1, R - v, R + v)
    return u, R


def decompose(x):
    """Vectorised e, p, j, s, a8, t for an array of x (int64)."""
    x = np.asarray(x, dtype=np.int64)
    e = np.floor(np.log2(x.astype(np.float64))).astype(np.int64)
    e += ((np.int64(1) << (e + 1)) <= x)          # guard against float error
    e -= ((np.int64(1) << e) > x)
    s = 21 - e
    xn = x << s
    p = e & 1
    return e, p, e >> 1, s, 128 * p + ((xn >> 14) & 127), xn & 0x3FFF


def zone1(p, a7):
    """Zone flag of the format Z9c for parity p and the 7 address bits a7 (ints or arrays): 1 where a7 >= 64
    (odd exponents) or a7 >= 76 (even exponents), i.e. a6 AND (p OR a5 OR a4 OR (a3 AND a2))."""
    return ((a7 >> 6) & 1) & (p | (a7 >> 5) | (a7 >> 4) | ((a7 >> 3) & (a7 >> 2))) & 1


def zone_used():
    """The addresses that occur in the zone format (an odd address inside zone 1 never does)."""
    return [a for a in range(256) if not (zone1(a >> 7, a & 127) and (a & 1))]


def zone_apply(a8, t):
    """(address, digit word) of the zone format for arrays a8, t of the plain format."""
    z = zone1(a8 >> 7, a8 & 127) == 1
    return np.where(z, a8 & ~np.int64(1), a8), np.where(z, ((a8 & 1) << 13) | (t >> 1), t)


def reference_all():
    """target and tol (uint64, units of 2^-64 turn) for x = 1..MAX_X, exactly as the official verifier."""
    path = os.path.join(WS, "scratch", "ref_all.npz")
    if os.path.exists(path):
        d = np.load(path)
        return d["target"], d["tol"]
    tg = np.empty(MAX_X, dtype=np.uint64)
    tl = np.empty(MAX_X, dtype=np.uint64)
    m64 = (1 << 64) - 1
    for x in range(1, MAX_X + 1):
        r = math.isqrt((1 << 130) // x)
        tg[x - 1] = ((r + 1) >> 1) & m64
        tl[x - 1] = min(r >> 25, 1 << 63)
    np.savez(path, target=tg, tol=tl)
    return tg, tl


def exact_reg(P: Params, a8, t):
    """Exact phase word 2^(26+glo) F in reg units (float64) for segment a8 and t."""
    p = a8 >> 7
    m = 1.0 + (a8 & 127) / 128.0 + np.asarray(t, dtype=np.float64) / 2.0 ** 21
    return 2.0 ** (26 + P.glo) * 2.0 ** (-p / 2.0) / np.sqrt(m)


def segment_sets(P: Params, a8):
    """Constraint sets of a segment: list of (t values, dropped low bits of lo, margin weight)."""
    p, a7 = a8 >> 7, a8 & 127
    base = (1 << 21) + (a7 << 14)
    out = []
    for dd in range(0, P.glo + 1):                 # level j = 10 - dd
        s = 1 + 2 * dd - p
        ts = np.arange(0, 1 << 14, 1 << s, dtype=np.int64)
        if dd == 0 and p == 1:
            ts = ts[base + ts <= MAX_X]
        drop = max(P.virt, 25 + P.glo + (10 - dd) - P.top_pos)
        if len(ts):
            out.append((ts, drop, float(1 << dd)))
    return out


def fit_segment(args):
    """Choose (H, C0, C1, C2) for one segment maximising the absolute margin
    min over constraint sets of weight * (tol - |err|), in reg units. Segments whose margin stays below
    P.want (in units of 1/64 reg) are searched again with doubled and quadrupled spans."""
    P, a8 = args
    best = None
    for mult in (1, 2, 4):
        Q = dataclasses.replace(P, c1_span=P.c1_span * mult, c2_span=P.c2_span * mult, w_span=P.w_span + mult - 1)
        r = _fit_segment(Q, a8)
        if best is None or r["margin"] > best["margin"]:
            best = r
        if best["margin"] >= P.want / 64.0:
            break
    return best


def _fit_segment(P, a8):
    sets = segment_sets(P, a8)
    tall = np.arange(0, 1 << 14, dtype=np.int64)
    tau = 2.0 * tall - (2 ** 14 - 1)
    ex = exact_reg(P, a8, tall)
    tol_all = ex / 2.0 ** 24
    q, sl, _ = np.polyfit(tau, ex, 2)
    c2c = int(round(q * 2.0 ** (27 - P.glo + P.fc2)))                 # c2 = Q * 2^(14+gu-glo) * 2^(13+fc2-gu)
    c1c = int(round(-sl * 2.0 ** (13 + P.gu - (P.glo - 1))))
    c2s = np.arange(c2c - (P.c2_span << P.fc2), c2c + (P.c2_span << P.fc2) + 1, dtype=np.int64)
    if P.c2_even:
        c2s = c2s[(c2s & 1) == 0]
    if C2_FIXED is not None:                 # experiments only (invsqrt_gen_c2exp.py): C2 given, fit C1 and C0
        c2s = np.array([C2_FIXED[a8]], dtype=np.int64)
    c1s = np.arange(c1c - (P.c1_span << P.gu), c1c + (P.c1_span << P.gu) + 1, dtype=np.int64)
    C2g, C1g = np.meshgrid(c2s, c1s, indexing="ij")
    res = []
    for ts, drop, wgt in sets:
        _, R = rows(P, C1g[:, :, None], C2g[:, :, None], ts[None, None, :])
        res.append((R, ex[ts], tol_all[ts], drop, wgt))
    R9, ex9 = res[-1][0], res[-1][1]
    mid = np.rint(0.5 * ((ex9 - R9).max(axis=2) + (ex9 - R9).min(axis=2))).astype(np.int64)
    best = None
    allowed = None
    if P.tie_c1:
        assert P.fc2 == P.gu, "tie_c1 needs the first A row unshifted"
        allowed = ((C1g ^ (C2g >> P.a_shift(P.a_lo))) & 1) == 0
    t0g = None
    if P.tie_t0:
        qmax = sum(rr(C2g, P.a_shift(i)) for i in P.a_rows)
        safe = (qmax <= (1 << (P.win - 1))) & (C1g >= (1 << P.win))
        allowed = safe if allowed is None else (allowed & safe)
        t0g = P.top_table(C1g)[1] & 1
    step = 2 if P.virt else 1
    if P.virt:
        mid = (mid >> 1) << 1
    for dw in range(-P.w_span * step, P.w_span * step + 1, step):
        W = mid + dw
        m = None
        for R, exs, tols, drop, wgt in res:
            v = W[:, :, None] + R
            if drop:
                # floor form: the dropped low bits of lo round the phase down; negated form: Hc half - floor(lo)
                v = -(((-v) >> drop) << drop) if P.neg else (v >> drop) << drop
            mm = wgt * (tols - np.abs(v - exs)).min(axis=2)
            m = mm if m is None else np.minimum(m, mm)
        if allowed is not None:
            m = np.where(allowed, m, -1e9)
        if t0g is not None:
            Wc = -W if P.neg else W          # the low wlo - 1 bits of C0 are those of W (floor form) or of -W
            par = np.zeros_like(W)
            for bit in range(1, P.wlo - 1):
                if (P.tie_mask >> bit) & 1:
                    par ^= (Wc >> bit) & 1
            assert P.tie_mask < (1 << (P.wlo - 1)), "the tie reads C0 below its top bit"
            m = np.where(par == t0g, m, -1e9)
        k = np.unravel_index(np.argmax(m), m.shape)
        if best is None or m[k] > best[0]:
            best = (float(m[k]), int(C2g[k]), int(C1g[k]), int(W[k]))
    marg, c2, c1, W = best
    u, R = rows(P, np.int64(c1), np.int64(c2), tall)
    valid = np.unique(np.concatenate([s[0] for s in sets]))
    assert u[valid].min() >= 0 and u[valid].max() < (1 << P.wu), (a8, int(u[valid].min()), int(u[valid].max()))
    assert 0 <= c2 < (1 << P.wc2) and 0 <= c1 < (1 << P.wu), (a8, c1, c2)
    lo_all = W + R[valid]
    half = 1 << (P.wlo - 1)
    if P.neg:
        Hc = -((-int(lo_all.max())) // half)
        assert Hc * half - int(lo_all.min()) < (1 << P.wlo), (a8, "lo range")
        N = (1 << P.wh) - Hc
        assert 0 <= N < (1 << (P.wh - 1)), (a8, Hc)
        return dict(a8=a8, H=N, C0=(Hc * half - W) % (1 << P.wlo), C1=c1, C2=c2, margin=marg)
    H = int(lo_all.min()) // half
    assert int(lo_all.max()) - H * half < (1 << P.wlo), (a8, "lo range")
    assert 0 <= H < (1 << P.wh), (a8, H)
    return dict(a8=a8, H=H, C0=(W - H * half) % (1 << P.wlo), C1=c1, C2=c2, margin=marg)


def fit_tables(P: Params, workers=3):
    assert not P.zone, "the zone format's tables come from a file (tabs=<name>), not from this fit"
    path = os.path.join(WS, "scratch", f"tables_{P.key()}.json")
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    import multiprocessing as mp
    with mp.get_context("fork").Pool(workers) as pool:
        segs = pool.map(fit_segment, [(P, a8) for a8 in range(256)], chunksize=8)
    tabs = dict(H=[s["H"] for s in segs], C0=[s["C0"] for s in segs], C1=[s["C1"] for s in segs],
                C2=[s["C2"] for s in segs], seg_margin=[s["margin"] for s in segs],
                params=dataclasses.asdict(P))
    with open(path, "w") as fh:
        json.dump(tabs, fh)
    return tabs


def load_tables(name):
    """Tables from a file instead of our own fit: ~/qcc/ws/invsqrt_gen/tabs/<name>.json with the lists H, C0,
    C1, C2 (256 entries each; "H" is N = 2^wh - Hc for neg=1), e.g. a tie set of the table-freedom workstream
    imported by invsqrt_gen_tabs.py. Nothing in such a file is trusted: the builder asserts the ties the
    circuit relies on and the exact check runs every input."""
    path = os.path.join(WS, "tabs", f"{name}.json")
    with open(path) as fh:
        d = json.load(fh)
    tabs = {k: [int(v) for v in d[k]] for k in ("H", "C0", "C1", "C2")}
    assert all(len(v) == 256 for v in tabs.values()), "a table has not 256 entries"
    tabs["sha"] = hashlib.sha256(json.dumps([tabs[k] for k in ("H", "C0", "C1", "C2")]).encode()).hexdigest()[:12]
    return tabs


def evaluate(P: Params, tabs, xs=None):
    """Bit-exact pipeline on xs (default all x). Returns dict of arrays."""
    if xs is None:
        xs = np.arange(1, MAX_X + 1, dtype=np.int64)
    e, p, j, s, a8, t = decompose(xs)
    if P.zone:                                      # zone format: address and digit word of the zone
        a8, t = zone_apply(a8, t)
    H = np.asarray(tabs["H"], dtype=np.int64)[a8]
    C0 = np.asarray(tabs["C0"], dtype=np.int64)[a8]
    C1 = np.asarray(tabs["C1"], dtype=np.int64)[a8]
    C2 = np.asarray(tabs["C2"], dtype=np.int64)[a8]
    u, R = rows(P, C1, C2, t)
    u &= (1 << P.wu) - 1
    lo = ((C0 - R) if P.neg else (C0 + R)) & ((1 << P.wlo) - 1)
    return dict(x=xs, e=e, j=j, a8=a8, t=t, H=H, u=u, lo=lo, c2=C2)


def phase_units(P: Params, ev, angles=None):
    """Delivered phase in units of 2^-64 turn (uint64, mod 1 turn) with nominal power-of-two angles, or with
    `angles` = dict position -> integer units of 2^-64 turn."""
    j = ev["j"]
    tot = np.zeros(len(j), dtype=np.uint64)

    def unit(pos):
        if angles is not None:
            return np.uint64(angles[pos])
        return np.uint64(1 << (63 - pos))           # 2^-pos half turns = 2^(63-pos) units

    for h in range(P.wh - 1 if P.neg else P.wh):    # H bit h at position 7 + j - h
        bit = ((ev["H"] >> h) & 1).astype(bool)
        for jj in range(11):
            sel = bit & (j == jj)
            tot[sel] += unit(7 + jj - h)
    if P.neg:                                       # sign extension of the negated word: ones at positions 0..j-1
        for q in range(10):
            tot[j > q] += unit(q)
    for k in range(P.virt, P.wlo):                  # lo bit k at position 25 + glo + j - k
        bit = ((ev["lo"] >> k) & 1).astype(bool)
        for jj in range(11):
            pos = 25 + P.glo + jj - k
            if pos > P.top_pos:
                continue
            sel = bit & (j == jj)
            tot[sel] += unit(pos)
    if P.neg:                                       # every angle is negative
        tot = (np.uint64(0) - tot).astype(np.uint64)
    return tot


def check(P: Params, tabs, angles=None, n_rot=None):
    target, tol = reference_all()
    ev = evaluate(P, tabs)
    th = phase_units(P, ev, angles)
    err = (th - target).view(np.int64)
    used = np.abs(err)
    ratio = used / tol.astype(np.float64)
    margin = tol.astype(np.int64) - used
    k = int(np.argmax(ratio))
    km = int(np.argmin(margin))
    out = dict(worst=float(ratio[k]), worst_x=int(ev["x"][k]), min_margin=int(margin[km]),
               min_margin_x=int(ev["x"][km]), margin_frac=float(margin[km] / tol[-1]),
               ok=bool((used <= tol).all()), n_fail=int((used > tol).sum()))
    out["by_level"] = {int(e): (float(ratio[ev["e"] == e].max()), float(np.sqrt((ratio[ev["e"] == e] ** 2).mean())))
                       for e in range(22)}
    if n_rot and out["min_margin"] > 0:
        eps = math.pi * out["min_margin"] / 2 ** 64 / n_rot
        out["t_per_rot"] = 0.56 * math.log2(1 / eps) + 5.3
    return out, ev, err


def main():
    import sys
    import time
    kw = {}
    for a in sys.argv[1:]:
        k, v = a.split("=")
        kw[k] = int(v)
    P = Params(**kw)
    t0 = time.time()
    tabs = fit_tables(P)
    t1 = time.time()
    out, ev, err = check(P, tabs, n_rot=44)
    lv = out.pop("by_level")
    print(f"params {kw} key {P.key()} fit {t1 - t0:.1f}s check {time.time() - t1:.1f}s")
    print(json.dumps(out))
    print("level max/rms: " + "  ".join(f"{e}:{lv[e][0]:.3f}/{lv[e][1]:.3f}" for e in range(12, 22)))
    sm = np.array(tabs["seg_margin"])
    print(f"segment margins (reg units): min {sm.min():.3f} at a8={int(sm.argmin())}, median {np.median(sm):.3f}, "
          f"negative {int((sm < 0).sum())}")
    print(f"H {min(tabs['H'])}..{max(tabs['H'])}  C1 {min(tabs['C1'])}..{max(tabs['C1'])}  C2 {min(tabs['C2'])}..{max(tabs['C2'])}")


if __name__ == "__main__":
    main()
