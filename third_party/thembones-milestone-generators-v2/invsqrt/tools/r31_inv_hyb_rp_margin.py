#!/usr/bin/env python3
"""r31_inv_hyb_rp_margin.py -- re-pricer of design INV-HYB (R31 panel): what the by-products of the design really
cost on the verifier's own terms.

Own code. Imports the designer's own model (tools/r31_inv-hyb_model.py: chain, divider recurrence, product rows)
and F1's own angle LP (tools/invsqrt_fk_angles.solve). The public model JSON (tables A, S, D, C, U0 of the public
46-wire circuit, Syed Asad Ali, Apache-2.0) is READ AS DATA. The rotation bits come from our own simulator run on
the public circuit FILE (tools/invsqrt_fk_angles.py extract pub.kmx OUT.npz): data, no public code is run.

Why: the verifier does not price rotations by the worst |err|/tol. It prices them by the ABSOLUTE minimum margin
m = min_x (tol(x) - |err(x)|) (units 2^-64 turn): T per rotation = 0.56 log2(n_rot 2^64 / (pi m)) + 5.3. The
circuit's 33 generic angles carry calibrated offsets (LP-optimal for the public error), so a variant must be read
with those angles and then with the angles re-tuned for the variant.

  rows    CACHE MODEL.json               the circuit's product rows (rounding-carry mux, parity wire) checked
                                         against U0 - F//2 on every input, public and with coarse rows
  lp      CACHE MODEL.json KMX [name ...]  margin with the circuit's angles and after the LP, per variant; names:
                                         base rows0 rows01 rows012 rows0123 rows01234 rows0to7 tail34,
                                         mask:0x3F (product rows of the mask), D64:MASK:TABLE.npy (64-row table)
  pick    CACHE MODEL.json KMX MASK      re-pick (A, S) of failing segments by absolute margin, then LP
  rebuild CACHE MODEL.json FORM_DIR      the fork's retire-and-rebuild programs (v9, v10, s12..s15) per variant
  fit64   CACHE MODEL.json KMX MASK      a 64-row coefficient table fitted for margin, then LP (RP_OUT saves it)
  coef    CACHE MODEL.json KMX ROWMASK hi|lo PRODMASK   coefficient rows rounded plainly, then LP
KMX = the public circuit as kmx text (e.g. ~/qcc/ws/invsqrt_fk/ang/pub.kmx); simulated once by our own simulator.
"""
import importlib.util
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("hyb", os.path.join(HERE, "r31_inv-hyb_model.py"))
H = importlib.util.module_from_spec(spec)
spec.loader.exec_module(H)
sys.path.insert(0, HERE)

MAX_X = 3133452
TOF_PER_COARSE_ROW = 2.0      # rounding-carry mux: 1 CCX forward + 1 CCX in the inverse (del_and is free) [P source]


def tper(m, n):
    return 0.56 * math.log2(n * 2 ** 64 / (math.pi * m)) + 5.3


# ---------------------------------------------------------------------------------------------------------------
def circuit_rows(C, q, coarse=0):
    """Accumulated product correction of the circuit's ten rows, as the generator emits them (lifecycle
    product_forward / _product_carry, parse-only reading): row i adds (1 - 2 q_i) (X_i + c_i), X_i = C >> (9 - i)
    (C for row 9), c_i = rounding carry: a mux (q_i == parity) ? C[s-1] : C[s] for s = 8 - i >= 1, nothing for
    s = 0, q_i ^ parity ^ 1 for s = -1. parity = C[8] in the public circuit (F's parity is 1 ^ C[8]).
    coarse rows (bit i set, i <= 7): X_i + C[s] with a plain carry C[s] (no mux), and the parity wire of the other
    rows moves to the lowest bit the coarse rows no longer telescope over."""
    bit = lambda j: (C >> j) & 1
    lo_s = min([8 - i for i in range(8) if not (coarse >> i) & 1] + [1])
    # F's parity: 1 + sum over non-coarse rows s=1..8 of C[s]^C[s-1]; telescopes when the non-coarse rows are a
    # contiguous run s = 1..smax; compute it generally instead
    par = np.zeros_like(C)
    for i in range(8):
        if not (coarse >> i) & 1:
            s = 8 - i
            par ^= bit(s) ^ bit(s - 1)
    p = par      # F mod 2 = 1 ^ par; the public circuit's parity wire C[8] equals par when no row is coarse
    acc = np.zeros_like(C)
    for i in range(10):
        s = 8 - i
        qi = (q >> i) & 1
        sg = 1 - 2 * qi
        if s >= 1:
            X = C >> (s + 1)
            if (coarse >> i) & 1:
                c = bit(s)
            else:
                c = np.where(qi == p, bit(s - 1), bit(s))
        elif s == 0:
            X, c = C >> 1, 0
        else:
            X, c = C, qi ^ p ^ 1
        acc += sg * (X + c)
    return acc, lo_s


def cmd_rows(c, model_path):
    import json
    m = json.load(open(model_path))
    Dt = np.array(m["D"], dtype=np.int64)
    v, rc = c["v"], c["rc"]
    C = 2 * Dt[v]
    q = c["h"] + 512
    U0 = H.seed_from_C(C, rc)
    # the designer's product_from_C reads mask 1 as "all nine rows"; row 0 alone is 0x201 (bit 9 has no effect)
    for mask in (0x0, 0x201, 0x3, 0x7, 0xF, 0x1F, 0xFF):
        U = H.product_from_C(c, C, U0, mask)
        acc, _ = circuit_rows(C, q, mask & 0xFF)
        d = (U - U0) - acc
        print(f"coarse mask {mask:#05x}: circuit rows minus model (U - U0): min {int(d.min())} max {int(d.max())} "
              f"-> {'EXACT' if d.min() == d.max() == 0 else 'constant offset' if d.min() == d.max() else 'MISMATCH'}; "
              f"rows with a mux removed {bin(mask & 0xFF).count('1')}: {TOF_PER_COARSE_ROW * bin(mask & 0xFF).count('1'):.0f} "
              f"Toffolis off (forward + inverse)", flush=True)


# ---------------------------------------------------------------------------------------------------------------
class Circuit:
    """The public circuit simulated once by our own simulator (tools/invsqrt_gen_circ.Sim, all inputs, HMR outcomes
    0): every Z_POW gate (generic, T, Clifford) and the half-turn sign part (Z, CZ, NEG) are described, per shift
    class k, as a bit of the model's word or a constant. A variant's phase is then rebuilt gate by gate from the
    variant's word with the circuit's own angles; gates that are not word bits (x-predicates) stay as they are.
    Self-check: the base rebuilt from the model's word equals the circuit's phase on every input."""

    def __init__(self, c, model_path, kmx):
        import json
        import invsqrt_lib as L
        import invsqrt_gen_circ as CC
        self.c = c
        self.m = json.load(open(model_path))
        t0 = time.time()
        header, ops = L.parse_kmx(kmx)
        ops = [op for op in ops if op[0] != "INFO_TAG"]
        nq = 1 + max(q for op in ops for q in op[1])
        sim = CC.Sim(c["x"].astype(np.uint64), nq, hmr="zero")
        sim.run(ops)
        self.angles = [a for a, _ in sim.rot]
        rows = [CC.unpack(b, MAX_X) if b else np.zeros(MAX_X, dtype=np.uint8) for _, b in sim.rot]
        rows.append(CC.unpack(sim.phase, MAX_X) if sim.phase else np.zeros(MAX_X, dtype=np.uint8))
        self.bits = np.array(rows, dtype=np.uint8)
        self.units = np.array([CC.half_turns_to_units(a) for a in self.angles] + [1 << 63], dtype=np.uint64)
        self.generic = [r for r, a in enumerate(self.angles) if (a * 4) % 1 != 0]
        self.K = len(self.generic)
        self.tol = c["tol"].astype(np.float64)
        th = self.phase_of(self.bits)
        self.err_c = ((th - c["target"]).view(np.int64)).astype(np.float64)
        print(f"circuit simulated in {time.time() - t0:.0f} s: Z_POW {len(self.angles)} (generic {self.K}) + sign part",
              flush=True)
        self.k = c["k"]
        self.classes = [np.flatnonzero(self.k == kk) for kk in range(11)]
        Dt = np.array(self.m["D"], dtype=np.int64)
        self.word0 = self.word_of(2 * Dt)
        self.err_m0 = self.err_of(self.word0)
        t0 = time.time()
        self.desc = self.describe(self.word0)
        self.fixed = [r for r, p in enumerate(self.desc) if p is None]
        print(f"phase rows described as word bits per class: {len(self.desc) - len(self.fixed)} of {len(self.desc)}; "
              f"kept fixed (not word bits): {[str(self.angles[r]) if r < len(self.angles) else 'sign' for r in self.fixed]} "
              f"({time.time() - t0:.0f} s)", flush=True)
        assert not [r for r in self.fixed if r in self.generic], "a generic rotation is not a word bit"
        chk = self.err_of_bits(self.rebuild(self.word0))
        print(f"self-check: base rebuilt from the model's word equals the circuit on every input: "
              f"{bool(np.array_equal(chk, self.err_c))}", flush=True)
        assert np.array_equal(chk, self.err_c)
        self.m_pub = float((self.tol - np.abs(self.err_c)).min())
        self.m_nom = float((self.tol - np.abs(self.err_m0)).min())
        print(f"public circuit min margin {self.m_pub:.0f} (T/rot {tper(self.m_pub, self.K):.4f}); the model's nominal "
              f"angles {self.m_nom:.0f}; worst |err|/tol circuit {(np.abs(self.err_c) / self.tol).max():.5f}, nominal "
              f"{(np.abs(self.err_m0) / self.tol).max():.5f}", flush=True)

    def phase_of(self, B):
        tot = np.zeros(B.shape[1], dtype=np.uint64)
        for r in range(B.shape[0]):
            tot += B[r].astype(np.uint64) * self.units[r]
        return tot

    def err_of_bits(self, B):
        return ((self.phase_of(B) - self.c["target"]).view(np.int64)).astype(np.float64)

    def rebuild(self, word):
        B = self.bits.copy()
        for r, per in enumerate(self.desc):
            if per is None:
                continue
            for kk, how in per.items():
                sel = self.classes[kk]
                if how[0] == "k":
                    B[r, sel] = how[1]
                else:
                    B[r, sel] = (((word[sel] >> how[1]) & 1).astype(np.uint8)) ^ how[2]
        return B

    def word_of(self, C, coarse=0, U0=None, tail_from=0):
        c = self.c
        if U0 is None:
            U0 = H.seed_from_C(C[c["v"]] if len(C) == 4096 else C, c["rc"])
        Cv = C[c["v"]] if len(C) == 4096 else C
        U = np.clip(H.product_from_C(c, Cv, U0, coarse), -16384, 16383)
        return H.public_word(c["d"], U, c["k"], 0, tail_from)

    def err_of(self, word):
        return ((word.astype(np.uint64) << np.uint64(27)) - self.c["target"]).view(np.int64).astype(np.float64)

    def describe(self, word):
        out = []
        for r in range(self.bits.shape[0]):
            per = {}
            for kk, sel in enumerate(self.classes):
                if len(sel) == 0:
                    continue
                b = self.bits[r][sel]
                if not b.any():
                    per[kk] = ("k", 0)
                    continue
                if b.all():
                    per[kk] = ("k", 1)
                    continue
                wk = word[sel]
                hit = None
                for t in range(37):
                    wt = ((wk >> t) & 1).astype(np.uint8)
                    if np.array_equal(b, wt):
                        hit = ("w", t, 0)
                        break
                    if np.array_equal(b, wt ^ 1):
                        hit = ("w", t, 1)
                        break
                if hit is None:
                    per = None
                    break
                per[kk] = hit
            out.append(per)
        return out

    def variant(self, word, lp=True):
        """(margin with the circuit's angles, failing inputs, margin after the LP over the generic angles, offsets)."""
        B = self.rebuild(word)
        e = self.err_of_bits(B)
        marg = self.tol - np.abs(e)
        m_pa, nfail = float(marg.min()), int((marg < 0).sum())
        if not lp:
            return m_pa, nfail, None, None
        import invsqrt_fk_angles as FA
        m_lp, dv = FA.solve(e, self.tol, B[self.generic])
        return m_pa, nfail, float(m_lp), dv

    def charge(self, m):
        return self.K * tper(m, self.K) / 4 if m > 0 else float("inf")


def variants_table(cir):
    Dt = np.array(cir.m["D"], dtype=np.int64)
    C = 2 * Dt
    return {
        "base": (C, 0, 0, 0),
        "rows01": (C, 0x3, 0, 2),
        "rows012": (C, 0x7, 0, 3),
        "rows0123": (C, 0xF, 0, 4),
        "rows01234": (C, 0x1F, 0, 5),
        "rows0to7": (C, 0xFF, 0, 8),
        "rows0": (C, 0x201, 0, 1),
        "tail34": (C, 0, 34, 0),
    }


def report(cir, name, word, coarse_rows, lp=True):
    t0 = time.time()
    m_pa, nfail, m_lp, _ = cir.variant(word, lp)
    base_charge = cir.charge(cir.m_pub)
    dtof = -TOF_PER_COARSE_ROW * coarse_rows
    if m_lp is not None and m_lp > 0:
        dch = cir.charge(m_lp) - base_charge
        net = dtof + dch
        tail = (f"after LP {m_lp:14.0f} ({m_lp / cir.m_pub:.3f}x): T/rot {tper(m_lp, cir.K):.4f}, rotation charge "
                f"{dch:+6.2f}; Toffolis {dtof:+5.1f}; net {net:+6.2f} Toffoli-eq = {46 * net:+7.1f} points at Q 46")
    else:
        tail = f"after LP {m_lp if m_lp is not None else float('nan'):14.0f}: INFEASIBLE (no angles pass every input)"
    print(f"{name:12s} circuit angles: margin {m_pa:14.0f}, failing {nfail:5d}; {tail} ({time.time() - t0:.0f} s)",
          flush=True)
    return m_pa, nfail, m_lp


def cmd_lp(c, model_path, kmx, names):
    cir = Circuit(c, model_path, kmx)
    table = variants_table(cir)
    for name in names or list(table):
        if name.startswith("mask:"):     # mask:0x3F = product rows 0..5 one bit shorter, public table
            mask = int(name.split(":", 1)[1], 0)
            C, tail, nrows = table["base"][0], 0, bin(mask & 0xFF).count("1")
            mask = mask if mask != 1 else 0x201
        elif name.startswith("D64:"):     # D64:MASK:PATH = a 64-row table from fit64, with product rows of MASK
            _, mk, path = name.split(":", 2)
            mask, tail = int(mk, 0), 0
            C, nrows = 2 * np.load(path), bin(int(mk, 0) & 0xFF).count("1")
        else:
            C, mask, tail, nrows = table[name]
        word = cir.word_of(C, mask, None, tail)
        report(cir, name, word, nrows)


def cmd_pick(c, model_path, kmx, mask, span=4):
    """Re-pick (A, S) per failing segment in the public format (A +- span, stored slope +- 2 span), choosing the
    candidate with the largest ABSOLUTE margin on that segment's inputs under the circuit's angles; then the LP."""
    cir = Circuit(c, model_path, kmx)
    A0 = np.array(cir.m["A"], dtype=np.int64)
    S0 = np.array(cir.m["S"], dtype=np.int64)
    Dt = np.array(cir.m["D"], dtype=np.int64)
    v, rc = c["v"], c["rc"]
    a = v >> 5
    nrows = bin(mask & 0xFF).count("1")
    word = cir.word_of(2 * Dt, mask)
    e = cir.err_of_bits(cir.rebuild(word))
    marg = cir.tol - np.abs(e)
    seg_m = np.full(128, np.inf)
    np.minimum.at(seg_m, a, marg)
    order = np.argsort(seg_m)
    print("tightest segments under the circuit's angles (segment: margin): " +
          ", ".join(f"{int(s)}: {seg_m[s]:.0f}" for s in order[:8]), flush=True)
    bad = sorted(set(a[marg < 0].tolist()))
    print(f"mask {mask:#x}: failing inputs {int((marg < 0).sum())} in segments {bad}", flush=True)
    newD = Dt.copy()
    changed = []
    for sgm in bad:
        idx = np.flatnonzero(a == sgm)
        w = v[idx] & 31
        sub = {key: c[key][idx] for key in ("d", "k", "h", "target", "tol")}
        best = None
        for dH in range(-span, span + 1):
            S1 = int(S0[sgm]) + 2 * dH
            if not 0 <= S1 < 64:
                continue
            for dA in range(-span, span + 1):
                A1 = int(A0[sgm]) + dA
                Cs = 2 * H.D_from_AS(A1, S1, w)
                Us = np.clip(H.product_from_C(sub, Cs, H.seed_from_C(Cs, rc[idx]), mask), -16384, 16383)
                ws = H.public_word(sub["d"], Us, sub["k"])
                # the circuit's own angles on this segment's inputs, every phase row rebuilt from the word
                bs = cir.bits[:, idx].copy()
                kk_arr = sub["k"]
                for r, per in enumerate(cir.desc):
                    if per is None:
                        continue
                    for kk, how in per.items():
                        sel = kk_arr == kk
                        if how[0] == "k":
                            bs[r, sel] = how[1]
                        else:
                            bs[r, sel] = (((ws[sel] >> how[1]) & 1).astype(np.uint8)) ^ how[2]
                es = ((cir.phase_of(bs) - sub["target"]).view(np.int64)).astype(np.float64)
                ms = float((cir.tol[idx] - np.abs(es)).min())
                if best is None or ms > best[0]:
                    best = (ms, A1, S1)
        print(f"  segment {sgm}: A {int(A0[sgm])} -> {best[1]}, S {int(S0[sgm])} -> {best[2]}: segment margin "
              f"{best[0]:.0f} under the circuit's angles", flush=True)
        newD[sgm << 5:(sgm + 1) << 5] = H.D_from_AS(best[1], best[2], np.arange(32))
        changed.append((sgm, int(A0[sgm]), best[1], int(S0[sgm]), best[2]))
    word2 = cir.word_of(2 * newD, mask)
    report(cir, f"m{mask:#x}+pick", word2, nrows)
    # which table words change: loader outputs are A (13 bits) + stored slope bits H0..H2; H3, H4 come from plans
    for sgm, a0, a1, s0, s1 in changed:
        print(f"  table change on row {sgm}: A {a0:013b} -> {a1:013b} ({bin(a0 ^ a1).count('1')} bits), stored "
              f"slope {s0:06b} -> {s1:06b} ({bin(s0 ^ s1).count('1')} bits)", flush=True)
    out_dir = os.environ.get("RP_OUT")
    if out_dir:
        np.save(os.path.join(out_dir, f"D_pick_{mask:#x}.npy"), newD)


# ---------------------------------------------------------------------------------------------------------------
def rows64(c0, S, lb=6):
    """D over one 64-root row in the circuit's coefficient encoding (F1's reading of the shipped rows, extended to
    six rows): D(j) = c0 - sum_i j_i round_half_up(S / 2^(lb-1-i))."""
    j = np.arange(1 << lb)
    out = np.full(1 << lb, c0, dtype=np.int64)
    for i in range(lb):
        den = 1 << (lb - 1 - i)
        out -= ((j >> i) & 1) * ((2 * S + den) // (2 * den))
    return out


def segment_margin(cir, idx, Cv_sub, mask):
    """Min absolute margin of the inputs idx under the circuit's own angles, coefficient values Cv_sub (= 2 D)."""
    c = cir.c
    sub = {key: c[key][idx] for key in ("d", "k", "h", "target", "tol")}
    Us = np.clip(H.product_from_C(sub, Cv_sub, H.seed_from_C(Cv_sub, c["rc"][idx]), mask), -16384, 16383)
    ws = H.public_word(sub["d"], Us, sub["k"])
    bs = cir.bits[:, idx].copy()
    for r, per in enumerate(cir.desc):
        if per is None:
            continue
        for kk, how in per.items():
            sel = sub["k"] == kk
            if how[0] == "k":
                bs[r, sel] = how[1]
            else:
                bs[r, sel] = (((ws[sel] >> how[1]) & 1).astype(np.uint8)) ^ how[2]
    es = ((cir.phase_of(bs) - sub["target"]).view(np.int64)).astype(np.float64)
    return float((cir.tol[idx] - np.abs(es)).min())


def cmd_fit64(c, model_path, kmx, mask, span=4):
    """A 64-row coefficient table (six coefficient rows, slope S of up to 7 bits, intercept c0) fitted per row
    for the largest ABSOLUTE margin under the circuit's angles (rows are independent: an input reads only its own
    root's coefficient), then the angle LP; with the product rows of MASK one bit shorter."""
    cir = Circuit(c, model_path, kmx)
    Dt = np.array(cir.m["D"], dtype=np.int64)
    v = c["v"]
    a6 = v >> 6
    order = np.argsort(a6, kind="stable")
    bounds = np.searchsorted(a6[order], np.arange(65))
    newD = Dt.copy()
    meta = []
    t0 = time.time()
    for a in range(64):
        idx = order[bounds[a]:bounds[a + 1]]
        roots = np.arange(a * 64, a * 64 + 64)
        jj = roots - roots[0]
        slope = -np.polyfit(jj, Dt[roots].astype(np.float64), 1)[0]
        S0 = int(round(slope * 32))
        best = None
        for S in range(S0 - span, S0 + span + 1):
            base = rows64(0, S)
            res = Dt[roots] - base
            cc0 = int(round((res.max() + res.min()) / 2))
            for cc in range(cc0 - span, cc0 + span + 1):
                Drow = cc + base
                if len(idx) == 0:
                    m = 0.0
                else:
                    Cv = 2 * Drow[v[idx] - a * 64]
                    m = segment_margin(cir, idx, Cv, mask)
                if best is None or m > best[0]:
                    best = (m, S, cc)
        m, S, cc = best
        newD[a * 64:(a + 1) * 64] = cc + rows64(0, S)
        meta.append((a, S, cc, m))
    Smax = max(t[1] for t in meta)
    cmax = max(t[2] for t in meta)
    worst = sorted(meta, key=lambda t: t[3])[:4]
    print(f"64-row fit (mask {mask:#x}) in {time.time() - t0:.0f} s: slope up to {Smax} ({Smax.bit_length()} bits), "
          f"intercept up to {cmax} ({cmax.bit_length()} bits); tightest rows {worst}", flush=True)
    word = cir.word_of(2 * newD, mask)
    report(cir, f"D64 m{mask:#x}", word, bin(mask & 0xFF).count("1"))
    if mask:
        report(cir, "D64 alone", cir.word_of(2 * newD, 0), 0)
    out_dir = os.environ.get("RP_OUT")
    if out_dir:
        np.save(os.path.join(out_dir, f"D64_m{mask:#x}.npy"), newD)


# ---------------------------------------------------------------------------------------------------------------
def D_plain(A, S, w, rowmask, opt):
    """The public coefficient rows (designer's D_from_AS, as the circuit's coefficient_forward) with the rows of
    rowmask rounded through ONE carry wire for both signs instead of the 1-CCX multiplexer of _coefficient_carry:
    opt 'hi' uses S[k] for both signs (dn := up), opt 'lo' uses S[k-1] (0 for k = 0) for both (up := dn)."""
    D = np.array(A, dtype=np.int64) + 0 * w
    for i in range(5):
        kk = 4 - i
        bit = (w >> i) & 1
        up = (S >> (kk + 1)) + ((S >> kk) & 1)
        dn = (S >> (kk + 1)) + (((S >> (kk - 1)) & 1) if kk else 0)
        if (rowmask >> i) & 1:
            if opt == "hi":
                dn = up
            else:
                up = dn
        D = D + np.where(bit == 0, up, -dn)
    return D


def cmd_coef(c, model_path, kmx, rowmask, opt, prodmask):
    cir = Circuit(c, model_path, kmx)
    A = np.array(cir.m["A"], dtype=np.int64)
    S = np.array(cir.m["S"], dtype=np.int64)
    Dt = np.array(cir.m["D"], dtype=np.int64)
    vv = np.arange(4096)
    D0 = D_plain(A[vv >> 5], S[vv >> 5], vv & 31, 0, opt)
    assert np.array_equal(D0, Dt), "coefficient rows as read"
    D1 = D_plain(A[vv >> 5], S[vv >> 5], vv & 31, rowmask, opt)
    dd = D1 - Dt
    hi_changed = int(((D1 >> 9) != (Dt >> 9)).sum())
    print(f"coefficient rows {rowmask:#x} plain ({opt}): D changes on {int((dd != 0).sum())} of 4096 roots "
          f"(range {int(dd.min())}..{int(dd.max())}); roots whose bits D9..D12 change: {hi_changed}; "
          f"CCX off: {3 * bin(rowmask).count('1')} (COEF1, COEF2, A_REST)", flush=True)
    ncoef = 3 * bin(rowmask).count("1")
    for pm in sorted({0, prodmask}):
        word = cir.word_of(2 * D1, pm)
        nrows = bin(pm & 0xFF).count("1") + ncoef / TOF_PER_COARSE_ROW
        report(cir, f"coef{rowmask:#x}{opt}+p{pm:#x}", word, nrows)


# ---------------------------------------------------------------------------------------------------------------
def divider_signs(d, U, k, tail_from=0):
    """Quotient signs of the divider, row by row (the designer's recurrence, r31_inv-hyb_model.public_word):
    S[j] = (r < 0) at the start of row j, j = 9..35, and S[36] = the final sign."""
    product = U & 32767
    r = 4096 - d
    S = {}
    for j in range(9, 36):
        sign = (r < 0).astype(np.int64)
        S[j] = sign
        if j <= 23:
            bit = (product >> (23 - j)) & 1
            if j == 9:
                bit = -bit
            raw = 2 * r + (2 * sign - 1) * d + bit
        elif j == 24:
            raw = 2 * r + (2 * sign - 1) * d
        else:
            if tail_from:
                kk = (j - 26 if j >= 28 else 0) + (1 if j >= tail_from else 0)
                if j < 28 and kk:
                    kk = 2
            else:
                kk = j - 26 if j >= 28 else 0
            V = d // 2
            A = ((((V >> (kk + 1)) + ((V >> kk) & 1)) << (kk + 1)) - (V & 1)) if kk else V
            u = (r - 1) // 2
            raw = 2 * (2 * u + sign + (2 * sign - 1) * A) + 1
        nextu = raw // 2
        r = 2 * nextu + (raw & 1)
    S[36] = (r < 0).astype(np.int64)
    return S


def eval_form(form, inputs):
    """A straight-line XOR-AND program given as data (inputs, gates with left / right index lists, output terms and
    constant): gate i = XOR(vals[left]) AND XOR(vals[right]); output = constant ^ XOR(vals[terms])."""
    vals = list(inputs)
    for g in form["gates"]:
        lft = np.zeros_like(vals[0])
        for j in g["left"]:
            lft ^= vals[j]
        rgt = np.zeros_like(vals[0])
        for j in g["right"]:
            rgt ^= vals[j]
        vals.append(lft & rgt)
    out = np.full_like(vals[0], form["output"]["constant"])
    for j in form["output"]["terms"]:
        out ^= vals[j]
    return out


def cmd_rebuild(c, model_path, form_dir, names):
    """The fork retires values across its width peaks and rebuilds them from exact AND programs that are correct
    on the CARE patterns of the public circuit (every x in the domain). A different product changes quotient
    signs on some inputs. For each variant: does each shipped program (or identity) still rebuild its value on
    every input?  v9 <- (v8, s12..s16) [v9_form_e2, 4 ANDs]; v10 <- v9 ^ g1 ^ (s13^g1)(v9^s14); s15 <- (v6..v11, s16)
    [hret_s15k7a, 7 ANDs]; s12 = v11 ^ (v10 & ~s13) [F3's identity]; s13 <- (v8..v11, s14) and s14 <- (v7..v11,
    s15) [patch formulas: care-pattern check]."""
    import json
    m = json.load(open(model_path))
    v9f = json.load(open(os.path.join(form_dir, "v9_form_e2.json")))
    s15f = json.load(open(os.path.join(form_dir, "hret_s15k7a_s15.json")))
    v = c["v"]
    vb = [(v >> i) & 1 for i in range(12)]
    Dt = np.array(m["D"], dtype=np.int64)
    C = 2 * Dt[v]
    U0 = H.seed_from_C(C, c["rc"])

    def signs_for(mask, Dtab=None, tail_from=0):
        Cv = C if Dtab is None else 2 * Dtab[v]
        U = np.clip(H.product_from_C(c, Cv, H.seed_from_C(Cv, c["rc"]), mask), -16384, 16383)
        return divider_signs(c["d"], U, c["k"], tail_from)

    S0 = signs_for(0)
    # alignment and polarity of the circuit's history signs against the model's row signs
    found = None
    for o in (-1, 0, 1):
        for pol in (0, 1):
            s = {j: S0[j + o] ^ pol for j in range(10, 35)}
            f = eval_form(v9f, [vb[8], s[12], s[13], s[14], s[15], s[16]])
            ok9 = int((f != vb[8 + 1]).sum())
            okid = int((s[12] != (vb[11] ^ (vb[10] & (1 - s[13])))).sum())
            print(f"alignment o={o:+d} polarity {pol}: v9 program wrong on {ok9} inputs; s12 identity wrong on {okid}",
                  flush=True)
            if ok9 == 0 and found is None:
                found = (o, pol)
    if found is None:
        print("no alignment makes the shipped v9 program exact on the base: inputs are in a frame; stop", flush=True)
        return
    o, pol = found
    print(f"using o={o:+d}, polarity {pol}", flush=True)

    def check(name, S):
        s = {j: S[j + o] ^ pol for j in range(10, 35)}
        f9 = eval_form(v9f, [vb[8], s[12], s[13], s[14], s[15], s[16]])
        bad9 = int((f9 != vb[9]).sum())
        # v10 from g1 of the v9 program
        vals = [vb[8], s[12], s[13], s[14], s[15], s[16]]
        g0 = ((vb[8] ^ s[13] ^ s[14] ^ s[15] ^ s[16]) & (vb[8] ^ s[12] ^ s[13] ^ s[14] ^ s[16]))
        g1 = (s[12] ^ s[14]) & (vb[8] ^ s[13] ^ s[15] ^ g0)
        f10 = vb[9] ^ g1 ^ ((s[13] ^ g1) & (vb[9] ^ s[14]))
        bad10 = int((f10 != vb[10]).sum())
        f15 = eval_form(s15f, [vb[6], vb[7], vb[8], vb[9], vb[10], vb[11], s[16]])
        bad15 = int((f15 != s[15]).sum())
        bad12 = int((s[12] != (vb[11] ^ (vb[10] & (1 - s[13])))).sum())
        out = []
        for tgt, sup in ((13, [8, 9, 10, 11]), (14, [7, 8, 9, 10, 11])):
            key0 = np.zeros(len(v), dtype=np.int64)
            for i, b in enumerate(sup):
                key0 |= vb[b] << i
            out.append((tgt, sup, key0))
        res = []
        for tgt, sup, key in out:
            nxt = s[tgt + 1]
            kk = key | (nxt << len(sup))
            res.append((tgt, kk, s[tgt]))
        return bad9, bad10, bad15, bad12, res, f9

    b9, b10, b15, b12, res0, _ = check("base", S0)
    base_sets = {}
    for tgt, kk, val in res0:
        pat = np.unique(kk * 2 + val)
        conflict = len(np.unique(kk)) != len(pat)
        base_sets[tgt] = set(pat.tolist())
        print(f"base: s{tgt} a function of its support on the care set: {not conflict} ({len(pat)} patterns)", flush=True)
    print(f"{'base':10s}: v9 program wrong on {b9}, v10 rule wrong on {b10}, s15 program wrong on {b15}, s12 identity "
          f"wrong on {b12}", flush=True)
    variants = {"rows0": (0x201, None, 0), "rows01": (0x3, None, 0), "rows012": (0x7, None, 0),
                "rows0123": (0xF, None, 0), "rows01234": (0x1F, None, 0), "rows0to7": (0xFF, None, 0),
                "tail34": (0, None, 34)}
    for name in names or list(variants):
        if name.startswith("pick"):       # pick:MASK:PATH_TO_D.npy (a re-picked table from the pick command)
            _, mk, path = name.split(":", 2)
            mask, Dtab, tail = int(mk, 0), np.load(path), 0
        else:
            mask, Dtab, tail = variants[name]
        S = signs_for(mask, Dtab, tail)
        nchg = int(np.any(np.array([S[j] != S0[j] for j in range(12, 18)]), axis=0).sum())
        b9, b10, b15, b12, res, _ = check(name, S)
        extra = []
        for tgt, kk, val in res:
            pat = set(np.unique(kk * 2 + val).tolist())
            new = pat - base_sets[tgt]
            extra.append(f"s{tgt}: {len(new)} care patterns outside the public set")
        print(f"{name:10s}: inputs whose signs s12..s17 change {nchg}; v9 program wrong on {b9}, v10 rule wrong on "
              f"{b10}, s15 program wrong on {b15}, s12 identity wrong on {b12}; " + "; ".join(extra), flush=True)


def main():
    cache, cmd = sys.argv[1], sys.argv[2]
    c = H.load_chain(cache)
    if cmd == "rows":
        cmd_rows(c, sys.argv[3])
    elif cmd == "rebuild":
        cmd_rebuild(c, sys.argv[3], sys.argv[4], sys.argv[5:])
    elif cmd == "fit64":
        cmd_fit64(c, sys.argv[3], sys.argv[4], int(sys.argv[5], 0))
    elif cmd == "coef":
        cmd_coef(c, sys.argv[3], sys.argv[4], int(sys.argv[5], 0), sys.argv[6], int(sys.argv[7], 0))
    elif cmd == "lp":
        cmd_lp(c, sys.argv[3], sys.argv[4], sys.argv[5:])
    elif cmd == "pick":
        cmd_pick(c, sys.argv[3], sys.argv[4], int(sys.argv[5], 0))
    else:
        raise SystemExit("unknown command")


if __name__ == "__main__":
    main()
