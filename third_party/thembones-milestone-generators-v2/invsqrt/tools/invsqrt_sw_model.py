#!/usr/bin/env python3
"""invsqrt_sw_model.py -- S1 sweep: exact integer model of the root-and-divide invsqrt family (our fork of the
public 46-wire circuit), parameterised, with the circuit's phase rebuilt from the model's word and a fast
max-min-margin angle LP.

Own code. Reads as DATA only: the chain cache of tools/r31_inv-hyb_model.py (k, v, d, h, rc, target, tol for all
3,133,452 inputs), coefficient tables (the public model JSON's D, or a D-per-root .npy), and the description of
the public circuit's phase gates, built once by `describe` with our own simulator (tools/invsqrt_gen_circ via the
re-pricer's Circuit class) from the public circuit file pub.kmx. No public code is imported or run.

Model chain (the public oracle's formulas as read by the R31 seats, here with every knob open):
  coefficient  C(v) = 2 D(v), D from the table (128 rows x 5 coefficient rows, or 64 rows x 6, ...)
  seed         U0 = -1 + sum_{i<3} (2 rc_i - 1) ((C + 2^(11-i)) >> (12-i))
  product      F = sum_{i<10} (2 q_i - 1) term_i(C), q = h + 512; row i < 8 rounds by a mux (public) or plainly
               (bit i of pmask: term = round(C / 2^(9-i)) * 2); U = U0 - F // 2   (15-bit signed)
  divider      27 non-restoring rows j = 9..35; rows j >= 25 add the root truncated by kk_j bits (public
               kk = 0,0,0,2,3,...,9 for j = 25..35); 14-bit remainder (overflow = infeasible)
  endpoint     p = quotient-sign word (37 bits), rotated by k
Units: 2^-64 turn as the verifier (target = round(2^64/sqrt(x)), tol = 2^40/sqrt(x)).

  invsqrt_sw_model.py describe CHAIN.npz MODEL.json PUB.kmx OUT.npz   one-time: phase-gate description
  invsqrt_sw_model.py selftest CHAIN.npz MODEL.json DESC.npz          base reproduces the circuit; LP timing
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

MAX_X = 3133452
MASK37 = (1 << 37) - 1
PUB_KK = (0, 0, 0, 2, 3, 4, 5, 6, 7, 8, 9)       # truncation of divider rows j = 25..35 (public schedule)


# ------------------------------------------------------------------------------------------------ chain + tables
def load_chain(path):
    z = np.load(path)
    keep = ("k", "v", "d", "h", "rc", "target", "tol")
    out = {}
    for key in keep:
        a = z[key]
        out[key] = a if key in ("target", "tol") else a.astype(np.int64)
    return out


def cmd_compact(chain_path, model_json, workdir):
    """write WORKDIR/chain.npz (small dtypes) and WORKDIR/tables/pub.npy (the public D per root, read as data)."""
    z = np.load(chain_path)
    os.makedirs(os.path.join(workdir, "tables"), exist_ok=True)
    np.savez_compressed(os.path.join(workdir, "chain.npz"), k=z["k"].astype(np.int8), v=z["v"].astype(np.int16),
                        d=z["d"].astype(np.int16), h=z["h"].astype(np.int16), rc=z["rc"].astype(np.int8),
                        target=z["target"], tol=z["tol"])
    np.save(os.path.join(workdir, "tables", "pub.npy"), table_pub(model_json))
    print("compact chain and pub table written", flush=True)


def table_pub(model_json):
    m = json.load(open(model_json))
    return np.array(m["D"], dtype=np.int64)


def table_npy(path):
    D = np.load(path).astype(np.int64)
    assert D.shape == (4096,), D.shape
    return D


# ------------------------------------------------------------------------------------------------ arithmetic
def seed_from_C(C, rc, offset=-1):
    """public seed word (three seed rows) as a function of the coefficient and the root class."""
    seed = np.full(len(C), offset, dtype=np.int64)
    for i in range(3):
        seed += (2 * ((rc >> i) & 1) - 1) * ((C + (1 << (11 - i))) >> (12 - i))
    return seed


def seed_range_bad(C, rc, offset=-1):
    """inputs whose seed leaves its register: 3-bit initial value, then 4 / 5 / 6 bits after rows 0 / 1 / 2."""
    seed = np.full(len(C), offset, dtype=np.int64)
    bad = np.full(len(C), offset < -1 or offset > 0)    # only -1 / 0 valid (fixed appended sign bit) [M F3]
    for i in range(3):
        seed += (2 * ((rc >> i) & 1) - 1) * ((C + (1 << (11 - i))) >> (12 - i))
        bad |= (seed < -(1 << (3 + i))) | (seed > (1 << (3 + i)) - 1)
    return bad


def product(C, U0, q, pmask=0):
    """U = U0 - F // 2. Row i (0..9) multiplies residual bit q_i: shift s = 8 - i. Rows i < 8 with bit i of pmask
    round plainly (one bit shorter: term = round_half_up(C / 2^(s+1)) * 2), else the public mux rounding
    (term = round_half_up(C / 2^s))."""
    F = np.zeros_like(q)
    for i in range(10):
        s = 8 - i
        if i < 8 and (pmask >> i) & 1:
            sh = s + 1
            term = ((C + (1 << (sh - 1))) >> sh) << 1
        elif s > 0:
            term = (C + (1 << (s - 1))) >> s
        elif s == 0:
            term = C
        else:
            term = 2 * C + 1
        F += (2 * ((q >> i) & 1) - 1) * term
    return U0 - F // 2


def divide(d, U, kk=PUB_KK):
    """public non-restoring recurrence with per-row truncation kk[j-25] for rows 25..35.
    Returns (word p before rotation, overflow mask of the 14-bit remainder register)."""
    product_bits = U & 32767
    r = 4096 - d
    phase = np.full(len(d), -((1 << 37) - (1 << 29)), dtype=np.int64)
    overflow = np.zeros(len(d), dtype=bool)
    V = d // 2
    for j in range(9, 36):
        s = (r < 0).astype(np.int64)
        phase -= s << (37 - j)
        if j <= 23:
            bit = (product_bits >> (23 - j)) & 1
            if j == 9:
                bit = -bit
            raw = 2 * r + (2 * s - 1) * d + bit
        elif j == 24:
            raw = 2 * r + (2 * s - 1) * d
        else:
            t = kk[j - 25]
            A = ((((V >> (t + 1)) + ((V >> t) & 1)) << (t + 1)) - (V & 1)) if t else V
            u = (r - 1) // 2
            raw = 2 * (2 * u + s + (2 * s - 1) * A) + 1
        nextu = raw // 2
        overflow |= (nextu < -8192) | (nextu >= 8192)
        nextu = ((nextu + 8192) & 16383) - 8192
        r = 2 * nextu + (raw & 1)
    phase -= (r < 0).astype(np.int64) << 1
    p = (phase - 1) & MASK37
    return p, overflow


def rotate(p, k):
    return ((p << k) + (p >> (27 - k))) & MASK37


def word_of(ch, D, pmask=0, kk=PUB_KK, seed_off=-1):
    """the model word (37 bits, rotated) and the overflow mask, for every input."""
    v = ch["v"]
    C = 2 * D[v]
    U0 = seed_from_C(C, ch["rc"], seed_off)
    U = product(C, U0, ch["h"] + 512, pmask)
    rng_bad = (U < -16384) | (U >= 16384) | seed_range_bad(C, ch["rc"], seed_off)
    U = np.clip(U, -16384, 16383)
    p, ovf = divide(ch["d"], U, kk)
    return rotate(p, ch["k"]), ovf | rng_bad


# ------------------------------------------------------------------------------------------------ product in place
def f3_parities(pmask):
    """F3's run parities (tools/invsqrt_fw_patch_prodmask.py): for each multiplexed row i < 8, (C indices, constant)."""
    rows = [i for i in range(8) if not (pmask >> i) & 1]
    runs = []
    for i in rows:
        if runs and runs[-1][1] == i - 1:
            runs[-1][1] = i
        else:
            runs.append([i, i])
    par, acc = {}, set()
    for L, U in runs:
        p = set(acc) ^ {8 - L}
        for i in range(L, U + 1):
            par[i] = (sorted(x for x in p if x != 0), 1)
        acc ^= {v for v in ({8 - L} ^ {8 - (U + 1)}) if v != 0}
    return par


def product_in_place_bad(C, h, U0, pmask, U):
    """the product word as the circuit grows it in place (lifecycle.product_forward with F3's mask patch; public
    parity C[8] for mask 0): row i adds +-(source + carry) into the (6+i)-bit target whose top bit is the previous
    residual bit; returns the inputs whose final 15-bit word differs from U mod 2^15 (a partial sum left its row)."""
    q = h + 512
    bit = lambda x, j: (x >> j) & 1
    par = f3_parities(pmask) if pmask else None
    T = U0 & 63
    for i in range(10):
        n = 6 + i
        qi = bit(q, i)
        if i > 0:
            T = T + (bit(q, i - 1) << (n - 1))
        S = (C >> (9 - i)) if i < 9 else C
        s = 8 - i
        if i < 8 and (pmask >> i) & 1:
            c = bit(C, s)
        elif i == 9:
            c = (qi ^ 1) if pmask else (qi ^ bit(C, 8) ^ 1)
        elif i == 8:
            c = np.zeros_like(C)
        else:
            if par is None:
                P = bit(C, 8)
            else:
                idx, const = par[i]
                P = np.full_like(C, const)
                for k in idx:
                    P ^= bit(C, k)
            c = np.where(qi == P, bit(C, s - 1), bit(C, s))
        T = np.where(qi == 0, T + S + c, T - S - c) & ((1 << n) - 1)
    return T != (U & 32767)


# ------------------------------------------------------------------------------------------------ endpoint frames
def divide_signs_u(d, U, kk):
    """all quotient signs s9..s36 and the 14 remainder wires u after the last row (R31-INV-DIV-reprice 3.2)."""
    product_bits = U & 32767
    r = 4096 - d
    S = {}
    V = d // 2
    for j in range(9, 36):
        s = (r < 0).astype(np.int64)
        S[j] = s
        if j <= 23:
            bit = (product_bits >> (23 - j)) & 1
            if j == 9:
                bit = -bit
            raw = 2 * r + (2 * s - 1) * d + bit
        elif j == 24:
            raw = 2 * r + (2 * s - 1) * d
        else:
            t = kk[j - 25]
            A = ((((V >> (t + 1)) + ((V >> t) & 1)) << (t + 1)) - (V & 1)) if t else V
            u = (r - 1) // 2
            raw = 2 * (2 * u + s + (2 * s - 1) * A) + 1
        nextu = raw // 2
        nextu = ((nextu + 8192) & 16383) - 8192
        r = 2 * nextu + (raw & 1)
    S[36] = (r < 0).astype(np.int64)
    return S, ((r - 1) // 2) & 16383


def affine_fit(t, v, nbits=12):
    """(mask, const) with t == parity(v & mask) ^ const on every input, or None (GF(2) elimination + full check)."""
    rng = np.random.default_rng(3)
    idx = rng.choice(len(t), min(4000, len(t)), replace=False)
    piv = {}
    for i in idx:
        a, b = int(v[i]) | (1 << nbits), int(t[i])
        for col in range(nbits + 1):
            if not (a >> col) & 1:
                continue
            if col in piv:
                pa, pb = piv[col]
                a ^= pa
                b ^= pb
            else:
                piv[col] = (a, b)
                break
        else:
            if b:
                return None
    coef = 0
    for col in sorted(piv, reverse=True):
        a, b = piv[col]
        rest = a & ~(1 << col)
        if b ^ (bin(rest & coef).count("1") & 1):
            coef |= 1 << col
    mask, const = coef & ((1 << nbits) - 1), (coef >> nbits) & 1
    par = np.zeros(len(t), dtype=np.int64)
    for i in range(nbits):
        if (mask >> i) & 1:
            par ^= (v >> i) & 1
    return (mask, const) if np.array_equal(par ^ const, t) else None


try:
    import numba as _nb

    @_nb.njit(cache=True)
    def _kernel(Dv, rc, h, d, k, pmask, kk, seed_off, word, ovf, uw, tail):
        """per input: seed, product (mask), 27 divider rows (truncation kk for rows 25..35), word, overflow, the
        14 remainder wires u and the signs of rows 25..36 (bit j-25 of tail)."""
        n = Dv.shape[0]
        for i in range(n):
            C = 2 * Dv[i]
            seed = seed_off
            sbad = seed < -1 or seed > 0          # only -1 / 0: the rows append a fixed sign bit [M F3 09:37Z]
            for t in range(3):
                sg = 2 * ((rc[i] >> t) & 1) - 1
                seed += sg * ((C + (1 << (11 - t))) >> (12 - t))
                if seed < -(1 << (3 + t)) or seed > (1 << (3 + t)) - 1:
                    sbad = True                   # the seed register after row t holds 4 + t bits
            q = h[i] + 512
            F = 0
            for r in range(10):
                s = 8 - r
                if r < 8 and (pmask >> r) & 1:
                    sh = s + 1
                    term = ((C + (1 << (sh - 1))) >> sh) << 1
                elif s > 0:
                    term = (C + (1 << (s - 1))) >> s
                elif s == 0:
                    term = C
                else:
                    term = 2 * C + 1
                F += (2 * ((q >> r) & 1) - 1) * term
            U = seed - (F // 2)
            bad = U < -16384 or U >= 16384 or sbad
            if U < -16384:
                U = -16384
            if U >= 16384:
                U = 16383
            pb = U & 32767
            di = d[i]
            V = di // 2
            rr = 4096 - di
            phase = -((1 << 37) - (1 << 29))
            tb = 0
            for j in range(9, 36):
                s = 1 if rr < 0 else 0
                if j >= 25:
                    tb |= s << (j - 25)
                phase -= s << (37 - j)
                if j <= 23:
                    bit = (pb >> (23 - j)) & 1
                    if j == 9:
                        bit = -bit
                    raw = 2 * rr + (2 * s - 1) * di + bit
                elif j == 24:
                    raw = 2 * rr + (2 * s - 1) * di
                else:
                    t = kk[j - 25]
                    if t:
                        A = (((V >> (t + 1)) + ((V >> t) & 1)) << (t + 1)) - (V & 1)
                    else:
                        A = V
                    u = (rr - 1) // 2
                    raw = 2 * (2 * u + s + (2 * s - 1) * A) + 1
                nextu = raw // 2
                if nextu < -8192 or nextu >= 8192:
                    bad = True
                nextu = ((nextu + 8192) & 16383) - 8192
                rr = 2 * nextu + (raw & 1)
            s36 = 1 if rr < 0 else 0
            tb |= s36 << 11
            phase -= s36 << 1
            p = (phase - 1) & ((1 << 37) - 1)
            kk_i = k[i]
            word[i] = ((p << kk_i) + (p >> (27 - kk_i))) & ((1 << 37) - 1)
            ovf[i] = bad
            uw[i] = ((rr - 1) // 2) & 16383
            tail[i] = tb
    @_nb.njit(cache=True)
    def _phase_kernel(word, k, fixedsum, units, tr, tb, ti, tg, Bg, target, e):
        """phase of every input from its word (rows that read a word bit, per shift class) + the precomputed
        fixed part; fills the generic rows' bits Bg and the signed error e (2^-64 turn, float)."""
        n = word.shape[0]
        L = tr.shape[1]
        for i in range(n):
            kk = k[i]
            tot = fixedsum[i]
            w = word[i]
            for m in range(L):
                r = tr[kk, m]
                if r < 0:
                    break
                b = ((w >> tb[kk, m]) & 1) ^ ti[kk, m]
                if b:
                    tot += units[r]
                g = tg[kk, m]
                if g >= 0:
                    Bg[g, i] = b
            diff = tot - target[i]
            e[i] = float(np.int64(diff))
    @_nb.njit(cache=True)
    def _margin_kernel(err0, tol, G, dv, marg):
        """marg = tol - |err0 + G^T dv| for every input."""
        K, n = G.shape
        for i in range(n):
            e = err0[i]
            for r in range(K):
                if G[r, i]:
                    e += dv[r]
            marg[i] = tol[i] - abs(e)
    HAVE_NUMBA = True
except Exception:          # local venv without numba: the numpy path below
    HAVE_NUMBA = False


def run_variant_fast(ch, D, pmask, kk, seed_off=-1):
    """numba path of run_variant: same outputs (S holds only rows 25..36, which the frames check reads)."""
    n = len(ch["v"])
    word = np.empty(n, dtype=np.int64)
    ovf = np.empty(n, dtype=np.bool_)
    uw = np.empty(n, dtype=np.int64)
    tail = np.empty(n, dtype=np.int64)
    _kernel(D[ch["v"]].astype(np.int64), ch["rc"], ch["h"], ch["d"], ch["k"], int(pmask),
            np.array(kk, dtype=np.int64), int(seed_off), word, ovf, uw, tail)
    S = {j: (tail >> (j - 25)) & 1 for j in range(25, 37)}
    return dict(word=word, overflow=ovf, S=S, u=uw)


def run_variant(ch, D, pmask, kk, seed_off=-1):
    """ONE divider pass: word (rotated), overflow mask (14-bit register or 15-bit product), signs, u wires."""
    if HAVE_NUMBA:
        return run_variant_fast(ch, D, pmask, kk, seed_off)
    v = ch["v"]
    C = 2 * D[v]
    U = product(C, seed_from_C(C, ch["rc"], seed_off), ch["h"] + 512, pmask)
    rng_bad = (U < -16384) | (U >= 16384) | seed_range_bad(C, ch["rc"], seed_off)
    U = np.clip(U, -16384, 16383)
    product_bits = U & 32767
    d = ch["d"]
    r = 4096 - d
    S = {}
    V = d // 2
    ovf = rng_bad
    phase = np.full(len(d), -((1 << 37) - (1 << 29)), dtype=np.int64)
    for j in range(9, 36):
        s = (r < 0).astype(np.int64)
        S[j] = s
        phase -= s << (37 - j)
        if j <= 23:
            bit = (product_bits >> (23 - j)) & 1
            if j == 9:
                bit = -bit
            raw = 2 * r + (2 * s - 1) * d + bit
        elif j == 24:
            raw = 2 * r + (2 * s - 1) * d
        else:
            t = kk[j - 25]
            A = ((((V >> (t + 1)) + ((V >> t) & 1)) << (t + 1)) - (V & 1)) if t else V
            u = (r - 1) // 2
            raw = 2 * (2 * u + s + (2 * s - 1) * A) + 1
        nextu = raw // 2
        ovf = ovf | (nextu < -8192) | (nextu >= 8192)
        nextu = ((nextu + 8192) & 16383) - 8192
        r = 2 * nextu + (raw & 1)
    S[36] = (r < 0).astype(np.int64)
    phase -= S[36] << 1
    p = (phase - 1) & MASK37
    return dict(word=rotate(p, ch["k"]), overflow=ovf, S=S, u=((r - 1) // 2) & 16383)


def frames_from(S, u, v):
    """endpoint reads from one divider pass (see frames)."""
    ub = [(u >> i) & 1 for i in range(14)]
    v0, v1 = v & 1, (v >> 1) & 1
    ok36 = bool(np.array_equal(S[36], ub[13]))
    oklow = all(np.array_equal(S[j], ub[35 - j] ^ v0) for j in range(27, 36))
    f26 = affine_fit(S[26] ^ ub[9], v)
    f25c = affine_fit(S[25] ^ ub[10], v)
    f25x = None if f25c is not None else affine_fit(S[25] ^ ub[10] ^ ((ub[9] ^ ub[8]) & (v0 ^ v1)), v)
    ok = ok36 and oklow and f26 is not None and (f25c is not None or f25x is not None)
    ftof = None if not ok else (0.0 if f25c is not None else 2.0)
    return dict(ok=ok, ok36=ok36, oklow=oklow, frame_toffolis=ftof, s26=f26, s25cnot=f25c, s25ccx=f25x)


def frames(ch, D, pmask, kk, seed_off=-1):
    """endpoint reads for schedule kk: s36 = u13; s_j = u_(35-j) ^ v0 (j = 27..35); s26 = u9 ^ affine(v);
    s25 = u10 ^ affine(v) (CNOTs only) or u10 ^ affine(v) ^ (u9 ^ u8)(v0 ^ v1) (the public CCX pair).
    Returns dict(ok, frame_toffolis, s26, s25cnot, s25ccx)."""
    v = ch["v"]
    C = 2 * D[v]
    U = np.clip(product(C, seed_from_C(C, ch["rc"], seed_off), ch["h"] + 512, pmask), -16384, 16383)
    S, u = divide_signs_u(ch["d"], U, kk)
    ub = [(u >> i) & 1 for i in range(14)]
    v0, v1 = v & 1, (v >> 1) & 1
    ok36 = bool(np.array_equal(S[36], ub[13]))
    oklow = all(np.array_equal(S[j], ub[35 - j] ^ v0) for j in range(27, 36))
    f26 = affine_fit(S[26] ^ ub[9], v)
    f25c = affine_fit(S[25] ^ ub[10], v)
    f25x = None if f25c is not None else affine_fit(S[25] ^ ub[10] ^ ((ub[9] ^ ub[8]) & (v0 ^ v1)), v)
    ok = ok36 and oklow and f26 is not None and (f25c is not None or f25x is not None)
    ftof = None if not ok else (0.0 if f25c is not None else 2.0)
    return dict(ok=ok, ok36=ok36, oklow=oklow, frame_toffolis=ftof, s26=f26, s25cnot=f25c, s25ccx=f25x)


# ------------------------------------------------------------------------------------------------ identities
def early_signs(d, U, upto=17):
    """quotient signs s_j at the start of rows j = 9..upto-1 (rows < 25 do not depend on the tail schedule)."""
    product_bits = U & 32767
    r = 4096 - d
    S = {}
    for j in range(9, upto):
        s = (r < 0).astype(np.int64)
        S[j] = s
        bit = (product_bits >> (23 - j)) & 1
        if j == 9:
            bit = -bit
        raw = 2 * r + (2 * s - 1) * d + bit
        nextu = raw // 2
        nextu = ((nextu + 8192) & 16383) - 8192
        r = 2 * nextu + (raw & 1)
    return S


def eval_form(form, inputs):
    """straight-line XOR-AND program (data): gate = XOR(left) AND XOR(right); out = const ^ XOR(terms)."""
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


class Identities:
    """The fork's retire-and-rebuild programs (v9 4 ANDs, v10 rule, s15 7 ANDs, F3's s12 identity, s13 / s14 care
    patterns) must stay exact for a changed product word (R31-INV-HYB-reprice `rebuild`). Alignment (o, pol) of the
    circuit's history signs is found on the public table once."""

    def __init__(self, ch, D_pub, form_dir):
        self.v9f = json.load(open(os.path.join(form_dir, "v9_form_e2.json")))
        self.s15f = json.load(open(os.path.join(form_dir, "hret_s15k7a_s15.json")))
        v = ch["v"]
        self.vb = [(v >> i) & 1 for i in range(12)]
        S0 = self.signs(ch, D_pub, 0, -1)
        self.o, self.pol = None, None
        for o in (-1, 0, 1):
            for pol in (0, 1):
                s = {j: S0[j + o] ^ pol for j in range(10, 17) if 9 <= j + o <= 16}
                if all(j in s for j in range(12, 17)):
                    f = eval_form(self.v9f, [self.vb[8], s[12], s[13], s[14], s[15], s[16]])
                    if int((f != self.vb[9]).sum()) == 0 and self.o is None:
                        self.o, self.pol = o, pol
        assert self.o is not None, "no alignment of the history signs"
        self.base = {}
        _, pats = self.check_signs(S0)
        self.base = pats

    def signs(self, ch, D, pmask, seed_off):
        v = ch["v"]
        C = 2 * D[v]
        U = np.clip(product(C, seed_from_C(C, ch["rc"], seed_off), ch["h"] + 512, pmask), -16384, 16383)
        return early_signs(ch["d"], U, 18)

    def check_signs(self, S0):
        vb = self.vb
        s = {j: S0[j + self.o] ^ self.pol for j in range(10, 17) if 9 <= j + self.o <= 17}
        f9 = eval_form(self.v9f, [vb[8], s[12], s[13], s[14], s[15], s[16]])
        g0 = ((vb[8] ^ s[13] ^ s[14] ^ s[15] ^ s[16]) & (vb[8] ^ s[12] ^ s[13] ^ s[14] ^ s[16]))
        g1 = (s[12] ^ s[14]) & (vb[8] ^ s[13] ^ s[15] ^ g0)
        f10 = vb[9] ^ g1 ^ ((s[13] ^ g1) & (vb[9] ^ s[14]))
        f15 = eval_form(self.s15f, [vb[6], vb[7], vb[8], vb[9], vb[10], vb[11], s[16]])
        bad = dict(v9=int((f9 != vb[9]).sum()), v10=int((f10 != vb[10]).sum()),
                   s15=int((f15 != s[15]).sum()), s12=int((s[12] != (vb[11] ^ (vb[10] & (1 - s[13])))).sum()))
        pats = {}
        for tgt, sup in ((13, [8, 9, 10, 11]), (14, [7, 8, 9, 10, 11])):
            key = np.zeros(len(vb[0]), dtype=np.int64)
            for i, b in enumerate(sup):
                key |= vb[b] << i
            key = key | (s[tgt + 1] << len(sup))
            pats[tgt] = set(np.unique(key * 2 + s[tgt]).tolist())
            if self.base:
                bad[f"s{tgt}"] = len(pats[tgt] - self.base[tgt])
        return bad, pats

    def check(self, ch, D, pmask, seed_off):
        bad, _ = self.check_signs(self.signs(ch, D, pmask, seed_off))
        return bad


# ------------------------------------------------------------------------------------------------ circuit phase
class Phase:
    """The public circuit's phase gates, each described per shift class k as a bit of the model word (or a
    constant); rows that are not word bits (x-predicate Cliffords) keep their simulated bits."""

    def __init__(self, desc_npz, ch):
        z = np.load(desc_npz, allow_pickle=False)
        self.units = z["units"].astype(np.uint64)
        self.generic = z["generic"].astype(np.int64)
        self.nrow = len(self.units)
        self.base = np.unpackbits(z["bits"], axis=1)[:, :MAX_X]
        self.desc = json.loads(str(z["desc"]))
        self.k = ch["k"]
        self.classes = [np.flatnonzero(self.k == kk) for kk in range(11)]
        self.tol = ch["tol"].astype(np.float64)
        self.target = ch["target"]
        self.m_pub = float(z["m_pub"])
        # per class: rows that read a word bit (row index, bit, invert) and rows that are a constant in the class
        self.wb, self.cst = [], []
        for kk in range(11):
            rows, bits, inv, crow, cval = [], [], [], [], []
            for r, per in enumerate(self.desc):
                if per is None or str(kk) not in per:
                    continue
                how = per[str(kk)]
                if how[0] == "k":
                    crow.append(r)
                    cval.append(how[1])
                else:
                    rows.append(r)
                    bits.append(how[1])
                    inv.append(how[2])
            self.wb.append((np.array(rows, dtype=np.int64), np.array(bits, dtype=np.int64),
                            np.array(inv, dtype=np.uint8)))
            self.cst.append((np.array(crow, dtype=np.int64), np.array(cval, dtype=np.uint8)))

    def _prep_fast(self):
        """precompute, for the numba phase kernel: per input the phase of every row that does not read the word
        (fixed rows, class constants), the generic-row bits that do not read the word, and per class the
        (row, bit, invert, generic index) of the word-reading rows."""
        n = self.base.shape[1]
        gidx = {int(g): i for i, g in enumerate(self.generic)}
        fixedsum = np.zeros(n, dtype=np.uint64)
        Bg0 = self.base[self.generic].copy()
        for r, per in enumerate(self.desc):
            if per is None and self.units[r]:
                fixedsum += self.base[r].astype(np.uint64) * self.units[r]
        tabs = []
        for kk in range(11):
            sel = self.classes[kk]
            crow, cval = self.cst[kk]
            for r, cv in zip(crow, cval):
                if cv and self.units[r]:
                    fixedsum[sel] += self.units[r]
                if int(r) in gidx:
                    Bg0[gidx[int(r)], sel] = cv
            rows, bits, inv = self.wb[kk]
            tabs.append([(int(r), int(b), int(iv), gidx.get(int(r), -1)) for r, b, iv in zip(rows, bits, inv)])
        L = max(len(t) for t in tabs)
        tr = np.full((11, L), -1, dtype=np.int64)
        tb = np.zeros((11, L), dtype=np.int64)
        ti = np.zeros((11, L), dtype=np.int64)
        tg = np.full((11, L), -1, dtype=np.int64)
        for kk, t in enumerate(tabs):
            for m, (r, b, iv, g) in enumerate(t):
                tr[kk, m], tb[kk, m], ti[kk, m], tg[kk, m] = r, b, iv, g
        self.fast = dict(fixedsum=fixedsum, Bg0=Bg0, tr=tr, tb=tb, ti=ti, tg=tg)

    def evaluate_fast(self, word, lp=True, **kw):
        if not hasattr(self, "fast"):
            self._prep_fast()
        f = self.fast
        Bg = f["Bg0"].copy()
        e = np.empty(len(word), dtype=np.float64)
        _phase_kernel(word, self.k, f["fixedsum"], self.units, f["tr"], f["tb"], f["ti"], f["tg"], Bg,
                      self.target, e)
        marg = self.tol - np.abs(e)
        m_pa, nfail = float(marg.min()), int((marg < 0).sum())
        if not lp:
            return m_pa, nfail, None, None
        m_lp, dv = lp_fast(e, self.tol, Bg, **kw)
        return m_pa, nfail, m_lp, dv

    def rebuild(self, word):
        B = self.base.copy()
        for kk in range(11):
            sel = self.classes[kk]
            if len(sel) == 0:
                continue
            rows, bits, inv = self.wb[kk]
            if len(rows):
                w = word[sel]
                B[np.ix_(rows, sel)] = (((w[None, :] >> bits[:, None]) & 1).astype(np.uint8)) ^ inv[:, None]
            crow, cval = self.cst[kk]
            if len(crow):
                B[np.ix_(crow, sel)] = cval[:, None]
        return B

    def err_of_bits(self, B):
        tot = np.zeros(B.shape[1], dtype=np.uint64)
        for r in range(B.shape[0]):
            if self.units[r]:
                tot += B[r].astype(np.uint64) * self.units[r]
        return ((tot - self.target).view(np.int64)).astype(np.float64)

    def evaluate(self, word, lp=True, drop=(), **kw):
        """(margin with the public circuit's angles, failing inputs, LP margin, offsets). drop: indices into the
        generic rotations to remove (their gates deleted, the LP runs over the others)."""
        if HAVE_NUMBA and not drop and not getattr(self, "slow", False):
            return self.evaluate_fast(word, lp, **kw)
        B = self.rebuild(word)
        keep = [g for i, g in enumerate(self.generic) if i not in set(drop)]
        for i in drop:
            B[self.generic[i]] = 0
        e = self.err_of_bits(B)
        marg = self.tol - np.abs(e)
        m_pa, nfail = float(marg.min()), int((marg < 0).sum())
        if not lp:
            return m_pa, nfail, None, None
        m_lp, dv = lp_fast(e, self.tol, B[keep], **kw)
        return m_pa, nfail, m_lp, dv

    def wrap_features(self):
        """generic rotations that are constant 0 in class k = 0 (they act only for larger shifts): the wrap
        features (1 AND + 1 rotation each in the public endpoint)."""
        out = []
        for i, g in enumerate(self.generic):
            per = self.desc[g]
            if per is not None and per.get("0", ["k", 0]) == ["k", 0]:
                out.append(i)
        return out


def tper(m, n=33):
    return 0.56 * math.log2(n * 2 ** 64 / (math.pi * m)) + 5.3


def lp_fast(err0, tol, G, start=4000, add=4000, rounds=40, verbose=False, stop_below=0.0):
    """max over integer offsets o (one per row of G) of min_x tol - |err0 + G^T o|: cutting planes, HiGHS.
    Small active sets (the optimum has at most K + 1 tight constraints). Returns (all-input margin, offsets).
    The LP over a subset of inputs bounds the true margin from ABOVE: once that bound is <= stop_below the set
    is infeasible (or useless) and the bound is returned at once."""
    from scipy.optimize import linprog
    K = G.shape[0]
    m0 = tol - np.abs(err0)
    act = np.unique(np.argpartition(m0, start)[:start])
    best = (float(m0.min()), np.zeros(K))
    Gf = None
    for rnd in range(rounds):
        Bm = G[:, act].T.astype(np.float64)
        A = np.vstack([np.hstack([Bm, np.ones((len(act), 1))]), np.hstack([-Bm, np.ones((len(act), 1))])])
        b = np.concatenate([tol[act] - err0[act], tol[act] + err0[act]])
        cvec = np.zeros(K + 1)
        cvec[-1] = -1.0
        res = linprog(cvec, A_ub=A, b_ub=b, bounds=[(-2.0 ** 44, 2.0 ** 44)] * K + [(None, None)],
                      method="highs")
        if res.status != 0:
            break
        if res.x[-1] <= stop_below:
            return (float(min(res.x[-1], best[0])), np.round(res.x[:K]))
        dv = np.round(res.x[:K])
        if HAVE_NUMBA:
            marg = np.empty(len(err0), dtype=np.float64)
            _margin_kernel(err0, tol, G, dv.astype(np.float64), marg)
        else:
            e = err0.copy()
            for r in np.flatnonzero(dv):
                e += dv[r] * G[r]
            marg = tol - np.abs(e)
        mmin = float(marg.min())
        if mmin > best[0]:
            best = (mmin, dv.copy())
        if verbose:
            print(f"    lp round {rnd}: {len(act)} rows, LP m {res.x[-1]:.0f}, all-input m {mmin:.0f}", flush=True)
        if mmin >= res.x[-1] - 2 * K:
            break
        cand = np.argpartition(marg, add)[:add]
        cand = cand[marg[cand] < res.x[-1] - 2 * K]
        new = np.setdiff1d(cand, act, assume_unique=False)
        if len(new) == 0:
            break
        act = np.concatenate([act, new])
    return best


# ------------------------------------------------------------------------------------------------ commands
def cmd_describe(chain_path, model_json, kmx, out):
    """simulate the public circuit once (our simulator, via the re-pricer's Circuit class) and save the phase-gate
    description: units per gate, generic gate rows, base bits, per-class word-bit map, public margin."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("rpm", os.path.join(HERE, "r31_inv_hyb_rp_margin.py"))
    RPM = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(RPM)
    z = np.load(chain_path)
    c = {key: z[key] for key in z.files}
    cir = RPM.Circuit(c, model_json, kmx)
    desc = []
    for per in cir.desc:
        if per is None:
            desc.append(None)
        else:
            desc.append({str(kk): [h[0], int(h[1])] + ([int(h[2])] if len(h) > 2 else []) for kk, h in per.items()})
    np.savez_compressed(out, units=cir.units, generic=np.array(cir.generic), bits=np.packbits(cir.bits, axis=1),
                        desc=json.dumps(desc), m_pub=cir.m_pub)
    print(f"saved {out}: rows {len(cir.units)}, generic {len(cir.generic)}, m_pub {cir.m_pub:.0f}", flush=True)


def cmd_selftest(chain_path, model_json, desc_npz):
    t0 = time.time()
    ch = load_chain(chain_path)
    ph = Phase(desc_npz, ch)
    D = table_pub(model_json)
    print(f"loaded in {time.time() - t0:.1f} s; m_pub {ph.m_pub:.0f}", flush=True)
    for name, pm, kk in (("pub", 0, PUB_KK), ("s27", 0, (0, 0, 1, 3, 4, 5, 6, 7, 8, 9, 10)),
                         ("p0to7", 0xFF, PUB_KK)):
        t0 = time.time()
        w, ovf = word_of(ch, D, pm, kk)
        t1 = time.time()
        m_pa, nf, m_lp, dv = ph.evaluate(w, verbose=True)
        print(f"{name}: overflow {int(ovf.sum())}; circuit angles m {m_pa:.0f} fail {nf}; LP m {m_lp:.0f} "
              f"T/rot {tper(m_lp):.4f}; word {t1 - t0:.1f} s, phase+LP {time.time() - t1:.1f} s", flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "describe":
        cmd_describe(*sys.argv[2:6])
    elif cmd == "selftest":
        cmd_selftest(*sys.argv[2:5])
    elif cmd == "compact":
        cmd_compact(*sys.argv[2:5])
