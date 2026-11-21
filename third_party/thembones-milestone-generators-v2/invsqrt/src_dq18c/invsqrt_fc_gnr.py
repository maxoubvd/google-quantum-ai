#!/usr/bin/env python3
"""invsqrt_fc_gnr.py -- FC: the grouped non-restoring HO stage (GNR) as the circuit computes it (own code).

Arithmetic (shared by the model and the gate builder):
  digits sigma_j in {0,1} (s_j = 1 - 2 sigma_j), j = 18 .. 0, on U's wires; held carries e_g in {0,1} at
  weight 2^(lo_g + 1) for every group g = (hi_g .. lo_g); quotient Phi = sum s_j 2^j + sum e_g 2^(lo_g+1).
  h(z) = 3 z^2 / 2^28 - z^3 / 2^53, z = Phi - B (B = biases added so far / 2^fq).
  The correction after group g is Hq(point_g) - Hq(point_{g-1}), Hq = h expanded in the circuit's basis
  {1, sigma_a, XOR(sigma set), AND(carry set), AND(carry set, sigma or XOR)} with every basis coefficient
  rounded to 2^-fq (same rounding at every point, so the sum telescopes exactly).
  remainder r (units 2^-fq):  r = y 2^fq;  per group: free NR steps (no gates beyond X); r += correction;
  r += bias c_g; e_g = [r >= 2^lo_g 2^fq]; r -= e_g 2^(lo_g+1) 2^fq.  After the last group: free steps to
  digit 0, then the final correction (deferred h-change of the last carry).  U = Phi + floor((r - sum c)/2^fq).
Commands:
  invsqrt_fc_gnr.py WD model '{"fq":3,"groups":[[18,14],[13,9],[8,6]],"keep":-1}'
     keep: drop basis terms whose |coefficient| < 2^keep (after rounding to 2^-fq).
"""
import importlib.util
import json
import os
import sys
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DEN = 1 << 60          # exact coefficient unit 2^-60


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------ polynomial in the circuit's basis
def h_basis(digits, carries, B_units, fq):
    """h(Phi - B) as a dict basis -> exact coefficient (int, units 2^-60).
    digits: list of j present; carries: list of (g, weight exponent l); B_units: bias in units 2^-fq.
    basis key: (frozenset of sigma indices XORed, frozenset of carry ids ANDed); () () = constant;
    the XOR set with one element is sigma_a itself."""
    atoms = [("s", j, 1 << j) for j in digits] + [("e", g, 1 << l) for g, l in carries]
    Bnum = Fraction(B_units, 1 << fq)
    # natural monomials: key (frozenset odd s indices, frozenset e ids) -> Fraction coefficient
    nat = {}

    def add(key, c):
        nat[key] = nat.get(key, Fraction(0)) + c

    lin = [(a[0], a[1], Fraction(a[2])) for a in atoms] + [("c", None, -Bnum)]

    def mono(items):
        s = set()
        e = set()
        for kind, idx, _ in items:
            if kind == "s":
                s ^= {idx}
            elif kind == "e":
                e.add(idx)
        return (frozenset(s), frozenset(e))

    n = len(lin)
    c2 = Fraction(3, 1 << 28)
    c3 = Fraction(-1, 1 << 53)
    for a in range(n):
        for b in range(n):
            add(mono([lin[a], lin[b]]), c2 * lin[a][2] * lin[b][2])
    for a in range(n):
        for b in range(n):
            ab = lin[a][2] * lin[b][2]
            for c in range(n):
                add(mono([lin[a], lin[b], lin[c]]), c3 * ab * lin[c][2])
    # to the circuit basis: prod s over S = 1 - 2 XOR(S) (S nonempty); times AND(E)
    out = {}

    def put(key, c):
        out[key] = out.get(key, Fraction(0)) + c

    for (S, E), c in nat.items():
        if c == 0:
            continue
        if not S:
            put((frozenset(), E), c)
        else:
            put((frozenset(), E), c)
            put((S, E), -2 * c)
    return {k: int(v * DEN) if (v * DEN).denominator == 1 else int(round(v * DEN)) for k, v in out.items() if v != 0}


def quant(coefs, fq, keep):
    """round every coefficient to 2^-fq (units of 2^-fq, int); drop |q| < 2^keep * 2^fq (non-constant terms)."""
    unit = DEN >> fq
    thr = (1 << fq) * (2.0 ** keep)
    q = {}
    for k, c in coefs.items():
        v = int(round(c / unit))
        if k != (frozenset(), frozenset()) and abs(v) < thr:
            continue
        if v:
            q[k] = v
    return q


def eval_basis(q, sig, car):
    """sum of q[key] * value(key) over all inputs; sig: dict j -> 0/1 array, car: dict g -> 0/1 array."""
    tot = None
    for (S, E), v in q.items():
        if S:
            x = None
            for j in S:
                x = sig[j].copy() if x is None else (x ^ sig[j])
        else:
            x = None
        if E:
            a = None
            for g in E:
                a = car[g].copy() if a is None else (a & car[g])
            x = a if x is None else (x & a)
        term = v * (x.astype(np.int64) if x is not None else 1)
        tot = term if tot is None else tot + term
    return tot


def run(y, prm):
    fq = int(prm.get("fq", 3))
    keep = float(prm.get("keep", -10))
    S = 1 << fq
    groups = prm.get("groups", [[18, 14], [13, 9], [8, 6]])
    biases = prm.get("biases")             # list of ints (units 2^-fq) or None: then chosen from the data
    n = len(y)
    r = y.astype(np.int64) * S
    sig = {}
    car = {}
    digits = []
    carries = []
    Bc = 0
    q_prev = {}
    val_prev = np.zeros(n, dtype=np.int64)
    info = []
    used_b = []
    for gi, (hi, lo) in enumerate(groups):
        for i in range(hi, lo - 1, -1):
            sg = (r < 0).astype(np.int64)
            sig[i] = sg
            digits.append(i)
            r = r - ((1 - 2 * sg) << i) * S
        q_now = quant(h_basis(digits, carries, Bc, fq), fq, keep)
        val_now = eval_basis(q_now, sig, car)
        r = r + (val_now - val_prev)
        lim = (1 << lo) * S
        c = (biases[gi] if biases else (int(-lim - r.min()) if r.min() < -lim else 0))
        used_b.append(int(c))
        r = r + c
        Bc += c
        e = (r >= lim).astype(np.int64)
        r = r - e * (2 * lim)
        ov = int(((r < -lim) | (r >= lim)).sum())
        car[gi] = e
        carries.append((gi, lo + 1))
        # the carry's h-change is deferred: the next correction subtracts Hq at the point WITHOUT the new
        # carry and with the old bias, i.e. q_prev must be Hq(point_g) as applied
        info.append(dict(g=gi, terms=len(q_now), new_terms=len(set(q_now) - set(q_prev)) if q_prev else len(q_now),
                         bias=int(c), ovf=ov, ones=int(e.sum())))
        q_prev, val_prev = q_now, val_now
    last = groups[-1][1]
    for i in range(last - 1, -1, -1):
        sg = (r < 0).astype(np.int64)
        sig[i] = sg
        digits.append(i)
        r = r - ((1 - 2 * sg) << i) * S
    q_now = quant(h_basis(digits, carries, Bc, fq), fq, keep)
    val_now = eval_basis(q_now, sig, car)
    r = r + (val_now - val_prev)
    info.append(dict(g="final", terms=len(q_now)))
    Phi = np.zeros(n, dtype=np.int64)
    for j in digits:
        Phi += (1 - 2 * sig[j]) << j
    for g, l in carries:
        Phi += car[g] << l
    U0 = Phi + np.floor_divide(r - sum(used_b), S)
    return U0, r, sig, car, info, used_b, q_now


def run_pd(y, prm):
    """per-digit variant (no carries): every digit is its own group, the correction is applied right after it,
    the remainder keeps one guard bit of headroom (positions -fq .. i+1 after digit i). Returns U0 (incl. the
    final remainder), the final remainder (units 2^-fq), the max |r| / 2^i ratio per step, per-step terms."""
    fq = int(prm.get("fq", 3))
    keep = float(prm.get("keep", -10))
    S = 1 << fq
    n = len(y)
    r = y.astype(np.int64) * S
    sig = {}
    digits = []
    q_prev = {}
    val_prev = np.zeros(n, dtype=np.int64)
    ratio = []
    terms = []
    top = int(prm.get("top", 18))
    lo = int(prm.get("lo", 0))
    for i in range(top, -1, -1):
        sg = (r < 0).astype(np.int64)
        sig[i] = sg
        digits.append(i)
        r = r - ((1 - 2 * sg) << i) * S
        if i >= lo:
            q_now = quant(h_basis(digits, [], 0, fq), fq, keep)
            val_now = eval_basis(q_now, sig, {})
            r = r + (val_now - val_prev)
            terms.append(len(set(k for k in set(q_now) | set(q_prev) if q_now.get(k, 0) != q_prev.get(k, 0))))
            q_prev, val_prev = q_now, val_now
        ratio.append(float(np.abs(r).max()) / float((1 << i) * S))
    Phi = np.zeros(n, dtype=np.int64)
    for j in digits:
        Phi += (1 - 2 * sig[j]) << j
    U0 = Phi + np.floor_divide(r, S)
    return U0, r, ratio, terms


def count(prm):
    """exact adder count of the construction from its term lists (no data needed except the biases):
    per correction the operand is sum over changed basis terms of dq * bit (+ constant + bias); positive and
    negative parts are added separately (r - N = NOT(NOT r + N)); each part is packed greedily into layers of
    distinct bit positions; a layer costs (top - lowest) ANDs with the two-carry ripple adder (positions from
    its lowest bit to the register top, top = lo + 1 for a group (guard / carry wire), 1 for the final);
    every term that ANDs a carry costs 1 AND to form (erased by measurement for free).
    Returns per-correction dicts and the total."""
    fq = int(prm.get("fq", 3))
    keep = float(prm.get("keep", -10))
    groups = prm.get("groups", [[18, 14], [13, 9], [8, 6]])
    biases = prm["biases"]
    digits, carries = [], []
    Bc = 0
    q_prev = {}
    out = []
    tot = 0
    points = []
    for gi, (hi, lo) in enumerate(groups):
        digits += list(range(hi, lo - 1, -1))
        points.append((list(digits), list(carries), Bc, lo + 1))
        Bc += biases[gi]
        carries.append((gi, lo + 1))
    digits += list(range(groups[-1][1] - 1, -1, -1))
    points.append((list(digits), list(carries), Bc, 1))
    for (dg, cr, B, top) in points:
        q_now = quant(h_basis(dg, cr, B, fq), fq, keep)
        keys = set(q_now) | set(q_prev)
        dq = {k: q_now.get(k, 0) - q_prev.get(k, 0) for k in keys}
        dq = {k: v for k, v in dq.items() if v}
        top_pos = top + fq                       # register top in units-of-2^-fq bit positions
        pos_bits, neg_bits = [], []
        n_and = 0
        for k, v in dq.items():
            if k == (frozenset(), frozenset()):
                continue                          # constant: merged into a layer (no extra cost)
            if k[1]:
                n_and += 1
            a = abs(v)
            b = 0
            while a:
                if a & 1:
                    (pos_bits if v > 0 else neg_bits).append(b)
                a >>= 1
                b += 1
        cost = 0
        layers = 0
        for bits in (pos_bits, neg_bits):
            cnt = {}
            for b in bits:
                cnt[b] = cnt.get(b, 0) + 1
            while cnt:
                lowest = min(cnt)
                cost += max(0, top_pos - lowest)
                layers += 1
                for b in list(cnt):
                    cnt[b] -= 1
                    if cnt[b] == 0:
                        del cnt[b]
        cost += n_and
        tot += cost
        out.append(dict(top=top, terms=len(dq), pos_bits=len(pos_bits), neg_bits=len(neg_bits), layers=layers,
                        and_terms=n_and, ands=cost))
        q_prev = q_now
    return out, tot


def main():
    wd, cmd = sys.argv[1], sys.argv[2]
    prm = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}
    if cmd == "count":
        out, tot = count(prm)
        for o in out:
            print(o)
        print(f"TOTAL ANDs one way (greedy layers, two-carry ripple): {tot}")
        return
    INJ = _load("fcinj", "invsqrt_fc_inj.py")
    CM = INJ.RP.CM
    best = dict(INJ.BEST)
    cm, U1, Ug, _ = INJ.words(wd, best)
    ideal = cm["Us"]
    if cmd == "pd":
        U0, r, ratio, terms = run_pd(U1, prm)
        e0 = U0 - ideal
        K = int(np.round(-(e0.max() + e0.min()) / 2.0))
        U = U0 + K
        e = U - ideal
        fq = int(prm.get("fq", 3))
        rint = np.floor_divide(r, 1 << fq)
        rv, rc = np.unique(rint, return_counts=True)
        lim = 1 << (26 - 7)
        Uc = np.clip(U, -lim, lim - 1)
        p = CM.RT.divide_m(cm, Uc, b_t=int(best.get("bt", 5)))
        res = CM.CORE.evaluate_word(CM._STATE["ph"], CM._STATE["ch"], p, lp=True)
        print(f"pd {prm}: max|r|/2^i per step {[round(x, 3) for x in ratio]}", flush=True)
        print(f"  changed terms per step {terms} (sum {sum(terms)}); K {K}; final int remainder "
              f"{list(zip(rv.tolist(), rc.tolist()))}; U - ideal {e.min():.3f} .. {e.max():.3f}", flush=True)
        print(f"  LP margin {res['m']:.4g}; clipped {int((U != Uc).sum())}", flush=True)
        print("RESULT " + json.dumps(dict(prm=prm, K=K, margin=float(res["m"]), e=[float(e.min()), float(e.max())],
                                          rvals=rv.tolist(), ratio=ratio, terms=terms)), flush=True)
        return
    U0, r, sig, car, info, used_b, qf = run(U1, prm)
    e0 = U0 - ideal
    K = int(np.round(-(e0.max() + e0.min()) / 2.0))
    U = U0 + K
    e = U - ideal
    fq = int(prm.get("fq", 3))
    rint = np.floor_divide(r, 1 << fq)
    rv, rc = np.unique(rint, return_counts=True)
    print(f"prm {prm}: groups {info}", flush=True)
    print(f"biases {used_b}; K {K}; final int remainder {list(zip(rv.tolist(), rc.tolist()))[:8]} ({len(rv)}); "
          f"U - ideal {e.min():.3f} .. {e.max():.3f}", flush=True)
    lim = 1 << (26 - 7)
    Uc = np.clip(U, -lim, lim - 1)
    p = CM.RT.divide_m(cm, Uc, b_t=int(best.get("bt", 5)))
    res = CM.CORE.evaluate_word(CM._STATE["ph"], CM._STATE["ch"], p, lp=True)
    st = np.stack([U] + [car[g] for g in sorted(car)] + [rint])
    print(f"LP margin {res['m']:.4g}; clipped {int((U != Uc).sum())}; distinct (U, carries, r_int) "
          f"{np.unique(st, axis=1).shape[1]} / distinct y {len(np.unique(U1))}", flush=True)
    print("RESULT " + json.dumps(dict(prm=prm, info=info, biases=used_b, K=K, e=[float(e.min()), float(e.max())],
                                      margin=float(res["m"]), rvals=rv.tolist(), final_terms=len(qf))), flush=True)


if __name__ == "__main__":
    main()
