#!/usr/bin/env python3
"""r33_inv_wire2_rp_tail.py -- re-pricer of R33 seat INV-WIRE2: what the tapered rows of the two correction
divisions cost at the precision they need (own code; the designer's tools and S1 / S2 imported read-only).

The design prices its tapered rows (divisor T_k = rhu(d / 2^k), error <= 2^(k-1) in the scaled frame) at the
line's truncated-tail price "m - k" (its reading -2). The line's tail row with kk dropped bits uses the divisor
2A + 1, A = V rounded to a multiple of 2^(kk+1) minus v0 (R31 rp_ovf, invsqrt_s2_root "ship"), whose error is up
to about 2^(kk+1), and costs n - 1 = m - kk [M F1]. This tool runs the two correction divisions with the line's
OWN tail form (scaled frame r <- 2 r - s dt, register u = floor(r / 2) in m + 2 wires, overflow counted per row,
then the final sign as the last digit, as the line's s36) and the design's downstream unchanged (higher-order
stage, main divider, angle LP on all 3,133,452 inputs), so that every price comes with the precision it buys.

  r33_inv_wire2_rp_tail.py WD one  '{"m":7, ..., "form":"line", "sched":[1,2,3,4,5,6,7,7]}'      one set
  r33_inv_wire2_rp_tail.py WD grid OUT.jsonl NPROC SETS.json                                       many sets
  r33_inv_wire2_rp_tail.py WD reads '{"m":7, ...}'      are the tail digits CNOT reads of the final u (frames)?
  r33_inv_wire2_rp_tail.py - row M_LO M_HI              AND count of a rounded tapered row with the REAL inputs
  r33_inv_wire2_rp_tail.py - price OUT.jsonl            re-price result lines (report law; tail at m - kk)
form "round" = the design's nr_tapered unchanged (its prices: reading -2 and middle are both printed).
"""
import importlib.util
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CM = _load("cm", "r33_inv-wire2_cm.py")          # the designer's model (read-only)
REP = _load("rep", "r33_inv-wire2_report.py")    # the designer's re-pricing (read-only)
ORIG_NR = CM.nr_tapered

CFG = dict(form="round", sched=None, ovf=[], tail_rr=[])


def nr_line(reg, d, a_top, r, F, mbits, adj0=0):
    """exact rows exactly as the design's nr_tapered, then len(sched) tail rows in the line's form, then the
    final sign as one more digit. Returns what nr_tapered returns (the final register is the window, which is
    entirely digits: zeros returned)."""
    m = mbits
    rows, digs = [], []
    for i, t in enumerate(range(a_top, -1, -1)):
        T = d << t
        b_reg = CM.bits_signed(reg) - (adj0 if i == 0 else 0)
        b_T = int(T.max()).bit_length() + 1
        virgin = min(max(t, 0), F)
        span = max(b_reg, b_T) - t
        live = max(b_reg, b_T) - virgin + i
        s = np.where(reg < 0, -1, 1).astype(np.int64)
        reg = reg - s * T
        digs.append(s)
        rows.append(dict(t=t, span=int(span), live=int(live)))
    n_ex = len(digs)
    rr = reg.copy()
    V = d >> 1
    lim = 1 << (m + 1)
    ovf = []
    rr_hist = []
    # division 1 is the first call of a word, division 2 the second ("sched1" / "sched2" override "sched")
    which = len(CFG["ovf"]) % 2
    sched = CFG.get("sched1" if which == 0 else "sched2") or CFG["sched"]
    CFG["used"].append(list(sched))
    for j, kk in enumerate(sched):
        s = np.where(rr < 0, -1, 1).astype(np.int64)
        if kk:
            A = (((V >> (kk + 1)) + ((V >> kk) & 1)) << (kk + 1)) - (V & 1)
            dt = 2 * A + 1
        else:
            dt = d
        rr = 2 * rr - s * dt
        u = rr >> 1
        ovf.append(int(((u < -lim) | (u >= lim)).sum()))
        digs.append(s)
        rr_hist.append(rr.copy())
        # price convention for the report's tapered-row law (PESS = -2: span - 4): span = m - kk + 4 -> m - kk
        rows.append(dict(t=-(j + 1), span=int(max(0, m - kk) + 4), live=int(n_ex + m + 2), kk=int(kk)))
    s = np.where(rr < 0, -1, 1).astype(np.int64)
    digs.append(s)
    rows.append(dict(t=-(len(sched) + 1), span=0, live=int(n_ex + m + 2), kk=None))
    CFG["ovf"].append(ovf)
    CFG["tail_rr"].append(rr_hist)
    n = len(digs)
    q = np.zeros(len(reg), dtype=np.int64)
    for j, sj in enumerate(digs):
        q += sj << (n - 1 - j)
    return q, np.zeros_like(reg), rows, digs[-1]


ORIG_WORD = CM.word_cm


def word_cm_2c(cm, prm):
    """ho "grpt2c": the design's grouped recurrence (mode grpt) with the square built from the TWO'S COMPLEMENT
    bits of F (top bit at negative weight), as a circuit holding F holds it, instead of the bits of |F| (which no
    stage computes). Same pair rule (e < f, e + f >= thr; diagonal if 2e >= thr - 1), same cubic and quartic."""
    if prm.get("ho") != "grpt2c":
        return ORIG_WORD(cm, prm)
    m = cm["m"]
    p0 = dict(prm)
    p0["ho"] = "none"
    p0["c0"] = 0
    U1, info = ORIG_WORD(cm, p0)
    assert info["n_clip"] == 0
    c0 = int(prm.get("c0", 0))
    kg = int(prm.get("kg", max(2, m - 2)))
    thr = int(prm.get("thr", 26))
    t_ho = int(prm.get("t", 13))
    Lb = 27 - m
    y = U1.astype(np.int64)
    F = np.zeros_like(y)
    R = y.copy()
    hi = Lb
    lo_end = max(0, Lb - t_ho)
    groups = []
    wf_max = 0
    while hi > lo_end:
        lo = max(hi - kg, lo_end)
        D = R >> lo
        F = F + (D << lo)
        WF = CM.bits_signed(F)
        wf_max = max(wf_max, WF)
        top = WF - 1
        w = [(1 << e) for e in range(top)] + [-(1 << top)]
        bits = [(F >> e) & 1 for e in range(top)] + [(F < 0).astype(np.int64)]
        S = np.zeros_like(F)
        for e in range(WF):
            if 2 * e >= thr - 1:
                S += bits[e] * (w[e] * w[e])
            for f in range(e + 1, WF):
                if e + f >= thr:
                    S += (bits[e] & bits[f]) * (2 * w[e] * w[f])
        Ff = F.astype(np.float64)
        Gt = Ff - 3.0 * S.astype(np.float64) / 2.0 ** 28 + Ff ** 3 / 2.0 ** 53 - 5.0 * Ff ** 4 / 2.0 ** 82
        R = y - np.round(Gt).astype(np.int64)
        groups.append([int(lo), int(hi), int(np.abs(D).max()).bit_length(), CM.bits_signed(R), int(WF)])
        hi = lo
    Ug = F + R
    npairs = sum(1 for e in range(lo_end, wf_max) for f in range(e + 1, wf_max) if e + f >= thr)
    q_int = np.round((Ug - U1).astype(np.float64)).astype(np.int64)
    U = c0 + U1 + q_int
    lim = 1 << (26 - m)
    n_clip = int(((U < -lim) | (U >= lim)).sum())
    U = np.clip(U, -lim, lim - 1)
    ue = U - cm["Us"]
    info.update(ho="grpt2c", t=t_ho, c0=c0, n_clip=n_clip, U_err=[float(ue.min()), float(ue.max())],
                grp=dict(kg=kg, groups=groups, pairs={thr: npairs}, lo_end=int(lo_end), thr=thr, wf=int(wf_max)),
                q_max=float((Ug - U1).max()), b_U=CM.bits_signed(U))
    return U, info


def set_form(prm):
    CFG["form"] = prm.get("form", "round")
    CFG["sched"] = prm.get("sched")
    CFG["sched1"] = prm.get("sched1")
    CFG["sched2"] = prm.get("sched2")
    CFG["ovf"] = []
    CFG["used"] = []
    CFG["tail_rr"] = []
    CM.nr_tapered = nr_line if CFG["form"] == "line" else ORIG_NR
    CM.word_cm = word_cm_2c


def price(out, pess):
    """the report's re-pricing with tapered rows read at PESS; for form 'line' only PESS = -2 is the line's law."""
    REP.PESS = pess
    st, tof, cost, tr = REP.reprice(out)
    return st, tof, cost, tr


def run(prm):
    set_form(prm)
    out = CM.run_one({k: v for k, v in prm.items() if k not in ("form", "sched", "sched1", "sched2", "tag")})
    out["prm"] = prm
    out["ovf"] = CFG["ovf"]
    out["used"] = CFG["used"]
    res = {}
    for pess in ((-2, 0) if CFG["form"] == "round" else (-2,)):
        st, tof, cost, tr = price(out, pess)
        res[str(pess)] = dict(tof=round(tof, 2), cost=(round(cost, 1) if math.isfinite(cost) else None),
                              gain=(round(100 * (1 - cost / REP.BASE), 2) if math.isfinite(cost) else None),
                              tr=(round(tr, 3) if tr == tr else None), DIV1=st["DIV1"], DIV1_INV=st["DIV1_INV"],
                              DIV2=st["DIV2"], DIV2_INV=st["DIV2_INV"], HO=st["HO"])
    out["price"] = res
    return out


def _worker(args):
    wd, prm = args
    if "ch" not in CM._STATE:
        CM.setup(wd)
    try:
        return run(prm)
    except Exception as e:  # a set that does not build is a result too
        return dict(prm=prm, error=repr(e))


# ------------------------------------------------------------------------------------- reads (frames)
def reads(wd, prm):
    """for the line form: is each tail digit (and the final sign) a CNOT read of the final window u (m + 2 bits)
    and at most one root bit (s = u[b] ^ v[c] ^ const)?  Division 1 only (division 2 is the same module)."""
    CM.setup(wd)
    set_form(prm)
    m = int(prm["m"])
    cm = CM.chain(CM._STATE["ch"], m)
    U, info = CM.word_cm(cm, {k: v for k, v in prm.items() if k not in ("form", "sched", "sched1", "sched2", "tag")})
    # recompute division 1 to get digits and the final window
    d = cm["d"]
    reg1 = cm["delta"] << int(prm.get("F1", 5))
    a1 = CM.a_top_of(reg1, d)
    CFG["tail_rr"] = []
    q, _, rows, _ = nr_line(reg1, d, a1, 0, int(prm.get("F1", 5)), m)
    rr = CFG["tail_rr"][-1][-1]
    u = (rr >> 1) & ((1 << (m + 2)) - 1)
    n = len(rows)
    digs = [((q >> (n - 1 - j)) & 1) for j in range(n)]   # not used: recompute signs directly below
    # signs of the tail rows and the final sign, recomputed from the history
    hist = CFG["tail_rr"][-1]
    n_ex = a1 + 1
    v = cm["v"]
    feats = [((u >> b) & 1) for b in range(m + 2)]
    vb = [np.zeros_like(v)] + [((v >> c) & 1) for c in range(m)]
    out = []
    prev = reg1.copy()
    # digit of tail row j = sign of the register before it
    befores = [None] * len(hist)
    # the register before tail row 0 is the register after the exact rows: rebuild it
    regx = reg1.copy()
    for t in range(a1, -1, -1):
        s = np.where(regx < 0, -1, 1)
        regx = regx - s * (d << t)
    befores[0] = regx
    for j in range(1, len(hist)):
        befores[j] = hist[j - 1]
    sig = [(b < 0).astype(np.int64) for b in befores] + [(hist[-1] < 0).astype(np.int64)]
    for j, sj in enumerate(sig):
        found = None
        for b in range(m + 2):
            for c in range(m + 1):
                x = feats[b] ^ vb[c]
                if np.array_equal(x, sj):
                    found = (b, c - 1, 0)
                elif np.array_equal(x ^ 1, sj):
                    found = (b, c - 1, 1)
                if found:
                    break
            if found:
                break
        out.append(found)
    print(f"m {m} sched {CFG['sched']}: tail digit j -> (u bit, v bit or -1, const) or None (needs a frame with ANDs)")
    for j, f in enumerate(out):
        print(f"  tail {j}{' (final sign)' if j == len(out) - 1 else ''}: {f}")
    print(f"  overflow per tail row: {CFG['ovf'][-1] if CFG['ovf'] else None}")


# ------------------------------------------------------------------------------------- row law, real inputs
def gf2_affine(X, y):
    """is y an affine function of the columns of X over GF(2)? (vectorised elimination; same test as the
    design's rowlaw.gf2_affine)"""
    A = np.concatenate([X, np.ones((X.shape[0], 1), dtype=np.uint8), y[:, None]], axis=1).astype(np.uint8)
    A = np.unique(A, axis=0)
    k = A.shape[1] - 1
    M = A.copy()
    r = 0
    for c in range(k):
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        piv = r + nz[0]
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        sel = M[:, c].astype(bool)
        sel[r] = False
        M[sel] ^= M[r]
        r += 1
        if r == M.shape[0]:
            break
    bad = (M[:, :k].sum(axis=1) == 0) & (M[:, k] == 1)
    return not bool(bad.any())


def row_check(m):
    """tapered row k of the design: res = reg - s T, T = rhu(d / 2^k) = T' + r, T' = d >> k, r = bit k-1 of d.
    On wires: reg, v, s. The adder adds (T' xor S) with carry-in (r xor sbit) -- T is NOT on wires. For each k:
    the smallest p such that, with ripple carries c_1..c_p computed, every result bit at or above p is affine in
    (register bits at or above p, c_p, sbit); and how many of c_1..c_p are themselves affine in the inputs
    (free). ANDs = p - free. Also the same with T precomputed (the design's rowlaw), for comparison."""
    W = m + 5
    mask = (1 << W) - 1
    res_line = []
    for k in range(1, m + 3):
        rows = []
        for Q in range(1 << m, 1 << (m + 1)):
            d = 2 * Q + 1
            T = max(1, (d + (1 << (k - 1))) >> k)
            Tp = max(1, (d + (1 << (k - 2))) >> (k - 1)) if k >= 2 else d
            reg = np.arange(-Tp, Tp + 1, dtype=np.int64)
            s = np.where(reg >= 0, 1, -1)
            res = reg - s * T
            rows.append(np.stack([np.full_like(reg, d), reg, s, res, np.full_like(reg, T)], axis=1))
        A = np.concatenate(rows)
        d, reg, s, res, T = A[:, 0], A[:, 1], A[:, 2], A[:, 3], A[:, 4]
        v = (d - 1) // 2 - (1 << m)
        ru = reg & mask
        sbit = (s > 0).astype(np.int64)                  # subtract when s = +1
        Tq = d >> k
        rbit = (d >> (k - 1)) & 1
        b = np.zeros((W, len(ru)), dtype=np.int64)
        for i in range(W):
            b[i] = ((Tq >> i) & 1) ^ sbit
        c = np.zeros((W + 1, len(ru)), dtype=np.int64)
        c[0] = rbit ^ sbit
        for i in range(W):
            a_i = (ru >> i) & 1
            c[i + 1] = (a_i & b[i]) | (a_i & c[i]) | (b[i] & c[i])
        chk = np.zeros_like(ru)
        for i in range(W):
            chk |= ((((ru >> i) & 1) ^ b[i] ^ c[i]) << i)
        assert np.array_equal(chk, res & mask), (m, k)
        inputs = [((ru >> i) & 1) for i in range(W)] + [((v >> i) & 1) for i in range(m)] + [sbit]
        Xin = np.stack(inputs, axis=1).astype(np.uint8)
        best = None
        for p in range(0, W):
            feats = [((ru >> i) & 1) for i in range(p, W)] + [c[p], sbit]
            X = np.stack(feats, axis=1).astype(np.uint8)
            ok = True
            for i in range(p, W):
                y = ((res & mask) >> i) & 1
                ok &= gf2_affine(X, y.astype(np.uint8))
                if not ok:
                    break
            if ok:
                best = p
                break
        free = sum(1 for i in range(1, best + 1) if gf2_affine(Xin, c[i].astype(np.uint8)))
        res_line.append((k, best, free, best - free, m + 2 - k, max(0, m - k)))
    return res_line


def main():
    wd, cmd = sys.argv[1], sys.argv[2]
    if cmd == "one":
        CM.setup(wd)
        prm = json.loads(sys.argv[3])
        print("RESULT " + json.dumps(run(prm)), flush=True)
    elif cmd == "grid":
        out_path, nproc, sets_path = sys.argv[3], int(sys.argv[4]), sys.argv[5]
        sets = json.load(open(sets_path))
        import multiprocessing as mp
        with mp.Pool(nproc) as pool, open(out_path, "a") as f:
            for res in pool.imap_unordered(_worker, [(wd, s) for s in sets]):
                f.write(json.dumps(res) + "\n")
                f.flush()
        print(f"done {len(sets)} sets", flush=True)
    elif cmd == "reads":
        reads(wd, json.loads(sys.argv[3]))
    elif cmd == "row":
        lo, hi = int(sys.argv[3]), int(sys.argv[4])
        print("k: carries needed p, of them free, ANDs, design middle (m + 2 - k), design reading -2 (m - k)")
        for m in range(lo, hi + 1):
            r = row_check(m)
            print(f"m {m}: " + "  ".join(f"k{k}: p{p} free{f} AND {a} (mid {mid}, -2 {lo2})" for k, p, f, a, mid, lo2 in r),
                  flush=True)
            print(f"   sum over k = 1..{m + 2}: ANDs {sum(x[3] for x in r)}, middle {sum(x[4] for x in r)}, "
                  f"reading -2 {sum(x[5] for x in r)}", flush=True)
    elif cmd == "price":
        for line in open(sys.argv[3]):
            if not line.strip():
                continue
            r = json.loads(line)
            if "error" in r:
                print("ERROR", r["prm"], r["error"])
                continue
            p = r["prm"]
            print(f"{json.dumps(p)}  q {r['q']} margin {r['margin']:.4g} ovf {r.get('ovf')} clip {r['n_clip']} "
                  f"price {json.dumps(r['price'])}")


if __name__ == "__main__":
    main()
