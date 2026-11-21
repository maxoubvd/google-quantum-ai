#!/usr/bin/env python3
"""r54_inv_tail_dc.py -- R54-INV-TAIL (designer): THE DON'T-CARE SET of a divider tail's last j rows, exact, on ALL
3,133,452 inputs, in the R37 / R49 integer model of k9f (read-only imports: tools/r37_inv_40_model.py; the chain reads
MY copy of FC's source, argv[1]).

The last j rows of a tail (tail rows 9-j .. 8 and the final sign) are a map of the window H = rr >> shift and the top
bits of v to the last j + 1 digits.  For every input x and every word w of those j + 1 digits (2^(j+1) words) the rest
of the chain is run exactly (DIV1: DIV2, stage table, MAIN, ring; DIV2: stage table, MAIN, ring) with TODAY's converged
angles (LP of the unchanged circuit): slack(x, w) = tol(x) - |err(x, w)|.  Per class c = (H, v top bits):
minslack[c, w] = min over the inputs of the class.  "w allowed on c at margin M0" iff minslack[c, w] >= M0.

  mode table:  r54_inv_tail_dc.py SRC WD HTAB.npz OUT.npz table DIV1|DIV2|MAIN j [NPROC] [extra=K]
                    extra=K: the class also reads K more low bits of v (finer classes: more freedom for a program)
  mode eval:   r54_inv_tail_dc.py SRC WD HTAB.npz OUT.json eval TABLE.npz PROG.json    (the program applied on all
                    inputs, the LP re-solved: the true margin of the candidate)
  mode sweep:  r54_inv_tail_dc.py SRC WD HTAB.npz OUT.jsonl sweep 'JSON list of [name, sched1, sched2, mt]'
                    the exact chain with other tail schedules (mt = MAIN's tapers {row b: taper}, null = k9f's), LP
                    re-solved: the price list of the schedule knobs.  DIV2's final sign follows the circuit's 2-bit rule
                    only when its schedule ends 6, 6 (other endings: G.ref_div's, which R49 checked on its level rows).
Safe to run: reads only; writes OUT under my own folder."""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r37_inv_40_model as RM          # noqa: E402  (reads sys.argv[1] = SRC at import)

F8, G, SM = RM.F8, RM.G, RM.SM
M_, BT, BLAST, SLOT, RING, A_TOP = 7, 6, 1, "zero", 26, 10
S1 = [0, 0, 0, 2, 3, 4, 5, 6]
S2 = [0, 0, 1, 2, 3, 4, 6, 6]


def div_trace(reg, d, a_top, sched):
    """G.ref_div with the remainder before every tail row kept."""
    reg = reg.astype(np.int64).copy()
    sg = []
    for t in range(a_top, -1, -1):
        neg = reg < 0
        sg.append(neg.astype(np.int64))
        reg = reg - np.where(neg, -1, 1) * (d << t)
    rr = reg
    rrs = []
    for kk in sched:
        rrs.append(rr.copy())
        neg = rr < 0
        sg.append(neg.astype(np.int64))
        rr = 2 * rr - np.where(neg, -1, 1) * F8.tail_dt(d, kk)
    sg.append((rr < 0).astype(np.int64))
    n = len(sg)
    q = np.zeros(len(d), dtype=np.int64)
    for j, g in enumerate(sg):
        q += (1 - 2 * g) << (n - 1 - j)
    return sg, q, rrs


def div2_any(q1, d, sched):
    """DIV2 with any schedule: the circuit's 2-bit final sign when it ends 6, 6, G.ref_div's otherwise."""
    if list(sched) == S2:
        return div2_true(q1, d)
    sg, P, rrs = div_trace(q1, d, A_TOP, sched)
    if list(sched[-2:]) == [6, 6]:
        M8 = rrs[7] >> 7
        A = 1 + ((d >> 7) & 1)
        res = np.where(M8 < 0, (M8 & 3) + A, (M8 & 3) - A)
        fin_t = (res >> 1) & 1
        P = P + 2 * (sg[-1] - fin_t)
        sg[-1] = fin_t.astype(np.int64)
    return sg, P, rrs


def sweep_one(spec, wd, hpath):
    """LP margin of the exact chain with the tail schedules of `spec` = [name, sched1, sched2, mt]."""
    global MT
    name, s1, s2, mt = spec
    z = np.load(os.path.join(wd, "chain.npz"))
    tgt, tol = z["target"].astype(np.uint64), z["tol"].astype(np.float64)
    x = np.arange(1, RM.MAXX + 1, dtype=np.int64)
    n, k, Q_, d, delta = F8.ref_chain(x, M_)
    _, q1, _ = div_trace(delta.astype(np.int64) << 5, d, A_TOP, s1)
    _, P, _ = div2_any(q1, d, s2)
    h, seen = htab(hpath)
    y = -((P + 1) >> 1)
    T = y >> 9
    out = dict(name=name, sched1=s1, sched2=s2, mt=mt)
    if (T + 1024).min() < 0 or (T + 1024).max() > 2047 or not seen[T + 1024].all():
        out.update(margin=None, note="T leaves the stage's table")
        return out
    U = y + h[T + 1024]
    keep = dict(MT)
    if mt:
        MT = {int(b): int(t) for b, t in mt.items()}
    try:
        if MT.get(1) == 6:
            Y = main_trace(U, d)[0]
        else:
            Y = RM.main_rows(U, d, M_, BT, BLAST, SLOT)       # the plain model (no XT): rows b 5..1 at tapers 1..5
    finally:
        MT = keep
    if Y.min() < (1 << 26) or Y.max() >= (1 << 27):
        out.update(margin=None, note="word leaves [2^26, 2^27)")
        return out
    fire, ang, fixed = ring_fast(Y, k, RING)
    tot = fixed.copy()
    for a_, f_ in zip(ang, fire):
        tot += np.where(f_ == 1, np.uint64(a_), np.uint64(0))
    err = (tot - tgt).view(np.int64).astype(np.float64)
    mg, conv, bound, dv = lp_dv(err, tol, fire)
    out.update(margin=mg, conv=conv)
    return out


def div2_true(q1, d):
    """DIV2 as the CIRCUIT computes it: G.ref_div with sched2 [.., 6, 6], except the final sign: the last row adds
    1 + v6 into TWO wires and the final sign is the top of that 2-bit result (the integer model keeps the full value:
    R49 measured 'the last digit differs from the model on 2.8%').  Checked in base(): y = the circuit's y."""
    assert S2[-2:] == [6, 6]
    sg, P, rrs = div_trace(q1, d, A_TOP, S2)
    M8 = rrs[7] >> 7
    su = M8 < 0
    A = 1 + ((d >> 7) & 1)
    L = M8 & 3
    res = np.where(su, L + A, L - A)
    fin_t = (res >> 1) & 1
    P_true = P + 2 * (sg[-1] - fin_t)
    sg[-1] = fin_t.astype(np.int64)
    return sg, P_true, rrs


MT = {5: 1, 4: 2, 3: 3, 2: 4, 1: 6}     # MAIN's tapered rows on k9f: bt 6, and the chain's R41_XT=MAIN runs the row
                                        # at taper 5 (b = 1) at taper 6 on TWO wires (its lowest target is skipped)


def main_trace(U, d):
    """MAIN as the CIRCUIT computes it (blast 1): RM.main_rows with the last row (b = 1) at taper 6 (R41_XT) and the
    final sign = the top of that row's 2-bit result; the remainder before the rows b <= 5 kept.  Checked against
    the circuit's word (tools/r54_inv_tail_dump.py): see the calibration line in the design file."""
    N = (np.int64(1) << (28 + M_)) + (U.astype(np.int64) << (M_ + 1))
    R = N.copy()
    Qw = np.zeros(len(d), dtype=np.int64)
    Rs, sg = {}, {}
    fin = None
    for b in range(26, BLAST - 1, -1):
        t = MT.get(b, 0)
        neg = R < 0
        if b <= 5:
            Rs[b] = R.copy()
            sg[b] = neg.astype(np.int64)
        s_ = np.where(neg, -1, 1)
        Qw = Qw + s_ * (1 << b)
        if b == 1:
            M = R >> 8                               # 4-bit window; the row adds 1 + v6 into its bits 1, 2
            A = 1 + ((d >> 7) & 1)
            L = (M >> 1) & 3
            res = np.where(neg, L + A, L - A)
            fin = ((res >> 1) & 1).astype(np.int64)
        else:
            dt = F8.tail_dt(d, t) if t else d
            R = R - s_ * (dt << b)
    Y = Qw - (fin << BLAST)
    Y = Y - (Y & ((1 << BLAST) - 1))
    return Y, Rs, sg, fin


def wvals(which, w, j):
    """value of a word of j + 1 bits (MSB = the first digit; int or array).  Dividers: sum (1 - 2 g_i) 2^(j-i);
    MAIN (rows b = j..1 and the final sign): sum_{b} (1 - 2 g_b) 2^b - 2 f."""
    if which == "MAIN":
        v = -2 * (w & 1)
        for i in range(j):
            v = v + ((1 - 2 * ((w >> (j - i)) & 1)) << (j - i))
        return v
    v = 0
    for i in range(j + 1):
        v = v + ((1 - 2 * ((w >> (j - i)) & 1)) << (j - i))
    return v


KEY = {"DIV1": "q1", "DIV2": "P", "MAIN": "Y"}
PROGS = {}          # programs of OTHER tails (joint fix rounds): set by tools/r54_inv_tail_alt.py before any run


def ring_fast(Y, k, N):
    """RM.ring_rows_n, vectorised (checked equal once in main)."""
    sh = 27 - N
    M = 1 << 64
    fire, ang = [], []
    for p in range(N):
        i = (p - k) % N
        fire.append(((Y >> (i + sh)) & 1).astype(np.uint8))
        ang.append(1 << (p + 27 + sh))
    fixed = np.zeros(len(Y), dtype=np.uint64)
    for j in range(10):
        on = (k >= j + 1)
        i = np.where(on, j - k + N, 0)
        bits = (on & (((Y >> (i + sh)) & 1) == 1)).astype(np.uint8)
        if j <= 6:
            fire.append(bits)
            ang.append(((1 << (j + 54)) - (1 << (j + 27 + sh))) % M)
        else:
            fixed += np.where(bits == 1, np.uint64((1 << (54 + j)) % M), np.uint64(0))
    return np.array(fire, dtype=np.uint8), ang, fixed


def lp_dv(err, tol, fire, rounds=400):
    """RM.lp_conv, returning the offsets of the best round too."""
    from scipy.optimize import linprog
    K = fire.shape[0]
    m0 = tol - np.abs(err)
    start = add = 4000
    act = np.unique(np.argpartition(m0, start)[:start])
    best, bestdv = float(m0.min()), np.zeros(K)
    bound, conv = float("inf"), False
    for rnd in range(rounds):
        Bm = fire[:, act].T.astype(np.float64)
        A = np.vstack([np.hstack([Bm, np.ones((len(act), 1))]), np.hstack([-Bm, np.ones((len(act), 1))])])
        b = np.concatenate([tol[act] - err[act], tol[act] + err[act]])
        cvec = np.zeros(K + 1)
        cvec[-1] = -1.0
        res = linprog(cvec, A_ub=A, b_ub=b, bounds=[(-2.0 ** 44, 2.0 ** 44)] * K + [(None, None)], method="highs")
        if res.status != 0:
            break
        bound = float(res.x[-1])
        dv = np.round(res.x[:K])
        marg = slack_of(err, tol, fire, dv)
        mmin = float(marg.min())
        if mmin > best:
            best, bestdv = mmin, dv.copy()
        if mmin >= bound - 2 * K:
            conv = True
            break
        cand = np.argpartition(marg, add)[:add]
        cand = cand[marg[cand] < bound - 2 * K]
        new = np.setdiff1d(cand, act, assume_unique=False)
        if len(new) == 0:
            conv = True
            break
        act = np.concatenate([act, new])
    return best, conv, bound, bestdv


def slack_of(err, tol, fire, dv):
    marg = np.empty(len(err), dtype=np.float64)
    if SM.HAVE_NUMBA:
        SM._margin_kernel(err, tol, fire, dv.astype(np.float64), marg)
    else:
        e = err.copy()
        for r in np.flatnonzero(dv):
            e += dv[r] * fire[r]
        marg = tol - np.abs(e)
    return marg


_H = {}


def htab(path):
    if path not in _H:
        z = np.load(path)
        h = z["h"].astype(np.int64).copy()
        seen = z["seen"].astype(bool)
        idx = np.arange(2048)
        good = idx[seen]
        near = good[np.argmin(np.abs(idx[:, None] - good[None, :]), axis=1)]
        h[~seen] = h[near[~seen]]
        _H[path] = (h, seen)
    return _H[path]


def downstream(stage, word, d, k, hpath, tgt):
    """stage 'q1': word = DIV1's quotient -> DIV2 -> y; stage 'P': word = DIV2's quotient -> y.
    Returns (err, fire, ok): ok False where T leaves the table or the word leaves [2^26, 2^27)."""
    h, seen = htab(hpath)
    if stage == "Y":
        Y = word
        bad = np.zeros(len(d), dtype=bool)
    else:
        if stage == "q1":
            sg2, P, rr2 = div2_true(word, d)
            if "DIV2" in PROGS:
                P = apply_to(dict(d=d, P=P, sg2=sg2, rr2=rr2), "DIV2", PROGS["DIV2"])
        else:
            P = word
        y = -((P + 1) >> 1)
        T = y >> 9
        bad = (T + 1024 < 0) | (T + 1024 > 2047)
        Tc = np.clip(T + 1024, 0, 2047)
        bad |= ~seen[Tc]
        U = y + h[Tc]
        if "MAIN" in PROGS:
            Y, Rs, sgm, fin = main_trace(U, d)
            Y = apply_to(dict(d=d, Y=Y, mR=Rs, msg=sgm, mfin=fin), "MAIN", PROGS["MAIN"])
        else:
            Y = main_trace(U, d)[0]
    bad |= (Y < (1 << 26)) | (Y >= (1 << 27))
    Y = np.clip(Y, 1 << 26, (1 << 27) - 1)
    fire, ang, fixed = ring_fast(Y, k, RING)
    tot = fixed.copy()
    for a_, f_ in zip(ang, fire):
        tot += np.where(f_ == 1, np.uint64(a_), np.uint64(0))
    err = (tot - tgt).view(np.int64).astype(np.float64)
    return err, fire, bad


def apply_prog(prog, H, vtop, nv, nH, j):
    """an XOR-AND program (inputs: nv v bits, then nH window bits; -1 = the constant 1) on all inputs -> word bits."""
    ins = [((vtop >> i) & 1) for i in range(nv)] + [((H >> i) & 1) for i in range(nH)]
    vals = list(ins)
    for g in prog["gates"]:
        a = np.zeros(len(H), dtype=np.int64)
        for i in g["left"]:
            a ^= 1 if i == -1 else vals[i]
        b = np.zeros(len(H), dtype=np.int64)
        for i in g["right"]:
            b ^= 1 if i == -1 else vals[i]
        vals.append(a & b)
    wnew = np.zeros(len(H), dtype=np.int64)
    for lst in prog["outputs"]:
        o = np.zeros(len(H), dtype=np.int64)
        for i in lst:
            o ^= 1 if i == -1 else vals[i]
        wnew = (wnew << 1) | o
    return wnew


def apply_to(B, which, prog):
    """the word of the tail `which` with the program (its "meta": j, extra, nH) in place of today's last j rows."""
    mt = prog["meta"]
    j, extra, nH = int(mt["j"]), int(mt.get("extra", 0)), int(mt["nH"])
    H, vtop, nv, wt, shift, kkmin = classes(B, which, j, extra)
    wnew = apply_prog(prog, H, vtop, nv, nH, j)
    return B[KEY[which]] - wvals(which, wt, j) + wvals(which, wnew, j)


def base(wd, hpath=None, progs=None):
    """progs: {tail: program} applied to the OTHER tails before the target tail is opened (joint fix rounds)."""
    progs = progs or {}
    z = np.load(os.path.join(wd, "chain.npz"))
    tgt, tol = z["target"].astype(np.uint64), z["tol"].astype(np.float64)
    x = np.arange(1, RM.MAXX + 1, dtype=np.int64)
    n, k, Q, d, delta = F8.ref_chain(x, M_)
    D = delta.astype(np.int64) << 5
    sg1, q1, rr1 = div_trace(D, d, A_TOP, S1)
    B = dict(tgt=tgt, tol=tol, k=k, d=d, q1=q1, sg1=sg1, rr1=rr1)
    if "DIV1" in progs:
        B["q1"] = q1 = apply_to(B, "DIV1", progs["DIV1"])
    sg2, P, rr2 = div2_true(q1, d)
    B.update(P=P, sg2=sg2, rr2=rr2)
    if "DIV2" in progs:
        B["P"] = P = apply_to(B, "DIV2", progs["DIV2"])
    if hpath:
        h, seen = htab(hpath)
        y = -((P + 1) >> 1)
        ycirc = np.load(hpath)["y"].astype(np.int64)
        B["y_differs_from_circuit"] = int((y != ycirc).sum())
        T = y >> 9
        B["T_outside_table"] = int(((T + 1024 < 0) | (T + 1024 > 2047) | ~seen[np.clip(T + 1024, 0, 2047)]).sum())
        U = y + h[np.clip(T + 1024, 0, 2047)]
        Y, Rs, sgm, fin = main_trace(U, d)
        B["Y_differs_from_R37_model"] = int((Y != RM.main_rows(U, d, M_, BT, BLAST, SLOT)).sum())
        B.update(Y=Y, mR=Rs, msg=sgm, mfin=fin)
        if "MAIN" in progs:
            B["Y"] = apply_to(B, "MAIN", progs["MAIN"])
    return B


def classes(B, which, j, extra):
    """class code per input for the last j rows of the tail `which`; returns (H, vtop, nH, nv, wtoday, word0)."""
    if which == "MAIN":
        kkmin = max(0, min(MT[b] for b in range(1, j + 1)) - extra)
        H = B["mR"][j] >> 8
        V = B["d"] >> 1
        vtop = (V >> kkmin) & ((1 << (7 - kkmin)) - 1)
        wt = np.zeros(len(H), dtype=np.int64)
        for b in range(j, 0, -1):
            wt = (wt << 1) | B["msg"][b]
        wt = (wt << 1) | B["mfin"]
        return H, vtop, 7 - kkmin, wt, 8, kkmin
    sched = S1 if which == "DIV1" else S2
    rrs = B["rr1"] if which == "DIV1" else B["rr2"]
    sg = B["sg1"] if which == "DIV1" else B["sg2"]
    k0 = 9 - j
    shift = min(sched[k - 1] + 1 - (k - k0) for k in range(k0, 9))
    kkmin = max(0, min(sched[k - 1] for k in range(k0, 9)) - extra)
    H = rrs[k0 - 1] >> shift
    V = B["d"] >> 1
    vtop = (V >> kkmin) & ((1 << (7 - kkmin)) - 1)
    n = len(sg)
    wt = np.zeros(len(H), dtype=np.int64)          # today's word as bits: bit i = sign bit g of digit n-1-j+i ... MSB first
    for i in range(j + 1):
        wt = (wt << 1) | sg[n - 1 - j + i]
    return H, vtop, 7 - kkmin, wt, shift, kkmin


def word_value(wbits, j):
    """digit word (bits g, MSB = the first of the j + 1 digits) -> its value sum (1 - 2 g) 2^(j - i)."""
    v = 0
    for i in range(j + 1):
        g = (wbits >> (j - i)) & 1
        v += (1 - 2 * g) << (j - i)
    return v


def one_word(args):
    which, j, w, wd, hpath, dv = args
    B = _B[0]
    wt = _B[1]
    valt = wvals(which, w, j)
    vt = wvals(which, wt, j)
    key = KEY[which]
    word = B[key] - vt + valt
    err, fire, bad = downstream(key, word, B["d"], B["k"], hpath, B["tgt"])
    sl = slack_of(err, B["tol"], fire, dv)
    sl[bad] = -1e18
    cid = _B[2]
    ms = np.full(_B[3], np.inf)
    np.minimum.at(ms, cid, sl)
    return w, ms


_B = [None, None, None, None]


def main():
    wd, hpath, outp, mode = (os.path.expanduser(sys.argv[i]) for i in (2, 3, 4, 5))
    t0 = time.time()
    if mode == "sweep":
        base_m = None
        with open(outp, "a") as f:
            for spec in json.loads(sys.argv[6]):
                o = sweep_one(spec, wd, hpath)
                if base_m is None and o.get("margin"):
                    base_m = o["margin"]                       # the first spec is the reference (k9f)
                if o.get("margin") and o["margin"] > 0:
                    o["ratio"] = round(o["margin"] / base_m, 4)
                    o["rot_toffoli_eq"] = round(-4.62 * float(np.log2(o["margin"] / base_m)), 2)
                o["t"] = round(time.time() - t0, 1)
                f.write(json.dumps(o) + "\n")
                f.flush()
                print(json.dumps(o), flush=True)
        return
    B = base(wd, hpath)
    err0, fire0, bad0 = downstream("Y", B["Y"], B["d"], B["k"], hpath, B["tgt"])
    assert not bad0.any()
    # the vectorised ring against the designer's, once, on a slice
    mg, conv, bound, dv = lp_dv(err0, B["tol"], fire0)
    print("base: LP margin %.0f conv %s bound %.0f; inputs whose y differs from the circuit's y: %d (%.0f s)"
          % (mg, conv, bound, B.get("y_differs_from_circuit", -1), time.time() - t0), flush=True)
    if mode == "table":
        which, j = sys.argv[6], int(sys.argv[7])
        nproc = int(sys.argv[8]) if len(sys.argv) > 8 and "=" not in sys.argv[8] else 8
        opt = dict(a.split("=", 1) for a in sys.argv[6:] if "=" in a)
        extra = int(opt.get("extra", 0))
        H, vtop, nv, wt, shift, kkmin = classes(B, which, j, extra)
        hmin, hmax = int(H.min()), int(H.max())
        nH = max(int(hmax).bit_length(), int(-hmin - 1).bit_length()) + 1
        code = ((H - hmin) << nv) | vtop
        ucode, cid = np.unique(code, return_inverse=True)
        ncls = len(ucode)
        # today's word must be a function of the class
        wmin = np.full(ncls, 1 << 30)
        wmax = np.full(ncls, -1)
        np.minimum.at(wmin, cid, wt)
        np.maximum.at(wmax, cid, wt)
        func = bool((wmin == wmax).all())
        cnt = np.bincount(cid, minlength=ncls)
        print("%s last %d rows: shift %d, v bits %d..6, H in [%d, %d] (%d bits), classes %d of %d patterns, today's "
              "word a function of the class: %s" % (which, j, shift, kkmin, hmin, hmax, nH, ncls, 1 << (nH + nv), func),
              flush=True)
        _B[0], _B[1], _B[2], _B[3] = B, wt, cid, ncls
        import multiprocessing as mp
        jobs = [(which, j, w, wd, hpath, dv) for w in range(1 << (j + 1))]
        ms = np.zeros((ncls, 1 << (j + 1)))
        with mp.get_context("fork").Pool(nproc) as pool:
            for w, col in pool.imap_unordered(one_word, jobs):
                ms[:, w] = col
                print("  word %d (value %+d): classes allowed at M0 = margin %d, margin/1.35 %d, 0 %d of %d (%.0f s)"
                      % (w, wvals(which, w, j), int((col >= mg - 1).sum()), int((col >= mg / 1.35).sum()),
                         int((col >= 0).sum()), ncls, time.time() - t0), flush=True)
        today = ms[np.arange(ncls), wmin]
        print("today's word: min over classes %.0f (= the margin %.0f)" % (today.min(), mg))
        np.savez_compressed(outp, ms=ms, Hval=(ucode >> nv) + hmin, vtop=ucode & ((1 << nv) - 1), nH=nH, nv=nv, j=j,
                            wtoday=wmin, cnt=cnt, margin=mg, dv=dv, which=which, shift=shift, kkmin=kkmin, func=func)
        for M0, nm in ((mg - 1, "margin"), (mg / 1.35, "margin/1.35"), (mg / 2, "margin/2"), (0.0, "0")):
            al = (ms >= M0).sum(axis=1)
            print("M0 = %s: allowed words per class: mean %.2f, classes with 1 / 2 / 3+ allowed words: %d / %d / %d"
                  % (nm, al.mean(), int((al == 1).sum()), int((al == 2).sum()), int((al >= 3).sum())))
        print("done %.0f s" % (time.time() - t0))
    elif mode == "eval":
        tab = np.load(sys.argv[6], allow_pickle=True)
        prog = json.load(open(sys.argv[7]))
        which, j = str(tab["which"]), int(tab["j"])
        opt = dict(a.split("=", 1) for a in sys.argv[8:] if "=" in a)
        H, vtop, nv, wt, shift, kkmin = classes(B, which, j, int(opt.get("extra", 0)))
        nH = int(tab["nH"])
        assert nv == int(tab["nv"])
        ins = [((vtop >> i) & 1) for i in range(nv)] + [((H >> i) & 1) for i in range(nH)]
        vals = list(ins)
        for g in prog["gates"]:
            a = np.zeros(len(H), dtype=np.int64)
            for i in g["left"]:
                a ^= 1 if i == -1 else vals[i]
            b = np.zeros(len(H), dtype=np.int64)
            for i in g["right"]:
                b ^= 1 if i == -1 else vals[i]
            vals.append(a & b)
        wnew = np.zeros(len(H), dtype=np.int64)
        for lst in prog["outputs"]:
            o = np.zeros(len(H), dtype=np.int64)
            for i in lst:
                o ^= 1 if i == -1 else vals[i]
            wnew = (wnew << 1) | o
        vt = wvals(which, wt, j)
        vn = wvals(which, wnew, j)
        key = KEY[which]
        word = B[key] - vt + vn
        err, fire, bad = downstream(key, word, B["d"], B["k"], hpath, B["tgt"])
        out = dict(prog=sys.argv[7], changed=int((vn != vt).sum()), bad=int(bad.sum()), base_margin=mg)
        sl = slack_of(err, B["tol"], fire, dv)
        out["margin_todays_angles"] = float(sl.min())
        if not bad.any():
            m2, c2, b2, _ = lp_dv(err, B["tol"], fire)
            out.update(margin=m2, conv=c2, lp_bound=b2)
        # injectivity of the whole 20-digit quotient given v (the leftover wire's clear needs it)
        full = word.astype(np.int64) + (1 << 21)
        codeq = (B["d"].astype(np.int64) << 23) | full
        coder = (B["d"].astype(np.int64) << 23) | (B[key] + (1 << 21))
        out["distinct_q_today"] = int(len(np.unique(coder)))
        out["distinct_q_new"] = int(len(np.unique(codeq)))
        json.dump(out, open(outp, "w"), indent=1)
        print(json.dumps(out))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
