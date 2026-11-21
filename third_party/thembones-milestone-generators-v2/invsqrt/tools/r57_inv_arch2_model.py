#!/usr/bin/env python3
"""r57_inv_arch2_model.py -- R57-INV-ARCH2 (designer): the invsqrt chain of family C as ONE exact integer model with
every length a parameter, on ALL 3,133,452 inputs, angle LP re-solved.  My own file; read-only imports of R54's exact
model of k9f (tools/r54_inv_tail_dc.py: reproduces k9f's generator, margin 51,917,278) and through it R37's model and
MY copy of FC's source (argv[1]).

  r57_inv_arch2_model.py SRC WD HTAB.npz OUT cal                  the base: k9f's chain, k9f's stage table -> 51,917,278
  r57_inv_arch2_model.py SRC WD HTAB.npz OUT terms                the table by function: which block moves which digit
  r57_inv_arch2_model.py SRC WD HTAB.npz OUT budget               where the tolerance goes: readout alone, + MAIN's
                                                                  tail, + the dividers' tails, + the stage's table
  r57_inv_arch2_model.py SRC WD HTAB.npz OUT sweep SPECS.json [NPROC]   any chain: spec keys
        m (7), F1 (5), a (exact rows - 1 of DIV1 / DIV2; default R37's), sched1, sched2 (tail tapers),
        ho: "k9f" (the circuit's table) | "ideal" (the whole series of the y that arrives) | "trunc" (ideal table of
        T = y >> s) | "quad" | "pade"; s, off; mt (MAIN's tapers {row: taper}; null = the circuit's at m 7), bt, blast (1),
        ring (26 = slots bits 1..26), name
Safe to run: reads only; writes OUT under my own folder."""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r54_inv_tail_dc as DC          # noqa: E402  (read-only; imports R37's model which reads sys.argv[1] = SRC)

RM, F8, G = DC.RM, DC.F8, DC.G
ATOP = {8: 8, 7: 10, 6: 12, 5: 14}


def targets(wd):
    z = np.load(os.path.join(wd, "chain.npz"))
    return z["target"].astype(np.uint64), z["tol"].astype(np.float64)


def rot_eq(mg, nrot):
    """Toffoli-equivalents of nrot rotations at LP margin mg (the verifier's law)."""
    if mg is None or mg <= 0:
        return None
    eps = math.pi * mg / 2.0 ** 64 / nrot
    return nrot * (0.56 * math.log2(1.0 / eps) + 5.3) / 4.0


def lp_word(Y, k, tgt, tol, N=26, extra=None):
    fire, ang, fixed = DC.ring_fast(Y, k, N)
    tot = fixed.copy()
    for a_, f_ in zip(ang, fire):
        tot += np.where(f_ == 1, np.uint64(a_), np.uint64(0))
    err = (tot - tgt).view(np.int64).astype(np.float64)
    if extra is not None and len(extra):
        fire = np.vstack([fire, np.array(extra, dtype=np.uint8)])
    mg, conv, bound, dv = DC.lp_dv(err, tol, fire)
    return mg, conv, bound, fire.shape[0]


def div2(q1, d, a, sched, m):
    sg, P, rrs = DC.div_trace(q1, d, a, sched)
    if m == 7 and list(sched[-2:]) == [6, 6] and len(sched) == 8:      # the circuit's 2-wire last row
        M8 = rrs[7] >> 7
        A = 1 + ((d >> 7) & 1)
        res = np.where(M8 < 0, (M8 & 3) + A, (M8 & 3) - A)
        fin_t = (res >> 1) & 1
        P = P + 2 * (sg[-1] - fin_t)
    return P


def main_word(U, d, spec):
    m = int(spec.get("m", 7))
    blast = int(spec.get("blast", 1))
    mt = spec.get("mt")
    if m == 7 and blast == 1 and (mt is None or mt == "k9f"):
        return DC.main_trace(U, d)[0]
    if m == 7 and blast == 1 and isinstance(mt, dict):
        keep = dict(DC.MT)
        DC.MT.clear()
        DC.MT.update({int(b): int(t) for b, t in mt.items()})
        try:
            if DC.MT.get(1) == 6:
                return DC.main_trace(U, d)[0]
            # generic rows with a per-row taper
            N = (np.int64(1) << (28 + m)) + (U.astype(np.int64) << (m + 1))
            R = N.copy()
            Qw = np.zeros(len(d), dtype=np.int64)
            for b in range(26, blast - 1, -1):
                t = DC.MT.get(b, 0)
                dt = F8.tail_dt(d, t) if t else d
                s_ = np.where(R >= 0, 1, -1)
                R = R - s_ * (dt << b)
                Qw = Qw + s_ * (1 << b)
            Y = Qw - ((R < 0).astype(np.int64) << blast)
            return Y - (Y & ((1 << blast) - 1))
        finally:
            DC.MT.clear()
            DC.MT.update(keep)
    bt = int(spec.get("bt", m - 1))
    return RM.main_rows(U, d, m, bt, blast, spec.get("slot", "zero"))


_X = {}


def chain_y(spec):
    m = int(spec.get("m", 7))
    F1 = int(spec.get("F1", 5))
    if "x" not in _X:
        _X["x"] = np.arange(1, RM.MAXX + 1, dtype=np.int64)
    if m not in _X:
        _X[m] = F8.ref_chain(_X["x"], m)
    n, k, Q, d, delta = _X[m]
    F1d = int(spec.get("F1d", F1))                  # DIV1's dividend shift (digits of q1); F1 is the stage's scale
    qs = int(spec.get("q1shift", 0))                # zeros DIV2 brings down under q1 (q1 one digit shorter each)
    a = spec.get("a")
    if a is None:
        a = ATOP[m] - (5 - F1d)
    a2 = spec.get("a2")
    if a2 is None:
        a2 = a + qs
    s1 = spec.get("sched1") or (DC.S1 if m == 7 else [0, 0] + list(range(1, m)))
    s2 = spec.get("sched2") or (DC.S2 if m == 7 else [0, 0] + list(range(1, m)))
    D = delta.astype(np.int64) << F1d
    _, q1, _ = DC.div_trace(D, d, a, s1)
    q1 = (q1 << qs) + int(spec.get("q1add", 0))
    P = div2(q1, d, a2, s2, m)
    y = -((P + 1) >> 1)
    return n, k, d, delta, y


def run_spec(spec, tgt, tol, hpath):
    t0 = time.time()
    out = dict(spec)
    m = int(spec.get("m", 7))
    n, k, d, delta, y = chain_y(spec)
    ho = spec.get("ho", "ideal")
    if ho == "k9f":
        h, seen = DC.htab(hpath)
        T = y >> 9
        if (T + 1024).min() < 0 or (T + 1024).max() > 2047 or not seen[T + 1024].all():
            out.update(margin=None, note="T leaves the stage's table")
            return out
        U = y + h[T + 1024]
    elif ho == "yex":                               # y itself exact (no divider error), the whole series
        eps = delta.astype(np.float64) / (4.0 ** (9 - m) * d.astype(np.float64) ** 2)
        U = np.round(2.0 ** 27 * ((1.0 + eps) ** -0.5 - 1.0)).astype(np.int64) + int(spec.get("off", 0))
    else:
        sp = dict(spec)
        sp["ho"] = ho
        sp.setdefault("off", -1)
        U = RM.ho(y, sp)
    cr = spec.get("cross")
    if cr:                                          # the T x L cross term the stage drops (h is read at L = 256)
        T_ = y >> 9
        L_ = y & 511
        yt = (T_.astype(np.float64) + 0.5) * 512.0
        if cr.get("mode") == "thr":                 # +-1 unit: (2 L8 - 1) sgn(y) [|T + 1/2| >= thr]
            c = np.where(np.abs(T_ + 0.5) >= float(cr["thr"]), 1, 0) * np.where(L_ >= 256, 1, -1) * np.where(yt < 0, -1, 1)
        else:
            nb = int(cr.get("bits", 1))             # top nb bits of L
            Lq = ((L_ >> (9 - nb)) << (9 - nb)) + (1 << (8 - nb)) if nb < 9 else L_
            c = np.round(6.0 * yt * (Lq - 256.0) / 2.0 ** 28).astype(np.int64)
        U = U + c.astype(np.int64)
        out["cross_abs_max"] = int(np.abs(c).max())
        out["cross_nonzero"] = float((c != 0).mean())
    Y = main_word(U, d, spec)
    out.update(y_min=int(y.min()), y_max=int(y.max()), U_min=int(U.min()), U_max=int(U.max()))
    if Y.min() < (1 << 26) or Y.max() >= (1 << 27):
        out.update(margin=None, note="word leaves [2^26, 2^27)", Y_min=int(Y.min()), Y_max=int(Y.max()))
        return out
    N = int(spec.get("ring", 26))
    extra = None
    if spec.get("feat") == "v":                     # (3): the 7 operand wires as extra rotation controls
        extra = [(((d >> 1) >> i) & 1).astype(np.uint8) for i in range(7)]
    elif spec.get("feat") == "vk":                  # the same, one angle per shift class (an upper bound: 77 controls)
        extra = [((((d >> 1) >> i) & 1) & (k == kk)).astype(np.uint8) for i in range(7) for kk in range(11)]
    mg, conv, bound, nrot = lp_word(Y, k, tgt, tol, N, extra)
    out.update(margin=mg, conv=bool(conv), lp_bound=bound, n_rot=nrot, rot_eq=rot_eq(mg, nrot),
               t=round(time.time() - t0, 1))
    return out


def worker(args):
    spec, wd, hpath = args
    tgt, tol = targets(wd)
    try:
        return run_spec(spec, tgt, tol, hpath)
    except Exception as e:                      # noqa: BLE001
        o = dict(spec)
        o["error"] = repr(e)
        return o


def mode_terms(wd, hpath, outp):
    """The table by function: which block moves which digit of the word, and by how much."""
    tgt, tol = targets(wd)
    n, k, d, delta, y = chain_y({})
    h, seen = DC.htab(hpath)
    T = y >> 9
    U = y + h[T + 1024]
    zero = np.zeros_like(y)
    words = dict(const=DC.main_trace(zero, d)[0],            # 2^35 / d alone
                 first=DC.main_trace(y, d)[0],               # + first order (DIV1, DIV2)
                 full=DC.main_trace(U, d)[0])                # + the stage
    quad = y + np.round(3.0 * ((T.astype(np.float64) + 0.5) * 512.0) ** 2 / 2.0 ** 28).astype(np.int64)
    words["quad"] = DC.main_trace(quad, d)[0]
    res = dict(inputs=int(len(y)))
    res["first_order_max_units_of_Y"] = int(np.abs(words["first"] - words["const"]).max())
    res["stage_max_units_of_Y"] = int(np.abs(words["full"] - words["first"]).max())
    res["cubic_max_units_of_Y"] = int(np.abs(words["full"] - words["quad"]).max())
    res["h_max_units_of_U"] = int(h[T + 1024].max())
    res["y_abs_max"] = int(np.abs(y).max())
    rows = []
    for b in range(26, 0, -1):
        f = lambda A, B: float((((A >> b) & 1) != ((B >> b) & 1)).mean())    # noqa: E731
        rows.append(dict(digit=b, moved_by_first_order=round(f(words["const"], words["first"]), 4),
                         moved_by_stage=round(f(words["first"], words["full"]), 4),
                         moved_by_cubic=round(f(words["quad"], words["full"]), 4)))
    res["digits"] = rows
    # tolerance in units of Y per class of k
    Yf = words["full"].astype(np.float64)
    res["tolerance_units_of_Y_min_max"] = [float((Yf * 2.0 ** -24).min()), float((Yf * 2.0 ** -24).max())]
    json.dump(res, open(outp, "w"), indent=1)
    print(json.dumps({kk: vv for kk, vv in res.items() if kk != "digits"}))
    for r in rows:
        print(r)


def mode_budget(wd, hpath, outp):
    """Where the tolerance goes.  Each line is the LP margin (26 slots + 7 features, all inputs) of a chain that is
    exact up to the named point."""
    tgt, tol = targets(wd)
    n, k, d, delta, y = chain_y({})
    res = []

    def put(name, Y):
        mg, conv, bound, nrot = lp_word(Y, k, tgt, tol, 26)
        r = dict(name=name, margin=mg, conv=bool(conv), rot_eq=rot_eq(mg, nrot))
        res.append(r)
        print(json.dumps(r), flush=True)

    print(json.dumps(dict(tol_min=float(tol.min()), tol_max=float(tol.max()))), flush=True)
    Yf = 2.0 ** 37 / np.sqrt(n.astype(np.float64))
    for nm, Yq in (("readout alone: the true word rounded to the nearest even", 2 * np.round(Yf / 2.0).astype(np.int64)),
                   ("readout alone: the true word floored to even", 2 * np.floor(Yf / 2.0).astype(np.int64)),
                   ("readout alone: the true word rounded to odd/2 centres (floor to even, +1 in the angles)",
                    2 * np.floor((Yf + 1.0) / 2.0).astype(np.int64))):
        Yq = np.clip(Yq, 1 << 26, (1 << 27) - 2)
        put(nm, Yq)
    eps = delta.astype(np.float64) / (16.0 * d.astype(np.float64) ** 2)
    for off in (0, -1, 1):
        Uex = np.round(2.0 ** 27 * ((1.0 + eps) ** -0.5 - 1.0)).astype(np.int64) + off
        put("U exact (rounded, off %d), MAIN with exact rows down to b = 1 (no taper)" % off,
            main_word(Uex, d, dict(mt={})))
        put("U exact (off %d), MAIN as the circuit (tapers 1, 2, 3, 4, 6)" % off, DC.main_trace(Uex, d)[0])
    for off in (-1, 0, 1):
        Ui = RM.ho(y, dict(ho="ideal", off=off))
        put("y as the circuit's dividers make it, stage = the whole series of that y (off %d), MAIN as the circuit" % off,
            DC.main_trace(Ui, d)[0])
    for s in (9, 8, 7):
        Ut = RM.ho(y, dict(ho="trunc", s=s, off=-1))
        put("the same with the stage an ideal table of T = y >> %d" % s, DC.main_trace(Ut, d)[0])
    h, seen = DC.htab(hpath)
    T = y >> 9
    put("k9f itself (the circuit's table)", DC.main_trace(y + h[T + 1024], d)[0])
    json.dump(res, open(outp, "w"), indent=1)


def main():
    wd, hpath, outp, mode = (os.path.expanduser(sys.argv[i]) for i in (2, 3, 4, 5))
    t0 = time.time()
    if mode == "cal":
        tgt, tol = targets(wd)
        o = run_spec(dict(name="k9f", ho="k9f"), tgt, tol, hpath)
        print(json.dumps(o))
        json.dump(o, open(outp, "w"))
    elif mode == "terms":
        mode_terms(wd, hpath, outp)
    elif mode == "budget":
        mode_budget(wd, hpath, outp)
    elif mode == "sweep":
        specs = json.load(open(sys.argv[6]))
        nproc = int(sys.argv[7]) if len(sys.argv) > 7 else 8
        import multiprocessing as mp
        with mp.get_context("fork").Pool(nproc) as pool, open(outp, "a") as f:
            for o in pool.imap(worker, [(s, wd, hpath) for s in specs]):
                f.write(json.dumps(o) + "\n")
                f.flush()
                print(json.dumps(o), flush=True)
    else:
        raise SystemExit(__doc__)
    print("done %.0f s" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
