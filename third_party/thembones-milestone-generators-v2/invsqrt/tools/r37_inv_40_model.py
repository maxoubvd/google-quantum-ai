#!/usr/bin/env python3
"""r37_inv_40_model.py -- R37 designer INV-40 (own code; imports a COPY of FC's frozen src read-only).

Exact integer model of family C at any root split m (8 = today's 41-wire circuit, 7 = 40 wires, 6 = 39 wires in FC's
layout: v m + dividend 27 + flags 4 + held carry 1 + 1 clean), on ALL 3,133,452 inputs, with a CONVERGED angle LP
(S1's lp_fast, stop_below = -inf, rounds 400; the printed margin is the all-input margin of the returned offsets and
`conv` says whether the cutting-plane loop met its own stop rule).

Chain (FC's references, plain DIV1 = cfg nocorr): n, k, Q, d = 2Q + 1, delta = n - 4^(9-m) d^2;
theta = G.ref_div(delta << F1, d, a, sched1); P = G.ref_div(theta, d, a, sched2); y = -((P + 1) >> 1);
U = HO(y); Y = F8.ref_main(U, d, m, bt); endpoint = F1's 27-ring (FC2's lp_margin copy: 27 slots + 7 tunable wrap
features, 3 fixed).

HO variants (spec["ho"]):
  "ideal"            U = round(f(y)) + off,   f(y) = 2^27 ((1 - y / 2^26)^(-1/2) - 1)   (the whole series)
  "pade"             U = round(2^28 y / (2^28 - 3 y)) + off                              (FC2's H1_exact)
  "trunc"            U = y + round(g(yt)) + off, T = y >> s, yt = (T + 1/2) 2^s, g = f(yt) - yt   (ideal table of T)
  "quad"             U = y + round(3 yt^2 / 2^28 + cub) + off, cub = 0 | "c3": 5 yt^3 / 2^55 exact | "st<q>": the
                     cubic rounded to steps of q units
  spec keys: m, s, bt, off, F1 (5), a (None = FC's a_top), sched1, sched2, name

  r37_inv_40_model.py SRC WD OUT.jsonl NPROC SPECS.json      (SPECS = a JSON list of spec dicts)
"""
import json
import math
import os
import sys
import time

import numpy as np

SRC = os.path.abspath(os.path.expanduser(sys.argv[1]))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, SRC)
import invsqrt_fc_g8 as F8   # noqa: E402
import invsqrt_fc_gen as G   # noqa: E402
import invsqrt_sw_model as SM  # noqa: E402

MAXX = F8.MAX_X


def ring_rows(Y, k):
    """FC2's ring model (invsqrt_fc2_h1lp.lp_margin, copied): fire rows, base angles, fixed phase."""
    N = 27
    w = [((Y >> i) & 1).astype(np.uint8) for i in range(N)]
    M = 1 << 64
    fire, ang = [], []
    for p in range(N):
        i = (p - k) % N
        bits = np.zeros(len(Y), dtype=np.uint8)
        for ii in range(N):
            sel = i == ii
            bits[sel] = w[ii][sel]
        fire.append(bits)
        ang.append(1 << (p + 27))
    fixed = np.zeros(len(Y), dtype=np.uint64)
    for j in range(10):
        on = (k >= j + 1)
        i = np.where(on, j - k + N, 0)
        bits = np.zeros(len(Y), dtype=np.uint8)
        for ii in range(N):
            sel = on & (i == ii)
            bits[sel] = w[ii][sel]
        if j <= 6:
            fire.append(bits)
            ang.append(((1 << (j + 54)) - (1 << (j + 27))) % M)
        else:
            fixed += np.where(bits == 1, np.uint64((1 << (54 + j - 7 + 7)) % M), np.uint64(0))
    return np.array(fire, dtype=np.uint8), ang, fixed


def ring_rows_n(Y, k, N):
    """a ring of N < 27 slots for the word bits sh .. 26 (sh = 27 - N; the low sh bits of Y are not slots): token j
    = word bit j + sh at slot (j + k) mod N, slot angle 2^(p + sh - 37) turns; a token with j + k >= N has wrapped to
    slot q = j + k - N and needs 2^(q - 10) turns: wrap features q = 0..9 as in the 27-ring (0..6 tunable; 7, 8, 9
    fixed at T, S, Z). Rotations: N + 7."""
    sh = 27 - N
    w = [((Y >> (i + sh)) & 1).astype(np.uint8) for i in range(N)]
    M = 1 << 64
    fire, ang = [], []
    for p in range(N):
        i = (p - k) % N
        bits = np.zeros(len(Y), dtype=np.uint8)
        for ii in range(N):
            sel = i == ii
            bits[sel] = w[ii][sel]
        fire.append(bits)
        ang.append(1 << (p + 27 + sh))
    fixed = np.zeros(len(Y), dtype=np.uint64)
    for j in range(10):
        on = (k >= j + 1)
        i = np.where(on, j - k + N, 0)
        bits = np.zeros(len(Y), dtype=np.uint8)
        for ii in range(N):
            sel = on & (i == ii)
            bits[sel] = w[ii][sel]
        if j <= 6:
            fire.append(bits)
            ang.append(((1 << (j + 54)) - (1 << (j + 27 + sh))) % M)
        else:
            fixed += np.where(bits == 1, np.uint64((1 << (54 + j)) % M), np.uint64(0))
    return np.array(fire, dtype=np.uint8), ang, fixed


def lp_conv(err, tol, fire, rounds=400):
    """lp_fast run to its own stop rule; returns (all-input margin, converged?, LP bound of the last round)."""
    from scipy.optimize import linprog
    K = fire.shape[0]
    m0 = tol - np.abs(err)
    start = add = 4000
    act = np.unique(np.argpartition(m0, start)[:start])
    best = float(m0.min())
    bound = float("inf")
    conv = False
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
        marg = np.empty(len(err), dtype=np.float64)
        if SM.HAVE_NUMBA:
            SM._margin_kernel(err, tol, fire, dv.astype(np.float64), marg)
        else:
            e = err.copy()
            for r in np.flatnonzero(dv):
                e += dv[r] * fire[r]
            marg = tol - np.abs(e)
        mmin = float(marg.min())
        best = max(best, mmin)
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
    return best, conv, bound


def f_ideal(yf):
    return 2.0 ** 27 * ((1.0 - yf / 2.0 ** 26) ** -0.5 - 1.0)


def chain(m, F1=5, a=None, sched1=None, sched2=None):
    x = np.arange(1, MAXX + 1, dtype=np.int64)
    n, k, Q, d, delta = F8.ref_chain(x, m)
    if a is None:
        a = {8: 8, 7: 10, 6: 12, 5: 14}[m] - (5 - F1)
    s1 = sched1 or ([0, 0] + list(range(1, m)))
    s2 = sched2 or ([0, 0] + list(range(1, m)))
    D = delta.astype(np.int64) << F1
    _, q1 = G.ref_div(D, d, a, s1)
    _, P = G.ref_div(q1, d, a, s2)
    y = -((P + 1) >> 1)
    return k, d, y, delta


def ho(y, spec):
    g = 5 - int(spec.get("F1", 5))
    if g:
        # U one (or g) bit(s) shorter: y is in units of 2^g; U = 2^g U_g + b9 (b9 = the constant low dividend bits)
        sc = float(1 << g)
        yy = sc * y.astype(np.float64) + float(spec.get("beta", (sc - 1) / 2))
        b9 = int(spec.get("b9", 0))
        kind = spec["ho"]
        if kind == "ideal":
            tgt = f_ideal(yy)
        elif kind == "trunc":
            s_ = int(spec["s"])
            yt = ((y >> s_).astype(np.float64) + 0.5) * 2.0 ** s_ * sc
            tgt = yy + (f_ideal(yt) - yt)
        else:
            raise SystemExit("ho kind not modelled for F1 < 5: " + kind)
        Ug = np.round((tgt + float(spec.get("off", 0)) - b9) / sc).astype(np.int64)
        return (Ug << g) + b9
    kind = spec["ho"]
    off = int(spec.get("off", -1))
    yf = y.astype(np.float64)
    if kind == "r37":                              # my gate-form stage (tools/r37_inv_40_ref.py), spec["st"]
        import r37_inv_40_ref as R37
        st = dict(R37.DEF7)
        st.update(spec.get("st", {}))
        if st.get("cub") == "exact":               # (model only) the cubic as the exact rounded value
            st2 = dict(st)
            st2["cub"] = []
            U, _ = R37.stage_ref(y, st2)
            s_ = int(spec.get("s", 9))
            yt = ((y >> s_).astype(np.float64) + 0.5) * 2.0 ** s_
            return U + np.round(5 * yt ** 3 / 2.0 ** 55).astype(np.int64) + off
        return R37.stage_ref(y, st)[0] + off
    if kind == "ideal":
        return np.round(f_ideal(yf)).astype(np.int64) + off
    if kind == "pade":
        return np.round(2.0 ** 28 * yf / (2.0 ** 28 - 3 * yf)).astype(np.int64) + off
    s = int(spec["s"])
    yt = ((y >> s).astype(np.float64) + 0.5) * 2.0 ** s
    if kind == "trunc":
        return y + np.round(f_ideal(yt) - yt).astype(np.int64) + off
    if kind == "quad":
        h = 3 * yt * yt / 2.0 ** 28
        cub = spec.get("cub", 0)
        c3 = 5 * yt ** 3 / 2.0 ** 55
        if cub == "c3":
            h = h + c3
        elif isinstance(cub, str) and cub.startswith("st"):
            q = float(cub[2:])
            h = h + np.round(c3 / q) * q
        return y + np.round(h).astype(np.int64) + off
    raise SystemExit("unknown ho " + kind)


def tper(mg, n=34):
    return 0.56 * math.log2(n * 2 ** 64 / (math.pi * mg)) + 5.3 if mg > 0 else float("nan")


def main_rows(U, d, m, bt, blast=0, slot="std"):
    """FC's ref_main with rows b = 26 .. blast; the final sign has weight 2^blast; slots below blast: `slot`
    'zero' | 'one' | 'rest' (the window bit under the sign: bit 9 + blast of the final remainder) | 'nrest'."""
    if blast == 0 and slot == "std":
        return F8.ref_main(U, d, m, bt)
    N = (np.int64(1) << (28 + m)) + (U.astype(np.int64) << (m + 1))
    R = N.copy()
    Qw = np.zeros(len(d), dtype=np.int64)
    for b in range(26, blast - 1, -1):
        t = max(0, bt - b)
        dt = F8.tail_dt(d, t) if t else d
        s_ = np.where(R >= 0, 1, -1)
        R = R - s_ * (dt << b)
        Qw = Qw + s_ * (1 << b)
    Y = Qw - ((R < 0).astype(np.int64) << blast)
    if blast:
        Y = Y - (Y & ((1 << blast) - 1))          # Qw's forced low ones are not slots
        rest = (R >> (m + 1 + blast)) & 1
        if slot == "one":
            Y = Y + (1 << (blast - 1))
        elif slot == "rest":
            Y = Y + (rest << (blast - 1))
        elif slot == "nrest":
            Y = Y + ((1 - rest) << (blast - 1))
    return Y


_CH = {}


def run_one(spec, tgt, tol):
    t0 = time.time()
    m = int(spec["m"])
    key = (m, spec.get("F1", 5), spec.get("a"), json.dumps(spec.get("sched1")), json.dumps(spec.get("sched2")))
    if key not in _CH:
        _CH[key] = chain(m, spec.get("F1", 5), spec.get("a"), spec.get("sched1"), spec.get("sched2"))
    k, d, y, delta = _CH[key]
    U = ho(y, spec)
    Y = main_rows(U, d, m, int(spec["bt"]), int(spec.get("blast", 0)), spec.get("slot", "std"))
    out = dict(spec)
    out.update(y_min=int(y.min()), y_max=int(y.max()), U_min=int(U.min()), U_max=int(U.max()),
               Y_min=int(Y.min()), Y_max=int(Y.max()))
    if Y.min() < (1 << 26) or Y.max() >= (1 << 27):
        out.update(margin=None, note="word leaves [2^26, 2^27)", t=round(time.time() - t0, 1))
        return out
    N = int(spec.get("ring", 27))
    fire, ang, fixed = ring_rows(Y, k) if N == 27 else ring_rows_n(Y, k, N)
    tot = fixed.copy()
    for a_, f_ in zip(ang, fire):
        tot += np.where(f_ == 1, np.uint64(a_), np.uint64(0))
    err = (tot - tgt).view(np.int64).astype(np.float64)
    mg, conv, bound = lp_conv(err, tol, fire)
    nrot = fire.shape[0]
    out.update(margin=mg, conv=conv, lp_bound=bound, n_rot=nrot, tper=round(tper(mg, nrot), 3),
               rot_cost=round(nrot * tper(mg, nrot) / 4.0, 2) if mg > 0 else None, t=round(time.time() - t0, 1))
    return out


def worker(args):
    spec, wd, cores = args
    if cores:
        try:
            os.sched_setaffinity(0, cores)
        except Exception:
            pass
    z = np.load(os.path.join(wd, "chain.npz"))
    tgt, tol = z["target"].astype(np.uint64), z["tol"].astype(np.float64)
    try:
        return run_one(spec, tgt, tol)
    except Exception as e:                      # noqa: BLE001
        o = dict(spec)
        o["error"] = repr(e)
        return o


def main():
    wd, outp, nproc, specp = sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
    specs = json.load(open(specp))
    cores = sorted(os.sched_getaffinity(0))
    import multiprocessing as mp
    with mp.Pool(nproc) as pool, open(outp, "a") as f:
        jobs = [(s, wd, None) for s in specs]
        for o in pool.imap_unordered(worker, jobs):
            f.write(json.dumps(o) + "\n")
            f.flush()
            print(json.dumps(o), flush=True)
    print("done; cores", cores, flush=True)


if __name__ == "__main__":
    main()
