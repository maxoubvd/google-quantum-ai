"""invsqrt_gen: assemble a complete invsqrt circuit from the task specification, check it exactly on every
input against the numeric model and the verifier's reference, and write .kmx / .kmb (own code).

usage: invsqrt_gen_build.py <variant> <outdir> [key=value ...] [--sample N] [--no-check] [--seeds N]
  key=value pairs that name a field of invsqrt_gen_model.Params select the numeric model, the others are options
  of the variant's builder.
variants:
  g0  first complete circuit: plain shifter, unary-iteration lookup, full-width rows, 4-tier streams with
      one-hot delay wires                                                       (model: glo=2 gu=1 fc2=1)
  g1  g0 + binary delay wires with monomial fix-ups + blocked carries           (blk_a, blk_b)
  g2  g1 + pool lookups with free erasure, address kept alive                   (n_early, n_scr, seed)
  g3  public register footprint (model: glo=2 gu=1 fc2=1 virt=1 tie_c1=1), replaceable lookup blocks, lent q21,
      options erase_addr=1 (address erased after the first lookup), refl=1 (reflection shifter), keep_c2=1,
      n_e, n_e1, s1, blk_a, blk_b, seed                                         (file invsqrt_gen_g3.py)
  g4  g3 with stage A as a product in a small measured register                 (file invsqrt_gen_g4.py)
      options deep (lookup solver 2), f1, f1b, fa (free fills)
  g5  g4 on 62 wires: exponent parked on t0..t4 + one flag, deferred stream fix-ups, block lists for the rows
      (files invsqrt_gen_g5.py, invsqrt_gen_code.py)
  g6  g5 with stage A as rows on the low window of u (model: add tie_t0=1), 13 lo wires before A, table wires
      measured with the repair deferred to the unlookup                         (file invsqrt_gen_g6.py)
  the best verified build and its options are in the state summary of research/ws_invsqrt_gen.md
The build is deterministic: tables come from the fit in invsqrt_gen_model.py (cached by parameter hash), the
lookup synthesiser draws from random.Random seeded by (seed, block name), gate order is fixed.
Sources of the ideas used are credited in the headers of invsqrt_gen_model.py, invsqrt_gen_parts.py and
invsqrt_gen_pool.py; no table, angle or gate list of any public circuit enters the generator.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_gen_model as M
import invsqrt_gen_parts as S
from invsqrt_gen_circ import Circ, Sim, kmx_to_kmb


def build_g0(P, tabs, angles=None):
    c = Circ(22)
    xq = list(range(22))
    dfun = lambda e: 10 - (e >> 1)
    # ---------------- NORM
    c.stage("NORM")
    f = S.thermo(c, xq)
    oh = [S.linear_of_e(c, f, lambda e, d=d: int(dfun(e) == d)) for d in range(11)]
    p = S.linear_of_e(c, f, lambda e: e & 1)
    db = [S.linear_of_e(c, f, lambda e, b=b: (dfun(e) >> b) & 1) for b in range(4)]
    S.thermo_erase(c, xq, f)
    stages = [(16, db[3], False), (8, db[2], False), (4, db[1], False), (2, db[0], False), (1, p, True)]
    S.shifter(c, xq, stages)
    # ---------------- LUT
    c.stage("LUT")
    H, u, lo, c2 = c.allocn(P.wh), c.allocn(P.wu), c.allocn(P.wlo), c.allocn(P.wc2)
    c.meta.update(H=H, u=u, lo=lo, c2=c2, p=p, oh=oh)
    data = [tabs["H"][a] | (tabs["C1"][a] << P.wh) | (tabs["C0"][a] << (P.wh + P.wu))
            | (tabs["C2"][a] << (P.wh + P.wu + P.wlo)) for a in range(256)]
    addr = [p] + [xq[k] for k in range(20, 13, -1)]
    S.qrom(c, addr, data, H + u + lo + c2)
    # ---------------- streams share the tier description
    tiers = [(1, lambda d: d & 1), (2, lambda d: d & 2), (4, lambda d: d & 4), (8, lambda d: d & 8)]
    cells = [[xq[0]], [xq[1], xq[2]], [xq[3], xq[4], xq[5], xq[6]], [xq[k] for k in range(3, 11)]]
    delays = list(range(11))
    c.stage("PH1")
    S.stream(c, H, tiers, db, cells, delays, oh, 7 + 10, P.top_pos, angles)
    c.stage("A")
    for i in P.a_rows:
        S.row(c, u, c2, P.a_shift(i), xq[i], 1)
    c.stage("B")
    for i in P.b_rows:
        S.row(c, lo, u, P.b_shift(i), xq[i], 1)
    c.stage("PH2")
    S.stream(c, lo, tiers, db, cells, delays, oh, 25 + P.glo + 10, P.top_pos, angles)
    c.stage("B^-1")
    for i in reversed(P.b_rows):
        S.row(c, lo, u, P.b_shift(i), xq[i], 0)
    c.stage("A^-1")
    for i in reversed(P.a_rows):
        S.row(c, u, c2, P.a_shift(i), xq[i], 0)
    c.stage("unLUT")
    S.qrom(c, addr, data, H + u + lo + c2)
    for w in H + u + lo + c2:
        c.free(w)
    # ---------------- UNNORM
    c.stage("UNNORM")
    S.shifter(c, xq, stages, inverse=True)
    f = S.thermo(c, xq)
    for d in range(11):
        S.linear_of_e(c, f, lambda e, d=d: int(dfun(e) == d), oh[d])
        c.free(oh[d])
    S.linear_of_e(c, f, lambda e: e & 1, p)
    c.free(p)
    for b in range(4):
        S.linear_of_e(c, f, lambda e, b=b: (dfun(e) >> b) & 1, db[b])
        c.free(db[b])
    S.thermo_erase(c, xq, f)
    c.stage("END")
    return c


def build_g1(P, tabs, angles=None, blk_a=32, blk_b=7):
    """G0 with (1) only p and the four binary delay wires kept (stream fix-ups through monomial wires) and
    (2) blocked row adders: forward rows defer the phase of their boundary carries to the mirrored reverse row,
    reverse rows repair theirs under the measurement outcome."""
    c = Circ(22)
    xq = list(range(22))
    dfun = lambda e: 10 - (e >> 1)
    c.stage("NORM")
    f = S.thermo(c, xq)
    p = S.linear_of_e(c, f, lambda e: e & 1)
    db = [S.linear_of_e(c, f, lambda e, b=b: (dfun(e) >> b) & 1) for b in range(4)]
    S.thermo_erase(c, xq, f)
    stages = [(16, db[3], False), (8, db[2], False), (4, db[1], False), (2, db[0], False), (1, p, True)]
    S.shifter(c, xq, stages)
    c.stage("LUT")
    H, u, lo, c2 = c.allocn(P.wh), c.allocn(P.wu), c.allocn(P.wlo), c.allocn(P.wc2)
    c.meta.update(H=H, u=u, lo=lo, c2=c2, p=p)
    data = [tabs["H"][a] | (tabs["C1"][a] << P.wh) | (tabs["C0"][a] << (P.wh + P.wu))
            | (tabs["C2"][a] << (P.wh + P.wu + P.wlo)) for a in range(256)]
    addr = [p] + [xq[k] for k in range(20, 13, -1)]
    S.qrom(c, addr, data, H + u + lo + c2)
    tier_bits = [0, 1, 2, 3]
    cells = [[xq[0]], [xq[1], xq[2]], [xq[3], xq[4], xq[5], xq[6]], [xq[k] for k in range(3, 11)]]
    delays = list(range(11))
    c.stage("PH1")
    S.stream_m(c, H, tier_bits, db, cells, delays, 7 + 10, P.top_pos, angles)
    c.stage("A")
    pa, pb = {}, {}
    for i in P.a_rows:
        pa[i] = S.row_blocked(c, u, c2, P.a_shift(i), xq[i], 1, blk_a, "fwd")
    c.stage("B")
    for i in P.b_rows:
        pb[i] = S.row_blocked(c, lo, u, P.b_shift(i), xq[i], 1, blk_b, "fwd")
    c.stage("PH2")
    S.stream_m(c, lo, tier_bits, db, cells, delays, 25 + P.glo + 10, P.top_pos, angles)
    c.stage("B^-1")
    for i in reversed(P.b_rows):
        S.row_blocked(c, lo, u, P.b_shift(i), xq[i], 0, blk_b, "rev", pb[i])
    c.stage("A^-1")
    for i in reversed(P.a_rows):
        S.row_blocked(c, u, c2, P.a_shift(i), xq[i], 0, blk_a, "rev", pa[i])
    c.stage("unLUT")
    S.qrom(c, addr, data, H + u + lo + c2)
    for w in H + u + lo + c2:
        c.free(w)
    c.stage("UNNORM")
    S.shifter(c, xq, stages, inverse=True)
    f = S.thermo(c, xq)
    S.linear_of_e(c, f, lambda e: e & 1, p)
    c.free(p)
    for b in range(4):
        S.linear_of_e(c, f, lambda e, b=b: (dfun(e) >> b) & 1, db[b])
        c.free(db[b])
    S.thermo_erase(c, xq, f)
    c.stage("END")
    return c


def addr_tt(bit):
    v = 0
    for a in range(256):
        if (a >> bit) & 1:
            v |= 1 << a
    return v


def build_g2(P, tabs, angles=None, blk_a=6, blk_b=7, n_early=13, n_scr=8, seed=1):
    """G1 with pool lookups (address kept alive): tables are built just in time and erased for free.
    LUT1: H, c1 (into u), the top n_early bits of c0 (into lo) + scratch; PH1; H erased; LUT2: c2; A;
    LUT3: the low bits of c0; c2 erased; B; PH2; B^-1; LUT2r: c2 again; low c0 bits erased; A^-1;
    unLUT: free erasures while the pool is rich, then 22 ANDs for the monomials of the two address nibbles."""
    import random
    from invsqrt_gen_pool import Pool, NotQuadratic, tt_bit
    c = Circ(22)
    xq = list(range(22))
    dfun = lambda e: 10 - (e >> 1)
    c.stage("NORM")
    f = S.thermo(c, xq)
    p = S.linear_of_e(c, f, lambda e: e & 1)
    db = [S.linear_of_e(c, f, lambda e, b=b: (dfun(e) >> b) & 1) for b in range(4)]
    S.thermo_erase(c, xq, f)
    stages = [(16, db[3], False), (8, db[2], False), (4, db[1], False), (2, db[0], False), (1, p, True)]
    S.shifter(c, xq, stages)
    pool = Pool(c, random.Random(seed))
    abits = [xq[14 + k] for k in range(7)] + [p]            # address bit b lives on abits[b]
    for b, w in enumerate(abits):
        pool.set(w, addr_tt(b))
    n_late = P.wlo - n_early
    c.stage("LUT1")
    H, u, lo_e, scr = c.allocn(P.wh), c.allocn(P.wu), c.allocn(n_early), c.allocn(n_scr)
    tg = {H[h]: tt_bit(tabs["H"], h) for h in range(P.wh)}
    tg.update({u[k]: tt_bit(tabs["C1"], k) for k in range(P.wu)})
    tg.update({lo_e[k]: tt_bit(tabs["C0"], n_late + k) for k in range(n_early)})
    per1 = pool.synth(H + u + lo_e + scr, tg)
    for w in scr:
        pool.erase(w)
        c.free(w)
    tier_bits = [0, 1, 2, 3]
    cells = [[xq[0]], [xq[1], xq[2]], [xq[3], xq[4], xq[5], xq[6]], [xq[k] for k in range(3, 11)]]
    delays = list(range(11))
    c.stage("PH1")
    c.meta["checks"] = [("PH1", "H", H)]
    S.stream_m(c, H, tier_bits, db, cells, delays, 7 + 10, P.top_pos, angles)
    for w in H:
        pool.erase(w)
        c.free(w)
    c.stage("LUT2")
    c2 = c.allocn(P.wc2)
    per2 = pool.synth(c2, {c2[k]: tt_bit(tabs["C2"], k) for k in range(P.wc2)})
    for w in u:
        pool.drop(w)
    c.stage("A")
    pa, pb = {}, {}
    for i in P.a_rows:
        pa[i] = S.row_blocked(c, u, c2, P.a_shift(i), xq[i], 1, blk_a, "fwd")
    c.stage("LUT3")
    lo_l = c.allocn(n_late)
    per3 = pool.synth(lo_l, {lo_l[k]: tt_bit(tabs["C0"], k) for k in range(n_late)})
    lo = lo_l + lo_e
    for w in c2:
        pool.erase(w)
        c.free(w)
    for w in lo:
        pool.drop(w)
    c.stage("B")
    c.meta["checks"] += [("B", "u", u), ("PH2", "lo", lo)]
    for i in P.b_rows:
        pb[i] = S.row_blocked(c, lo, u, P.b_shift(i), xq[i], 1, blk_b, "fwd")
    c.stage("PH2")
    S.stream_m(c, lo, tier_bits, db, cells, delays, 25 + P.glo + 10, P.top_pos, angles)
    c.stage("B^-1")
    for i in reversed(P.b_rows):
        S.row_blocked(c, lo, u, P.b_shift(i), xq[i], 0, blk_b, "rev", pb[i])
    for k, w in enumerate(lo):
        pool.set(w, tt_bit(tabs["C0"], k))
    c.stage("LUT2r")
    c2 = c.allocn(P.wc2)
    per4 = pool.synth(c2, {c2[k]: tt_bit(tabs["C2"], k) for k in range(P.wc2)})
    for w in lo_l:
        pool.erase(w)
        c.free(w)
    c.stage("A^-1")
    for i in reversed(P.a_rows):
        S.row_blocked(c, u, c2, P.a_shift(i), xq[i], 0, blk_a, "rev", pa[i])
    for k, w in enumerate(u):
        pool.set(w, tt_bit(tabs["C1"], k))
    c.stage("unLUT")
    rest = c2 + lo_e + u
    mono = []
    while rest:
        w = rest[0]
        try:
            pool.erase(w)
        except NotQuadratic:
            assert not mono, "monomials did not suffice"
            for nib in (abits[:4], abits[4:]):
                mw = S.Monos(c, nib, [mk for mk in range(16) if bin(mk).count("1") >= 2])
                for mk, wire in mw.w.items():
                    if bin(mk).count("1") >= 2:
                        v = (1 << 256) - 1
                        for b in range(4):
                            if (mk >> b) & 1:
                                v &= pool.val[nib[b]]
                        pool.set(wire, v)
                mono.append(mw)
            continue
        c.free(w)
        rest.pop(0)
    for mw in mono:
        for mk, wire in mw.w.items():
            if bin(mk).count("1") >= 2:
                pool.drop(wire)
        mw.erase()
    c.meta["lut_ands"] = dict(LUT1=per1, LUT2=per2, LUT3=per3, LUT2r=per4)
    c.stage("UNNORM")
    S.shifter(c, xq, stages, inverse=True)
    f = S.thermo(c, xq)
    S.linear_of_e(c, f, lambda e: e & 1, p)
    c.free(p)
    for b in range(4):
        S.linear_of_e(c, f, lambda e, b=b: (dfun(e) >> b) & 1, db[b])
        c.free(db[b])
    S.thermo_erase(c, xq, f)
    c.stage("END")
    return c


def build_g3(P, tabs, **kw):
    from invsqrt_gen_g3 import build_g3 as b
    return b(P, tabs, **kw)


def build_g4(P, tabs, **kw):
    from invsqrt_gen_g4 import build_g4 as b
    return b(P, tabs, **kw)


def build_g5(P, tabs, **kw):
    from invsqrt_gen_g5 import build_g5 as b
    return b(P, tabs, **kw)


def build_g6(P, tabs, **kw):
    from invsqrt_gen_g6 import build_g6 as b
    return b(P, tabs, **kw)


VARIANTS = {"g0": build_g0, "g1": build_g1, "g2": build_g2, "g3": build_g3, "g4": build_g4, "g5": build_g5,
            "g6": build_g6}


def stage_index(c, name):
    for i, o in enumerate(c.ops):
        if o[0] == "@stage" and o[1] == name:
            return i
    raise KeyError(name)


def check(c, P, tabs, xs=None, hmr="rand", seed=1, verbose=True):
    """Exact run on xs (default: all x). Compares u and lo with the model after B, and the phase with the
    verifier's reference. Returns a dict."""
    if xs is None:
        xs = np.arange(1, M.MAX_X + 1, dtype=np.int64)
        target, tol = M.reference_all()
    else:
        xs = np.asarray(xs, dtype=np.int64)
        ta, tl = M.reference_all()
        target, tol = ta[xs - 1], tl[xs - 1]
    nq = c.num_qubits()
    sim = Sim(xs.astype(np.uint64), nq, hmr=hmr, seed=seed)
    ev = M.evaluate(P, tabs, xs)
    checks = c.meta.get("checks") or [("PH2", n, c.meta[n]) for n in ("u", "lo", "H")]
    checks = [(ch + (0,))[:4] for ch in checks]          # (stage, register, wires, dropped low bits)
    reg_ok, pos = {}, 0
    for idx, name, wires, sh in sorted((stage_index(c, st), n, w, s) for st, n, w, s in checks):
        sim.run(c.ops, pos, idx)
        pos = idx
        got = sim.value(wires).astype(np.int64)
        want = (ev[name] >> sh) & ((1 << len(wires)) - 1)
        reg_ok[name] = int((got != want).sum())
    sim.run(c.ops, pos)
    th = sim.phase_units()
    err = (th - target).view(np.int64)
    used = np.abs(err)
    ratio = used / tol.astype(np.float64)
    margin = tol.astype(np.int64) - used
    k, km = int(np.argmax(ratio)), int(np.argmin(margin))
    xr = sim.value(range(22)).astype(np.int64)
    dirty = [w for w in range(22, nq) if sim.q[w]]
    tof = sim.tof / sim.n
    rots = [a for a, _ in sim.rot]
    n_rot = sum(1 for a in rots if (a * 4) % 1 != 0)
    n_t = sum(1 for a in rots if (a * 4) % 1 == 0 and (a * 2) % 1 != 0)
    out = dict(ok=bool((used <= tol).all()) and not dirty and bool((xr == xs).all()) and not sim.bad_free,
               phase_ok=bool((used <= tol).all()), n_fail=int((used > tol).sum()), x_restored=bool((xr == xs).all()),
               dirty=dirty[:8], bad_free=sim.bad_free[:8], reg_mismatch=reg_ok, Q=nq, tof=tof, n_rot=n_rot, n_t=n_t,
               worst=float(ratio[k]), worst_x=int(xs[k]), min_margin=int(margin[km]), min_margin_x=int(xs[km]),
               margin_frac=float(margin[km] / M.reference_all()[1][-1]), n_ops=len(c.real_ops()))
    if out["min_margin"] > 0 and n_rot:
        eps = math.pi * out["min_margin"] / 2 ** 64 / n_rot
        tper = 0.56 * math.log2(1 / eps) + 5.3
        out["t_per_rot"] = tper
        out["toffoli_eq"] = tof + (n_t + n_rot * tper) / 4
        out["cost"] = nq * out["toffoli_eq"]
    out["stages"] = sim.stage_table()
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = sys.argv[1:]
    variant, outdir = args[0], args[1]
    rawkw = dict(a.split("=") for a in args[2:] if "=" in a)
    tabname = rawkw.pop("tabs", None)                               # tables from a file instead of our fit
    allkw = {k: int(v) for k, v in rawkw.items()}
    pfields = set(M.Params.__dataclass_fields__)
    kw = {k: v for k, v in allkw.items() if k in pfields}          # model parameters
    bkw = {k: v for k, v in allkw.items() if k not in pfields}     # builder options
    P = M.Params(**kw)
    t0 = time.time()
    tabs = M.load_tables(tabname) if tabname else M.fit_tables(P)
    if tabname:
        allkw["tabs"] = tabname
    c = VARIANTS[variant](P, tabs, **bkw)
    kw = allkw
    os.makedirs(outdir, exist_ok=True)
    kmx = os.path.join(outdir, "invsqrt.kmx")
    kmb = os.path.join(outdir, "invsqrt.kmb")
    c.write_kmx(kmx)
    kc = kmx_to_kmb(kmx, kmb)
    sha = hashlib.sha256(open(kmb, "rb").read()).hexdigest()
    st = c.static_counts()
    tabid = f"{tabname} ({tabs['sha']})" if tabname else P.key()
    print(f"variant {variant} params {kw} tables {tabid}  built in {time.time() - t0:.1f}s")
    print(f"  kmb sha256 {sha}  Q {kc.num_qubits}  bits {kc.num_bits}  ops {len(kc)}  "
          f"static CCX {st.get('CCX', 0)} CSWAP {st.get('CSWAP', 0)} CCZ {st.get('CCZ', 0)} Z_POW {st.get('Z_POW', 0)}")
    print("  peak live wires per stage: " + "  ".join(f"{n}:{p}" for n, p in c.stage_peaks))
    if c.meta.get("lut"):
        print("  lookup blocks: " + "  ".join(f"{k}: {v['ands']}/{v['outputs']}+{v['scratch']} ({v['source']})"
                                               for k, v in c.meta["lut"].items()))
    if "--no-check" in flags:
        return
    xs = None
    if "--sample" in flags:
        n = int(flags[flags.index("--sample") + 1])
        rng = np.random.default_rng(0)
        xs = np.unique(np.concatenate([rng.integers(1, M.MAX_X + 1, n), np.arange(1, 2049),
                                       np.arange(M.MAX_X - 2048, M.MAX_X + 1)]))
    seeds = int(flags[flags.index("--seeds") + 1]) if "--seeds" in flags else 1
    for sd in range(seeds):
        t1 = time.time()
        r = check(c, P, tabs, xs, seed=sd + 1)
        stg = r.pop("stages")
        print(f"  seed {sd + 1} [{time.time() - t1:.1f}s] " + json.dumps(r))
        if sd == 0:
            print("  stage table (mean executed Toffolis, live wires at entry):")
            print("   " + "  ".join(f"{n}:{t:.2f}/{l}" for n, t, l in stg))
            if c.meta.get("lut"):
                print("  lookup blocks: " + "  ".join(f"{k}: {v['ands']} ANDs / {v['outputs']} out + {v['scratch']} scr"
                                                       f" ({v['source']})" for k, v in c.meta["lut"].items()))
    with open(os.path.join(outdir, "build.json"), "w") as fh:
        json.dump(dict(variant=variant, params=kw, tables=tabid, sha256=sha, check=r, stages=stg), fh, indent=1)


if __name__ == "__main__":
    main()
