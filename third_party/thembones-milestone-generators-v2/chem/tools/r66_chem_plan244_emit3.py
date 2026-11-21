"""R66-CHEM-PLAN244: a copy of tools/r62_chem_bilin_emit2.py with ONE addition (marked TP): a helper's point may be
a LIST of points of D (plan entries of tools/r66_chem_plan244_mask3.py --tp): its low form is the sum of the
indicators (one AND as before). Everything else is unchanged; the pre-pass checks every wire on every run address.
Below: the header of the copied file.
R62-CHEM-BILIN: a copy of tools/r58_chem_kept_rp_emit68.py with ONE addition (marked S5): a point product may
read helpers and take the spill of the first point of D (plan entries of tools/r62_chem_bilin_mask2.py --s5 with
nine fields); the helpers are therefore erased AFTER the point stage. Plans with five-field entries emit as before.
Below: the header of the copied file.
R58-CHEM-KEPT re-pricer's copy of tools/r58_chem_kept_emit68.py. ONE change (marked RP): the pre-pass checks
the rows the verifier runs (rows < --nreal, 1,219), not the 13 padding rows 1,219-1,231 (never run); the pure
corrections are solved on those rows. Everything else is the designer's emitter unchanged.
Original docstring:
R58-CHEM-KEPT: EMITTER for a kept-product loader plan on builder C2's 68-AND HIGH CACHE (eight spare wires).

A copy of chem_ld2_c4.emit_c4 (C2's emitter for the 68-AND cache: the cache of H = {f : sum_V f = 0}, the high
parts of the controls and the pure corrections in the cache basis; imported read-only: build_hcache, TSpan,
hmask_truth) with the three additions of tools/r58_chem_kept_emit.py: helpers on free wires, reads, HMR + CZ.
Needs a plan that meets C2's conditions (tools/r58_chem_kept_c4chk.py --V ... --plan): every block form in H and
every wire's pure block part in H (tau == kappa); it fails loudly otherwise.

  CHEM_BROT=7 CHEM_HB=77 python -I r58_chem_kept_emit68.py PAD.npz PLAN.json --V 0x... --out k1.kmx
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_wwidth  # noqa: E402,F401
import chem_ld2_c4 as c4  # noqa: E402

bits = core.bits
NREAL = 1219                                      # RP: set from --nreal in main()


def emit_kept68(an, V, helpers, prodsA, pprods):
    """helpers: (low mask, block truth, fan-out); prodsA: (y1, a1, y2, a2, s, rF, rG) with a = high monomial mask;
    pprods: (low mask 1, block truth 1, low mask 2 or 0, block truth 2 or None, s): features (truth1, low1),
    (truth2, low2)."""
    BROT, J, QUBIT_CAP = core.BROT, core.J, core.QUBIT_CAP
    HB = c4.HB
    A, lb, L, nl, n_sel = an.A, an.lb, an.L, an.nl, an.n_sel
    assert an.H == HB and an.hbits == 7
    assert len(an.excl) == 1 and not an.excl[0], "single-phase low cache only"
    a = np.arange(A)
    lidx, hidx = a & (L - 1), a >> lb
    ALL = (1 << A) - 1
    LBf = []
    for bs in an.basis:
        tt = np.zeros(L, np.uint8)
        for m in bits(bs["anf"]):
            tt[m] = 1
        tt = core.mobius(tt, lb)
        LBf.append(core.to_int(tt[lidx]))

    def addr_truth(g):
        bv = np.array([(g >> p) & 1 for p in range(HB)], np.uint8)
        return core.to_int(bv[hidx])

    ops_c, truths = c4.build_hcache(V)
    nwc = len(truths)
    CT = c4.TSpan()
    CT.add((1 << HB) - 1, 1 << nwc)
    for w, f in enumerate(truths):
        CT.add(f, 1 << w)

    def hw(truth):
        r, tag = CT.reduce(truth)
        assert r == 0, "a block form is outside H"
        return tag

    def wire_q(w):
        return n_sel + BROT * w[0] + w[1]

    ops = []
    free = sorted(wire_q(c["wire"]) for c in an.carriers if c["clean"]) + list(range(n_sel + J * BROT, QUBIT_CAP))
    used_max = n_sel + J * BROT - 1

    def take():
        nonlocal used_max
        if not free:
            raise ValueError("out of qubits")
        q = free.pop(0)
        used_max = max(used_max, q)
        return q

    def give(q):
        free.append(q)
        free.sort()

    lowq = [None] * nl
    for t, bs in enumerate(an.basis):
        if bs["kind"] == "lin":
            lowq[t] = bs["bit"]
        elif bs["kind"] == "car":
            lowq[t] = wire_q(bs["wire"])
    prep = []
    for c in an.carriers:
        q = wire_q(c["wire"])
        if c["conj"]:
            prep.append(f"CX q{q - 1} q{q}")
        if c["const"]:
            prep.append(f"X q{q}")
        if c["clean"]:
            for t in c["combo"]:
                prep.append(f"CX q{lowq[t]} q{q}")
    ops += prep

    def controls(s1, s2, const1=False, const2=False):
        s1, s2 = set(s1), set(s2)
        assert s1 and s2 and s1 != s2, "degenerate controls"
        p1 = min(s1 - s2) if s1 - s2 else min(s1)
        pre = [f"CX q{q} q{p1}" for q in sorted(s1) if q != p1]
        if p1 in s2:
            s2 = (s1 ^ s2) | {p1}
            const2 = const2 ^ const1
        if const1:
            pre.append(f"X q{p1}")
        p2 = min(q for q in s2 if q != p1)
        pre += [f"CX q{q} q{p2}" for q in sorted(s2) if q != p2]
        if const2:
            pre.append(f"X q{p2}")
        return pre, p1, p2

    cq = [lb + i for i in range(7)]
    made_high = []
    for (t1, t2) in ops_c:
        pre, p1, p2 = controls([cq[w] for w in t1], [cq[w] for w in t2])
        t = take()
        ops.extend(pre)
        ops.append(f"CCX q{p1} q{p2} q{t}")
        ops.extend(reversed(pre))
        cq.append(t)
        made_high.append((pre, p1, p2, t))

    def low_and(t):
        par, i = an.recipe[t]
        qs = [lowq[x] for x in par]
        piv = next(q for q in qs if q != i)
        pre = [f"CX q{q} q{piv}" for q in qs if q != piv]
        tq = take()
        ops.extend(pre)
        ops.append(f"CCX q{piv} q{i} q{tq}")
        ops.extend(reversed(pre))
        lowq[t] = tq
        return (pre, piv, i, tq, t)

    def low_del(rec):
        pre, piv, i, tq, t = rec
        ops.extend(pre)
        ops.append(f"HMR q{tq} b0")
        ops.append(f"CZ q{piv} q{i} if b0")
        ops.extend(reversed(pre))
        lowq[t] = None
        give(tq)

    nw = len(an.out_wires)
    wq = [wire_q(w) for w in an.out_wires]

    def lowf(mask):
        f = 0
        for t in bits(mask):
            f ^= LBf[t]
        return f

    def targets(s):
        return [q for q in range(nw) if (an.G[q] & s).bit_count() & 1]

    # ---- classical pre-pass
    def hform(lo, cl):                            # TP: a helper form = low form + block form (either may be absent)
        return lowf(lo) ^ (addr_truth(cl) if cl is not None else 0)

    hv = [hform(lo1, cl1) & hform(lo2, cl2) for (lo1, cl1, lo2, cl2, _) in helpers]
    acc = [0] * nw
    for j, (_, _, _, _, fan) in enumerate(helpers):
        for q in targets(fan):
            acc[q] ^= hv[j]
    for (y1, a1, y2, a2, s, rF, rG) in prodsA:
        f1 = lowf(y1) ^ addr_truth(c4.hmask_truth(a2))
        f2 = lowf(y2) ^ addr_truth(c4.hmask_truth(a1))
        for j in rF:
            f1 ^= hv[j]
        for j in rG:
            f2 ^= hv[j]
        P = f1 & f2
        tg = targets(s)
        assert tg
        for q in tg:
            acc[q] ^= P
    for (l1, t1, l2, t2, s, r1, r2) in pprods:
        f1 = lowf(l1) ^ (addr_truth(t2) if t2 is not None else 0)
        f2 = lowf(l2) ^ addr_truth(t1)
        for j in r1:                                  # S5
            f1 ^= hv[j]
        for j in r2:
            f2 ^= hv[j]
        P = f1 & f2
        tg = targets(s)
        assert tg
        for q in tg:
            acc[q] ^= P
    REAL = (1 << NREAL) - 1                       # RP: the verifier runs rows < NREAL only
    pure = core.Span()
    pure.add(ALL & REAL, 1)
    for t in range(nl):
        pure.add(LBf[t] & REAL, 1 << (1 + t))
    for w, f in enumerate(truths):
        pure.add(addr_truth(f) & REAL, 1 << (1 + nl + w))
    corr = []
    for q in range(nw):
        r, c = pure.reduce((an.delta_int[q] ^ acc[q]) & REAL)
        assert r == 0, f"wire {an.out_wires[q]}: not loaded, or its pure block part is outside H (tau != kappa)"
        corr.append(c)

    core_made = [low_and(t) for t in an.core]
    n_and = 0
    cmask = (1 << nwc) - 1

    def hset(tag):
        return {cq[w] for w in bits(tag & cmask)}, bool(tag >> nwc & 1)

    # ---- helpers
    hrec = []
    for (lo1, cl1, lo2, cl2, fan) in helpers:     # TP: (low 1 + block 1) AND (low 2 + block 2)
        h1, k1 = hset(hw(cl1)) if cl1 is not None else (set(), False)
        h2, k2 = hset(hw(cl2)) if cl2 is not None else (set(), False)
        s1 = {lowq[t] for t in bits(lo1)} | h1
        s2 = {lowq[t] for t in bits(lo2)} | h2
        assert None not in s1 and None not in s2
        pre, p1, p2 = controls(s1, s2, k1, k2)
        t = take()
        ops.extend(pre)
        ops.append(f"CCX q{p1} q{p2} q{t}")
        ops.extend(reversed(pre))
        for q in targets(fan):
            ops.append(f"CX q{t} q{wq[q]}")
        hrec.append((pre, p1, p2, t))
        n_and += 1

    def do_product(low1, tag1, low2, tag2, s, rF=(), rG=()):
        nonlocal n_and
        h1, k1 = hset(tag1)
        h2, k2 = hset(tag2)
        s1 = {lowq[t] for t in bits(low1)} | h1 | {hrec[j][3] for j in rF}
        s2 = {lowq[t] for t in bits(low2)} | h2 | {hrec[j][3] for j in rG}
        assert None not in s1 and None not in s2
        pre, p1, p2 = controls(s1, s2, k1, k2)
        tg = targets(s)
        first = wq[tg[0]]
        fan = [f"CX q{first} q{wq[q]}" for q in tg[1:]]
        ops.extend(pre)
        ops.extend(fan)
        ops.append(f"CCX q{p1} q{p2} q{first}")
        ops.extend(reversed(fan))
        ops.extend(reversed(pre))
        n_and += 1

    for (y1, a1, y2, a2, s, rF, rG) in prodsA:
        do_product(y1, hw(c4.hmask_truth(a2)), y2, hw(c4.hmask_truth(a1)), s, rF, rG)
    for (l1, t1, l2, t2, s, r1, r2) in pprods:        # S5: the point stage runs while the helpers are alive
        do_product(l1, hw(t2) if t2 is not None else 0, l2, hw(t1), s, r1, r2)
    for (pre, p1, p2, t) in reversed(hrec):
        ops.extend(pre)
        ops.append(f"HMR q{t} b0")
        ops.append(f"CZ q{p1} q{p2} if b0")
        ops.extend(reversed(pre))
        give(t)
    for q in range(nw):
        for t in range(nl):
            if corr[q] >> (1 + t) & 1:
                ops.append(f"CX q{lowq[t]} q{wq[q]}")
    for q in range(nw):
        c = corr[q]
        if c & 1:
            ops.append(f"X q{wq[q]}")
        for w in range(nwc):
            if c >> (1 + nl + w) & 1:
                ops.append(f"CX q{cq[w]} q{wq[q]}")
    for rec in reversed(core_made):
        low_del(rec)
    for pre, p1, p2, t in reversed(made_high):
        ops.extend(pre)
        ops.append(f"HMR q{t} b0")
        ops.append(f"CZ q{p1} q{p2} if b0")
        ops.extend(reversed(pre))
    ops.extend(reversed(prep))
    head = []
    for i in range(n_sel):
        head.append(f"APPEND_TO_REGISTER q{i} r0")
    head.append('REGISTER r0 "selection"')
    for j in range(J):
        for b in range(BROT):
            head.append(f"APPEND_TO_REGISTER q{n_sel + BROT * j + b} r{j + 1}")
        head.append(f'REGISTER r{j + 1} "word_{j}"')
    infra = len(ops_c) + len(an.comp)
    stats = dict(k=an.k, infra=infra, products=n_and, toffoli=infra + n_and, qubits=used_max + 1,
                 helpers=len(helpers), ops=len(ops))
    return "\n".join(head + ops) + "\n", stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pad")
    ap.add_argument("plan")
    ap.add_argument("--V", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--nreal", type=int, default=1219, help="RP: rows the verifier runs (the pre-pass checks these)")
    a = ap.parse_args()
    global NREAL
    NREAL = a.nreal
    core.J = 5
    words = np.asarray(np.load(a.pad)["angle_words"], np.int64)
    zt = np.load(a.pad)                         # RP guard: the run rows must match the table's N + R B
    if "W" in zt.files and "U" in zt.files:
        nt = int(zt["U"].shape[-1] + zt["W"].shape[0] * zt["W"].shape[1])
        assert NREAL == nt, f"--nreal {NREAL} but the table has {nt} rows (N + R B)"
    assert NREAL <= words.shape[0]
    pl = json.load(open(a.plan))
    an = core.analyze(words, pl["k"], 4)
    nl = an.nl
    V = int(a.V, 16)
    LB = []
    for bs in an.basis:
        tt = 0
        for l in range(16):
            v = 0
            for m in bits(bs["anf"]):
                if (l & m) == m:
                    v ^= 1
            tt |= v << l
        LB.append(tt)
    sp = core.Span()
    for t in range(nl):
        sp.add(LB[t], 1 << t)

    def dmask(p):
        r, c = sp.reduce(1 << p)
        assert r == 0
        return c

    full = (1 << c4.HB) - 1
    def dmask_any(p):                             # TP: a helper on several points of D
        if isinstance(p, list):
            m = 0
            for x in p:
                m ^= dmask(x)
            return m
        return dmask(p)

    helpers = []
    for (p, cl, fan) in pl["helpers"]:
        if isinstance(p, list):                   # TP: pair helper = (1_A + delta_p3)(1_B + delta_p2), A, B disjoint
            assert len(p) == 2 and int(cl[0]) & int(cl[1]) == 0
            helpers.append((dmask(p[1]), int(cl[0]), dmask(p[0]), int(cl[1]), int(fan)))
        else:
            helpers.append((dmask(p), None, 0, int(cl), int(fan)))
    prodsA = [(int(y1), int(a1), int(y2), int(a2), int(s), list(rF), list(rG))
              for (y1, a1, y2, a2, s, ph, rF, rG) in pl["stageA"]]
    pprods = []
    for ent in pl["point"]:
        p1, c1, p2, c2, s = ent[:5]
        r1, r2, sp1, sp2 = (list(ent[5]), list(ent[6]), int(ent[7]), int(ent[8])) if len(ent) > 5 else ([], [], 0, 0)
        t1, t2 = int(c1), (int(c2) if p2 is not None else None)
        ok1 = t1 not in (0, full)
        ok2 = t2 is not None and t2 not in (0, full)
        d0 = dmask(pl["D0"]) if (sp1 or sp2) else 0     # S5: the spill of the first point of D
        if ok1 and ok2:
            pprods.append((dmask(p1) ^ (d0 if sp1 else 0), t1, dmask(p2) ^ (d0 if sp2 else 0), t2, int(s), r1, r2))
        else:
            assert not (r1 or r2 or sp1 or sp2), "S5: a degenerate point product with reads"
            if ok1:
                pprods.append((dmask(p1), t1, 0, None, int(s), [], []))
            elif ok2:
                pprods.append((dmask(p2), t2, 0, None, int(s), [], []))
    text, st = emit_kept68(an, V, helpers, prodsA, pprods)
    open(a.out, "w").write(text)
    print(f"EMITTED k {pl['k']} on the 68-AND cache: {st['toffoli']} Toffolis = infra {st['infra']} + "
          f"{st['products']} ANDs ({len(prodsA)} stage A + {st['helpers']} helpers + {len(pprods)} point), qubits "
          f"{st['qubits']}; every wire checked on every RUN address (rows < {NREAL}) by the pre-pass", flush=True)


if __name__ == "__main__":
    main()
