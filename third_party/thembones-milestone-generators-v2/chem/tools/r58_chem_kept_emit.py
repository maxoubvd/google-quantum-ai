"""R58-CHEM-KEPT: EMITTER for a loader plan with KEPT PRODUCTS (helpers on spare wires, read by later forms).

A copy of chem_loaders_core.emit (the builders' emitter: caches, products, corrections, every wire checked on every
address) with three additions:
  * a HELPER is a product whose CCX lands on a FREE wire (not on an output wire); its fan-out is plain CNOTs from
    that wire; it is erased after its last reader by HMR + CZ on its two forms (no Toffoli);
  * a READ: the form of a later product adds the helper's wire (one more CNOT into the form's pivot);
  * the classical pre-pass evaluates the helpers and the reads.
Input: the padded table and the plan written by tools/r58_chem_kept_mask.py --plan.

  CHEM_BROT=7 python -I r58_chem_kept_emit.py PAD.npz PLAN.json --out k1.kmx
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402

bits, to_int, mobius, Span = core.bits, core.to_int, core.mobius, core.Span


def emit_kept(an, helpers, prodsA, pprods):
    """helpers: (low mask, high mask, fan-out); prodsA: (y1, a1, y2, a2, s, rF, rG) with form1 = low y1 + high a2
    + helpers rF, form2 = low y2 + high a1 + helpers rG; pprods: (y1, a1, y2, a2, s) after the helpers are erased."""
    BROT, J, QUBIT_CAP = core.BROT, core.J, core.QUBIT_CAP
    A, lb, L, nl, n_sel = an.A, an.lb, an.L, an.nl, an.n_sel
    assert len(an.excl) == 1 and not an.excl[0], "this emitter handles the single-phase low cache only"
    a = np.arange(A)
    lidx, hidx = a & (L - 1), a >> lb
    ALL = (1 << A) - 1
    LBf = []
    for bs in an.basis:
        tt = np.zeros(L, np.uint8)
        for m in bits(bs["anf"]):
            tt[m] = 1
        tt = mobius(tt, lb)
        LBf.append(to_int(tt[lidx]))
    HMf = [None] + [to_int((hidx & M) == M) for M in range(1, an.H)]

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
    highq = {}
    made_high = []
    for i in range(an.hbits):
        if (1 << i) < an.H:
            highq[1 << i] = lb + i
    for M in an.hi_list:
        top = M.bit_length() - 1
        t = take()
        c1, c2 = highq[M ^ (1 << top)], lb + top
        ops.append(f"CCX q{c1} q{c2} q{t}")
        highq[M] = t
        made_high.append((c1, c2, t))

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

    def func(low_mask, high_mask):
        f = 0
        for t in bits(low_mask):
            f ^= LBf[t]
        for M in bits(high_mask):
            f ^= HMf[M + 1]
        return f

    def targets(s):
        return [q for q in range(nw) if (an.G[q] & s).bit_count() & 1]

    # ---- classical pre-pass (helpers, reads, products), then the pure correction every wire still needs
    hv = [func(lo, 0) & func(0, hi) for (lo, hi, _) in helpers]
    acc = [0] * nw
    for j, (_, _, fan) in enumerate(helpers):
        for q in targets(fan):
            acc[q] ^= hv[j]
    for (y1, a1, y2, a2, s, rF, rG) in prodsA:
        f1, f2 = func(y1, a2), func(y2, a1)
        for j in rF:
            f1 ^= hv[j]
        for j in rG:
            f2 ^= hv[j]
        P = f1 & f2
        tg = targets(s)
        assert tg, "product without a target"
        for q in tg:
            acc[q] ^= P
    for (y1, a1, y2, a2, s) in pprods:
        P = func(y1, a2) & func(y2, a1)
        tg = targets(s)
        assert tg, "point product without a target"
        for q in tg:
            acc[q] ^= P
    pure = Span()
    pure.add(ALL, 1)
    for t in range(nl):
        pure.add(LBf[t], 1 << (1 + t))
    for M in range(1, an.H):
        pure.add(HMf[M], 1 << (nl + M))
    corr = []
    for q in range(nw):
        r, c = pure.reduce(an.delta_int[q] ^ acc[q])
        assert r == 0, f"wire {an.out_wires[q]} is not loaded correctly (circuit {an.k})"
        corr.append(c)

    core_made = [low_and(t) for t in an.core]
    n_and = 0

    def fold(s1, s2):
        """CNOTs that leave form 1 on wire p1 and form 2 on wire p2 (the builders' pivot rule)."""
        assert s1 and s2 and s1 != s2, "degenerate product"
        p1 = min(s1 - s2) if s1 - s2 else min(s1)
        pre = [f"CX q{q} q{p1}" for q in sorted(s1) if q != p1]
        if p1 in s2:
            s2 = (s1 ^ s2) | {p1}
        p2 = min(q for q in s2 if q != p1)
        pre += [f"CX q{q} q{p2}" for q in sorted(s2) if q != p2]
        return pre, p1, p2

    # ---- helpers: CCX onto a free wire, fan-out by CNOTs from it
    hrec = []
    for (lo, hi, fan) in helpers:
        s1 = {lowq[t] for t in bits(lo)}
        s2 = {highq[M + 1] for M in bits(hi)}
        assert None not in s1
        pre, p1, p2 = fold(s1, s2)
        t = take()
        ops.extend(pre)
        ops.append(f"CCX q{p1} q{p2} q{t}")
        ops.extend(reversed(pre))
        for q in targets(fan):
            ops.append(f"CX q{t} q{wq[q]}")
        hrec.append((pre, p1, p2, t))
        n_and += 1

    def do_product(y1, a1, y2, a2, s, rF=(), rG=()):
        nonlocal n_and
        s1 = {lowq[t] for t in bits(y1)} | {highq[M + 1] for M in bits(a2)} | {hrec[j][3] for j in rF}
        s2 = {lowq[t] for t in bits(y2)} | {highq[M + 1] for M in bits(a1)} | {hrec[j][3] for j in rG}
        assert None not in s1 and None not in s2
        pre, p1, p2 = fold(s1, s2)
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
        do_product(y1, a1, y2, a2, s, rF, rG)
    # ---- erase the helpers (HMR + CZ on the two forms: no Toffoli)
    for (pre, p1, p2, t) in reversed(hrec):
        ops.extend(pre)
        ops.append(f"HMR q{t} b0")
        ops.append(f"CZ q{p1} q{p2} if b0")
        ops.extend(reversed(pre))
        give(t)
    for (y1, a1, y2, a2, s) in pprods:
        do_product(y1, a1, y2, a2, s)
    # ---- pure corrections, cache erasure (as the builders' emitter)
    for q in range(nw):
        for t in range(nl):
            if corr[q] >> (1 + t) & 1:
                ops.append(f"CX q{lowq[t]} q{wq[q]}")
    for q in range(nw):
        c = corr[q]
        if c & 1:
            ops.append(f"X q{wq[q]}")
        for M in range(1, an.H):
            if c >> (nl + M) & 1:
                ops.append(f"CX q{highq[M]} q{wq[q]}")
    for rec in reversed(core_made):
        low_del(rec)
    for c1, c2, t in reversed(made_high):
        ops.append(f"HMR q{t} b0")
        ops.append(f"CZ q{c1} q{c2} if b0")
    ops.extend(reversed(prep))
    head = []
    for i in range(n_sel):
        head.append(f"APPEND_TO_REGISTER q{i} r0")
    head.append('REGISTER r0 "selection"')
    for j in range(J):
        for b in range(BROT):
            head.append(f"APPEND_TO_REGISTER q{n_sel + BROT * j + b} r{j + 1}")
        head.append(f'REGISTER r{j + 1} "word_{j}"')
    stats = dict(k=an.k, infra=an.infra, products=n_and, toffoli=an.infra + n_and, qubits=used_max + 1,
                 helpers=len(helpers), ops=len(ops))
    return "\n".join(head + ops) + "\n", stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pad")
    ap.add_argument("plan")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    core.BROT, core.J = int(os.environ.get("CHEM_BROT", "7")), 5
    words = np.asarray(np.load(a.pad)["angle_words"], np.int64)
    pl = json.load(open(a.plan))
    an = core.analyze(words, pl["k"], 4)
    nl, NH = an.nl, an.H
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
    sp = Span()
    for t in range(nl):
        sp.add(LB[t], 1 << t)

    def dmask(p):
        r, c = sp.reduce(1 << p)
        assert r == 0
        return c

    def hmask(c):
        """block function (bit h = value on block h, zero on block 0) -> mask over the high monomials M - 1."""
        if c & 1:
            c ^= (1 << NH) - 1                    # the constant is a pure low part of the product: free
        m = 0
        for M in range(1, NH):
            v = 0
            sub = M
            while True:
                v ^= (c >> sub) & 1
                if sub == 0:
                    break
                sub = (sub - 1) & M
            if v:
                m |= 1 << (M - 1)
        return m

    helpers = [(dmask(p), hmask(int(cl)), int(fan)) for (p, cl, fan) in pl["helpers"]]
    prodsA = [(int(y1), int(a1), int(y2), int(a2), int(s), list(rF), list(rG))
              for (y1, a1, y2, a2, s, ph, rF, rG) in pl["stageA"]]
    pprods = []
    for (p1, c1, p2, c2, s) in pl["point"]:
        h1 = hmask(int(c1))
        h2 = hmask(int(c2)) if p2 is not None else 0
        if h1 and h2:
            pprods.append((dmask(p1), h1, dmask(p2), h2, int(s)))
        elif h1:
            pprods.append((dmask(p1), h1, 0, 0, int(s)))
        elif h2:
            pprods.append((dmask(p2), h2, 0, 0, int(s)))
    text, st = emit_kept(an, helpers, prodsA, pprods)
    open(a.out, "w").write(text)
    print(f"EMITTED k {pl['k']}: {st['toffoli']} Toffolis = infra {st['infra']} + {st['products']} ANDs "
          f"({len(prodsA)} stage A + {st['helpers']} helpers + {len(pprods)} point), qubits {st['qubits']}; every wire "
          f"checked on every address by the pre-pass", flush=True)


if __name__ == "__main__":
    main()
