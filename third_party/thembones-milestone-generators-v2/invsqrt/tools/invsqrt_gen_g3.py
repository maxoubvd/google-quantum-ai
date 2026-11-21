"""invsqrt_gen: variant g3 (own code). The register footprint of the public file with our own numerics:
lo on 20 wires (its second guard bit is virtual), u on 17 wires during B (its guard bit exists during A only),
c2 9 bits, H 8 bits; every lookup goes through the replaceable block interface of invsqrt_gen_pool.py; the
always-one wire q21 is lent out after normalisation; optionally the address wires are erased after the first
lookup and rebuilt before the final unlookup (erase_addr=1), as the public file does.

Schedule (erase_addr = 0):  NORM, LUT1 [H, c1, top n_e bits of c0, scratch], PH1, erase H, LUT2 [c2], A,
  LUT3 [low bits of c0], erase c2, B, PH2, B^-1, LUT2r [c2], erase low c0 bits, A^-1, unLUT, UNNORM.
Schedule (erase_addr = 1):  NORM, LUT1 [H, c1, n_e1 bits of c0, scratch], erase address, PH1, erase H,
  LUT1b [rest of the top c0 bits, c2], erase scratch, A, LUT3, erase c2, B, PH2, B^-1, LUT2r, erase low c0 bits,
  A^-1, ADDR [address rebuilt from the table wires], unLUT, UNNORM.
"""
from __future__ import annotations

import sys

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_gen_parts as S
from invsqrt_gen_circ import Circ
from invsqrt_gen_pool import Pool, NotQuadratic, tt_bit, lookup_block


def addr_tt(bit):
    v = 0
    for a in range(256):
        if (a >> bit) & 1:
            v |= 1 << a
    return v


def build_g3(P, tabs, angles=None, blk_a=6, blk_b=7, n_e=11, n_e1=6, s1=10, erase_addr=0, seed=1, refl=0,
             keep_c2=0, s3=0, s2r=0):
    assert P.virt == 1 and P.tie_c1 == 1 and P.gu == 1 and P.fc2 == 1, "g3 needs the virtual-guard model"
    import random
    c = Circ(22)
    xq = list(range(22))
    dfun = lambda e: 10 - (e >> 1)
    nlo = P.wlo - 1                      # lo wires: bits 1..wlo-1 of the model's lo
    nu = P.wu - 1                        # u wires during B: bits 1..wu-1 of the model's u
    kg = P.a_shift(P.a_lo)               # the guard bit of u starts as bit kg of c2
    # ---------------- NORM
    c.stage("NORM")
    f = S.thermo(c, xq)
    p = S.linear_of_e(c, f, lambda e: e & 1)
    db = [S.linear_of_e(c, f, lambda e, b=b: (dfun(e) >> b) & 1) for b in range(4)]
    S.thermo_erase(c, xq, f)
    stages = [(16, db[3], False), (8, db[2], False), (4, db[1], False), (2, db[0], False), (1, p, True)]
    if refl:
        S.shifter_refl(c, xq, p, db)
    else:
        S.shifter(c, xq, stages)
    c.x(xq[21])                          # xn bit 21 is always 1: clear it and lend the wire
    c.lend(xq[21])
    pool = Pool(c, random.Random(seed))
    names, lutlog = {}, {}
    abits = [xq[14 + k] for k in range(7)] + [p]
    for b, w in enumerate(abits):
        pool.set(w, addr_tt(b))
        names[w] = f"a{b}"

    def lut(block, outs, max_scr):
        """outs: list of (name, truth table); allocates the wires, runs the block, returns (wires, scratch)."""
        wires = c.allocn(len(outs))
        for (n, _), w in zip(outs, wires):
            names[w] = n
        inputs = sorted(((names[w], w) for w in pool.val if w not in wires), key=lambda t: t[0])
        scr = lookup_block(c, pool, block, inputs, [(n, w, v) for (n, v), w in zip(outs, wires)], max_scr,
                           seed=seed, log=lutlog)
        for i, w in enumerate(scr):
            names[w] = f"{block}.s{i}"
        return wires, scr

    def erase(wires, keep_failures=False, stable=None):
        """Free erasure of table wires. If a wire is not a quadratic form over the live pool and `stable` is
        given (pool wires that stay alive), helper wires are added: each is one AND of two random forms over the
        stable wires and earlier helpers, so it enlarges the span by products of degree 3 and more and can
        itself be erased for free afterwards."""
        left, helpers = [], []
        hrng = random.Random(f"{seed}:helper:{len(c.ops)}")
        for w in wires:
            while True:
                try:
                    pool.erase(w)
                    break
                except NotQuadratic:
                    if stable is None or len(helpers) >= 16:
                        if not keep_failures:
                            raise
                        left.append(w)
                        w = None
                        break
                    base = sorted(stable) + helpers
                    A = {x for x in base if hrng.random() < 0.5} or {base[0]}
                    B = {x for x in base if hrng.random() < 0.5} or {base[-1]}
                    h = c.alloc()
                    pool.val[h] = 0
                    names[h] = f"helper{len(helpers)}"
                    pool.and_into(h, A, 0, B, 0)
                    helpers.append(h)
            if w is not None:
                names.pop(w, None)
                c.free(w)
        for h in reversed(helpers):
            pool.erase(h)
            names.pop(h)
            c.free(h)
        c.meta.setdefault("helpers", []).append(len(helpers))
        return left

    def stream_ctx():
        tier_bits = [0, 1, 2, 3]
        cells = [[xq[0]], [xq[1], xq[2]], [xq[3], xq[4], xq[5], xq[6]], [xq[k] for k in range(3, 11)]]
        return tier_bits, cells, list(range(11))

    n_late = nlo - n_e
    tH = [(f"H{h}", tt_bit(tabs["H"], h)) for h in range(P.wh)]
    tU = [(f"u{k}", tt_bit(tabs["C1"], k)) for k in range(1, P.wu)]
    tLo = [(f"lo{k}", tt_bit(tabs["C0"], k)) for k in range(1, P.wlo)]      # index k-1 in the lo register
    tC2 = [(f"c2_{k}", tt_bit(tabs["C2"], k)) for k in range(P.wc2)]
    lo_late_t, lo_early_t = tLo[:n_late], tLo[n_late:]
    # ---------------- LUT1
    c.stage("LUT1")
    if erase_addr:
        first = tH + tU + lo_early_t[:n_e1]
    else:
        first = tH + tU + lo_early_t
    w1, scr = lut("LUT1", first, s1)
    H, u = w1[:P.wh], w1[P.wh:P.wh + nu]
    lo_e = w1[P.wh + nu:]
    if erase_addr:
        for w in abits[:7]:
            pool.erase(w)
            names.pop(w)
            c.lend(w)
        pool.erase(p)
        names.pop(p)
        c.free(p)
    else:
        erase(scr)
        scr = []
    # ---------------- PH1
    tier_bits, cells, delays = stream_ctx()
    c.stage("PH1")
    c.meta["checks"] = [("PH1", "H", H, 0)]
    S.stream_m(c, H, tier_bits, db, cells, delays, 7 + 10, P.top_pos, angles)
    erase(H)
    # ---------------- c2 (and, in the address-free schedule, the rest of the early lo bits)
    c.stage("LUT2")
    if erase_addr:
        w2, scr2 = lut("LUT1b", lo_early_t[n_e1:] + tC2, 0)
        lo_e = lo_e + w2[:n_e - n_e1]
        c2 = w2[n_e - n_e1:]
        erase(scr)
        scr = []
    else:
        c2, _ = lut("LUT2", tC2, 0)
    for w in u:
        pool.drop(w)
    # ---------------- A
    c.stage("A")
    g = c.alloc()
    c.cx(c2[kg], g)
    pa, pb = {}, {}
    for i in P.a_rows:
        pa[i] = S.row_blocked(c, [g] + u, c2, P.a_shift(i), xq[i], 1, blk_a, "fwd")
    c.free(g)                            # the guard bit of u is |0> after A for every input
    # ---------------- LUT3: low bits of c0
    c.stage("LUT3")
    lo_l, scr3 = lut("LUT3", lo_late_t, s3)
    lo = lo_l + lo_e
    erase(scr3)
    erase(c2, stable=lo)
    for w in lo:
        pool.drop(w)
    # ---------------- B, PH2, B^-1
    c.stage("B")
    c.meta["checks"] += [("B", "u", u, 1), ("PH2", "lo", lo, 1)]
    for i in P.b_rows:
        pb[i] = S.row_virt(c, lo, u, P.b_shift(i), xq[i], 1, blk_b, "fwd")
    c.stage("PH2")
    S.stream_m(c, lo, tier_bits, db, cells, delays, 25 + P.glo + 10 - 1, P.top_pos, angles)
    c.stage("B^-1")
    for i in reversed(P.b_rows):
        S.row_virt(c, lo, u, P.b_shift(i), xq[i], 0, blk_b, "rev", pb[i])
    for k, w in enumerate(lo):
        pool.set(w, tLo[k][1])
    # ---------------- c2 again, low lo bits erased, A^-1
    c.stage("LUT2r")
    c2, scr2r = lut("LUT2r", tC2, s2r)
    erase(scr2r)
    left = erase(lo_l, keep_failures=True, stable=lo_e + c2)
    c.stage("A^-1")
    g = c.alloc()
    for i in reversed(P.a_rows):
        S.row_blocked(c, [g] + u, c2, P.a_shift(i), xq[i], 0, blk_a, "rev", pa[i])
    c.cx(c2[kg], g)
    c.free(g)
    for k, w in enumerate(u):
        pool.set(w, tU[k][1])
    erase(left)
    # ---------------- address back (address-free schedule)
    if erase_addr:
        c.stage("ADDR")
        if not keep_c2:
            erase(c2)
            c2 = []
        wa, _ = lut("ADDR", [(f"a{b}", addr_tt(b)) for b in range(8)], 0)
        abits = wa
        p = wa[7]
    # ---------------- unLUT
    c.stage("unLUT")
    rest = list(c2) + list(lo_e) + list(u)
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
        names.pop(w, None)
        c.free(w)
        rest.pop(0)
    for mw in mono:
        for mk, wire in mw.w.items():
            if bin(mk).count("1") >= 2:
                pool.drop(wire)
        mw.erase()
    c.meta["lut"] = lutlog
    # ---------------- UNNORM
    c.stage("UNNORM")
    if erase_addr:
        # move the rebuilt address bits onto q14..q20 (SWAPs are free) and the parity onto a non-I/O wire
        cur = list(abits[:7])
        tgt = [xq[14 + k] for k in range(7)]
        if p < 22:
            q = c.alloc_high()
            c.swap(p, q)
            c.free(p)
            p = q
        for k in range(7):
            s_, t_ = cur[k], tgt[k]
            if s_ == t_:
                continue
            if c.is_free(t_):
                c.reclaim(t_)
            else:
                assert t_ in cur, "target wire is busy"
            c.swap(s_, t_)
            for j in range(7):
                if cur[j] == t_ and j != k:
                    cur[j] = s_
            cur[k] = t_
            if s_ not in cur:
                c.free(s_)
        for t_ in tgt:
            c.lendable.discard(t_)
    c.reclaim(xq[21])
    c.x(xq[21])
    stages = stages[:4] + [(1, p, True)]         # p may live on a new wire after the address was rebuilt
    if refl:
        S.shifter_refl(c, xq, p, db, inverse=True)
    else:
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
