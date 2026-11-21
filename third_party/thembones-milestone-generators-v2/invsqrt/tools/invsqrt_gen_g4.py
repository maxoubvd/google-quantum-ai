"""invsqrt_gen: variant g4 (own code) = g3 (address erased, reflection shifter) with stage A rebuilt as a
product in a small measured register.

Stage A of g3 added the eight signed c2 rows straight into the 18-wire u register (17 ANDs per row). Here
(structure of the public file's stage A as read from a functional probe, own construction):
  q = sum_i sigma_i rr(c2, 13 - i) is accumulated in its own register, least significant row first. The partial
  sum after the rows with shifts 7..k is below 2^(10 - k) in magnitude, so the register grows by one
  sign-extension wire (a CNOT) per row and a row costs 10 - k ANDs: 52 for the eight rows.
  u17 = T1 - S' in ONE 17-bit subtraction (16 ANDs), with T1 = (C1 >> 1) + 512 looked up into u and
  S' = (q >> 1) + 512: bit 0 of q equals c2 bit 7 = bit 0 of C1 for every t (the tie of the model), so the
  model's u18 = C1 - q is even and u17 = C1/2 - (q >> 1) exactly [P].
  Then q is measured out; the phase is repaired for free in A^-1, which recomputes q (same values), adds S'
  back and uncomputes q by the mirrored rows: 52 + 16 + 52 plus block repairs.
The q register (11 wires) and its carries have to fit next to u, c2 and the early lo wires, so only n_e = 7 lo
bits are loaded before A and all but lo_rem = 3 lo wires are erased before A^-1.

Schedule: NORM, LUT1 [H, T1 -> u, n_e1 lo bits, scratch], erase address, PH1, erase H, LUT1b [rest of the early
lo bits, c2], A, LUT3 [late lo bits], erase c2, B, PH2, B^-1, LUT2r [c2], erase all lo wires but lo_rem, A^-1,
ADDR [address from u, c2, lo_rem wires], unLUT, UNNORM.
"""
from __future__ import annotations

import random
import sys

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_gen_parts as S
from invsqrt_gen_circ import Circ
from invsqrt_gen_pool import Pool, NotQuadratic, tt_bit, lookup_block
from invsqrt_gen_g3 import addr_tt


def build_g4(P, tabs, angles=None, blk_b=7, n_e=7, n_e1=6, s1=8, s1b=6, s_addr=10, lo_rem=3, blk_q=2, blk_s=3,
             blk_sr=6, cu_max=6, seed=1, deep=0, f1=0, f1b=0, fa=0):
    assert P.virt == 1 and P.tie_c1 == 1 and P.gu == 1 and P.fc2 == 1, "g4 needs the virtual-guard model"
    c = Circ(22)
    xq = list(range(22))
    dfun = lambda e: 10 - (e >> 1)
    nlo = P.wlo - 1
    nu = P.wu - 1
    kg = P.a_shift(P.a_lo)               # bit 0 of q equals bit kg of c2
    kmax = kg                            # largest A shift (row of t_{a_lo})
    wq0 = P.wc2 - kmax + 2               # q register width after its first row
    off = 1 << (P.wc2 + 2 - 2)           # 2^9: S' = (q >> 1) + off is non-negative
    # ---------------- NORM
    c.stage("NORM")
    f = S.thermo(c, xq)
    p = S.linear_of_e(c, f, lambda e: e & 1)
    db = [S.linear_of_e(c, f, lambda e, b=b: (dfun(e) >> b) & 1) for b in range(4)]
    S.thermo_erase(c, xq, f)
    S.shifter_refl(c, xq, p, db)
    c.x(xq[21])
    c.lend(xq[21])
    pool = Pool(c, random.Random(seed))
    names, lutlog = {}, {}
    abits = [xq[14 + k] for k in range(7)] + [p]
    for b, w in enumerate(abits):
        pool.set(w, addr_tt(b))
        names[w] = f"a{b}"

    def lut(block, outs, max_scr, empty=0):
        wires = c.allocn(len(outs))
        for (n, _), w in zip(outs, wires):
            names[w] = n
        inputs = sorted(((names[w], w) for w in pool.val if w not in wires), key=lambda t: t[0])
        scr = lookup_block(c, pool, block, inputs, [(n, w, v) for (n, v), w in zip(outs, wires)], max_scr,
                           seed=seed, log=lutlog, deep=deep, empty=empty)
        for i, w in enumerate(scr):
            names[w] = f"{block}.s{i}"
        return wires, scr

    def erase(wires, stable=None):
        """Free erasure; helper wires (one AND of two forms over `stable` and earlier helpers) when the pool is
        too small for a quadratic form. See invsqrt_gen_g3.py."""
        helpers = []
        hrng = random.Random(f"{seed}:helper:{len(c.ops)}")
        for w in wires:
            while True:
                try:
                    pool.erase(w)
                    break
                except NotQuadratic:
                    if stable is None or len(helpers) >= 24:
                        raise
                    base = sorted(stable) + helpers
                    A = {x for x in base if hrng.random() < 0.5} or {base[0]}
                    B = {x for x in base if hrng.random() < 0.5} or {base[-1]}
                    h = c.alloc()
                    pool.val[h] = 0
                    names[h] = f"helper{len(helpers)}"
                    pool.and_into(h, A, 0, B, 0)
                    helpers.append(h)
            names.pop(w, None)
            c.free(w)
        for h in reversed(helpers):
            pool.erase(h)
            names.pop(h)
            c.free(h)
        c.meta.setdefault("helpers", []).append(len(helpers))

    tier_bits = [0, 1, 2, 3]
    cells = [[xq[0]], [xq[1], xq[2]], [xq[3], xq[4], xq[5], xq[6]], [xq[k] for k in range(3, 11)]]
    delays = list(range(11))

    T1 = [(v >> 1) + off for v in tabs["C1"]]
    assert max(T1) < (1 << nu)
    n_late = nlo - n_e
    tH = [(f"H{h}", tt_bit(tabs["H"], h)) for h in range(P.wh)]
    tU = [(f"u{k}", tt_bit(T1, k - 1)) for k in range(1, P.wu)]
    tLo = [(f"lo{k}", tt_bit(tabs["C0"], k)) for k in range(1, P.wlo)]
    tC2 = [(f"c2_{k}", tt_bit(tabs["C2"], k)) for k in range(P.wc2)]
    lo_late_t, lo_early_t = tLo[:n_late], tLo[n_late:]
    # ---------------- LUT1, address erased
    c.stage("LUT1")
    w1, scr = lut("LUT1", tH + tU + lo_early_t[:n_e1], s1, f1)
    H, u = w1[:P.wh], w1[P.wh:P.wh + nu]
    lo_e = w1[P.wh + nu:]
    for w in abits[:7]:
        pool.erase(w)
        names.pop(w)
        c.lend(w)
    pool.erase(p)
    names.pop(p)
    c.free(p)
    # ---------------- PH1
    c.stage("PH1")
    c.meta["checks"] = [("PH1", "H", H, 0)]
    S.stream_m(c, H, tier_bits, db, cells, delays, 7 + 10, P.top_pos, angles)
    erase(H)
    # ---------------- rest of the early lo bits and c2
    c.stage("LUT2")
    w2, scr1b = lut("LUT1b", lo_early_t[n_e1:] + tC2, s1b, f1b)
    lo_e = lo_e + w2[:n_e - n_e1]
    c2 = w2[n_e - n_e1:]
    erase(scr + scr1b)

    # ---------------- stage A pieces
    def q_rows(pend_in):
        """q = sum_i sigma_i rr(c2, k_i) on fresh wires, least significant row first; boundary carries are
        measured with the repair left to a later run of the same rows. Returns (wires, {k: pending bits})."""
        q = c.allocn(wq0)
        out = {}
        for k in range(kmax, -1, -1):
            if k < kmax:
                w = c.alloc()
                c.cx(q[-1], w)                   # sign extension
                q.append(w)
            out[k] = S.row_blocked(c, q, c2, k, xq[P.nt - 1 - k], 0, blk_q, "fwd",
                                   pend_in[k] if pend_in else None)
        return q, out

    def u_add_s(q, sign, blk, mode, pend=None):
        """u += sign * S' with S' = (q >> 1) + off: operand = q[1..] with its top (sign) bit flipped."""
        opnd = q[1:]
        c.x(opnd[-1])
        if sign < 0:
            for a in u:
                c.x(a)
        out = S.add_blocked(c, u, opnd, None, blk, mode, pend)
        if sign < 0:
            for a in u:
                c.x(a)
        c.x(opnd[-1])
        return out

    def blk_u(n):
        return n if n - 1 <= cu_max else (n + 1) // 2

    # ---------------- A
    c.stage("A")
    for w in u:
        pool.drop(w)
    q, pq = q_rows(None)
    c.cx(c2[kg], q[0])                           # bit 0 of q is c2 bit kg for every input
    c.free(q[0])
    psub = u_add_s(q, -1, blk_s, "fwd")
    qm = {}
    for b in range(1, len(q)):
        qm[b] = c.hmr(q[b])                      # repaired in A^-1 when q is back on wires
        c.free(q[b])
    # ---------------- LUT3: late lo bits; c2 erased
    c.stage("LUT3")
    lo_l, _ = lut("LUT3", lo_late_t, 0)
    lo = lo_l + lo_e
    erase(c2, stable=lo)
    for w in lo:
        pool.drop(w)
    # ---------------- B, PH2, B^-1
    c.stage("B")
    c.meta["checks"] += [("B", "u", u, 1), ("PH2", "lo", lo, 1)]
    pb = {}
    for i in P.b_rows:
        pb[i] = S.row_virt(c, lo, u, P.b_shift(i), xq[i], 1, blk_b, "fwd")
    c.stage("PH2")
    S.stream_m(c, lo, tier_bits, db, cells, delays, 25 + P.glo + 10 - 1, P.top_pos, angles)
    c.stage("B^-1")
    for i in reversed(P.b_rows):
        S.row_virt(c, lo, u, P.b_shift(i), xq[i], 0, blk_b, "rev", pb[i])
    for k, w in enumerate(lo):
        pool.set(w, tLo[k][1])
    # ---------------- c2 again; lo erased down to lo_rem wires
    c.stage("LUT2r")
    c2, _ = lut("LUT2r", tC2, 0)
    # the lo wires kept through A^-1 must, together with c2, still tell all 256 addresses apart: otherwise the
    # other lo wires are not functions of what stays and cannot be erased at all
    # (c2 takes only 169 distinct values). No three single lo wires do it; three XOR forms of lo wires do
    # (invsqrt_gen_sep.py), and a XOR form costs only CNOTs.
    from invsqrt_gen_sep import separating_forms
    forms = separating_forms(tabs["C2"], tabs["C0"], nlo, lo_rem, seed=seed)
    assert forms is not None, "no lo_rem linear forms of lo separate the addresses together with c2"
    red = []                                     # reduced echelon: every pivot bit occurs in one form only
    for fm in forms:
        for pv, g in red:
            if (fm >> pv) & 1:
                fm ^= g
        pv = fm.bit_length() - 1
        red = [(p_, g ^ fm if (g >> pv) & 1 else g) for p_, g in red] + [(pv, fm)]
    lo_keep = []
    for i, (pv, fm) in enumerate(red):
        pool.xor_into(lo[pv], {lo[j] for j in range(nlo) if (fm >> j) & 1 and j != pv})
        names[lo[pv]] = f"lf{i}"
        lo_keep.append(lo[pv])
    erase([w for w in lo if w not in lo_keep], stable=lo_keep + c2)
    # ---------------- A^-1
    c.stage("A^-1")
    q, pq2 = q_rows(pq)
    c.cx(c2[kg], q[0])
    c.free(q[0])
    for b in range(1, len(q)):
        c.z(q[b], cond=qm[b])
        c.bit_free(qm[b])
    u_add_s(q, +1, blk_sr, "rev", psub)
    q[0] = c.alloc()
    c.cx(c2[kg], q[0])
    for k in range(0, kmax + 1):
        S.row_blocked(c, q, c2, k, xq[P.nt - 1 - k], 1, blk_u(len(q)), "rev", pq2[k])
        if k < kmax:
            c.cx(q[-2], q[-1])                   # the top wire is a sign-extension copy again
            c.free(q.pop())
    for w in q:
        c.free(w)
    for k, w in enumerate(u):
        pool.set(w, tU[k][1])
    # ---------------- address back
    c.stage("ADDR")
    wa, scr_a = lut("ADDR", [(f"a{b}", addr_tt(b)) for b in range(8)], s_addr, fa)
    abits = wa
    p = wa[7]
    # ---------------- unLUT
    c.stage("unLUT")
    rest = list(scr_a) + list(c2) + list(lo_keep) + list(u)
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
    cur = list(abits[:7])
    tgt = [xq[14 + k] for k in range(7)]
    if p < 22:
        w = c.alloc_high()
        c.swap(p, w)
        c.free(p)
        p = w
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
    S.shifter_refl(c, xq, p, db, inverse=True)
    f = S.thermo(c, xq)
    S.linear_of_e(c, f, lambda e: e & 1, p)
    c.free(p)
    for b in range(4):
        S.linear_of_e(c, f, lambda e, b=b: (dfun(e) >> b) & 1, db[b])
        c.free(db[b])
    S.thermo_erase(c, xq, f)
    c.stage("END")
    return c
