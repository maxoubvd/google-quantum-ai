"""invsqrt_gen: variant g6 (own code) = g5 (62 wires, parked exponent) with stage A as rows on the low window
of u, 13 lo wires loaded before A, and table wires measured with their repair deferred to the unlookup.

What is taken from the public file (structure only, read in research/ws_invsqrt_decode.md 2, 2b, 4): stage A as
add-subtract rows on a low window of u with one pending wire that takes the carries out of the window; the top
of u left as table wires until the third lookup is done; the pending wire added into the top of u before B,
measured, and repaired in the mirrored rows; table wires measured early with the outcome kept for the unlookup.
The numbers, the rounding rule of the rows, the way the top of u is restored and every gate are ours.

Numerics (unchanged model, tables with tie_t0=1): u18 = C1 - sum_k sigma rr(c2, k), k = 0..7, bit 0 virtual.
  window  = bits 1..10 of u18 (10 wires, table = bits of C1) + pending wire P at weight 2^11
  top     = bits 11..17 (7 wires, table T = (C1 >> 11) - P0), P starts as P0 = 1 - bit 10 of C1 (one CNOT)
  safety  : sum_k rr(c2, k) <= 2^10, so the 12-bit window value never leaves [0, 2^12) [P, in the fit]
  merge   : top += P (6 ANDs), P measured; its phase is repaired when P is on a wire again (below)
  un-merge: after B^-1 the lo wires hold C0 again, and the fit ties bit 1 of C0 to bit 0 of T, so
            P = (top bit 0) xor (lo bit 1) by two CNOTs; top -= P (6 ANDs); P measured again, and this time its
            phase is repaired in the reverse rows: P = P0 xor (XOR of the carries out of the window), each of
            which is (carry into the top window wire) AND (that wire) in the mirrored row [P] -> one conditioned
            CZ per row, plus Z / NEG for P0.
  rows    : rows=0/1: our placeholder rows (invsqrt_gen_parts.row_win): every row ripples to the top of the
            window, 11 ANDs forward (into P), 10 in reverse. rows=2: the rows module's emit_A / emit_A_inv
            (same arithmetic, 60 forward, 65 in reverse with 21 clean wires).

Schedule: NORM, park, LUT1 [H, u, n_e1 lo, n_c1 c2; scratch], erase address, PH1, erase H, LUT1b [rest of the
n_e early lo wires and of c2], A (window rows), LUT3 [late lo wires; pool: early lo, top of u, c2], c2 measured,
MERGE, B, PH2, B^-1, UNMERGE, LUT2r [c2; pool: lo, top of u], lo measured, A^-1, ADDR [address from u, c2],
unLUT (all measured table wires repaired against the address and a chain of helper products), UNNORM.
With late=1 (phase module version 2): no H in LUT1 and no PH1 stage; ADDR also builds H; the hi stream runs
inside UNNORM.

Modules of other workstreams of our team, imported from pinned verbatim copies (see the comments at their
import): rows=1/2 (workstream invsqrt_rows: B, B^-1, the A rows), ph=1/2 (workstream invsqrt_ph: normaliser,
exponent code, streams, small-x path). Lookup programs: our own solver (invsqrt_gen_pool.py) or, per block, a
program file of the lookup workstream replayed through the same bookkeeping (lut/<block>.prog.json).

Builds published from this file (entry lines in research/ws_invsqrt_gen.md and in best/<tag>/ENTRY):
g6a, g6b (rows=1), g6c, g6d, g6e (rows=2 ph=1); they rebuild bit-identically from their entry lines.
"""
from __future__ import annotations

import hashlib
import json
import os
import random
import sys

import numpy as np

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_gen_model as M
import invsqrt_gen_parts as S
from invsqrt_gen_circ import Circ
from invsqrt_gen_code import ExpCode
from invsqrt_gen_pool import Pool, NotQuadratic, tt_bit, lookup_block, Span, bits_of, ONE
from invsqrt_gen_g3 import addr_tt


def digits(v):
    return [int(ch) for ch in str(v)]


def rows_b_inv(c, RB, t, u, lo, ctx, pool_first, pool_rest, after_first):
    """The reverse pass of the rows module (invsqrt_rows_b.emit_B_inv, same steps in the same order, built from
    its row and merge emitters) with the two things the integration needs: the rows of the last group (digits
    t4..t0, which run while the exponent code is off the t wires and four control wires are alive) get their own
    pool size, and `after_first()` is called when they are done (the code is parked there)."""
    nt, wlo = len(t), len(lo)
    lay, pend, ghosts = ctx["lay"], ctx["pend"], ctx["ghosts"]
    for gi in range(1, len(lay)):
        pw, tf, w = pend[gi]
        if gi == 1:
            RB._merge(c, pw, tf, lo[w:wlo - 1], lo[wlo - 1], True)
        else:
            pup, _, wup = pend[gi - 1]
            RB._merge(c, pw, tf, lo[w:wup], pup, True)
    m_pend = {}
    for gi in range(1, len(lay)):
        pw, tf, w = pend[gi]
        m_pend[gi] = c.hmr(pw)
        c.z(tf, cond=m_pend[gi])
        c.free(pw)
    qghost = ctx.get("qghost", {})
    for k, (bq, a) in qghost.items():
        c.z(t[nt - 1 - k], cond=bq)
    for gi in range(len(lay) - 1, -1, -1):
        grp, w, has_p = lay[gi]
        pool = pool_first if gi == len(lay) - 1 else pool_rest
        for k in reversed(grp):
            zfix = {}
            for k2 in grp:
                if k2 <= k and k2 in qghost:
                    zfix.setdefault(qghost[k2][1], []).append(qghost[k2][0])
            if has_p:
                RB._row_rev(c, lo[:w], None, u, k, t[nt - 1 - k], ghosts, m_pend[gi], pool, zfix)
            else:
                RB._row_rev(c, lo[:w], lo[w], u, k, t[nt - 1 - k], ghosts, None, pool, zfix)
        if has_p:
            c.bit_free(m_pend[gi])
        for k in grp:
            if k in qghost:
                c.bit_free(qghost.pop(k)[0])
        if gi == len(lay) - 1:
            after_first()
    assert not ghosts, "unrepaired forward outcomes"


def build_g6(P, tabs, angles=None, blk_b=7, blk_bn=5, blk_r=668, blk_rn=65432, n_e=13, n_e1=6, s1=8, s1b=0, s3=0,
             s2r=0, s_addr=12, blk_a=56, late_top=1, seed=1, deep=0, f1=0, f1b=0, f3=0, f2r=0, fa=0, rows=0,
             pool_n=7, pool_w=10, h7=0, unl=1, n_c1=0, ph=0, pool_a=0, c2_pre=0, qmax=62, late=0, hdrop=0,
             park=0, grp=0, lo_pre=0, c2_hi=0, useed=0, c2_pm=0, lo_pm=0, hlate=0, hsel=0, s1e=0, reuse=0,
             r2g=23324, pool_un=0, alazy=0, tin=0, rpin=1, akeep=0, akm=0, acut=0, atm=0, apre=0, apc=0,
             zcol=0, zh=0, ztp=0, f1z=2, orda=0, ordl=0, ztq=0, aon=0, lon=0, lonh=6543210, l1bh=0, l1bo=0,
             l1bw=65427, akp=0, m1b=0, l3h=0, l2h=0, ainv=0, angf=0):
    """rows=1: B and B^-1 from the rows module. h7=1: the top bit of the hi word is 1 for every address and gets
    no wire (as in the public file). unl=1: the unlookup uses a chain of products of dense forms over the
    address (about 15) instead of the 22 products of the two address nibbles. n_e1, n_c1: lo wires and c2 wires
    that LUT1 already builds next to H and u (the rest of the 13 early lo wires and of c2 comes in LUT1b).
    ph: phase-delivery module (pinned copy number). late=1 (needs ph >= 2): the hi stream runs inside UNNORM
    (the module's "late PH1"): the hi word is no output of LUT1; it is looked up together with the address at
    the end, streamed by pd.unnorm and erased in its hook against the address and the unlookup's helper wires,
    which also are the pool of the small-x special bits.
    Model option neg=1 (needs ph >= 3, late=1): negated-word tables, 40 rotations; the B rows then run on
    complemented digit wires (X gates around B and B^-1, taken off around every call into the phase module).
    park=1 (ph >= 3): the exponent code comes off the t wires only before the last group of the forward rows,
    with the control E6 parked on the wire of xn bit 6: three control wires next to the rows of t4..t0 instead
    of four. grp=1: the rows module's partition 2,3,3,6 instead of the public 3,3,3,5. lo_pre=K: K lo wires are
    measured before the block that reloads c2 (room for it on fewer wires)."""
    assert P.virt == 1 and P.tie_c1 == 1 and P.gu == 1 and P.fc2 == 1 and P.tie_t0 == 1, "g6 needs the window model"
    assert not P.neg or (ph >= 3 and late and rows >= 2), "the negated-word tables need phase module 3, late=1"
    assert not park or (ph >= 3 and rows >= 2)
    c = Circ(22)
    xq = list(range(22))
    dfun = lambda e: 10 - (e >> 1)
    nlo = P.wlo - 1
    nu = P.wu - 1
    win = P.win                          # window = bits 0..win-1 of u18, wires for bits 1..win-1
    kg = P.a_shift(P.a_lo)               # bit 0 of u18 is c2 bit kg before the first A row
    kmax = kg
    # ---------------- NORM
    c.stage("NORM")
    pd = code = None
    if ph:
        # Phase-delivery module of the workstream invsqrt_ph (normaliser, exponent code on t0..t3 + one flag,
        # both streams, small-x path: 41 rotations). The module is under development, so the build imports a
        # pinned verbatim copy of tools/invsqrt_ph_{mod,net,stream,sx}.py:
        # ph=1 -> tools/invsqrt_gen_pin/ph1/, files of 2026-11-11 12:17Z (sha256 a07dbfb2..., b3273f39...,
        #         f331c2d7..., 77997c42...): version 1, hi stream in the place of the public file;
        # ph=2 -> tools/invsqrt_gen_pin/ph2/, files of 12:42Z (sha256 c148df5f..., 9aeac6e6..., f331c2d7...,
        #         77997c42...): version 2 (tested by its author at 12:30Z), used here for the late hi stream.
        # ph=3 -> tools/invsqrt_gen_pin/ph3/, files of 13:27Z (sha256 56ceeb32..., 9aeac6e6..., 1a378bff...,
        #         77997c42...): version 3: negated-word form (PhaseDelivery(neg=True), 40 rotations) and
        #         code_off(park=True); credit for both ideas: workstream invsqrt_ph.
        pin = os.path.join(os.path.dirname(os.path.abspath(__file__)), "invsqrt_gen_pin", f"ph{ph}")
        for old in [k for k in sys.modules if k.startswith("invsqrt_ph_")]:
            del sys.modules[old]                 # one pinned copy per build, also inside a sweep
        sys.path[:] = [q for q in sys.path if "invsqrt_gen_pin" not in q]
        sys.path.insert(0, pin)
        from invsqrt_ph_mod import PhaseDelivery
        assert h7 == 1, "the phase module is used with the constant top bit of H"
        assert not late or ph >= 2
        if ph == 1:
            pd = PhaseDelivery(c, xq, angles, P.top_pos)
        elif ph == 2:
            pd = PhaseDelivery(c, xq, None, None, P.top_pos)
        elif ph == 3:
            pd = PhaseDelivery(c, xq, None, None, P.top_pos, neg=bool(P.neg))
        else:
            # ph=4 -> tools/invsqrt_gen_pin/ph4/, files of 2026-11-11 18:12Z (mod sha256 d81ee01a17ee1f93...,
            #         net 2c04b8342d14922d... = the network file its author tested at 17:12Z as "v4b"; the module
            #         file is a later state of the same version): forward network of 46 swaps (NORM 65) and
            #         UNNORM on the half thermometer (64.5 = 9 + 51 + 4.5). Same calls as ph=3.
            # ph=5 -> tools/invsqrt_gen_pin/ph5/, files of 19:18Z = "v4d" as tested by its author (mod sha256
            #         13b6ebd7f14bc348..., net f6dd1b9c083e7cb7..., taken together): backward network of 47
            #         swaps + 2 ANDs (UNNORM 62.5); the late hi stream holds 13 module wires instead of 15.
            #         The same files have the zone option (model option zone=1): zone flag and the rotation of
            #         the wires 0..14 at the end of NORM and at the start of UNNORM (+19 each), exponent code on
            #         xq[0..4], special bits for the odd x of zone 1 in 192..255.
            # ph=6 -> tools/invsqrt_gen_pin/ph6/, pinned 2026-11-12 02:30Z: the module file is the copy made by the
            #         workstream invsqrt_zpol (tools/invsqrt_zpol_phmod.py, sha256 2f912832dc9852d7...; v4d gate
            #         for gate at its defaults) with v4d's net / stream / sx files (f6dd1b9c..., 1a378bff...,
            #         77997c42...). Called with zsb="gate" (zone special bits with one shared gate product
            #         [e < 8] AND xq13: 1 + positions Toffolis instead of 2 per position) and sx_tries=2000 (more
            #         random multipliers in the small-x synthesis; never more products). Credit: invsqrt_zpol.
            # ph >= 6 in general: a pin folder ph<n>/ in this form (made by `tools/invsqrt_gen_hand.py pin`): the four
            #         module files and PIN.json with their sha256 and the keyword options that the constructor
            #         gets for this pin, so a further copy of the module with the same calls needs no edit here.
            #         A pinned file that is no longer the file that was pinned stops the build.
            opt6 = {}
            if ph >= 6:
                rec = json.load(open(os.path.join(pin, "PIN.json")))
                for fn, h in rec["sha256"].items():
                    assert hashlib.sha256(open(os.path.join(pin, fn), "rb").read()).hexdigest() == h,                         f"pinned file ph{ph}/{fn} is not the file that was pinned"
                opt6 = dict(rec["options"])
            ang_hi = ang_lo = None
            if angf:
                # angf=<k> (owner G2, 2026-11-12): rotation angles tuned by a linear programme for the largest
                # minimum margin over all inputs (tools/invsqrt_int_g2angles.py): integer offsets of 2^-64 turn on
                # the nominal -2^-q half turns, so every angle stays dyadic; Z, S, T (positions 0, 1, 2) stay
                # exact, the rotation count is unchanged. Pinned file tools/invsqrt_gen_pin/ang<k>.json, sha256
                # checked here. Valid for the circuit whose rotation bits it was tuned on (same tables, rows and
                # streams; lookup programs do not change those bits). ang1: tuned on g8x15, margin 56,710,366 ->
                # 75,502,023 units.
                from fractions import Fraction
                sha_a = {1: "3b0be9a43126266e4c4355b3484f416adf58ed335f373462b7e84ebccfdde45e"}[angf]
                f_a = os.path.join(os.path.dirname(os.path.abspath(__file__)), "invsqrt_gen_pin", f"ang{angf}.json")
                assert hashlib.sha256(open(f_a, "rb").read()).hexdigest() == sha_a, \
                    f"invsqrt_gen_pin/ang{angf}.json is not the file that was pinned"
                rec_a = json.load(open(f_a))
                assert P.neg and P.zone and ph >= 6
                ang_lo = {int(q): Fraction(n_, d_) for q, (n_, d_) in rec_a["lo"].items()}
                ang_hi = {q: Fraction(-1, 1 << q) for q in range(0, 18)}
                ang_hi.update({int(q): Fraction(n_, d_) for q, (n_, d_) in rec_a["hi"].items()})
                assert sorted(ang_lo) == list(range(11, P.top_pos + 1)) and all(int(q) >= 3 for q in rec_a["hi"])
            if P.zone:
                assert ph >= 5, "the zone format needs phase module 5"
                pd = PhaseDelivery(c, xq, ang_hi, ang_lo, P.top_pos, neg=bool(P.neg), half=True, zone=True, **opt6)
            else:
                pd = PhaseDelivery(c, xq, None, None, P.top_pos, neg=bool(P.neg), half=True, **opt6)
        p = pd.norm(keep_tree=True) if (ztp or ztq) else pd.norm()   # parity on wire 21; code on xq[0..4]; flag pd.F
        code_off, code_on = pd.code_off, pd.code_on
    else:
        f = S.thermo(c, xq)
        p = S.linear_of_e(c, f, lambda e: e & 1)
        db = [S.linear_of_e(c, f, lambda e, b=b: (dfun(e) >> b) & 1) for b in range(4)]
        fF = S.linear_of_e(c, f, lambda e: int(dfun(e) >= 1))
        fG = S.linear_of_e(c, f, lambda e: int(dfun(e) >= 2))
        S.thermo_erase(c, xq, f)
        S.shifter_refl(c, xq, p, db)
        c.x(xq[21])
        c.lend(xq[21])
        code = ExpCode(c, xq, fF, fG, db)
        code.park()
        code_off, code_on = (lambda: code.unpark(flags=False)), code.park
    pool = Pool(c, random.Random(seed))
    names, lutlog = {}, {}
    abits = [xq[14 + k] for k in range(7)] + [p]
    # Zone format (model option zone=1; format and tables: workstream invsqrt_tab). Only 198 of the 256 formal
    # addresses occur: inside zone 1 the wire of address bit 0 is zero (that bit went to the top digit wire).
    # Every function of the address is kept on all 256 formal addresses with the value of an address that does
    # not occur copied from the even address below it. The address wires and the tables have this property,
    # AND and XOR keep it, so two such functions are equal where it matters exactly if they are equal as
    # 256-bit vectors, and nothing else in the bookkeeping changes; only "every function is reached" means
    # rank 198 instead of 256.
    src_addr = [(a & ~1) if (P.zone and M.zone1(a >> 7, a & 127)) else a for a in range(256)]
    nused = len(set(src_addr))
    pool.n = nused
    pool.redundant = {a: s_ for a, s_ in enumerate(src_addr) if s_ != a}

    def att(b):
        """Truth table of address wire b."""
        return sum(1 << a for a in range(256) if (src_addr[a] >> b) & 1) if P.zone else addr_tt(b)

    if P.zone:
        assert nused == 198 and P.a_lo == 5 and late and rows == 3 and P.neg, "zone format: see the options of g8"
        for k_ in ("H", "C0", "C1", "C2"):
            assert all(tabs[k_][a] == tabs[k_][src_addr[a]] for a in range(256)), \
                f"zone tables: {k_} of an address that does not occur is no copy of the address below it"
    for b, w in enumerate(abits):
        pool.set(w, att(b))
        names[w] = f"a{b}"
    ztree, ztparts = {}, {}
    if ztp or ztq:
        # the forward zone tree handed over by the phase module (ph >= 9, norm(keep_tree=True)); a0 is not in it,
        # so its values are the same with the 198-class convention of the address
        av = {b: att(b) for b in range(8)}
        tv = {"t1": av[3] & av[2], "t2": ONE & ~av[5] & ~av[4]}
        tv["t3"] = tv["t2"] & ~tv["t1"]
        tv["t4"] = tv["t3"] & ~av[7]
        tv["z"] = av[6] & ~tv["t4"]
        for nm, v in tv.items():
            ztree[nm] = pd.ztree[nm]
            pool.set(ztree[nm], v)
            names[ztree[nm]] = f"zt_{nm}"
        ztparts = {"t1": [abits[3], abits[2]], "t2": [abits[5], abits[4]], "t3": [ztree["t2"], ztree["t1"]],
                   "t4": [ztree["t3"], abits[7]], "z": [abits[6], ztree["t4"]]}

    consumed = set()                     # input wires that a tied output was written onto (tin)

    def lut(block, outs, max_scr, empty=0, consume=(), needs=(), mid=None, onto=None):
        # onto = {output name: live pool wire} (owner B1): that output is made ON this wire, starting from what it
        # holds (pre-filled output of invsqrt_gen_pool.lookup_block); the wire is pool until it is finished
        # reuse: sum of 1 (LUT1), 2 (LUT1b), 4 (LUT3), 8 (LUT2r), 16 (ADDR): blocks built with wire reuse.
        # tin: same bits: tied outputs are left out of the lookup WITHOUT an extra scratch wire and written
        #      afterwards, in place onto an input wire of `consume` that the caller would release anyway
        #      (one wire less at the block's peak: the 59-wire layout), else onto a fresh wire.
        bit = 16 if block.startswith("ADDR") else {"LUT1": 1, "LUT1z": 1, "LUT1b": 2, "LUT3": 4, "LUT2r": 8}[block]
        if (reuse | tin) & bit or onto:
            return lut_reuse(block, outs, max_scr, empty, bool(reuse & bit), set(consume) if tin & bit else set(),
                             set(needs), mid, onto)
        assert mid is None and not onto
        wires = c.allocn(len(outs))
        for (n, _), w in zip(outs, wires):
            names[w] = n
        inputs = sorted(((names[w], w) for w in pool.val if w not in wires), key=lambda t: t[0])
        scr = lookup_block(c, pool, block, inputs, [(n, w, v) for (n, v), w in zip(outs, wires)], max_scr,
                           seed=seed, log=lutlog, deep=deep, empty=empty)
        for i, w in enumerate(scr):
            names[w] = f"{block}.s{i}"
        return wires, scr

    def lut_reuse(block, outs, max_scr, empty, extra, consume, needs, mid, onto_map=None):
        """"Wire reuse" (proposed by the architect workstream for the tie sets of the table-freedom
        workstream): an output that is an XOR of pool wires and of other outputs of the block is no output of
        the lookup at all; it is written by CNOTs after the block.
        extra (option reuse): its wire is one more scratch wire of the block; after the block the scratch
        wires beyond max_scr are released (erased if they are quadratic forms over the pool, measured
        otherwise). Same wire count, one more scratch wire per tie.
        not extra (option tin): no scratch wire is added. A tied output whose XOR contains an input wire of
        `consume` is written in place onto that wire (which the caller would have released: the block is one
        wire narrower); the others get fresh wires, those that need wires of `needs` before `mid()` is called
        (where the caller releases them), the rest after it."""
        onto_map = dict(onto_map or {})
        base, gens = Span(), [None]
        base.add(ONE, 1)
        for w in sorted(pool.val):
            if pool.val[w] and w not in onto_map.values():   # a pre-filled wire is no tie material
                gens.append(("w", w))
                base.add(pool.val[w], 1 << (len(gens) - 1))
        real, tied = [], []
        for n, v in outs:
            rem, tag = base.reduce(v)
            if rem == 0 and n not in onto_map:
                tied.append((n, v, [gens[i] for i in bits_of(tag)]))
                if os.environ.get("INVSQRT_INT_DEBUG"):
                    print(f"  [{block}] tie {n} = " + " ^ ".join("1" if g is None else (names.get(g[1], g[1]) if
                          g[0] == "w" else g[1]) for g in tied[-1][2]))
            else:
                gens.append(("o", n))
                base.add(v, 1 << (len(gens) - 1))
                real.append((n, v))
        # pre-filled outputs (onto) are finished last: until then they are pool
        real = [r for r in real if r[0] not in onto_map] + [r for r in real if r[0] in onto_map]
        fresh = iter(c.allocn(sum(1 for n, _ in real if n not in onto_map)))
        wires = [onto_map[n] if n in onto_map else next(fresh) for n, _ in real]
        for (n, _), w in zip(real, wires):
            names[w] = n
        inputs = sorted(((names[w], w) for w in pool.val if w not in wires), key=lambda t: t[0])
        scr = lookup_block(c, pool, block, inputs, [(n, w, v) for (n, v), w in zip(real, wires)],
                           max_scr + (len(tied) if extra else 0), seed=seed, log=lutlog, deep=deep, empty=empty)
        for i, w in enumerate(scr):
            names[w] = f"{block}.s{i}"
        consumed.update(onto_map.values())           # pre-filled hosts now hold outputs (not released as inputs)
        by_name = {n: w for (n, _), w in zip(real, wires)}

        def write(n, v, comb):
            ws = {g[1] if g[0] == "w" else by_name[g[1]] for g in comb if g is not None}
            const = sum(1 for g in comb if g is None) & 1
            onto = sorted(w for w in ws if w in consume and w not in consumed)
            if onto:                                 # in place: the output takes over an input wire
                w = onto[0]
                pool.xor_into(w, ws - {w}, const)
                consumed.add(w)
            else:
                w = c.alloc()
                pool.val[w] = 0
                pool.xor_into(w, ws, const)
            names[w] = n
            assert pool.val[w] == v, f"{block}: tied output {n} is wrong"
            by_name[n] = w

        if tied and extra:
            erase_or_defer(scr[max_scr:])            # the wires the tied outputs are written on now
            scr = scr[:max_scr]
        early = [t for t in tied if any(g is not None and g[0] == "w" and g[1] in needs for g in t[2])]
        for n, v, comb in early:
            write(n, v, comb)
        if mid:
            mid()
        for n, v, comb in tied:
            if n not in by_name:
                write(n, v, comb)
        if tied:
            lutlog[block]["tied"] = len(tied)
        return [by_name[n] for n, _ in outs], scr

    deferred = []                        # (outcome bit, truth table) of table wires measured early

    def defer(wires):
        for w in wires:
            m = c.hmr(w)
            deferred.append((m, pool.val[w]))
            pool.drop(w)
            names.pop(w, None)
            c.free(w)

    def erase_or_defer(wires):
        """Free erasure now if the wire is a quadratic form over the pool, else measure and repair at the end."""
        n = 0
        for w in wires:
            try:
                pool.erase(w)
            except NotQuadratic:
                defer([w])
                n += 1
                continue
            names.pop(w, None)
            c.free(w)
        c.meta.setdefault("deferred", []).append(n)

    tier_bits = [0, 1, 2, 3]
    cells = [[xq[0]], [xq[1], xq[2]], [xq[3], xq[4], xq[5], xq[6]], [xq[k] for k in range(3, 11)]]
    delays = list(range(11))

    C1 = tabs["C1"]
    T = [P.top_table(v)[1] for v in C1]
    assert min(T) >= 0 and max(T) < (1 << (P.wu - win))
    assert all((bin(tabs["C0"][a] & P.tie_mask).count("1") ^ T[a]) & 1 == 0 for a in range(256)), \
        "tables lack the tie of C0 to T"
    nh = P.wh
    if P.neg:                                    # the table holds N = 2^wh - Hc: wh - 1 wires, no constant bit
        assert h7 and max(tabs["H"]) < (1 << (P.wh - 1)) and min(tabs["H"]) >= 0
        nh = P.wh - 1
    elif h7:
        assert min(tabs["H"]) >= (1 << (P.wh - 1)), "the top bit of H is not constant"
        nh = P.wh - 1
    tH = [(f"H{h}", tt_bit(tabs["H"], h)) for h in range(nh)]
    tU = [(f"u{k}", tt_bit(C1, k)) for k in range(1, win)] + \
         [(f"u{win + b}", tt_bit(T, b)) for b in range(P.wu - win)]
    tLo = [(f"lo{k}", tt_bit(tabs["C0"], k)) for k in range(1, P.wlo)]
    tC2 = [(f"c2_{k}", tt_bit(tabs["C2"], k)) for k in range(P.wc2)]
    n_late = nlo - n_e
    if late_top:
        lo_early_t, lo_late_t = tLo[:n_e], tLo[n_e:]
    else:
        lo_late_t, lo_early_t = tLo[:n_late], tLo[n_late:]
    # ---------------- LUT1, address erased
    c.stage("LUT1")
    nhf = 0 if late else nh                      # hi wires built by LUT1 (none with the late hi stream)
    # c2 bits built by LUT1 / by LUT1b. c2_hi=1: LUT1 takes the TOP n_c1 bits (the top bits of c2 are smooth
    # functions of the address, which the lookup workstream's design A builds for 1 to 3 ANDs each from the
    # bare address and nobody builds cheaply later).
    c2_first = tC2[P.wc2 - n_c1:] if c2_hi else tC2[:n_c1]
    c2_rest = tC2[:P.wc2 - n_c1] if c2_hi else tC2[n_c1:]
    outs1 = tH[:nhf] + tU + lo_early_t[:n_e1] + c2_first

    def lut_ord(block, outs, rev, last, *a, **kw):
        """rev=1 (owner B1): the block synthesises its outputs in reversed order (the ones named in `last` at the
        end, so that a tied output still comes after the outputs its XOR uses); wires returned in the order of
        `outs`."""
        if not rev:
            return lut(block, outs, *a, **kw)
        o2 = [o for o in reversed(outs) if o[0] not in last] + [o for o in outs if o[0] in last]
        w2, s2 = lut(block, o2, *a, **kw)
        byn = {n: w for (n, _), w in zip(o2, w2)}
        return [byn[n] for n, _ in outs], s2

    def onto_last(outs, k, hosts, skip):
        """The last k real outputs of `outs` (names in `skip` are ties written by the builder) made in place on
        the wires `hosts` (pool wires the caller releases right after the block)."""
        real_n = [n for n, _ in outs if n not in skip]
        assert k <= len(hosts), "not enough host wires"
        return dict(zip(real_n[len(real_n) - k:], hosts[:k]))

    w1b_merged = None
    if ztq:
        # ztq=1 (owner B1): ONE block; the forward zone tree's five wires are LUT1's last five real output wires,
        # pre-filled: pool for every earlier output, then finished from what they hold (no split, no erasure)
        real_names = [n for n, _ in outs1 if n != "c2_6"]          # c2_6: the table set's tie, written in place
        onto1 = dict(zip(real_names[-5:], [ztree[k] for k in ("t1", "t2", "t3", "t4", "z")]))
        if m1b:
            # m1b=1 (owner B1): LUT1b's outputs (lo18..lo20, c2_8) join LUT1 as its very last outputs, made on
            # address wires (order l1bw; erased after LUT1 anyway): no LUT1b block
            outs_m = lo_early_t[n_e1:] + c2_rest
            hosts_m = [w for w in (abits[b] for b in digits(l1bw))][:len(outs_m)]
            onto1.update(dict(zip([n for n, _ in outs_m], hosts_m)))
            w_all, scr = lut_ord("LUT1", outs1 + outs_m, ordl, ("c2_6",), s1, f1, consume=abits, onto=onto1)
            w1, w1b_merged = list(w_all[:len(outs1)]), list(w_all[len(outs1):])
        else:
            w1, scr = lut("LUT1", outs1, s1, f1, consume=abits, onto=onto1)
    elif lon:
        # lon=K (owner B1): LUT1's last K real outputs are made in place on address wires (pre-filled with a_k;
        # erased after LUT1 anyway), so K fresh wires more are free during the block: scratch for free fills
        onto1 = onto_last(outs1, lon, [abits[b] for b in digits(lonh)], ("c2_6",))   # lonh: host order
        w1, scr = lut("LUT1", outs1, s1 + lon, f1, consume=abits, onto=onto1)
    elif not ztp:
        w1, scr = lut_ord("LUT1", outs1, ordl, ("c2_6",), s1, f1, consume=abits)
    else:
        # ztp=K (owner B1; the pool mechanism of DESIGN3 item I1 where the wires exist): the forward zone flag's
        # AND tree (t1..t4, z: five functions of the address that NORM makes anyway) is pool of LUT1's first K
        # outputs; then it is erased against its own factors (free) and LUT1z makes the rest on its wires.
        k1 = ztp
        w1a, scra = lut("LUT1", outs1[:k1], s1, f1, consume=abits)
        for nm in ("z", "t4", "t3", "t2", "t1"):
            w = ztree[nm]
            pool.erase(w, others=ztparts[nm])
            names.pop(w, None)
            c.free(w)
        w1b, scrb = lut("LUT1z", outs1[k1:], 0, f1z, consume=abits)
        w1, scr = list(w1a) + list(w1b), list(scra) + list(scrb)

    def small_x(pool_wires):
        """Data of the small-x path for the phase module: for x < 64 the lo bits above position 10 go out through
        the hi stream. The module needs the hi and lo words of these 63 inputs (functions of the address: all
        their t bits are 0) and wires that stay unchanged during the hi stream with their values on them."""
        if P.zone:
            # Zone format (phase module's interface, research/ws_invsqrt_ph.md, section Z9c): the words of
            # x = 1..255 from the zone model (for the odd x of zone 1 in 128..255 the digit word is 2^13, not 0,
            # and their address bit 0 sits alone on the wire of digit 13, which must be in the pool).
            xs = np.arange(1, 256, dtype=np.int64)
            ev = M.evaluate(P, tabs, xs)
            a8s, ts = ev["a8"], ev["t"]
            assert not (ts & 0x1FFF).any()
            return dict(H={int(x): int(tabs["H"][int(a)]) for x, a in zip(xs, a8s)},
                        lo={int(x): int(v) >> 1 for x, v in zip(xs, ev["lo"])},
                        pool=[(w, sum(((pool.val[w] >> int(a)) & 1) << int(x) for x, a in zip(xs, a8s)))
                              for w in pool_wires]
                        + [(xq[13], sum(((int(t_) >> 13) & 1) << int(x) for x, t_ in zip(xs, ts)))])
        xs = np.arange(1, 64, dtype=np.int64)
        _, _, _, _, a8s, ts = M.decompose(xs)
        assert not ts.any()
        ev = M.evaluate(P, tabs, xs)
        return dict(H={int(x): int(tabs["H"][int(a)]) for x, a in zip(xs, a8s)},
                    lo={int(x): int(v) >> 1 for x, v in zip(xs, ev["lo"])},
                    pool=[(w, sum(((pool.val[w] >> int(a)) & 1) << int(x) for x, a in zip(xs, a8s)))
                          for w in pool_wires])
    # l1bh=K (owner B1): K address wires (order l1bw, digits = address bit, 7 = parity) stay alive through LUT1b
    # as its pool; the first l1bo of them carry LUT1b's last outputs (pre-filled); the rest are erased after it
    keep_a = [w for w in (abits[b] for b in digits(l1bw)) if w not in consumed][:l1bh] if l1bh else []
    for w in abits[:7]:
        if w in consumed or w in keep_a:
            continue                             # tin: a tied LUT1 output was written onto this address wire
        pool.erase(w)
        names.pop(w)
        c.lend(w)
    if p not in consumed and p not in keep_a:
        pool.erase(p)
        names.pop(p)
        if p < 22:
            c.lend(p)                            # phase module: the parity sat on the wire of bit 21
        else:
            c.free(p)
    for cw in [w for w in abits if w in consumed]:
        nw = c.alloc()                           # the table wire leaves the input wire (one SWAP, free)
        c.swap(cw, nw)
        pool.val[nw] = pool.val.pop(cw)
        names[nw] = names.pop(cw)
        if cw < 22:
            c.lend(cw)
        else:
            c.free(cw)
        w1 = [nw if w == cw else w for w in w1]
        if w1b_merged is not None:
            w1b_merged = [nw if w == cw else w for w in w1b_merged]
    H, u = w1[:nhf], w1[nhf:nhf + nu]
    lo_e = w1[nhf + nu:nhf + nu + n_e1]
    c2 = w1[nhf + nu + n_e1:]
    uw, top = u[:win - 1], u[win - 1:]
    # ---------------- PH1
    c.meta["checks"] = []
    if not late:
        c.stage("PH1")
        c.meta["checks"] = [("PH1", "H", H, 0)]
        if ph:
            pd.ph1(H + [None] * (P.wh - nh), small_x([w for w in sorted(pool.val) if w not in H]))
        else:
            code.unpark()
            S.stream_d(c, H + [None] * (P.wh - nh), tier_bits, code.b, cells, delays, 7 + 10, P.top_pos, angles)
            code.park()
        erase_or_defer(H)
    # ---------------- rest of the early lo wires and c2
    c.stage("LUT1b")
    if s1e:                                      # LUT1's scratch wires are released before LUT1b instead of after
        erase_or_defer(scr)                      # it: the spec of LUT1b then does not depend on what the LUT1
        scr = []                                 # program leaves on them, and LUT1b gets the wires as its own
    outs_b = lo_early_t[n_e1:] + c2_rest
    onto_b = onto_last(outs_b, l1bo, keep_a, ()) if l1bo else None
    if w1b_merged is not None:                   # m1b: LUT1 made these outputs already
        assert not l1bh
        w2, scr1b = w1b_merged, []
    else:
        w2, scr1b = lut("LUT1b", outs_b, s1b, f1b, onto=onto_b)
    if keep_a:
        for w in keep_a:                         # address wires kept for LUT1b and not hosting an output: erased
            if w in consumed:
                continue
            pool.erase(w)
            names.pop(w)
            c.lend(w) if w < 22 else c.free(w)
        erase_or_defer(scr + scr1b)              # scratch first, so that the moves below add no wire
        scr, scr1b = [], []
        for cw in [w for w in keep_a if w in consumed]:
            nw = c.alloc()                       # a LUT1b output on an input wire moves off it (one SWAP, free)
            c.swap(cw, nw)
            pool.val[nw] = pool.val.pop(cw)
            names[nw] = names.pop(cw)
            c.lend(cw) if cw < 22 else c.free(cw)
            w2 = [nw if w == cw else w for w in w2]
    lo_e = lo_e + w2[:n_e - n_e1]
    c2 = (w2[n_e - n_e1:] + c2) if c2_hi else (c2 + w2[n_e - n_e1:])      # c2[k] = wire of bit k
    erase_or_defer(scr + scr1b)
    # ---------------- A: rows on the window, pending wire on top
    c.stage("A")
    for w in uw:
        pool.drop(w)

    def make_pw():
        w = c.alloc()
        c.cx(uw[-1], w)
        c.x(w)                                   # P0 = 1 - bit (win-1) of C1
        return w

    pw = None if (rows == 3 and alazy) else make_pw()
    tw = [xq[i] for i in range(P.nt)]
    RB = R2 = None
    if rows == 3:
        # Second rows module (workstream invsqrt_rows2, tools/invsqrt_rows2_b.py: another partition of the B
        # rows with four pending wires; B 195, B^-1 219 at pools 8 / 6 where the first module has 232.3; same
        # arithmetic). Pinned verbatim copies: invsqrt_gen_rows2p1.py = its file v1d of 2026-11-11 17:46Z,
        # sha256 49528455643a9d06c12a75f7f3e11c9180fe3681302ea7127c3511207f34c451, and the first module's file
        # it takes its row primitives from, invsqrt_gen_rowsb3.py = tools/invsqrt_rows_b.py of 12:35Z,
        # sha256 4a063c2c58de5917a6dcf7fddd98d97949158ec1b07c75e46353cc0de0d0d6e9. The copy binds itself to a
        # module named invsqrt_rows_b when it is imported; that name is given to the pinned copy first.
        import importlib
        RB = importlib.import_module("invsqrt_gen_rowsb3")
        sys.modules["invsqrt_rows_b"] = RB
        # rpin=<k>: a later pinned copy invsqrt_gen_rows2p<k>.py of the same module with the same calls (made by
        # `tools/invsqrt_gen_hand.py rpin`; its sha256 in invsqrt_gen_pin/rows2p<k>.json is checked here).
        if rpin != 1:
            rec = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "invsqrt_gen_pin",
                                              f"rows2p{rpin}.json")))
            here = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"invsqrt_gen_rows2p{rpin}.py")
            assert hashlib.sha256(open(here, "rb").read()).hexdigest() == rec["sha256"],                 f"invsqrt_gen_rows2p{rpin}.py is not the file that was pinned"
        R2 = importlib.import_module(f"invsqrt_gen_rows2p{rpin}")
        R2.bind("invsqrt_gen_rowsb3")
    elif rows:
        # Rows module of the workstream invsqrt_rows (tools/invsqrt_rows_b.py). The build imports a pinned
        # verbatim copy, so that published circuits rebuild bit-identically while the module moves on:
        # rows=1: invsqrt_gen_rowsb1.py = the module as of 2026-11-11 11:45Z, sha256 188310f9...24fa. B / B^-1 only,
        #         flat rows (B 206, B^-1 218.5), reverse loop copied here (rows_b_inv); stage A = parts.row_win.
        # rows=2: invsqrt_gen_rowsb2.py = the module's final file of its first session, 12:08Z,
        #         sha256 13e029fe419971b042d5701d57a9b99a5ffef06e3e628fb35d4f4318a0123c2e. B with its tree device
        #         (195), B^-1 through its own emit_B_inv with hooks and the tail variant, stage A rows through
        #         emit_A / emit_A_inv (60 forward, 63 to 80 in reverse depending on the clean wires).
        import importlib
        RB = importlib.import_module(f"invsqrt_gen_rowsb{rows}")
    pa = {}
    if rows == 3:
        # alazy=1: the pending wire is made by the module when the rows are done (one clean wire less in A)
        actx = R2.emit_A(c, tw, uw, make_pw if alazy else pw, c2, nrows=kmax + 1, kg=kg)
        pw = actx["pend"] if alazy else pw
    elif rows >= 2:
        actx = RB.emit_A(c, tw, uw, pw, c2, nrows=kmax + 1, kg=kg)
    else:
        for k in range(0, kmax + 1):
            pa[k] = S.row_win(c, uw + [pw], c2, k, kg, xq[P.nt - 1 - k], 1, digits(blk_a), "fwd")
    # ---------------- LUT3: late lo wires (pool: early lo, top of u, c2); c2 measured
    c.stage("LUT3")
    onto3 = None
    if l3h:
        # l3h=1 (owner B1): the c2 wires that c2_pre would measure here host LUT3's last real outputs instead
        # (pre-filled): pool until those outputs are finished, no measurement
        assert c2_pre > 0 and tin & 4
        onto3 = onto_last(lo_late_t, c2_pre, list(c2[:c2_pre]), ("lo2", "lo5"))
        c2 = c2[c2_pre:]
    elif c2_pre > 0:                             # room for the late lo wires on one wire less (61-wire builds)
        defer(c2[:c2_pre])
        c2 = c2[c2_pre:]
    elif c2_pre < 0:                             # the same with the top wires of c2 (smooth functions of the
        defer(c2[c2_pre:])                       # address: the weakest members of the pool of this lookup)
        c2 = c2[:c2_pre]
    if c2_pm:                                    # ... or with exactly the c2 bits of this mask (for tie sets
        assert not c2_pre                        # that need other c2 columns in the pool of LUT3)
        gone = [w for k, w in enumerate(c2) if (c2_pm >> k) & 1]
        defer(gone)
        c2 = [w for w in c2 if w not in gone]
    if tin & 4:
        # tied late lo wires get no wire during the lookup: those whose XOR needs c2 wires are written right
        # after it, then c2 is released, then the others are written (the stage is one wire narrower)
        c2_now = list(c2)
        lo_l, scr3 = lut("LUT3", lo_late_t, s3, f3, needs=c2_now, mid=lambda: erase_or_defer(c2_now), onto=onto3)
        c2 = []
    else:
        lo_l, scr3 = lut("LUT3", lo_late_t, s3, f3)
    lo = (lo_e + lo_l) if late_top else (lo_l + lo_e)
    erase_or_defer(list(scr3) + list(c2))
    for w in lo + top:
        pool.drop(w)
    c.stage("MERGE")
    S.add_into(c, top, [], pw)                   # top += P
    m_p = c.hmr(pw)                              # repaired in UNMERGE
    c.free(pw)
    # ---------------- B, PH2, B^-1
    c.stage("B")
    c.meta["checks"] += [("B", "u", u, 1), ("PH2", "lo", lo, 1)]
    if rows:
        # rows module: same arithmetic as row_virt on windows that shrink in groups, three pending wires across
        # PH2. Its forward pass holds at most 6 clean wires (pending wires included), so the whole of it runs
        # with the code off (four control wires: 14 + 4 + 37 + 6 = 61); its reverse pass runs the rows of the
        # last group (t4..t0) first, with pool_n clean wires, then the code is parked and the other groups get
        # pool_w.
        assert P.b_shift(P.nt - 1) == 0 and len(lo) == 20 and len(u) == 17
        if rows == 3:
            # r2g: the group sizes as digits, 23324 = the module's G60 (four pending wires, 6 clean wires
            # forward), 2444 = its G59 (three pending wires, 5 clean wires forward)
            groups = R2.groups_from_sizes(tuple(digits(r2g)))
            # zone format: the exponent code sits on the wires of digits 0..4 (not 0..3), so the rows of
            # t4..t0 must all be in the last group, the only one that runs with the code off (2435, 12335)
            assert not P.zone or (park and digits(r2g)[-1] >= 5), "zone format: last group = rows of t4..t0"
            emit_B = lambda **kw: R2.emit_B(c, tw, u, lo, groups=groups, hook=kw.get("hook"))
            emit_B_inv = lambda ctx, pool, hook: R2.emit_B_inv(c, tw, u, lo, ctx, pool=pool, hook=hook, tail=True,
                                                               pool_un=pool_un or None)
        else:
            groups = RB.BEST_GROUPS if grp else RB.PUBLIC_GROUPS
            emit_B = lambda **kw: RB.emit_B(c, tw, u, lo, groups=groups, **kw)
            emit_B_inv = lambda ctx, pool, hook: RB.emit_B_inv(c, tw, u, lo, ctx, pool=pool, hook=hook, tail=True)
        last = len(groups) - 1

        def flip_t():
            """Negated-word tables: the rows see every digit complemented, so they add where they subtract
            (lo = C0 + sum sigma rr(u, k)). The phase module reads the t wires as they are: the complement is
            taken off around each of its calls. X gates only."""
            if P.neg:
                for w in tw:
                    c.x(w)

        if park:
            def f_hook(event, gi):
                if event == "group_start" and gi == last:
                    flip_t()
                    pd.code_off(park=True)       # three control wires; E6 sits on the wire of xn bit 6
                    flip_t()

            flip_t()
            bctx = emit_B(tree=True, hook=f_hook)
        else:
            code_off()
            flip_t()
            bctx = emit_B(tree=(rows >= 2))
        flip_t()
        c.stage("PH2")
        if ph:
            pd.ph2(lo)
        else:
            S.stream_d(c, lo, tier_bits, code.b, cells, delays, 25 + P.glo + 10 - 1, P.top_pos, angles)
        c.stage("B^-1")
        if rows >= 2:
            def b_hook(event, gi):
                if event == "group_end" and gi == last:
                    flip_t()
                    code_on()
                    flip_t()

            flip_t()
            emit_B_inv(bctx, {"default": pool_w, last: pool_n}, b_hook)
            flip_t()
        else:
            rows_b_inv(c, RB, tw, u, lo, bctx, pool_n, pool_w, code_on)
    else:
        pb = {}
        for i in P.b_rows:
            if i == 4:
                code_off()                       # rows of t4..t0 and PH2 need t0..t4 clean
            pb[i] = S.row_virt(c, lo, u, P.b_shift(i), xq[i], 1, blk_b if i > 4 else blk_bn, "fwd")
        c.stage("PH2")
        if ph:
            pd.ph2(lo)
        else:
            S.stream_d(c, lo, tier_bits, code.b, cells, delays, 25 + P.glo + 10 - 1, P.top_pos, angles)
        c.stage("B^-1")
        for i in reversed(P.b_rows):
            S.row_virt(c, lo, u, P.b_shift(i), xq[i], 0, digits(blk_rn) if i <= 4 else digits(blk_r), "rev", pb[i])
            if i == 4:
                code_on()
    for k, w in enumerate(lo):
        pool.set(w, tLo[k][1])
    # ---------------- top of u back to the table
    c.stage("UNMERGE")
    pw = c.alloc()
    c.cx(top[0], pw)
    for bit in range(1, P.wlo):                  # P = (T0 xor P) xor (tied C0 bits = T0)
        if (P.tie_mask >> bit) & 1:
            c.cx(lo[bit - 1], pw)
    c.z(pw, cond=m_p)
    c.bit_free(m_p)
    for w in top:
        c.x(w)
    S.add_into(c, top, [], pw)                   # top -= P
    for w in top:
        c.x(w)
    m_x = c.hmr(pw)                              # repaired by the reverse rows of A
    c.free(pw)
    for b, w in enumerate(top):
        pool.set(w, tU[win - 1 + b][1])
    # ---------------- c2 again; lo measured
    c.stage("LUT2r")
    lo_first = lo[:lo_pre] if lo_pre > 0 else (lo[lo_pre:] if lo_pre < 0 else [])
    if lo_pm:                                    # exactly the lo wires of this mask (bit k = wire lo_k, k = 1..20)
        assert not lo_pre
        lo_first = [w for k, w in enumerate(lo) if (lo_pm >> (k + 1)) & 1]
    if lo_first:                                 # 60-wire builds: this block is 15 + 10 + 27 inputs + 9 outputs wide,
        defer(lo_first)                          # so lo wires are measured before it instead of after it (K > 0:
    had = set(consumed)
    onto2 = None
    if l2h:
        # l2h=K (owner B1): LUT2r's last K real outputs are made on lo wires (pre-filled; measured after LUT2r
        # anyway) outside the tie of c2_8; K more scratch wires for free fills
        lon_ = {f"lo{k + 1}": w for k, w in enumerate(lo)}
        hosts2 = [lon_[n] for n in ("lo3", "lo4", "lo5", "lo7", "lo8", "lo9", "lo10", "lo12", "lo14", "lo15",
                                    "lo16", "lo20") if n in lon_ and lon_[n] not in lo_first]
        onto2 = onto_last(tC2, l2h, hosts2, ("c2_8",))
    c2, scr2r = lut("LUT2r", tC2, s2r + l2h, f2r,      # the lowest K, K < 0: the highest -K)
                    consume=[w for w in lo if w not in lo_first], onto=onto2)
    erase_or_defer(list(scr2r))
    lo_left = [w for w in lo if w not in lo_first and w not in consumed - had]
    # apre=n (owner B1, item I1): the first n outputs of ADDR (a0, a1, ...) are looked up HERE, before A^-1, on the
    # pool top of u + lo + c2 (34 wires), each on the wire of a lo wire measured just before (highest lo first,
    # the kept ones of akm excepted). apc=1: one sub-block per output (pool constant), else one block of n.
    # They stay alive through A^-1 (which then has n clean wires fewer) and are pool of the rest of ADDR.
    a_all = [(f"a{b}", att(b)) for b in range(8)] + (tH if late else [])
    if zcol:
        # zcol=1 (owner B1, DESIGN3 item I2, S3-INV-E 1.3): the zone flag z is one more (last) output of ADDR; the
        # phase module (ph >= 7, tools/invsqrt_int_phmod.py) then uses this wire for the backward zone rotation
        # instead of making z with 5 ANDs. z stays alive through the hi stream and is measured after the rotation,
        # its phase repaired with the helper chain, which is erased only then.
        assert P.zone and late and hlate and ph >= 7, "zcol needs the zone format, the late hi stream and ph >= 7"
        raise NotImplementedError("zcol: not finished; after the backward rotation wire xq[14] holds the true a0 "
                                  "in zone 1, so z cannot be repaired with the unlookup chain there. Use zh=1.")
        ztt = sum(1 << a for a in range(256) if M.zone1(a >> 7, a & 127))
        zf = sum(1 << a for a in range(256)
                 if (a >> 6) & 1 and ((a >> 7) & 1 or (a >> 5) & 1 or (a >> 4) & 1 or ((a >> 3) & 1 and (a >> 2) & 1)))
        assert ztt == zf, "zcol: the zone model and the module's flag formula differ"
        a_all = a_all + [("z", ztt)]
    pre_w = []
    if apre:
        keep_m = [w for k, w in enumerate(lo) if (akm >> (k + 1)) & 1]
        drop = [w for w in lo_left[::-1] if w not in keep_m][:apre]
        if apc:
            for g in range(apre):
                defer([drop[g]])
                w_g, _ = lut(f"ADDRe{g}", a_all[g:g + 1], 0, 0, consume=[])
                pre_w += list(w_g)
        else:
            defer(drop)
            w_g, _ = lut("ADDRe", a_all[:apre], 0, min(fa, apre - 1), consume=[])
            pre_w += list(w_g)
        lo_left = [w for w in lo_left if w not in drop]      # (a measured wire's id may now hold an output)
    # akeep=K (owner B1 since 2026-11-12, DESIGN3 item I1, design S3-INV-T 1.2): K of these lo wires stay alive
    # through A^-1 as pool of the address lookup (they hold bits of C0, functions of the address); A^-1 then runs
    # with K clean wires fewer. Which ones: akm = mask (bit k = wire lo_k), else the highest K.
    lo_kept = []
    if akeep or akm:
        if akm:
            lo_kept = [w for k, w in enumerate(lo) if (akm >> (k + 1)) & 1 and w in lo_left]
        else:
            lo_kept = lo_left[::-1][:akeep]
        assert len(lo_kept) == (bin(akm).count("1") if akm else akeep), "akeep / akm: wire not available"
    defer([w for w in lo_left if w not in lo_kept])
    if atm:
        # atm = mask (bit b = wire u_(win+b), the top of u holding T): these wires are measured here, before A^-1,
        # so that A^-1 keeps its clean wires while lo wires stay alive in their place (pool trade for ADDR)
        gone = [w for b, w in enumerate(top) if (atm >> b) & 1]
        defer(gone)
        u = [w for w in u if w not in gone]
    # ---------------- A^-1: mirrored rows on the window
    c.stage("A^-1")
    if rows == 3 and ainv:
        # ainv=<k> (owner G2, 2026-11-12): A^-1 by the pinned copy tools/invsqrt_gen_ainv<k>.py of G2's
        # tools/invsqrt_int_ainv.py: the same joint blocks and repairs, but a block's top chain may run with fewer
        # carry wires (bounded blocks, boundary repaired under its outcome), so the planner affords the joint
        # block of rows 5..8 at 17 clean wires (A^-1 71 -> 70). Copies are never edited; the sha256 is checked.
        import importlib
        # 2: the same plus the carry-in of a joint block's widest row shed and re-made under the outcome of the
        #    carry above it, with erasures that cancel AND terms by live wires (A^-1 70 -> 69.5 at 17 wires)
        sha_k = {1: "8b86938eed3e4e1cdc24152583784ea2755aba3922d1d960a035de7bb7d57639",
                 2: "b24d1d7b9c667ca47339e3830cd6b3c85c4534d2388e5ee7028b14a4a575e47b"}[ainv]
        here = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"invsqrt_gen_ainv{ainv}.py")
        assert hashlib.sha256(open(here, "rb").read()).hexdigest() == sha_k, \
            f"invsqrt_gen_ainv{ainv}.py is not the file that was pinned"
        AI = importlib.import_module(f"invsqrt_gen_ainv{ainv}")
        AI.bind("invsqrt_gen_rowsb3")
        AI.emit_A_inv_ext(c, tw, uw, c2, actx, m_x, pool=pool_a if pool_a else qmax - len(c.live))
    elif rows == 3:
        # the second module's reverse stage A builds its rows from the c2 wires passed here
        R2.emit_A_inv(c, tw, uw, c2, actx, m_x, pool=pool_a if pool_a else qmax - len(c.live))
    elif rows >= 2:
        # c2 sits on other wires than in the forward pass (it was measured after LUT3 and rebuilt by LUT2r), and
        # the pinned emit_A_inv takes its row descriptions (with the wire lists) from the forward context: they
        # are rebuilt here for the wires of now. The outcome bits in the context are keyed by row and bit.
        actx["specs"] = [RB._spec_A(c2, k, tw[P.nt - 1 - k], kg) for k in range(kmax + 1)]
        # the module may hold every wire that is free here: narrow rows are then undone jointly
        RB.emit_A_inv(c, tw, uw, c2, actx, m_x, pool=pool_a if pool_a else qmax - len(c.live))
    else:
        for k in range(kmax, -1, -1):
            S.row_win(c, uw, c2, k, kg, xq[P.nt - 1 - k], 0, [len(uw)], "rev", pa[k], [m_x])
    c.z(uw[-1], cond=m_x)                        # P0 = 1 - bit (win-1) of C1
    c.neg(cond=m_x)
    c.bit_free(m_x)
    for k, w in enumerate(uw):
        pool.set(w, tU[k][1])
    # ---------------- address back
    c.stage("ADDR")
    sep = {tuple((pool.val[w] >> a) & 1 for w in list(u) + list(c2) + lo_kept) for a in range(256)}
    assert len(sep) == nused, "u and c2 do not tell the addresses apart"
    had = set(consumed)
    a_outs = a_all[len(pre_w):]                  # apre: the first outputs exist already (before A^-1)
    if aon or akp:
        # aon=K (owner B1, item I1): ADDR's last K real outputs are made in place on c2 / u wires that the tie of
        # a3 does not use (measured after ADDR anyway): K fresh wires more for scratch and free fills.
        # akp=1: the lo wires kept through A^-1 (akeep / akm) host the last outputs (pre-filled) instead of
        # being measured after ADDR: the chain of S3-INV-T 1.2 without a measured pool wire
        # atm with akp (G2, 2026-11-12): allowed when aon = 0, because then only the kept lo wires host outputs
        # (the u names below shift when atm has removed u wires, so the u hosts of aon would be wrong)
        assert not (atm and aon), "atm with aon: the u host names would shift"
        cw = {n: w for n, w in zip([f"c2_{k}" for k in range(P.wc2)], c2)}
        uw_ = {f"u{k + 1}": w for k, w in enumerate(u)}
        hosts = [cw[n] for n in ("c2_0", "c2_2", "c2_7", "c2_8") if n in cw] + \
                [uw_[n] for n in ("u2", "u3", "u4", "u6", "u9", "u12", "u14", "u16") if n in uw_]
        hosts = (list(lo_kept) if akp else []) + [w for w in hosts if w not in consumed]
        ontoa = onto_last(a_outs, aon + (len(lo_kept) if akp else 0), hosts, ("a3",))
        wa, scr_a = lut("ADDR", a_outs, s_addr + aon, fa, consume=list(u) + list(c2), onto=ontoa)
        scr_a = list(scr_a) + [w for w in lo_kept if w not in consumed]
    elif not acut:
        wa, scr_a = lut_ord("ADDR", a_outs, orda, ("a3",), s_addr, fa, consume=list(u) + list(c2))
        scr_a = list(scr_a) + lo_kept            # kept lo wires (akeep): measured with the table wires below
    else:
        # acut = digits n0 n1 ... (owner B1, item I1): ADDR as a chain of sub-blocks ADDRp0, ADDRp1, ... with n_g
        # outputs each, in the order a0..a7, H0.. (the tied output counted). Before a sub-block, kept lo wires are
        # measured, one for each output that has no clean wire; the last sub-block first measures every kept lo
        # wire left and takes the clean wires beyond its outputs as scratch (at most s_addr); earlier sub-blocks
        # have none. Each sub-block is a lookup block of its own (spec and program files ADDRp<g>.*).
        cuts = digits(acut)
        assert sum(cuts) == len(a_outs), "acut: the digits must add up to the outputs of ADDR"
        wa, scr_a, left, pos = [], [], list(lo_kept), 0
        for g, n in enumerate(cuts):
            outs_g = a_outs[pos:pos + n]
            pos += n
            last = g == len(cuts) - 1
            take = len(left) if last else max(0, n - (qmax - len(c.live)))
            defer(left[:take])
            left = left[take:]
            scr_g = max(0, min(s_addr, qmax - len(c.live) - n)) if last else 0
            w_g, s_g = lut(f"ADDRp{g}", outs_g, scr_g, min(fa, max(0, n + scr_g - 1)), consume=list(u) + list(c2))
            wa += list(w_g)
            scr_a += list(s_g)
        assert not left
    wa = pre_w + list(wa)
    u = [w for w in u if w not in consumed - had]     # tin: a tied address or hi bit took over one of these wires
    c2 = [w for w in c2 if w not in consumed - had]
    abits = wa[:8]
    p = wa[7]
    if late:
        H = wa[8:8 + len(tH)]                    # the hi word, looked up here for the late hi stream
    zw = wa[8 + len(tH)] if zcol else None       # zcol: the zone flag's wire
    # ---------------- unLUT: everything that is left is measured, then repaired against the address
    c.stage("unLUT")
    kept = []
    if hlate:
        # hlate=K: the helper chain is built only after the hi stream (in the phase module's hook). K table
        # wires (u first, then c2) stay alive through the stream as the pool of the small-x special bits, which
        # the helpers were before; they are measured in the hook like every other table wire. No helper has to
        # be erased and made again around the stream (hdrop), which is where this saves its two ANDs.
        assert late and unl and not hdrop
        order = [list(u) + list(c2), list(c2) + list(u), list(u)[::-1] + list(c2), list(c2)[::-1] + list(u)[::-1],
                 list(u)[:10] + list(c2) + list(u)[10:], list(c2)[:5] + list(u)][hsel]
        kept = order[:hlate]                     # hsel: which table wires are kept (the pool decides whether a
                                                 # special bit costs 2 ANDs or 3)
        if P.zone:
            # Zone format: the odd x of zone 1 in 192..255 (32 inputs: e = 7, address 192 + 2k, digit word 2^13)
            # can have the top lo wire set, which the lo stream does not deliver at their level. The phase
            # module adds the bits of "hi word + top lo bit" to the hi stream as [e < 8] AND xq[13] AND B, with
            # B an XOR of pool wires; so every bit that the carry of that +1 reaches (a "target") must be, on
            # these 32 inputs, an XOR of wires that are alive during the stream. The address, the hi word and
            # as many table wires as fit next to the stream give at most 29 of the 32 dimensions; ALL table
            # wires together give 32 [checked below]. So before the table wires are measured, each target is
            # written by CNOTs onto one of the table wires that would be measured anyway (in place: no new
            # wire, no Toffoli). These wires stay alive through the stream with the kept ones and are measured
            # after it like every table wire (they are functions of the address).
            zin = [192 + 2 * k_ for k_ in range(32)]
            assert all(src_addr[a] == a for a in zin)
            ev_z = M.evaluate(P, tabs, np.array([a - 128 + 129 for a in zin], dtype=np.int64))   # x = 193 + 2k
            assert [int(a) for a in ev_z["a8"]] == zin and all(int(t_) == 1 << 13 for t_ in ev_z["t"])
            ztarget = {}
            for k_, (a, lo21) in enumerate(zip(zin, ev_z["lo"])):
                have = int(tabs["H"][a]) + 0x700             # hi word N and the ones above it (level j = 3)
                diff = ((have + ((int(lo21) >> 1) >> 19)) ^ have) & 0x7FF
                for b_ in range(11):
                    if (diff >> b_) & 1:
                        ztarget[10 - b_] = ztarget.get(10 - b_, 0) | (1 << k_)

            def zvec(w):
                return sum(((pool.val[w] >> a) & 1) << k_ for k_, a in enumerate(zin))

            cand = [w for w in list(u) + list(c2) + list(scr_a)]
            zwires = []
            for q_ in sorted(ztarget, reverse=True):
                # elimination over {1, address, hi word, wires that already hold a target} (no tag) and the
                # other table wires as they are now (tag = which of them)
                zlead = {}

                def zred(v_, tag_):
                    while v_:
                        h_ = v_.bit_length() - 1
                        if h_ not in zlead:
                            return v_, tag_
                        v_, tag_ = v_ ^ zlead[h_][0], tag_ ^ zlead[h_][1]
                    return 0, tag_

                for v_ in [(1 << 32) - 1] + [zvec(w) for w in list(abits) + list(H) + zwires]:
                    v_, _ = zred(v_, 0)
                    if v_:
                        zlead[v_.bit_length() - 1] = (v_, 0)
                rest = [w for w in cand if w not in zwires]
                for j_, w in enumerate(rest):
                    v_, tag_ = zred(zvec(w), 1 << j_)
                    if v_:
                        zlead[v_.bit_length() - 1] = (v_, tag_)
                v_, tag_ = zred(ztarget[q_], 0)
                assert v_ == 0, "zone special bit: not an XOR of the address, the hi word and the table wires"
                if tag_ == 0:
                    continue                                 # already an XOR of what stays alive anyway
                ws = [w for j_, w in enumerate(rest) if (tag_ >> j_) & 1]
                host = ws[-1]                                # the target is written onto this table wire
                pool.xor_into(host, set(ws[:-1]))
                names[host] = f"z{q_}"
                zwires.append(host)
            sp_ = Span()
            sp_.add((1 << 32) - 1, 1)
            for w in list(abits) + list(H) + zwires:
                sp_.add(zvec(w), 1)
            assert all(sp_.reduce(v_)[0] == 0 for v_ in ztarget.values())
            kept = zwires + [w for w in order if w not in zwires][:max(0, hlate - len(zwires))]
            c.meta["zone_positions"] = sorted(ztarget)
        defer([w for w in list(scr_a) + list(c2) + list(u) if w not in kept])
    else:
        defer(list(scr_a) + list(c2) + list(u))
    mono, helpers = [], []
    hform = {}                                   # helper wire -> its two factors, as lists of truth tables
    held = {w: pool.val.pop(w) for w in H} if (late and not hlate) else {}     # helpers: over the address only

    def build_chain():
        # every function of the address must be a quadratic form over the address wires and the helper wires:
        # 1 + n + n (n - 1) / 2 >= 256 needs n >= 23 wires; a helper is a product of two dense forms over the
        # wires below it (so it is erased for free afterwards) and adds at most n - 2 dimensions: 17 helpers.
        hrng = random.Random(f"{useed or seed}:unlut")      # useed: seed of this chain alone
        while True:
            basis = pool.phase_basis(list(pool.val))
            if len(basis[0].lead) == nused:      # 256, or the 198 addresses of the zone format
                break
            base = sorted(pool.val)
            best = None
            for _ in range(8 if unl == 1 else unl):          # keep the product that adds most to the span
                A = {x for x in base if hrng.random() < 0.5} or {base[0]}
                B = {x for x in base if hrng.random() < 0.5} or {base[-1]}
                fv = pool.form_val(A) & pool.form_val(B)
                if A == B or fv in (0, pool.form_val(A), pool.form_val(B)):
                    continue
                pool.val[-1] = fv                # probe entry, never a wire
                gain = len(pool.phase_basis(list(pool.val))[0].lead)
                del pool.val[-1]
                if best is None or gain > best[0]:
                    best = (gain, A, B)
            if best is None:
                continue
            h = c.alloc()
            pool.val[h] = 0
            hform[h] = ([pool.val[x] for x in sorted(best[1])], [pool.val[x] for x in sorted(best[2])])
            pool.and_into(h, best[1], 0, best[2], 0)
            helpers.append(h)
        return basis

    basis = None
    if unl and not hlate:
        basis = build_chain()
    elif not unl:
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
        basis = pool.phase_basis(list(pool.val))
    if not hlate:
        for m, fv in deferred:
            pool.fix_phase(fv, m, basis)
            c.bit_free(m)
    c.meta["n_deferred"] = len(deferred)
    c.meta["n_helpers"] = len(helpers)
    c.meta["lut"] = lutlog
    if late:
        # ---------------- UNNORM with the hi stream inside (phase module). The helper wires stay: they are the
        # pool of the small-x special bits and, in the module's hook, what the hi wires are erased against.
        assert unl and ph >= 2
        pool.val.update(held)
        c.stage("UNNORM")

        def move_high(q):
            """Take a wire off a lent input wire (free: one SWAP)."""
            w = c.alloc_high()
            c.swap(q, w)
            if q in pool.val:
                pool.val[w] = pool.val.pop(q)
            c.free(q)
            return w

        # hdrop: the last helpers are erased here (free) and made again in the hook (one AND each): that many
        # wires fewer during the hi stream, which is the widest part of UNNORM
        dropped = []
        for _ in range(hdrop):
            h = helpers.pop()
            dropped.append(hform[h])
            pool.erase(h)
            c.free(h)
        H = [move_high(w) if w < 22 else w for w in H]
        helpers = [move_high(w) if w < 22 else w for w in helpers]
        kept = [move_high(w) if w < 22 else w for w in kept]
        if zw is not None and zw < 22:
            zw = move_high(zw)
        if pd.F < 22:
            pd.F = move_high(pd.F)
        if p < 22:
            p = move_high(p)
        cur = list(abits[:7])
        tgt = [xq[14 + k] for k in range(7)]
        for k in range(7):
            s_, t_ = cur[k], tgt[k]
            if s_ == t_:
                continue
            if c.is_free(t_):
                c.reclaim(t_)
            else:
                assert t_ in cur, "target wire is busy"
            c.swap(s_, t_)
            vs, vt = pool.val.pop(s_, None), pool.val.pop(t_, None)
            if vs is not None:
                pool.val[t_] = vs
            if vt is not None:
                pool.val[s_] = vt
            for j in range(7):
                if cur[j] == t_ and j != k:
                    cur[j] = s_
            cur[k] = t_
            if s_ not in cur:
                c.free(s_)
        for t_ in tgt:
            c.lendable.discard(t_)
        c.reclaim(xq[21])
        # pool of the special bits: address wires, parity, helpers and the hi wires themselves (the module
        # changes pool wires only inside one special bit and restores them before the next stream step);
        # with fewer than about 30 wires a special bit needs two products instead of one
        small = small_x(sorted(pool.val))

        zt = {}

        def zone_helpers():
            """zh=1 (owner B1, DESIGN3 item I2): the zone flag's own AND tree as the first four helpers of the
            unlookup chain, t1 = a3 a2, t2 = (1 + a5)(1 + a4), t3 = t2 (1 + t1), t4 = t3 (1 + p): the backward zone
            rotation then needs only z = a6 AND NOT t4 (1 AND, in the phase module, ph >= 8). They depend on
            a2..a6 and p only, which the rotation does not move, so they are erased after it (z_done)."""
            a2, a3, a4, a5 = (xq[14 + k] for k in (2, 3, 4, 5))

            def mk(A, ca, B, cb):
                h = c.alloc()
                pool.val[h] = 0
                pool.and_into(h, A, ca, B, cb)
                helpers.append(h)
                return h
            t1 = mk({a3}, 0, {a2}, 0)
            t2 = mk({a5}, 1, {a4}, 1)
            t3 = mk({t2}, 0, {t1}, 1)
            t4 = mk({t3}, 0, {p}, 1)
            zt.update(t1=t1, t2=t2, t3=t3, t4=t4, parts={t4: [t3, p], t3: [t2, t1], t2: [a5, a4], t1: [a3, a2]})
            if zh >= 2:
                # zh=2: z = a6 (1 + t4) itself is the fifth helper (a product of two forms): the backward flag
                # then costs nothing beyond the chain (module ph >= 10, zext)
                a6 = xq[20]
                zz = mk({a6}, 0, {t4}, 1)
                zt["z"] = zz
                zt["parts"][zz] = [a6, t4]

        def z_done():
            for h in ([zt["z"]] if "z" in zt else []) + [zt["t4"], zt["t3"], zt["t2"], zt["t1"]]:
                pool.erase(h, others=zt["parts"][h])
                c.free(h)

        def hi_done():
            if hlate:
                # the kept table wires and the hi wires are measured like every table wire; only now the helper
                # chain is built (over the address wires alone), every kept outcome is repaired, and the
                # helpers are erased
                defer(list(kept) + list(H))
                if zh:
                    zone_helpers()
                basis3 = build_chain()
                for m, fv in deferred:
                    pool.fix_phase(fv, m, basis3)
                    c.bit_free(m)
                c.meta["n_deferred"] = len(deferred)
                c.meta["n_helpers"] = len(helpers)
                for h in reversed(helpers[(5 if zh >= 2 else 4):] if zh else helpers):
                    pool.erase(h)
                    c.free(h)
                return
            for fa_, fb_ in reversed(dropped):   # helpers dropped before the stream, in creation order again
                inv = {v: w for w, v in pool.val.items() if w not in H}
                h = c.alloc()
                pool.val[h] = 0
                pool.and_into(h, {inv[v] for v in fa_}, 0, {inv[v] for v in fb_}, 0)
                helpers.append(h)
            basis2 = pool.phase_basis([w for w in pool.val if w not in H])
            for w in H:
                m = c.hmr(w)
                pool.fix_phase(pool.val[w], m, basis2)
                c.bit_free(m)
                pool.drop(w)
                c.free(w)
            for h in reversed(helpers):
                pool.erase(h)
                c.free(h)

        if zh >= 2:
            assert P.zone and hlate and ph >= 10 and not zcol, "zh=2 needs the zone format, hlate and ph >= 10"
            pd.unnorm(pw=p, H=H + [None] * (P.wh - nh), small=small, hook=hi_done,
                      zext=lambda: zt["z"], hook2=z_done)
        elif zh:
            assert P.zone and hlate and ph >= 8 and not zcol, "zh needs the zone format, hlate and ph >= 8"
            pd.unnorm(pw=p, H=H + [None] * (P.wh - nh), small=small, hook=hi_done,
                      zt4=lambda: zt["t4"], hook2=z_done)
        else:
            pd.unnorm(pw=p, H=H + [None] * (P.wh - nh), small=small, hook=hi_done)
        c.meta["ph_info"] = {k: v for k, v in pd.info.items()}
        c.stage("END")
        return c
    for h in reversed(helpers):
        pool.erase(h)
        c.free(h)
    for mw in mono:
        for mk, wire in mw.w.items():
            if bin(mk).count("1") >= 2:
                pool.drop(wire)
        mw.erase()
    # ---------------- UNNORM
    c.stage("UNNORM")
    cur = list(abits[:7])
    tgt = [xq[14 + k] for k in range(7)]
    if p < 22:
        w = c.alloc_high()
        c.swap(p, w)
        c.free(p)
        p = w
    fobj = pd if ph else code
    if fobj.F < 22:                              # the flag was raised again on a lent I/O wire
        w = c.alloc_high()
        c.swap(fobj.F, w)
        c.free(fobj.F)
        fobj.F = w
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
    if ph:
        pd.unnorm(pw=p)                          # parity on a wire of its own, the wire of bit 21 is |0>
        c.meta["ph_info"] = {k: v for k, v in pd.info.items()}
        c.stage("END")
        return c
    c.x(xq[21])
    code.unpark()
    db = code.b
    S.shifter_refl(c, xq, p, db, inverse=True)
    f = S.thermo(c, xq)
    S.linear_of_e(c, f, lambda e: e & 1, p)
    c.free(p)
    for b in range(4):
        S.linear_of_e(c, f, lambda e, b=b: (dfun(e) >> b) & 1, db[b])
        c.free(db[b])
    S.linear_of_e(c, f, lambda e: int(dfun(e) >= 1), code.F)
    c.free(code.F)
    S.linear_of_e(c, f, lambda e: int(dfun(e) >= 2), code.G)
    c.free(code.G)
    S.thermo_erase(c, xq, f)
    c.meta["code_moves"] = dict(code.moves)
    c.stage("END")
    return c
