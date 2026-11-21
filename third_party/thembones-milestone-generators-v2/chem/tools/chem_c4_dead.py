#!/usr/bin/env python3
"""chem_c4_dead.py: builder C4 (holder 4). Dead words of the LAST batch as free unknowns in a cells-bound group's column
relations (audit R85-CHEM-LOADERS candidate 2's first test).

A dead word (exact don't-care, tools/r62_chem_cache_dead.dead_mask, read-only) of batch K (read by loader K alone when
K = 14) lets every bit of that word take any value. Toggling bit b of slot j at address a = 16 h + l toggles wire (j, b)'s
delta at one address: a rank-one change u_h (x) v_l of that wire's mixed table (u_h: the monomials M containing h; v_l: the
low coordinates of the point function at l). Usable cells (so that nothing else moves): bit b >= 1 (bit 0 is the wire),
the wire is a basis output of the group and appears in no relation (no dependent wire reads it), l != 0 (kappa, the C-4
V-sum at l = 0) and block h not in V (the V-sums of the other rows), block 76 excluded; a conj carrier's bit 5 toggles
with bit 6 (bit 6 XOR bit 5 kept).

Test: unknowns = the group's 43 shared kernel directions an output (from chem_c4_drop's dump, C4_DROPDUMP) + the usable
dead cells; equations = the killed column kept at zero + m random column relations. Prints the consistency over draws.
  python -I chem_c4_dead.py <table.npz> --k 14 --dump k14_g4kill.npz --draws 20
"""
import argparse
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_c4_wb  # noqa: E402
from r62_chem_cache_dead import dead_mask  # noqa: E402

bits = core.bits


def parity(x):
    return bin(x).count("1") & 1


def usable_cells(words0, edges, k, an, wires, outs, V, batches=None):
    """[(output o, h, l, wire (j, b), inV, batch)] for the usable dead cells on the group's outputs: dead words of
    batch k (default) or of the given batches (holder 4, 19:13Z: batch k - 1 too: a toggle of W_{k-1} moves delta_k the
    same way, and moves delta_{k-1}: the loader k - 1 is then rebuilt; a cell dead in both is taken once, own batch first)."""
    # holder 7 (Nov 19 10:35Z, the 18 x 63 shape): the table's own row count and block count (chem_ld2_hb: CHEM_HB, default
    # 77 = the 1,219 / 1,220-row tables): the last (padded) block is LAST = HB - 1 and its REAL rows are l < A - 16 LAST
    # (3 at 1,219 rows: l 0-2 as before; 10 at 1,210 rows); every number is as before at HB 77
    import chem_ld2_hb as hbm
    A = min(words0.shape[0], hbm.NFULL)
    LAST = hbm.LAST
    LREAL = A - 16 * LAST
    dm = dead_mask(words0[:A], edges, 7, 0)
    dep = 0                                         # basis outputs read by a dependent (relation) wire
    for q, g in enumerate(an.G):
        if g and (g & (g - 1)):
            dep |= g
    conj = {c["wire"][0] for c in an.carriers if c["conj"]}
    # holder 5 (Nov 18 22:36Z): C4_DEAD_XV="bt:0xV,..." = the block set V of the OTHER reader of batch bt (a kept loader,
    # V-padded): a toggle of batch bt at a block of that V and a low point 0-2 moves its kappa / block-76 real rows, which
    # its padding cannot absorb (measured: k11's plain S5T failed the C-4 closure after k12's form toggled such cells)
    xv = {}
    for sp in os.environ.get("C4_DEAD_XV", "").split(","):
        if ":" in sp:
            xb, xvv = sp.split(":")
            xv[int(xb)] = int(xvv, 16)
    # holder 5 (Nov 19 00:50Z): C4_DEAD_LPTS="bt:lo-hi,..." = the only low points at which batch bt's dead words may be
    # toggled. A kept reader's S5T+ reads its batches exactly at low points 0-12 (stage A) and 13 (phase 2), its dead
    # unknowns sit at 14-15 (builder C3's KT-J, 00:46Z): a joint build in its 'order A' lets the form touch the kept
    # reader's batch only at points 13-15, the kept reader's phases 1-2 then re-solved on the new words.
    lp = {}
    for sp in os.environ.get("C4_DEAD_LPTS", "").split(","):
        if ":" in sp:
            xb, rg = sp.split(":")
            lo, hi = (int(x) for x in rg.split("-"))
            lp[int(xb)] = set(range(lo, hi + 1))
    # holder 5 (Nov 19 01:06Z): C4_DEAD_WONLY="bt:j.b+j.b+...,..." = the only wires (slot j, bit b) of batch bt whose dead
    # words may be toggled (a joint build with a WIRED neighbour: its injected pool-free groups 0-1 sit on the other wires)
    wo = {}
    for sp in os.environ.get("C4_DEAD_WONLY", "").split(","):
        if ":" in sp:
            xb, ws = sp.split(":")
            wo[int(xb)] = {tuple(int(x) for x in w.split(".")) for w in ws.split("+") if w}
    out = []
    for o in outs:
        if dep >> o & 1:
            continue
        q = wires[o]
        j, b = an.out_wires[q]
        if b == 0:
            continue
        if b == 6 and j in conj:
            continue
        for a in range(A):
            h, l = a >> 4, a & 15
            if h >= LAST:
                continue
            inV = (V >> h) & 1
            if inV and l < LREAL:                   # kappa (l = 0) and the last block's real rows cannot absorb it
                continue
            for bt in (batches if batches else [k]):
                if bt in xv and l < LREAL and (xv[bt] >> h) & 1:
                    continue
                if bt in lp and l not in lp[bt]:
                    continue
                if bt in wo and (j, b) not in wo[bt]:
                    continue
                if 0 <= bt < 15 and dm[a, 5 * bt + j]:
                    out.append((o, h, l, (j, b), bool(inV), bt))   # inV: block 76's padding row (76, l) toggles with it
                    break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, default=14)
    ap.add_argument("--dump", required=True)
    ap.add_argument("--draws", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    import chem_ld2_e19grp as eg
    z = np.load(a.table)
    w0 = np.asarray(z["angle_words"], np.int64)
    edges = np.asarray(z["tree_edges"], np.int64)
    V, words = eg.block_set(w0, a.k)
    an = core.analyze(words, a.k, 4)
    wires = chem_c4_wb.wire_basis(an)
    d = np.load(a.dump)
    outs = [int(x) for x in d["outs"]]
    live = int(d["live"])
    t = int(d["killed"][0])
    R0 = [d["R0"][i].astype(np.uint32) for i in range(len(outs))]
    N = [d["N"][i].astype(np.uint32) for i in range(d["N"].shape[0])]
    nr, K, b = int(d["nr"]), len(N), len(outs)
    cells = usable_cells(w0, edges, a.k, an, wires, outs, V)
    by_o = {o: [c for c in cells if c[0] == o] for o in outs}
    print(f"k{a.k} group {outs}: wires {[an.out_wires[wires[o]] for o in outs]}, |V| {bin(V).count('1')}, "
          f"usable dead cells by output {[len(by_o[o]) for o in outs]}, blocks {len({c[1] for c in cells})}")
    # effects: E rows (nr ints, coordinate masks) of a toggle at (h, l)
    def v_low(l):
        f = 0
        for m in range(1, 16):
            if (m & l) == l:
                f ^= an.mono_coords[m]
        return f

    def effect(h, l, inV=False):
        vl = v_low(l)
        e = [vl if ((M & h) == h) else 0 for M in range(1, nr + 1)]
        if inV:                                     # the padding row (LAST, l) of the same wire toggles too
            import chem_ld2_hb as hbm
            e = [x ^ (vl if ((M & hbm.LAST) == hbm.LAST) else 0) for x, M in zip(e, range(1, nr + 1))]
        return e

    eff = {o: [effect(c[1], c[2], c[4]) for c in by_o[o]] for o in outs}
    # fix round 2: a toggle's content on the killed column t is compensated by pool fan-outs (the pool's t-columns span
    # all but one direction): each cell's effect becomes E + (pool combination with the same t-column); a cell outside
    # the span is paired with the previous such cell of the same output (their XOR is inside)
    comp = os.environ.get("C4_DEAD_COMP", "1") == "1" and "P" in d.files
    if comp:
        P = d["P"].astype(np.uint32)
        tcols = []
        for pi in range(P.shape[0]):
            v = 0
            for r in range(nr):
                if (int(P[pi, r]) >> t) & 1:
                    v |= 1 << r
            tcols.append(v)
        piv = {}
        for pi, v in enumerate(tcols):
            tag = 1 << pi
            while v:
                hh = v.bit_length() - 1
                if hh not in piv:
                    piv[hh] = (v, tag)
                    break
                v ^= piv[hh][0]
                tag ^= piv[hh][1]
        def express(c):
            tag = 0
            while c:
                hh = c.bit_length() - 1
                if hh not in piv:
                    return None
                c ^= piv[hh][0]
                tag ^= piv[hh][1]
            return tag
        for o in outs:
            new_e, pend = [], None
            for E in eff[o]:
                c = 0
                for r in range(nr):
                    if (E[r] >> t) & 1:
                        c |= 1 << r
                tag = express(c)
                if tag is None:
                    if pend is None:
                        pend = (E, c)
                        continue
                    E = [x ^ y for x, y in zip(E, pend[0])]
                    tag = express(c ^ pend[1])
                    pend = None
                    if tag is None:
                        continue
                Ec = list(E)
                for pi in bits(tag):
                    for r in range(nr):
                        Ec[r] ^= int(P[pi, r])
                assert not any((x >> t) & 1 for x in Ec)
                new_e.append(Ec)
            eff[o] = new_e
        print(f"  compensated through the pool: unknowns by output {[len(eff[o]) for o in outs]}", flush=True)
    lc = bits(live)
    rng = random.Random(a.seed)
    # unknown index: per output o (gi): K kernel bits, then its dead cells
    offs, n = {}, 0
    for gi, o in enumerate(outs):  # (unknown counts after the compensation)
        offs[o] = n
        n += K + len(eff[o])
    for m in [int(x) for x in os.environ.get("C4_DEAD_MS", "2,3,4,5").split(",")]:
        ok = 0
        short = []
        for draw in range(a.draws):
            rows = []                               # equations as (mask over n unknowns, rhs)
            # the killed column t stays zero for the dead toggles (kernel directions are zero there already)
            for o in outs:
                for r in range(nr):
                    mask = 0
                    for i, E in enumerate(eff[o]):
                        if (E[r] >> t) & 1:
                            mask |= 1 << (offs[o] + K + i)
                    if mask:
                        rows.append((mask, 0))   # empty once compensated
            for _ in range(m):
                vs = {o: rng.getrandbits(15) & live for o in outs}
                for r in range(nr):
                    mask, rhs = 0, 0
                    for o in outs:
                        rhs ^= parity(int(R0[outs.index(o)][r]) & vs[o])
                        for j in range(K):
                            if parity(int(N[j][r]) & vs[o]):
                                mask |= 1 << (offs[o] + j)
                        for i, E in enumerate(eff[o]):
                            if parity(E[r] & vs[o]):
                                mask |= 1 << (offs[o] + K + i)
                    rows.append((mask, rhs))
            # Gaussian elimination with rhs
            piv = {}
            bad = 0
            for mask, rhs in rows:
                v = mask | (rhs << n)
                while v & ((1 << n) - 1):
                    hh = (v & ((1 << n) - 1)).bit_length() - 1
                    e = piv.get(hh)
                    if e is None:
                        piv[hh] = v
                        break
                    v ^= e
                else:
                    if v >> n & 1:
                        bad += 1
            if bad == 0:
                ok += 1
            short.append(bad)
        print(f"  m = {m} relations: consistent on {ok} of {a.draws} draws; inconsistent equations per draw "
              f"min {min(short)} median {sorted(short)[len(short) // 2]}; unknowns {n}, equations about "
              f"{len(rows)}", flush=True)


if __name__ == "__main__":
    main()


# ------------------------------------------------------------------------------------- the engine hook (build step)
HOOK = {"cells": None, "m": 6, "conj": set(), "stats": {"calls": 0, "solved": 0, "fail": 0}}


def _effect(an, nr, h, l, inV):
    vl = 0
    for mm in range(1, 16):
        if (mm & l) == l:
            vl ^= an.mono_coords[mm]
    import chem_ld2_hb as hbm
    LAST = hbm.LAST                                 # holder 7: the padded last block (76 at 77 blocks, 75 at 76)
    e = np.zeros(nr, np.uint32)
    for M in range(1, nr + 1):
        bit = ((M & h) == h)
        if inV:
            bit ^= ((M & LAST) == LAST)
        if bit:
            e[M - 1] = vl
    return e


def solve_dead(dec, outs, live, killed, R0, N, K, P, usable, sols, kern, nr, b, rows_of, ncell):
    """Column relations with the dead toggles (compensated through the pool) as extra unknowns. Applies nothing;
    returns (chosen pool combinations per output, list of chosen cells with their raw effect) or None."""
    import chem_c4_drop as dr
    HOOK["stats"]["calls"] += 1
    cells_all = HOOK["cells"]
    an = dec.an
    t = killed[0]
    liveu = np.uint32(live)
    gi_of = {o: gi for gi, o in enumerate(outs)}
    cells = [c for c in cells_all if c[0] in gi_of]
    if not cells:
        return None
    # pool t-columns and their span (tags over usable indices)
    piv = {}
    for pi in range(len(usable)):
        v = 0
        for r in range(nr):
            if (int(P[pi][r]) >> t) & 1:
                v |= 1 << r
        tag = 1 << pi
        while v:
            hh = v.bit_length() - 1
            if hh not in piv:
                piv[hh] = (v, tag)
                break
            v ^= piv[hh][0]
            tag ^= piv[hh][1]

    def express(c):
        tag = 0
        while c:
            hh = c.bit_length() - 1
            if hh not in piv:
                return None
            c ^= piv[hh][0]
            tag ^= piv[hh][1]
        return tag
    unk = {gi: [] for gi in range(b)}               # per output: (raw effect sum, compensation tag, cell list)
    pend = {}
    for c in cells:
        gi = gi_of[c[0]]
        E = _effect(an, nr, c[1], c[2], c[4])
        cvec = 0
        for r in range(nr):
            if (int(E[r]) >> t) & 1:
                cvec |= 1 << r
        tag = express(cvec)
        group = [c]
        if tag is None:
            if gi not in pend:
                pend[gi] = (E, cvec, c)
                continue
            E0, c0, cc0 = pend.pop(gi)
            E = E ^ E0
            tag = express(cvec ^ c0)
            group = [cc0, c]
            if tag is None:
                continue
        Ec = E.copy()
        for pi in bits(tag):
            Ec ^= P[pi]
        assert not ((Ec >> np.uint32(t)) & 1).any()
        unk[gi].append((E, tag, group, Ec & liveu))
    rng = dec.rng
    m = HOOK["m"]
    offs, n = [], 0
    for gi in range(b):
        offs.append(n)
        n += K + len(unk[gi])
    for _trial in range(6):
        cons = []
        for _ in range(m):
            vs = [rng.getrandbits(32) & live for _ in range(b)]
            for r in range(nr):
                mask, rhs = 0, 0
                for gi in range(b):
                    rhs ^= parity(int(R0[gi][r]) & vs[gi])
                    for j in range(K):
                        if parity(int(N[j][r]) & vs[gi]):
                            mask |= 1 << (offs[gi] + j)
                    for i, u in enumerate(unk[gi]):
                        if parity(int(u[3][r]) & vs[gi]):
                            mask |= 1 << (offs[gi] + K + i)
                cons.append((mask, rhs))
        x = _solve(cons, n, rng)
        if x is None:
            continue
        cols = []
        chosen, picked = [], []
        for gi in range(b):
            col = R0[gi].copy()
            comb = sols[gi]
            xs = (x >> offs[gi]) & ((1 << (K + len(unk[gi]))) - 1)
            for j in range(K):
                if xs >> j & 1:
                    col ^= N[j]
                    comb ^= kern[j]
            for i, u in enumerate(unk[gi]):
                if xs >> (K + i) & 1:
                    col ^= u[3]
                    comb ^= u[1]
                    picked.append((outs[gi], u[0], u[2]))
            cols.append(col)
            chosen.append(comb)
        if dr._rank(rows_of(cols)) <= ncell - m:
            HOOK["stats"]["solved"] += 1
            return chosen, picked
    HOOK["stats"]["fail"] += 1
    return None


def _solve(cons, n, rng):
    piv = {}
    for mask, rhs in cons:
        v = mask | (rhs << n)
        while v & ((1 << n) - 1):
            h = (v & ((1 << n) - 1)).bit_length() - 1
            e = piv.get(h)
            if e is None:
                piv[h] = v
                break
            v ^= e
        else:
            if v >> n & 1:
                return None
    x = 0
    for i in range(n):
        if i not in piv and rng.random() < 0.5:
            x |= 1 << i
    for h in sorted(piv):
        v = piv[h]
        rest = v & ((1 << n) - 1) & ~(1 << h)
        if (v >> n & 1) ^ (bin(rest & x).count("1") & 1):
            x |= 1 << h
    return x


def install(cells, m, conj):
    import chem_c4_drop as dr
    HOOK["cells"], HOOK["m"], HOOK["conj"] = cells, m, conj
    dr.DEAD_HOOK = solve_dead


def write_toggles(words0, k, an, wires, toggles, conj):
    """Words with every chosen dead toggle applied (bit b of slot j at address 16 h + l of the cell's batch (batch k for a
    5-tuple); with bit 6 too when b = 5 in a conj-carrier slot). toggles: list of (o, h, l, (j, b), inV[, batch]) cells."""
    w = words0.copy()
    for c in toggles:
        o, h, l, (j, b), inV = c[:5]
        bt = c[5] if len(c) > 5 else k
        a = 16 * h + l
        w[a, 5 * bt + j] ^= 1 << b
        if b == 5 and j in conj:
            w[a, 5 * bt + j] ^= 1 << 6
    return w


# ------------------------------------------------------------------- row-bound groups: dead-paid left-kernel vectors
HOOK_ROW = {"cells": None, "m": 4, "target": None, "stats": {"calls": 0, "solved": 0, "fail": 0}}


def solve_dead_row(dec, outs, live, killed, R0, N, K, P, usable, sols, kern, nr, b, rows_of, base):
    """A row-bound group (rank = rows - 1): the other outputs fixed, their left kernel L holds >= 6 dependencies; the
    output with the most dead cells takes m of them (14 m equations on its kernel bits + pool-compensated toggles).
    Returns (chosen pool combinations per output, picked cells) or None."""
    import chem_c4_drop as dr
    gi_now = getattr(dec, "_gi", -1)
    specs = HOOK_ROW.get("specs")                   # holder 4, 18:42Z: several groups, each (vectors, rank asked)
    if specs:
        if gi_now not in specs:
            return None
        m, rk_spec = specs[gi_now]
    else:
        if HOOK_ROW["target"] is not None and gi_now != HOOK_ROW["target"]:
            return None
        m, rk_spec = HOOK_ROW["m"], HOOK_ROW.get("rk")
    HOOK_ROW["stats"]["calls"] += 1
    HOOK_ROW["stats"][f"calls_g{gi_now}"] = HOOK_ROW["stats"].get(f"calls_g{gi_now}", 0) + 1
    an = dec.an
    t = killed[0] if killed else None
    if t is None:
        return None
    liveu = np.uint32(live)
    cells = [c for c in HOOK_ROW["cells"] if c[0] in outs]
    if not cells:
        return None
    cnt = {o: sum(1 for c in cells if c[0] == o) for o in outs}
    ostar = max(outs, key=lambda o: cnt[o])
    gs = outs.index(ostar)
    piv = {}
    for pi in range(len(usable)):
        v = 0
        for r in range(nr):
            if (int(P[pi][r]) >> t) & 1:
                v |= 1 << r
        tag = 1 << pi
        while v:
            hh = v.bit_length() - 1
            if hh not in piv:
                piv[hh] = (v, tag)
                break
            v ^= piv[hh][0]
            tag ^= piv[hh][1]

    def express(c):
        tag = 0
        while c:
            hh = c.bit_length() - 1
            if hh not in piv:
                return None
            c ^= piv[hh][0]
            tag ^= piv[hh][1]
        return tag
    unk, pend = [], None
    for c in cells:
        if c[0] != ostar:
            continue
        E = _effect(an, nr, c[1], c[2], c[4])
        cv = 0
        for r in range(nr):
            if (int(E[r]) >> t) & 1:
                cv |= 1 << r
        tag = express(cv)
        group = [c]
        if tag is None:
            if pend is None:
                pend = (E, cv, c)
                continue
            E = E ^ pend[0]
            tag = express(cv ^ pend[1])
            group = [pend[2], c]
            pend = None
            if tag is None:
                continue
        Ec = E.copy()
        for pi in bits(tag):
            Ec ^= P[pi]
        unk.append((E, tag, group, Ec & liveu))
    rk_need = rk_spec or (nr - m)                   # holder 4, 18:33Z: the rank asked, apart from the vectors' number
    st = HOOK_ROW["stats"]
    st["cells_ostar"] = max(st.get("cells_ostar", 0), cnt[ostar])
    st["unk"] = max(st.get("unk", 0), K + len(unk))
    st["base"] = base
    rng = dec.rng
    lc = bits(live)
    for _trial in range(int(os.environ.get("C4_DEAD_ROW_TRIALS", "40"))):
        xs = [rng.getrandbits(K) if K else 0 for _ in range(b)]
        cols = []
        for gi in range(b):
            v = R0[gi].copy()
            for j in bits(xs[gi]):
                v ^= N[j]
            cols.append(v)
        others = [cols[gi] for gi in range(b) if gi != gs]
        A = rows_of(others)
        _, lk = dr._elim_with_comb(A)            # combinations of rows that vanish on the other outputs
        if len(lk) < m:
            continue
        lams = []
        ind = {}
        for _ in range(200):
            lam = 0
            for q in lk:
                if rng.random() < 0.5:
                    lam ^= q
            if not lam:
                continue
            v, ok = lam, True
            while v:
                hh = v.bit_length() - 1
                if hh not in ind:
                    ind[hh] = v
                    break
                v ^= ind[hh]
            else:
                ok = False
            if ok:
                lams.append(lam)
            if len(lams) == m:
                break
        if len(lams) < m:
            continue
        nU = K + len(unk)
        cons = []
        for lam in lams:
            rows_l = bits(lam)
            base_v = int(np.bitwise_xor.reduce(R0[gs][rows_l])) if rows_l else 0
            nv = [int(np.bitwise_xor.reduce(N[j][rows_l])) if rows_l else 0 for j in range(K)]
            uv = [int(np.bitwise_xor.reduce(u[3][rows_l])) if rows_l else 0 for u in unk]
            for tt in lc:
                mask = 0
                for j in range(K):
                    if (nv[j] >> tt) & 1:
                        mask |= 1 << j
                for i in range(len(unk)):
                    if (uv[i] >> tt) & 1:
                        mask |= 1 << (K + i)
                cons.append((mask, (base_v >> tt) & 1))
        x = _solve(cons, nU, rng)
        if x is None:
            continue
        newc = list(cols)
        v = R0[gs].copy()
        comb = sols[gs]
        picked = []
        for j in range(K):
            if x >> j & 1:
                v ^= N[j]
                comb ^= kern[j]
        for i, u in enumerate(unk):
            if x >> (K + i) & 1:
                v ^= u[3]
                comb ^= u[1]
                picked.append((ostar, u[0], u[2]))
        newc[gs] = v
        rk = dr._rank(rows_of(newc))
        st["consistent"] = st.get("consistent", 0) + 1
        if rk > rk_need:
            st["rank_short"] = st.get("rank_short", 0) + 1
            continue
        chosen = []
        for gi in range(b):
            if gi == gs:
                chosen.append(comb)
            else:
                cc = sols[gi]
                for j in bits(xs[gi]):
                    cc ^= kern[j]
                chosen.append(cc)
        HOOK_ROW["stats"]["solved"] += 1
        HOOK_ROW["stats"][f"solved_g{gi_now}"] = HOOK_ROW["stats"].get(f"solved_g{gi_now}", 0) + 1
        HOOK_ROW["stats"][f"rank_g{gi_now}"] = (base, rk)
        HOOK_ROW["last_rank"] = (base, rk)
        return chosen, picked
    HOOK_ROW["stats"]["fail"] += 1
    return None


def install_row(cells, m, target_gi, conj, rk=None, specs=None):
    import chem_c4_drop as dr
    HOOK_ROW["cells"], HOOK_ROW["m"], HOOK_ROW["target"], HOOK_ROW["rk"] = cells, m, target_gi, rk
    HOOK_ROW["specs"] = specs
    HOOK["conj"] = conj
    dr.DEAD_ROW_HOOK = solve_dead_row
