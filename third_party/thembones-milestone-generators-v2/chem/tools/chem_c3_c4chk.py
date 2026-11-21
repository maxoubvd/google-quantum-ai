"""Builder C3 copy of R58-CHEM-KEPT r58_chem_kept_c4chk.py (fix: the class parity swap never raises StopIteration;
a partition without a swap is skipped). Original docstring follows.

R58-CHEM-KEPT: THE POINT STAGE WITH CLASS-MASK HELPERS (kept products read on BOTH sides, linear and exact).

  low points = 16 = (16 - d) for stage A + a set D of d points for the point stage;
  stage A    = the plain blocks engine on the table restricted to the points outside D (15 - d coordinates): exact
               there, FREE on D (a stage-A form may take any value at a point of D);
  helper j   = ONE AND: [block in class C_j] AND [low point = p_j]  (a block form times a low form), kept on a spare
               wire while stage A runs; the classes of the helpers of one point are DISJOINT sets of blocks;
  unknowns   = per stage-A form (two per product): SPILL at a point p of D (flip the form at p on all blocks) and
               READ of helper j (flip it at p_j on the blocks of C_j): the other form's value at the point, on those
               blocks, goes to the product's fan-out. Cross terms of the two sides of one product are multiples of
               the helpers themselves (classes disjoint) and are cancelled by the helpers' own free fan-out;
  target     = empty as many (output, point) cells of D as the unknowns allow (76 equations a cell); the systems of
               the d points are independent;
  point stage = the live cells, two per AND, by point products (c1 + delta_p2)(c2 + delta_p1).
The plan is built and re-evaluated on all 77 x 16 addresses: every output must be pure (block-only + low-only).

  CHEM_BROT=7 python -I r58_chem_kept_mask.py TABLE.npz --k 1 [--D 13,14,15] [--K 8] [--seed 0] [--seeds 2]
        [--effort 2] [--json out.json]
"""
import argparse
import json
import os
import random
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402

bits = core.bits
NH, NL = 77, 16
ONES76 = (1 << 76) - 1
FULL = (1 << NH) - 1


def red(v):
    return (v >> 1) ^ (ONES76 if v & 1 else 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, default=1)
    ap.add_argument("--D", default="13,14,15")
    ap.add_argument("--K", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--effort", type=int, default=2)
    ap.add_argument("--json", default=None)
    ap.add_argument("--pseeds", type=int, default=4, help="partitions of the blocks tried")
    ap.add_argument("--etries", type=int, default=6, help="random sets of emptied outputs tried per level")
    ap.add_argument("--plan", default=None, help="write the checked plan (for tools/r58_chem_kept_emit.py)")
    ap.add_argument("--onlyK", action="store_true")
    ap.add_argument("--dist", default=None, help="helpers per point, e.g. 6,1,1 (default: every split of K)")
    ap.add_argument("--V", default=None, help="builder C2 block set V of the loader (hex, bit h = block h)")
    ap.add_argument("--planA", default=None, help="stage A from tools/r58_chem_kept_c2a.py (C2's engine)")
    a = ap.parse_args()
    t0 = time.time()
    rng = random.Random(1000 + a.seed)
    planA = json.load(open(a.planA)) if a.planA else None
    if planA:
        a.D = ",".join(str(x) for x in planA["D"])
        a.k = planA["k"]
    core.BROT, core.J = int(os.environ.get("CHEM_BROT", "7")), 5
    w = np.asarray(np.load(a.table)["angle_words"], np.int64)
    an = core.analyze(w, a.k, 4)
    X0 = an.X0.copy()
    nr, nb = X0.shape
    nl = an.nl
    assert an.H == NH and nl == 15, (an.H, nl)
    D = [int(x) for x in a.D.split(",")]
    nD = len(D)
    assert 0 not in D
    LB = []
    for bs in an.basis:
        tt = 0
        for l in range(NL):
            v = 0
            for m in bits(bs["anf"]):
                if (l & m) == m:
                    v ^= 1
            tt |= v << l
        LB.append(tt)
    sp = core.Span()
    for t in range(nl):
        sp.add(LB[t], 1 << t)
    ech = []                                    # echelon of the point indicators: (pivot coordinate, vector)
    mDraw = []
    for p in D:
        r, c = sp.reduce(1 << p)
        assert r == 0
        mDraw.append(c)
        v = c
        for piv, e in ech:
            if v >> piv & 1:
                v ^= e
        assert v
        piv = v.bit_length() - 1
        ech = [(pv, (e ^ v) if e >> piv & 1 else e) for pv, e in ech] + [(piv, v)]
    drop = sorted(pv for pv, _ in ech)
    XA = np.zeros_like(X0)
    for M in range(nr):
        for o in range(nb):
            v = int(X0[M, o])
            for piv, e in ech:
                if v >> piv & 1:
                    v ^= e
            XA[M, o] = v
    ctl = None
    for sd in range(a.seeds):
        an.X0 = X0.copy()
        try:
            p = core.BlockDecomp(an, sd, a.effort).run()
        except (RuntimeError, AssertionError):
            continue
        ctl = len(p) if ctl is None else min(ctl, len(p))
    bestA = None
    if planA:
        bestA = ([tuple(int(x) for x in q) for q in planA["prods"]], f"C2 engine, {planA['pfile']} rep {planA['rep']}",
                 planA["codim"])
    for sd in range(0 if planA else a.seeds):
        an.X0 = XA.copy()
        bd = core.BlockDecomp(an, sd, a.effort)
        try:
            p = bd.run()
        except (RuntimeError, AssertionError):
            continue
        if bestA is None or len(p) < len(bestA[0]):
            bestA = (p, [list(b["ts"]) for b in bd.binfo], list(bd.evec))
    an.X0 = X0
    if bestA is None:
        print("stage A: the plain engine found no plan")
        return
    prodsA, blocksA, codA = bestA
    nA = len(prodsA)
    print(f"k {a.k}: outputs {nb}, control (plain engine, whole table) {ctl}; D {D} (dropped coordinates {drop}); "
          f"stage A on {nl - nD} coordinates, blocks {blocksA}, codim {codA}: {nA} products", flush=True)

    rep = sum(1 << (h * NL) for h in range(NH))
    lowg = [LB[t] * rep for t in range(nl)]
    highg = [None]
    for Mm in range(1, NH):
        g = 0
        for h in range(NH):
            if (h & Mm) == Mm:
                g |= ((1 << NL) - 1) << (h * NL)
        highg.append(g)

    def Lg(mask):
        g = 0
        for t in bits(mask):
            g ^= lowg[t]
        return g

    def Hg(mask):
        g = 0
        for i in bits(mask):
            g ^= highg[i + 1]
        return g

    def blockg(c):
        g = 0
        for h in bits(c):
            g |= ((1 << NL) - 1) << (h * NL)
        return g

    def col(g, p):
        v = 0
        for h in range(NH):
            if g >> (h * NL + p) & 1:
                v |= 1 << h
        return v

    fb = [0] * nb
    for M in range(nr):
        hg = highg[M + 1]
        for o in range(nb):
            v = int(X0[M, o])
            if v:
                fb[o] ^= hg & Lg(v)
    F, G, S = [], [], []
    for (y1, a1, y2, a2, s, _) in prodsA:
        F.append(Lg(y1) ^ Hg(a2))
        G.append(Lg(y2) ^ Hg(a1))
        S.append(s)

    def residual(Fs, Gs, extra=()):
        gb = list(fb)
        for f, g, s in zip(Fs, Gs, S):
            pr = f & g
            for o in bits(s):
                gb[o] ^= pr
        for pr, s in extra:
            for o in bits(s):
                gb[o] ^= pr
        return gb

    def pure_outside(gb, pts):
        for o in range(nb):
            c0 = col(gb[o], 0)
            for l in range(1, NL):
                if l in pts:
                    continue
                if (col(gb[o], l) ^ c0) not in (0, FULL):
                    return False
        return True

    def rho(gb):
        out = {}
        for o in range(nb):
            c0 = col(gb[o], 0)
            for ip, p in enumerate(D):
                out[(o, ip)] = col(gb[o], p) ^ c0
        return out

    gb = residual(F, G)
    assert pure_outside(gb, D), "stage A is not exact outside D"
    R0 = rho(gb)
    Vset = int(a.V, 16) if a.V else 0
    assert bin(Vset).count("1") % 2 == 0

    def parV(x):
        return bin(x & Vset).count("1") & 1

    # the table's kappa: V-sum at low point 0 of every wire's delta; in the output basis: G kap = kappa
    kq = []
    for q in range(len(an.out_wires)):
        d = an.delta_int[q]
        kq.append(sum((d >> (16 * h)) & 1 for h in range(NH) if Vset >> h & 1) & 1)
    rows = [(an.G[q] | (kq[q] << nb)) for q in range(len(an.out_wires))]
    piv = {}
    for r in rows:
        while r & ((1 << nb) - 1):
            hb = (r & ((1 << nb) - 1)).bit_length() - 1
            if hb in piv:
                r ^= piv[hb]
            else:
                piv[hb] = r
                break
        else:
            assert r == 0, "kappa is not consistent with the output basis"
    kap = [0] * nb
    for hb in sorted(piv):
        r = piv[hb]
        v = (r >> nb) & 1
        for b in range(hb):
            if r >> b & 1:
                v ^= kap[b]
        kap[hb] = v
    assert all((bin(an.G[q] & sum(kap[b] << b for b in range(nb))).count("1") & 1) == kq[q]
               for q in range(len(an.out_wires)))
    pi0 = [parV(col(gb[o], 0)) ^ kap[o] for o in range(nb)]  # V-parity still owed by each output after stage A
    oddA = sum(parV(col(F[q], 0)) + parV(col(G[q], 0)) for q in range(nA))
    print(f"C-4: |V| {bin(Vset).count(chr(49))}; stage-A block forms outside H_V: {oddA} of {2 * nA}; outputs whose "
          f"pure block part is outside H_V after stage A: {sum(pi0)} of {nb}", flush=True)
    Fc = [[col(F[q], p) for p in D] for q in range(nA)]
    Gc = [[col(G[q], p) for p in D] for q in range(nA)]

    def attempt(K, pseed, rot, verbose):
        rng = random.Random(1000 + pseed)
        vprint = print if verbose else (lambda *x, **y: None)
        # helpers: K over the points, classes = a random partition of the blocks (class 0 = no helper)
        Kp = list(rot)                              # helpers per point
        helpers = []                                 # (ip, class mask)
        for ip in range(nD):
            k = Kp[ip]
            wts = [2.0 * nA] + [2.0 * nA + nb] * k
            tot = sum(wts)
            sizes = [int(NH * x / tot) for x in wts]
            i = 0
            while sum(sizes) < NH:
                sizes[i % len(sizes)] += 1
                i += 1
            perm = list(range(1, NH))
            rng.shuffle(perm)
            perm = [0] + perm                       # block 0 stays outside every helper class (no constant term)
            posn = sizes[0]
            spare = [h for h in perm[1:sizes[0]]]       # blocks of the no-helper class (block 0 stays)
            for i in range(1, k + 1):
                m = 0
                for h in perm[posn:posn + sizes[i]]:
                    m |= 1 << h
                posn += sizes[i]
                if parV(m):                             # the class indicator must lie in H_V: swap one block
                    # C3 fix: any (block of the class, spare block) pair of opposite V-membership will do; the
                    # design's version raised StopIteration when the class lay wholly inside V (or outside it)
                    sw = next(((hi, ho) for hi in range(NH) if m >> hi & 1 for ho in spare
                               if (Vset >> hi & 1) != (Vset >> ho & 1)), None)
                    if sw is None:
                        raise ValueError("no parity swap for this partition")
                    hin, hout = sw
                    spare.remove(hout)
                    spare.append(hin)
                    m = (m & ~(1 << hin)) | (1 << hout)
                helpers.append((ip, m))
        hl = [[j for j, (ip, _) in enumerate(helpers) if ip == p] for p in range(nD)]

        def unknowns(ip, with_helpers=True):
            U = []
            for q in range(nA):
                for sd in (0, 1):
                    other = Gc[q][ip] if sd == 0 else Fc[q][ip]
                    U.append((("spill", q, sd, ip), other, S[q]))
                    if with_helpers:
                        for j in hl[ip]:
                            U.append((("read", q, sd, j), helpers[j][1] & other, S[q]))
            if with_helpers:
                for j in hl[ip]:
                    for o in range(nb):
                        U.append((("fan", j, 0, o), helpers[j][1], 1 << o))
            return U

        def solve(U, ip, E):
            pos = {o: i for i, o in enumerate(E)}
            lpos = {o: i for i, o in enumerate([o for o in range(nb) if o not in pos])}
            base = 76 * len(E)
            spn = core.Span()
            for ui, (_, vec, s) in enumerate(U):
                rv = red(vec)
                pv = parV(vec)
                c = 0
                if rv or pv:
                    for o in bits(s):
                        if o in pos:
                            c ^= rv << (pos[o] * 76)
                        elif pv:                        # a live cell block form must stay in H_V
                            c ^= 1 << (base + lpos[o])
                if c:
                    spn.add(c, 1 << ui)
            tgt = 0
            for o in E:
                tgt ^= red(R0[(o, ip)]) << (pos[o] * 76)
            for o, i in lpos.items():
                if parV(R0[(o, ip)]):
                    tgt ^= 1 << (base + i)
            r, c = spn.reduce(tgt)
            return (None if r else [U[i] for i in bits(c)]), len(spn)

        def capacity(U, ip, order):
            """the largest number of cells of this point that can be emptied (several random output sets a level)."""
            e = min(nb, len(U) // 76)
            while e > 0:
                for trial in range(a.etries):
                    oo = list(order)
                    if trial:
                        rng.shuffle(oo)
                    sol, rk = solve(U, ip, oo[:e])
                    if sol is not None:
                        return e, rk
                e -= 1
            return 0, 0

        ev = [o for o in range(nb - 1, -1, -1) if not pi0[o]]
        od = [o for o in range(nb - 1, -1, -1) if pi0[o]]
        order = ev + od                              # outputs emptied at every point need an even pure part
        cap0, cap1, info = [], [], []
        for ip in range(nD):
            U0 = unknowns(ip, False)
            e0, _ = capacity(U0, ip, order)
            U1 = unknowns(ip, True)
            e1, rk = capacity(U1, ip, order)
            cap0.append(e0)
            cap1.append(e1)
            info.append((len(U0), len(U1), rk))
            vprint(f"   point {D[ip]}: spill only {len(U0)} unknowns -> {e0} cells emptied; with {Kp[ip]} helpers "
                  f"{len(U1)} unknowns -> {e1} cells emptied (rank {rk} = {76 * e1})", flush=True)

        def stage_cost(live_sets):
            """ANDs of the point stage for the live point set of each output (pairs; triangles for odd counts)."""
            n, odd = 0, 0
            for P in live_sets:
                n += len(P) // 2
                odd += len(P) % 2
            return n + (odd + 1) // 2 if nD >= 3 else n + odd

        def layout(cap):
            """emptied outputs per point: a common tail + disjoint extras (an output never keeps exactly one cell,
            when it can be avoided)."""
            e_min = min(min(cap), len(ev))
            E = [order[:e_min] for _ in range(nD)]
            nxt = e_min
            for ip in range(nD):
                x = cap[ip] - e_min
                E[ip] = E[ip] + order[nxt:nxt + x]
                nxt += x
            return E

        def run(cap, with_helpers):
            cap = list(cap)
            Us = [unknowns(ip, with_helpers) for ip in range(nD)]
            while True:
                fails = [0] * nD
                for trial in range(a.etries):
                    if trial:
                        rng.shuffle(ev)
                        rng.shuffle(od)
                        order[:] = ev + od
                    E = layout(cap)
                    sols, ok = [], True
                    for ip in sorted(range(nD), key=lambda i: cap[i]):
                        sol, _ = solve(Us[ip], ip, E[ip])
                        if sol is None:
                            fails[ip] += 1
                            ok = False
                            break
                        sols.append((ip, sol))
                    if ok:
                        return E, [sol for _, sol in sorted(sols)]
                cap[max(range(nD), key=lambda i: fails[i])] -= 1

        res = {}
        for label, cap, wh in (("spill only", cap0, False), (f"K {K} helpers", cap1, True)):
            E, sols = run(cap, wh)
            F2, G2 = list(F), list(G)
            dg = [sum(1 << (h * NL + p) for h in range(NH)) for p in D]
            hg = [blockg(m) & dg[ip] for (ip, m) in helpers]
            hfan = [0] * len(helpers)
            sF = [[0] * nD for _ in range(nA)]
            sG = [[0] * nD for _ in range(nA)]
            rF = [[0] * len(helpers) for _ in range(nA)]
            rG = [[0] * len(helpers) for _ in range(nA)]
            nread = 0
            for sol in sols:
                for (kind, q, sd, idx), _, _ in sol:
                    if kind == "spill":
                        (sF if sd == 0 else sG)[q][idx] ^= 1
                    elif kind == "read":
                        (rF if sd == 0 else rG)[q][idx] ^= 1
                        nread += 1
                    else:
                        hfan[q] ^= 1 << idx
            for q in range(nA):
                for ip in range(nD):
                    if sF[q][ip]:
                        F2[q] ^= dg[ip]
                    if sG[q][ip]:
                        G2[q] ^= dg[ip]
                for j, (ip, _) in enumerate(helpers):
                    x, y = rF[q][j], rG[q][j]
                    if x:
                        F2[q] ^= hg[j]
                    if y:
                        G2[q] ^= hg[j]
                    if (sF[q][ip] & y) ^ (sG[q][ip] & x) ^ (x & y):      # cross term = the helper itself
                        hfan[j] ^= S[q]
            used = [j for j in range(len(helpers)) if wh]
            extra = [(hg[j], hfan[j]) for j in used]
            gb2 = residual(F2, G2, extra)
            assert pure_outside(gb2, D), "reads disturbed the points outside D"
            R2 = rho(gb2)
            for ip in range(nD):
                assert all(R2[(o, ip)] in (0, FULL) for o in E[ip]), "emptied cells are not empty"
            # point stage: every wire pure block part must end in H_V (tau == kappa of builder C2)
            live = [[ip for ip in range(nD) if o not in E[ip]] for o in range(nb)]
            pi = [parV(col(gb2[o], 0)) ^ kap[o] for o in range(nb)]
            pp, fails = [], 0
            good, bad = {}, []                       # two-cell outputs by their point pair
            for o in range(nb):
                P = list(live[o])
                if len(P) == 0:
                    fails += pi[o]
                elif len(P) == 1:
                    fails += pi[o]
                    pp.append((P[0], None, R2[(o, P[0])], 0, 1 << o))
                elif len(P) == 2:
                    sg = parV(R2[(o, P[0])] & R2[(o, P[1])])
                    if pi[o] ^ sg:
                        bad.append(o)
                    else:
                        good.setdefault(tuple(P), []).append([o, 1 << o])
                else:
                    done = False
                    for z in range(3):
                        x, y = [P[i] for i in range(3) if i != z]
                        if not pi[o] ^ parV(R2[(o, x)] & R2[(o, y)]):
                            pp.append((x, y, R2[(o, x)], R2[(o, y)], 1 << o))
                            pp.append((P[z], None, R2[(o, P[z])], 0, 1 << o))
                            done = True
                            break
                    if not done:
                        fails += 1
                        pp.append((P[0], P[1], R2[(o, P[0])], R2[(o, P[1])], 1 << o))
                        pp.append((P[2], None, R2[(o, P[2])], 0, 1 << o))
            for o in bad:                            # the partner product also feeds o; o loads the difference
                P = tuple(live[o])
                done = False
                for rec in good.get(P, []):
                    o2 = rec[0]
                    c1 = R2[(o, P[0])] ^ R2[(o2, P[0])]
                    c2 = R2[(o, P[1])] ^ R2[(o2, P[1])]
                    if not pi[o] ^ parV(R2[(o2, P[0])] & R2[(o2, P[1])]) ^ parV(c1 & c2):
                        rec[1] |= 1 << o
                        pp.append((P[0], P[1], c1, c2, 1 << o))
                        done = True
                        break
                if not done:
                    fails += 1
                    pp.append((P[0], P[1], R2[(o, P[0])], R2[(o, P[1])], 1 << o))
            for P, recs in good.items():
                for o2, fan in recs:
                    pp.append((P[0], P[1], R2[(o2, P[0])], R2[(o2, P[1])], fan))
            extra2 = list(extra)
            for (i1, i2, c1, c2, s) in pp:
                g = blockg(c1) & dg[i1]
                if i2 is not None:
                    g ^= blockg(c2) & dg[i2]
                extra2.append((g, s))
            gb3 = residual(F2, G2, extra2)
            ok = pure_outside(gb3, [])
            assert ok, "the plan is not exact"
            pure_pp = [0] * nb
            for (i1, i2, c1, c2, sx) in pp:
                if i2 is not None:
                    for o in bits(sx):
                        pure_pp[o] ^= c1 & c2
            oddw = sum(parV(col(gb3[o], 0) ^ pure_pp[o]) ^ kap[o] for o in range(nb))
            oddp = sum(parV(c1) + (parV(c2) if i2 is not None else 0) for (i1, i2, c1, c2, sx) in pp)
            oddh = sum(parV(m) for (_, m) in helpers) if wh else 0
            vprint(f"   C-4 conditions for the 68-AND cache: stage-A forms outside H_V {oddA}; helper classes outside "
                   f"{oddh}; point-stage forms outside {oddp}; WIRES whose pure block part is outside H_V {oddw} "
                   f"(not fixed by the point stage: {fails})", flush=True)
            nh = len(used)
            tot = nA + nh + len(pp)
            ncell = sum(len(e) for e in E)
            vprint(f"{label}: cells emptied per point {[len(e) for e in E]} = {ncell} of {nD * nb}; point stage "
                  f"{len(pp)} ANDs for {nD * nb - ncell} live cells; helpers {nh}; TOTAL {tot} = {nA} + {nh} + {len(pp)}; "
                  f"kept reads {nread}; PLAN CHECKED on all 77 x 16 addresses (every output pure): {ok}", flush=True)
            res[label] = dict(emptied=[len(e) for e in E], point=len(pp), helpers=nh,
                              total=tot + 1000 * (oddw + oddp + oddh), reads=nread, checked=bool(ok),
                              c4=[oddA, oddh, oddp, oddw])
            if verbose and wh and a.plan and K == a.K:
                stA = []
                for q, pr in enumerate(prodsA):
                    y1, a1, y2, a2, sq, ph = (int(x) for x in pr[:6])
                    for ip in range(nD):
                        if sF[q][ip]:
                            y1 ^= mDraw[ip]
                        if sG[q][ip]:
                            y2 ^= mDraw[ip]
                    stA.append([str(y1), str(a1), str(y2), str(a2), str(sq), ph,
                                [j for j in range(len(helpers)) if rF[q][j]],
                                [j for j in range(len(helpers)) if rG[q][j]]])
                json.dump(dict(k=a.k, D=D, nb=nb, total=tot,
                               helpers=[[D[ip], str(m), str(hfan[j])] for j, (ip, m) in enumerate(helpers)],
                               stageA=stA,
                               point=[[D[i1], str(c1), (None if i2 is None else D[i2]), str(c2), str(s)]
                                      for (i1, i2, c1, c2, s) in pp]), open(a.plan, "w"))
                vprint(f"plan written: {a.plan}")
        return res

    out = {}
    for K in ([a.K] if a.onlyK else sorted({a.K, a.K - 1}, reverse=True)):
        best = None
        def comps(n, parts, cap):
            if parts == 1:
                return [(n,)] if n <= cap else []
            return [(x,) + r for x in range(min(n, cap), -1, -1) for r in comps(n - x, parts - 1, x)]
        dists = [tuple(int(x) for x in a.dist.split(","))] if a.dist else comps(K, nD, K)
        for pseed in range(a.seed, a.seed + a.pseeds):
            for rot in dists:
                try:
                    r = attempt(K, pseed, rot, False)
                except ValueError as ex:                # C3: a partition without a parity swap is skipped
                    print(f"   skipped: partition seed {pseed}, helpers per point {rot}: {ex}", flush=True)
                    continue
                t = r[f"K {K} helpers"]["total"]
                if best is None or t < best[0]:
                    best = (t, pseed, rot, r)
        t, pseed, rot, r = best
        attempt(K, pseed, rot, True)
        out[K] = dict(total=t, pseed=pseed, rot=rot, res=r)
        print(f"SUMMARY k {a.k} D {D} K {K}: control {ctl}; stage A {nA}; spill only "
              f"{r['spill only']['total']}; with helpers {t} (partition seed {pseed}, helpers per point {rot}); "
              f"{time.time() - t0:.0f} s", flush=True)
    if a.json:
        json.dump(dict(k=a.k, D=D, control=ctl, stageA=nA, out=out), open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
