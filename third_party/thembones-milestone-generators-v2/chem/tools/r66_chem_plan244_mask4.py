"""R66-CHEM-PLAN244 (second holder): tools/r66_chem_plan244_mask3.py with ONE addition, marked SL:
  SL --slots : SLOT FILLS. A stage-A product with ONE feature (a "single": (low form) x (block form), as the
            staircase stage A of tools/r66_chem_plan244_stairA.py --lucky 1 leaves 12 to 18 of them) has an empty
            side. The tool puts the indicator of ONE point p of phase 1 on its block-form side (fixed before the
            solve) and a FREE BLOCK FORM c on its low-form side: the product then adds c(h) at point p to its own
            fan-out and a pure block part elsewhere: 75 linear unknowns a single in the joint system of phase 1
            (the price of half a point product), so F rises by one for every two singles. The V-parity an emptied
            output owes is one more linear equation a output (the pure part c AND (block form) moves it); c is put
            into H_V with the don't-care block. The block-form side's own spills and reads are switched off in
            phase 1 (they would multiply c) and come back in phase 2, where c is known. Needs --tp.
Below: the header of tools/r66_chem_plan244_mask3.py.
R66-CHEM-PLAN244: a copy of tools/r62_chem_bilin_mask2.py with ONE addition, marked TP:
  TP --tp : A TWO-POINT HELPER ON THE 8TH WIRE. W = [block in C3] AND [low point in {p2, p3}] (one AND: a block form
            times the low form delta_p2 + delta_p3). Phase 1 becomes ONE joint linear system of the points p2, p3:
            classes at p2: rest / C1 (helper of p2) / C3, at p3: rest / C2 (helper of p3) / C3, C3 disjoint from C1
            and C2; unknown vectors: spill at p2, spill at p3, reads of the three helpers (5 x 2 nA) against
            2 x 75 equations an output emptied at both points: --F 10 outputs (S5: 8).
            THE CROSS TERMS at C3 of the two sides of a product differ at the two points by
            (s2_F + s3_F) y_W + (s2_G + s3_G) x_W; summed over the products that feed an emptied output this must be
            0 (W's own fan-out cancels the common part): F quadratic bit conditions, met by sampling the solution
            space. C1 and C2 take --c12 blocks (31: 2 lucky bits each, chosen class by class) so that the coupled
            part (spills + reads of W) keeps about 24 free dimensions.
            Phase 2 (the first point last) is S5's, with 2 (30 - F) point-product forms: a 13-block class costs
            22 lucky bits at stage A 149. --lshare DIR lets several processes search the same lucky class and
            share the first find (same --seed, different --lseed). --dist must be 5,1,1 and --K 8.
            Plan: 8 helpers, 30 - F point products: 149 + 8 + 20 = 177.
Below: the header of tools/r62_chem_bilin_mask2.py.
R62-CHEM-BILIN, SECOND TOOL (tools/r62_chem_bilin_mask.py + --s5, marked S5):
  S5 --s5 : THE POINT PRODUCTS READ THE HELPERS OF THE FIRST POINT. Order: (1) points 2 and 3 of D with one helper
            each (8 cells, the same 8 outputs); (2) the 22 point products of the live outputs, as before;
            (3) the first point LAST: its 30 cells are emptied by the stage-A forms AND by the 44 forms of the
            point products (a point product is pure block at the first point: flipping one of its forms on a class
            sends the other form's column there to its own output). A class then has 2 nA + 44 + 30 unknowns:
            12 blocks free, 13 blocks at 18 lucky bits (stage A 149), so FIVE helpers hold the 75 blocks
            (class 0: 11, helpers 13, 13, 13, 13, 12). The lucky classes are found one at a time (an 11-block
            base in reduced echelon form, then every pair of further blocks). Plan: 7 helpers, 22 point products.
            --dist must be 5,1,1. The plan file's point entries get four more fields (reads and spill of each form).
Below: the header of tools/r62_chem_bilin_mask.py.
R62-CHEM-BILIN: a copy of tools/r58_chem_kept_rp_mask.py with ONE addition, marked R62:
  R62 --cls N : CLASS-BY-CLASS CHOICE OF THE HELPER CLASSES. The linear system of a point decouples by class
              (class j: z = spill + read_j and the helper's fan-out are its own unknowns), so a helper class C is
              USABLE for a set of emptied outputs iff its own small system is consistent:
                  sum_f z_f g_f(h) [f feeds o] + phi_o = R0_o(h)   for h in C, o emptied,
              (|C| x cells equations, 2 nA + cells unknowns). When 2 nA + cells < |C| x cells (stage A 149: 328
              against 330) a random partition is consistent with probability 2^-2 a class, 2^-12 for six classes;
              chosen class by class it costs about 4 tries a class. Up to N tries a class; then the old random
              class. Everything else (solve, layout, point stage, C-4 conditions, plan file) is unchanged.
R58-CHEM-KEPT re-pricer's copy (tools/r58_chem_kept_c4chk.py, unchanged except where marked RP):
  RP --dc   : rows >= --nreal (1,219) are never run by the verifier; at a point p of D every block h with 16 h + p
              >= nreal is a DON'T-CARE: its equation is dropped from every emptied cell, and when such a block lies
              in V the live-cell parity row is dropped too (the free value fixes the parity); every check runs on
              the real addresses only.
  RP --eqw  : equal class sizes (each class: (size - 1) equations a cell against 2 nA unknowns); the original
              weights 2 nA + nb for helper classes make them larger than the no-helper class.
Original docstring:
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
    ap.add_argument("--dc", action="store_true", help="RP: blocks of rows >= nreal are don't-cares at D")
    ap.add_argument("--eqw", action="store_true", help="RP: equal class weights")
    ap.add_argument("--nreal", type=int, default=1219, help="RP: rows the verifier runs")
    ap.add_argument("--cls", type=int, default=0, help="R62: tries a class for the class-by-class choice")
    ap.add_argument("--s5", action="store_true", help="S5: the first point last, read by the point products")
    ap.add_argument("--s5bases", type=int, default=20000, help="S5: bases tried for one lucky class")
    ap.add_argument("--tp", action="store_true", help="TP: a two-point helper on the 8th wire (needs --s5)")
    ap.add_argument("--F", type=int, default=10, help="TP: outputs emptied at both points of phase 1")
    ap.add_argument("--c12", type=int, default=31, help="TP: blocks of the one-point helper classes C1, C2")
    ap.add_argument("--c3", type=int, default=20, help="TP: blocks of the two-point class C3")
    ap.add_argument("--tptries", type=int, default=300, help="TP: draws of (output set, classes)")
    ap.add_argument("--qsamples", type=int, default=200000, help="TP: samples of the solution space a draw")
    ap.add_argument("--phase1only", action="store_true", help="TP: stop after the exact check of phase 1")
    ap.add_argument("--samecls", action="store_true", help="TP: one common class at both points (dead: rank 806)")
    ap.add_argument("--lshare", default=None, help="TP: folder shared by the processes that search a lucky class")
    ap.add_argument("--lseed", type=int, default=0, help="TP: this process's seed for the lucky-class search")
    ap.add_argument("--slots", action="store_true", help="SL: slot fills on the singles of stage A (needs --tp)")
    a = ap.parse_args()
    if a.tp:
        assert a.s5 and a.dc, "--tp needs --s5 and --dc"
    lrng = random.Random(77000 + a.lseed)           # TP: the lucky-class search's own generator
    t0 = time.time()
    rng = random.Random(1000 + a.seed)
    planA = json.load(open(a.planA)) if a.planA else None
    if planA:
        a.D = ",".join(str(x) for x in planA["D"])
        a.k = planA["k"]
    core.BROT, core.J = int(os.environ.get("CHEM_BROT", "7")), 5
    w = np.asarray(np.load(a.table)["angle_words"], np.int64)
    zt = np.load(a.table)                       # RP guard: the run rows must match the table's N + R B
    if "W" in zt.files and "U" in zt.files:
        nt = int(zt["U"].shape[-1] + zt["W"].shape[0] * zt["W"].shape[1])
        assert a.nreal == nt, f"--nreal {a.nreal} but the table has {nt} rows (N + R B)"
    assert a.nreal <= w.shape[0]
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
    slot_of = {}                                # SL: product -> index of its point in D (1 or 2)
    if a.slots:
        assert a.tp and nD == 3
        ns_ = 0
        for q, pr in enumerate(prodsA):
            y1, a1, y2, a2 = (int(x) for x in pr[:4])
            if a2 == 0 and a1 != 0 and y1 != 0:    # (low form y1) x (block form a1): the low side is empty
                ip = 1 + (ns_ & 1)
                slot_of[q] = ip
                ns_ += 1
                prodsA[q] = (y1, a1, mDraw[ip], 0) + tuple(pr[4:])   # the block-form side takes the indicator
        print(f"SL: {len(slot_of)} singles of stage A take a slot (points "
              f"{[D[ip] for ip in sorted(set(slot_of.values()))]}): {75 * len(slot_of)} more unknowns in phase 1",
              flush=True)
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

    # RP: real blocks at each low point (all blocks unless --dc), don't-care blocks at each point of D
    realL = [FULL] * NL
    if a.dc:
        for l in range(NL):
            realL[l] = sum(1 << h for h in range(NH) if 16 * h + l < a.nreal)
    assert realL[0] == FULL
    dcb = [[h for h in range(NH) if not realL[p] >> h & 1] for p in D]
    m76 = []
    for ip in range(nD):
        m = ONES76
        for h in dcb[ip]:
            m &= ~(1 << (h - 1))
        m76.append(m)
    print(f"RP: dc {a.dc} (nreal {a.nreal}), equal weights {a.eqw}; don't-care blocks at D {D}: {dcb}", flush=True)

    def pure_outside(gb, pts):
        for o in range(nb):
            c0 = col(gb[o], 0)
            for l in range(1, NL):
                if l in pts:
                    continue
                if ((col(gb[o], l) ^ c0) & realL[l]) not in (0, realL[l]):
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

    dcV = [next((h for h in dcb[ip] if Vset >> h & 1), None) for ip in range(nD)]   # RP: a free block in V

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

    def class_ok(ip, m, E):
        """R62: is the helper class m (a block mask without block 0) consistent by itself for the emptied set E?"""
        spn = core.Span()
        pos = {o: i for i, o in enumerate(E)}
        for q in range(nA):
            so = [pos[o] for o in bits(S[q]) if o in pos]
            if not so:
                continue
            for other in (Gc[q][ip], Fc[q][ip]):
                v = other & m
                if v:
                    c = 0
                    for i in so:
                        c ^= v << (i * NH)
                    spn.add(c, 0)
        for i in range(len(E)):
            spn.add(m << (i * NH), 0)
        tgt = 0
        for o, i in pos.items():
            tgt ^= (R0[(o, ip)] & m) << (i * NH)
        r, _ = spn.reduce(tgt)
        return r == 0

    cls_stat = {"tries": 0, "found": 0, "fallback": 0}

    def attempt(K, pseed, rot, verbose):
        rng = random.Random(1000 + pseed)
        vprint = print if verbose else (lambda *x, **y: None)
        # helpers: K over the points, classes = a random partition of the blocks (class 0 = no helper)
        Kp = list(rot)                              # helpers per point
        K0 = Kp[0]
        if a.s5:
            assert nD == 3 and Kp[1] == 1 and Kp[2] == 1 and K0 >= 2, "--s5 needs --dist k,1,1"
            Kp[0] = 0                               # S5: the first point gets its helpers after the point products
            if a.tp:                                # TP: the three helpers of points 2, 3 are made in phase 1
                Kp[1] = Kp[2] = 0
        helpers = []                                 # (ip, class mask)
        for ip in range(nD):
            k = Kp[ip]
            wts = [2.0 * nA] + [2.0 * nA + (0 if a.eqw else nb)] * k        # RP: --eqw
            tot = sum(wts)
            NHr = NH - len(dcb[ip])                                          # RP: real blocks at this point
            sizes = [int(NHr * x / tot) for x in wts]
            i = 0
            while sum(sizes) < NHr:
                sizes[i % len(sizes)] += 1
                i += 1
            perm = [h for h in range(1, NH) if h not in dcb[ip]]
            rng.shuffle(perm)
            perm = [0] + perm                       # block 0 stays outside every helper class (no constant term)
            # RP: the don't-care blocks stay in the no-helper class (they are outside perm, never swapped)
            posn = sizes[0]
            spare = [h for h in perm[1:sizes[0]]]       # blocks of the no-helper class (block 0 stays)
            chosen = None                                # R62: class-by-class choice of the helper classes
            if a.cls and k >= 2:
                Eall = list(range(nb))
                while (sizes[0] - 1) * nb > 2 * nA and sizes[0] > 1:   # class 0 has no fan-out bits: keep it free
                    sizes[0] -= 1
                    sizes[min(range(1, k + 1), key=lambda j: sizes[j])] += 1
                pool = list(perm[1:])
                rng.shuffle(pool)
                pool = pool[sizes[0] - 1:]               # class 0 = block 0 + the first sizes[0] - 1 blocks

                def mask_of(cand):
                    mm = 0
                    for h in cand:
                        mm |= 1 << h
                    return mm

                def pick(pool, n, also_rest):
                    for _ in range(a.cls):
                        cand = rng.sample(pool, n)
                        mm = mask_of(cand)
                        if parV(mm):
                            continue
                        cls_stat["tries"] += 1
                        if not class_ok(ip, mm, Eall):
                            continue
                        if also_rest:
                            rest = mask_of([h for h in pool if h not in cand])
                            if parV(rest):
                                continue
                            cls_stat["tries"] += 1
                            if not class_ok(ip, rest, Eall):
                                continue
                        return cand, mm
                    return None
                chosen = []
                for i in range(1, k):
                    got = pick(pool, sizes[i], i == k - 1)
                    if got is None:
                        cls_stat["fallback"] += 1
                        raise ValueError("R62: no consistent class in --cls tries")
                    cls_stat["found"] += 1
                    pool = [h for h in pool if h not in got[0]]
                    chosen.append(got[1])
                chosen.append(mask_of(pool))
                cls_stat["found"] += 1
            for i in range(1, k + 1):
                m = 0
                for h in perm[posn:posn + sizes[i]]:
                    m |= 1 << h
                posn += sizes[i]
                if chosen is not None:
                    helpers.append((ip, chosen[i - 1]))
                    continue
                if parV(m):                             # the class indicator must lie in H_V: swap one block
                    # builder C3's fix (tools/chem_c3_c4chk.py), ported: any pair of opposite V-membership
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
                rv = red(vec) & m76[ip]                 # RP: no equation at a don't-care block
                pv = parV(vec) if dcV[ip] is None else 0  # RP: a free block in V fixes a live cell's parity
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
                tgt ^= (red(R0[(o, ip)]) & m76[ip]) << (pos[o] * 76)
            for o, i in lpos.items():
                if dcV[ip] is None and parV(R0[(o, ip)]):
                    tgt ^= 1 << (base + i)
            r, c = spn.reduce(tgt)
            return (None if r else [U[i] for i in bits(c)]), len(spn)

        def capacity(U, ip, order):
            """the largest number of cells of this point that can be emptied (several random output sets a level)."""
            e = min(nb, len(U) // (NH - 1 - len(dcb[ip])))   # R62: 75 equations a cell under --dc (was 76)
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

        def tp_phase1():
            """TP: points 2, 3 of D as ONE linear system: helpers C1 (point 2), C2 (point 3), C3 (both points);
            then the F quadratic bit conditions of the cross terms at C3 by sampling the solution space."""
            assert nD == 3 and dcb[1] == dcb[2] and dcV[1] is not None and dcV[2] is not None and not helpers
            Fn = a.F
            if len(ev) < Fn and not slot_of:
                raise ValueError(f"TP: only {len(ev)} outputs with an even pure part, {Fn} needed")
            realb12 = [h for h in range(1, NH) if h not in dcb[1]]
            j1, j2, jW = 0, 1, 2

            def mask_of(c):
                mm = 0
                for h in c:
                    mm |= 1 << h
                return mm

            def even_class(pool, n, ipk, E10, tries):
                for _ in range(tries):
                    cand = rng.sample(pool, n)
                    mm = mask_of(cand)
                    if parV(mm):
                        continue
                    if ipk is None or class_ok(ipk, mm, E10):
                        return mm
                return None

            stat = dict(draws=0, noclass=0, incons=0, noquad=0)
            for draw in range(a.tptries):
                stat["draws"] += 1
                rng.shuffle(ev)
                order[:] = ev + od
                E10 = list(ev[:Fn])
                if slot_of:                          # SL: the parity rows let any output be emptied
                    allo = list(range(nb))
                    rng.shuffle(allo)
                    E10 = allo[:Fn]
                # TP: the PAIR helper: class mA at point 2, class mB at point 3, DISJOINT (with one common class the
                # two points see the same partner columns up to a constant: the reads then act alike at both points:
                # measured rank 806 for 880 equations, 16:17Z)
                mA = even_class(realb12, a.c3, None, E10, 50)
                mB = even_class([h for h in realb12 if not mA >> h & 1], a.c3, None, E10, 50) if mA is not None else None
                if mA is None or mB is None:
                    stat["noclass"] += 1
                    continue
                if a.samecls:                        # the measured dead variant, kept for the record
                    mB = mA
                rest1 = [h for h in realb12 if not mA >> h & 1]
                rest2 = [h for h in realb12 if not mB >> h & 1]
                m1 = even_class(rest1, a.c12, None if slot_of else 1, E10, 400)
                m2 = even_class(rest2, a.c12, None if slot_of else 2, E10, 400) if m1 is not None else None
                if m1 is None or m2 is None:
                    stat["noclass"] += 1
                    continue
                pos = {o: i for i, o in enumerate(E10)}

                def lift(v1, v2, s):
                    r1 = red(v1) & m76[1]
                    r2 = red(v2) & m76[2]
                    c = 0
                    if r1 or r2:
                        for o in bits(s):
                            i = pos.get(o)
                            if i is not None:
                                c ^= (r1 << (i * 76)) ^ (r2 << ((Fn + i) * 76))
                    return c

                U = []                               # five unknowns a form: spill 2, spill 3, read C1, C2, C3
                for q in range(nA):
                    for sd in (0, 1):
                        o1 = Gc[q][1] if sd == 0 else Fc[q][1]
                        o2 = Gc[q][2] if sd == 0 else Fc[q][2]
                        U.append((("spill", q, sd, 1), lift(o1, 0, S[q])))
                        U.append((("spill", q, sd, 2), lift(0, o2, S[q])))
                        U.append((("read", q, sd, j1), lift(m1 & o1, 0, S[q])))
                        U.append((("read", q, sd, j2), lift(0, m2 & o2, S[q])))
                        U.append((("read", q, sd, jW), lift(mA & o1, mB & o2, S[q])))
                for o in E10:
                    U.append((("fan", j1, 0, o), lift(m1, 0, 1 << o)))
                    U.append((("fan", j2, 0, o), lift(0, m2, 1 << o)))
                    U.append((("fan", jW, 0, o), lift(mA, mB, 1 << o)))
                tgt = 0
                for o in E10:
                    i = pos[o]
                    tgt ^= ((red(R0[(o, 1)]) & m76[1]) << (i * 76)) ^ ((red(R0[(o, 2)]) & m76[2]) << ((Fn + i) * 76))
                dead = 0
                if slot_of:                          # SL
                    pbase = 2 * Fn * 76              # one parity row an emptied output
                    hfree = dcV[1]
                    for q, ip in slot_of.items():
                        for u in range(10 * q + 5, 10 * q + 10):     # the block-form side: off in phase 1
                            U[u] = (U[u][0], 0)
                            dead |= 1 << u
                        g0 = col(G[q], 0)
                        gfree = g0 >> hfree & 1
                        for h in realb12:
                            v = lift(1 << h, 0, S[q]) if ip == 1 else lift(0, 1 << h, S[q])
                            if Vset >> h & 1 and ((g0 >> h & 1) ^ gfree):
                                for o in bits(S[q]):
                                    i = pos.get(o)
                                    if i is not None:
                                        v ^= 1 << (pbase + i)
                            U.append((("slot", q, h, ip), v))
                    for o in E10:
                        if pi0[o]:
                            tgt ^= 1 << (pbase + pos[o])
                if os.environ.get("R66_DEBUG") and draw == 0:   # TP diagnostic: ranks class by class
                    def rank_of(vs):
                        pv = {}
                        for v in vs:
                            while v:
                                b = v.bit_length() - 1
                                e = pv.get(b)
                                if e is None:
                                    pv[b] = v
                                    break
                                v ^= e
                        return len(pv)

                    def rowmask(ipi, m):
                        mm = 0
                        for i in range(Fn):
                            mm |= ((m >> 1) & m76[ipi]) << (((ipi - 1) * Fn + i) * 76)
                        return mm
                    mr1, mr2, mall = mask_of(rest1), mask_of(rest2), mask_of(realb12)
                    allv = [v for _, v in U]
                    for name, ipi, m in (("R@2", 1, mr1 & ~m1), ("C1@2", 1, m1), ("A@2", 1, mA),
                                         ("R@3", 2, mr2 & ~m2), ("C2@3", 2, m2), ("B@3", 2, mB),
                                         ("all@2", 1, mall), ("all@3", 2, mall)):
                        rm = rowmask(ipi, m)
                        print(f"   DEBUG rows {name}: {bin(rm).count('1')} equations, rank {rank_of([v & rm for v in allv])}",
                              flush=True)
                    for name, sel in (("spills", lambda t: t[0] == "spill"), ("reads", lambda t: t[0] == "read"),
                                      ("spill+readW", lambda t: t[0] == "spill" or (t[0] == "read" and t[3] == jW)),
                                      ("fans", lambda t: t[0] == "fan")):
                        vs = [v for t, v in U if sel(t)]
                        print(f"   DEBUG unknowns {name}: {len(vs)} vectors, rank {rank_of(vs)}", flush=True)
                    print(f"   DEBUG |ev| {len(ev)}, E {E10}, fan-out weights on E: "
                          f"{sorted(bin(S[q] & mask_of(E10)).count('1') for q in range(nA))[:20]} ...", flush=True)
                piv, null = {}, []
                for ui, (_, v) in enumerate(U):
                    c = 1 << ui
                    while v:
                        e = piv.get(v.bit_length() - 1)
                        if e is None:
                            break
                        v ^= e[0]
                        c ^= e[1]
                    if v:
                        piv[v.bit_length() - 1] = (v, c)
                    else:
                        null.append(c)
                v, c0 = tgt, 0
                while v:
                    e = piv.get(v.bit_length() - 1)
                    if e is None:
                        break
                    v ^= e[0]
                    c0 ^= e[1]
                if v:
                    stat["incons"] += 1
                    vprint(f"   TP draw {draw}: inconsistent (rank {len(piv)} of {2 * Fn * 75} equations, "
                           f"{len(U)} unknowns)", flush=True)
                    continue
                Em = mask_of(E10)

                def kappa(c):
                    kk = 0
                    for q in range(nA):
                        b = (c >> (10 * q)) & 1023       # F: bits 0-4, G: bits 5-9 (spill 2, spill 3, C1, C2, C3)
                        if ((b ^ (b >> 1)) & (b >> 9) & 1) ^ (((b >> 5) ^ (b >> 6)) & (b >> 4) & 1):
                            kk ^= S[q]
                    return kk & Em

                # the null vectors that move the coupled part (spills, reads of C3) come first
                cm = 0
                for q in range(nA):
                    cm |= 0b1001110011 << (10 * q)
                null = [n for n in null if not n & dead]     # SL: the switched-off unknowns stay at 0
                null.sort(key=lambda n: 0 if n & cm else 1)
                nc = sum(1 for n in null if n & cm)
                found, ns = None, 0
                lim = min(a.qsamples, 1 << min(nc, 40))
                for ns in range(lim):
                    c = c0
                    r = rng.getrandbits(nc) if ns else 0
                    for i in bits(r):
                        c ^= null[i]
                    c &= ~dead                       # SL: zero vectors, kept at 0
                    if not kappa(c):
                        found = c
                        break
                vprint(f"   TP draw {draw}: rank {len(piv)}, {len(U)} unknowns, null space {len(null)} "
                       f"({nc} move the coupled part); cross-term bits met after {ns + 1} samples: "
                       f"{found is not None}", flush=True)
                if found is None:
                    stat["noquad"] += 1
                    continue
                helpers.append((1, m1))
                helpers.append((2, m2))
                helpers.append(((1, 2), (mA, mB)))
                vprint(f"   TP PHASE 1: {Fn} outputs emptied at points {D[1]}, {D[2]} with helpers C1 {a.c12} blocks, "
                       f"C2 {a.c12}, pair helper A {a.c3} at {D[1]} / B {a.c3} at {D[2]}; draws {stat}; {time.time() - t0:.0f} s", flush=True)
                return [[], list(E10), list(E10)], [[(U[i][0], None, None) for i in bits(found)]]
            raise ValueError(f"TP: phase 1 found no plan at F {Fn}: {stat}")

        res = {}
        for label, cap, wh in (("spill only", cap0, False), (f"K {K} helpers", cap1, True)):
            if a.s5 and wh and a.tp:                 # TP phase 1: one joint system of points 2, 3
                E, sols = tp_phase1()
            elif a.s5 and wh:                        # S5 phase 1: points 2, 3 with the same emptied outputs
                Us = [unknowns(ip, True) for ip in range(nD)]
                E = None
                e8 = min(cap1[1], cap1[2], len(ev))
                for trial in range(max(a.etries, 60)):   # cheap: two small solves a trial
                    if trial:
                        rng.shuffle(ev)
                        rng.shuffle(od)
                        order[:] = ev + od
                    E8 = order[:e8]
                    s1, _ = solve(Us[1], 1, E8)
                    s2, _ = solve(Us[2], 2, E8) if s1 is not None else (None, 0)
                    if s1 is not None and s2 is not None:
                        E, sols = [[], list(E8), list(E8)], [s1, s2]
                        break
                if E is None:
                    raise ValueError("S5: phase 1 found no common output set")
            else:
                E, sols = run(cap, wh)
            F2, G2 = list(F), list(G)
            dg = [sum(1 << (h * NL + p) for h in range(NH)) for p in D]
            hg = [(blockg(m) & dg[ip]) if isinstance(ip, int) else       # TP: the pair helper
                  ((blockg(m[0]) & dg[ip[0]]) | (blockg(m[1]) & dg[ip[1]])) for (ip, m) in helpers]
            hfan = [0] * len(helpers)
            sF = [[0] * nD for _ in range(nA)]
            sG = [[0] * nD for _ in range(nA)]
            rF = [[0] * len(helpers) for _ in range(nA)]
            rG = [[0] * len(helpers) for _ in range(nA)]
            nread = 0
            slotc = {}                               # SL: product -> block set of its slot fill
            for sol in sols:
                for (kind, q, sd, idx), _, _ in sol:
                    if kind == "spill":
                        (sF if sd == 0 else sG)[q][idx] ^= 1
                    elif kind == "read":
                        (rF if sd == 0 else rG)[q][idx] ^= 1
                        nread += 1
                    elif kind == "slot":
                        slotc[q] = slotc.get(q, 0) ^ (1 << sd)
                    else:
                        hfan[q] ^= 1 << idx
            for q, c in slotc.items():               # SL: c into H_V by the don't-care block, then into the form
                if parV(c):
                    c ^= 1 << dcV[1]
                F2[q] ^= blockg(c)
                F[q] ^= blockg(c)
                for ip in range(nD):
                    Fc[q][ip] = col(F[q], D[ip])
                amask = 0                            # c in the basis of the high monomials (triangular solve)
                coef = {}
                for Mm in range(1, NH):
                    v = c >> Mm & 1
                    for M2, cf in coef.items():
                        if cf and (Mm & M2) == M2:
                            v ^= 1
                    coef[Mm] = v
                    if v:
                        amask |= 1 << (Mm - 1)
                assert Hg(amask) == blockg(c), "SL: the monomial form of a slot fill"
                pr = prodsA[q]
                prodsA[q] = (pr[0], pr[1], pr[2], amask) + tuple(pr[4:])
            if slotc:
                vprint(f"   SL: {len(slotc)} slot fills placed (block sets of "
                       f"{sorted(bin(c).count('1') for c in slotc.values())} blocks)", flush=True)
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
                    if not isinstance(ip, int):          # TP: the cross term at the helper's first point; the
                        ip = ip[0]                       # difference at its second point sums to 0 on emptied outputs
                    if (sF[q][ip] & y) ^ (sG[q][ip] & x) ^ (x & y):      # cross term = the helper itself
                        hfan[j] ^= S[q]
            used = [j for j in range(len(helpers)) if wh]
            extra = [(hg[j], hfan[j]) for j in used]
            gb2 = residual(F2, G2, extra)
            assert pure_outside(gb2, D), "reads disturbed the points outside D"
            R2 = rho(gb2)
            for ip in range(nD):
                rm = realL[D[ip]]
                assert all((R2[(o, ip)] & rm) in (0, rm) for o in E[ip]), "emptied cells are not empty"
                if dcV[ip] is not None:                  # RP: the free block puts every live column into H_V
                    for o in range(nb):
                        if o not in E[ip] and parV(R2[(o, ip)]):
                            R2[(o, ip)] ^= 1 << dcV[ip]
            if a.tp and wh and a.phase1only:         # TP: the exact check of phase 1 has passed (asserts above)
                print(f"TP PHASE 1 CHECKED EXACT on every real address: {len(E[1])} outputs emptied at points "
                      f"{D[1]} and {D[2]} with 3 helpers (one of them two-point); kept reads {nread}; "
                      f"{time.time() - t0:.0f} s", flush=True)
                sys.exit(0)
            # point stage: every wire pure block part must end in H_V (tau == kappa of builder C2)
            live = [[ip for ip in range(nD) if o not in E[ip] and not (a.s5 and wh and ip == 0)]
                    for o in range(nb)]            # S5: the first point is emptied in phase 2
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
            pp13 = [0] * len(pp)                     # S5: what a point product adds at the first point
            ppr = [([], [], 0, 0) for _ in pp]       # S5: (reads of form 1, reads of form 2, spill 1, spill 2)
            if a.s5 and wh:
                if any(i2 is None for (_, i2, _, _, _) in pp):
                    raise ValueError("S5: a single-cell point product")
                npp = len(pp)
                nU = 2 * nA + 2 * npp + nb
                MU = (1 << nU) - 1
                oth, fanm = [], []
                for q in range(nA):
                    oth += [Gc[q][0], Fc[q][0]]
                    fanm += [S[q], S[q]]
                for (_, _, c1, c2, sx) in pp:
                    oth += [c1, c2]                  # flipping form 1 (block c2) sends c1; form 2 sends c2
                    fanm += [sx, sx]
                fo = [0] * nb
                for u, sx in enumerate(fanm):
                    for o in bits(sx):
                        fo[o] |= 1 << u
                realb = [h for h in range(1, NH) if h not in dcb[0]]
                gmh = {h: sum(1 << u for u, v in enumerate(oth) if v >> h & 1) for h in realb}
                rowc = {}

                def rows_of(h):
                    if h not in rowc:
                        rowc[h] = [(gmh[h] & fo[o]) | (1 << (2 * nA + 2 * npp + o)) | (((R0[(o, 0)] >> h) & 1) << nU)
                                   for o in range(nb)]
                    return rowc[h]

                def rref(blocks):
                    """reduced echelon form of the class system of these blocks: pivot bit -> row; None if 0 = 1."""
                    piv = {}
                    for h in blocks:
                        for r in rows_of(h):
                            while True:
                                low = r & MU
                                if not low:
                                    break
                                e = piv.get(low.bit_length() - 1)
                                if e is None:
                                    break
                                r ^= e
                            if r & MU:
                                piv[(r & MU).bit_length() - 1] = r
                            elif r:
                                return None, 0
                    pm = 0
                    for hb in sorted(piv):               # back substitution, lowest pivot first
                        r = piv[hb]
                        for b in bits(r & pm):
                            r ^= piv[b]
                        piv[hb] = r
                        pm |= 1 << hb
                    return piv, pm

                def fullred(r, piv, pm):
                    for b in bits(r & pm):
                        r ^= piv[b]
                    return r

                def shared_file(idx):
                    return os.path.join(a.lshare, f"cls_s{pseed}_{idx}.json")

                def shared_get(idx):                 # TP: another process found this class already
                    if a.lshare and os.path.exists(shared_file(idx)):
                        try:
                            return [int(h) for h in json.load(open(shared_file(idx)))]
                        except (ValueError, OSError):
                            return None
                    return None

                def shared_put(idx, got):            # TP: the first find wins (hard link = atomic)
                    if not a.lshare:
                        return got
                    tmp = shared_file(idx) + f".{os.getpid()}.tmp"
                    json.dump([int(h) for h in got], open(tmp, "w"))
                    try:
                        os.link(tmp, shared_file(idx))
                    except OSError:
                        pass
                    os.remove(tmp)
                    return [int(h) for h in json.load(open(shared_file(idx)))]

                def lucky(pool, n, idx=0):
                    """a class of n blocks from the pool, consistent by itself, even V parity: (n - 2)-base + a pair."""
                    for it in range(a.s5bases):
                        if a.lshare and it % 8 == 0:
                            g = shared_get(idx)
                            if g is not None:
                                return g
                        base = (lrng if a.tp else rng).sample(pool, n - 2)
                        piv, pm = rref(base)
                        if piv is None:
                            continue
                        cls_stat["tries"] += 1
                        pb = parV(sum(1 << h for h in base))
                        cands = [h for h in pool if h not in base]
                        red = {h: [fullred(r, piv, pm) for r in rows_of(h)] for h in cands}
                        for i, h1 in enumerate(cands):
                            p1 = {}
                            bad = False
                            for r in red[h1]:
                                while True:
                                    low = r & MU
                                    if not low:
                                        break
                                    e = p1.get(low.bit_length() - 1)
                                    if e is None:
                                        break
                                    r ^= e
                                if r & MU:
                                    p1[(r & MU).bit_length() - 1] = r
                                elif r:
                                    bad = True
                                    break
                            if bad:
                                continue
                            v1 = pb ^ (Vset >> h1 & 1)
                            for h2 in cands[i + 1:]:
                                if v1 ^ (Vset >> h2 & 1):
                                    continue
                                p2 = {}
                                good2 = True
                                for r in red[h2]:
                                    while True:
                                        low = r & MU
                                        if not low:
                                            break
                                        hb = low.bit_length() - 1
                                        e = p1.get(hb)
                                        if e is None:
                                            e = p2.get(hb)
                                            if e is None:
                                                break
                                        r ^= e
                                    if r & MU:
                                        p2[(r & MU).bit_length() - 1] = r
                                    elif r:
                                        good2 = False
                                        break
                                if good2:
                                    return shared_put(idx, base + [h1, h2])
                    return None

                # sizes: class 0 = block 0 + n0 blocks (free), K0 helper classes
                cap_h = (2 * nA + 2 * npp + nb) // nb                 # free size of a helper class
                n0 = (2 * nA + 2 * npp) // nb                         # free size of class 0 (blocks besides block 0)
                need = len(realb) - n0 - K0 * cap_h                   # classes that take one block more
                if need > K0 - 1:
                    raise ValueError(f"S5: {need} lucky classes needed for {K0} helpers")
                pool = list(realb)
                rng.shuffle(pool)
                newc = []
                for i in range(max(0, need)):
                    got = lucky(pool, cap_h + 1, i)
                    if got is None:
                        raise ValueError("S5: no lucky class")
                    cls_stat["found"] += 1
                    pool = [h for h in pool if h not in got]
                    newc.append(got)
                    vprint(f"   S5: lucky class {i + 1} of {need} ({cap_h + 1} blocks) after {cls_stat['tries']} bases, "
                           f"{time.time() - t0:.0f} s", flush=True)
                okrest = False
                for _ in range(200):                                   # the free classes: even V parity each
                    rng.shuffle(pool)
                    rest, p0 = [], 0
                    for i in range(K0 - len(newc)):
                        take = pool[p0:p0 + cap_h]
                        p0 += len(take)
                        rest.append(take)
                    if len(pool) - p0 > n0:
                        continue
                    if all(c and not parV(sum(1 << h for h in c)) for c in rest):
                        okrest = True
                        break
                if not okrest:
                    raise ValueError("S5: no even split of the free classes")
                j0 = len(helpers)
                for c in newc + rest:
                    helpers.append((0, sum(1 << h for h in c)))
                hl[0][:] = list(range(j0, len(helpers)))
                hg += [blockg(m) & dg[0] for (_, m) in helpers[j0:]]
                hfan += [0] * (len(helpers) - j0)
                for q in range(nA):
                    rF[q] += [0] * (len(helpers) - j0)
                    rG[q] += [0] * (len(helpers) - j0)
                U0 = unknowns(0, True)
                for kq, (_, _, c1, c2, sx) in enumerate(pp):
                    for sd, colv in ((0, c1), (1, c2)):
                        U0.append((("pspill", kq, sd, 0), colv, sx))
                        for j in hl[0]:
                            U0.append((("pread", kq, sd, j), helpers[j][1] & colv, sx))
                E0 = list(range(nb))
                sol0, rk0 = solve(U0, 0, E0)
                if sol0 is None:
                    raise ValueError(f"S5: phase 2 not solvable (rank {rk0})")
                psp = [[0, 0] for _ in pp]
                prd = [[[], []] for _ in pp]
                for (kind, q, sd, idx), _, _ in sol0:
                    if kind == "spill":
                        (sF if sd == 0 else sG)[q][0] ^= 1
                    elif kind == "read":
                        (rF if sd == 0 else rG)[q][idx] ^= 1
                        nread += 1
                    elif kind == "fan":
                        hfan[q] ^= 1 << idx
                    elif kind == "pspill":
                        psp[q][sd] ^= 1
                    else:
                        prd[q][sd].append(idx)
                        nread += 1
                for q in range(nA):
                    if sF[q][0]:
                        F2[q] ^= dg[0]
                    if sG[q][0]:
                        G2[q] ^= dg[0]
                    for j in hl[0]:
                        x, y = rF[q][j], rG[q][j]
                        if x:
                            F2[q] ^= hg[j]
                        if y:
                            G2[q] ^= hg[j]
                        if (sF[q][0] & y) ^ (sG[q][0] & x) ^ (x & y):
                            hfan[j] ^= S[q]
                allb = sum(1 << h for h in range(NH))
                for kq, (_, _, c1, c2, sx) in enumerate(pp):
                    u1 = allb if psp[kq][0] else 0
                    u2 = allb if psp[kq][1] else 0
                    for j in prd[kq][0]:
                        u1 ^= helpers[j][1]
                    for j in prd[kq][1]:
                        u2 ^= helpers[j][1]
                    pp13[kq] = ((c2 ^ u1) & (c1 ^ u2)) ^ (c1 & c2)
                    for j in hl[0]:
                        x, y = int(j in prd[kq][0]), int(j in prd[kq][1])
                        if (psp[kq][0] & y) ^ (psp[kq][1] & x) ^ (x & y):
                            hfan[j] ^= sx
                    ppr[kq] = (list(prd[kq][0]), list(prd[kq][1]), psp[kq][0], psp[kq][1])
                used = list(range(len(helpers)))
                extra = [(hg[j], hfan[j]) for j in used]
                E[0] = E0
            extra2 = list(extra)
            for kq, (i1, i2, c1, c2, s) in enumerate(pp):
                g = blockg(c1) & dg[i1]
                if i2 is not None:
                    g ^= blockg(c2) & dg[i2]
                g ^= blockg(pp13[kq]) & dg[0]
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
            oddh = sum((parV(m) if isinstance(m, int) else parV(m[0]) + parV(m[1])) for (_, m) in helpers) if wh else 0
            vprint(f"   C-4 conditions for the 68-AND cache: stage-A forms outside H_V {oddA}; helper classes outside "
                   f"{oddh}; point-stage forms outside {oddp}; WIRES whose pure block part is outside H_V {oddw} "
                   f"(not fixed by the point stage: {fails})", flush=True)
            nh = len(used)
            tot = nA + nh + len(pp)
            ncell = sum(len(e) for e in E)
            vprint(f"{label}: cells emptied per point {[len(e) for e in E]} = {ncell} of {nD * nb}; point stage "
                  f"{len(pp)} ANDs for {nD * nb - ncell} live cells; helpers {nh}; TOTAL {tot} = {nA} + {nh} + {len(pp)}; "
                  f"kept reads {nread}; PLAN CHECKED on "
                  f"{'the real addresses (rows < ' + str(a.nreal) + ')' if a.dc else 'all 77 x 16 addresses'} "
                  f"(every output pure): {ok}", flush=True)
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
                               helpers=[[(D[ip] if isinstance(ip, int) else [D[i] for i in ip]),
                                         (str(m) if isinstance(m, int) else [str(x) for x in m]), str(hfan[j])]
                                        for j, (ip, m) in enumerate(helpers)],      # TP: the pair helper
                               stageA=stA,
                               point=[[D[i1], str(c1), (None if i2 is None else D[i2]), str(c2), str(s),
                                       ppr[kq][0], ppr[kq][1], ppr[kq][2], ppr[kq][3]]
                                      for kq, (i1, i2, c1, c2, s) in enumerate(pp)], s5=bool(a.s5), D0=D[0]),
                          open(a.plan, "w"))
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
                    r = attempt(K, pseed, rot, bool(a.s5))   # S5: one verbose pass (the class search is the cost)
                except ValueError as ex:                # C3's fix: a partition without a parity swap is skipped
                    print(f"   skipped: partition seed {pseed}, helpers per point {rot}: {ex}", flush=True)
                    continue
                t = r[f"K {K} helpers"]["total"]
                if best is None or t < best[0]:
                    best = (t, pseed, rot, r)
            if a.s5 and best is not None:               # S5: the first seed with a plan is the plan
                break
        if best is None:                                # S5: every seed failed: say so and stop cleanly
            print(f"NO PLAN k {a.k} D {D} K {K}: every partition seed {a.seed} .. {a.seed + a.pseeds - 1} was skipped; "
                  f"give other seeds", flush=True)
            sys.exit(3)
        t, pseed, rot, r = best
        if not a.s5:
            attempt(K, pseed, rot, True)
        out[K] = dict(total=t, pseed=pseed, rot=rot, res=r)
        if a.cls:
            print(f"R62 class-by-class choice: classes found {cls_stat['found']}, consistency tests {cls_stat['tries']}, "
                  f"partitions given up {cls_stat['fallback']}", flush=True)
        print(f"SUMMARY k {a.k} D {D} K {K}: control {ctl}; stage A {nA}; spill only "
              f"{r['spill only']['total']}; with helpers {t} (partition seed {pseed}, helpers per point {rot}); "
              f"{time.time() - t0:.0f} s", flush=True)
    if a.json:
        json.dump(dict(k=a.k, D=D, control=ctl, stageA=nA, out=out), open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
