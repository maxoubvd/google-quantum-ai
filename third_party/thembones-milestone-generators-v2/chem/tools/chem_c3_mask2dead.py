"""chem_c3_mask2dead.py: builder C3 (Nov 17 13:45Z): a THREE-WAY MERGE (diff3) of R62-CHEM-BILIN's
r62_chem_bilin_mask2.py (S5, marked S5 / R62) and R62-CHEM-CACHE's re-pricer r62_chem_cache_rp_mask.py
(dead words, marked R62RP), both copies of r58_chem_kept_rp_mask.py; no line of either is changed. For the
CACHE re-pricer's STACKING TEST (S5 + --dead). Their two headers follow.
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
R62-CHEM-CACHE re-pricer's copy of tools/r58_chem_kept_rp_mask.py, unchanged except where marked R62RP:
  R62RP --dead : the table's own EXACT DON'T-CARE WORDS (tools/r62_chem_cache_dead.dead_mask: words the decoder never
              reads) become unknowns of the point stage. At point p of D and block h, the wires (j, b) of loader k
              whose word (k, j) or (k-1, j) of row 16 h + p is dead are free there (bits --deadbits; bit 0 only when
              the neighbour loader that shares the word is not wired, --wired). The free output vectors are the kernel
              of the LIVE wires' rows of G (outputs are basis outputs); each kernel vector x is one unknown with block
              vector 1 << h and fan-out x. A chosen unknown is a TABLE CHANGE at that address on dead wires only; it
              is applied to the residual, the live cells and the point stage are built on it, and the whole plan is
              re-checked on every run address. The changed (row, word, bit) list is printed (and --chg writes it).
  The solve, classes, C-4 conditions, check: as in the R58 re-pricer's tool.
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
    ap.add_argument("--dead", action="store_true", help="R62RP: the table's dead words are free at the points of D")
    ap.add_argument("--deadbits", default="1,2,3,4,5", help="R62RP: word bits a dead word frees (bit 0: see --wired)")
    ap.add_argument("--wired", default="4,6,8,9,10,12,13,14",
                    help="R62RP: wired loaders (bit 0 of a word they share is NOT free); G8F's eight by default")
    ap.add_argument("--deaduse", default="either", choices=["either", "own", "prev"],
                    help="R62RP: dead in batch k or k-1 (either), batch k only (own), batch k-1 only (prev)")
    ap.add_argument("--chg", default=None, help="R62RP: write the table changes of the best plan (json)")
    a = ap.parse_args()
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

    # R62RP: dead wires at the points of D -> kernel of the LIVE wires' rows of G = free output vectors there
    deadU = [[] for _ in range(nD)]             # per point: (h, x) with x a mask over the basis outputs
    deadfree = {}                               # (ip, h) -> set of free wire indices
    nwires = len(an.out_wires)
    if a.dead:
        from r62_chem_cache_dead import dead_mask
        edges = np.asarray(zt["tree_edges"], np.int64)
        dm = dead_mask(w[:a.nreal], edges, core.BROT, 0)
        wired = {int(x) for x in a.wired.split(",") if x != ""}
        dbits = [int(x) for x in a.deadbits.split(",") if x != ""]
        kk = a.k
        ncolw = w.shape[1]
        wpos = {wb: q for q, wb in enumerate(an.out_wires)}
        npairs, nfree = 0, 0
        for ip, p in enumerate(D):
            for h in range(NH):
                row = 16 * h + p
                if row >= a.nreal:
                    continue
                free = set()
                for j in range(5):
                    nbrs = []
                    if a.deaduse in ("either", "own") and 5 * kk + j < ncolw and dm[row, 5 * kk + j]:
                        nbrs.append(kk + 1)              # word (k, j) is also the input of loader k + 1
                    if a.deaduse in ("either", "prev") and kk >= 1 and dm[row, 5 * (kk - 1) + j]:
                        nbrs.append(kk - 1)              # word (k - 1, j) is also the output of loader k - 1
                    if not nbrs:
                        continue
                    npairs += 1
                    for b in dbits:
                        if (j, b) in wpos:
                            free.add(wpos[(j, b)])
                    if any(x not in wired for x in nbrs) and (j, 0) in wpos:
                        free.add(wpos[(j, 0)])          # bit 0 only when a sharing neighbour is not wired
                if not free:
                    continue
                deadfree[(ip, h)] = free
                nfree += len(free)
                ker = core.nullspace([an.G[q] for q in range(nwires) if q not in free], nb)
                for x in ker:
                    deadU[ip].append((h, x))
        print(f"R62RP: loader {kk}, D {D}: dead (row, word) pairs at D {npairs}, free wires {nfree}; dead unknowns per "
              f"point {[len(u) for u in deadU]} (use {a.deaduse}, bits {dbits}, wired {sorted(wired)}, outputs {nb}, "
              f"wires {nwires})", flush=True)

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
            for (h, x) in deadU[ip]:                    # R62RP: a table change at (h, D[ip]) along x
                U.append((("dead", h, x, ip), 1 << h, x))
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

        res = {}
        for label, cap, wh in (("spill only", cap0, False), (f"K {K} helpers", cap1, True)):
            if a.s5 and wh:                          # S5 phase 1: points 2, 3 with the same emptied outputs
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
            hg = [blockg(m) & dg[ip] for (ip, m) in helpers]
            hfan = [0] * len(helpers)
            sF = [[0] * nD for _ in range(nA)]
            sG = [[0] * nD for _ in range(nA)]
            rF = [[0] * len(helpers) for _ in range(nA)]
            rG = [[0] * len(helpers) for _ in range(nA)]
            nread = 0
            dsel = {}                                   # R62RP: (ip, h) -> output vector of the table change
            for sol in sols:
                for (kind, q, sd, idx), _, _ in sol:
                    if kind == "spill":
                        (sF if sd == 0 else sG)[q][idx] ^= 1
                    elif kind == "read":
                        (rF if sd == 0 else rG)[q][idx] ^= 1
                        nread += 1
                    elif kind == "dead":
                        dsel[(idx, q)] = dsel.get((idx, q), 0) ^ sd
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
            chg = []                                    # R62RP: the table changes (row, word slot j, bit b)
            dextra = []                                 # C3 (14:03Z): kept apart: S5's phase 2 rebuilds `extra`
            for (ip, h), x in sorted(dsel.items()):
                if not x:
                    continue
                extra.append((1 << (h * NL + D[ip]), x))
                dextra.append((1 << (h * NL + D[ip]), x))
                for q in range(nwires):
                    if bin(an.G[q] & x).count("1") & 1:
                        assert q in deadfree[(ip, h)], "R62RP: a table change on a LIVE wire"
                        chg.append((16 * h + D[ip], an.out_wires[q][0], an.out_wires[q][1]))
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

                def lucky(pool, n):
                    """a class of n blocks from the pool, consistent by itself, even V parity: (n - 2)-base + a pair."""
                    for _ in range(a.s5bases):
                        base = rng.sample(pool, n - 2)
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
                                    return base + [h1, h2]
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
                    got = lucky(pool, cap_h + 1)
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
                # C3 (13:57Z): phase 2 (the first point) WITHOUT the dead-word unknowns, as the CACHE re-pricer's note
                # says it may: S5's solution loop below reads every unknown that is not spill / read / fan / pspill as a
                # point product's read, so a ("dead", h, x, ip) unknown indexed prd[h] (IndexError on every seed)
                U0 = [u for u in unknowns(0, True) if u[0][0] != "dead"]
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
                # C3 (14:03Z): S5 rebuilds the helper corrections; the dead-word table changes of phase 1 stay in
                extra = [(hg[j], hfan[j]) for j in used] + dextra
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
            oddh = sum(parV(m) for (_, m) in helpers) if wh else 0
            vprint(f"   C-4 conditions for the 68-AND cache: stage-A forms outside H_V {oddA}; helper classes outside "
                   f"{oddh}; point-stage forms outside {oddp}; WIRES whose pure block part is outside H_V {oddw} "
                   f"(not fixed by the point stage: {fails})", flush=True)
            nh = len(used)
            tot = nA + nh + len(pp)
            ncell = sum(len(e) for e in E)
            vprint(f"{label}: cells emptied per point {[len(e) for e in E]} = {ncell} of {nD * nb}; point stage "
                  f"{len(pp)} ANDs for {nD * nb - ncell} live cells; helpers {nh}; TOTAL {tot} = {nA} + {nh} + {len(pp)}; "
                  f"kept reads {nread}; R62RP table changes {len(chg)} bits on dead wires in "
                  f"{len(set((r, j) for r, j, _ in chg))} (row, word) pairs; PLAN CHECKED on "
                  f"{'the real addresses (rows < ' + str(a.nreal) + ')' if a.dc else 'all 77 x 16 addresses'} "
                  f"(every output pure): {ok}", flush=True)
            res[label] = dict(emptied=[len(e) for e in E], point=len(pp), helpers=nh,
                              total=tot + 1000 * (oddw + oddp + oddh), reads=nread, checked=bool(ok),
                              c4=[oddA, oddh, oddp, oddw], changes=len(chg))
            if verbose and wh and a.chg:
                json.dump(dict(k=a.k, D=D, total=tot, changes=[list(map(int, c)) for c in chg]), open(a.chg, "w"))
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
