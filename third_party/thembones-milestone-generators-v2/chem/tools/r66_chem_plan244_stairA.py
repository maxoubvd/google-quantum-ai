"""R66-CHEM-PLAN244 (second holder): STAGE A OF THE KEPT LOADER AS A STAIRCASE ("SA2").

= the set-up of tools/chem_c3_c2a.py (builder C2's basis of the loader with ONE block replaced by the indicators of
the three low points D; the padded table; C2's engine, imported read-only) + the stage loop of
tools/r55_chem_engine2_stairb.py (R55's staircase: the features of the earlier stage carry spill bits in the later
stage's coordinates; one linear solve empties outputs and kills ranks there) on the FOUR blocks outside D.
The kept stage's stage A has always been planned as four blocks at once (149 = codim [16, 15, 15, 16]); a
staircase was never run inside it (R58-CHEM-KEPT-reprice row 9).

Output: a stage-A plan in the format of tools/r58_chem_kept_c2a.py (for tools/r66_chem_plan244_mask3.py --planA):
exact on the 13 low points outside D, free on D.

  CHEM_BROT=7 CHEM_HB=77 python -I r66_chem_plan244_stairA.py PAD.npz --k 1 --pfile k1bs1_D13-14-15.json --rep 4
        --stages "0,1,2|3" --kills 2:2 [--lucky 1] [--seed 0] [--zseed 1] [--effort 1] [--tries 12] --out planA.json
  --stages: indices into the four blocks outside D (in the basis file's order, the replaced block left out);
  --kills per later stage: emptied outputs : rank kills : extra (as in r55_chem_engine2_stairb.py).
  --lucky N : LUCKY RANK KILLS (new here). R55's budget law stops the linear kills where the spill unknowns of a
              coordinate (one per earlier feature: 224) are used up: two emptied outputs + two kills = 208-210.
              The unknowns a coordinate that are left are a free affine space of spills; a block's flattening
              (rows x live cells) of rank r drops to r - 1 on a random point of it with probability
              2^-(cells - r + 1) (2^-12 at 84 cells, rank 73). The tool walks that space (one free vector flipped a
              step) until every block of the stage has lost N more ranks. No equation is added: the spill stays a
              solution of the linear system.
Needs numba (box: ~/qcc/.venv_work/bin/python).
"""
import argparse
import copy
import json
import os
import random
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_wwidth  # noqa: E402,F401
import chem_ld2_basis as lb2  # noqa: E402
import chem_ld2_engine as e2  # noqa: E402
import r55_chem_engine2_lib as L  # noqa: E402

bits = core.bits


def plan_stage(an2, X, blocks, seed, effort, fill=0, double=False):
    import chem_ld2_prot as pr
    a = copy.copy(an2)
    a.X0 = X.copy()
    a.blocks_override = [list(b) for b in blocks]
    if double:
        E, price, _ = pr.search_protected_double(a, seed, effort, fill=fill)
    elif fill:
        E, price, _ = pr.search_protected(a, seed, effort, fill=fill)
    else:
        E, price, _ = pr.search_protected(a, seed, effort)
    return [list(p) for p in E.plan()], int(price), list(E.codim)


def spill_null(X, feats, coords, killed, Z, nr, nb):
    """L.solve_spill without extra relations, coordinate by coordinate (the coordinates do not share unknowns);
    returns per coordinate (particular solution, basis of the free solutions, rank), masks over the features."""
    live_mask = ((1 << nb) - 1) & ~sum(1 << o for o in killed)
    nk, nz = len(killed), len(Z)
    zpar = [[L.par(z & H) for z in Z] for (_, _, H, _) in feats]
    out = []
    for j in coords:
        sp, null = core.Span(), []
        for f, (_, _, H, s) in enumerate(feats):
            col = 0
            for ki, o in enumerate(killed):
                if s >> o & 1:
                    col |= H << (ki * nr)
            sl = s & live_mask
            if sl:
                for zi in range(nz):
                    if zpar[f][zi]:
                        col |= sl << (nk * nr + zi * nb)
            if not col:
                null.append(1 << f)
                continue
            r, c = sp.reduce(col, 1 << f)
            if r:
                sp.piv[r.bit_length() - 1] = (r, c)
            else:
                null.append(c)
        cj = [L.colint(X, o, j) for o in range(nb)]
        tgt = 0
        for ki, o in enumerate(killed):
            tgt |= cj[o] << (ki * nr)
        for zi, z in enumerate(Z):
            m = 0
            for o in range(nb):
                if live_mask >> o & 1 and L.par(z & cj[o]):
                    m |= 1 << o
            tgt |= m << (nk * nr + zi * nb)
        r, c = sp.reduce(tgt)
        if r:
            return None
        out.append((c, null, len(sp)))
    return out


def rank_rows(rows):
    piv = {}
    for v in rows:
        while v:
            b = v.bit_length() - 1
            e = piv.get(b)
            if e is None:
                piv[b] = v
                break
            v ^= e
    return len(piv)


def lucky_block(Xres, feats, bcoords, nulls, killed, nr, nb, want, rng, nsamples):
    """Walk the free spills of one block's coordinates until its flattening has lost `want` ranks.
    nulls: {coordinate: basis of the free spills (masks over the features)}; the particular solutions are ALREADY
    applied to Xres. Returns ({coordinate: mask over features to add} or None, rank before, best rank, samples,
    free vectors)."""
    live = [o for o in range(nb) if o not in killed]
    cell = {(o, j): i for i, (o, j) in enumerate((o, j) for j in bcoords for o in live)}
    rows = [0] * nr
    for M in range(nr):
        v = 0
        for (o, j), i in cell.items():
            if int(Xres[M, o]) >> j & 1:
                v |= 1 << i
        rows[M] = v
    eff, who = [], []
    for j in bcoords:
        for n in nulls[j]:
            e = [0] * nr
            for f in bits(n):
                _, _, H, s = feats[f]
                cm = 0
                for o in bits(s):
                    i = cell.get((o, j))
                    if i is not None:
                        cm |= 1 << i
                if cm:
                    for M in bits(H):
                        e[M] ^= cm
            if any(e):
                eff.append(e)
                who.append((j, n))
    r0 = rank_rows(rows)
    if want <= 0:
        return {j: 0 for j in bcoords}, r0, r0, 0, len(eff)
    state = [0] * len(eff)
    best = r0
    for it in range(nsamples):
        i = rng.randrange(len(eff))
        state[i] ^= 1
        e = eff[i]
        rows = [x ^ y for x, y in zip(rows, e)]
        r = rank_rows(rows)
        best = min(best, r)
        if r <= r0 - want:
            add = {j: 0 for j in bcoords}
            for q, on in enumerate(state):
                if on:
                    add[who[q][0]] ^= who[q][1]
            return add, r0, r, it + 1, len(eff)
    return None, r0, best, nsamples, len(eff)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pad")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--pfile", required=True)
    ap.add_argument("--rep", type=int, required=True, help="the block of the basis replaced by the point block")
    ap.add_argument("--D", default="13,14,15")
    ap.add_argument("--stages", default="0,1,2|3")
    ap.add_argument("--kills", default="2:2")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--zseed", type=int, default=1)
    ap.add_argument("--effort", type=int, default=1)
    ap.add_argument("--tries", type=int, default=12)
    ap.add_argument("--fill", type=int, default=0, help="first stage: start every block at the usual codim + fill")
    ap.add_argument("--double", action="store_true", help="first stage: search_protected_double")
    ap.add_argument("--lucky", type=int, default=0, help="extra rank kills a block by sampling the free spills")
    ap.add_argument("--lsamples", type=int, default=200000, help="samples a block for the lucky kills")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    t0 = time.time()
    core.J = 5
    words = np.asarray(np.load(a.pad)["angle_words"], np.int64)
    an = core.analyze(words, a.k, 4)
    nl = an.nl
    D = [int(x) for x in a.D.split(",")]
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
    sp0 = core.Span()
    for t in range(nl):
        sp0.add(LB[t], 1 << t)
    mD = []
    for p in D:
        r, c = sp0.reduce(1 << p)
        assert r == 0
        mD.append(c)
    p0 = json.load(open(a.pfile))
    P = [int(x, 16) for x in p0["P"]]
    five = [list(b) for b in p0["blocks"]]
    P2 = list(P)
    for i, t in enumerate(five[a.rep]):
        P2[t] = mD[i]
    try:
        lb2.mat_inv(P2, nl)
    except ValueError:
        print(f"k {a.k} rep {a.rep} D {D}: the point block is not a complement of the other four blocks")
        sys.exit(3)
    e2.install()
    an2 = lb2.transform(an, P2, free_top=False, blocks=five)
    four = [five[b] for b in range(5) if b != a.rep]
    m12 = np.uint32(sum(1 << t for b in four for t in b))
    X2 = an2.X0.copy()
    nr, nb = X2.shape
    stages = [[int(x) for x in st.split(",")] for st in a.stages.split("|")]
    assert sorted(i for st in stages for i in st) == [0, 1, 2, 3], "the stages must partition the four blocks"
    kills = [tuple(int(x) for x in kk.split(":")) for kk in a.kills.split("|")] if len(stages) > 1 else []
    rng = random.Random(a.zseed)
    prods, feats, report = [], [], []
    Xres = (X2 & m12).copy()
    tag = f"k {a.k} stages {a.stages} kills {a.kills} lucky {a.lucky} seed {a.seed} zseed {a.zseed}"
    for si, st in enumerate(stages):
        blks = [four[i] for i in st]
        coords = [t for b in blks for t in b]
        m = np.uint32(sum(1 << t for t in coords))
        note = ""
        if si:
            kk = kills[min(si - 1, len(kills) - 1)]
            kout, krank, nextra = kk[0], kk[1], (kk[2] if len(kk) > 2 else 0)
            sp = None
            neq = rk = 0
            nf = len(feats)
            for attempt in range(a.tries):
                pick = rng.sample(range(nb), kout + 2 + 2 * len(blks))
                killed = sorted(pick[:kout])
                if attempt < 3 and kout:                     # first tries: the outputs held by the most fan-outs
                    cover = sorted(range(nb), key=lambda o: -sum(1 for f in feats if f[3] >> o & 1))
                    killed = sorted(cover[attempt:attempt + kout])
                    pick = killed + [o for o in pick if o not in killed]
                oB, oC = pick[kout], pick[kout + 1]
                extra = []
                for bi, b in enumerate(blks):
                    if nextra in (1, 2):
                        extra.append([(b[0], pick[kout + 2 + 2 * bi]), (b[1], pick[kout + 2 + 2 * bi])])
                    if nextra == 2:
                        extra.append([(b[-2], pick[kout + 3 + 2 * bi]), (b[-1], pick[kout + 3 + 2 * bi])])
                    if 11 <= nextra <= 13:
                        extra.append([(b[0], oB), (b[1], oC)])
                    if 12 <= nextra <= 13:
                        extra.append([(b[1], oB), (b[2], oC)])
                    if nextra in (13,):
                        extra.append([(b[2], oB), (b[0], oC)])
                Z, zs = [], core.Span()
                while len(Z) < krank:
                    z = rng.getrandbits(nr)
                    if z and zs.add(z):
                        Z.append(z)
                if a.lucky:
                    assert not extra, "--lucky runs without extra relations"
                    sols = spill_null(Xres, feats, coords, killed, Z, nr, nb)
                    neq = len(coords) * (kout * nr + krank * nb)
                    if sols is None:
                        sp = None
                        continue
                    rk = sum(x[2] for x in sols)
                    sp = [0] * nf
                    for ci, (c, _, _) in enumerate(sols):
                        for f in bits(c):
                            sp[f] |= 1 << coords[ci]
                    Xtry = Xres.copy()
                    ptry = [list(q) for q in prods]
                    L.apply_spill(Xtry, ptry, feats, sp, nr)
                    okl = True
                    for b in blks:
                        nulls = {j: sols[coords.index(j)][1] for j in b}
                        add, r0, r1, ns, nfree = lucky_block(Xtry, feats, b, nulls, killed, nr, nb, a.lucky, rng,
                                                             a.lsamples)
                        print(f"   lucky kills (try {attempt}, emptied {killed}), block {b}: {nfree} free spill "
                              f"vectors; rank {r0} -> {r1} after {ns} samples ({time.time() - t0:.0f} s)",
                              flush=True)
                        if add is None:
                            okl = False
                            break
                        sp2 = [0] * nf
                        for j, msk in add.items():
                            for f in bits(msk):
                                sp2[f] |= 1 << j
                        L.apply_spill(Xtry, ptry, feats, sp2, nr)
                        for f in range(nf):
                            sp[f] ^= sp2[f]
                    if not okl:
                        sp = None
                        continue
                    break
                sp, neq, rk = L.solve_spill(Xres, feats, coords, killed, Z, extra, nr, nb)
                if sp is not None:
                    break
            if sp is None:
                print(f"stage {si}: SPILL SYSTEM HAS NO SOLUTION in {a.tries} tries ({neq} equations, "
                      f"{nf * len(coords)} unknowns, rank {rk})", flush=True)
                print(f"RESULT {tag}: NO SPILL", flush=True)
                return
            L.apply_spill(Xres, prods, feats, sp, nr)
            L.check_spill(Xres, coords, killed, Z, extra, nr)
            note = (f"spill {neq} equations rank {rk} on {nf * len(coords)} unknowns (try {attempt}; emptied "
                    f"{killed}); ")
        try:
            p, price, cod = plan_stage(an2, Xres & m, blks, a.seed, a.effort,
                                       fill=a.fill if si == 0 else 0, double=a.double and si == 0)
        except Exception as exc:  # noqa: BLE001
            print(f"stage {si}: {note}the engine failed: {type(exc).__name__}: {exc}", flush=True)
            print(f"RESULT {tag}: ENGINE FAILED at stage {si} (after "
                  f"{' + '.join(str(r['products']) for r in report) or '-'})", flush=True)
            return
        assert np.array_equal(L.tensor_of(p, nr, nb), Xres & m), f"stage {si} plan is not exact"
        n0 = len(prods)
        prods += p
        feats += L.feats_of(p, n0)
        Xres &= ~m
        report.append(dict(stage=st, products=len(p), codim=[int(x) for x in cod]))
        print(f"stage {si} blocks {st}: {note}{len(p)} products, codim {list(cod)} ({time.time() - t0:.0f} s)",
              flush=True)
    ten = L.tensor_of(prods, nr, nb)
    ok = bool(np.array_equal(ten & m12, X2 & m12))
    print(f"RESULT {tag} pfile {os.path.basename(a.pfile)} rep {a.rep} D {D}: STAGE A {len(prods)} = "
          f"{' + '.join(str(r['products']) for r in report)}, codim {[r['codim'] for r in report]}; exact on the "
          f"12 coordinates outside D: {ok}; {time.time() - t0:.0f} s", flush=True)
    if ok and a.out:
        old = [[int(x) for x in q] for q in lb2.map_back(prods, P2)]
        json.dump(dict(k=a.k, D=D, rep=a.rep, pfile=os.path.basename(a.pfile), seed=a.seed, effort=a.effort,
                       products=len(old), codim=[int(x) for r in report for x in r["codim"]],
                       stages=a.stages, kills=a.kills, lucky=a.lucky, zseed=a.zseed,
                       prods=[[str(x) for x in q] for q in old]), open(a.out, "w"))


if __name__ == "__main__":
    main()
