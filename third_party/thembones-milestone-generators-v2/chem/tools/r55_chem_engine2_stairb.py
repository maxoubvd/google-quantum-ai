"""R55-CHEM-ENGINE2: the STAIRCASE on the BUILDERS' chain (searched basis, protected search), exact, emitted.

Every stage is planned by C2's engine, read-only (chem_ld2_basis.transform with the entry's basis file +
chem_ld2_prot.search_protected + Engine.plan); between stages the spill system of r55_chem_engine2_lib is solved.
The don't-cares of block 76 are fixed FIRST by this tool (the top row of the mixed tensor is set to a sum of other
rows that agrees with it on the real addresses: one global dependency, which every stage and the spill keep for
free), so the engines run with free_top off on a fully specified table.

  CHEM_BROT=7 python -I r55_chem_engine2_stairb.py TABLE.npz --k 7 --pfile k7_basis.json
        --stages "0,1,2|3,4" --kills 1:3:2 [--effort 1] [--seed 0] [--emit out.kmx] [--json out.json] [--ctl]
  --kills per later stage: emptied outputs : rank kills : extra designed relations per block: 0 none; 1, 2 = rank-one
  relations (one more output on a sum of two coordinates); 11, 12, 13 = one, two, three PLANAR relations
  e_j1 (x) e_oB + e_j2 (x) e_oC on two outputs oB, oC shared by all blocks of the stage; 21..26 = that many
  relations of a rank-distance-2 code (no rank-one combination) on THREE shared outputs.
  --engine plain uses chem_loaders_core.BlockDecomp on the same transformed table instead of C2's search.
Needs numba for C2's engine (box: ~/qcc/.venv_work/bin/python).
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
import chem_own_compile as cc  # noqa: E402
import chem_ld2_wwidth  # noqa: E402,F401
import chem_ld2_basis as lb2  # noqa: E402
import r55_chem_engine2_lib as L  # noqa: E402

bits = core.bits
par = L.par


def pad_top(an):
    """X0 with the top row replaced by a sum of other rows that agrees with it outside the don't-care space."""
    X = an.X0.copy()
    nr, nb = X.shape
    dl = lb2.dont_care_low(an)
    if not dl:
        return X, None
    phis = core.nullspace(dl, an.nl)              # functionals that see only the fixed part of the top row
    sp = core.Span()
    for M in range(nr - 1):
        v = 0
        for o in range(nb):
            for i, ph in enumerate(phis):
                if par(ph & int(X[M, o])):
                    v |= 1 << (len(phis) * o + i)
        sp.add(v, 1 << M)
    t = 0
    for o in range(nb):
        for i, ph in enumerate(phis):
            if par(ph & int(X[nr - 1, o])):
                t |= 1 << (len(phis) * o + i)
    r, c = sp.reduce(t)
    if r:
        return X, None
    new = np.zeros(nb, np.uint32)
    for M in bits(c):
        new ^= X[M]
    dsp = core.Span()
    for v in dl:
        dsp.add(v)
    for o in range(nb):
        assert dsp.reduce(int(new[o]) ^ int(X[nr - 1, o]))[0] == 0, "padding touches a real address"
    X[nr - 1] = new
    return X, c


def plan_stage(an2, X, blocks, seed, effort, engine, warm=None, keep=None):
    """warm = (prot, start_R) per block of this stage, taken from a joint plan (C2's engine only);
    keep: a dict that receives the engine's spans (R, prot) for later warm starts."""
    a = copy.copy(an2)
    a.X0 = X.copy()
    a.blocks_override = [list(b) for b in blocks]
    if engine == "plain":
        best = None
        for sd in range(seed, seed + 3):
            bd = core.BlockDecomp(a, sd, 3)
            try:
                p = bd.run()
            except (RuntimeError, AssertionError):
                continue
            if best is None or len(p) < len(best[0]):
                best = ([list(x) for x in p], len(p), list(bd.evec))
        if best is None:
            raise RuntimeError("plain engine found no plan")
        return best
    import chem_ld2_prot as pr
    if warm is not None:
        E, price, hist = pr.search_protected(a, seed, effort, prot=[list(x) for x in warm[0]],
                                             start_R=[list(x) for x in warm[1]])
    else:
        E, price, hist = pr.search_protected(a, seed, effort)
    if keep is not None:
        keep["R"] = [list(x) for x in E.R]
        keep["prot"] = [list(x) for x in E.prot]
    prods = E.plan()
    return [list(p) for p in prods], int(price), list(E.codim)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, default=7)
    ap.add_argument("--pfile", required=True)
    ap.add_argument("--stages", default="0,1,2|3,4")
    ap.add_argument("--kills", default="1:3:0")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--effort", type=int, default=1)
    ap.add_argument("--zseed", type=int, default=1)
    ap.add_argument("--tries", type=int, default=12)
    ap.add_argument("--engine", default="c2", help="c2 | plain")
    ap.add_argument("--ctl", action="store_true", help="also run the builders' own command (free_top) as control")
    ap.add_argument("--noctl", action="store_true", help="skip the control on the padded table")
    ap.add_argument("--warmA", action="store_true", help="first stage starts from the joint plan's relations of "
                    "its blocks (needs the control run, C2's engine)")
    ap.add_argument("--emit", default=None)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    t0 = time.time()
    if a.engine == "c2":
        import chem_ld2_engine as e2
        e2.install()
    words, _ = cc.load_words(a.table)
    lbs = cc.feasible_lbs(words)
    an, lb = cc.analyze_first(words, a.k, lbs, False)
    bf = json.load(open(a.pfile))
    P = [int(x, 16) for x in bf["P"]]
    five = bf.get("blocks", [[3 * i, 3 * i + 1, 3 * i + 2] for i in range(5)])
    nr, nb = an.X0.shape
    out = dict(k=a.k, stages=a.stages, kills=a.kills, seed=a.seed, effort=a.effort, engine=a.engine)
    if a.ctl and a.engine == "c2":
        import chem_ld2_prot as pr
        an2c = lb2.transform(an, P, free_top=True, blocks=five)
        E, price, hist = pr.search_protected(an2c, a.seed, a.effort)
        print(f"control (builders' command: free_top, protected): {price} products, codim {E.codim} "
              f"({time.time() - t0:.0f} s)", flush=True)
        out["control_free_top"] = int(price)
    Xp, V = pad_top(an)
    an_p = copy.copy(an)
    an_p.X0 = Xp
    an2 = lb2.transform(an_p, P, free_top=False, blocks=five)
    X2 = an2.X0.copy()
    print(f"k {a.k}: outputs {nb}, rows {nr}, top row padded: {V is not None}", flush=True)
    nctl = None
    if not a.noctl:
        joint = {}
        pc, pricec, codc = plan_stage(an2, X2, five, a.seed, a.effort, a.engine, keep=joint)
        nctl = len(pc)
        print(f"control on the padded table, all five blocks: {len(pc)} products, codim {codc} "
              f"sum {sum(codc)} ({time.time() - t0:.0f} s)", flush=True)
        out["control_padded"] = len(pc)
    stages = [[int(x) for x in st.split(",")] for st in a.stages.split("|")]
    kills = [tuple(int(x) for x in kk.split(":")) for kk in a.kills.split("|")]
    rng = random.Random(a.zseed)
    prods, feats, report = [], [], []
    Xres = X2.copy()
    for si, st in enumerate(stages):
        blks = [five[i] for i in st]
        coords = [t for b in blks for t in b]
        m = np.uint32(sum(1 << t for t in coords))
        note = ""
        if si:
            kk = kills[si - 1]
            kout, krank, nextra = kk[0], kk[1], (kk[2] if len(kk) > 2 else 0)
            sp = None
            for attempt in range(a.tries):
                pick = rng.sample(range(nb), kout + 2 + 2 * len(blks))
                killed = sorted(pick[:kout])
                if attempt < 3 and kout:                     # first tries: the outputs held by the most fan-outs
                    cover = sorted(range(nb), key=lambda o: -sum(1 for f in feats if f[3] >> o & 1))
                    killed = sorted(cover[attempt:attempt + kout])
                    pick = killed + [o for o in pick if o not in killed]
                oB, oC = pick[kout], pick[kout + 1]            # the two outputs shared by the planar relations
                extra = []
                for bi, b in enumerate(blks):
                    if nextra in (1, 2):                       # rank-one: one more output on a sum of coordinates
                        extra.append([(b[0], pick[kout + 2 + 2 * bi]), (b[1], pick[kout + 2 + 2 * bi])])
                    if nextra == 2:
                        extra.append([(b[-2], pick[kout + 3 + 2 * bi]), (b[-1], pick[kout + 3 + 2 * bi])])
                    if 11 <= nextra <= 13:                     # planar relations on the common outputs oB, oC
                        extra.append([(b[0], oB), (b[1], oC)])
                    if 12 <= nextra <= 13:
                        extra.append([(b[1], oB), (b[2], oC)])
                    if nextra in (13,):
                        extra.append([(b[2], oB), (b[0], oC)])
                    if 21 <= nextra <= 26:                     # a rank-distance-2 code on three shared outputs
                        outs3 = [o for o in pick if o not in killed][:3]
                        extra += L.mrd_relations(b, outs3, nextra - 20)
                Z, zs = [], core.Span()
                while len(Z) < krank:
                    z = rng.getrandbits(nr)
                    if z and zs.add(z):
                        Z.append(z)
                sp, neq, rk = L.solve_spill(Xres, feats, coords, killed, Z, extra, nr, nb)
                if sp is not None:
                    break
            if sp is None:
                print(f"stage {si}: SPILL SYSTEM HAS NO SOLUTION in {a.tries} tries ({neq} equations, "
                      f"{len(feats) * len(coords)} unknowns, rank {rk})", flush=True)
                return
            L.apply_spill(Xres, prods, feats, sp, nr)
            L.check_spill(Xres, coords, killed, Z, extra, nr)
            note = (f"spill {neq} equations rank {rk} on {len(feats) * len(coords)} unknowns (try {attempt}; emptied "
                    f"{killed}, extra {extra}); ")
        warm = None
        if si == 0 and a.warmA and a.engine == "c2" and not a.noctl:
            warm = ([joint["prot"][i] for i in st], [joint["R"][i] for i in st])
        try:
            p, price, cod = plan_stage(an2, Xres & m, blks, a.seed, a.effort, a.engine, warm=warm)
        except Exception as exc:  # noqa: BLE001
            print(f"stage {si}: {note}the engine failed: {type(exc).__name__}: {exc}", flush=True)
            return
        assert np.array_equal(L.tensor_of(p, nr, nb), Xres & m), f"stage {si} plan is not exact"
        n0 = len(prods)
        prods += p
        feats += L.feats_of(p, n0)
        Xres &= ~m
        report.append(dict(stage=st, products=len(p), codim=cod))
        print(f"stage {si} blocks {st}: {note}{len(p)} products, codim {cod} sum {sum(cod)} "
              f"({time.time() - t0:.0f} s)", flush=True)
    ok = np.array_equal(L.tensor_of(prods, nr, nb), X2)
    print(f"stages {a.stages} kills {a.kills} engine {a.engine} seed {a.seed}: TOTAL {len(prods)} = "
          f"{' + '.join(str(r['products']) for r in report)} (control padded {nctl}); "
          f"mixed tensor identity: {ok}; {time.time() - t0:.0f} s", flush=True)
    out.update(total=len(prods), tensor_ok=bool(ok), report=report)
    if ok and a.emit:
        old = lb2.map_back(prods, P)
        text, st = core.emit(an, old)
        with open(a.emit, "w") as fh:
            fh.write(text)
        print("EMITTED with chem_loaders_core.emit (every wire on every address):",
              dict(toffoli=st["toffoli"], qubits=st["qubits"], products=st["products"], infra=st.get("infra")),
              flush=True)
        out["emit"] = dict(toffoli=st["toffoli"], qubits=st["qubits"], products=st["products"])
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
