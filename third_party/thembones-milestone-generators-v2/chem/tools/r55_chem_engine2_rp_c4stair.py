"""R55-CHEM-ENGINE2 re-pricer: does the STAIRCASE stack with C2's C-4 sector closure (the 68-AND high cache)?

On C2's V-padded table (c4/<build>k13_pad.npz: block 76's don't-care rows already hold the V-consistent padding,
so no don't-cares are left and the designer's pad_top is a no-op), loader k, C2's own plan file for the basis:
  1. control: C2's engine on all five blocks (search_protected), extraction E.plan() (or sectored), closure seeds;
  2. the staircase exactly as tools/r55_chem_engine2_stairb.py (stage A by C2's engine, spill system of
     tools/r55_chem_engine2_lib.py, stage B by C2's engine), same extraction as 1;
  3. every plan: chem_loaders_core.emit (69-AND cache) and, for each closure seed, chem_ld2_c4.closure3(an, prods,
     V, seed); the first solved one is emitted by chem_ld2_c4.emit_c4 (68-AND cache; both emitters check every wire
     on every address). Reads only; C2's and the designer's files are imported, never written.
Usage: CHEM_BROT=7 CHEM_HB=77 python -I r55_chem_engine2_rp_c4stair.py PAD.npz PROJ.json PLAN.json --k K
       --stages "0,1,2|3,4" --kills 1:2:0 [--seed 0] [--zseed 1] [--sect 0] [--cseeds 8] [--noctl] [--json f]
       [--emit prefix]
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
import chem_ld2_c4 as c4  # noqa: E402
import r55_chem_engine2_lib as L  # noqa: E402


def plan_stage(an2, X, blocks, seed, effort, sect, xseed):
    import chem_ld2_prot as pr
    a = copy.copy(an2)
    a.X0 = X.copy()
    a.blocks_override = [list(b) for b in blocks]
    E, price, hist = pr.search_protected(a, seed, effort)
    if sect:
        prods, _ = c4.extract_sectored(E, rounds=sect, seed=xseed)
    else:
        prods = E.plan()
    return [list(p) for p in prods], int(price), list(E.codim)


def close(an, old, V, cseed0, ncs, tag, out, emit):
    """closure3 over ncs seeds; returns the 68-AND toffoli count or None; records the diagnostics."""
    rec = []
    for cs in range(cseed0, cseed0 + ncs):
        new, info = c4.closure3(an, [tuple(p) for p in old], V, seed=cs)
        rec.append(dict(seed=cs, mismatch0=info["mismatch0"], multi=info["multi"],
                        multi_products=info["multi_products"], mutable_rank=info.get("mutable_rank"),
                        mismatch=info["mismatch"], solved=bool(info["solved"])))
        if info["solved"]:
            text, st = c4.emit_c4(an, new, V)
            out[tag + "_c4"] = dict(seed=cs, toffoli=st["toffoli"], products=st["products"], qubits=st["qubits"])
            if emit:
                with open(f"{emit}_{tag}_c4.kmx", "w") as fh:
                    fh.write(text)
            out[tag + "_closure"] = rec
            print(f"{tag}: CLOSED at closure seed {cs}: emit_c4 {st['toffoli']} Toffolis ({st['products']} "
                  f"products, qubits {st['qubits']})", flush=True)
            return st["toffoli"]
    out[tag + "_closure"] = rec
    mr = [r["mutable_rank"] for r in rec]
    print(f"{tag}: closure not solved in {ncs} seeds (mismatch0 {rec[0]['mismatch0']}, multi sectors "
          f"{rec[0]['multi']} / {rec[0]['multi_products']} products, mutable rank {min(mr)}-{max(mr)}, "
          f"residual {min(r['mismatch'] for r in rec)}-{max(r['mismatch'] for r in rec)} wires)", flush=True)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pad")
    ap.add_argument("proj")
    ap.add_argument("plan")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--stages", default="0,1,2|3,4")
    ap.add_argument("--kills", default="1:2:0")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--effort", type=int, default=1)
    ap.add_argument("--zseed", type=int, default=1)
    ap.add_argument("--tries", type=int, default=12)
    ap.add_argument("--sect", type=int, default=0, help="0 = E.plan(); R = chem_ld2_c4.extract_sectored rounds R")
    ap.add_argument("--xseed", type=int, default=0)
    ap.add_argument("--cseeds", type=int, default=8)
    ap.add_argument("--cseed0", type=int, default=0)
    ap.add_argument("--noctl", action="store_true")
    ap.add_argument("--json", default=None)
    ap.add_argument("--emit", default=None)
    a = ap.parse_args()
    t0 = time.time()
    import chem_ld2_engine as e2
    e2.install()
    words, _ = cc.load_words(a.pad)
    lbs = cc.feasible_lbs(words)
    an, lb = cc.analyze_first(words, a.k, lbs, False)
    assert not lb2.dont_care_low(an), "table has don't-cares: give the V-padded table"
    plan0 = json.load(open(a.plan))
    P = [int(x, 16) for x in plan0["P"]]
    five = plan0.get("blocks") or [[3 * i, 3 * i + 1, 3 * i + 2] for i in range(5)]
    V = int(json.load(open(a.proj))["V"][str(a.k)], 16)
    an2 = lb2.transform(an, P, free_top=False, blocks=five)
    X2 = an2.X0.copy()
    nr, nb = X2.shape
    out = dict(k=a.k, stages=a.stages, kills=a.kills, seed=a.seed, zseed=a.zseed, sect=a.sect, xseed=a.xseed,
               rows=int(words.shape[0]), V=hex(V))
    print(f"k {a.k}: rows {words.shape[0]}, mixed tensor {nr} x {nb}, |V| {bin(V).count('1')}", flush=True)
    if not a.noctl:
        pc, _, codc = plan_stage(an2, X2, five, a.seed, a.effort, a.sect, a.xseed)
        assert np.array_equal(L.tensor_of(pc, nr, nb), X2), "control plan is not exact"
        oldc = lb2.map_back(pc, P)
        _, stc = core.emit(an, oldc)
        out["control"] = dict(products=len(pc), codim=codc, toffoli69=stc["toffoli"])
        print(f"control (C2's engine, five blocks): {len(pc)} products, codim {codc}, emit {stc['toffoli']} "
              f"({time.time() - t0:.0f} s)", flush=True)
        close(an, oldc, V, a.cseed0, a.cseeds, "control", out, a.emit)
    stages = [[int(x) for x in st.split(",")] for st in a.stages.split("|")]
    kills = [tuple(int(x) for x in kk.split(":")) for kk in a.kills.split("|")]
    rng = random.Random(a.zseed)
    prods, feats, report = [], [], []
    Xres = X2.copy()
    for si, st in enumerate(stages):
        blks = [five[i] for i in st]
        coords = [t for b in blks for t in b]
        m = np.uint32(sum(1 << t for t in coords))
        if si:
            kout, krank = kills[si - 1][0], kills[si - 1][1]
            sp = None
            for attempt in range(a.tries):
                pick = rng.sample(range(nb), kout + 2 + 2 * len(blks))     # same rng use as stairb
                killed = sorted(pick[:kout])
                if attempt < 3 and kout:
                    cover = sorted(range(nb), key=lambda o: -sum(1 for f in feats if f[3] >> o & 1))
                    killed = sorted(cover[attempt:attempt + kout])
                Z, zs = [], core.Span()
                while len(Z) < krank:
                    z = rng.getrandbits(nr)
                    if z and zs.add(z):
                        Z.append(z)
                sp, neq, rk = L.solve_spill(Xres, feats, coords, killed, Z, [], nr, nb)
                if sp is not None:
                    break
            if sp is None:
                print(f"stage {si}: spill system has no solution", flush=True)
                out["fail"] = "spill"
                break
            L.apply_spill(Xres, prods, feats, sp, nr)
            L.check_spill(Xres, coords, killed, Z, [], nr)
        try:
            p, _, cod = plan_stage(an2, Xres & m, blks, a.seed, a.effort, a.sect, a.xseed)
        except Exception as exc:  # noqa: BLE001
            print(f"stage {si}: the engine failed: {type(exc).__name__}: {exc}", flush=True)
            out["fail"] = f"stage {si}: {exc}"
            break
        assert np.array_equal(L.tensor_of(p, nr, nb), Xres & m), f"stage {si} plan is not exact"
        n0 = len(prods)
        prods += p
        feats += L.feats_of(p, n0)
        Xres &= ~m
        report.append(dict(stage=st, products=len(p), codim=cod))
    if "fail" not in out:
        ok = np.array_equal(L.tensor_of(prods, nr, nb), X2)
        old = lb2.map_back(prods, P)
        text, st = core.emit(an, old)
        out.update(total=len(prods), tensor_ok=bool(ok), report=report, toffoli69=st["toffoli"])
        print(f"staircase {a.stages} kills {a.kills}: TOTAL {len(prods)} = "
              f"{' + '.join(str(r['products']) for r in report)}; identity {ok}; emit {st['toffoli']} "
              f"({time.time() - t0:.0f} s)", flush=True)
        if a.emit:
            with open(f"{a.emit}_stair69.kmx", "w") as fh:
                fh.write(text)
        if ok:
            close(an, old, V, a.cseed0, a.cseeds, "stair", out, a.emit)
    out["seconds"] = round(time.time() - t0, 1)
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
