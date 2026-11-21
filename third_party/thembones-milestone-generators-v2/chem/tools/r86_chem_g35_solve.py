#!/usr/bin/env python3
"""r86_chem_g35_solve.py: seat R86-CHEM-G35. Group 4 of a 24-output wired loader at 35 products by a STRUCTURED relation
space: the column kill (as built) + C (x) psi = three forced relations sharing ONE low vector psi (C = 3 combinations of
the five non-payer outputs), delivered by the pool's shared kernel (K directions an output) and the outputs' usable dead
cells (pool-compensated on the killed column, C4's rule), then rank 70 by five one-output left-kernel vectors paid by
the payer (C4's C4_DEAD_ROW route). Price: builder C2's exact matroid-parity pricer (chem_ld2_e19grp.Group, read-only)
on ALL 14 relations. Plan: C2's extract_plan; high vectors solved row by row; multiplied out against the dumped
pre-kill residual on every row x coordinate x output (killed column included).

  python -I r86_chem_g35_solve.py k9_g4.npz --payer 18 --rel 19+20,21+22,23+19 --psi auto --draws 20 --seed 1
"""
import argparse
import json
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_ld2_e19grp as eg  # noqa: E402  (Span, nullspace, Group, extract_plan; read-only)

bits = eg.bits


def par(x):
    return bin(x).count("1") & 1


def solve_lin(eqs, n, rng):
    """eqs: list of (mask over n unknowns, rhs). Returns a random solution (int) or None; also the rank."""
    piv = {}
    for mask, rhs in eqs:
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
                return None, len(piv)
    x = 0
    free = [j for j in range(n) if j not in piv]
    for j in free:
        if rng.random() < 0.5:
            x |= 1 << j
    for h in sorted(piv):
        v = piv[h]
        rhs = v >> n & 1
        rest = v & ((1 << n) - 1) & ~(1 << h)
        if rhs ^ par(rest & x):
            x |= 1 << h
        else:
            x &= ~(1 << h)
    for mask, rhs in eqs:
        assert par(mask & x) == rhs
    return x, len(piv)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--payer", type=int, required=True)
    ap.add_argument("--rel", required=True, help="three combinations of outputs, e.g. 19+20,21+22,23+19")
    ap.add_argument("--psi", default="auto", help="hex mask over the live coordinates, or auto / rand")
    ap.add_argument("--draws", type=int, default=20)
    ap.add_argument("--lamtries", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--rank", type=int, default=70)
    ap.add_argument("--out", default="")
    ap.add_argument("--nplans", type=int, default=400)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    z = np.load(a.dump)
    meta = json.load(open(a.dump.replace(".npz", ".json")))
    outs, live, nr, K = meta["outs"], meta["live"], meta["nr"], meta["K"]
    t = meta["killed"][0]
    lc = bits(live)
    A = len(lc)
    b = len(outs)
    gi = {o: i for i, o in enumerate(outs)}
    R0 = [[int(v) for v in row] for row in z["R0"]]
    N = [[int(v) for v in row] for row in z["N"]]
    X = [[int(v) for v in row] for row in z["X"]]          # pre-kill residual (b x nr), all 15 coordinates
    P = [[int(v) for v in row] for row in z["P"]]
    sols = [int(s) for s in meta["sols"]]
    kern = [int(s) for s in meta["kern"]]
    togg = {o: [] for o in outs}
    for c in meta["comp"]:
        togg[c["o"]].append(c)
    payer = a.payer
    rels = [[int(x) for x in r.split("+")] for r in a.rel.split(",")]
    nonp = [o for o in outs if o != payer]
    assert all(o in nonp for r in rels for o in r), "relations must avoid the payer"
    # unknown layout: k_o (K bits) for each non-payer output, then each non-payer toggle
    off, n = {}, 0
    for o in nonp:
        off[o] = n
        n += K
    toff = {}
    for o in nonp:
        for i, c in enumerate(togg[o]):
            toff[(o, i)] = n
            n += 1

    def psi_auto():
        # psi with v.psi = 1 for as many usable toggles as possible: random search over masks on the live coordinates
        best, bs = None, -1
        for _ in range(4000):
            m = rng.getrandbits(15) & live
            if not m:
                continue
            sc = 0
            for o in nonp:
                for c in togg[o]:
                    if any(par(int(v) & m) for v in c["Ec"]):
                        sc += 1
            if sc > bs:
                best, bs = m, sc
        return best
    summary = []
    for draw in range(a.draws):
        if a.psi == "auto":
            psi = psi_auto() if draw == 0 or rng.random() < 0.5 else (rng.getrandbits(15) & live or live)
        elif a.psi == "rand":
            psi = rng.getrandbits(15) & live or live
        else:
            psi = int(a.psi, 16) & live
        # (1) C (x) psi: for each relation, each row: sum_{o in c} parity((R0_o + sum k N + sum toggles)[r] & psi) = 0
        eqs = []
        for r_ in rels:
            for r in range(nr):
                mask, rhs = 0, 0
                for o in r_:
                    rhs ^= par(R0[gi[o]][r] & psi)
                    for j in range(K):
                        if par(N[j][r] & psi):
                            mask ^= 1 << (off[o] + j)
                    for i, c in enumerate(togg[o]):
                        if par(c["Ec"][r] & psi):
                            mask ^= 1 << toff[(o, i)]
                eqs.append((mask, rhs))
        x, rk1 = solve_lin(eqs, n, rng)
        line = dict(draw=draw, psi=hex(psi), eq1=len(eqs), unk1=n, rank1=rk1, ok1=x is not None)
        if x is None:
            summary.append(line)
            print(json.dumps(line), flush=True)
            continue
        # residuals of the non-payer outputs (live columns), and their pool combinations
        res = {}
        comb = {}
        used = {}
        for o in outs:
            res[o] = list(R0[gi[o]])
            comb[o] = sols[gi[o]]
            used[o] = []
        for o in nonp:
            for j in range(K):
                if x >> (off[o] + j) & 1:
                    res[o] = [u ^ v for u, v in zip(res[o], N[j])]
                    comb[o] ^= kern[j]
            for i, c in enumerate(togg[o]):
                if x >> toff[(o, i)] & 1:
                    res[o] = [u ^ v for u, v in zip(res[o], c["Ec"])]
                    comb[o] ^= int(c["tag"])
                    used[o].append(i)
        for r_ in rels:                                   # exact re-check of the forced relations
            for r in range(nr):
                acc_ = 0
                for o in r_:
                    acc_ ^= res[o][r]
                assert par(acc_ & psi) == 0
        # (2) the payer: lambdas in the left kernel of the other five outputs' residual rows
        def rowvec(os_, r):
            v = 0
            for q, o in enumerate(os_):
                v |= res[o][r] << (16 * q)
            return v
        import chem_c4_drop as dr
        _, lk = dr._elim_with_comb([rowvec(nonp, r) for r in range(nr)])
        got = None
        pt = togg[payer]
        npay = K + len(pt)
        def lam_eqs(lam):
            rows_l = bits(lam)
            base_v = 0
            for r in rows_l:
                base_v ^= R0[gi[payer]][r]
            nv = []
            for j in range(K):
                v = 0
                for r in rows_l:
                    v ^= N[j][r]
                nv.append(v)
            uv = []
            for c in pt:
                v = 0
                for r in rows_l:
                    v ^= c["Ec"][r]
                uv.append(v)
            out_ = []
            for tt in lc:
                mask = 0
                for j in range(K):
                    if nv[j] >> tt & 1:
                        mask |= 1 << j
                for i in range(len(pt)):
                    if uv[i] >> tt & 1:
                        mask |= 1 << (K + i)
                out_.append((mask, base_v >> tt & 1))
            return out_

        def apply_payer(y):
            resp = list(R0[gi[payer]])
            cp = sols[gi[payer]]
            up = []
            for j in range(K):
                if y >> j & 1:
                    resp = [u ^ v for u, v in zip(resp, N[j])]
                    cp ^= kern[j]
            for i, c in enumerate(pt):
                if y >> (K + i) & 1:
                    resp = [u ^ v for u, v in zip(resp, c["Ec"])]
                    cp ^= int(c["tag"])
                    up.append(i)
            return resp, cp, up
        for _try in range(a.lamtries):
            eqs2, nl_ = [], 0
            for _att in range(80):
                lam = 0
                for q in lk:
                    if rng.random() < 0.5:
                        lam ^= q
                if not lam:
                    continue
                y, _ = solve_lin(eqs2 + lam_eqs(lam), npay, rng)
                if y is None:
                    continue
                eqs2 += lam_eqs(lam)
                nl_ += 1
                res[payer], comb[payer], used[payer] = apply_payer(y)
                rk = dr._rank([rowvec(outs, r) for r in range(nr)])
                if rk <= a.rank:
                    got = (rk, len(lk), nl_)
                    break
            if got:
                break
        line.update(lk=len(lk), rank=(got[0] if got else dr._rank([rowvec(outs, r) for r in range(nr)])), ok2=got is not None)
        if got is None:
            summary.append(line)
            print(json.dumps(line), flush=True)
            continue
        # (3) price with C2's exact pricer on ALL relations of the final W
        crow = []
        for r in range(nr):
            v = 0
            for q, o in enumerate(outs):
                xx = res[o][r]
                for ci, tt in enumerate(lc):
                    if xx >> tt & 1:
                        v |= 1 << (q * A + ci)
            crow.append(v)
        W = eg.Span()
        wb = [v for v in crow if W.add(v)]
        phis = eg.nullspace(wb, A * b)
        # structured part: (sum over c of e_o) (x) psi in compact coordinates
        psic = sum(1 << ci for ci, tt in enumerate(lc) if psi >> tt & 1)
        Ws = eg.Span()
        for v in phis:
            Ws.add(v)
        struct_in = all(Ws.reduce(sum(psic << (gi[o] * A) for o in r_)) == 0 for r_ in rels)
        g = dict(a=A, b=b, wb=[hex(v) for v in wb], prods=[])
        G = eg.Group(g, seed=a.seed + draw)
        pr = G.price(phis, reps=2)
        line.update(dimW=len(wb), c=len(phis), struct_in_E=struct_in, price=pr["prod"], nu=pr["nu"], need=pr["need"],
                    miss=pr["miss"], prof={str(k_): v for k_, v in sorted(pr["prof"].items())},
                    toggles={str(o): len(used[o]) for o in outs})
        print(json.dumps(line), flush=True)
        summary.append(line)
        if pr["prod"] is not None and pr["prod"] <= 35:
            plan, info = eg.extract_plan(G, phis, rng)
            line["plan_info"] = info
            # (4) high vectors row by row, products, multiplied out against the PRE-KILL residual X
            fs = eg.Span()
            feats = []
            for (y1, y2, s) in plan:
                for y in (y1, y2):
                    feats.append((y, s, G.feat(y, s)))
            piv = {}
            for idx, (_, _, f) in enumerate(feats):
                v, c = f, 1 << idx
                while v:
                    h = v.bit_length() - 1
                    if h not in piv:
                        piv[h] = (v, c)
                        break
                    v ^= piv[h][0]
                    c ^= piv[h][1]
                assert v, "features dependent"
            Acoef = [0] * len(feats)
            for r, v in enumerate(crow):
                c = 0
                while v:
                    h = v.bit_length() - 1
                    assert h in piv, "row outside the plan's span"
                    v ^= piv[h][0]
                    c ^= piv[h][1]
                for idx in bits(c):
                    Acoef[idx] |= 1 << r
            prods = []
            for i in range(0, len(feats), 2):
                (y1c, s, _), (y2c, _, _) = feats[i], feats[i + 1]
                y1 = sum(1 << lc[ci] for ci in bits(y1c))
                y2 = sum(1 << lc[ci] for ci in bits(y2c))
                prods.append((y1, Acoef[i], y2, Acoef[i + 1], sum(1 << outs[q] for q in bits(s))))
            # multiply out: X_o + pool(comb_o) + dead effects E + products = 0 on every row and all 15 coordinates
            bad = 0
            for o in outs:
                acc = list(X[gi[o]])
                for pi in bits(comb[o]):
                    acc = [u ^ v for u, v in zip(acc, P[pi])]
                for i in used[o]:
                    c = togg[o][i]
                    acc = [u ^ v for u, v in zip(acc, c["E"])]
                for (y1, A1, y2, A2, s) in prods:
                    if s >> o & 1:
                        for r in range(nr):
                            m = (y1 if A1 >> r & 1 else 0) ^ (y2 if A2 >> r & 1 else 0)
                            acc[r] ^= m
                bad += sum(1 for v in acc if v)
            line.update(products=len(prods), mismatched_rows=bad)
            print(json.dumps(dict(draw=draw, products=len(prods), mismatched_rows=bad, info=info)), flush=True)
            if a.out and bad == 0:
                # many distinct 35-plans of the same W (greedy random over the fan-out kernels, each re-checked):
                # the engine's C-4 closure retry (tau == kappa, about 2^-6 a plan) cycles through them
                fans = [(s, Y) for s, Y in G.kernels(phis) if len(Y) >= 2]
                plans, keys = [], set()
                for _g in range(a.nplans * 40):
                    if len(plans) >= a.nplans:
                        break
                    order = list(range(len(fans)))
                    rng.shuffle(order)
                    sp, pl = eg.Span(), []
                    for fi in order:
                        s, Y = fans[fi]
                        cand = []
                        d = len(Y)
                        pts = [0]
                        for q in range(d):
                            pts = pts + [p_ ^ Y[q] for p_ in pts]
                        pts = [p_ for p_ in pts if p_]
                        rng.shuffle(pts)
                        y1 = pts[0]
                        y2 = next(p_ for p_ in pts if p_ != y1)
                        f1, f2 = G.feat(y1, s), G.feat(y2, s)
                        save = dict(sp.piv)
                        if sp.add(f1) and sp.add(f2):
                            pl.append((y1, y2, s))
                        else:
                            sp.piv = save
                        if len(pl) == 35:
                            break
                    if len(pl) != 35 or any(sp.reduce(w) for w in wb):
                        continue
                    key = frozenset((frozenset((y1, y2, y1 ^ y2)), s) for (y1, y2, s) in pl)
                    if key in keys:
                        continue
                    keys.add(key)
                    plans.append(dict(planes=[[y1, y2, s] for (y1, y2, s) in pl]))
                json.dump(dict(dump=a.dump, payer=payer, rels=rels, psi=psi, comb={str(o): str(comb[o]) for o in outs},
                               toggle_idx={str(o): used[o] for o in outs},
                               toggles={str(o): [togg[o][i]["cells"] for i in used[o]] for o in outs},
                               prods=[[p[0], hex(p[1]), p[2], hex(p[3]), p[4]] for p in prods], line=line,
                               wb=[hex(v) for v in wb]),
                          open(a.out, "w"))
                json.dump(plans, open(a.out.replace(".json", "_plans.json"), "w"))
                print(f"PLAN WRITTEN {a.out}; {len(plans)} distinct 35-plans -> {a.out.replace('.json', '_plans.json')}",
                      flush=True)
                break
    ok = [s_ for s_ in summary if s_.get("price") is not None]
    print("SUMMARY", json.dumps(dict(draws=len(summary), consistent1=sum(1 for s_ in summary if s_["ok1"]),
                                     rank_met=sum(1 for s_ in summary if s_.get("ok2")),
                                     prices=sorted(s_["price"] for s_ in ok))), flush=True)


if __name__ == "__main__":
    main()
