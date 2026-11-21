#!/usr/bin/env python3
"""r65_chem_wired221_dc.py: designer R65-CHEM-WIRED221. LIVE KILL UNITS AT THEIR EXACT PRICE.

R63-CHEM-WIRED's owner engine (tools/r63_chem_wired_dc.py, imported read-only: DcKill, make_tie, rebase, ...) pays
the kill systems of a wired loader with re-fan-out bits and with the free bits of dead words. Here a third payer
joins: bit b (b <= --maxbit) of a LIVE two-body word, priced by the exact single-flip change of the verifier's
error_2^2 (tools/r65_chem_wired221_exact.py: formula checked against the verifier to 1e-6 mHa^2). A unit joins when
its price is at most --tau (mHa^2); units are offered cheapest first after the pool and the dead bits, so the solve
takes a cheapest basis. The changed table is evaluated by the verifier's own evaluate_words (weights held).

    CHEM_BROT=7 CHEM_HB=77 python -I tools/r65_chem_wired221_dc.py TABLE.npz --k 14 --plan "6,6:1:3:0,6:1:3:1,7:2:0:2"
        [--modes reg,tie3,-,-] [--avoid 0,1] [--seeds 0-7] [--pad] [--own cur,prev] [--worder auto]
        [--tau 1e-3] [--maxbit 1] [--emit DIR] [--stop N]
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("CHEM_BROT", "7")
os.environ.setdefault("CHEM_HB", "77")
import chem_loaders_core as core  # noqa: E402
import chem_ld2_wwidth  # noqa: E402,F401
import chem_ld2_e17kill as ek  # noqa: E402
import r63_chem_wired_dc as r63  # noqa: E402
from r62_chem_cache_dead import dead_mask  # noqa: E402

bits = core.bits


def free_bits(words, edges, k, an, V, own, lmin, live=None):
    """dict basis output -> list of free bits. live = (tau, maxbit, price[bit][col] arrays over all rows, n1)."""
    A0, ncol = words.shape
    J = core.J
    nr = an.nr
    dead = dead_mask(words, edges, core.BROT)
    qidx = {w: q for q, w in enumerate(an.out_wires)}
    ginv = r63.ginv_masks(an)
    L = an.L
    dl = []
    for l in range(L):
        v = 0
        for m in range(1, L):
            if (m & l) == l:
                v ^= an.mono_coords[m]
        dl.append(v)
    last = an.H - 1

    def unit_vec(h, l):
        vec = np.zeros(nr, np.uint32)
        for M in range(1, nr + 1):
            if (M & h) == h:
                vec[M - 1] = dl[l]
        if (V is not None) and ((V >> h) & 1) and last <= nr:
            vec[last - 1] ^= np.uint32(dl[l])
        return vec

    free = {}
    seen = set()
    cols = []
    if "cur" in own:
        cols += [(J * k + j, j) for j in range(J)]
    if "prev" in own and k > 0:
        cols += [(J * (k - 1) + j, j) for j in range(J)]
    nfree = 0
    for col, j in cols:
        for a in np.nonzero(dead[:, col])[0].tolist():
            h, l = a >> 4, a & 15
            if l < lmin or h >= last or (a, j) in seen:
                continue
            seen.add((a, j))
            vec = unit_vec(h, l)
            if not vec.any():
                continue
            for b in range(1, core.BROT - 1):
                q = qidx.get((j, b))
                if q is None or q not in ginv:
                    continue
                o, gm = ginv[q]
                free.setdefault(o, []).append(dict(vec=vec, g=gm, a=a, col=col, bit=b, h=h, l=l, cost=0.0))
                nfree += 1
    nlive = 0
    if live is not None:
        tau, maxbit, price, n1 = live
        for j in range(J):
            cj = [c for c, jj in cols if jj == j]
            for a in range(n1, A0):
                h, l = a >> 4, a & 15
                if l < lmin or h >= last or (a, j) in seen:
                    continue
                vec = None
                for b in range(1, maxbit + 1):
                    q = qidx.get((j, b))
                    if q is None or q not in ginv:
                        continue
                    best = None
                    for c in cj:
                        p = float(price[b][c][a - n1])
                        if best is None or p < best[0]:
                            best = (p, c)
                    if best is None or best[0] > tau:
                        continue
                    if vec is None:
                        vec = unit_vec(h, l)
                        if not vec.any():
                            break
                    o, gm = ginv[q]
                    free.setdefault(o, []).append(dict(vec=vec, g=gm, a=a, col=best[1], bit=b, h=h, l=l, cost=best[0]))
                    nlive += 1
        for o in free:
            free[o].sort(key=lambda f: f["cost"])
    return free, nfree, nlive, int(dead[:, [c for c, _ in cols]].sum())


def load_prices(path, tau, maxbit):
    """live tuple (tau, maxbit, price[bit][col], n1) from a --pricesonly file."""
    pz = np.load(path)
    pb, pc, pd = [int(x) for x in pz["bits"]], [int(x) for x in pz["cols"]], pz["dE2"]
    price = {b: {c: pd[bi][ci] for ci, c in enumerate(pc)} for bi, b in enumerate(pb)}
    assert all(b in price for b in range(1, maxbit + 1)), "price file lacks a bit"
    return (tau, maxbit, price, int(pz["n1"]))


def supply(an, words, edges, k, own, lmin, live):
    """SECOND HOLDER: table units an output can pay kills with: dead addresses of its word + live words of its bit
    whose exact flip price is at most tau. dict (word, bit) -> (dead, live)."""
    J = core.J
    dpw = r63.dead_per_word(words, edges, k, own, lmin)
    cols = {j: ([J * k + j] if "cur" in own else []) + ([J * (k - 1) + j] if ("prev" in own and k > 0) else [])
            for j in range(J)}
    A0 = words.shape[0]
    last = an.H - 1
    dead = dead_mask(words, edges, core.BROT)
    out = {}
    for q, (j, b) in enumerate(an.out_wires):
        if not (1 <= b <= core.BROT - 2):
            continue
        nl = 0
        if live is not None and b <= live[1]:
            tau, maxbit, price, n1 = live
            a = np.arange(n1, A0)
            ok = ((a & 15) >= lmin) & ((a >> 4) < last) & ~dead[n1:, cols[j]].any(axis=1)
            best = np.min(np.stack([np.asarray(price[b][c], float) for c in cols[j]]), axis=0)
            nl = int((ok & (best <= tau)).sum())
        out[(j, b)] = (dpw[j], nl)
    return out


def placed_order(an, sup, sizes, pattern="", explicit=""):
    """wire order for r63.rebase: OUTPUTS (not words) placed into the plan's groups.
    explicit = "w:b,w:b,.." (the generic outputs in group order), or pattern = the supply class of every group
    (0 = poorest), e.g. 0,3,2,1: group 1 takes the poorest outputs, group 2 the richest, group 3 the next, group 4
    the second poorest. Wired (bit 0) wires first, bit-6 copies last, as analyze() needs no order among them."""
    wires = {w: q for q, w in enumerate(an.out_wires)}
    gens = [w for w in an.out_wires if 1 <= w[1] <= core.BROT - 2]
    if explicit:
        seq = [tuple(int(x) for x in t.split(":")) for t in explicit.split(",")]
        assert sorted(seq) == sorted(gens), "explicit order must list every generic output once"
    else:
        asc = sorted(gens, key=lambda w: (sup[w][0] + sup[w][1], sup[w][0], -w[1], w[0]))
        cls = [int(x) for x in pattern.split(",")]
        assert sorted(cls) == list(range(len(sizes))), "pattern must be a permutation of the groups"
        take, pos = {}, 0
        for c in range(len(sizes)):
            g = cls.index(c)
            take[g] = sorted(asc[pos:pos + sizes[g]])
            pos += sizes[g]
        seq = [w for g in range(len(sizes)) for w in take[g]]
    first = [q for q, w in enumerate(an.out_wires) if w[1] == 0]
    lastw = [q for q, w in enumerate(an.out_wires) if w[1] > core.BROT - 2]
    return first + [wires[w] for w in seq] + lastw, seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--plan", required=True)
    ap.add_argument("--avoid", default="")
    ap.add_argument("--seeds", default="0-3")
    ap.add_argument("--effort", type=int, default=3)
    ap.add_argument("--pad", action="store_true")
    ap.add_argument("--tries", type=int, default=400)
    ap.add_argument("--own", default="cur,prev")
    ap.add_argument("--lmin", type=int, default=3)
    ap.add_argument("--emit", default=None)
    ap.add_argument("--stop", type=int, default=0)
    ap.add_argument("--modes", default="")
    ap.add_argument("--tau", type=float, default=0.0, help="live units: exact price of the flip at most this (mHa^2)")
    ap.add_argument("--maxbit", type=int, default=1)
    ap.add_argument("--worder", default="")
    ap.add_argument("--wpat", default="", help="word order as a pattern over the words ranked by dead addresses (0 = poorest), e.g. 0,3,1,2,4")
    ap.add_argument("--border", action="store_true", help="outputs ordered by bit plane first (bit 1 of every word, then bit 2 ...)")
    ap.add_argument("--noeval", action="store_true")
    ap.add_argument("--oauto", default="", help="SECOND HOLDER: outputs placed by supply: class of every group (0 = poorest), e.g. 0,3,2,1")
    ap.add_argument("--oorder", default="", help="SECOND HOLDER: the generic outputs in group order, w:b,w:b,..")
    ap.add_argument("--prices", default=None, help="npz of r65_chem_wired221_exact.py --pricesonly (bits, cols, dE2 in mHa^2): no pricer in this process")
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    z = np.load(a.table)
    d = {kk: np.array(z[kk]) for kk in z.files}
    words = np.asarray(d["angle_words"], np.int64)
    edges = np.asarray(d["tree_edges"], np.int64)
    V, kap = None, 0
    wp = words
    if a.pad:
        import chem_ld2_c4 as c4
        V = c4.vset(words, a.k)
        wp = c4.pad_table(words, {a.k: V})
    an = core.analyze(wp, a.k, 4)
    if a.pad:
        kap = c4.kappa_wires(an, V)
    own = a.own.split(",")
    if a.worder or a.border or a.wpat:
        dpw = r63.dead_per_word(words, edges, a.k, own, a.lmin)
        if a.wpat:
            rk = sorted(range(core.J), key=lambda j: (dpw[j], j))
            wo = [rk[int(x)] for x in a.wpat.split(",")]
        elif a.worder in ("", "auto"):
            wo = sorted(range(core.J), key=lambda j: (dpw[j], j))
        else:
            wo = [int(x) for x in a.worder.split(",")]
        if a.border:
            order = sorted(range(len(an.out_wires)), key=lambda q: (min(an.out_wires[q][1], 2), wo.index(an.out_wires[q][0]), an.out_wires[q][1]))
        else:
            order = sorted(range(len(an.out_wires)), key=lambda q: (wo.index(an.out_wires[q][0]), an.out_wires[q][1]))
        r63.rebase(an, order)
        print(f"k{a.k}: dead addresses a word {dpw}; word order {wo}; bit-plane first {a.border}", flush=True)
    live = None
    P = None
    if a.prices:
        a.noeval = True
        if a.tau > 0:
            live = load_prices(a.prices, a.tau, a.maxbit)
    if a.oauto or a.oorder:
        sup = supply(an, words, edges, a.k, own, a.lmin, live)
        sizes = [g[0] for g in ek.parse_plan(a.plan)]
        order, seq = placed_order(an, sup, sizes, a.oauto, a.oorder)
        r63.rebase(an, order)
        pos, txt = 0, []
        for b in sizes:
            txt.append(" ".join(f"{w[0]}:{w[1]}({sup[w][0]}+{sup[w][1]})" for w in seq[pos:pos + b]))
            pos += b
        print(f"k{a.k}: outputs placed (word:bit(dead+live units)): " + " | ".join(txt), flush=True)
    elif a.tau > 0 or not a.noeval:
        import r65_chem_wired221_exact as ex
        P = ex.Pricer(d)
    if a.tau > 0 and not a.prices:
        J = core.J
        cols = [J * a.k + j for j in range(J)] + ([J * (a.k - 1) + j for j in range(J)] if a.k > 0 else [])
        price = {}
        for b in range(1, a.maxbit + 1):
            price[b] = {}
            for c in cols:
                dE2, _, _ = P.flip_prices(c, xormask=(1 << b) | ((1 << 6) if b == 5 else 0))
                price[b][c] = dE2 * 1e6
        live = (a.tau, a.maxbit, price, P.n1)
    free, nfree, nlive, ndead = free_bits(words, edges, a.k, an, V, own, a.lmin, live)
    print(f"k{a.k}: generic outputs {an.X0.shape[1]}, infra {an.infra}; dead words owned {ndead}; dead bits {nfree}; "
          f"live units {nlive} (tau {a.tau}, bits 1..{a.maxbit})", flush=True)
    plan = ek.parse_plan(a.plan)
    avoid = [int(x) for x in a.avoid.split(",") if x.strip()]
    best = None
    for s in ek.rng_list(a.seeds):
        t0 = time.time()
        if a.modes:
            kd = r63.make_tie(an, s, max(2, a.effort), plan, avoid, a.modes.split(","), V, kap, a.tries if a.pad else 1, free)
        else:
            kd = r63.DcKill(an, s, max(2, a.effort), plan, avoid, V, kap, a.tries if a.pad else 1, free)
        try:
            prods = kd.run()
        except (RuntimeError, AssertionError) as exc:
            print(f"seed {s}: FAILED {exc}", flush=True)
            continue
        n = len(prods)
        used = {o: list(v) for o, v in kd.used.items()}
        nflip = sum(len(v) for v in used.values())
        lcost = sum(f["cost"] for lst in used.values() for f in lst)
        nl = sum(1 for lst in used.values() for f in lst if f["cost"] != 0.0)
        line = " | ".join(f"{b}->{m}" + (f" [kills {st[0]},{st[1]} pool {st[2]} free {st[3]}]" if st else "")
                          for b, m, st in kd.gcount)
        print(f"seed {s}: {n} products = {an.infra + n}{' (-1 with C-4)' if a.pad else ''}; {line}; bits flipped "
              f"{nflip}, live {nl}, formula price {lcost:+.5f} mHa^2 ({time.time() - t0:.1f} s)", flush=True)
        key = (n, lcost)
        if best is None or key < best[0]:
            best = (key, s, prods, used)
        if a.stop and n <= a.stop:
            break
    if best is None:
        return
    (n, lcost), s, prods, used = best
    print(f"BEST: {n} products, seed {s}, formula price of the live flips {lcost:+.5f} mHa^2", flush=True)
    X0 = an.X0.copy()
    for o, lst in used.items():
        for f in lst:
            for o2 in bits(f["g"]):
                X0[:, o2] ^= f["vec"]
    w2 = r63.modified_words(words, used)
    nch = int((w2 != words).sum())
    dead = dead_mask(words, edges, core.BROT)
    live_changed = int(((w2 != words) & ~dead).sum())
    print(f"words changed {nch}, of them live {live_changed}", flush=True)
    if P is not None and not a.noeval:
        e2b, e1b, lb = P.exact(words)
        e2, e1, lam = P.exact(w2)
        print(f"VERIFIER'S EVALUATION (weights held): d(error_2^2) {(e2 - e2b) * 1e6:+.5f} mHa^2 (base "
              f"{e2b * 1e6:.4f}); error_1 {e1 * 1e3:.5f} mHa (base {e1b * 1e3:.5f}); lambda_eff {lam:.6f} "
              f"(base {lb:.6f})", flush=True)
    w2p = w2
    if a.pad:
        V2 = c4.vset(w2, a.k)
        if V2 != V:
            print("the block set V moved: no emission", flush=True)
            return
        w2p = c4.pad_table(w2, {a.k: V})
    an2 = core.analyze(w2p, a.k, 4)
    t_new = r63.wire_tensor(an2, an2.X0)
    an2.G, an2.X0, an2.out_level = an.G, X0, an.out_level
    t_mod = r63.wire_tensor(an2, X0)
    bad = sum(1 for x, y in zip(t_new, t_mod) if not np.array_equal(x, y))
    print(f"wires whose tensor differs from the model: {bad} of {len(t_new)}", flush=True)
    if bad:
        return
    if a.pad and c4.kappa_wires(an2, V) != kap:
        print("kappa moved: no emission", flush=True)
        return
    if a.emit:
        os.makedirs(a.emit, exist_ok=True)
        if a.pad:
            text, st = c4.emit_c4(an2, prods, V)
        else:
            text, st = core.emit(an2, prods)
        with open(os.path.join(a.emit, f"angle_k_{a.k}.kmx"), "w") as fh:
            fh.write(text)
        np.savez(os.path.join(a.emit, f"table_k{a.k}.npz"), angle_words=w2, tree_edges=edges)
        print(f"emitted {a.emit}/angle_k_{a.k}.kmx: toffoli {st['toffoli']}, qubits {st.get('qubits')} (checked on "
              f"every address by the emitter); new words in {a.emit}/table_k{a.k}.npz", flush=True)


if __name__ == "__main__":
    main()
