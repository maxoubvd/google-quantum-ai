#!/usr/bin/env python3
"""R32 re-pricer of INV-CLEAN (2026-11-13): the designer's undo-adder game (tools/r32_inv-clean_pebble.py and
tools/r32_inv-clean_sign.py, copied here, not imported for the sweep) with ONE MOVE ADDED, to test whether the
game's "optimum over every adaptive schedule" is an optimum over the circuits a builder can write.

The added move (flag ccz=1), "repay at distance two":
  a pending phase (-1)^(c_i) is repaid by ONE CCZ and Cliffords, with NO extra wire, while c_(i-2) is on a wire.
  Identity: c_i = MAJ(x_i, y_i, c_(i-1)) = x_i y_i + p_i c_(i-1) (p_i = x_i + y_i over GF(2)) and
            MAJ(x, y, c) = (x + c)(y + c) + c,  so
            (-1)^(c_i) = CZ(x_i, y_i) * CCZ(p_i, x_(i-1) + c_(i-2), y_(i-1) + c_(i-2)) * CZ(p_i, c_(i-2)).
  The three CCZ arguments are put on wires by CNOTs (x_i into y_i; c_(i-2) into x_(i-1) and y_(i-1)), undone after.
  Anchors: c_1 (i = 0): the carry-in itself is one AND of affine forms (feed rows: of v0, low, cur, sign), so
  (-1)^(c_1) = CZ * CCZ(p_1, l1, l2) * CZ(p_1, l3): always available (sign game: needs the sign on its wire);
  c_2 (i = 1): the carry-in on its wire (tail rows: the root wire, always).
  In the designer's game the same repayment is "re-place c_(i-1) (1 AND, one free wire), repay, remove": the
  same Toffoli count, so the move can only matter where a wire is short at the moment of the repayment.

  r32_inv_clean_rp_game.py check                 base game here == designer's game (row prices) on k2 rows
  r32_inv_clean_rp_game.py sched [extra]         k2 clean-up row by row: designer's game / game + CCZ move / rule
  r32_inv_clean_rp_game.py table [nmax]          n x clean wires: game, game + CCZ (feed and tail rows)
  r32_inv_clean_rp_game.py sign N CLEAN GCOST    the sign-pebble game (rows 12, 13), without and with the move
  r32_inv_clean_rp_game.py row N CLEAN FEED      the extended game's schedule (main line + repairs)
"""
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import numba as nb
    HAVE_NUMBA = True
except Exception:
    HAVE_NUMBA = False

INF = 1e18


def _sweep_py(V, L, C, feed, pow3, off, d, ccz):
    delta = 0.0
    for m in range(L, -1, -1):
        for cfg in range(pow3[m]):
            x = cfg
            present = 0
            for i in range(m):
                d[i] = x % 3
                x //= 3
                if d[i] == 1:
                    present += 1
            base = off[m] + cfg
            for f in range(2):
                if not feed and f == 0:
                    continue
                for t in range(2):
                    if t == 1 and m < L:
                        continue
                    sid = base * 4 + f * 2 + t
                    if m == L and t == 1 and cfg == 0 and (f == 0 or not feed):
                        continue
                    wires = present + (f if feed else 0)
                    best = INF
                    if m == L and t == 0 and L > 0 and d[L - 1] == 1:
                        c = 1.0 + V[base * 4 + f * 2 + 1]
                        if c < best:
                            best = c
                    if m == L and t == 0 and L == 0 and f == 1:
                        c = 1.0 + V[base * 4 + f * 2 + 1]
                        if c < best:
                            best = c
                    if m < L and wires < C:
                        pred = (f == 1) if m == 0 else (d[m - 1] == 1)
                        if pred:
                            c = 1.0 + V[(off[m + 1] + cfg + pow3[m]) * 4 + f * 2]
                            if c < best:
                                best = c
                    for i in range(m):
                        pred = (f == 1) if i == 0 else (d[i - 1] == 1)
                        di = d[i]
                        if di == 1:
                            if pred:
                                c = V[(base - pow3[i]) * 4 + f * 2 + t]
                            else:
                                c = 0.5 * V[(base - pow3[i]) * 4 + f * 2 + t] + 0.5 * V[(base + pow3[i]) * 4 + f * 2 + t]
                            if c < best:
                                best = c
                        elif di == 2:
                            if pred:
                                c = V[(base - 2 * pow3[i]) * 4 + f * 2 + t]
                                if c < best:
                                    best = c
                            elif ccz:
                                if i == 0:
                                    anc = True
                                elif i == 1:
                                    anc = f == 1
                                else:
                                    anc = d[i - 2] == 1
                                if anc:
                                    c = 1.0 + V[(base - 2 * pow3[i]) * 4 + f * 2 + t]
                                    if c < best:
                                        best = c
                        else:
                            if pred and wires < C:
                                c = 1.0 + V[(base + pow3[i]) * 4 + f * 2 + t]
                                if c < best:
                                    best = c
                    if feed:
                        if f == 1:
                            c = V[base * 4 + t]
                            if c < best:
                                best = c
                        elif wires < C:
                            c = 1.0 + V[base * 4 + 2 + t]
                            if c < best:
                                best = c
                    if best > V[sid]:
                        if best - V[sid] > delta:
                            delta = best - V[sid]
                        V[sid] = best
    return delta


def _sweep_sign_py(V, L, C, gcost, pow3, off, d, ccz):
    delta = 0.0
    for m in range(L, -1, -1):
        for cfg in range(pow3[m]):
            x = cfg
            present = 0
            for i in range(m):
                d[i] = x % 3
                x //= 3
                if d[i] == 1:
                    present += 1
            base = off[m] + cfg
            for g in range(2):
                for f in range(2):
                    for t in range(2):
                        if t == 1 and m < L:
                            continue
                        sid = ((base * 2 + g) * 2 + f) * 2 + t
                        if m == L and t == 1 and cfg == 0 and f == 0 and g == 1:
                            continue
                        wires = present + f + g
                        best = INF
                        if m == L and t == 0 and d[L - 1] == 1:
                            c = 1.0 + V[((base * 2 + g) * 2 + f) * 2 + 1]
                            if c < best:
                                best = c
                        if m < L and wires < C:
                            pred = (f == 1) if m == 0 else (d[m - 1] == 1)
                            if pred:
                                c = 1.0 + V[(((off[m + 1] + cfg + pow3[m]) * 2 + g) * 2 + f) * 2]
                                if c < best:
                                    best = c
                        for i in range(m):
                            pred = (f == 1) if i == 0 else (d[i - 1] == 1)
                            di = d[i]
                            if di == 1:
                                lo = V[(((base - pow3[i]) * 2 + g) * 2 + f) * 2 + t]
                                if pred:
                                    c = lo
                                else:
                                    hi = V[(((base + pow3[i]) * 2 + g) * 2 + f) * 2 + t]
                                    c = 0.5 * lo + 0.5 * hi
                                if c < best:
                                    best = c
                            elif di == 2:
                                if pred:
                                    c = V[(((base - 2 * pow3[i]) * 2 + g) * 2 + f) * 2 + t]
                                    if c < best:
                                        best = c
                                elif ccz:
                                    if i == 0:
                                        anc = g == 1
                                    elif i == 1:
                                        anc = f == 1
                                    else:
                                        anc = d[i - 2] == 1
                                    if anc:
                                        c = 1.0 + V[(((base - 2 * pow3[i]) * 2 + g) * 2 + f) * 2 + t]
                                        if c < best:
                                            best = c
                            else:
                                if pred and wires < C:
                                    c = 1.0 + V[(((base + pow3[i]) * 2 + g) * 2 + f) * 2 + t]
                                    if c < best:
                                        best = c
                        if g == 1:
                            if f == 1:
                                c = V[((base * 2 + g) * 2 + 0) * 2 + t]
                                if c < best:
                                    best = c
                            elif wires < C:
                                c = 1.0 + V[((base * 2 + g) * 2 + 1) * 2 + t]
                                if c < best:
                                    best = c
                            c = V[((base * 2 + 0) * 2 + f) * 2 + t]
                            if c < best:
                                best = c
                        elif wires < C:
                            c = gcost + V[((base * 2 + 1) * 2 + f) * 2 + t]
                            if c < best:
                                best = c
                        if best > V[sid]:
                            if best - V[sid] > delta:
                                delta = best - V[sid]
                            V[sid] = best
    return delta


if HAVE_NUMBA:
    _sweep = nb.njit(cache=True)(_sweep_py)
    _sweep_sign = nb.njit(cache=True)(_sweep_sign_py)
else:
    _sweep = _sweep_py
    _sweep_sign = _sweep_sign_py


def _tabs(L):
    pow3 = np.array([3 ** i for i in range(L + 2)], dtype=np.int64)
    off = np.zeros(L + 2, dtype=np.int64)
    for m in range(L + 1):
        off[m + 1] = off[m] + pow3[m]
    return pow3, off


def solve(n, C, feed, ccz, tol=1e-10, max_sweeps=20000):
    L = max(0, n - 2)
    pow3, off = _tabs(L)
    V = np.zeros(int(off[L + 1]) * 4, dtype=np.float64)
    d = np.zeros(L + 1, dtype=np.int64)
    for _ in range(max_sweeps):
        delta = _sweep(V, L, C, feed, pow3, off, d, ccz)
        if delta < tol or V[0 if feed else 2] > 1e6:
            break
    start = V[0 if feed else 2]
    return (float(start) if start < 1e6 else float("inf")), V, (L, pow3, off)


def game_price(n, clean, feed, ccz):
    """designer's units: feed rows get clean + 1 wires and their first carry-in AND is priced elsewhere."""
    if n <= 1:
        return 0.0
    if n == 2:
        return 1.0
    val, _, _ = solve(n, clean + (1 if feed else 0), feed, ccz)
    return val - (1.0 if feed else 0.0)


def solve_sign(n, clean, gcost, ccz, tol=1e-10, max_sweeps=20000):
    L = n - 2
    C = clean + 2
    pow3, off = _tabs(L)
    V = np.zeros(int(off[L + 1]) * 8, dtype=np.float64)
    d = np.zeros(L + 1, dtype=np.int64)
    for _ in range(max_sweeps):
        delta = _sweep_sign(V, L, C, float(gcost), pow3, off, d, ccz)
        if delta < tol:
            break
    start = V[((0 * 2 + 1) * 2 + 0) * 2 + 0]
    return float(start) - 1.0, V, (L, pow3, off)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    keep = sys.argv
    sys.argv = [keep[0], "none"]
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.argv = keep
    return mod


K2_ROWS = ([dict(j=j, n=13, feed=True, clean=c) for j, c in ((12, 2), (13, 3), (14, 4), (15, 5))]
           + [dict(j=j, n=13, feed=True, clean=6) for j in range(16, 25)]
           + [dict(j=j, n=13 - t, feed=False, clean=7)
              for j, t in zip(range(25, 36), (0, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10))])


def trace(n, C, feed, V, tabs, pick, limit=400):
    """moves of the extended game's optimal policy (ccz = 1) for a given outcome picker."""
    L, pow3, off = tabs
    m, cfg, f, t = 0, 0, (0 if feed else 1), 0
    out = []
    k = 0
    for _ in range(limit):
        d = [(cfg // int(pow3[i])) % 3 for i in range(m)]
        wires = sum(1 for x in d if x == 1) + (f if feed else 0)
        base = int(off[m]) + cfg
        if m == L and t == 1 and cfg == 0 and (f == 0 or not feed):
            break
        cands = []
        if m == L and t == 0 and L > 0 and d[L - 1] == 1:
            cands.append((1.0 + V[base * 4 + f * 2 + 1], "top", (m, cfg, f, 1)))
        if m < L and wires < C and ((f == 1) if m == 0 else d[m - 1] == 1):
            cands.append((1.0 + V[(int(off[m + 1]) + cfg + int(pow3[m])) * 4 + f * 2], f"place c{m + 1}",
                          (m + 1, cfg + int(pow3[m]), f, 0)))
        for i in range(m):
            pred = (f == 1) if i == 0 else (d[i - 1] == 1)
            p3 = int(pow3[i])
            if d[i] == 1:
                if pred:
                    cands.append((V[(base - p3) * 4 + f * 2 + t], f"remove c{i + 1}", (m, cfg - p3, f, t)))
                else:
                    e = 0.5 * V[(base - p3) * 4 + f * 2 + t] + 0.5 * V[(base + p3) * 4 + f * 2 + t]
                    cands.append((e, f"GHOST c{i + 1}", ((m, cfg - p3, f, t), (m, cfg + p3, f, t))))
            elif d[i] == 2:
                if pred:
                    cands.append((V[(base - 2 * p3) * 4 + f * 2 + t], f"repay c{i + 1} (free)", (m, cfg - 2 * p3, f, t)))
                else:
                    anc = True if i == 0 else ((f == 1) if i == 1 else d[i - 2] == 1)
                    if anc:
                        cands.append((1.0 + V[(base - 2 * p3) * 4 + f * 2 + t],
                                      f"CCZ-repay c{i + 1} from c{i - 1} (1 CCZ, no wire)", (m, cfg - 2 * p3, f, t)))
            else:
                if pred and wires < C:
                    cands.append((1.0 + V[(base + p3) * 4 + f * 2 + t], f"re-place c{i + 1}", (m, cfg + p3, f, t)))
        if feed:
            if f == 1:
                cands.append((V[base * 4 + t], "erase carry-in", (m, cfg, 0, t)))
            elif wires < C:
                cands.append((1.0 + V[base * 4 + 2 + t], "make carry-in", (m, cfg, 1, t)))
        cands.sort(key=lambda z: (round(z[0], 9), 0 if ("remove" in z[1] or "repay" in z[1] and "CCZ" not in z[1]
                                                          or "erase" in z[1]) else 1))
        val, name, nxt = cands[0]
        if name.startswith("GHOST"):
            nxt = nxt[pick(k)]
            k += 1
        out.append((name, wires, val, (m, cfg, f, t)))
        m, cfg, f, t = nxt
    return out


def trace_sign(n, clean, gcost, V, tabs, pick, ccz=1, limit=400):
    """moves of the sign-pebble game's optimal policy (with the CCZ move when ccz = 1)."""
    L, pow3, off = tabs
    C = clean + 2
    m, cfg, g, f, t = 0, 0, 1, 0, 0
    out = []
    k = 0

    def sid(m_, cfg_, g_, f_, t_):
        return (((int(off[m_]) + cfg_) * 2 + g_) * 2 + f_) * 2 + t_

    for _ in range(limit):
        d = [(cfg // int(pow3[i])) % 3 for i in range(m)]
        wires = sum(1 for x in d if x == 1) + f + g
        if m == L and t == 1 and cfg == 0 and f == 0 and g == 1:
            break
        cands = []
        if m == L and t == 0 and d[L - 1] == 1:
            cands.append((1.0 + V[sid(m, cfg, g, f, 1)], "top", (m, cfg, g, f, 1)))
        if m < L and wires < C and ((f == 1) if m == 0 else d[m - 1] == 1):
            cands.append((1.0 + V[sid(m + 1, cfg + int(pow3[m]), g, f, 0)], f"place c{m + 1}",
                          (m + 1, cfg + int(pow3[m]), g, f, 0)))
        for i in range(m):
            pred = (f == 1) if i == 0 else (d[i - 1] == 1)
            p3 = int(pow3[i])
            if d[i] == 1:
                if pred:
                    cands.append((V[sid(m, cfg - p3, g, f, t)], f"remove c{i + 1}", (m, cfg - p3, g, f, t)))
                else:
                    e = 0.5 * V[sid(m, cfg - p3, g, f, t)] + 0.5 * V[sid(m, cfg + p3, g, f, t)]
                    cands.append((e, f"GHOST c{i + 1}", ((m, cfg - p3, g, f, t), (m, cfg + p3, g, f, t))))
            elif d[i] == 2:
                if pred:
                    cands.append((V[sid(m, cfg - 2 * p3, g, f, t)], f"repay c{i + 1} (free)", (m, cfg - 2 * p3, g, f, t)))
                elif ccz:
                    anc = (g == 1) if i == 0 else ((f == 1) if i == 1 else d[i - 2] == 1)
                    if anc:
                        src = "the carry-in's inputs (sign on its wire)" if i == 0 else (
                            "the carry-in" if i == 1 else f"c{i - 1}")
                        cands.append((1.0 + V[sid(m, cfg - 2 * p3, g, f, t)],
                                      f"CCZ-repay c{i + 1} from {src} (1 CCZ, no wire)", (m, cfg - 2 * p3, g, f, t)))
            else:
                if pred and wires < C:
                    cands.append((1.0 + V[sid(m, cfg + p3, g, f, t)], f"re-place c{i + 1}", (m, cfg + p3, g, f, t)))
        if g == 1:
            if f == 1:
                cands.append((V[sid(m, cfg, g, 0, t)], "erase carry-in", (m, cfg, g, 0, t)))
            elif wires < C:
                cands.append((1.0 + V[sid(m, cfg, g, 1, t)], "make carry-in", (m, cfg, g, 1, t)))
            cands.append((V[sid(m, cfg, 0, f, t)], "MEASURE SIGN", (m, cfg, 0, f, t)))
        elif wires < C:
            cands.append((gcost + V[sid(m, cfg, 1, f, t)], f"REBUILD SIGN ({gcost:g} AND)", (m, cfg, 1, f, t)))
        free = ("remove", "repay c", "erase", "MEASURE")
        cands.sort(key=lambda z: (round(z[0], 9), 0 if any(z[1].startswith(w) for w in free) else 1))
        val, name, nxt = cands[0]
        if name.startswith("GHOST"):
            nxt = nxt[pick(k)]
            k += 1
        out.append((name, wires, val, (m, cfg, g, f, t)))
        m, cfg, g, f, t = nxt
    return out


def cost_of(moves, gcost):
    c = 0.0
    for name, _, _, _ in moves:
        if name.startswith(("place", "re-place", "make", "top", "CCZ")):
            c += 1.0
        elif name.startswith("REBUILD"):
            c += gcost
    return c


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "signrow":
        n, clean, gcost = int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
        ccz = int(sys.argv[5]) if len(sys.argv) > 5 else 1
        val, V, tabs = solve_sign(n, clean, gcost, ccz)
        print(f"n {n} clean {clean} gcost {gcost:g} ccz {ccz}: optimum {val:.4f} (row price; first carry-in AND "
              f"and today's sign rebuild priced elsewhere)")
        main = trace_sign(n, clean, gcost, V, tabs, lambda q: 0, ccz)
        gpos = [i for i, z in enumerate(main) if z[0].startswith("GHOST")]
        print(" main line (every ghost returns 0):")
        for i, (name, w, v, st) in enumerate(main):
            print(f"   {i:3d} wires {w}: {name}")
        print(f" ANDs on the main line: {cost_of(main, gcost):g} (incl. the first carry-in AND)")
        G = len(gpos)
        tot = 0.0
        for mask in range(1 << G):
            mv = trace_sign(n, clean, gcost, V, tabs, lambda q, mask=mask: (mask >> q) & 1 if q < G else 0, ccz)
            tot += cost_of(mv, gcost) / (1 << G)
        print(f" expected ANDs over the 2^{G} outcome paths of the main-line ghosts (later ghosts 0): {tot:.4f} "
              f"(row price {tot - 1:.4f})")
        for g_, idx in enumerate(gpos):
            alt = trace_sign(n, clean, gcost, V, tabs, lambda q, g_=g_: 1 if q == g_ else 0, ccz)
            later = {st: i for i, (_, _, _, st) in enumerate(main) if i > idx}
            print(f" repair when ghost #{g_} ({main[idx][0]}) returns 1:")
            for name, w, v, st in alt[idx + 1:]:
                if st in later:
                    print(f"      back on the main line at step {later[st]}")
                    break
                print(f"      wires {w}: {name}")
        sys.exit(0)
    if cmd == "check":
        D = _load("pebble", "r32_inv-clean_pebble.py")
        seen = {}
        bad = 0
        for r in K2_ROWS:
            key = (r["n"], r["clean"], r["feed"])
            if key in seen:
                continue
            a = D.game_price(r["n"], r["clean"], r["feed"])
            b = game_price(r["n"], r["clean"], r["feed"], 0)
            seen[key] = (a, b)
            bad += abs(a - b) > 1e-9
            print(f"n {r['n']:2d} clean {r['clean']} feed {int(r['feed'])}: designer {a:.4f}  copy {b:.4f}", flush=True)
        print("check:", "EQUAL" if bad == 0 else f"{bad} DIFFER")
    elif cmd == "sched":
        P = _load("invadd", "r31_inv-wire_invadd.py")
        extra = int(sys.argv[2]) if len(sys.argv) > 2 else 0
        cache = {}
        tf = tr = tg = tc = 0.0
        print("row  n feed clean  forward  rule    game   game+CCZ  best(rule,game)  best(rule,game+CCZ)")
        for r in K2_ROWS:
            clean = r["clean"] + extra
            key = (r["n"], clean, r["feed"])
            if key not in cache:
                g0 = game_price(r["n"], clean, r["feed"], 0)
                g1 = game_price(r["n"], clean, r["feed"], 1)
                pr = P.inv_cost(r["n"], r["n"] - 2, clean, r["feed"])
                cache[key] = (g0, g1, pr[0] if pr is not None else float("inf"))
            g0, g1, pub = cache[key]
            b0, b1 = min(g0, pub), min(g1, pub)
            tf += r["n"] - 1
            tr += pub
            tg += b0
            tc += b1
            print(f"{r['j']:3d} {r['n']:2d} {int(r['feed']):4d} {clean:5d}  {r['n'] - 1:7d}  {pub:5.1f}  {g0:7.3f}  {g1:7.3f}"
                  f"  {b0:8.3f}  {b1:8.3f}", flush=True)
        print(f"forward {tf:.1f}; rule {tr:.2f} (penalty {tr - tf:.2f}); best(rule, game) {tg:.2f} (penalty {tg - tf:.2f}); "
              f"best(rule, game+CCZ) {tc:.2f} (penalty {tc - tf:.2f})")
    elif cmd == "table":
        nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 13
        for feed in (True, False):
            print(f"feed {int(feed)}: n, clean 0..9: game / game+CCZ")
            for n in range(4, nmax + 1):
                line = []
                for clean in range(0, 10):
                    a = game_price(n, clean, feed, 0)
                    b = game_price(n, clean, feed, 1)
                    line.append(f"{a:6.3f}/{b:6.3f}")
                print(f"  n {n:2d}: " + " ".join(line), flush=True)
    elif cmd == "sign":
        n, clean, gcost = int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
        a, _, _ = solve_sign(n, clean, gcost, 0)
        b, _, _ = solve_sign(n, clean, gcost, 1)
        h0 = game_price(n, clean, True, 0)
        h1 = game_price(n, clean, True, 1)
        print(f"n {n} clean {clean} gcost {gcost}: sign held: game {h0:.4f}, +CCZ {h1:.4f}; "
              f"sign pebble: game {a:.4f}, +CCZ {b:.4f}", flush=True)
    elif cmd == "row":
        n, clean, feed = int(sys.argv[2]), int(sys.argv[3]), bool(int(sys.argv[4]))
        C = clean + (1 if feed else 0)
        val0, _, _ = solve(n, C, feed, 0)
        val, V, tabs = solve(n, C, feed, 1)
        sub = 1.0 if feed else 0.0
        print(f"n {n} clean {clean} feed {int(feed)}: game {val0 - sub:.4f}; game+CCZ {val - sub:.4f}")
        main = trace(n, C, feed, V, tabs, lambda k: 0)
        gpos = [i for i, z in enumerate(main) if z[0].startswith("GHOST")]
        print(" main line (every ghost returns 0):")
        for i, (name, w, v, st) in enumerate(main):
            print(f"   {i:3d} wires {w}: {name}")
        for g, idx in enumerate(gpos):
            alt = trace(n, C, feed, V, tabs, lambda q, g=g: 1 if q == g else 0)
            later = {st: i for i, (_, _, _, st) in enumerate(main) if i > idx}
            print(f" repair when ghost #{g} ({main[idx][0]}) returns 1:")
            for name, w, v, st in alt[idx + 1:]:
                if st in later:
                    print(f"      back on the main line at step {later[st]}")
                    break
                print(f"      wires {w}: {name}")
