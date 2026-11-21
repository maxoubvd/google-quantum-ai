#!/usr/bin/env python3
"""R32 seat INV-CLEAN (2026-11-13): the width-limited undo adder as an exact game, solved to optimality. Own code.

One undo row of n target bits (n - 2 real source bits, the public rows) holds L = n - 2 carries c_1..c_L on
wires; the last carry goes into the top bit with one more AND ("top"). Rules, each read in the generator's
inverse adder (paired_arithmetic.inverse_add, fusion.arithmetic) and in kickmix's HMR:
  place c_i     1 AND, needs c_(i-1) on a wire (c_0 = the row's carry-in) and a free wire; the first placement
                of c_i also writes target bit i and pays the forward row's debt on c_i (a Z); a re-placement
                pays a pending phase of c_i (a Z)
  remove c_i    free while c_(i-1) is on a wire (X-measurement, CZ repair on the bit below)
  ghost c_i     free at any time (X-measurement); with probability 1/2 it leaves a phase (-1)^(c_i) to repay
  repay c_i     free while c_(i-1) is on a wire (majority is quadratic), or by re-placing c_i
  top           1 AND once all L carries have been placed and c_L is on a wire
  carry-in      tail rows: a root wire, always there, no wire of ours. Feed rows: one AND to make, free to
                erase at any time (a quadratic function of v0, low, cur, sign), counted in the wire budget.
A row is done when top is done, no carry is on a wire and no phase is pending. The solver minimises the EXPECTED
number of ANDs over every adaptive schedule with at most C wires at any time (value iteration on all states;
states: per placed carry {absent, on a wire, absent with a pending phase}, the placed carries are a prefix).
The generator's rule (segments whose boundary carries stay to the end, then measured and repaired by a comparator
with probability 1/2) is one family of schedules of this game; its prices come from tools/r31_inv-wire_invadd.py
(the first-round model, exact against the public file).

  r32_inv-clean_pebble.py table [nmax]          optimum against the generator's rule, n and wires
  r32_inv-clean_pebble.py row N C FEED          one row: optimum, and the schedule on the all-zero outcome path
  r32_inv-clean_pebble.py sched ROWS.json|k2    the divider clean-up row by row (ROWS from r31_inv-wire_rows.py)
"""
import importlib.util
import json
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


def _sweep_py(V, L, C, feed, pow3, off, d):
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
                            elif wires < C and False:
                                pass
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


if HAVE_NUMBA:
    _sweep = nb.njit(cache=True)(_sweep_py)
else:
    _sweep = _sweep_py


def solve(n, C, feed, tol=1e-10, max_sweeps=5000):
    """optimal expected ANDs of one undo row: n target bits, C wires for carries (feed rows: the carry-in wire
    counts and its first making is included in the result). Returns (value, V, tables) or (inf, ...)."""
    L = max(0, n - 2)
    pow3 = np.array([3 ** i for i in range(L + 2)], dtype=np.int64)
    off = np.zeros(L + 2, dtype=np.int64)
    for m in range(L + 1):
        off[m + 1] = off[m] + pow3[m]
    V = np.zeros(int(off[L + 1]) * 4, dtype=np.float64)
    d = np.zeros(L + 1, dtype=np.int64)
    prev = -1.0
    for it in range(max_sweeps):
        delta = _sweep(V, L, C, feed, pow3, off, d)
        start = V[(0) * 4 + (0 if feed else 2)]
        if delta < tol:
            break
        if start > 1e6:
            break
        prev = start
    start = V[0 * 4 + (0 if feed else 2)]
    return (float(start) if start < 1e6 else float("inf")), V, (L, pow3, off)


def trace(n, C, feed, V, tabs, outcome=0, limit=400):
    """moves of the optimal policy when every ghost measurement returns `outcome`."""
    L, pow3, off = tabs
    m, cfg, f, t = 0, 0, (0 if feed else 1), 0
    out = []
    ghosts = 0
    pick = outcome if callable(outcome) else (lambda k: outcome)
    for _ in range(limit):
        d = [(cfg // int(pow3[i])) % 3 for i in range(m)]
        present = sum(1 for x in d if x == 1)
        wires = present + (f if feed else 0)
        base = int(off[m]) + cfg
        if m == L and t == 1 and cfg == 0 and (f == 0 or not feed):
            break
        cands = []
        if m == L and t == 0 and ((L > 0 and d[L - 1] == 1) or (L == 0 and f == 1)):
            cands.append((1.0 + V[base * 4 + f * 2 + 1], "top", (m, cfg, f, 1)))
        if m < L and wires < C and ((f == 1) if m == 0 else d[m - 1] == 1):
            cands.append((1.0 + V[(int(off[m + 1]) + cfg + int(pow3[m])) * 4 + f * 2], f"place c{m + 1} (first: writes bit {m + 1})",
                          (m + 1, cfg + int(pow3[m]), f, 0)))
        for i in range(m):
            pred = (f == 1) if i == 0 else (d[i - 1] == 1)
            p3 = int(pow3[i])
            if d[i] == 1:
                if pred:
                    cands.append((V[(base - p3) * 4 + f * 2 + t], f"remove c{i + 1} (free)", (m, cfg - p3, f, t)))
                else:
                    e = 0.5 * V[(base - p3) * 4 + f * 2 + t] + 0.5 * V[(base + p3) * 4 + f * 2 + t]
                    extra = V[(base + p3) * 4 + f * 2 + t] - V[(base - p3) * 4 + f * 2 + t]
                    cands.append((e, f"GHOST c{i + 1} (outcome 1 costs {extra:.3f} more)",
                                  ((m, cfg - p3, f, t), (m, cfg + p3, f, t))))
            elif d[i] == 2:
                if pred:
                    cands.append((V[(base - 2 * p3) * 4 + f * 2 + t], f"repay c{i + 1} (free)", (m, cfg - 2 * p3, f, t)))
            else:
                if pred and wires < C:
                    cands.append((1.0 + V[(base + p3) * 4 + f * 2 + t], f"re-place c{i + 1}", (m, cfg + p3, f, t)))
        if feed:
            if f == 1:
                cands.append((V[base * 4 + t], "erase carry-in (free)", (m, cfg, 0, t)))
            elif wires < C:
                cands.append((1.0 + V[base * 4 + 2 + t], "make carry-in", (m, cfg, 1, t)))
        cands.sort(key=lambda z: (z[0], 0 if "remove" in z[1] or "repay" in z[1] or "erase" in z[1] else 1))
        val, name, nxt = cands[0]
        if name.startswith("GHOST"):
            nxt = nxt[pick(ghosts)]
            ghosts += 1
        out.append((name, wires, val, (m, cfg, f, t)))
        m, cfg, f, t = nxt
    return out


def branches(n, C, feed, V, tabs):
    """the optimal schedule as a main line (every ghost returns 0) and, for each ghost, the repair when that one
    ghost returns 1 (later ghosts 0), printed until the schedule is back on the main line."""
    main = trace(n, C, feed, V, tabs, 0)
    print(" main line (every ghost measurement returns 0):")
    k = 0
    gpos = []
    for idx, (name, wires, val, st) in enumerate(main):
        tag = ""
        if name.startswith("GHOST"):
            tag = f"   <- ghost #{k}"
            gpos.append(idx)
            k += 1
        print(f"   {idx:3d} wires {wires}: {name}{tag}")
    for g, idx in enumerate(gpos):
        alt = trace(n, C, feed, V, tabs, lambda q, g=g: 1 if q == g else 0)
        later = {st: i for i, (_, _, _, st) in enumerate(main) if i > idx}
        print(f" repair when ghost #{g} ({main[idx][0].split(' (')[0]}) returns 1:")
        cost = 0
        for name, wires, val, st in alt[idx + 1:]:
            if st in later:
                print(f"      back on the main line at step {later[st]}; ANDs spent in this repair: {cost}")
                break
            if name.startswith("place") or name.startswith("re-place") or name.startswith("make") or name == "top":
                cost += 1
            print(f"      wires {wires}: {name}")
        else:
            print(f"      (ends without rejoining; ANDs spent after the ghost: {cost})")


def public_law():
    spec = importlib.util.spec_from_file_location("invadd", os.path.join(HERE, "r31_inv-wire_invadd.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.argv_backup = sys.argv
    sys.argv = [sys.argv[0], "none"]
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.argv = sys.argv_backup
    return mod


def game_price(n, clean, feed):
    """price of one row in the units of the generator's rule: feed rows have clean + 1 wires (the carry-in wire
    is allocated before the adder) and their first carry-in AND is priced with the feeds, not with the row."""
    if n <= 1:
        return 0.0
    if n == 2:
        return 1.0
    val, _, _ = solve(n, clean + (1 if feed else 0), feed)
    return val - (1.0 if feed else 0.0)


K2_ROWS = ([dict(j=j, n=13, feed=True, clean=c) for j, c in ((12, 2), (13, 3), (14, 4), (15, 5))]
           + [dict(j=j, n=13, feed=True, clean=6) for j in range(16, 25)]
           + [dict(j=j, n=13 - t, feed=False, clean=7)
              for j, t in zip(range(25, 36), (0, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10))])


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "table":
        nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 13
        P = public_law()
        for feed in (False, True):
            print(f"feed {int(feed)}: rows n, columns clean wires 0..9: optimum of the game / generator's rule")
            for n in range(4, nmax + 1):
                line = []
                for clean in range(0, 10):
                    g = game_price(n, clean, feed)
                    r = P.inv_cost(n, n - 2, clean, feed)
                    pub = r[0] if r is not None else float("inf")
                    line.append(f"{g:6.3f}/{pub:5.1f}")
                print(f"  n {n:2d}: " + "  ".join(line), flush=True)
    elif cmd == "row":
        n, C, feed = int(sys.argv[2]), int(sys.argv[3]), bool(int(sys.argv[4]))
        val, V, tabs = solve(n, C + (1 if feed else 0), feed)
        print(f"n {n} clean {C} feed {int(feed)}: optimum {val - (1.0 if feed else 0.0):.4f} (forward n - 1 = {n - 1})")
        branches(n, C + (1 if feed else 0), feed, V, tabs)
    elif cmd == "sched":
        P = public_law()
        rows = K2_ROWS if sys.argv[2] == "k2" else None
        if rows is None:
            raw = json.load(open(sys.argv[2]))
            rows = [dict(j=r.get("j", i), n=r["n"], feed=r["feed"], clean=46 - r["alloc"], tof=r["tof"])
                    for i, r in enumerate(raw)]
        extra = int(sys.argv[3]) if len(sys.argv) > 3 else 0
        tg = tp = tf = 0.0
        print("row  n feed clean  forward  rule   game   measured")
        for r in rows:
            clean = r["clean"] + extra
            g = game_price(r["n"], clean, r["feed"])
            pr = P.inv_cost(r["n"], r["n"] - 2, clean, r["feed"])
            pub = pr[0] if pr is not None else float("inf")
            best = min(g, pub)
            tg += best
            tp += pub
            tf += r["n"] - 1
            print(f"{r['j']:3d} {r['n']:2d} {int(r['feed']):4d} {clean:5d}  {r['n'] - 1:7d}  {pub:5.1f}  {g:6.3f}  "
                  f"{r.get('tof', '')}", flush=True)
        print(f"forward adders {tf:.1f}; generator's rule {tp:.1f} (penalty {tp - tf:.1f}); "
              f"best of game and rule {tg:.2f} (penalty {tg - tf:.2f})")
