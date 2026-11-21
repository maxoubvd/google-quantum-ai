#!/usr/bin/env python3
"""R32 seat INV-CLEAN (2026-11-13): the undo-adder game of tools/r32_inv-clean_pebble.py with one more pebble,
the row's own sign wire. Own code.

Variant "sign inside its row" (feed rows 12..15 of the divider clean-up): the row's sign wire s_j is needed for
the first frame (before any carry), for making / erasing the carry-in (the feed formula reads it) and for the
last frame (after the last write). Between those it idles. If it has a small exact rebuild program (gcost ANDs;
s12: 1 AND from v10, v11 and the top target bit, whose two frames cancel so that bit can stay unframed) it can be
X-measured after the first frame (its repair is that program as a phase: for s12 a quadratic form, free) and
rebuilt when the carry-in or the last frame needs it. Rules added to the game:
  rebuild sign   gcost ANDs, needs a free wire          measure sign   free
  make / erase carry-in need the sign on a wire
  start: sign on a wire (today's rebuild, priced with the sign rebuilds); end: sign on a wire (top bit of U_j)
Wires: carries + carry-in + sign <= clean + 2 (today both the sign and the carry-in are allocated before the adder
counts its clean wires).

  r32_inv-clean_sign.py one N CLEAN GCOST        optimum and the schedule on the all-0 and all-1 outcome paths
  r32_inv-clean_sign.py many "N CLEAN GCOST" ...  several, one process each (logs sign_N_CLEAN_GCOST.log)
"""
import importlib.util
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import numba as nb
    HAVE_NUMBA = True
except Exception:
    HAVE_NUMBA = False

INF = 1e18


def _sweep_sign_py(V, L, C, gcost, pow3, off, d):
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


_sweep_sign = nb.njit(cache=True)(_sweep_sign_py) if HAVE_NUMBA else _sweep_sign_py


def solve_sign(n, clean, gcost, tol=1e-10, max_sweeps=5000):
    """feed row with the sign pebble. Returns the row's price in the generator's units (the first carry-in AND and
    today's sign rebuild are priced elsewhere), the value array and the index tables."""
    L = n - 2
    C = clean + 2
    pow3 = np.array([3 ** i for i in range(L + 2)], dtype=np.int64)
    off = np.zeros(L + 2, dtype=np.int64)
    for m in range(L + 1):
        off[m + 1] = off[m] + pow3[m]
    V = np.zeros(int(off[L + 1]) * 8, dtype=np.float64)
    d = np.zeros(L + 1, dtype=np.int64)
    for _ in range(max_sweeps):
        delta = _sweep_sign(V, L, C, float(gcost), pow3, off, d)
        if delta < tol:
            break
    start = V[((0 * 2 + 1) * 2 + 0) * 2 + 0]       # m = 0, sign on a wire, no carry-in, top not done
    return float(start) - 1.0, V, (L, pow3, off)


def trace_sign(n, clean, gcost, V, tabs, outcome=0, limit=400):
    L, pow3, off = tabs
    C = clean + 2
    m, cfg, g, f, t = 0, 0, 1, 0, 0
    out = []
    ghosts = 0
    pick = outcome if callable(outcome) else (lambda k: outcome)

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
            cands.append((1.0 + V[sid(m + 1, cfg + int(pow3[m]), g, f, 0)],
                          f"place c{m + 1} (first: writes bit {m + 1})", (m + 1, cfg + int(pow3[m]), g, f, 0)))
        for i in range(m):
            pred = (f == 1) if i == 0 else (d[i - 1] == 1)
            p3 = int(pow3[i])
            if d[i] == 1:
                if pred:
                    cands.append((V[sid(m, cfg - p3, g, f, t)], f"remove c{i + 1} (free)", (m, cfg - p3, g, f, t)))
                else:
                    e = 0.5 * V[sid(m, cfg - p3, g, f, t)] + 0.5 * V[sid(m, cfg + p3, g, f, t)]
                    extra = V[sid(m, cfg + p3, g, f, t)] - V[sid(m, cfg - p3, g, f, t)]
                    cands.append((e, f"GHOST c{i + 1} (outcome 1 costs {extra:.3f} more)",
                                  ((m, cfg - p3, g, f, t), (m, cfg + p3, g, f, t))))
            elif d[i] == 2:
                if pred:
                    cands.append((V[sid(m, cfg - 2 * p3, g, f, t)], f"repay c{i + 1} (free)",
                                  (m, cfg - 2 * p3, g, f, t)))
            else:
                if pred and wires < C:
                    cands.append((1.0 + V[sid(m, cfg + p3, g, f, t)], f"re-place c{i + 1}", (m, cfg + p3, g, f, t)))
        if g == 1:
            if f == 1:
                cands.append((V[sid(m, cfg, g, 0, t)], "erase carry-in (free)", (m, cfg, g, 0, t)))
            elif wires < C:
                cands.append((1.0 + V[sid(m, cfg, g, 1, t)], "make carry-in", (m, cfg, g, 1, t)))
            cands.append((V[sid(m, cfg, 0, f, t)], "MEASURE SIGN (free)", (m, cfg, 0, f, t)))
        elif wires < C:
            cands.append((gcost + V[sid(m, cfg, 1, f, t)], f"REBUILD SIGN ({gcost} AND)", (m, cfg, 1, f, t)))
        cands.sort(key=lambda z: (z[0], 0 if "free" in z[1] else 1))
        val, name, nxt = cands[0]
        if name.startswith("GHOST"):
            nxt = nxt[pick(ghosts)]
            ghosts += 1
        out.append((name, wires, val, (m, cfg, g, f, t)))
        m, cfg, g, f, t = nxt
    return out


def branches(n, clean, gcost, V, tabs):
    main = trace_sign(n, clean, gcost, V, tabs, 0)
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
    for gi, idx in enumerate(gpos):
        alt = trace_sign(n, clean, gcost, V, tabs, lambda q, gi=gi: 1 if q == gi else 0)
        later = {st: i for i, (_, _, _, st) in enumerate(main) if i > idx}
        print(f" repair when ghost #{gi} ({main[idx][0].split(' (')[0]}) returns 1:")
        cost = 0.0
        for name, wires, val, st in alt[idx + 1:]:
            if st in later:
                print(f"      back on the main line at step {later[st]}; ANDs spent in this repair: {cost}")
                break
            if name.startswith("place") or name.startswith("re-place") or name.startswith("make") or name == "top":
                cost += 1
            if name.startswith("REBUILD"):
                cost += gcost
            print(f"      wires {wires}: {name}")
        else:
            print(f"      (ends without rejoining; ANDs spent after the ghost: {cost})")


def base_price(n, clean):
    spec = importlib.util.spec_from_file_location("pebble", os.path.join(HERE, "r32_inv-clean_pebble.py"))
    mod = importlib.util.module_from_spec(spec)
    keep = sys.argv
    sys.argv = [keep[0], "none"]
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.argv = keep
    return mod.game_price(n, clean, True)


if __name__ == "__main__":
    if sys.argv[1] == "one":
        n, clean, gcost = int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
        val, V, tabs = solve_sign(n, clean, gcost)
        base = base_price(n, clean)
        print(f"n {n} clean {clean}, sign rebuild {gcost} AND: with the sign pebble {val:.4f}; "
              f"sign held (game) {base:.4f}; change {val - base:+.4f}", flush=True)
        branches(n, clean, gcost, V, tabs)
    elif sys.argv[1] == "many":
        procs = []
        for spec_ in sys.argv[2:]:
            a = spec_.split()
            name = "sign_" + "_".join(a) + ".log"
            procs.append((name, subprocess.Popen([sys.executable, os.path.abspath(__file__), "one"] + a,
                                                 stdout=open(name, "w"), stderr=subprocess.STDOUT)))
        for name, p in procs:
            p.wait()
            print(name, "->", open(name).readline().strip(), flush=True)
