"""chem_loaders: exact gauge moves on an angle-word table (the four word forms) and a search for cheaper loaders.

The verifier decodes row words on the Givens tree: v[p], v[c] = cos(th) v[p], sin(th) v[p] for every edge (p, c).
Flipping the sign of cos on edge g = (x, c) (word g -> word g + 128 after the compensation below, i.e. the
bit7 xor bit6 "carrier" bit of that word flips) is exact when the sign change of v[x] it causes is undone:
negate (word -> -word mod 256) the edge that created v[x] and every earlier edge out of x. For an edge out of the
root there is nothing to negate before the first edge and the row's vector changes sign as a whole, which no
metric sees. These 75 moves per row generate every exact form change (theta, theta + pi, -theta, pi - theta).

flip(words, edges, g, rows)   apply the move for word g on the given rows (boolean mask or None = all rows)
check(words0, words1, edges)  largest deviation of the decoded vectors, up to the sign of each row
search(...)                   greedy: uniform moves (all rows) and carrier-restoring moves (rows where a low
                              monomial is 1), accepted when the exact price of the loaders drops.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402


def tree(edges):
    edges = [tuple(int(x) for x in e) for e in edges]
    made_by = {c: t for t, (p, c) in enumerate(edges)}          # node -> edge that created it
    out = {}
    for t, (p, c) in enumerate(edges):
        out.setdefault(p, []).append(t)
    return edges, made_by, out


def flip(words, edges, g, rows=None):
    """New table with the cos sign of word g flipped on `rows` and the exact compensation applied."""
    edges, made_by, out = tree(edges)
    w = np.array(words, dtype=np.int64, copy=True)
    sel = slice(None) if rows is None else np.asarray(rows, bool)
    x = edges[g][0]
    w[sel, g] ^= 0x80
    neg = [t for t in out[x] if t < g]
    if x in made_by:
        neg.append(made_by[x])
    for t in neg:
        w[sel, t] = (-w[sel, t]) & 0xFF
    return w


def decode(words, edges, brot=8):
    words = np.asarray(words, np.int64)
    rows = words.shape[0]
    n = int(np.max(edges)) + 1
    v = np.zeros((rows, n))
    v[:, 0] = 1.0
    th = np.pi * words / (1 << (brot - 1))
    for t, (p, c) in enumerate(edges):
        pv, cv = v[:, p].copy(), v[:, c].copy()
        v[:, p] = np.cos(th[:, t]) * pv - np.sin(th[:, t]) * cv
        v[:, c] = np.sin(th[:, t]) * pv + np.cos(th[:, t]) * cv
    return v


def check(words0, words1, edges):
    a, b = decode(words0, edges), decode(words1, edges)
    return float(np.max(np.minimum(np.abs(a - b).max(axis=1), np.abs(a + b).max(axis=1))))


def _price(words, ks, effort, workers, lb, cache):
    """Toffolis of the circuits ks (cached on the bytes of the two batches each circuit depends on)."""
    todo, keys = [], {}
    for k in ks:
        lo = max(0, 5 * k - 5)
        key = (k, words[:, lo:5 * k + 5].tobytes())
        keys[k] = key
        if key not in cache:
            todo.append(k)
    if todo:
        jobs = [(words, k, lb, 0, effort, False) for k in todo]
        if workers > 1 and len(jobs) > 1:
            from concurrent.futures import ProcessPoolExecutor
            with ProcessPoolExecutor(min(3, workers)) as ex:
                res = list(ex.map(core.compile_circuit, jobs, chunksize=1))
        else:
            res = [core.compile_circuit(j) for j in jobs]
        for k, r in zip(todo, res):
            cache[keys[k]] = r["toffoli"]
    return {k: cache[keys[k]] for k in ks}


def search(words, edges, out=None, rounds=1, effort=2, workers=3, lb=4, log=print):
    """Greedy search over exact form changes. Returns (new words, per-pass Toffolis before, after)."""
    words = np.asarray(words, np.int64)
    A, ncol = words.shape
    K = -(-ncol // 5)
    cache = {}
    cur = _price(words, range(K), effort, workers, lb, cache)
    start = sum(cur.values())
    log(f"start: {start} per pass")
    low = np.arange(A) & 15
    monos = [m for m in range(1, 16) if m & (m - 1)]
    edges_l, made_by, outs = tree(edges)
    for rnd in range(rounds):
        gained = 0
        for g in range(ncol):
            x = edges_l[g][0]
            touched = {g} | {t for t in outs[x] if t < g} | ({made_by[x]} if x in made_by else set())
            ks = sorted({k for t in touched for k in (t // 5, t // 5 + 1) if k < K})
            cands = [None]
            # carrier-restoring rows: the low monomials other words of this slot carry
            car = ((words[:, g] >> 7) ^ (words[:, g] >> 6)) & 1
            if not car.any():
                for m in monos:
                    cands.append((low & m) == m)
            best = None
            for rows_mask in cands:
                w2 = flip(words, edges, g, rows_mask)
                new = _price(w2, ks, effort, workers, lb, cache)
                d = sum(new[k] - cur[k] for k in ks)
                if d < 0 and (best is None or d < best[0]):
                    best = (d, w2, new)
            if best is not None:
                d, w2, new = best
                dev = check(words, w2, edges)
                assert dev < 1e-9, dev
                words = w2
                cur.update(new)
                gained -= d
                log(f"round {rnd} word {g}: {d} -> {sum(cur.values())} per pass")
        if not gained:
            break
    total = sum(cur.values())
    if out:
        np.save(out, words)
    log(f"done: {start} -> {total} per pass" + (f"; words written to {out}" if out else ""))
    return words, start, total
