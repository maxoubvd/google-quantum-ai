#!/usr/bin/env python3
"""r46_dqi_scan_rank.py -- R46-DQI-SCAN: redundant wires of the DQI decoder at its peaks, on sampled valid shots.

P1's bit-parallel value model (tools/dqi_pp_act.py: valid verifier samples, the verifier's weight law, random HMR
outcomes; phases ignored), with a snapshot of every wire at chosen op boundaries (the state just BEFORE op index oi).
Per site: L = wires that are nonzero on at least one sample, the GF(2) rank of [1 | those wires] over the samples and
every affine relation among them (free columns of the reduced echelon form; column order: constant, wires 0..901, so
a dependent wire is written on wires of lower index). A site with nullity 0 is TIGHT on the samples (exact: no affine
relation exists on the reachable set); a relation found on the samples is a candidate that needs a proof.
Usage:
  python r46_dqi_scan_rank.py sweep <dqi.kmb> <out.json> [stride ops=35000] [runs=20] [seed=5] [offset=0]
  python r46_dqi_scan_rank.py sites <dqi.kmb> <out.json> <sites.npy> [runs=20] [seed=5] [keep.npz]
"""
import json
import os
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
from dqi_pp_act import pack, v, km, _rng  # noqa: E402

ALL = np.uint64(0xFFFFFFFFFFFFFFFF)


@njit(cache=True)
def simulate_snapop(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, site_oi, snap):
    vv = vinit.copy()
    bits = binit.copy()
    s = np.uint64(seed) * np.uint64(0x9E3779B97F4A7C15) + np.uint64(1)
    stack = np.zeros(256, np.uint64)
    sp = 0
    active = ALL
    i3 = 0
    i2 = 0
    i1 = 0
    ib = 0
    ic = 0
    js = 0
    ns = site_oi.shape[0]
    for oi in range(ops.shape[0]):
        while js < ns and site_oi[js] == oi:
            snap[js, :] = vv
            js += 1
        op = ops[oi]
        if op == 0x81:
            c = cd[ic]
            ic += 1
            stack[sp] = active
            sp += 1
            active = active & bits[c]
            continue
        base = op & 0x7F
        if base == 0x01:
            if sp > 0:
                sp -= 1
                active = stack[sp]
            continue
        nqk = (op & 0x30) >> 4
        a = 0
        b = 0
        t = 0
        if nqk == 3:
            a = q3a[i3]
            b = q3b[i3]
            t = q3c[i3]
            i3 += 1
        elif nqk == 2:
            a = q2a[i2]
            b = q2b[i2]
            i2 += 1
        elif nqk == 1:
            a = q1[i1]
            i1 += 1
        bit = 0
        if (op & 0x40) != 0 and base != 0x06 and base != 0x07:
            bit = bt[ib]
            ib += 1
        m = active
        if (op & 0x80) != 0:
            m = m & bits[cd[ic]]
            ic += 1
        if base == 0x11:
            vv[a] ^= m
        elif base == 0x12:
            vv[a] &= ~m
        elif base == 0x21:
            vv[b] ^= vv[a] & m
        elif base == 0x22:
            d = (vv[a] ^ vv[b]) & m
            vv[a] ^= d
            vv[b] ^= d
        elif base == 0x31:
            vv[t] ^= vv[a] & vv[b] & m
        elif base == 0x32:
            d = (vv[b] ^ vv[t]) & vv[a] & m
            vv[b] ^= d
            vv[t] ^= d
        elif base == 0x50:
            s = _rng(s)
            bits[bit] = (bits[bit] & ~m) | (s & m)
            vv[a] &= ~m
        elif base == 0x40:
            bits[bit] &= ~m
        elif base == 0x41:
            bits[bit] |= m
        elif base == 0x42:
            bits[bit] ^= m
    return js


@njit(cache=True)
def rref(M, ncols):
    """Reduced row echelon form in place; M (n, W) uint64, bit c of a row = column c. Returns pivot columns."""
    n = M.shape[0]
    W = M.shape[1]
    piv = np.full(ncols, -1, np.int64)
    row = 0
    one = np.uint64(1)
    for col in range(ncols):
        w = col >> 6
        bb = np.uint64(col & 63)
        p = -1
        for r in range(row, n):
            if (M[r, w] >> bb) & one:
                p = r
                break
        if p < 0:
            continue
        if p != row:
            for x in range(W):
                tmp = M[row, x]
                M[row, x] = M[p, x]
                M[p, x] = tmp
        for r in range(n):
            if r != row and ((M[r, w] >> bb) & one):
                for x in range(W):
                    M[r, x] ^= M[row, x]
        piv[row] = col
        row += 1
        if row >= n:
            break
    return piv[:row]


def unpack_site(S):
    """S (R, nq) uint64 -> (64 R, nq) uint8 bits."""
    R, nq = S.shape
    W = np.unpackbits(S.copy().view(np.uint8).reshape(R, nq * 8), axis=1, bitorder="little").reshape(R, nq, 64)
    return W.transpose(0, 2, 1).reshape(R * 64, nq)


def relations(Wb):
    """Wb (N, nq) uint8. Returns live wires, rank, list of relations (each a list: -1 = constant, else wire)."""
    N, nq = Wb.shape
    live = np.flatnonzero(Wb.any(axis=0))
    A = np.concatenate([np.ones((N, 1), np.uint8), Wb[:, live]], axis=1)
    nc = A.shape[1]
    pad = (-nc) % 64
    if pad:
        A = np.concatenate([A, np.zeros((N, pad), np.uint8)], axis=1)
    M = np.ascontiguousarray(np.packbits(A, axis=1, bitorder="little")).view(np.uint64).copy()
    piv = rref(M, nc)
    rank = len(piv)
    pivset = set(piv.tolist())
    rels = []
    if rank < nc:
        cols = [-1] + live.tolist()
        for f in range(nc):
            if f in pivset:
                continue
            w, bb = f >> 6, np.uint64(f & 63)
            rows = np.flatnonzero((M[:rank, w] >> bb) & np.uint64(1))
            rels.append([cols[f]] + [cols[int(piv[i])] for i in rows])
    return live, rank, rels


def site_labels(c, sites):
    tp = c.tag_positions()
    tags = c.tags
    out = []
    stack = []
    k = 0
    for s in sites:
        while k < len(tp) and tp[k] <= s:
            t = tags[k]
            if t.startswith("@enter:"):
                stack.append(t[7:])
            elif t.startswith("@exit:") and stack:
                nm = t[6:]
                if nm in stack:
                    while stack and stack[-1] != nm:
                        stack.pop()
                    stack.pop()
            k += 1
        out.append((stack[0] if stack else "-", stack[-1] if stack else "-", len(stack)))
    return out


def main():
    t0 = time.time()
    mode = sys.argv[1]
    c = kmb.read(sys.argv[2])
    out = sys.argv[3]
    ops = np.ascontiguousarray(np.asarray(c.ops))
    n = len(ops)
    keep = None
    if mode == "sweep":
        S = int(sys.argv[4]) if len(sys.argv) > 4 else 35000
        R = int(sys.argv[5]) if len(sys.argv) > 5 else 20
        seed = int(sys.argv[6]) if len(sys.argv) > 6 else 5
        off = int(sys.argv[7]) if len(sys.argv) > 7 else 0
        sites = np.arange(off % S, n, S, dtype=np.int64)
    else:
        sites = np.sort(np.load(sys.argv[4]).astype(np.int64))
        R = int(sys.argv[5]) if len(sys.argv) > 5 else 20
        seed = int(sys.argv[6]) if len(sys.argv) > 6 else 5
        keep = sys.argv[7] if len(sys.argv) > 7 else None
    ns = len(sites)
    nq = c.num_qubits
    print("read %.0fs: %d ops, %d qubits; %d sites, %d samples" % (time.time() - t0, n, nq, ns, 64 * R), flush=True)
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    wprob = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    rng = np.random.default_rng(seed)
    snaps = np.zeros((R, ns, nq), np.uint64)
    for r in range(R):
        vinit = np.zeros(nq, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        syn, mbuc, _ = v.sample_shots(spec, tables, rng, 64, wprob)
        vinit[:840] = pack(v.to_bits(syn, spec.degree))
        mb = pack(v.to_bits(mbuc, spec.degree))
        binit[:len(mb)] = mb
        js = simulate_snapop(ops, *arr, vinit, binit, seed * 1000 + r, sites, snaps[r])
        assert js == ns, (js, ns)
        if r % 4 == 3 or r == R - 1:
            print("run %d done (%.0fs)" % (r, time.time() - t0), flush=True)
    if keep:
        np.savez(keep, sites=sites, snaps=snaps)
    labels = site_labels(c, sites.tolist())
    res = []
    for j in range(ns):
        Wb = unpack_site(snaps[:, j, :])
        live, rank, rels = relations(Wb)
        res.append({"op": int(sites[j]), "stage": labels[j][0], "mark": labels[j][1], "L": int(len(live)),
                    "rank": int(rank), "null": int(len(live) + 1 - rank), "zero": [int(x) for x in
                    np.setdiff1d(np.arange(nq), live)] if len(live) >= nq - 40 else [],
                    "rels": [r[:40] for r in rels[:60]], "rel_len": [len(r) for r in rels[:60]]})
        if j % 200 == 199:
            print("site %d / %d (%.0fs)" % (j + 1, ns, time.time() - t0), flush=True)
    json.dump({"kmb": sys.argv[2], "runs": R, "seed": seed, "nq": nq, "sites": res}, open(out, "w"))
    # summary by stage
    by = {}
    for r in res:
        d = by.setdefault(r["stage"], [0, 0, 0, 0, 0, 0])
        d[0] += 1
        d[1] = max(d[1], r["L"])
        d[2] = max(d[2], r["L"] - r["null"])
        if r["L"] >= nq - 2:
            d[3] += 1
            if r["null"] == 0:
                d[4] += 1
        if r["null"] == 0:
            d[5] += 1
    print("%-16s %6s %6s %10s %12s %10s %10s" % ("stage", "sites", "max L", "max L-null", "sites L>=Q-2", "tight", "null=0"))
    for k, d in by.items():
        print("%-16s %6d %6d %10d %12d %10d %10d" % (k[:16], d[0], d[1], d[2], d[3], d[4], d[5]))
    print("done %.0fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
