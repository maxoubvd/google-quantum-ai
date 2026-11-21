#!/usr/bin/env python3
"""R91-INV-STEP part 3: width profile of fc18_M by value, on ALL inputs x = 1..3,133,452.

Own bit-packed value simulator of the kmx ops (X, CX, CCX, CSWAP, HMR, BIT_INVERT, PUSH/POP_CONDITION; phase ops change no
value). HMR outcomes in four modes: zero, one, random 1, random 2 (one coin a shot an HMR). A wire is LIVE after op i if
it is 1 in some shot of some mode. Width at op i = |live before op i UNION live after op i| (a fresh target counts at
its op, a measured wire at its HMR). Checks at the end: x restored and every scratch wire 0 in every shot and mode.
Executed Toffoli-like a shot per op: the exact coin probabilities of r91_inv_step_exec.py (per_op_p.json).
Op numbering = kickmix's (instructions after the 23 header lines), as R80 / R85 count.
Usage: python r91_inv_step_width.py FILE.kmx PER_OP_P.json OUTDIR
"""
import json
import os
import sys
import time

import numpy as np

TOF = ("CCX", "CCZ", "CSWAP")
MAXX = 3 * (2 ** 10 - 2) ** 2


def parse(path):
    ops = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("APPEND_TO_REGISTER") or s.startswith("REGISTER"):
            continue
        cond = None
        if " if " in s:
            s, c = s.split(" if ")
            cond = int(c.strip()[1:])
        p = s.split()
        ops.append((p[0], p[1:], cond))
    return ops


def qi(a):
    return int(a[1:])


def run_mode(ops, mode, nq, N, W, valid, x_words, lt64, last_use):
    q = np.zeros((nq, W), dtype=np.uint64)
    q[:22] = x_words
    bits = {}
    stack = []
    rng = np.random.default_rng(1000 + mode)
    ones = valid.copy()
    live = np.zeros(nq, dtype=bool)
    for k in range(nq):
        live[k] = bool(q[k].any())
    live_after = np.zeros((len(ops), nq), dtype=bool)
    flag = np.zeros(nq, dtype=np.int8)   # +1: wire == [x < 64], -1: wire == NOT [x < 64] (on valid shots)
    flag_after = np.zeros((len(ops), nq), dtype=np.int8)
    nlt = valid ^ lt64
    zero_w = np.zeros(W, dtype=np.uint64)
    for i, (nm, args, cond) in enumerate(ops):
        for b in last_use.get(i - 1, ()):
            bits.pop(b, None)
        m = ones
        if stack:
            m = stack[-1]
        if cond is not None:
            m = m & bits.get(cond, zero_w)
        written = ()
        if nm == "X":
            a = qi(args[0]); q[a] ^= m; written = (a,)
        elif nm == "CX":
            a, b = qi(args[0]), qi(args[1]); q[b] ^= q[a] & m; written = (b,)
        elif nm == "CCX":
            a, b, c = qi(args[0]), qi(args[1]), qi(args[2]); q[c] ^= q[a] & q[b] & m; written = (c,)
        elif nm == "CSWAP":
            c, a, b = qi(args[0]), qi(args[1]), qi(args[2])
            t = (q[a] ^ q[b]) & q[c] & m
            q[a] ^= t; q[b] ^= t; written = (a, b)
        elif nm == "HMR":
            a, b = qi(args[0]), int(args[1][1:])
            if mode == 0:
                coin = zero_w
            elif mode == 1:
                coin = valid
            else:
                coin = rng.integers(0, 2 ** 63, size=W, dtype=np.uint64) ^ (
                    rng.integers(0, 2, size=W, dtype=np.uint64) << np.uint64(63))
                coin &= valid
            old = bits.get(b, zero_w)
            bits[b] = (coin & m) | (old & ~m)
            q[a] &= ~m
            written = (a,)
        elif nm == "BIT_INVERT":
            b = int(args[0][1:])
            bits[b] = bits.get(b, zero_w) ^ m
        elif nm == "PUSH_CONDITION":
            stack.append(m)   # m = (enclosing mask) AND the pushed bit
        elif nm == "POP_CONDITION":
            stack.pop()
        for a in written:
            live[a] = bool(q[a].any())
            flag[a] = 1 if np.array_equal(q[a], lt64) else (-1 if np.array_equal(q[a], nlt) else 0)
        live_after[i] = live
        flag_after[i] = flag
    ok_x = bool(np.array_equal(q[:22], x_words))
    ok_s = not bool(q[22:].any())
    return live_after, flag_after, ok_x, ok_s


def main():
    kmx, pp, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    ops = parse(kmx)
    nq = 1 + max(qi(a) for nm, args, c in ops for a in args if a.startswith("q"))
    per = json.load(open(pp))
    per = [r for r in per if r["name"] not in ("APPEND_TO_REGISTER", "REGISTER")]
    assert len(per) == len(ops), (len(per), len(ops))
    pex = np.array([(r["p"] if (r["name"] in TOF and r["p"] is not None) else 0.0) for r in per])
    print("ops", len(ops), "qubits", nq, "exact executed Toffoli-like", round(pex.sum(), 4), flush=True)
    N = MAXX
    W = (N + 63) // 64
    xs = np.arange(1, N + 1, dtype=np.uint64)
    valid = np.zeros(W, dtype=np.uint64)
    full = N // 64
    valid[:full] = np.uint64(0xFFFFFFFFFFFFFFFF)
    if N % 64:
        valid[full] = np.uint64((1 << (N % 64)) - 1)
    x_words = np.zeros((22, W), dtype=np.uint64)
    pad = np.zeros(W * 64, dtype=np.uint64)
    pad[:N] = xs
    for j in range(22):
        bitsj = ((pad >> np.uint64(j)) & np.uint64(1)).astype(np.uint8)
        x_words[j] = np.packbits(bitsj.reshape(W, 64)[:, ::-1], axis=1).view(">u8").astype(np.uint64).ravel()
    # sanity: bit k of word w = shot 64 w + k
    assert ((x_words[0][0] >> np.uint64(0)) & np.uint64(1)) == 1  # x = 1 has bit 0 set
    # indicator [x < 64] packed, for the f3 search
    lt64 = np.zeros(W, dtype=np.uint64)
    small = ((pad >= 1) & (pad < 64)).astype(np.uint8)
    lt64[:] = np.packbits(small.reshape(W, 64)[:, ::-1], axis=1).view(">u8").astype(np.uint64).ravel()
    last = {}
    for i, (nm, args, cond) in enumerate(ops):
        if cond is not None:
            last[cond] = i
        if nm == "HMR":
            last[int(args[1][1:])] = max(last.get(int(args[1][1:]), i), i)
        if nm == "BIT_INVERT":
            last[int(args[0][1:])] = i
    last_use = {}
    for b, i in last.items():
        last_use.setdefault(i, []).append(b)
    lives = []
    flags = []
    for mode in range(4):
        la, fa, okx, oks = run_mode(ops, mode, nq, N, W, valid, x_words, lt64, last_use)
        print("mode", mode, "x restored", okx, "scratch clean", oks, "t", round(time.time() - t0, 1), flush=True)
        lives.append(la)
        flags.append(fa)
    L = lives[0] | lives[1] | lives[2] | lives[3]
    # the class flag [x < 64] (route N's removed wire): a wire that equals it (or its complement) in every mode
    F = np.where((flags[0] == flags[1]) & (flags[1] == flags[2]) & (flags[2] == flags[3]), flags[0], 0)
    np.save(os.path.join(outdir, "flag_lt64.npy"), F)
    L0 = np.zeros(nq, dtype=bool)
    L0[:22] = True
    before = np.vstack([L0[None, :], L[:-1]])
    U = before | L
    width = U.sum(axis=1)
    np.save(os.path.join(outdir, "live_after.npy"), L)
    np.save(os.path.join(outdir, "width.npy"), width)
    hist = {int(w): int((width == w).sum()) for w in np.unique(width)}
    print("width histogram", hist)
    for w in (39, 38, 37):
        sel = width == w
        print("width %d: ops %d, executed Toffoli-like a shot %.4f, static Toffoli-like %d" % (
            w, int(sel.sum()), float(pex[sel].sum()),
            int(sum(1 for i in np.flatnonzero(sel) if ops[i][0] in TOF))))
    # runs at the peak
    runs = []
    i = 0
    n = len(ops)
    while i < n:
        if width[i] == 39:
            j = i
            while j + 1 < n and width[j + 1] == 39:
                j += 1
            touched = set()
            for k in range(i, j + 1):
                for a in ops[k][1]:
                    if a.startswith("q"):
                        touched.add(qi(a))
            liv = set(np.flatnonzero(U[i:j + 1].all(axis=0)).tolist())
            opening = sorted(set(np.flatnonzero(U[i] & ~(before[i] if i > 0 else L0)).tolist())) if i > 0 else []
            # the wire newly live at op i (fresh target) or, if none, the wire whose life started latest before i
            closing = sorted(set(np.flatnonzero(L[j] & ~L[j + 1]).tolist())) if j + 1 < n else []
            passengers = sorted(liv - touched)
            f3 = sorted(np.flatnonzero((F[i:j + 1] != 0).all(axis=0)).tolist())
            runs.append({"start": i, "end": j, "n_ops": j - i + 1, "exec_tof": float(pex[i:j + 1].sum()),
                         "f3_wires": f3,
                         "static_tof": int(sum(1 for k in range(i, j + 1) if ops[k][0] in TOF)),
                         "opening": opening, "closing": closing, "passengers": passengers,
                         "first_op": " ".join([ops[i][0]] + ops[i][1]), "last_op": " ".join([ops[j][0]] + ops[j][1])})
            i = j + 1
        else:
            i += 1
    print("runs at 39:", len(runs), "ops in them", sum(r["n_ops"] for r in runs))
    json.dump({"hist": hist, "runs": runs, "n_ops": n}, open(os.path.join(outdir, "runs39.json"), "w"), indent=1)
    # birth op of each wire's current life (for 'latest arrival') per run start
    birth = np.zeros(nq, dtype=np.int64)
    births_at = []
    for k in range(n):
        if k == 0:
            newl = L[0] & ~L0
        else:
            newl = L[k] & ~L[k - 1]
        birth[newl] = k
        births_at.append(birth.copy())
    for r in runs:
        b = births_at[r["start"]]
        livs = np.flatnonzero(U[r["start"]])
        order = sorted(livs.tolist(), key=lambda w: -b[w])
        r["latest_arrivals"] = [(int(w), int(b[w])) for w in order[:3]]
    json.dump({"hist": hist, "runs": runs, "n_ops": n}, open(os.path.join(outdir, "runs39.json"), "w"), indent=1)
    print("time", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
