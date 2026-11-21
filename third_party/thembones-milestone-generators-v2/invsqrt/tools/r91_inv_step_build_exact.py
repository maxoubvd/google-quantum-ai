#!/usr/bin/env python3
"""R91-INV-C2-SPLITDEBT: exact re-check of a split lead on ALL inputs x = 1..3,133,452.

For each pair (debt b, point t): simulate every input (the census simulator, bit-packed, one shot an input), record
f_b (the measured wire's value just before b's HMR), the wires live right after the HMR (measured wire excluded) and the
wires right after op t. Solve f_b = g(H wires) + h(t wires), g and h of degree <= 2, EXACTLY on all inputs:
  sample rows -> Gauss-Jordan -> a particular solution -> evaluated on ALL inputs -> failing inputs added -> repeat,
until the solution matches f_b on every input (ALIVE, the solution is printed) or the system is inconsistent on the
enlarged sample (DEAD: a proof on all inputs). Coins: COINS=one (all 1) or a seed (random coins, one a shot an HMR);
run both: at an unconditioned point outside every block the state must not depend on the coins.
Usage: [COINS=one|SEED=k] python r91_inv_step_build_exact.py FILE.kmx EXEC_BLOCKS.json OUT.json b:t [b:t ...]
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r91_inv_step_width as Wd  # noqa: E402

TOF = ("CCX", "CCZ", "CSWAP")
COINS_ONE = os.environ.get("COINS", "") == "one"
SEED = int(os.environ.get("SEED", "11"))


def bit_of(words, idx):
    """bits of packed words (bit k of word w = shot 64 w + k) at shot indices idx -> bool array"""
    return ((words[idx >> 6] >> (idx & 63).astype(np.uint64)) & np.uint64(1)).astype(bool)


def simulate_all(ops, want_h, want_t):
    """want_h: set of coin ids (record f and H-state at their HMR); want_t: set of op indices (record state after op).
    Returns dict f[b], H[b] (list of (wire, words)), T[t] (39 x W words), selT[t] (W words), N."""
    N = Wd.MAXX
    W = (N + 63) // 64
    xs = np.arange(1, N + 1, dtype=np.uint64)
    pad = np.zeros(W * 64, dtype=np.uint64)
    pad[:N] = xs
    valid = np.zeros(W, dtype=np.uint64)
    valid[: N // 64] = np.uint64(0xFFFFFFFFFFFFFFFF)
    if N % 64:
        valid[N // 64] = np.uint64((1 << (N % 64)) - 1)

    def packw(bits):
        return np.packbits(bits.astype(np.uint8).reshape(W, 64)[:, ::-1], axis=1).view(">u8").astype(np.uint64).ravel()

    q = np.zeros((39, W), dtype=np.uint64)
    x0 = []
    for j in range(22):
        q[j] = packw(((pad >> np.uint64(j)) & np.uint64(1)).astype(bool))
        x0.append(q[j].copy())
    bits = {}
    stack = []
    zero = np.zeros(W, dtype=np.uint64)
    rng = np.random.default_rng(SEED)
    F, H, T, S = {}, {}, {}, {}
    for i, (nm, args, cond) in enumerate(ops):
        mm = stack[-1] if stack else valid
        if cond is not None:
            mm = mm & bits.get(cond, zero)
        if nm == "X":
            q[Wd.qi(args[0])] ^= mm
        elif nm == "CX":
            q[Wd.qi(args[1])] ^= q[Wd.qi(args[0])] & mm
        elif nm in TOF:
            if nm == "CCX":
                q[Wd.qi(args[2])] ^= q[Wd.qi(args[0])] & q[Wd.qi(args[1])] & mm
            elif nm == "CSWAP":
                c, a, b2 = Wd.qi(args[0]), Wd.qi(args[1]), Wd.qi(args[2])
                t = (q[a] ^ q[b2]) & q[c] & mm
                q[a] ^= t
                q[b2] ^= t
            if i in want_t:
                T[i] = q.copy()
                S[i] = mm.copy()
        elif nm == "HMR":
            a, b = Wd.qi(args[0]), int(args[1][1:])
            if b in want_h:
                F[b] = q[a].copy()
                H[b] = [(w, q[w].copy()) for w in range(39) if w != a]
            if COINS_ONE:
                coin = valid.copy()
            else:
                coin = rng.integers(0, 2 ** 63, size=W, dtype=np.uint64) ^ (
                    rng.integers(0, 2, size=W, dtype=np.uint64) << np.uint64(63))
                coin &= valid
            old = bits.get(b, zero)
            bits[b] = (coin & mm) | (old & ~mm)
            q[a] &= ~mm
        elif nm == "BIT_INVERT":
            b = int(args[0][1:])
            bits[b] = bits.get(b, zero) ^ mm
        elif nm == "PUSH_CONDITION":
            stack.append(mm)
        elif nm == "POP_CONDITION":
            stack.pop()
    ok_x = all(np.array_equal(q[j] & valid, x0[j] & valid) for j in range(22))
    clean = not (q[22:] & valid).any()
    return F, H, T, S, N, W, valid, ok_x, clean


def solve(cols_s, f_s):
    """Gauss-Jordan over GF(2): cols_s (K, M) bool, f_s (K,) bool. Returns (consistent, coef bool[M] or None, rank)."""
    K, M = cols_s.shape
    A = np.concatenate([cols_s, f_s[:, None]], axis=1)
    nb = M + 1
    nw = (nb + 63) // 64
    padc = np.zeros((K, nw * 64), dtype=bool)
    padc[:, :nb] = A
    R = np.packbits(padc.reshape(K, nw, 64)[:, :, ::-1], axis=2).view(">u8").astype(np.uint64).reshape(K, nw)
    used = np.zeros(K, dtype=bool)
    piv = {}
    for c in range(M):
        w, s = divmod(c, 64)
        colbits = ((R[:, w] >> np.uint64(s)) & np.uint64(1)).astype(bool)
        cand = np.flatnonzero(colbits & ~used)
        if len(cand) == 0:
            continue
        p = cand[0]
        used[p] = True
        piv[c] = p
        others = np.flatnonzero(colbits)
        others = others[others != p]
        if len(others):
            R[others] ^= R[p]
    wt, st = divmod(M, 64)
    tb = ((R[:, wt] >> np.uint64(st)) & np.uint64(1)).astype(bool)
    if tb[~used].any():
        return False, None, len(piv)
    coef = np.zeros(M, dtype=bool)
    for c, p in piv.items():
        coef[c] = tb[p]
    return True, coef, len(piv)


def main():
    kmx, blocksp, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    pairs = [(int(s.split(":")[0]), int(s.split(":")[1])) for s in sys.argv[4:]]
    t0 = time.time()
    ops = Wd.parse(kmx)
    F, H, T, S, N, W, valid, ok_x, clean = simulate_all(ops, {b for b, _ in pairs}, {t for _, t in pairs})
    print("ALL inputs simulated: N", N, "x restored", ok_x, "scratch clean", clean, "coins", "one" if COINS_ONE else
          "seed %d" % SEED, "t", round(time.time() - t0, 1), flush=True)
    rng = np.random.default_rng(12345)
    out = []
    for b, t in pairs:
        f = F[b] & valid
        selc = int(np.unpackbits(S[t].view(np.uint8)).sum())
        hl = [(w, v & valid) for w, v in H[b] if (v & valid).any()]
        tl = [(w, T[t][w] & valid) for w in range(39) if (T[t][w] & valid).any()]
        # monomial list: H part then t part; each as (kind, wires)
        mons = [("H", ())] + [("H", (w,)) for w, _ in hl] + [("H", (hl[a][0], hl[c][0])) for a in range(len(hl))
                                                            for c in range(a + 1, len(hl))]
        mons += [("T", (w,)) for w, _ in tl] + [("T", (tl[a][0], tl[c][0])) for a in range(len(tl))
                                               for c in range(a + 1, len(tl))]
        Hd = {w: v for w, v in hl}
        Td = {w: v for w, v in tl}

        def mon_words(m):
            kind, ws = m
            D = Hd if kind == "H" else Td
            if not ws:
                return valid
            if len(ws) == 1:
                return D[ws[0]]
            return D[ws[0]] & D[ws[1]]

        idx = np.unique(np.concatenate([np.arange(0, 4096), rng.integers(0, N, 16384)])).astype(np.int64)
        verdict, coef, rk = None, None, None
        for rnd in range(12):
            cols = np.stack([bit_of(mon_words(m), idx) for m in mons], axis=1)
            fs = bit_of(f, idx)
            ok, coef, rk = solve(cols, fs)
            if not ok:
                verdict = "DEAD (inconsistent on %d inputs, rank %d of %d monomials, round %d)" % (len(idx), rk, len(mons), rnd)
                break
            acc = np.zeros(W, dtype=np.uint64)
            for m, cbit in zip(mons, coef):
                if cbit:
                    acc ^= mon_words(m)
            diff = (acc ^ f) & valid
            nbad = int(np.unpackbits(diff.view(np.uint8)).sum())
            print("b%d t%d round %d: sample %d rank %d of %d; solution wrong on %d of %d inputs" % (
                b, t, rnd, len(idx), rk, len(mons), nbad, N), flush=True)
            if nbad == 0:
                verdict = "ALIVE (exact on all %d inputs)" % N
                break
            bad = np.flatnonzero(np.unpackbits(diff.view(np.uint8), bitorder="little")[:N])
            add = bad if len(bad) <= 4096 else rng.choice(bad, 4096, replace=False)
            idx = np.unique(np.concatenate([idx, add]))
        terms_g = [m[1] for m, cb in zip(mons, coef) if cb and m[0] == "H"] if coef is not None else []
        terms_h = [m[1] for m, cb in zip(mons, coef) if cb and m[0] == "T"] if coef is not None else []
        rec = {"b": b, "t": t, "gate_t": " ".join([ops[t][0]] + ops[t][1]) + (" if b%d" % ops[t][2] if ops[t][2] else ""),
               "sel_shots_t": selc, "verdict": verdict, "rank": rk, "n_mons": len(mons),
               "g_terms": [list(x) for x in terms_g] if verdict and verdict.startswith("ALIVE") else None,
               "h_terms": [list(x) for x in terms_h] if verdict and verdict.startswith("ALIVE") else None}
        print("RESULT b%d at op %d (%s): %s; g %d terms, h %d terms" % (b, t, rec["gate_t"], verdict, len(terms_g),
                                                                     len(terms_h)), flush=True)
        out.append(rec)
    json.dump({"coins": "one" if COINS_ONE else SEED, "x_restored": bool(ok_x), "scratch_clean": bool(clean),
               "results": out}, open(outp, "w"), indent=1)
    print("time", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
