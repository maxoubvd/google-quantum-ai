#!/usr/bin/env python
"""chem_wr7_roundg (WR7) = r50_chem_lbl_round.py + the last refit block checked (chem_wr7_mkround.py).
r50_chem_lbl_round = r48_chem_q7c_round.py + ONE option, --wires (R50-CHEM-LBL, construction 9 / 9b): bit-0 delta
wires on address bits with the address of a row decided at freeze time (DynChain of r50_chem_lbl_chain.py). Three
hooks, each marked "r50_chem_lbl": the objects after the init / resume block (a resumed run replays its frozen wire
words), the freeze (a wire word whose parity its row's earlier words already decide takes the nearest word of that
parity), and the end (rows permuted inside their groups to the addresses chosen; rounded.npz is in address order,
`order` composed). Needs --no-relations --no-pairs --carriers none (R7g's 7-bit recipe). Everything else unchanged.

r48_chem_q7c_round: copy of chem_own_round7a.py (2026-11-12 06:20Z, workstream chem_own2) with additions:
  --log-every N    log the error every N iterations inside a stage (convergence curves of the refits; 0 = off)
Everything below is chem_own_round7a.py unchanged unless a comment says "chem_own2".

chem_own_round7a: chem_own_round7.py plus the ROUNDING-AWARE DESCENT ON EVERY ANGLE (panel item C2):
  --aware EPS      every free angle is evaluated at +EPS*d and -EPS*d (d uniform in half a rounding step, fresh at
                   every recentering) and the objective is the mean of the two errors: to third order
                   error + EPS^2/2 d'Hd, a penalty on the curvature of the error along the rounding still to come;
  --aware-pre N    a first stage of N iterations with nothing frozen (the descent proper);
  --aware-until F  the term stays in the stage refits while the frozen share is below F (default 1.0);
  --aware-seed S   seed of the noise. With --aware 0 this is chem_own_round7.py.
Made by make_round7a.py (scratchpad of the session of 2026-11-11) from chem_own_round7.py as of 22:40Z.

chem_own_round7: chem_own_round.py with two additions for the panel design CHEM-1 (research/panel), step 1:
  --bits N     word width (default 8): scale 2^(N-1) steps per half turn, signed range -2^(N-2) .. 2^(N-2)-1,
               carrier = bit N-1 xor bit N-2;
  twin pairs   if INIT.npz holds `pairs` (flat row indices master, twin; made by tools/chem_own_pairs.py): the
               two rows have the same vector and the weights p w and (1 - p) w. Their free angles stay tied,
               and an entry of the pair is frozen JOINTLY: on an uncoded column the floor / ceil combination
               whose weighted mean p q1 + (1 - p) q2 is nearest to the angle (first-order spacing of a third
               of a step at p = 1/3: "antithetic split copies", the panel's mechanism); on a coded column each
               twin takes its own nearest coded word, both at the same moment. A pair whose gauge moves
               differ after a stage is untied from then on (counted in the log). --no-pairs ignores them.
               Since 19:20Z (fix round 1): the joint choice also on coded columns (--pair-coded), with the
               current weights of the two rows (--pair-weights) and a cost for the step between the two
               words (--pair-step); --pair-mode same is the control (the unsplit vector); --pairs-first
               freezes the paired vectors alone in stage 0 (probe of the mechanism on its own vectors).
The log also gives, per stage, the share of the frozen impact that belongs to the 1 / 2 / 4 / 6 % strongest
vectors, and at the end the damage attributed to them.
Everything else is chem_own_round.py as of 2026-11-11 17:30Z (kept separate so that running long roundings,
which reload that file at every piece, are not touched).

chem_own_round: OUR progressive rounding of a continuous DF-THC table to 8-bit Givens words, with codes on the
low bits, the rows dealt over the address blocks, and every word kept in the signed range without a wall.

    chem_own_round.py INIT.npz OUT_DIR --codes "0:1:3-14,1:1:3-14" [--xcodes "2:2:8-14"] [--schedule full|long]
                      [--iters 3000] [--final 6000] [--max-minutes 25] [--no-deal] [--order ORDER.npy]
                      [--stop-at 0.02] [--no-relations] [--carriers 3,5,6,9,10|none] [--direct]
    chem_own_round.py --selftest          (CPU only: codes, residues, gauge moves, dealing)

INIT.npz: U (R, B, 76), W (R, B, C) of a continuous fit (tools/chem_factor_fit.py). OUT_DIR gets log.txt,
state.npz after every stage (rerun the same command to resume), and at the end rounded.npz:
    words (R*B, 75) in 0..255 on the log-depth tree, W (R, B, C), U (decoded vectors), order (row order
    applied: order[address - 76] = row index of INIT), codes, xcodes, trace.
--max-minutes M stops after the stage that passes M minutes (GPU jobs here are short pieces under a lock).

Codes ('bit:degree:first-last' over the fifteen five-word batches, comma separated):
  --codes   per column: that bit of every word of the column is a Reed-Muller codeword of that degree in the
            four low address bits on every 16-address block (degree 1 = affine: the public bit-0 code);
  --xcodes  the same for the XOR of the five words of a batch at that bit (marculera's XOR-of-five codes).
  The coded bits of a column must be contiguous from bit 0.

Method and credits.
  * Staged freezing with a refit of everything still free after each stage: ahmedberatozer's progressive rounding
    as used in marculera's public pipeline (quantize_rm.py --open: "open codewords", the cheapest codeword that
    agrees with the frozen entries is re-chosen at every stage); slack entries for the XOR codes: the public
    rmcode.slack idea in the form of our chem_lens tool. Nothing of the public code is imported or copied: the
    Reed-Muller codewords, the residue rounding, the codeword choice, the checks and the gauge are written here.
  * Engine: float32 L-BFGS on the float64-difference objective and row-scaled angles (our chem_factor_fit
    dlbfgs / chem_lens_codes), impact of an entry = weight norm^2 x mass entering the rotation x squared offset.
  * Signed range without a wall: the result of our chemistry-lens workstream (research/ws_chem_lens.md section 6,
    tools/chem_lens_gauge.py): a free angle may cross +-pi/2 when every entry the crossing negates is negatable;
    the rows are re-gauged exactly once at the end of every stage. New here: the negatable rule for PER-COLUMN
    codes above bit 0 (a frozen word w may be negated iff -w = w modulo 2^(coded bits), i.e. w = 0 modulo
    2^(coded bits - 1)), and the row dealing over the TRUE address blocks (a group starts inside a block).
  * Carriers (--carriers, default the pair monomials 3,5,6,9,10): bit7 xor bit6 of word position j is the j-th
    monomial of the low address bits on every row and in every batch, as in the public words (the "word forms"
    of marculera's table, read off its words with our chem_loaders_core.analyze). Here the pattern is a gauge
    fixed BEFORE the rounding: the entry's angle is theta + pi where the monomial is 1 and theta stays in the
    signed range, so the pattern costs no accuracy and never meets the codes (it only touches bit 7).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

N = 76
NB_BATCH = 15
HALF = 64                 # signed range of a word: -64 .. 63 (theta in [-pi/2 - pi/256, pi/2 - pi/256))
FULL = tuple(round(1.0 - f, 4) for f in (0.9995, 0.999, 0.9975, 0.995, 0.99, 0.98, 0.965, 0.95, 0.925, 0.9, 0.86,
                                         0.82, 0.77, 0.72, 0.66, 0.6, 0.54, 0.48, 0.42, 0.36, 0.3, 0.25, 0.2, 0.16,
                                         0.12, 0.09, 0.065, 0.045, 0.03, 0.02, 0.012, 0.006, 0.0))
LONG = tuple(round(1.0 - f, 5) for f in (0.9995, 0.999, 0.998, 0.997, 0.9955, 0.994, 0.992, 0.99, 0.9875, 0.985,
                                         0.98, 0.975, 0.97, 0.965, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.9, 0.885,
                                         0.87, 0.855, 0.84, 0.82, 0.8, 0.78, 0.76, 0.74, 0.72, 0.69, 0.66, 0.63, 0.6,
                                         0.57, 0.54, 0.51, 0.48, 0.45, 0.42, 0.39, 0.36, 0.33, 0.3, 0.27, 0.24, 0.21,
                                         0.18, 0.15, 0.125, 0.1, 0.08, 0.06, 0.045, 0.03, 0.02, 0.012, 0.006, 0.0))


# ------------------------------------------------------------------------------------------------ codes
def rm_code(d: int) -> np.ndarray:
    """All codewords of the Reed-Muller code RM(d, 4): (2^k, 16) array of 0/1, position = low address."""
    pos = np.arange(16)
    b = [(pos >> i) & 1 for i in range(4)]
    mons = [np.ones(16, np.int64)]
    if d >= 1:
        mons += b
    if d >= 2:
        mons += [b[i] & b[j] for i in range(4) for j in range(i + 1, 4)]
    if d >= 3:
        mons += [b[i] & b[j] & b[k] for i in range(4) for j in range(i + 1, 4) for k in range(j + 1, 4)]
    G = np.array(mons, np.int64)
    k = len(G)
    coef = (np.arange(1 << k)[:, None] >> np.arange(k)[None, :]) & 1
    return (coef @ G) % 2


def parse_codes(spec: str, ncols: int = N - 1, bits: int = 8) -> np.ndarray:
    """deg[column, bit] = degree of the per-column code on that bit (0 = none; -1 = RM(0,4), the word bit constant
    on every 16-row block: rm_code(-1) gives the two constant codewords; chem_own2, 2026-11-14)."""
    deg = np.zeros((ncols, bits), np.int64)
    for item in filter(None, (spec or "").split(",")):
        bit, d, span = item.split(":")
        lo, hi = (int(x) for x in span.split("-")) if "-" in span else (int(span),) * 2
        deg[5 * lo:5 * hi + 5, int(bit)] = int(d)
    return deg


def parse_xcodes(spec: str, nbatch: int = NB_BATCH, bits: int = 8) -> np.ndarray:
    """xdeg[batch, bit] = degree of the code on the XOR of the batch's five words at that bit (0 = none)."""
    xdeg = np.zeros((nbatch, bits), np.int64)
    for item in filter(None, (spec or "").split(",")):
        bit, d, span = item.split(":")
        lo, hi = (int(x) for x in span.split("-")) if "-" in span else (int(span),) * 2
        xdeg[lo:hi + 1, int(bit)] = int(d)
    return xdeg


def residue_round(t, r, m):
    """Nearest integer to t that is congruent to r modulo m (arrays; m per column or scalar)."""
    return r + m * np.round((t - r) / m)


def into_range(v, m):
    """The word of the same residue modulo m inside the signed range -64..63 that is nearest to v."""
    hi, lo = HALF - 1, -HALF
    v = np.where(v > hi, v - m * np.ceil((v - hi) / m), v)
    return np.where(v < lo, v + m * np.ceil((lo - v) / m), v)


def check_vec(p: int, d: int) -> int:
    """Column p of the parity-check matrix of RM(d, 4) as a bit mask (monomials of degree <= 3 - d)."""
    b = [(p >> i) & 1 for i in range(4)]
    mons = [1]
    if d <= 2:
        mons += b
    if d <= 1:
        mons += [b[i] & b[j] for i in range(4) for j in range(i + 1, 4)]
    return sum(v << i for i, v in enumerate(mons))


def xor_slack(c, blk, pos, xdeg):
    """Entries frozen last: per block and XOR-coded batch, the weakest word of each of the weakest rows whose
    parity checks are independent (11 rows for an affine code, 5 for degree 2, 1 for degree 3)."""
    res = np.zeros(c.shape, bool)
    for k in np.flatnonzero(xdeg.any(axis=1)):
        d = int(xdeg[k][xdeg[k] > 0].min())
        need = {1: 11, 2: 5, 3: 1}[d]
        cols = np.arange(5 * k, 5 * k + 5)
        cmin, jmin = c[:, cols].min(axis=1), c[:, cols].argmin(axis=1)
        for b in np.unique(blk):
            idx = np.flatnonzero(blk == b)
            if len(idx) <= 4:
                continue
            piv = {}
            for i in idx[np.argsort(cmin[idx])]:
                v = check_vec(int(pos[i]), d)
                while v and v.bit_length() in piv:
                    v ^= piv[v.bit_length()]
                if v:
                    piv[v.bit_length()] = v
                    res[i, cols[jmin[i]]] = True
                    if len(piv) == need:
                        break
    return res


_SUB5 = ((np.arange(32)[:, None] >> np.arange(5)[None, :]) & 1).astype(np.float64)   # subset of the 5 words
_PAR5 = (_SUB5.sum(axis=1).astype(np.int64)) & 1
_XOR5 = np.arange(32)[:, None] ^ np.arange(32)[None, :]                             # [state after, F] -> before


def _row_dp(m, delta, target):
    """Rows i of one block (or half block): flip a subset F_i of the five words with |F_i| = m[cw, i] (mod 2) at cost
    sum_{j in F_i} delta[i, j], the flips of word j summing to target bit j (mod 2).  Exact DP (state = the five
    flip parities).  Returns (cost per codeword (ncw,), subsets (ncw, nrows))."""
    ncw, nr = m.shape
    big = 1e300
    cF = delta @ _SUB5.T                                                           # (nr, 32)
    V = np.full((ncw, 32), big)
    V[:, 0] = 0.0
    choice = np.zeros((nr, ncw, 32), np.int64)
    for i in range(nr):
        cand = V[:, _XOR5] + cF[i][None, None, :]
        cand = np.where(_PAR5[None, None, :] != m[:, i][:, None, None], big, cand)
        choice[i] = cand.argmin(axis=2)
        V = np.take_along_axis(cand, choice[i][:, :, None], axis=2)[:, :, 0]
    pick = np.zeros((ncw, nr), np.int64)
    st = np.full(ncw, int(target), np.int64)
    ac = np.arange(ncw)
    for i in range(nr - 1, -1, -1):
        F = choice[i][ac, st]
        pick[:, i] = F
        st = st ^ F
    return V[:, int(target)], pick


def joint_parity_xor(t0, t1, d, xd):
    """Five words of one XOR-coded batch whose bit carries BOTH a per-column block parity (d = 3: even on the 16 rows;
    d = 4: even on both 8-row halves) and the XOR-of-five code RM(xd, 4) (chem_own2, Nov 14: the L2 rung on a
    low-degree-delta table).  t0, t1 (blocks, 16, 5): costs of the bit values.  Exact per block: every codeword of
    the XOR code (pruned by a lower bound), DP over the rows for the flips.  Returns the bits (blocks, 16, 5)."""
    beta = (t1 < t0).astype(np.int64)
    delta = np.abs(t1 - t0)
    cws = rm_code(int(xd))
    halves = [np.arange(16)] if d == 3 else [np.arange(8), np.arange(8, 16)]
    out = beta.copy()
    for b in range(beta.shape[0]):
        x = np.bitwise_xor.reduce(beta[b], axis=1)                                 # (16,) XOR of the five bits
        m = cws ^ x[None, :]                                                       # (ncw, 16) row corrections
        tgt = [int(sum(((int(beta[b][h, j].sum()) & 1) << j) for j in range(5))) for h in halves]
        # lower bound per codeword: every row's cheapest subset of the right parity, the parities ignored
        cF = delta[b] @ _SUB5.T
        ce = np.where(_PAR5[None, :] == 0, cF, np.inf).min(axis=1)
        co = np.where(_PAR5[None, :] == 1, cF, np.inf).min(axis=1)
        lb = np.where(m == 1, co[None, :], ce[None, :]).sum(axis=1)
        order = np.argsort(lb, kind="stable")
        best, bestF = np.inf, None
        for c0 in range(0, len(order), 64):
            idx = order[c0:c0 + 64]
            if lb[idx[0]] >= best:
                break
            tot = np.zeros(len(idx))
            Fs = np.zeros((len(idx), 16), np.int64)
            for h, tg in zip(halves, tgt):
                cst, pk = _row_dp(m[idx][:, h], delta[b][h], tg)
                tot += cst
                Fs[:, h] = pk
            j = int(np.argmin(tot))
            if tot[j] < best:
                best, bestF = float(tot[j]), Fs[j]
        flips = ((bestF[:, None] >> np.arange(5)[None, :]) & 1)
        out[b] = beta[b] ^ flips
    return out


def _shared_checks(d: int) -> np.ndarray:
    """Parity checks of the SHARED codes (chem_own2, Nov 14): the syndrome a word must share with the words of the
    same slot in every batch of the code, on every block (L2S's condition: the high ANF coefficients the same for
    all the batches, a free pattern per block, slot and bit).  13: the 16-row parity (code 3 relaxed); 14: the two
    8-row parities (code 4 relaxed); 12: the RM(2,4) coset = the 5 checks of RM(1,4) (code 2 relaxed)."""
    p = np.arange(16)
    if d == 13:
        return np.ones((1, 16), np.int64)
    if d == 14:
        return np.stack([(p < 8).astype(np.int64), (p >= 8).astype(np.int64)])
    if d == 12:
        return np.stack([np.ones(16, np.int64)] + [(p >> i) & 1 for i in range(4)])
    raise ValueError(d)


def shared_coset_bits(t0, t1, cols, d, xor0=False):
    """Cheapest bits of the columns `cols` (one code-d group per slot = cols % 5) on every block, every column's
    syndrome equal to its slot group's common one, chosen to minimise the group's total (exact: a DP over the 16
    positions per word for every syndrome, then the best common syndrome per block and slot).  With xor0 the five
    slots' syndromes must XOR to 0 on every block (exact DP over the slots), so that the L2 rung's XOR-of-five codes
    can still be put on these bits later (route C).  t0, t1 (blocks, 16, ncols).  Returns the bits (blocks, 16,
    ncols)."""
    H = _shared_checks(d)
    r = H.shape[0]
    S = 1 << r
    hm = (H * (1 << np.arange(r))[:, None]).sum(axis=0)                            # syndrome mask of each position
    nb, _, nc = t0.shape
    V = np.full((nb, nc, S), np.inf)
    V[:, :, 0] = 0.0
    choice = np.zeros((16, nb, nc, S), np.int8)
    for i in range(16):
        perm = np.arange(S) ^ int(hm[i])
        a0 = V + t0[:, i, :][:, :, None]                                           # bit 0 here: syndrome unchanged
        a1 = V[:, :, perm] + t1[:, i, :][:, :, None]                               # bit 1 here: came from s ^ h_i
        choice[i] = (a1 < a0)
        V = np.minimum(a0, a1)
    target = np.zeros((nb, nc), np.int64)
    groups = [np.flatnonzero(np.asarray(cols) % 5 == j) for j in range(5)]
    if not xor0:
        for g in groups:
            if len(g):
                target[:, g] = V[:, g, :].sum(axis=1).argmin(axis=1)[:, None]
    else:
        tot = np.stack([V[:, g, :].sum(axis=1) if len(g) else np.where(np.arange(S)[None, :] == 0, 0.0, np.inf)
                        * np.ones((nb, 1)) for g in groups], axis=1)               # (nb, 5, S)
        W_ = np.full((nb, S), np.inf)
        W_[:, 0] = 0.0
        back = np.zeros((5, nb, S), np.int64)
        idx = np.arange(S)[:, None] ^ np.arange(S)[None, :]                        # [state after, p] -> before
        for j in range(5):
            cand = W_[:, idx] + tot[:, j, :][:, None, :]                          # (nb, state after, p)
            back[j] = cand.argmin(axis=2)
            W_ = cand.min(axis=2)
        st = np.zeros(nb, np.int64)
        for j in range(4, -1, -1):
            p_ = back[j][np.arange(nb), st]
            if len(groups[j]):
                target[:, groups[j]] = p_[:, None]
            st = st ^ p_
        assert not st.any()
    out = np.zeros((nb, 16, nc), np.int64)
    s = target.copy()
    ib, ic = np.meshgrid(np.arange(nb), np.arange(nc), indexing="ij")
    for i in range(15, -1, -1):
        bit = choice[i][ib, ic, s].astype(np.int64)
        out[:, i, :] = bit
        s = np.where(bit == 1, s ^ int(hm[i]), s)
    assert not s.any()
    return out


_RM2W = None


def joint_rm2_affine(t0, t1):
    """Five words of one XOR-coded batch whose bit is per-column RM(2,4) (code 2) AND whose XOR must be affine
    (RM(1,4), the L2 rung's bit-1 code) (chem_own2, Nov 14).  Exact per block: the cheapest RM(2,4) word of every
    quadratic class for each word, then a DP over the five words (state = the XOR of the quadratic parts, which
    must end at 0).  t0, t1 (blocks, 16, 5).  Returns the bits (blocks, 16, 5)."""
    global _RM2W
    if _RM2W is None:
        _RM2W = rm_code(2)                     # index bits: 1, x0..x3 (affine), then the six products (quadratic)
    W = _RM2W.astype(np.float64)
    nb = t0.shape[0]
    out = np.zeros(t0.shape, np.int64)
    for b in range(nb):
        cost = t0[b].sum(axis=0)[None, :] + W @ (t1[b] - t0[b])                   # (2048, 5)
        cls = cost.reshape(64, 32, 5)
        ba, bc = cls.argmin(axis=1), cls.min(axis=1)                               # (64, 5)
        V = np.full(64, np.inf)
        V[0] = 0.0
        back = np.zeros((5, 64), np.int64)
        for j in range(5):
            cand = V[np.arange(64)[:, None] ^ np.arange(64)[None, :]] + bc[:, j][None, :]
            back[j] = cand.argmin(axis=1)
            V = cand.min(axis=1)
        st = 0
        for j in range(4, -1, -1):
            q = int(back[j][st])
            out[b, :, j] = _RM2W[32 * q + int(ba[q, j])]
            st ^= q
    return out


_PROP = (2, 3, 4, 12, 13, 14)        # the low-degree-delta codes (on bits >= 1)


def prop_columns(deg):
    """Columns whose bits >= 1 carry a low-degree-delta code (chem_own2)."""
    return np.array([any(int(deg[c_, b]) in _PROP for b in range(1, deg.shape[1])) for c_ in range(deg.shape[0])])


def column_slack(c, blk, pos, deg):
    """--col-slack (chem_own2, Nov 14): entries frozen last for the low-degree-delta codes.  Per full block and
    property-coded column, the weakest entries that can carry its checks: one for a 16-row parity (codes 3 / 13), the
    weakest of each 8-row half (4 / 14), the weakest five rows with independent RM(1,4) checks (2 / 12).  They absorb
    every parity fix of bits >= 1 until they are frozen themselves, at the end."""
    res = np.zeros(c.shape, bool)
    for col in np.flatnonzero(prop_columns(deg)):
        ds = {int(deg[col, b]) for b in range(1, deg.shape[1])}
        kind = 2 if ds & {2, 12} else (4 if ds & {4, 14} else 3)
        for b in np.unique(blk):
            idx = np.flatnonzero(blk == b)
            if len(idx) < 16:
                continue
            order = idx[np.argsort(c[idx, col], kind="stable")]
            if kind == 3:
                res[order[0], col] = True
            elif kind == 4:
                res[[i for i in order if pos[i] < 8][0], col] = True
                res[[i for i in order if pos[i] >= 8][0], col] = True
            else:
                piv = {}
                for i in order:
                    v = check_vec(int(pos[i]), 2)
                    while v and v.bit_length() in piv:
                        v ^= piv[v.bit_length()]
                    if v:
                        piv[v.bit_length()] = v
                        res[i, col] = True
                        if len(piv) == 5:
                            break
    return res


def residues(t, c, blk, pos, deg, xdeg, frozen, value, reserved=None):
    """Low-bit residue of every entry on the cheapest codewords that agree with the frozen entries, bit by bit
    from bit 0. t (rows, 75) in steps, c the weights, blk / pos the address block and low address of each row,
    value the frozen words (signed integers). Returns (r, modulus per column)."""
    rows, ncols = t.shape
    nblocks = int(blk.max()) + 1
    r = np.zeros((rows, ncols), np.int64)
    nbits = np.zeros(ncols, np.int64)
    for bit in range(deg.shape[1]):
        m = 1 << bit
        pc = np.flatnonzero(deg[:, bit])
        xb = np.flatnonzero(xdeg[:, bit])
        # chem_own2, Nov 14: an XOR-coded batch whose five columns ALSO carry a per-column block code on this bit is
        # solved jointly (codes 3 / 4: joint_parity_xor); code 2 under an RM(2,4) XOR code needs nothing (the XOR of
        # five RM(2,4) words is one); other mixtures are refused
        jt, xkeep = [], []
        for k in xb:
            dk = deg[5 * k:5 * k + 5, bit]
            if not np.any(dk != 0):
                xkeep.append(int(k))
                continue
            assert np.all(dk == dk[0]), f"batch {k}, bit {bit}: per-column codes on only some of an XOR batch's words"
            if int(dk[0]) in (3, 4) or (int(dk[0]) == 2 and int(xdeg[k, bit]) == 1):
                jt.append((int(k), int(dk[0])))
            elif not (int(dk[0]) == 2 and int(xdeg[k, bit]) >= 2):
                raise NotImplementedError(f"per-column code {int(dk[0])} with XOR code {int(xdeg[k, bit])} (batch {k}, bit {bit})")
        xb = np.array(xkeep, np.int64)
        jcols = np.array([5 * k + j for k, _ in jt for j in range(5)], np.int64)
        act = np.concatenate([pc] + [np.arange(5 * k, 5 * k + 5) for k in xb]).astype(np.int64)
        if not len(act):
            continue
        assert np.all(nbits[act] == bit), "coded bits of a column must be contiguous from bit 0"
        cost = [c[:, act] * (into_range(residue_round(t[:, act], r[:, act] + v * m, 2 * m), 2 * m) - t[:, act]) ** 2
                for v in (0, 1)]
        if reserved is not None:
            rs = reserved[:, act].copy()
            # only the words of XOR-coded batches absorb relations, and (chem_own2, --col-slack) the slack entries of
            # columns whose bit `bit` >= 1 carries a low-degree-delta code
            rs[:, :len(pc)] &= (np.isin(deg[act[:len(pc)], bit], _PROP) & (bit >= 1))[None, :]
            if rs.any():
                cost = [np.where(rs, 1e-9 * x, x) for x in cost]
        have = (value[:, act] >> bit) & 1
        big = 1e4 * max(float(x.max()) for x in cost) + 1e-300
        cost = [np.where(frozen[:, act], np.where(have == v, 0.0, big), cost[v]) for v in (0, 1)]
        table = [np.zeros((nblocks, 16, len(act))) for _ in (0, 1)]
        for v in (0, 1):
            table[v][blk, pos] = cost[v]
        rbit = np.zeros((rows, len(act)), np.int64)
        for k, dd in jt:                                        # joint: block parity + XOR-of-five (chem_own2)
            sidx = np.flatnonzero(np.isin(act[:len(pc)], np.arange(5 * k, 5 * k + 5)))
            assert len(sidx) == 5
            if dd == 2:
                rbit[:, sidx] = joint_rm2_affine(table[0][:, :, sidx], table[1][:, :, sidx])[blk, pos]
            else:
                rbit[:, sidx] = joint_parity_xor(table[0][:, :, sidx], table[1][:, :, sidx], dd, int(xdeg[k, bit]))[blk, pos]
        for d in (np.unique(deg[pc, bit]) if len(pc) else []):
            sidx = np.flatnonzero((deg[act[:len(pc)], bit] == d) & ~np.isin(act[:len(pc)], jcols))
            if not len(sidx):
                continue
            if int(d) in (12, 13, 14):          # shared codes: one syndrome per block and slot, common to the batches
                rbit[:, sidx] = shared_coset_bits(table[0][:, :, sidx], table[1][:, :, sidx],
                                                  act[:len(pc)][sidx], int(d), xor0=bit in (1, 2))[blk, pos]
                continue
            if int(d) == 3:
                # RM(3,4) = the 2^15 even-weight words (one parity per block): the cheapest codeword is every entry's
                # cheaper bit, with the cheapest single flip where that leaves the block's parity odd (exact; frozen
                # entries carry the prohibitive cost and are flipped only if nothing else can); chem_own2, Nov 14
                t0, t1 = table[0][:, :, sidx], table[1][:, :, sidx]           # (blocks, 16, columns)
                beta = (t1 < t0).astype(np.int64)
                delta = np.abs(t1 - t0)
                bi, ci = np.nonzero(beta.sum(1) % 2)
                jmin = delta.argmin(1)
                beta[bi, jmin[bi, ci], ci] ^= 1
                rbit[:, sidx] = beta[blk, pos]
                continue
            if int(d) == 4:
                # code 4 (chem_own2, Nov 14; NOT RM(4,4)): the ANF coefficients of x0x1x2x3 AND x0x1x2 both zero =
                # even weight on both 8-row halves of the block (low addresses 0-7 and 8-15); the exact rule of code 3
                # applied per half
                t0, t1 = table[0][:, :, sidx], table[1][:, :, sidx]
                nb_, nc_ = t0.shape[0], t0.shape[2]
                t0, t1 = t0.reshape(nb_, 2, 8, nc_), t1.reshape(nb_, 2, 8, nc_)
                beta = (t1 < t0).astype(np.int64)
                delta = np.abs(t1 - t0)
                bi, hi_, ci = np.nonzero(beta.sum(2) % 2)
                jmin = delta.argmin(2)
                beta[bi, hi_, jmin[bi, hi_, ci], ci] ^= 1
                rbit[:, sidx] = beta.reshape(nb_, 16, nc_)[blk, pos]
                continue
            words = rm_code(int(d))
            gain = (table[1] - table[0])[:, :, sidx].transpose(0, 2, 1) @ words.T.astype(np.float64)
            pick = words[gain.argmin(axis=2)]                                   # (blocks, columns, 16)
            rbit[:, sidx] = pick[blk, :, pos]
        for q_, k in enumerate(xb):
            sidx = len(pc) + 5 * q_ + np.arange(5)
            t0, t1 = table[0][:, :, sidx], table[1][:, :, sidx]                 # (blocks, 16, 5)
            beta = (t1 < t0).astype(np.int64)
            delta = np.abs(t1 - t0)
            x = beta.sum(-1) % 2
            jmin, dmin = delta.argmin(-1), delta.min(-1)
            words = rm_code(int(xdeg[k, bit]))
            costcw = ((words[None, :, :] != x[:, None, :]) * dmin[:, None, :]).sum(-1)
            cw = words[costcw.argmin(1)]
            flip = (cw != x).astype(np.int64)
            bi, pi = np.nonzero(flip)
            beta[bi, pi, jmin[bi, pi]] ^= 1
            rbit[:, sidx] = beta[blk, pos]
        r[:, act] += rbit * m
        nbits[act] += 1
    return r, 1 << nbits


def check_codes(words, addr, deg, xdeg):
    """Raises if a per-column or XOR-of-five code fails on any address block. words: integers (any range)."""
    blk, pos = addr // 16, addr % 16
    w = np.asarray(words).astype(np.int64)
    for d_ in (12, 13, 14):                  # shared codes: the slot group's syndromes agree on every full block
        H = _shared_checks(d_)
        for bit in range(deg.shape[1]):
            cs = np.flatnonzero(deg[:, bit] == d_)
            for j in range(5):
                g = cs[cs % 5 == j]
                if not len(g):
                    continue
                for b in np.unique(blk):
                    rows = np.flatnonzero(blk == b)
                    if len(rows) < 16:
                        continue
                    v = (w[rows][:, g] >> bit) & 1                                  # (16, len(g)), rows in address order
                    syn = (H[:, pos[rows]] @ v) % 2                                 # (r, len(g))
                    if np.any(syn != syn[:, :1]):
                        raise AssertionError(f"shared code {d_} fails: bit {bit}, slot {j}, block {b}")
            if bit in (1, 2) and len(cs) and np.all(np.isin(np.arange(5), cs % 5)):    # xor0: slot syndromes XOR to 0
                for b in np.unique(blk):
                    rows = np.flatnonzero(blk == b)
                    if len(rows) < 16:
                        continue
                    acc = np.zeros(H.shape[0], np.int64)
                    for j in range(5):
                        c0 = cs[cs % 5 == j][0]
                        acc ^= (H[:, pos[rows]] @ ((w[rows, c0] >> bit) & 1)) % 2
                    if acc.any():
                        raise AssertionError(f"shared code {d_}: slot syndromes do not XOR to 0, bit {bit}, block {b}")
    for col, bit in zip(*np.nonzero(deg)):
        if int(deg[col, bit]) in (12, 13, 14):
            continue
        v = (w[:, col] >> bit) & 1
        if int(deg[col, bit]) == 3:          # even weight on every full block (a partial block always extends)
            for b in np.unique(blk):
                rows = np.flatnonzero(blk == b)
                if len(rows) == 16 and int(v[rows].sum()) % 2:
                    raise AssertionError(f"code fails: bit {bit} of column {col}, block {b}")
            continue
        if int(deg[col, bit]) == 4:          # even weight on both 8-row halves of every full block
            for b in np.unique(blk):
                rows = np.flatnonzero(blk == b)
                if len(rows) == 16:
                    h = pos[rows] >= 8
                    if int(v[rows][~h].sum()) % 2 or int(v[rows][h].sum()) % 2:
                        raise AssertionError(f"code fails (halves): bit {bit} of column {col}, block {b}")
            continue
        cw = rm_code(int(deg[col, bit]))
        for b in np.unique(blk):
            rows = np.flatnonzero(blk == b)
            if not np.any(np.all(cw[:, pos[rows]] == v[rows], axis=1)):
                raise AssertionError(f"code fails: bit {bit} of column {col}, block {b}")
    for k, bit in zip(*np.nonzero(xdeg)):
        cw = rm_code(int(xdeg[k, bit]))
        x = np.bitwise_xor.reduce((w[:, 5 * k:5 * k + 5] >> bit) & 1, axis=1)
        for b in np.unique(blk):
            rows = np.flatnonzero(blk == b)
            if not np.any(np.all(cw[:, pos[rows]] == x[rows], axis=1)):
                raise AssertionError(f"XOR code fails: bit {bit} of batch {k}, block {b}")


# ------------------------------------------------------------------------------------------- the target
DATA = os.environ.get("QCC_DATA", os.path.expanduser("~/qcc/data"))


class Target:
    """The two-body target in the packed space of orbital pairs (i >= j), from the hosts' reference loader only:
        error_2 = || sqrt(m m^T) * (A^T A - T) ||_F,   A[(r, c), p] = sum_b w[r, b, c] u[r, b, i_p] u[r, b, j_p],
        T = tei + (S (x) I + I (x) S) / 2,   m = 1 on diagonal pairs and 2 elsewhere.
    (Our formulation from tools/chem_factor_lib.py, where it is checked against the verifier's error_2; repeated
    here so that this tool loads nothing but numpy, JAX and the hosts' verifier package.)"""

    def __init__(self):
        for p in (os.path.join(DATA, "verifiers"), os.path.join(DATA, "kickmix")):
            if p not in sys.path:
                sys.path.insert(0, p)
        os.environ.setdefault("DFTHC_DATA", DATA)
        from dfthcsos.scoreinfra.reference import check_fcidump, load_reference
        from dfthcsos.scoreinfra.dfthc_math import symmetric_from_lower
        ref = load_reference(check_fcidump())
        with np.load(os.path.join(DATA, "verifiers", "dfthcsos", "FeMoCo76o113e.fixed_shifts.npz")) as z:
            S = np.asarray(symmetric_from_lower(np.array(z["Hsym_lower"]), N), np.float64)
        tei = np.asarray(ref["tei"], np.float64)
        i, j = np.tril_indices(N)
        T = tei[i, j][:, i, j].copy()
        eye = (i == j).astype(np.float64)
        T += 0.5 * (np.outer(S[i, j], eye) + np.outer(eye, S[i, j]))
        m = np.where(i == j, 1.0, 2.0)
        self.i, self.j, self.T, self.M = i, j, T, np.outer(m, m)

    def error2_np(self, u, w):
        R, B, C = w.shape
        u = u.reshape(R, B, N)
        u = u / np.linalg.norm(u, axis=2, keepdims=True)
        A = np.einsum("rbc,rbp->rcp", w, u[:, :, self.i] * u[:, :, self.j]).reshape(R * C, -1)
        H = A.T @ A - self.T
        return float(np.sqrt(np.sum(self.M * H * H)))


# ---------------------------------------------------------------------------------------- tree and gauge
def levels(edges):
    """Runs of consecutive edges touching disjoint orbitals (one vectorised update each)."""
    runs, current, used = [], [], set()
    for t, (p, q) in enumerate(edges):
        if p in used or q in used:
            runs.append(np.array(current))
            current, used = [], set()
        current.append(t)
        used.update((int(p), int(q)))
    return runs + [np.array(current)]


def log_depth_tree(n: int = N) -> np.ndarray:
    edges, level, depth = [], [0], 0
    while len(edges) < n - 1:
        nxt = []
        for label in level:
            nxt.append(label)
            right = label + (1 << depth)
            if right < n:
                edges.append((label, right))
                nxt.append(right)
        level, depth = nxt, depth + 1
    return np.asarray(edges, dtype=np.int64)


def decode_np(theta_over_pi, edges):
    th = np.pi * np.asarray(theta_over_pi, np.float64)
    v = np.zeros((th.shape[0], edges.shape[0] + 1))
    v[:, 0] = 1.0
    for t, (p, c) in enumerate(edges):
        a, b = v[:, p].copy(), v[:, c].copy()
        v[:, p] = np.cos(th[:, t]) * a - np.sin(th[:, t]) * b
        v[:, c] = np.sin(th[:, t]) * a + np.cos(th[:, t]) * b
    return v


def carrier_pattern(addr, masks, ncols=N - 1):
    """K[row, column] = 1 where the word carries a half turn: the monomial masks[column % 5] of the four low
    address bits (bit i of the mask = address bit i) is 1 at the row's address. The word is then w + 128 with w
    in the signed range, i.e. bit7 xor bit6 of the stored word equals that monomial: the same pure function of
    the low address in every batch, which our compiler uses as a ready-made low monomial ("carrier")."""
    K = np.zeros((len(addr), ncols))
    if len(masks):
        assert len(masks) == 5 and all(0 < m < 16 for m in masks)
        low = np.asarray(addr) % 16
        for j, m in enumerate(masks):
            K[:, j::5] = ((low & m) == m)[:, None]
    return K


def angles_and_mass(vectors, edges, K=None):
    """theta' / pi in [-1/2, 1/2] such that the angles theta' + K (K = 0 / 1 half turns, default none: the
    signed gauge) reproduce each row up to a global sign, and the squared amplitude that enters each rotation.
    A half turn on an entry is exact because the sign it gives the merged amplitude is absorbed by the angles
    above it in the tree (they are extracted afterwards)."""
    a = np.array(vectors, dtype=np.float64)
    a /= np.linalg.norm(a, axis=1, keepdims=True)
    theta = np.zeros((len(a), len(edges)))
    mass = np.zeros((len(a), len(edges)))
    for t in range(len(edges) - 1, -1, -1):
        p, q = edges[t]
        s = np.where(a[:, p] > 0, 1.0, np.where(a[:, p] < 0, -1.0, np.where(a[:, q] > 0, -1.0, 1.0)))
        if K is not None:
            s = s * (1.0 - 2.0 * K[:, t])
        theta[:, t] = np.arctan2(s * a[:, q], s * a[:, p])
        h = np.hypot(a[:, p], a[:, q])
        mass[:, t] = h * h
        a[:, p] = s * h
        a[:, q] = 0.0
    theta /= np.pi
    if K is not None:
        theta = np.mod(theta - K + 1.0, 2.0) - 1.0
    return theta, mass


def neighbour_matrix(edges):
    """NB[t, k] = 1 when a crossing (half turn) at edge t negates edge k: the edge into t's parent and t's
    earlier siblings."""
    n = len(edges)
    NB = np.zeros((n, n), np.int32)
    for t, (p, _) in enumerate(edges):
        for k, (pp, qq) in enumerate(edges):
            if (qq == p) or (pp == p and k < t):
                NB[t, k] = 1
    return NB


def negatable(frozen, wv, protected, NB):
    """(negat, flip_ok): which entries may change sign and which crossings are therefore allowed.
    A frozen word -64 is special: its negation is +64 = -64 plus a half turn, so the stored word stays and the
    sign change passes on to that entry's own neighbours; it is negatable exactly when that crossing is allowed
    (least fixed point; the neighbours of an edge have smaller indices, so the iteration ends)."""
    edge = frozen & (wv == -HALF)
    negat = ~edge & ~(protected & ~edge)
    for _ in range(NB.shape[0] + 1):
        flip_ok = ((~negat).astype(np.int32) @ NB.T) == 0
        new = np.where(edge, flip_ok, negat)
        if np.array_equal(new, negat):
            break
        negat = new
    flip_ok = ((~negat).astype(np.int32) @ NB.T) == 0
    return negat, flip_ok


def regauge(th, frozen, flip_ok, negat, edges, lo, hi):
    """Exact re-gauging of every row into [lo, hi) (units of pi). Free entries whose crossing is allowed get a
    half turn and the entries it negates change sign; free entries whose crossing is not allowed are clipped;
    frozen entries only ever change sign, and only where `negat` allows. Returns (angles, crossings, negated
    frozen entries, clipped entries)."""
    out = np.where(frozen, th, np.mod(th + 1.0, 2.0) - 1.0)
    free = ~frozen
    sign = np.ones((th.shape[0], edges.shape[0] + 1), np.int8)
    nflip = nneg = nclip = 0
    for k in range(edges.shape[0] - 1, -1, -1):
        p, q = edges[k]
        neg = sign[:, p] != sign[:, q]
        assert not np.any(neg & ~negat[:, k]), "a protected entry would be negated"
        a = np.where(neg, -out[:, k], out[:, k])
        nneg += int(np.sum(neg & frozen[:, k]))
        oor = (a >= hi) | (a < lo)
        fl = oor & flip_ok[:, k]
        a = np.where(fl, np.where(a >= hi, a - 1.0, a + 1.0), a)
        cl = oor & free[:, k] & ~flip_ok[:, k]
        a = np.where(cl, np.clip(a, lo, hi - 1e-9), a)
        assert not np.any(oor & frozen[:, k] & ~flip_ok[:, k]), "a frozen entry left the signed range"
        out[:, k] = a
        nflip += int(fl.sum())
        nclip += int(cl.sum())
        sign[fl, p] *= -1
        sign[:, q] = 1
    return out, nflip, nneg, nclip


# ------------------------------------------------------------------------------------------- row dealing
def deal_order(strength, R, B, n_one=N):
    """Row order by dealing: inside every group the rows go, strongest first, round-robin over the TRUE
    16-address blocks the group touches (a group starts and ends inside a block), so that every block gets the
    same mix of strong and weak rows. strength: (R, B). Returns order with order[address - n_one] = flat row."""
    order = np.empty(R * B, np.int64)
    for r in range(R):
        a0 = n_one + r * B
        addrs = np.arange(a0, a0 + B)
        blocks = np.unique(addrs // 16)
        slots = {int(b): list(addrs[addrs // 16 == b]) for b in blocks}
        ranks = np.argsort(-strength[r], kind="stable")
        j = 0
        for row in ranks:
            for _ in range(len(blocks)):
                b = int(blocks[j % len(blocks)])
                j += 1
                if slots[b]:
                    order[slots[b].pop(0) - n_one] = r * B + row
                    break
    assert len(set(order.tolist())) == R * B
    return order


# ------------------------------------------------------------------------------------------------ main
def selftest():
    rng = np.random.default_rng(1)
    for d, k in ((0, 1), (1, 5), (2, 11), (3, 15)):
        cw = rm_code(d)
        assert cw.shape == (1 << k, 16) and len({tuple(x) for x in cw.tolist()}) == 1 << k
    # every affine codeword is orthogonal to the parity checks of degree <= 2 (and so on)
    for d in (1, 2):
        H = np.array([[(check_vec(p, d) >> i) & 1 for p in range(16)] for i in range(check_vec(15, d).bit_length())])
        assert not ((rm_code(d) @ H.T) % 2).any()
    edges = log_depth_tree(N)
    NB = neighbour_matrix(edges)
    rows = 300
    addr = N + np.arange(rows)
    blk, pos = addr // 16, addr % 16
    t = rng.normal(0, 30, (rows, N - 1))
    t = np.clip(t, -64.4, 63.4)
    c = rng.random((rows, N - 1)) ** 4
    for codes, xcodes in (("0:1:3-14", ""), ("0:1:3-14,1:1:3-14", ""), ("0:1:3-14,1:1:8-14", "1:2:3-7,2:2:8-14"),
                          ("0:1:0-14,1:2:0-14,2:1:5-9", ""), ("0:-1:0-2,0:1:3-14", ""),
                          ("0:1:3-14,1:3:3-14,2:3:3-14,3:3:3-14,4:3:3-14,5:3:3-14,6:3:3-14", ""),
                          ("0:4:0-2,0:1:3-14,1:4:0-14,2:4:0-14,3:4:0-14,4:4:0-14,5:4:0-14,6:4:0-14", ""),
                          ("0:1:3-14,1:3:3-14,2:3:3-14,3:3:3-14", "1:1:3-14,2:2:3-14"),
                          ("0:1:3-14,1:4:3-14,2:4:3-14", "1:1:3-14,2:2:3-14"),
                          ("0:1:3-14,1:2:3-14,2:2:3-14", "2:2:3-14"),
                          ("0:1:3-14,1:2:3-14,2:2:3-14,3:2:3-14", "1:1:3-14,2:2:3-14"),
                          ("0:2:0-2,0:1:3-14,1:12:0-14,2:12:0-14,3:13:0-14,4:14:0-14", "")):
        deg, xdeg = parse_codes(codes), parse_xcodes(xcodes)
        frozen = np.zeros(t.shape, bool)
        value = np.zeros(t.shape, np.int64)
        res = xor_slack(c, blk, pos, xdeg) if xdeg.any() else None
        if prop_columns(deg).any():                          # chem_own2: the column slack of the property codes
            cs_ = column_slack(c, blk, pos, deg)
            res = cs_ if res is None else (res | cs_)
        for frac in (0.0, 0.3, 0.7, 1.0):                 # staged: frozen entries must be respected by later choices
            r_, m_ = residues(t, c, blk, pos, deg, xdeg, frozen, value, res)
            q = into_range(residue_round(t, r_, m_), m_).astype(np.int64)
            assert q.min() >= -HALF and q.max() <= HALF - 1
            assert np.array_equal(q[frozen], value[frozen]), "a frozen entry changed"
            assert np.all((q - r_) % m_ == 0)
            check_codes(q, addr, deg, xdeg)
            pri = rng.random(t.shape)
            if res is not None:
                pri = np.where(res, 2.0 + pri, pri)     # slack entries last
            new = (~frozen) & (pri <= np.quantile(pri, frac))
            frozen |= new
            value = np.where(new, q, value)
            t = np.where(frozen, t, np.clip(t + rng.normal(0, 0.7, t.shape), -64.4, 63.4))   # free angles drift
        check_codes(value, addr, deg, xdeg)
    # gauge: with nothing frozen every row is re-gauged exactly
    th = rng.uniform(-1, 1, (200, N - 1))
    v0 = decode_np(th, edges)
    fr = np.zeros(th.shape, bool)
    negat, flip_ok = negatable(fr, np.zeros(th.shape, np.int64), fr, NB)
    lo, hi = -0.5 - 0.5 / 128, 0.5 - 0.5 / 128
    out, nf, nn, ncl = regauge(th, fr, flip_ok, negat, edges, lo, hi)
    v1 = decode_np(out, edges)
    sgn = np.sign(np.sum(v0 * v1, axis=1))[:, None]
    assert np.abs(v0 - sgn * v1).max() < 1e-9 and out.min() >= lo and out.max() < hi and ncl == 0 and nf > 0
    # gauge with frozen, partly protected entries: protected ones never change, frozen ones only change sign,
    # and a row that needed no clip is reproduced exactly (crossings and negated frozen entries included)
    exact = 0
    for sigma, pfr, ppr in ((0.2, 0.5, 0.5), (0.03, 0.5, 0.5), (0.03, 0.3, 0.2), (0.02, 0.9, 0.1)):
        th2 = out + rng.normal(0, sigma, out.shape)
        fr = rng.random(th.shape) < pfr
        wv = np.round(np.where(fr, out, 0.0) * 128).astype(np.int64)
        th2 = np.where(fr, wv / 128.0, th2)
        prot = fr & (rng.random(th.shape) < ppr)
        negat, flip_ok = negatable(fr, wv, prot, NB)
        out2, nf, nn, ncl = regauge(th2, fr, flip_ok, negat, edges, lo, hi)
        assert np.array_equal(out2[prot & (wv != -HALF)], th2[prot & (wv != -HALF)])
        assert np.allclose(np.abs(out2[fr]), np.abs(th2[fr])) and out2.min() >= lo and out2.max() < hi
        assert np.array_equal(np.round(out2[fr] * 128), out2[fr] * 128), "a frozen entry left the word grid"
        va, vb = decode_np(th2, edges), decode_np(out2, edges)
        moved = np.abs(va - np.sign(np.sum(va * vb, axis=1))[:, None] * vb).max(axis=1) > 1e-9
        clip_rows = np.any(~fr & ((out2 == lo) | (out2 == hi - 1e-9)), axis=1)
        assert not np.any(moved & ~clip_rows), "a row that needed no clip was changed by the re-gauge"
        exact += int(np.sum(np.any(out2 != th2, axis=1) & ~clip_rows))
    assert exact > 20, exact
    # carrier pattern: any row carries the half turns exactly, and bit7 xor bit6 of the words is the pattern
    vec = rng.normal(size=(300, N))
    Kp = carrier_pattern(N + np.arange(300), [3, 5, 6, 9, 10])
    thp, _ = angles_and_mass(vec, edges, Kp)
    vv, vn = decode_np(thp + Kp, edges), vec / np.linalg.norm(vec, axis=1, keepdims=True)
    assert np.abs(vv - np.sign(np.sum(vv * vn, axis=1))[:, None] * vn).max() < 1e-9
    assert thp.min() >= -0.5 and thp.max() <= 0.5 and 0.15 < Kp.mean() < 0.3
    wds = np.mod(np.round(np.clip(thp, lo, hi - 1e-9) * 128).astype(np.int64) + 128 * Kp.astype(np.int64), 256)
    assert np.array_equal(((wds >> 7) ^ (wds >> 6)) & 1, Kp.astype(np.int64))
    # dealing over the true blocks
    R, B = 8, 143
    strength = rng.random((R, B)) ** 3
    order = deal_order(strength, R, B)
    a = N + np.arange(R * B)
    top = np.zeros(R * B, bool)
    for r in range(R):
        top[r * B + np.argsort(-strength[r])[:9]] = True
    per_block = np.bincount((a // 16)[top[order]], minlength=a.max() // 16 + 1)
    assert per_block.max() <= 2, per_block
    print(f"selftest ok: codes, staged residues with frozen entries (4 code sets), exact re-gauge "
          f"(last case: {nf} crossings, {nn} frozen negations, {ncl} clipped; {exact} re-gauged rows without a clip "
          f"reproduced exactly), carrier pattern exact ({Kp.mean():.3f} of the entries), dealing")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("init")
    ap.add_argument("out")
    ap.add_argument("--codes", default="0:1:3-14")
    ap.add_argument("--xcodes", default="")
    ap.add_argument("--no-deal", action="store_true")
    ap.add_argument("--order", default=None, help="a row order (npy, order[address - 76] = row) instead of the deal")
    ap.add_argument("--wires", default="", help="r50_chem_lbl: bit-0 delta wires k:j:0:m,... or @file (tools/r50_chem_lbl_chain.py): "
                    "the address of a row decided at freeze time; needs --no-relations --no-pairs --carriers none")
    ap.add_argument("--no-relations", action="store_true", help="control: round to nearest, codes only define the structure")
    ap.add_argument("--slack-at", type=float, default=1.0,
                    help="chem_own2 (R8, Nov 15): once the schedule passes this frozen share, the slack entries of the "
                         "XOR / column codes are frozen first (default 1.0: last, as before)")
    ap.add_argument("--col-slack", action="store_true",
                    help="chem_own2: slack entries for the low-degree-delta codes (the weakest entry per block-column and "
                         "check absorbs the parity fixes, frozen last; a frozen word of such a column may be negated "
                         "while its block-column's slack is free)")
    ap.add_argument("--carriers", default="3,5,6,9,10",
                    help="five low-address monomials (masks) carried on bit7 xor bit6 of the five word positions; "
                         "'none' = every word in the signed range. Default: the five lowest pair monomials, the "
                         "first ones our compiler would otherwise have to compute")
    ap.add_argument("--schedule", default="full", choices=["full", "long", "fine2", "fine4"],
                    help="chem_own2: fine2 / fine4 = the long schedule with 2 / 4 stages in place of each of its stages")
    ap.add_argument("--iters", type=int, default=3000)
    ap.add_argument("--late-iters", type=int, default=0,
                    help="iterations of the stages after --late-from (0 = --iters). Half of the rounding price is "
                         "paid before 4%% of the entries are frozen and 86%% before 20%%, so the early stages get the time")
    ap.add_argument("--late-from", type=float, default=0.25, help="frozen share after which a stage runs --late-iters")
    ap.add_argument("--gate", default="", help="'share:mHa[,share:mHa]': stop (exit 5) if the error after the stage "
                                               "that reaches a frozen share is above its value")
    ap.add_argument("--final", type=int, default=6000)
    ap.add_argument("--kappa", type=float, default=5e-5)
    ap.add_argument("--chunk", type=int, default=200)
    ap.add_argument("--recenter", type=int, default=1000)
    ap.add_argument("--max-minutes", type=float, default=25.0)
    ap.add_argument("--stop-at", type=float, default=2.0, help="probe: stop after the stage that reaches this frozen share")
    ap.add_argument("--x64", action="store_true",
                    help="stage optimiser in float64 (several times slower on a consumer card; float32 gradients stop "
                         "the L-BFGS of a long stage early, see the log)")
    ap.add_argument("--bits", type=int, default=8, help="word width")
    ap.add_argument("--no-pairs", action="store_true", help="ignore the twin pairs of INIT (they are then plain rows)")
    ap.add_argument("--pair-mode", choices=("joint", "same"), default="joint",
                    help="joint: floor / ceil combination with the least first-order error (the mechanism); same: both "
                         "rows take the same nearest word, i.e. the unsplit vector (the control)")
    ap.add_argument("--pair-coded", type=int, default=1, help="1: joint choice on coded columns too (each row in its own "
                    "residue class); 0: there each row takes its own nearest coded word (the run of 18:37Z)")
    ap.add_argument("--pair-weights", choices=("now", "start"), default="now",
                    help="weights in the joint choice: the current ones of the two rows, or the p of the start")
    ap.add_argument("--pair-step", type=float, default=3.4e-4,
                    help="cost of the step between master and twin: this times (q_m - q_s)^4 is added to the squared "
                         "first-order offset (steps^2); 3.4e-4 is the raw second-order term measured at 7 bits")
    ap.add_argument("--pmass", type=float, default=0.0, help="scale every angle variable also by 1/sqrt(max(mass, FLOOR)) (0 = off)")
    ap.add_argument("--pull", type=float, default=0.0, help="pull of every free angle to its word: final strength MU (1 = raw rounding energy; 0 = off)")
    ap.add_argument("--pull-start", type=float, default=0.0, help="strength at the start of the first stage (default MU / 1000)")
    ap.add_argument("--pull-pre", type=int, default=0, help="iterations of the first stage (nothing frozen) in which the pull rises")
    ap.add_argument("--pull-until", type=float, default=1.0, help="keep the pull in the refits while the frozen share is below this")
    ap.add_argument("--tail-pull", type=int, default=0, help="r47: with --pull: only the free words within this many steps of a "
                                                             "trivial word (0 or a quarter turn) are pulled, and to that trivial "
                                                             "word (the tail-aware refit: small tails go to zero or grow)")
    ap.add_argument("--hold", type=int, default=0, help="diagnostic: refit N iterations at the frozen set of state.npz, report, stop")
    ap.add_argument("--memory", type=int, default=30, help="pairs kept by L-BFGS")
    ap.add_argument("--aware", type=float, default=0.0, help="rounding-aware term: noise scale EPS (1 = the full rounding step)")
    ap.add_argument("--aware-pre", type=int, default=0, help="iterations of a first stage with nothing frozen")
    ap.add_argument("--aware-until", type=float, default=1.0, help="keep the term in the refits while the frozen share is below this")
    ap.add_argument("--aware-seed", type=int, default=0)
    ap.add_argument("--pairs-first", choices=("none", "uncoded", "all"), default="none",
                    help="probe: stage 0 freezes exactly the entries of the paired vectors on these columns")
    ap.add_argument("--log-every", type=int, default=0, help="chem_own2: log the error every N iterations of a stage")
    ap.add_argument("--rank-raw", default=None, help="chem_own2: (R*B, 75) raw diagonal H_kk (r48_chem_q7c_sdiag.py *_raw.npy)")
    ap.add_argument("--blend-lambda", type=float, default=0.25, help="chem_own2: --rank blend = S_kk + lambda H_kk")
    ap.add_argument("--absorb-ratio", type=float, default=0.3, help="chem_own2: --rank twophase: S_kk / H_kk below this = absorbable")
    ap.add_argument("--absorb-from", type=float, default=0.5, help="chem_own2: --rank twophase: absorbable entries held free "
                                                                     "until this frozen share, then frozen first")
    ap.add_argument("--switch-at", type=float, default=0.5, help="chem_own2: --rank sdiag_stiff: S_kk order below this "
                                                                 "frozen share, the stiffness order (present masses) above")
    ap.add_argument("--rank", default="impact", choices=("impact", "stiff", "stiffdyn", "sdiag", "blend", "twophase",
                                                         "sdiag_stiff", "sratio", "sdyn"),
                    help="chem_own2: order of freezing; impact = |w|^2 mass offset^2 (the old rule), stiff = |w|^2 mass "
                         "(offset ignored), sdiag = the weight-compensated stiffness S_kk of --rank-file "
                         "(tools/r48_chem_q7c_sdiag.py); the linear model prices stiff at 0.46 and sdiag at 0.33 of impact; "
                         "sdyn = S_kk recomputed from the PRESENT state in a background CPU thread (the start table is "
                         "--rank-file if given; every finished recompute is used from the next stage on)")
    ap.add_argument("--rank-file", default=None, help="chem_own2: (R*B, 75) npy in the row order of INIT (for --rank sdiag)")
    ap.add_argument("--sdyn-cpu", type=int, default=None, help="chem_own2: core of the --rank sdyn worker process (default: "
                                                                "the job's own cores)")
    ap.add_argument("--labels", default=None, help="chem_own2: orbital labels of the log-depth tree, perm[orbital] = label "
                                                   "(panel S3-CHEM-R, tools/panel3_chemr_labels.py)")
    ap.add_argument("--direct", action="store_true",
                    help="CPU only, no refit: round every entry at once (a coded table for pipeline tests, not a candidate)")
    ap.add_argument("--edges", default=None, help="r47: a general tree: npz with `edges` (75 x 2, parents first, root 0) "
                                                  "instead of the log-depth tree and --labels")
    ap.add_argument("--keep-states", default="", help="r47: frozen shares at which state.npz is also kept as state_<share>.npz")
    ap.add_argument("--fine", default=None, help="r47 (a bound, not a table): bool (R*B, 75) npy in the row order of INIT: "
                                                 "these words take a grid of HALF the step (one more bit)")
    ap.add_argument("--rank-pow", type=float, default=0.0, help="r47: the ranking value is divided by |w_row|^(2 * POW) "
                                                                "(0 = the builders' order; 1 = the row weight taken out)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    logf = open(os.path.join(a.out, "log.txt"), "a")

    def log(*x):
        s = " ".join(str(v) for v in x)
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    z = np.load(a.init)
    U, w = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
    R, B, C = w.shape
    U = U / np.linalg.norm(U, axis=2, keepdims=True)
    if a.order:
        order = np.load(a.order).astype(np.int64)
    elif a.no_deal:
        order = np.arange(R * B)
    else:
        order = deal_order(np.sum(w ** 2, axis=2), R, B)
    global HALF
    HALF = 1 << (a.bits - 2)
    DYN = None                                      # r50_chem_lbl: one DynChain per group when --wires is given
    pairs0 = np.array(z["pairs"], np.int64).reshape(-1, 2) if ("pairs" in z.files and not a.no_pairs) else np.zeros((0, 2), np.int64)
    if len(pairs0) and not a.order:
        # the twin sits on the address next to its master: swap it with whatever the deal put there
        inv = np.empty(R * B, np.int64)
        inv[order] = np.arange(R * B)
        taken = set(inv[pairs0.ravel()].tolist())
        for m_, s_ in pairs0:
            pm, ps = int(inv[m_]), int(inv[s_])
            for cand in (pm + 1, pm - 1):
                if cand // B == pm // B and 0 <= cand < R * B and (cand == ps or cand not in taken):
                    other = order[cand]
                    order[cand], order[ps] = s_, other
                    inv[s_], inv[other] = cand, ps
                    taken.discard(ps)
                    taken.add(cand)
                    break
    U = U.reshape(-1, N)[order].reshape(R, B, N)
    rank_sk = rank_hk = None
    if a.rank in ("sdiag", "blend", "twophase", "sdiag_stiff", "sratio") or (a.rank == "sdyn" and a.rank_file):
        # chem_own2: the static S_kk, dealt order
        rank_sk = np.load(a.rank_file).astype(np.float64).reshape(R * B, N - 1)[order]
    if a.rank in ("blend", "twophase", "sratio"):
        rf_ = a.rank_raw or (a.rank_file[:-4] + "_raw.npy")
        rank_hk = np.load(rf_).astype(np.float64).reshape(R * B, N - 1)[order]
    w = w.reshape(-1, C)[order].reshape(R, B, C)
    inv = np.empty(R * B, np.int64)
    inv[order] = np.arange(R * B)
    pm_rows, ps_rows = inv[pairs0[:, 0]], inv[pairs0[:, 1]]          # rows of the ordered table: masters, twins
    if len(pm_rows) and a.carriers not in ("", "none"):
        raise SystemExit("twin pairs need --carriers none for now: the two addresses of a pair carry different half "
                         "turns, so their angle tables differ by signs and 'tied' is not 'equal' (to be added)")
    edges = log_depth_tree(N)
    if a.labels:        # chem_own2: orbital labels of the tree (perm[orbital] = label; orbital inv[label] at a label)
        perm_ = np.load(a.labels).astype(np.int64)
        inv_ = np.empty(N, np.int64)
        inv_[perm_] = np.arange(N)
        assert inv_[0] == 0, "orbital 0 must keep label 0 (the tree's root)"
        edges = inv_[edges]
    if a.edges:         # r47: any spanning tree (edge order = column order)
        edges = np.array(np.load(a.edges)["edges"], np.int64)
        assert edges.shape == (N - 1, 2) and edges[0, 0] == 0
        seen_ = {0}
        for p_, q_ in edges:
            assert int(p_) in seen_ and int(q_) not in seen_, "edges must list parents first and form a tree"
            seen_.add(int(q_))
    SUB = np.ones((R * B, N - 1))
    if a.fine:          # r47: one more bit on a set of words (bound)
        SUB = np.where(np.load(a.fine).reshape(R * B, N - 1)[order], 2.0, 1.0)
    KEEP = [float(x) for x in a.keep_states.split(",") if x]
    NB = neighbour_matrix(edges)
    scale = float(1 << (a.bits - 1))
    deg, xdeg = parse_codes(a.codes), parse_xcodes(a.xcodes)
    pcbits = (deg != 0).sum(axis=1)                                  # per-column coded bits of every column
    PCOL = prop_columns(deg)                                         # chem_own2: columns with low-degree-delta codes
    xcol = np.repeat(xdeg[:, 1:].any(axis=1), 5)                     # columns of batches with an XOR code above bit 0
    ccols = np.flatnonzero(deg.any(axis=1) | np.repeat(xdeg.any(axis=1), 5))
    early = 0.5 / scale
    lo, hi = -0.5 - early, 0.5 - early
    addr = N + np.arange(R * B)
    blk, pos = addr // 16, addr % 16
    nblocks = int(blk.max()) + 1
    masks = [] if a.carriers in ("", "none") else [int(x) for x in a.carriers.split(",")]
    K = carrier_pattern(addr, masks)                                 # half turns: the angle of an entry is theta + K
    sched = FULL if a.schedule == "full" else LONG
    if a.schedule.startswith("fine"):       # chem_own2: every long stage cut into k equal steps of the frozen share
        kf, out_, prev_ = int(a.schedule[4:]), [], 0.0
        for f_ in LONG:
            out_ += [round(prev_ + (f_ - prev_) * j_ / kf, 7) for j_ in range(1, kf + 1)]
            prev_ = f_
        sched = tuple(out_)
    PRE = 1 if ((a.aware > 0 and a.aware_pre > 0) or (a.pull > 0 and a.pull_pre > 0)) else 0
    if PRE:
        sched = (0.0,) + tuple(sched)             # stage 0 = the rounding-aware descent, nothing frozen
    AWARE = a.aware > 0
    state_path = os.path.join(a.out, "state.npz")
    refc = Target()

    def finish(wv, w_, trace_, seconds):
        """Checks and writes the all-frozen table: signed words wv, weights w_."""
        assert wv.min() >= -HALF and wv.max() <= HALF - 1, "a word is outside the signed range"
        if not a.no_relations:
            check_codes(wv, addr, deg, xdeg)
        order_ = order
        if DYN is not None:     # r50_chem_lbl: every row to the address its relations chose (inside its group)
            w9perm = np.concatenate([w9g * B + np.argsort(DYN[w9g].mrow[:B]) for w9g in range(R)])   # w9perm[new place] = old row
            wv, w_, order_ = wv[w9perm], w_.reshape(R * B, -1)[w9perm].reshape(w_.shape), order[w9perm]
            w9bad = 0
            for (w9k, w9j, _, w9mask) in WIRES:
                w9bad += int(np.sum(((wv[:, 5 * w9k + w9j] ^ wv[:, 5 * w9k - 5 + w9j]) & 1) != _wlin(addr, w9mask)))
            log(f"wires: rows placed; relations broken {w9bad} of {len(WIRES) * R * B}" + ("" if w9bad == 0 else
                " (REPAIR: tools/r50_chem_lbl_wtab2.py --no-assign --floor 1e-4 on the finished table)"))
        top = 1 << (a.bits - 1)
        words = np.mod(wv + top * K.astype(np.int64), 2 * top)
        carried = ((words >> (a.bits - 1)) ^ (words >> (a.bits - 2))) & 1
        assert np.array_equal(carried, K.astype(np.int64)), "the two top bits do not carry the carrier pattern"
        u = decode_np(words / scale, edges).reshape(R, B, N)
        e_fin = refc.error2_np(u, w_)
        np.savez(os.path.join(a.out, "rounded.npz"), words=words, W=w_, U=u, order=order_, codes=deg, xcodes=xdeg,
                 carriers=np.array(masks, np.int64), tree_edges=edges, trace=np.array(trace_), bits=a.bits,
                 pairs=np.stack([pm_rows, ps_rows], axis=1))
        out = {"codes": a.codes, "xcodes": a.xcodes, "relations": not a.no_relations, "schedule": a.schedule,
               "iters": a.iters, "carriers": a.carriers, "direct": bool(a.direct), "error_2_mHa": 1e3 * e_fin,
               "shape": [R, B, C], "init": os.path.abspath(a.init), "bits": a.bits, "pairs": int(len(pm_rows))}
        json.dump(out, open(os.path.join(a.out, "result.json"), "w"))
        log(f"all frozen, every word in the signed range{', codes verified on every block' if not a.no_relations else ''}"
            f", carriers {masks}: error_2 {1e3 * e_fin:.4f} mHa  ({seconds:.0f} s)")

    if a.direct:        # CPU only: every entry rounded at once and nothing refitted (a table for pipeline tests)
        theta, mass = angles_and_mass(U.reshape(-1, N), edges, K)
        fr0 = np.zeros(theta.shape, bool)
        ng, fo = negatable(fr0, np.zeros(theta.shape, np.int64), fr0, NB)
        theta = regauge(theta, fr0, fo, ng, edges, lo, hi)[0]
        c = (np.linalg.norm(w, axis=2).reshape(-1)[:, None] ** 2) * mass + 1e-30
        res = xor_slack(c, blk, pos, xdeg) if xdeg.any() else None
        r_, m_ = residues(theta * scale, c, blk, pos, deg, xdeg, fr0, np.zeros(theta.shape, np.int64), res)
        log(f"direct rounding of {a.init} {(R, B, C)} codes [{a.codes}] xcodes [{a.xcodes}] carriers {masks}: "
            f"the fit itself {1e3 * refc.error2_np(U, w):.4f} mHa")
        finish(into_range(residue_round(theta * scale, r_, m_), m_).astype(np.int64), w, [], 0.0)
        return

    import jax
    jax.config.update("jax_enable_x64", bool(a.x64))
    jax.config.update("jax_default_matmul_precision", "highest")
    import jax.numpy as jnp
    import optax
    DT = jnp.float64 if a.x64 else jnp.float32
    lv = levels(edges)

    def d_penalty(v):      # sign mixing 1/2 sum_rc ((sum_b |w|)^2 - (sum_b w)^2): the part of Lambda the fit controls
        return 0.5 * jnp.sum(jnp.sum(jnp.abs(v), axis=1) ** 2 - jnp.sum(v, axis=1) ** 2)
    i_, j_ = jnp.asarray(refc.i), jnp.asarray(refc.j)
    sM32 = jnp.asarray(np.sqrt(refc.M), DT)
    sM64 = np.sqrt(refc.M)
    sw = float(np.sqrt(np.mean(w ** 2)))
    nbw = np.linalg.norm(w, axis=2).reshape(-1)
    nbw_pair = nbw.copy()                           # weight norm of the vector a row belongs to (pair: both rows)
    nbw_pair[pm_rows] = nbw_pair[ps_rows] = nbw[pm_rows] + nbw[ps_rows]
    pshare = nbw[pm_rows] / np.maximum(nbw[pm_rows] + nbw[ps_rows], 1e-300)     # p of every pair (the master's share)
    srow = 1.0 / np.maximum(nbw_pair / np.sqrt(np.mean(nbw ** 2)), 3e-2)
    mst = np.arange(R * B)                          # the row whose angles a row follows (itself unless a tied twin)
    mst[ps_rows] = pm_rows
    ucol = np.ones(N - 1, bool)                     # columns without any code: joint floor / ceil rounding of a pair
    ucol[ccols] = False
    rank = np.argsort(-nbw_pair ** 2 + 1e-30 * np.arange(R * B), kind="stable")
    vec_first = np.ones(R * B, bool)
    vec_first[ps_rows] = False                      # one entry per vector in the strength ranking
    ranked_rows = [r_ for r_ in rank if vec_first[r_]]
    nvec = len(ranked_rows)
    cls = []                                        # rows of the 1 / 2 / 4 / 6 % strongest vectors
    for share_ in (0.01, 0.02, 0.04, 0.06):
        m_ = np.zeros(R * B, bool)
        m_[ranked_rows[:int(round(share_ * nvec))]] = True
        m_[ps_rows] = m_[pm_rows]
        cls.append(m_)
    attributed, total_inc = np.zeros(4), 0.0

    def stage_scale(th_):
        """Scale of the angle variables for one stage: the row scale, with --pmass also 1/sqrt(mass of the rotation)."""
        sc_ = np.repeat(srow[:, None], N - 1, axis=1)
        if a.pmass > 0:
            m_ = angles_and_mass(decode_np(th_ + K, edges), edges)[1]
            sc_ = sc_ / np.sqrt(np.maximum(m_, a.pmass))
            sc_[ps_rows] = np.where((mst[ps_rows] == pm_rows)[:, None], sc_[pm_rows], sc_[ps_rows])   # a tied twin: its master's scale
        return sc_

    def reference(th, w_, full=False):
        """Float64 reference at the angles th (+ half turns K) and the weights w_: A0, the weighted residual r0 and
        the error; with full=True also what the difference decode needs (cos and sin of every rotation, the
        amplitude that enters it, the vectors and their pairs)."""
        ang = np.pi * (th + K)
        c0, s0 = np.cos(ang), np.sin(ang)
        v = np.zeros((th.shape[0], N))
        v[:, 0] = 1.0
        ain = np.zeros_like(th)
        for t, (p_, q_) in enumerate(edges):
            a_ = v[:, p_].copy()
            ain[:, t] = a_
            v[:, p_] = c0[:, t] * a_
            v[:, q_] = s0[:, t] * a_
        u = v.reshape(R, B, N)
        p0 = u[:, :, refc.i] * u[:, :, refc.j]
        A0 = np.einsum("rbc,rbp->rcp", w_, p0).reshape(R * C, -1)
        r0 = sM64 * (A0.T @ A0 - refc.T)
        eps = float(np.sqrt(np.sum(r0 * r0)))
        return (A0, r0, eps, (c0, s0, ain, u, p0)) if full else (A0, r0, eps)

    opt = optax.lbfgs(memory_size=a.memory)

    @jax.jit
    def run(params, st, A0, r0, eps0, off, fmask, wall, own0, th0, c0, s0, ain0, u0, p0, v0, mst_j, nz, sc_j, pw, qt):
        def f(p):
            # The variables are CHANGES since the float64 reference: p["tau"] of the scaled angles, p["v"] of the
            # scaled weights; both are zero at every recentering and the float64 master copies are kept on the
            # CPU. (Until 2026-11-11 15:40Z the variables were the angles and weights themselves in float32: a
            # late L-BFGS step is smaller than the last digit of such a number and was rounded away.)
            dfree = (p["tau"] * sc_j)[mst_j]                  # a tied twin moves with its master
            dwall = jnp.clip(own0 + dfree, lo, hi - 1e-7) - th0
            dth = (1.0 - fmask) * jnp.where(wall > 0, dwall, dfree)
            # The change dA of the table since the reference is computed AS A CHANGE. (Until 2026-11-11 15:30Z it
            # was A(p) - A0 with A(p) evaluated in float32: an error of 1e-7 |A| in every entry whatever the size
            # of the true change, i.e. about 5e-9 Ha of noise in the value, the size of one iteration's progress
            # late in a long stage. In float64 the same stage reached in 12,000 iterations what float32 missed
            # in 40,000.) Rotation k changes by (dc, ds) = (cos, sin)(a0 + d) - (cos, sin)(a0), exact to float32
            # relative accuracy through the half-angle forms; the vectors by du, the pairs by dp.
            def core(dth_, pv):
                half = (0.5 * jnp.pi) * dth_
                sh, ch = jnp.sin(half), jnp.cos(half)
                dc = -2.0 * (s0 * ch + c0 * sh) * sh
                ds = 2.0 * (c0 * ch - s0 * sh) * sh
                du = jnp.zeros((dth_.shape[0], N), dth_.dtype)
                for ts in lv:
                    pp, qq = edges[ts, 0], edges[ts, 1]
                    da, a0 = du[:, pp], ain0[:, ts]
                    du = du.at[:, pp].set(dc[:, ts] * a0 + (c0[:, ts] + dc[:, ts]) * da) \
                           .at[:, qq].set(ds[:, ts] * a0 + (s0[:, ts] + ds[:, ts]) * da)
                du = du.reshape(R, B, N)
                ui, uj, dui, duj = u0[:, :, i_], u0[:, :, j_], du[:, :, i_], du[:, :, j_]
                dp = dui * uj + (ui + dui) * duj
                wv_ = (v0 + pv) * sw
                dA = (jnp.einsum("rbc,rbp->rcp", pv * sw, p0) + jnp.einsum("rbc,rbp->rcp", wv_, dp)).reshape(R * C, -1)
                S = A0.T @ dA
                dr = sM32 * (S + S.T + dA.T @ dA)
                delta = jnp.sum(dr * (2.0 * r0 + dr))
                return delta / (eps0 + jnp.sqrt(jnp.maximum(eps0 * eps0 + delta, 0.0)))
            wv = (v0 + p["v"]) * sw
            if AWARE:       # the two-point form: the rounding noise still to come on every free angle (nz is zero elsewhere)
                deps = 0.5 * (core(dth + nz, p["v"]) + core(dth - nz, p["v"]))
            else:
                deps = core(dth, p["v"])
            # the pull to the words, as its CHANGE since the reference (pw is zero when it is off or the entry is frozen)
            pen = jnp.sum(pw * dth * (2.0 * (th0 - qt) + dth))
            return off + deps + pen + a.kappa * d_penalty(wv)
        vg = optax.value_and_grad_from_state(f)

        def body(_, carry):
            p, s = carry
            value, grad = vg(p, state=s)
            upd, s = opt.update(grad, s, p, value=value, grad=grad, value_fn=f)
            return optax.apply_updates(p, upd), s
        return jax.lax.fori_loop(0, a.chunk, body, (params, st))

    if os.path.exists(state_path):
        with np.load(state_path) as s:
            theta, w, frozen, qv, stage0, spent = s["theta"], s["w"], s["frozen"], s["qv"], int(s["stage"]), float(s["spent"])
            RES = s["res"] if "res" in s.files else None
            trace = [tuple(x) for x in s["trace"]] if "trace" in s.files else []
            if "mst" in s.files:
                mst = np.array(s["mst"])
        log(f"resume at stage {stage0} (frozen {frozen.mean():.4f}, {spent / 60:.1f} minutes spent so far)")
    else:
        theta, _ = angles_and_mass(U.reshape(-1, N), edges, K)
        fr0 = np.zeros(theta.shape, bool)
        ng, fo = negatable(fr0, np.zeros(theta.shape, np.int64), fr0, NB)
        theta, _, _, _ = regauge(theta, fr0, fo, ng, edges, lo, hi)
        frozen = np.zeros(theta.shape, bool)
        qv = np.zeros_like(theta)
        stage0, spent, RES, trace = 0, 0.0, None, []
        log(f"init {a.init} {(R, B, C)} codes [{a.codes}] ({int((deg != 0).sum())} coded column-bits) xcodes [{a.xcodes}] "
            f"relations {not a.no_relations} dealt {not a.no_deal and not a.order} carriers {masks} "
            f"({K.mean():.4f} of the entries carry a half turn) schedule {a.schedule} x {a.iters}")
        if AWARE:
            log(f"rounding-aware term: EPS {a.aware} (noise +-{0.5 * a.aware:.4f} of a step on every free angle, two-point form), "
                f"first stage {a.aware_pre} iterations with nothing frozen, kept in the refits below frozen {a.aware_until}")
        log(f"start {1e3 * reference(theta, w)[2]:.4f} mHa (the fit itself: {1e3 * refc.error2_np(U, w):.4f}); word width {a.bits} "
            f"bits; twin pairs {len(pm_rows)}" + (f" (p {pshare.min():.3f}-{pshare.max():.3f}, largest angle difference inside a "
            f"pair {np.abs(theta[pm_rows] - theta[ps_rows]).max():.1e})" if len(pm_rows) else ""))
    if a.wires:     # r50_chem_lbl: bit-0 delta wires on address bits (construction 9 / 9b)
        from r50_chem_lbl_chain import DynChain, lin as _wlin, parse as _wparse
        assert a.no_relations and len(pm_rows) == 0 and not np.any(K), "--wires needs --no-relations --no-pairs --carriers none"
        WIRES = _wparse(a.wires)
        DYN = [DynChain(WIRES, addr[w9g * B:(w9g + 1) * B], B) for w9g in range(R)]
        WCOL = np.zeros(N - 1, bool)
        WCOL[list(DYN[0].col)] = True
        w9zq = np.round(qv * scale).astype(np.int64)
        w9rep = 0
        for w9i, w9c in zip(*np.nonzero(frozen & WCOL[None, :])):      # a resumed run: the frozen wire words again
            DYN[int(w9i) // B].commit(int(w9i) % B, int(w9c), int(w9zq[w9i, w9c]))
            w9rep += 1
        log(f"wires: {len(WIRES)} in {len(DYN[0].chains)} chains on {int(WCOL.sum())} columns; the address of a row is decided "
            f"at freeze time; {w9rep} frozen wire words replayed")
    SDYN = a.rank == "sdyn"
    if SDYN:        # chem_own2: S_kk from the present state, recomputed in a SPAWNED worker process while the GPU
        # refits (a thread held the GIL: the main loop then got a third of a core and the card idled, 2x slower [M])
        import multiprocessing as mp_
        from r48_chem_q7c_sdiag import state_table as _state_table
        sdyn_pool = mp_.get_context("spawn").Pool(1)

        def sdyn_start(th_, w_):
            return sdyn_pool.apply_async(_state_table, (th_, w_, K, edges, a.sdyn_cpu)), time.time()

        sdyn_path = os.path.join(a.out, "sdyn.npy")
        if stage0 > 0 and os.path.exists(sdyn_path):
            rank_sk = np.load(sdyn_path)
            log("sdyn: the last S_kk table of this run reloaded")
        elif rank_sk is None or stage0 > 0:
            t_ = time.time()
            rank_sk = _state_table(np.where(frozen, qv, theta), w, K, edges)
            np.save(sdyn_path, rank_sk)
            log(f"sdyn: S_kk computed from the present state ({time.time() - t_:.0f} s)")
        sdyn_job, sdyn_from = None, None

    def rules():
        """Frozen entries that must keep their word: per-column codes on bits 0..nb-1 survive a negation only if
        w = 0 modulo 2^(nb-1); in a batch with an XOR code above bit 0 a frozen word is protected once the slack
        entries of its block and batch are frozen (before that they re-absorb the changed relation)."""
        wv = np.round(qv * scale).astype(np.int64)
        modh = np.where(pcbits >= 2, 1 << np.maximum(pcbits - 1, 0), 1)[None, :]
        protected = frozen & (pcbits[None, :] >= 2) & (np.mod(wv, modh) != 0)
        if a.col_slack and RES is not None and PCOL.any():
            # chem_own2: a negation keeps bit 0 and changes bits >= 1, whose low-degree-delta parities the block-column's
            # slack re-absorbs while it is free: protect only once a slack entry of the block-column is frozen
            fbc = np.zeros((nblocks, N - 1), bool)
            rr_, cl_ = np.nonzero(frozen & RES & PCOL[None, :])
            fbc[blk[rr_], cl_] = True
            protected = np.where(PCOL[None, :], protected & fbc[blk], protected)
        if xcol.any():
            if RES is None:
                protected |= frozen & xcol[None, :] & (wv != 0)
            else:
                fb = np.zeros((nblocks, NB_BATCH), bool)
                rr, cl = np.nonzero(frozen & RES)
                fb[blk[rr], cl // 5] = True
                protected |= frozen & xcol[None, :] & fb[blk][:, np.arange(N - 1) // 5]
        return negatable(frozen, wv, protected, NB)

    def coded_round(now, c):
        """Nearest allowed word of every entry (units of pi), its residue and the modulus of its column."""
        t = now * scale
        q = into_range(np.round(t), 1.0)
        if a.fine:      # r47: half steps on the fine words (no codes with it)
            return np.round(t * SUB) / SUB / scale, np.zeros(t.shape, np.int64), np.ones(N - 1, np.int64)
        if a.no_relations or len(ccols) == 0:
            return q / scale, np.zeros(t.shape, np.int64), np.ones(N - 1, np.int64)
        r_, m_ = residues(t, c, blk, pos, deg, xdeg, frozen, np.round(qv * scale).astype(np.int64), RES)
        v = into_range(residue_round(t, r_, m_), m_)
        q[:, ccols] = v[:, ccols]
        return q / scale, r_, m_

    jcols = (np.ones(N - 1, bool) if a.pair_coded else ucol.copy()) & ~xcol    # columns with a joint choice
    prow = np.zeros(R * B, bool)
    prow[pm_rows] = prow[ps_rows] = True

    def joint_round(now, q, rr, mm):
        """Twin pairs still tied, one decision per entry. Master and twin each take the word below or above the angle
        in their own residue class (step 1 on an uncoded column, 2 under a bit-0 code); the combination is the one
        with the least first-order error sum_c (w_mc d_m + w_sc d_s)^2 / (|w_m| + |w_s|)^2 (d in steps) plus
        --pair-step times the fourth power of the step between the two words. Mode `same`: both rows take the same
        nearest word wherever their residues agree. Returns the cost (units of pi^2) of every tied master entry,
        the squared weight norm of every pair, and which pairs are tied."""
        cost_pi2 = np.zeros(now.shape)
        tied = mst[ps_rows] == pm_rows
        if not tied.any():
            return cost_pi2, np.zeros(0), tied
        m_, s_ = pm_rows[tied], ps_rows[tied]
        if a.pair_weights == "now":
            wf = w.reshape(-1, C)
            qa, qb, qd = np.sum(wf[m_] ** 2, axis=1), np.sum(wf[m_] * wf[s_], axis=1), np.sum(wf[s_] ** 2, axis=1)
        else:
            qa, qb, qd = nbw[m_] ** 2, nbw[m_] * nbw[s_], nbw[s_] ** 2
        tot = (np.sqrt(qa) + np.sqrt(qd)) ** 2
        qa, qb, qd = (qa / tot)[:, None], (qb / tot)[:, None], (qd / tot)[:, None]
        t = now[m_] * scale
        hi_, lo_ = HALF - 1, -HALF

        def first(qm, qs):
            dm, ds = qm - t, qs - t
            return qa * dm * dm + 2.0 * qb * dm * ds + qd * ds * ds

        if a.pair_mode == "same":
            cm, cs = q[m_] * scale, q[s_] * scale
            eq = rr[m_] == rr[s_]
            cs = np.where(eq & jcols[None, :], cm, cs)
            use = jcols[None, :] & ~frozen[m_]
            q[s_] = np.where(use, cs / scale, q[s_])
            cost_pi2[m_] = np.where(use, first(cm, cs) / scale ** 2, 0.0)
            return cost_pi2, tot, tied
        mod = np.broadcast_to(np.asarray(mm, np.float64)[None, :], t.shape)
        lom = rr[m_] + mod * np.floor((t - rr[m_]) / mod)
        los = rr[s_] + mod * np.floor((t - rr[s_]) / mod)
        best, bm, bs = np.full(t.shape, np.inf), q[m_] * scale, q[s_] * scale
        for am in (0.0, 1.0):
            for as_ in (0.0, 1.0):
                qm, qs = lom + am * mod, los + as_ * mod
                cost = first(qm, qs) + a.pair_step * (qm - qs) ** 4
                ok = (qm <= hi_) & (qs <= hi_) & (qm >= lo_) & (qs >= lo_)
                better = ok & (cost < best)
                best, bm, bs = np.where(better, cost, best), np.where(better, qm, bm), np.where(better, qs, bs)
        use = jcols[None, :] & np.isfinite(best) & ~frozen[m_]
        q[m_] = np.where(use, bm / scale, q[m_])
        q[s_] = np.where(use, bs / scale, q[s_])
        cost_pi2[m_] = np.where(use, best / scale ** 2, 0.0)
        return cost_pi2, tot, tied

    t0 = time.time()
    sc = stage_scale(theta)
    hnorm = None
    tau64, v64 = theta / sc, w / sw                       # float64 master copies of the scaled angles and weights
    zero = {"tau": jnp.zeros(theta.shape, DT), "v": jnp.zeros(w.shape, DT)}
    HOLD = a.hold > 0
    if HOLD:        # an extra refit at the present frozen set, before the stage the state would continue with
        sched = tuple(sched[:stage0]) + (float(frozen.mean()),) + tuple(sched[stage0:])
    for k in range(stage0, len(sched)):
        frac = sched[k]
        # Zero of the float32 objective: the error and the sign mixing BEFORE this stage's freeze. The stage ends
        # close to them, so the value L-BFGS sees is small there and its last float32 digit is far below the
        # decrease of one iteration. (Until 2026-11-11 15:15Z the zero was the error right after the freeze,
        # 50-60 mHa higher: the value then sat near -0.05 Ha, one float32 step was 4e-9 Ha, the decrease per
        # iteration late in a long stage 2e-9 Ha, and the line search ran on rounding. Measured on the public
        # pipeline's point, stage 0: 14.365 mHa after 40,000 iterations; float64: 14.321 after 12,000.)
        e_base = reference(theta, w)[2]
        d_base = 0.5 * float(np.sum(np.abs(w).sum(axis=1) ** 2 - w.sum(axis=1) ** 2))
        need = int(round(frac * theta.size)) - int(frozen.sum())
        if HOLD and k == stage0:
            need = 0
        if need > 0:
            now = np.where(frozen, qv, theta)
            mass = angles_and_mass(decode_np(now + K, edges), edges)[1]
            c = (nbw[:, None] ** 2) * mass + 1e-30
            if RES is None and (xdeg.any() or a.col_slack):
                RES = xor_slack(c, blk, pos, xdeg) if xdeg.any() else np.zeros(c.shape, bool)
                if a.col_slack:
                    CS_ = column_slack(c, blk, pos, deg)
                    RES = RES | CS_
                    log(f"column slack (low-degree-delta codes): {int(CS_.sum())} entries in "
                        f"{int(prop_columns(deg).sum())} columns")
                log(f"slack entries: {int(RES.sum())} ({RES.mean():.4f} of all entries)")
            q, res_r, res_m = coded_round(now, c)
            clash = int(np.sum(frozen & (np.round(q * scale) != np.round(qv * scale))))
            if clash:
                log(f"WARNING stage {k}: {clash} frozen entries disagree with every codeword (the codes cannot hold)")
            jcost, jtot, tied = joint_round(now, q, res_r, res_m)
            impact = c * (q - now) ** 2
            if a.rank == "stiff":       # chem_own2: decide in order of stiffness, whatever the present offset
                impact = c.copy()
            elif a.rank == "sdiag":     # chem_own2: decide in order of the weight-compensated stiffness S_kk
                impact = rank_sk + 1e-30 * c
            elif a.rank == "sdyn":      # chem_own2: S_kk of the present state (the latest finished recompute)
                if sdyn_job is not None and sdyn_job[0].ready():
                    try:
                        new_ = sdyn_job[0].get()
                    except Exception as ex_:        # noqa: BLE001 - logged, the run goes on with the last table
                        log(f"sdyn: WARNING the recompute failed ({ex_!r}); started again")
                        new_ = None
                    if new_ is not None:
                        prev_ = rank_sk
                        fr_ = ~frozen
                        rho_ = np.corrcoef(np.log(new_[fr_] + 1e-30), np.log(prev_[fr_] + 1e-30))[0, 1] if fr_.any() else 1.0
                        grew_ = int(np.sum(fr_ & (new_ > 10.0 * prev_ + 1e-20)))
                        rank_sk = new_
                        np.save(sdyn_path, rank_sk)
                        log(f"sdyn: S_kk of the state at the start of stage {sdyn_from} ({time.time() - sdyn_job[1]:.0f} s "
                            f"after it was asked for) in use from stage {k}: corr(log) with the table before {rho_:.3f} "
                            f"over the free entries; {grew_} free entries more than 10x stiffer")
                    sdyn_job = None
                if sdyn_job is None:
                    sdyn_job, sdyn_from = sdyn_start(now.copy(), w.copy()), k
                impact = rank_sk + 1e-30 * c
            elif a.rank == "sratio":    # chem_own2: absorbable fraction S_kk / H_kk of the start x the PRESENT stiffness
                ratio_ = np.where(rank_hk > 1e-12 * float(rank_hk.max()), np.clip(rank_sk / np.maximum(rank_hk, 1e-300), 0.0, 1.0), 1.0)
                impact = ratio_ * c
            elif a.rank == "stiffdyn":  # chem_own2: stiffness with the PRESENT weights and masses (ranking only)
                impact = (np.linalg.norm(w, axis=2).reshape(-1)[:, None] ** 2) * mass + 1e-30
            elif a.rank == "sdiag_stiff":   # chem_own2: S_kk order to --switch-at, then the stiffness order
                impact = (rank_sk + 1e-30 * c) if float(frozen.mean()) < a.switch_at else c.copy()
            elif a.rank == "blend":     # chem_own2: S_kk + lambda H_kk (absorbable stiff entries move to the middle)
                impact = rank_sk + a.blend_lambda * rank_hk + 1e-30 * c
            elif a.rank == "twophase":  # chem_own2: absorbable entries held free to --absorb-from, then frozen first
                absorb_ = rank_sk < a.absorb_ratio * rank_hk
                if float(frozen.mean()) < a.absorb_from:
                    impact = np.where(absorb_, 1e-30 * c, rank_sk + 1e-30 * c)
                else:
                    impact = np.where(absorb_, 1e3 * (1.0 + rank_hk / max(float(rank_hk.max()), 1e-300)),
                                      rank_sk + 1e-30 * c)
            if a.rank_pow:      # r47: the row weight taken out of the ranking value
                impact = impact / np.maximum(np.linalg.norm(w, axis=2).reshape(-1)[:, None] ** (2.0 * a.rank_pow), 1e-300)
            if tied.any():      # a tied pair is one decision: its impact sits on the master, the twin follows
                m_, s_ = pm_rows[tied], ps_rows[tied]
                cpair = jtot[:, None] * mass[m_]
                impact[m_] = np.where(jcols[None, :], cpair * jcost[m_], impact[m_] + impact[s_])
                impact[s_] = -1e300
            if RES is not None:
                if frac > a.slack_at:   # chem_own2 (R8, Nov 15): past --slack-at the slack entries go FIRST, while
                    # free entries are left to compensate their parity-forced moves (the end-of-run "slack wall")
                    impact = np.where(RES, 1e300, impact)
                else:
                    impact = np.where(RES, -1.0 / (1.0 + impact), impact)      # slack entries are frozen last
            impact[frozen] = -np.inf
            was = frozen.copy()
            if a.pairs_first != "none" and k == 0 and tied.any():      # probe: the paired vectors, and only they
                sel = np.zeros(theta.shape, bool)
                sel[pm_rows[tied]] = (ucol if a.pairs_first == "uncoded" else np.ones(N - 1, bool))[None, :]
                pick = np.flatnonzero(sel & ~frozen)
            else:
                pick = np.argsort(-impact, axis=None)[:need]
            if DYN is not None:     # r50_chem_lbl: a wire word takes the parity its row's relations already decide
                w9t = now * scale
                w9n = w9f = w9m = 0
                for w9flat in pick:
                    w9i, w9c = divmod(int(w9flat), N - 1)
                    if not WCOL[w9c]:
                        continue
                    w9z = int(round(float(q[w9i, w9c]) * scale))
                    w9p = DYN[w9i // B].ask(w9i % B, w9c)
                    w9n += 1
                    if w9p is not None:
                        w9f += 1
                        if (w9z & 1) != w9p:
                            w9lo = int(np.floor(w9t[w9i, w9c]))
                            w9cand = [w9v for w9v in (w9lo - 1, w9lo, w9lo + 1, w9lo + 2) if (w9v & 1) == w9p and -HALF <= w9v <= HALF - 1]
                            w9z = min(w9cand, key=lambda w9v: abs(w9v - w9t[w9i, w9c]))
                            q[w9i, w9c] = w9z / scale
                            w9m += 1
                    DYN[w9i // B].commit(w9i % B, w9c, w9z)
                if w9n:
                    log(f"wires: {w9n} wire words frozen in this stage, {w9f} with a decided parity, {w9m} of them moved off their nearest word")
            frozen.flat[pick] = True
            qv.flat[pick] = q.flat[pick]
            pr_, pc_ = np.unravel_index(pick, theta.shape)
            if tied.any():
                twin_of = np.full(R * B, -1)
                twin_of[pm_rows[tied]] = ps_rows[tied]
                has = twin_of[pr_] >= 0
                frozen[twin_of[pr_[has]], pc_[has]] = True
                qv[twin_of[pr_[has]], pc_[has]] = q[twin_of[pr_[has]], pc_[has]]
            imp_ = np.maximum(impact.flat[pick], 0.0)
            stage_share = np.array([imp_[m_[pr_]].sum() for m_ in cls]) / max(imp_.sum(), 1e-300)
            pair_note = ""
            if len(pm_rows):    # raw kick of the pairs' own entries (everything else as before the freeze)
                newp = frozen & ~was & prow[:, None]
                e_p = reference(np.where(newp, qv, now), w)[2]
                stepped = newp[pm_rows] & (np.round(qv[pm_rows] * scale) != np.round(qv[ps_rows] * scale))
                pair_note = (f"; pairs' entries frozen {int(newp[pm_rows].sum())} ({int(stepped.sum())} with a step between "
                             f"master and twin), their raw kick {1e6 * (e_p * e_p - e_base * e_base):.1f} mHa^2")
        else:
            stage_share = np.zeros(4)
            pair_note = ""
        negat, flip_ok = rules()
        tied = mst[ps_rows] == pm_rows
        if tied.any():          # a tied pair takes the same gauge moves: a crossing only where both rows allow it
            m_, s_ = pm_rows[tied], ps_rows[tied]
            both = flip_ok[m_] & flip_ok[s_]
            flip_ok[m_], flip_ok[s_] = both, both
            theta[s_] = np.where(frozen[s_], theta[s_], theta[m_])
        wall = (~frozen) & ~flip_ok
        sc = stage_scale(np.where(frozen, qv, theta))
        sc_dev = jnp.asarray(sc, DT)
        tau64 = np.where(frozen, qv, theta) / sc
        params = zero
        st = opt.init(params)
        iters = a.final if frac >= 1.0 else (a.late_iters if a.late_iters and frac > a.late_from else a.iters)
        if HOLD and k == stage0:
            iters = a.hold
        if PRE and k == 0:
            iters = max(a.aware_pre if a.aware > 0 else 0, a.pull_pre if a.pull > 0 else 0)
        pull_on = a.pull > 0 and frac < 1.0 and (frac < a.pull_until or (PRE and k == 0))
        aware_on = AWARE and frac < 1.0 and (frac < a.aware_until or (PRE and k == 0))
        rng_n = np.random.default_rng([a.aware_seed, k])
        done, e_first, unstuck, e_blk = 0, None, 0, None
        blk_start = None        # chem_own2 (R8, Nov 15): float64 state at the start of the present refit block
        fm, wl = jnp.asarray(frozen, DT), jnp.asarray(wall, DT)
        while done < iters:
            if tied.any():
                tau64[ps_rows[tied]] = np.where(frozen[ps_rows[tied]], tau64[ps_rows[tied]], tau64[pm_rows[tied]])
            own = tau64 * sc
            stuck = wall & ((own < lo) | (own > hi - 1e-9))     # zero gradient past the clip: put it on the boundary
            if stuck.any():
                own = np.where(stuck, np.clip(own, lo, hi - 1e-9), own)
                tau64 = own / sc
                unstuck += int(stuck.sum())
            th = np.where(frozen, qv, np.where(wall, np.clip(own, lo, hi - 1e-9), own))
            A0, r0, eps0, (c0, s0, ain0, u0, p0) = reference(th, v64 * sw, full=True)
            if e_first is None:
                e_first = eps0
            # chem_own2 (R8, Nov 15): a block that RAISED the error (or made it non-finite) is undone. Until now the
            # test below also ended the stage on a rise and kept the bad state: request G's stage 87 (Nov 14) went
            # from 19.41 to 110.9 mHa in one block, the stage ended after 51 s, and the run carried that state on.
            if e_blk is not None and blk_start is not None and not aware_on and not pull_on and \
                    (not np.isfinite(eps0) or eps0 - e_blk > 1e-7):
                log(f"  stage {k}: a refit block raised error_2 from {1e3 * e_blk:.4f} to {1e3 * eps0:.4f} mHa: "
                    f"undone, the stage ends at the block's start")
                tau64, v64 = blk_start
                own = tau64 * sc
                th = np.where(frozen, qv, np.where(wall, np.clip(own, lo, hi - 1e-9), own))
                break
            blk_start = (tau64.copy(), v64.copy())
            if e_blk is not None and ((e_blk - eps0 < 2e-9 and not aware_on and not pull_on)
                                      or (abs(e_blk - eps0) < 2e-9 and not (PRE and k == 0))):
                break       # the stage has stopped moving (below 0.000002 mHa per 1,000 iterations): near the end
                            # of a run the line search of every iteration then runs to its last step, 20 times the cost
            e_blk = eps0
            if a.log_every and done and done % a.log_every == 0:      # chem_own2: convergence curve of the refit
                log(f"  stage {k} iteration {done}: error_2 {1e3 * eps0:.5f} mHa")
            ref = (jnp.asarray(A0, DT), jnp.asarray(r0, DT), DT(eps0), DT((eps0 - e_base) - a.kappa * d_base))
            dref = (jnp.asarray(own, DT), jnp.asarray(th, DT), jnp.asarray(c0, DT), jnp.asarray(s0, DT),
                    jnp.asarray(ain0, DT), jnp.asarray(u0, DT), jnp.asarray(p0, DT), jnp.asarray(v64, DT))
            if aware_on:        # fresh rounding noise for this block: uniform in half a step, times EPS, free angles only
                nz_j = jnp.asarray((a.aware / scale) * rng_n.uniform(-0.5, 0.5, theta.shape) * (~frozen), DT)
                if done:        # the objective changed with the noise: L-BFGS starts afresh (30 iterations rebuild its memory)
                    st = opt.init(zero)
            else:
                nz_j = jnp.zeros(theta.shape, DT)
            if pull_on:         # the words the free angles would get now, and the strength of this block
                m_now = angles_and_mass(decode_np(th + K, edges), edges)[1]
                c_now = (nbw[:, None] ** 2) * m_now + 1e-30
                q_now = coded_round(th, c_now)[0]
                if hnorm is None:
                    e_all = reference(np.where(frozen, qv, q_now), v64 * sw)[2]
                    hnorm = (e_all * e_all - eps0 * eps0) / max(float(np.sum(c_now * (q_now - th) ** 2 * (~frozen))), 1e-300)
                    log(f"pull: raw energy of rounding every free angle now {1e6 * (e_all * e_all - eps0 * eps0):.0f} mHa^2 (error {1e3 * e_all:.1f} mHa)")
                nb_ = max(1, iters // a.recenter)
                mu0_ = a.pull_start if a.pull_start > 0 else a.pull / 1000.0
                mu = mu0_ * (a.pull / mu0_) ** (min(done // a.recenter, nb_ - 1) / max(nb_ - 1, 1)) if (PRE and k == 0) else a.pull
                pmask_ = ~frozen
                if a.tail_pull:     # r47: the small tails only, to "keep" / "pass on"
                    q_now = np.round(2.0 * th) / 2.0
                    pmask_ = (~frozen) & (np.abs(th - q_now) * scale <= a.tail_pull + 0.5)
                pw_j = jnp.asarray(mu * hnorm * c_now * pmask_ / (2.0 * eps0), DT)
                qt_j = jnp.asarray(q_now, DT)
                if done and not aware_on:
                    st = opt.init(zero)
                off_now = float(np.sqrt(np.sum(c_now * (q_now - th) ** 2 * pmask_) * hnorm))
                if a.tail_pull and done % (10 * a.recenter) == 0:
                    d_ = np.abs(th - q_now) * scale
                    log(f"  tail pull: words within {a.tail_pull} steps of a trivial word {int(pmask_.sum())}, of them under half a "
                        f"step {int((pmask_ & (d_ <= 0.5)).sum())}, under 0.1 step {int((pmask_ & (d_ <= 0.1)).sum())}")
                if done % (10 * a.recenter) == 0:
                    log(f"  pull {mu:.2e} at iteration {done}: error {1e3 * eps0:.3f} mHa, raw energy still to round {1e6 * off_now ** 2:.0f} mHa^2")
            else:
                pw_j, qt_j = jnp.zeros(theta.shape, DT), jnp.zeros(theta.shape, DT)
            for _ in range(max(1, min(a.recenter, iters - done) // a.chunk)):
                params, st = run(params, st, *ref, fm, wl, *dref, jnp.asarray(mst), nz_j, sc_dev, pw_j, qt_j)
                done += a.chunk
            # fold the changes into the float64 copies and move the origin to the new point; the L-BFGS memory is
            # made of differences and survives, only the last point it keeps has to move with the origin
            tau64 = tau64 + np.where(frozen, 0.0, np.asarray(params["tau"], np.float64)[mst])
            v64 = v64 + np.asarray(params["v"], np.float64)
            prev = optax.tree_utils.tree_get(st, "params")
            st = optax.tree_utils.tree_set(st, params={"tau": prev["tau"] - params["tau"], "v": prev["v"] - params["v"]})
            params = zero
        # chem_wr7: the LAST block of the refit is checked too (the guard above sees a block only when the next starts)
        if blk_start is not None and e_blk is not None and not aware_on and not pull_on:
            w7own = tau64 * sc
            w7th = np.where(frozen, qv, np.where(wall, np.clip(w7own, lo, hi - 1e-9), w7own))
            w7e = reference(w7th, v64 * sw)[2]
            if not np.isfinite(w7e) or w7e - e_blk > 1e-7:
                log(f"  stage {k}: the last refit block raised error_2 from {1e3 * e_blk:.4f} to {1e3 * w7e:.4f} mHa: "
                    f"undone (chem_wr7 guard)")
                tau64, v64 = blk_start
        own = tau64 * sc
        theta = np.where(frozen, qv, np.where(wall, np.clip(own, lo, hi - 1e-9), own))
        theta, nf, nn, ncl = regauge(theta, frozen, flip_ok, negat, edges, lo, hi)
        qv = np.where(frozen, theta, qv)
        untied = 0
        for m_, s_ in zip(pm_rows[tied], ps_rows[tied]):      # the gauge moves of the two rows must have been the same
            fr_ = ~frozen[m_]
            if not np.allclose(theta[m_][fr_], theta[s_][fr_], atol=1e-12):
                mst[s_] = s_
                untied += 1
        w = v64 * sw
        eps = reference(theta, w)[2]
        inc = eps * eps - e_base * e_base
        attributed += 1e6 * inc * stage_share
        total_inc += 1e6 * inc
        now_spent = spent + time.time() - t0
        trace.append((float(frozen.mean()), 1e3 * e_first, 1e3 * eps))
        log(f"stage {k:2d}: frozen {frozen.mean():.4f}  error_2 {1e3 * e_first:8.3f} -> {1e3 * eps:8.3f} mHa  crossings {nf} "
            f"negated frozen {nn} walled {int(wall.sum())} clipped {ncl} put back {unstuck}  ({now_spent:.0f} s)  damage "
            f"{1e6 * inc:+.2f} mHa^2 this stage, share of the top 1/2/4/6% vectors "
            f"{' '.join(f'{x:.2f}' for x in stage_share)}" + (f"; pairs untied {untied}, still tied {int((mst[ps_rows] == pm_rows).sum())}" if len(pm_rows) else "")
            + pair_note + (f", whole raw kick {1e6 * (e_first * e_first - e_base * e_base):.1f}" if pair_note else ""))
        if HOLD and k == stage0:
            log(f"hold: {done} iterations at frozen {frozen.mean():.4f}: error_2 {1e3 * e_first:.4f} -> {1e3 * eps:.4f} mHa, "
                f"error^2 {1e6 * (eps * eps - e_first * e_first):+.2f} mHa^2 (the state is not written)")
            return
        extra = {} if RES is None else {"res": RES}
        np.savez(state_path, theta=theta, w=w, frozen=frozen, qv=qv, stage=k + 1, spent=now_spent, trace=np.array(trace), mst=mst, **extra)
        for ks_ in KEEP:        # r47
            if (sched[k - 1] if k else 0.0) < ks_ <= frac:
                np.savez(os.path.join(a.out, f"state_{ks_:.4f}.npz"), theta=theta, w=w, frozen=frozen, qv=qv, stage=k + 1,
                         spent=now_spent, trace=np.array(trace), mst=mst, **extra)
        if frozen.mean() >= a.stop_at and not frozen.all():
            log(f"stopped at frozen {frozen.mean():.4f} (--stop-at): a probe, no words written")
            return
        for g in filter(None, a.gate.split(",")):
            gshare, gval = (float(x) for x in g.split(":"))
            prev = sched[k - 1] if k else 0.0
            if prev < gshare <= frac:
                if 1e3 * eps > gval:
                    log(f"GATE FAILED at frozen {frac:.4f}: error_2 {1e3 * eps:.3f} mHa is above {gval} (--gate): stopped")
                    sys.exit(5)
                log(f"gate at frozen {frac:.4f}: error_2 {1e3 * eps:.3f} mHa <= {gval}: passed")
        if (time.time() - t0) / 60.0 > a.max_minutes and k + 1 < len(sched):
            log(f"stopped after stage {k} by --max-minutes; rerun to resume")
            return
    assert frozen.all()
    if a.fine:
        log(f"all frozen on the mixed grid (--fine: a bound, no words written): error_2 {1e3 * reference(theta, w)[2]:.4f} mHa")
        return
    if total_inc > 0:
        log("damage attributed to the 1 / 2 / 4 / 6 % strongest vectors (stage damage x their share of the impact frozen "
            "in the stage; stages run in this process only): " + " ".join(f"{x / total_inc:.3f}" for x in attributed)
            + f" of {total_inc:.1f} mHa^2")
    finish(np.round(qv * scale).astype(np.int64), w, trace, spent + time.time() - t0)


if __name__ == "__main__":
    main()
