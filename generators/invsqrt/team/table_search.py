#!/usr/bin/env python3
# New file added by the team (Maxime Bendavid), 2026-11-18, for idea INV-A02 (plan step S2.7).
# The integer phase model below is a parameterised re-statement of Syed Asad Ali's `whole_model.py`
# (src/work/c10/narrow/native/candidate01/, Apache-2.0); see ../NOTICE. Syed Asad Ali is the original builder.
"""Coefficient-table search for the invsqrt generator (64-row variant), certified on every input.

The coefficient for root v (12 bits, v = Q - 4096) is built by the circuit as

    D(v) = A[a] + sum_{i < L} (1 - 2 v_i) * ((S[a] >> (L - i)) + q_i),   a = v >> L,  C = 2 D,

with q_i the circuit's rounding carry (k = L - 1 - i: q = v_i (S_{k-1} xor S_k) xor S_k for k > 0, and
q = (1 - v_i) S_0 for k = 0). L = 5 for the 128-row table (K = 7) and L = 6 for the 64-row table (K = 6).
This formula reproduces the frozen K = 7 model's 4,096 coefficients exactly (checked by `--check-k7`).

For each x in [1, 3,133,452], `phase_words` replays the generator's exact integer model (normalisation, integer
root, seed, 10 coefficient rows, 27 signed-digit divider rows), and the margin is the verifier's tolerance minus
|phase error| in 2^-64 turns, with the power-of-two ("base") angles. A table passes when every margin is
positive and no physical register overflows. Rows are independent for the base angles, so each row is searched
on its own inputs only; the 34 angles are refitted afterwards (step S2.10) on the whole domain.

    python table_search.py --check-k7                      # regression: K = 7 model margin must be 31,966,138
    python table_search.py --search --out model_k6.json    # K = 6 search (12 processes)

Complexity: preparing the domain is O(N) with N = 3,133,452; evaluating one (A, S) candidate on a row costs
O(N / 64 x 37) vectorized operations (~10 ms); the search evaluates a window of candidates per row, rows in
parallel. Memory ~ 0.5 GB.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from multiprocessing import get_context

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(HERE), "src")
NARROW = os.path.join(SRC, "work", "c10", "narrow")
K7_MODEL = os.path.join(NARROW, "research", "recompute_root", "models", "phasefit_z1_K7_G4_m37_virtual_quadS2_36_12c0.json")
K7_SHA = "24fff823faa911afe93d02afc3ded0bbd156a210a85f1265c989a087f8fd1219"
MAX_X = 3133452
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
VERIFIER = os.path.join(REPO, "google-verifiable-quantum-circuit-compilation", "verifiers", "invsqrt")


# ---------------------------------------------------------------- the domain (computed once)
def domain():
    """Per-input quantities that do not depend on the table, plus the verifier's targets and tolerances."""
    sys.path.insert(0, VERIFIER)
    import verify_submission_invsqrt as V
    spec = V.ProblemSpecInvSqrt.for_n_p(10)
    assert spec.max_x == MAX_X
    x = np.arange(1, MAX_X + 1, dtype=np.int64)
    leading = np.searchsorted(np.left_shift(np.ones(22, dtype=np.int64), np.arange(22)), x, side="right") - 1
    exponent = (21 - leading) // 2
    norm = x << (2 * exponent)
    Q = np.sqrt(norm << 4).astype(np.int64)
    assert np.all(Q * Q <= norm << 4) and np.all((Q + 1) * (Q + 1) > norm << 4)
    v = Q - 4096
    d = 2 * Q + 1
    h = ((norm << 6) - d * d) // 64
    target = np.empty(MAX_X, dtype=np.uint64)
    tol = np.empty(MAX_X, dtype=np.int64)
    for s in range(0, MAX_X, 1 << 18):
        t, tl = V.reference_values(spec, x[s:s + (1 << 18)].astype(np.uint64))
        target[s:s + len(t)] = t
        tol[s:s + len(t)] = np.asarray(tl, dtype=np.int64)
    return dict(x=x, exponent=exponent, v=v, d=d, q=h + 512, target=target, tol=tol)


# ---------------------------------------------------------------- the integer model
def coefficients(A, S, L):
    """D(v) for all 4,096 roots from row words A, S (arrays of length 4096 >> L)."""
    v = np.arange(4096, dtype=np.int64)
    a = v >> L
    Sv = np.asarray(S, dtype=np.int64)[a]
    D = np.asarray(A, dtype=np.int64)[a].copy()
    for i in range(L):
        k = L - 1 - i
        s = (v >> i) & 1
        if k:
            q = (s & (((Sv >> (k - 1)) ^ (Sv >> k)) & 1)) ^ ((Sv >> k) & 1)
        else:
            q = (1 - s) & (Sv & 1)
        D = D + (1 - 2 * s) * ((Sv >> (k + 1)) + q)
    return D


def seeds(C):
    v = np.arange(4096, dtype=np.int64)
    rootclass = (v * (v + 1) // 2) & 7
    seed = -np.ones(4096, dtype=np.int64)
    for i in range(3):
        seed += (2 * ((rootclass >> i) & 1) - 1) * ((C + (1 << (11 - i))) >> (12 - i))
    return seed


def phase_words(C, U0, dom, sel=None):
    """Exact model phase words (2^-37 turns) and an overflow flag per input. Mirrors whole_model.phase_units."""
    v, d, q, exponent = dom["v"], dom["d"], dom["q"], dom["exponent"]
    if sel is not None:
        v, d, q, exponent = v[sel], d[sel], q[sel], exponent[sel]
    Cx, Ux = C[v], U0[v]
    F = np.zeros(len(v), dtype=np.int64)
    for i in range(10):
        shift = 8 - i
        term = (Cx + (1 << (shift - 1))) >> shift if shift > 0 else Cx if shift == 0 else 2 * Cx + 1
        F += (2 * ((q >> i) & 1) - 1) * term
    U = Ux - F // 2
    overflow = (U < -16384) | (U >= 16384)
    product = U & 32767
    r = 4096 - d
    phase = np.full(len(v), -((1 << 37) - (1 << 29)), dtype=np.int64)
    for j in range(9, 36):
        sign = (r < 0).astype(np.int64)
        phase -= sign << (37 - j)
        if j <= 23:
            bit = (product >> (23 - j)) & 1
            if j == 9:
                bit = -bit
            raw = 2 * r + (2 * sign - 1) * d + bit
        elif j == 24:
            raw = 2 * r + (2 * sign - 1) * d
        else:
            k = j - 26 if j >= 28 else 0
            V = d // 2
            Aa = ((((V >> (k + 1)) + ((V >> k) & 1)) << (k + 1)) - (V & 1)) if k else V
            u = (r - 1) // 2
            raw = 2 * (2 * u + sign + (2 * sign - 1) * Aa) + 1
        nextu = raw // 2
        overflow |= (nextu < -8192) | (nextu >= 8192)
        nextu = ((nextu + 8192) & 16383) - 8192
        r = 2 * nextu + (raw & 1)
    phase -= (r < 0).astype(np.int64) << 1
    p = (phase - 1) & ((1 << 37) - 1)
    word = ((p << exponent) + (p >> (27 - exponent))) & ((1 << 37) - 1)
    return word, overflow


def margins(C, U0, dom, sel=None):
    """tol - |error| in 2^-64 turns with the base angles; overflowed inputs get a margin of -2^62."""
    word, overflow = phase_words(C, U0, dom, sel)
    target = dom["target"] if sel is None else dom["target"][sel]
    tol = dom["tol"] if sel is None else dom["tol"][sel]
    err = ((word.astype(np.uint64) << np.uint64(27)) - target).view(np.int64)
    m = tol - np.abs(err)
    m[overflow] = -(1 << 62)
    return m


# ---------------------------------------------------------------- K = 7 regression
def check_k7(dom) -> dict:
    raw = open(K7_MODEL, "rb").read()
    assert hashlib.sha256(raw).hexdigest() == K7_SHA
    model = json.loads(raw)
    D = coefficients(model["A"], model["S"], 5)
    assert np.array_equal(D, np.asarray(model["D"])), "coefficient formula"
    C = 2 * D
    U0 = seeds(C)
    assert np.array_equal(U0, np.asarray(model["U0"])), "seed formula"
    m = margins(C, U0, dom)
    return dict(min_margin=int(m.min()), failures=int((m <= 0).sum()), expected=model["minimum_margin_units"])


# ---------------------------------------------------------------- K = 6 search
_DOM = None


def _row_inputs(a, L):
    return np.flatnonzero((_DOM["v"] >> L) == a)


def _eval_row(a, L, A, S, sel):
    """Minimum margin of row a for words (A, S): build that row's coefficients only."""
    rows = 4096 >> L
    Aw = np.zeros(rows, dtype=np.int64)
    Sw = np.zeros(rows, dtype=np.int64)
    Aw[a], Sw[a] = A, S
    C = 2 * coefficients(Aw, Sw, L)
    U0 = seeds(C)
    return int(margins(C, U0, _DOM, sel).min())


def search_row(args):
    a, L, A0, S0, s_bits, a_bits, radius_a, radius_s = args
    sel = _row_inputs(a, L)
    best = (_eval_row(a, L, A0, S0, sel), A0, S0)
    # Coarse-to-fine coordinate search on (A, S) inside the register limits.
    for step_a, step_s in ((8, 4), (4, 2), (2, 1), (1, 1)):
        improved = True
        while improved:
            improved = False
            _, Ab, Sb = best
            for dA in range(-radius_a, radius_a + 1, step_a):
                for dS in range(-radius_s, radius_s + 1, step_s):
                    A, S = Ab + dA, Sb + dS
                    if not (0 <= A < (1 << a_bits) and 0 <= S < (1 << s_bits)):
                        continue
                    m = _eval_row(a, L, A, S, sel)
                    if m > best[0]:
                        best, improved = (m, A, S), True
            radius_a, radius_s = max(radius_a // 2, 2), max(radius_s // 2, 2)
    return dict(a=a, A=int(best[1]), S=int(best[2]), min_margin=int(best[0]), inputs=int(len(sel)))


def initial_k6(model_k7):
    """Start from the K = 7 coefficients: A = value at the row centre, S = 32 x average decrement per root."""
    D7 = np.asarray(model_k7["D"], dtype=np.int64)
    A0, S0 = [], []
    for a in range(64):
        seg = D7[a * 64:(a + 1) * 64]
        A0.append(int(round((seg[31] + seg[32]) / 2)))
        S0.append(int(round(32 * (seg[0] - seg[63]) / 63)))
    return A0, S0


# ---------------------------------------------------------------- angle refit (feasibility; S2.10 does the final one)
UNIT = 1 << 27
COEFF = np.array([(1 << i) + (1 << 27 if i == 0 else 0) for i in range(27)] + [-(1 << 27)]
                 + [1 << (27 + j) for j in range(1, 7)], dtype=np.int64)


def feature_masks(word, exponent):
    """34 feature bits per input, as in Syed Asad Ali's research/postgate_angle_refit/refit.py."""
    z = word & ((1 << 27) - 1)
    packed = z.astype(np.uint64)
    packed |= (exponent == 0).astype(np.uint64) << np.uint64(27)
    for j in range(1, 7):
        packed |= (((z >> j) & 1) * (exponent > j)).astype(np.uint64) << np.uint64(27 + j)
    return packed


def _features(mask):
    return ((mask[:, None] >> np.arange(34, dtype=np.uint64)[None, :]) & np.uint64(1)).astype(np.float64)


def refit(model, dom, rounds_max=24):
    """Active-set LP over 34 angle deltas maximising the minimum margin, then an exact screen on all inputs.

    Same method as refit.py (HiGHS LP on an active set grown with the worst inputs, deltas bounded to +-1 unit
    of 2^-37 turns, quantised to 2^-64 turns). Returns the exact minimum margin, failures and the deltas."""
    from scipy.optimize import linprog
    C = np.asarray(model["C"], dtype=np.int64)
    U0 = np.asarray(model["U0"], dtype=np.int64)
    word, overflow = phase_words(C, U0, dom)
    assert not overflow.any(), "physical overflow: no angle can repair it"
    errors = ((word.astype(np.uint64) << np.uint64(27)) - dom["target"]).view(np.int64)
    tols = dom["tol"]
    masks = feature_masks(word, dom["exponent"])
    margin0 = (tols - np.abs(errors)).astype(np.float64) / UNIT
    active = set(np.argpartition(margin0, 512)[:512].tolist())
    active.update(np.linspace(0, len(tols) - 1, 512, dtype=np.int64).tolist())
    best, best_min, log = np.zeros(34), float(margin0.min()), []

    def screen(delta):
        out = np.empty(len(tols))
        for s in range(0, len(tols), 1 << 16):
            pred = _features(masks[s:s + (1 << 16)]) @ delta
            out[s:s + len(pred)] = tols[s:s + len(pred)] / UNIT - np.abs(errors[s:s + len(pred)] / UNIT + pred)
        return out

    for it in range(rounds_max):
        ix = np.array(sorted(active), dtype=np.int64)
        F = _features(masks[ix])
        ones = np.ones((len(ix), 1))
        A = np.concatenate((np.concatenate((F, ones), 1), np.concatenate((-F, ones), 1)), 0)
        b = np.concatenate(((tols[ix] - errors[ix]) / UNIT, (tols[ix] + errors[ix]) / UNIT))
        obj = np.zeros(35)
        obj[-1] = -1
        res = linprog(obj, A_ub=A, b_ub=b, bounds=[(-1.0, 1.0)] * 34 + [(0, None)], method="highs")
        if not res.success:
            log.append(dict(iteration=it, status=str(res.message)))
            break
        delta = res.x[:34]
        sc = screen(delta)
        log.append(dict(iteration=it, active=len(active), lp_margin37=float(res.x[-1]), full_margin37=float(sc.min())))
        if sc.min() > best_min:
            best, best_min = delta.copy(), float(sc.min())
        if res.x[-1] - sc.min() <= 2e-8:
            break
        active.update(np.argpartition(sc, 256)[:256].tolist())
    quant = np.rint(best * UNIT).astype(np.int64)
    new_err = errors.copy()
    for s in range(0, len(tols), 1 << 16):
        new_err[s:s + (1 << 16)] += _features(masks[s:s + (1 << 16)]).astype(np.int64) @ quant
    m = tols - np.abs(new_err)
    return dict(min_margin_units=int(m.min()), failures=int((m <= 0).sum()),
                worst_error_over_tolerance=float((np.abs(new_err) / tols).max()),
                delta_turn_units64=quant.tolist(), rounds=log)


def main() -> int:
    global _DOM
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check-k7", action="store_true")
    ap.add_argument("--search", action="store_true")
    ap.add_argument("--s-bits", type=int, default=7)
    ap.add_argument("--a-bits", type=int, default=13)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--out")
    ap.add_argument("--refit", help="model JSON to refit (angle feasibility)")
    a = ap.parse_args()
    t0 = time.perf_counter()
    _DOM = domain()
    print(f"domain ready in {time.perf_counter() - t0:.1f} s", flush=True)
    if a.check_k7:
        r = check_k7(_DOM)
        print(json.dumps(r), f"[{time.perf_counter() - t0:.1f} s]")
        if r["min_margin"] != r["expected"] or r["failures"]:
            return 1
    if a.refit:
        r = refit(json.load(open(a.refit)), _DOM)
        print(json.dumps({k: v for k, v in r.items() if k != "delta_turn_units64"}), f"[{time.perf_counter() - t0:.1f} s]")
        if a.out:
            with open(a.out, "w") as f:
                json.dump(r, f)
    if a.search:
        model_k7 = json.load(open(K7_MODEL))
        A0, S0 = initial_k6(model_k7)
        tasks = [(r, 6, A0[r], S0[r], a.s_bits, a.a_bits, 32, 16) for r in range(64)]
        with get_context("fork").Pool(a.workers) as pool:
            rows = pool.map(search_row, tasks)
        A = [r["A"] for r in rows]
        S = [r["S"] for r in rows]
        D = coefficients(A, S, 6)
        C = 2 * D
        U0 = seeds(C)
        m = margins(C, U0, _DOM)
        out = dict(name="team_K6_rows64", K=6, low_bits=6, A_bits=a.a_bits, S_bits=a.s_bits, row_count=64,
                   status="PASS_EXACT_BASE_ANGLES" if m.min() > 0 else "FAIL", minimum_margin_units=int(m.min()),
                   failures=int((m <= 0).sum()), A=A, S=S, D=D.tolist(), C=C.tolist(), U0=U0.tolist(),
                   row_details=rows, worst_row=min(rows, key=lambda r: r["min_margin"]),
                   seconds=round(time.perf_counter() - t0, 1))
        print(json.dumps({k: v for k, v in out.items() if k not in ("A", "S", "D", "C", "U0", "row_details")}))
        if a.out:
            with open(a.out, "w") as f:
                json.dump(out, f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
