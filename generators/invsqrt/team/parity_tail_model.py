#!/usr/bin/env python3
# New file added by the team (Maxime Bendavid), 2026-11-19, plan step S2.19 (experiment E1 of idea INV-B10).
# Builds on team/table_search.py (our re-statement of Syed Asad Ali's integer model, Apache-2.0; see ../NOTICE).
"""E1: error budget of a history-free ("parity") divider tail, in the exact integer phase model.

Incumbent divider rows (j = 9..35 in `table_search.phase_words`) are non-restoring: the digit is s = [r < 0] and
r' = 2r -/+ d (+ a fed bit in rows 9..23). A parity row (INV-B10) uses p = [r >= 0] and
r' = 2r - d if p else 2r + d + 1 on the interval [-d-1, d] (d = 2Q + 1 is odd), so p = r' mod 2 and the row can be
undone without a stored digit. Its digits are the binary quotient of (r0 + d + 1) / (2d + 1) instead of r0 / d:
written in the incumbent's place (s = 1 - p, same bit positions of the phase word), the tail differs from the
incumbent by a factor (1 + 1/(2d)) on its own value, which the 34-angle refit can only absorb on average.

For each tail start j_tail, rows j >= j_tail (and the final rounding digit) become parity rows on the full divisor d
(version A of the S2.18 report); the head rows are unchanged, including the incumbent's own feed-free rows below
j_tail. The table is the frozen K = 7 model (v12), so the reference is the incumbent model under the same LP refit.

    python parity_tail_model.py [--tails 24 26 28 30]

Cost: one exact model evaluation is O(27 N) vectorised (N = 3,133,452, about 1 s); one refit about 5-10 s.
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np

import table_search as T


def make_phase_words(j_tail: int):
    """phase_words with parity rows for j >= j_tail (and for the final rounding digit at weight 2)."""

    def phase_words(C, U0, dom, sel=None):
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
            if j >= j_tail:  # parity row on the full divisor: digit p = [r >= 0], s = 1 - p
                p = (r >= 0).astype(np.int64)
                phase -= (1 - p) << (37 - j)
                r = np.where(p == 1, 2 * r - d, 2 * r + d + 1)
                overflow |= (r < -d - 1) | (r > d) | (r < -16384) | (r > 16383)
                continue
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
        if j_tail <= 35:  # final rounding digit as one more parity digit
            phase -= (1 - (r >= 0).astype(np.int64)) << 1
        else:
            phase -= (r < 0).astype(np.int64) << 1
        pw = (phase - 1) & ((1 << 37) - 1)
        word = ((pw << exponent) + (pw >> (27 - exponent))) & ((1 << 37) - 1)
        return word, overflow

    return phase_words


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tails", type=int, nargs="+", default=[24, 26, 28, 30])
    ap.add_argument("--json")
    a = ap.parse_args()
    t0 = time.perf_counter()
    dom = T.domain()
    model = json.load(open(T.K7_MODEL))
    original = T.phase_words
    rows = []
    # Reference: the incumbent model, base angles and refit.
    m0 = T.margins(np.asarray(model["C"]), np.asarray(model["U0"]), dom)
    ref = T.refit(model, dom)
    rows.append(dict(j_tail="none (incumbent)", freed_rows=0, base_min_margin=int(m0.min()), base_failures=int((m0 <= 0).sum()),
                     refit_min_margin=ref["min_margin_units"], refit_failures=ref["failures"],
                     worst_ratio=round(ref["worst_error_over_tolerance"], 5)))
    print(json.dumps(rows[-1]), f"[{time.perf_counter() - t0:.1f} s]", flush=True)
    for jt in a.tails:
        T.phase_words = make_phase_words(jt)
        try:
            m = T.margins(np.asarray(model["C"]), np.asarray(model["U0"]), dom)
            word, overflow = T.phase_words(np.asarray(model["C"]), np.asarray(model["U0"]), dom)
            row = dict(j_tail=jt, freed_rows=36 - jt, base_min_margin=int(m.min()), base_failures=int((m <= 0).sum()),
                       overflows=int(overflow.sum()))
            if not overflow.any():
                r = T.refit(model, dom)
                row.update(refit_min_margin=r["min_margin_units"], refit_failures=r["failures"],
                           worst_ratio=round(r["worst_error_over_tolerance"], 5),
                           margin_vs_reference=round(r["min_margin_units"] / ref["min_margin_units"], 4))
        finally:
            T.phase_words = original
        rows.append(row)
        print(json.dumps(row), f"[{time.perf_counter() - t0:.1f} s]", flush=True)
    if a.json:
        json.dump(rows, open(a.json, "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
