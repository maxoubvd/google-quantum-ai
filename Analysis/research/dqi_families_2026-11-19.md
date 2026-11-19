# dqi: divergent decoder families on paper, with red team (S2.21, GEN-B11)

Written 2026-11-19 for milestone M2.C. Scope: qualitatively different decoder families or stage-level
redesigns for the `dqi` task, each priced on paper (qubit peak register by register, Toffolis forward +
uncompute, failure rate, cost `R x Q^2 x T`), red-teamed against `experiments/NULLS.md` (dqi, 75 + 24 rows),
`Analysis/research/dqi.md` section 4 and the public-frontier notes, and ranked. Nothing here was built or run;
all numbers are paper estimates from exact formulas or from the v12 stage ledger
(`generators/dqi/src/controller_work/after_v3/integration_v6/SOURCE_MODEL_LEDGER.json`).

## Summary

1. **Reference:** public frontier `e39b1396`, 890 Q, 365,893 mean Toffolis, cost 3.4526e11 (R = 1.19126).
   The M2.C target (1.34x on dqi) is cost <= 2.577e11, i.e. T <= ~274k at 890 Q (<= 5,462 GF(2^12)
   multiplication-equivalents instead of ~7,668).
2. **Eight families priced** (F1-F8): explicit-convergent extractor, re-based extractor arithmetic,
   convert-free per-position Forney by a Boolean quotient phase, reciprocal (Horiguchi-Koetter / Bezout-reuse)
   convert, decode-re-encode reverse, narrow-width line, cross-stage retention with measurement-based
   uncomputation, quartic-node CRT.
3. **No family reaches the 1.34x target on paper.** Seven are structural nulls with a number attached; one
   (F3) survives only conditionally, as a cheap mathematical screen with a large ceiling.
4. **Why (six laws, section 1):** a phase oracle costs the multiplicative complexity (MC) of its Boolean
   function; an in-place map and its inverse cost the same; one O(n)-word temporary (35 words) more than
   doubles the cost through Q^2; each extra Toffoli per scan position costs +1.1%; the implicit Thiele chain
   needs ~595 GF(2^24) products per direction at 1x width; the weight law forbids truncating the weight.
5. **F3 (convert-free Forney, new DQI-B13)** removes the convert/unconvert pair and the projective finish
   (-122k static T, 33%) by computing at each of the 4,025 scan positions the phase
   `[sigma(X)=0] * Tr(nu_X * Omega(X)/O(X)^2)` directly from (E, O, Omega), which fit the 840-qubit core.
   It pays only if that quotient phase costs <= ~30 Toffolis above the 11-T zero test. Our best construction
   (GF(64)/GF(8) tower) needs ~69; the degree bound only proves >= 11. Ceiling -20% (cost 2.75e11), P ~0.03.
6. **F2 (arithmetic) and F6 (width)** leave only small residues already in the register: the GF(2^24)
   rank-below-99 lottery (DQI-B08, <= -1%) and one qubit at 888 by norm-based inverse helpers (<= -0.4%).
7. **Hard nulls recorded:** explicit convergents (+840 Q at slot 35), subfield-normalized pivots (168 T
   instead of 99 per product), shared-pivot CRT batching (break-even at 3 simultaneous heads, +24 Q each),
   reciprocal convert (in-place multiply = in-place divide), decode-re-encode (in-place encode is a Hankel
   solve, out-of-place needs +840 Q), held products across the mirror (needs >= 267 products per held word),
   quartic nodes (-52k T of queries eaten by +48 Q and GF(2^48) inverses).
8. **Recommendation:** (1) run the F3 MC screen (exact/heuristic MC of `Tr(a/b)` over GF(2^n), n = 3..6, then
   a tower-plus-search bound at n = 12; 1 day; kills or opens the only family with a -20% ceiling);
   (2) fold F6-lite (888 Q via 12-bit-norm inverse helpers) into the A-track (DQI-A09). No family deserves a
   full-scale January build on current evidence.
9. **Plan consequence:** M2.C go/no-go ("1-2 families priced below the frontier") is **not met**, except
   conditionally through F3. The leaders' gap is more likely on chemistry (or on a dqi idea outside the
   in-place CRT/Thiele class that we cannot see); January family effort should go to chemistry first.

## 0. Cost model used throughout

Cost = R x Q^2 x T with R = 1.19126 at zero failures. At the frontier (Q = 890, T = 365,893):

| Lever | Task cost | Toffoli equivalent |
|---|---:|---:|
| -1,000 mean Toffolis | -0.273% | 1,000 |
| -1 qubit | -0.224% | ~820 |
| +12 qubits (one GF(2^12) word) | +2.72% | ~+9.9k |
| +24 qubits (one GF(2^24) word) | +5.47% | ~+20.0k |
| +420 qubits (35 words: one polynomial) | +117% | — |
| +1 Toffoli at each of the 4,025 scan positions | +1.10% | 4,025 |
| +1 expected failed shot per 102,400-shot draw | ~+0.46% | ~1,680 |

Stage anatomy (v12 static, used as the proxy for the frontier's shares): key = extractor + synthesis 216,928
(forward 101,087, reverse 115,841), of which 1,190 GF(2^24) query MACs at 99 T = 117,810 and synthesis
65,076; convert (Forney preimage) 114,423 (57,034 + ~57,389); scan 44,286 (one pass, 11 T x ~4,025
positions); projective finish 7,856. Width 890 = 840 core + 24-bit carry + up to 26 helpers/flags.

## 1. Structural laws used by the red team

These are the facts that decided most verdicts. L1-L3 are standard; L4-L6 are arithmetic on the task.

- **L1 (phase oracle = multiplicative complexity).** A phase `(-1)^f(x)` costs at most MC(f) Toffolis: compute
  f with an XOR-AND circuit (one temporary AND per gate), kick the phase, and uncompute every AND by
  X-measurement plus a CZ fix-up (0 Toffolis). A phase of algebraic degree d needs at least d - 2
  Toffoli-likes. The scan cell is at this floor: `[w = 0] * l(R)` has degree 13 and costs 11.
- **L2 (in-place symmetry).** An in-place bijection and its inverse have the same optimal cost (C and its
  adjoint). In-place multiplication by a unit of `A_sigma = GF(2^12)[x]/sigma` costs exactly what in-place
  division by it costs. Out-of-place versions are cheap (evaluation-based products at 42 T per point, MBU
  uncompute at the bilinear-form rank) but need a fresh polynomial register (L3).
- **L3 (width price).** Every fast algorithm (half-gcd, product trees, Newton inversion, explicit block
  matrices, out-of-place products mod sigma) needs at least one temporary polynomial of ~35 words (420 Q),
  which multiplies the cost by >= 2.17. The core is forced in place with O(1) words of scratch (NR-0265).
- **L4 (per-position budget).** Anything executed per scan position is paid ~4,025 times: a family that
  deletes D static Toffolis elsewhere may add at most D / 4,025 Toffolis per position. Classical skipping of a
  position by a quantum predicate is impossible (no qubit-to-bit copy, varianceofx B1e).
- **L5 (implicit chain floor).** At 1x width the extractor keeps one projective head per unprocessed node and
  evaluates each new node through all earlier pivots: sum (i - 1) = 595 GF(2^24) products per direction.
  Evaluating an explicit convergent at a public node would be free (CNOT), but explicit convergents need the
  current and the previous pair, about 2x the information (the Proos-Zalka factor).
- **L6 (failure budget).** One failed shot per draw costs ~0.46%; a class above ~1e-5 per shot is
  unaffordable. P(weight = 35) = 2.77% and P(weight <= 17) = 7.7e-7: the decoder must solve full weight 35,
  and only branches below ~1e-5 can be dropped (A07 territory, not a family).

## 2. Comparison table

Cost ratio = (Q/890)^2 x T/365,893 at equal R. "Best case" is the most favourable assumption that survives
the red team; P is the probability that a full build reaches at least break-even.

| Family | Stages replaced | Q peak | T (or #mult) best case | Cost best case | vs 3.4526e11 | P | Effort to full build | Red-team verdict |
|---|---|---:|---:|---:|---:|---:|---|---|
| F1 Explicit / blocked convergents ("free queries") | extractor queries + synthesis | 890 + 840 at slot 35 (blocks of 12: +624) | ~-50k T if width were free | ~1.13e12 | +226% | 0.01 | n/a | **Dead** (L3, L5): width |
| F2 Re-based extractor MAC arithmetic | 1,190 GF(2^24) MACs | 890 (normalization), 962 (batch m = 4) | rank 96: -3.6k | 3.418e11 | -1.0% | 0.10 | 40 h (rank search, = DQI-B08) | **Dead** except the rank lottery |
| F3 Convert-free Forney by Boolean quotient phase | convert + unconvert + finish (122k) | 890 (core = sigma, Omega; scan not at peak) | 4,025 x C_q - 122k: C_q = 12 -> 292k; C_q = 69 (our tower) -> 521k | 2.75e11 (C_q = 12) / 4.92e11 (C_q = 69) | -20% / +42% | 0.03 | screen 8 h; build ~120 h | **Conditional**: needs MC(Tr(a/b)) <= ~30 at n = 12 |
| F4 Reciprocal convert (Horiguchi-Koetter, Bezout reuse) | convert | 890 | >= +40k (per-position inverse-trace) | >= 3.83e11 | >= +11% | 0.01 | n/a | **Dead** (L2): dominated by F3 |
| F5 Decode-re-encode reverse | reverse key + unconvert | 890 in place / ~1,730 out of place | in place: no change; out of place: ~-60k at +840 Q | ~1.1e12 | > +200% | 0.00 | n/a | **Dead**: in-place encode = Hankel solve |
| F6 Narrow-width line (12-bit carry, norm inverses) | carry bank, inverse helpers | 888 (lite) / 878 (full, if T flat) | lite: +0.4 to +1.5k at -2 Q | 3.44e11 (lite) | -0.05 to -0.35% | 0.35 (lite) | 20 h (lite) | **Full dead** (+118k T for the subfield carry); lite survives as an A-item |
| F7 Cross-stage retention with MBU (99 -> 24 T on the reverse copy) | reverse recomputations | 890 + 24 per held word | -75 T per held product | > 3.6e11 | >= +5% | 0.01 | n/a | **Dead**: needs >= 267 products per held word |
| F8 Quartic-node CRT (GF(2^48) heads) | extractor | ~938 | queries -51.8k, overheads +10k | ~3.40e11 to 3.6e11 | -1.6% to +5% | 0.03 | 150 h | **Dead / marginal** (cubic precedent NR-0280: +48k) |

## 3. Families

Each family lists: (1) datapath, (2) qubit peak, (3) Toffolis and failure rate, (4) cost, (5) red team,
(6) scaled prototype, (7) P(success) and effort.

### F1. Explicit or blocked convergents ("free queries")

1. **Datapath.** Replaces the 1,190 query MACs and the synthesis. Keep the convergent polynomials of the
   continued fraction explicitly, in GF(2^12) coefficient form, during extraction: (A_k, B_k) and
   (A_{k-1}, B_{k-1}). Slot k+1 evaluates the four polynomials at the public node beta_{k+1} (a CNOT map),
   forms the new pivot from the head and these four GF(2^24) values (2 products, 1 inverse, 1 product,
   ~660 T), then updates A_{k+1} = P_{k+1}(x) A_k + q_{k+1}(x) A_{k-1} and the same for B (P linear with
   quantum coefficients, q a public quadratic). At the end (A_35, B_35) = (Omega, sigma): no synthesis.
   Blocked hybrid: materialize the 2x2 matrix of a block of b pivots only while it is applied to every later
   head, then unmaterialize it.
2. **Qubit peak.** At slot k: four convergents of degree ~k (4k + 4 words) plus 35 - k unprocessed heads
   (2 words each) = 70 + 2k + 4 words; at k = 35, 144 words = 1,728 qubits, plus carry and helpers:
   Q ~ 1,730 (L5: the Proos-Zalka factor 2). Blocks: 4(b + 1) words of scratch, 624 qubits for b = 12, and a
   right-looking application needs every later head normalized after each block.
3. **Toffolis.** Explicit updates by evaluation-based products: 2 x sum_{k=1..35} (k + 1) = 1,330 GF(2^12)
   products = 55.9k per direction, plus 35 x ~660 = 23.1k for pivots: ~79k per direction against ~109k today
   (queries 58.9k + synthesis 32.5k + per-slot overheads ~18k), so about -30k per direction, -60k for the pair,
   if width were free. Blocks: applying a materialized block to one head costs ~560 T (2 GF(2^24) products,
   normalization) against b x 99 for the chain: break-even at b ~ 6; for b = 12 it saves ~628 T per
   (block, head) but materializing costs ~12k per block (make + unmake): +2.4k net for the first block of 12,
   a loss for the later ones. Failure rate unchanged.
4. **Cost.** Q 1,730, T ~306k: ratio 3.78 x 0.86 = 3.26, cost ~1.13e12 (+226%).
5. **Red team.** L3 and L5: the previous convergent cannot be dropped, and recovering it from (A_k, B_k) is a
   Euclid of cost O(k^2). Same obstacle as DQI-B09 (written there as "48k bits at slot k"); NR-0045 (block
   continued fractions b in {2, 3, 4, 6}: marginal) is the small-block end of the same trade. **Dead.**
   Reopen only for a representation that stores one convergent pair plus O(1) words and recovers the
   previous pair with O(k) products.
6. **Prototype.** None. If reopened: RS over GF(2^6) (m = 63, 8 syndromes, 4 quadratic nodes), exhaustive.
7. **P = 0.01**; no build planned.

### F2. Re-based extractor arithmetic (normalized pivots, shared-pivot batching, rank, inversion-free pivots)

1. **Datapath.** Keeps the in-place Thiele chain and changes the arithmetic of its query MAC
   `d + P_l(alpha_i) n`, where P_l = a_l + b_l x has GF(2^12) coefficients and n is a GF(2^24) head:
   (a) **monic pivots** by an equivalence transformation of the continued fraction (P'_l = P_l / b_l), so
   that P'_l(alpha_i) n = a'_l n + (alpha_i n) is one GF(2^12) x GF(2^24) product (84 T) plus a public
   scaling; (b) **CRT batching** of the m later heads that share pivot l: P_l(x) N(x) mod prod q_i with
   N the CRT combination of the m heads (degree 2m - 1), by evaluation at 2m + 1 public points;
   (c) **bilinear rank** of GF(2^24) below 99; (d) **inversion-free (projective) pivots**, dropping the 35
   deferred pivot inverses per direction.
2. **Qubit peak.** (a), (c), (d): 890. (b): the left-looking chain must run m slots at once, one 24-bit carry
   each: 890 + 24(m - 1), 962 for m = 4.
3. **Toffolis.** (a) The transformation moves the scale into the partial numerators, q'_l = q_l /
   (b_{l-1} b_l), a quantum GF(2^12) scalar that multiplies d at every step: two 84-T products per step instead
   of one 99-T product, +69 x 1,190 = **+82k**. (b) (2m + 1) x 42 = 84m + 42 against 99m: break-even at
   m = 2.8; m = 4 saves 18 T per 4 MACs (-4.5% of the queries, -5.3k). (c) Every multiplicative construction
   gives at least 99: F8 elliptic curve 14 x 6 + 15 = 99 (Chih-Kai Wang); GF(16)^6 11 x mu_2(4) = 11 x 9 = 99;
   GF(64)^4 7 x 15 = 105; GF(256)^3 5 x 24 = 120; GF(4096)^2 3 x 42 = 126; GF(4)^12 needs mu_4(12) <= 32,
   below what curve point counts over F4 give (~34); a genus-2 curve over F8 needs 2n + g - 1 = 17 units
   (102). Rank 96 would save 3 x 1,190 = -3.6k. (d) Each later MAC then carries a scaling product (+99 x
   595 per direction) to save ~6.8k of inverses: **+110k**. Failure rate unchanged in all four.
4. **Cost.** (a) ~4.2e11 (+22%); (b) (962/890)^2 x 0.986 = 1.15, ~3.97e11 (+15%); (c) 3.418e11 (-1.0%);
   (d) ~4.6e11.
5. **Red team.** (a), (b), (d) die on arithmetic; (b) is also the generic (ab, ac) rank-additivity
   objection of DQI-B10. (c) is DQI-B08 unchanged: NR-0021/0022/0033 found nothing below 42 for GF(2^12); 99
   is reached by two independent constructions. **Dead except the rank lottery**, which stays B08 (no new ID).
6. **Prototype.** For (c) only: flip-graph or SAT search for a rank-98 decomposition of GF(2^24)
   multiplication seeded with the 99-term scheme; the 576-pair bilinear check validates a candidate in 1 s.
7. **P = 0.10** for (c) at <= -1%; 40 h. (a), (b), (d): 0.

### F3. Convert-free Forney by a Boolean quotient phase (new: DQI-B13)

1. **Datapath.** Removes the convert (Forney preimage), the projective finish and the unconvert; keeps the
   extractor, the synthesis and the scan. After synthesis the core holds sigma = E^2 + x O^2 and Omega
   (70 words, exactly like (sigma, R) today). The additive FFT (CNOT only) gives at each scan position X the
   three words w = E(X) + sqrt(X) O(X) (sigma(X) = w^2, so sigma(X) = 0 iff w = 0), b = O(X) and
   a = sqrt(nu_X Omega(X)), where nu_X is the public trace-dual mask times the public Forney normalization
   (powers of X, the gauge factor g(X)). The error value is Omega(X)/sigma'(X) = Omega(X)/O(X)^2, and
   Tr(nu Omega/O^2) = Tr(sqrt(nu Omega)/O) (Frobenius), so the scan cell applies the phase
   `[w = 0] * Tr(a/b)`. By L1 it costs the 11 ANDs of the zero test plus MC(Tr(a/b)) ANDs, all uncomputed by
   measurement. The ratio Omega/sigma' is invariant under the common projective scale of (sigma, Omega), so no
   normalizing finish is needed; at non-roots b may be 0, which the form Tr(a b^4094) handles (value 0, and
   z = 0 there anyway); at roots b != 0 because sigma is squarefree.
2. **Qubit peak.** Core 840 + FFT and cell helpers. A tower circuit for Tr(a/b) holds GF(64) intermediates
   N, T (6 + 6), N', T'' (3 + 3) and at most ~30 live AND temporaries, so the scan stage reaches ~880-890,
   below or at the extractor peak. **Q = 890** (set by the extractor, unchanged).
3. **Toffolis.** Deleted: convert pair 114,423 + finish 7,856 = **-122,279** static T. Added: 4,025 x C_q,
   with C_q = MC(Tr(a/b)) on GF(2^12)^2 (+ at most 1 for folding the zero test). Break-even C_q = 30.4.
   - Lower bound: Tr(a b^4094) has degree 12 (the inverse has degree 11), and it is the Dillon PS_ap bent
     function on 24 variables (MC >= 12 for bent functions on 2n variables): **C_q >= 11-12**.
   - Our construction (tower GF(2^12) = GF(64)[t], t^2 + t + g): b^-1 = conj(b)/N(b), so
     Tr_12(a/b) = Tr_6(T/N) with T = Tr_{12/6}(a conj(b)) = a0 b1 + a1 b0 (two GF(64) products) and
     N = b0^2 + b0 b1 + g b1^2 (one GF(64) product); recursing, C(2n) <= C(n) + 3 mu_2(n), C(3) <= 6 (every
     6-variable function has MC <= 6): C(6) <= 24 and **C(12) <= 69**.
   - Failures: the convert's degree-gap handling (part of the dominant ~3e-6/shot half-Euclid class) and its
     gap-2 coverage disappear: a possible -0.1 to -0.2% bonus, not counted below.
4. **Cost.** T = 365.9k - 122.3k + 4,025 C_q. C_q = 12: 292k, **2.75e11 (-20%)**; C_q = 20: 324k, 3.06e11
   (-11%); C_q = 30: break-even; C_q = 69: 521k, 4.92e11 (+42%). Reaching the 1.34x target needs C_q <= 7,
   below the lower bound: F3 alone cannot reach it, but it is the only family with a double-digit ceiling.
5. **Red team.**
   - *Balram Srinivasan's null* (per-position Forney needs an inverse <= 16 T; current 66) priced the
     12-bit inverse with outputs and an uncompute (~66 + 66). F3 needs only one trace bit as a phase (L1), so
     the right price is MC(Tr(a/b)), between 11 and 69, and the budget is 30 because the whole convert pair and
     finish go. It escapes the recorded null in form, not yet in number.
   - *How likely is 30?* It needs a circuit 2.3x better than the tower. The AES precedent (GF(2^8) inverse:
     tower ~36-40 ANDs, best known 32) shows optimization gaining 15-20%, not 2x. Joint sharing among the
     three GF(64) products {a0 b1, a1 b0, b0 b1} might give ~40 instead of 45, so ~60-64 overall. **P ~0.03.**
   - NR-0083 and NR-0199 (convert replacements >= 180k per direction) priced coefficient-domain chains, not
     a per-position phase: consistent, not contradicting.
   - Width and failure risk: none found (the scan is not at the peak; root positions never divide by zero).
   - **Survives conditionally**: a cheap mathematical screen decides it.
6. **Prototype (the screen).** Exact or heuristic MC of Tr(a/b) on GF(2^n)^2 for n = 2..6 (4 to 12
   variables): SAT with an AND budget for n <= 4, heuristic XOR-AND synthesis for n = 5, 6, all compared with
   the tower values C(2) <= 3, C(4) <= 12, C(6) <= 24. Decision rule: if the best found at n = 4 and 6 is
   >= 75% of the tower, record F3 as a null (C(12) ~50-69 > 30); if it is <= 45%, spend a second day on a
   GF(2^12) search seeded with the tower circuit, and verify any candidate exhaustively over all 2^24 inputs
   (bit-sliced, seconds). Only a verified C_q <= 25 justifies the build: generator change in the scan cell
   (12-bit inputs a, b, w), removal of the convert/finish/unconvert stages, official verifier gate.
7. **P = 0.03.** Effort: screen 8 h; full build ~120 h (scan cell, pipeline surgery in the v12 generator,
   gates).

### F4. Reciprocal convert (Horiguchi-Koetter, Bezout reuse from the extractor)

1. **Datapath.** Uses the continued-fraction determinant identity Omega B' - A' sigma = +-G (G the product of
   the public node polynomials; the constant is public if sigma is kept as the unnormalized last
   denominator). At roots Omega(X) = +-G(X)/B'(X), so e = +-G(X)/(B'(X) sigma'(X)): 1/e is a polynomial value.
   (a) Synthesize (sigma, B') instead of (sigma, Omega) (the other column of the same 2x2 product), form
   W = B' O^2 mod sigma in place (no Euclid), and scan with the single-variable phase `[w = 0] * Tr(nu'/W(X))`.
   (b) Skip W and scan with `[w = 0] * Tr(nu'/(B'(X) O(X)^2))`.
2. **Qubit peak.** 890 (sigma, B' or W: 70 words).
3. **Toffolis.** (a) In-place multiplication by O^2 modulo sigma is the adjoint of today's in-place division
   by O^2 modulo sigma (L2): the convert pair costs the same, nothing is saved; on top, 4,025 x C_r with
   C_r = MC(Tr(nu x^-1)) >= 10 (degree 11; with the zero test the phase has degree 23) and ~39-42 by the tower
   (the relative trace of nu conj(W) is linear, so only N (15) and the GF(64) recursion (24) remain):
   **+40k to +170k**. (b) Deletes the convert like F3 but its cell function Tr(a/(b c^2)) contains an extra
   GF(2^12) product: strictly more expensive than F3. Failure rate: as F3 for (b), unchanged for (a).
4. **Cost.** (a) >= 3.83e11 (+11%), likely ~5e11; (b) dominated by F3.
5. **Red team.** L2 settles (a): trading the Forney division for a modular multiplication saves nothing in
   place. Bezout reuse (DQI-B05, NR-0027) gives Omega^-1 or 1/e as polynomials while the trace phase needs e,
   so one inversion always remains, either modular (convert) or per position (F3). **Dead; dominated by F3.**
6. **Prototype.** None (if F3's screen opens, check the identity on a GF(2^6) toy for completeness).
7. **P = 0.01.**

### F5. Decode-re-encode reverse (a cheaper algorithm for the inverse map)

1. **Datapath.** Replaces the reverse key and the unconvert (~173k static T) by an encoder from the decoded
   state: Omega = R sigma' mod sigma, then the syndromes as the first 70 Laurent coefficients of
   Omega/sigma at infinity (a power-series division), instead of replaying the extractor backwards.
2. **Qubit peak.** In place, the map (sigma, Omega) -> s splits into a triangular half (c_0..c_34 from Omega,
   in place, cheap) and a Hankel half (sigma -> c_35..c_69 given c_0..c_34), which is the application of a
   quantum Hankel matrix: doing it in place needs its LU factorization, i.e. the Schur/Thiele key stage
   itself. Out of place the series needs a fresh 840-qubit register: Q ~ 1,730.
3. **Toffolis.** In place: triangular half 595 products (~25k) + Hankel half >= the forward key (~101k) +
   R -> Omega (= the unconvert by L2): no saving. Out of place: Omega step = unconvert (L2), series division
   sum_{k<70} min(k, 35) ~ 1,820 products (~76k); and the input register, which now holds (sigma, R), must
   still be erased, which needs a decode: the out-of-place encode does not even solve the problem.
4. **Cost.** In place >= today; out of place > 1e12.
5. **Red team.** Self-defeating ("in-place encode = decode"); NR-0084 adds that no measurement record can
   replace the recomputation (kickmix HMR outcomes are independent of the qubit). **Dead.**
6. **Prototype.** None.
7. **P = 0.00.**

### F6. Narrow-width line (12-bit carry, norm-based inverse helpers)

1. **Datapath.** Full version: keep the extractor's projective denominator in the subfield GF(2^12) (scale
   the chain pair (n, d) by conj(d) after each step) so the carry bank needs 12 qubits, and compute every
   GF(2^24) inverse through its 12-bit norm (Balram Srinivasan's closing generalized: x^-1 = conj(x) N(x)^-1,
   12 helpers instead of 24). Lite version: only the norm-based helpers, at the peak sites that are not
   out-of-place GF(2^24) multiply windows.
2. **Qubit peak.** Full: 840 + 12 carry + ~12 helpers + flags = 866-878. Lite: 888, because the out-of-place
   multiply windows hold 840 + 24 carry + 24 product (NR-0082, NR-0120, NR-0138).
3. **Toffolis.** Full: the chain step (n, d) -> ((d + P n)/q, n) makes the new d the old n, a full GF(2^24)
   value; re-normalizing it to the subfield costs one GF(2^24) product plus a norm per step:
   ~+118k for the pair. Lite: one or two qubits at +0.4k to +1.5k (NR-0081 measured 889 -> 888 at +424 T on
   another line; the norm path is even cheaper per call: 340 vs 360 T). Failures unchanged if the rare
   branches (norm in GF(64), zero source) keep their handlers.
4. **Cost.** Full at 866 Q and ~484k T: (866/890)^2 x 1.32 = 1.25, ~4.3e11 (+25%). Lite: 3.43-3.45e11
   (-0.05 to -0.35%).
5. **Red team.** NR-0260 (the 25-wire budget is spent by two independent families), NR-0082 (multiply
   windows), NR-0101/0104 (a narrower division costs more than the qubit is worth). **Full: dead. Lite:
   survives as a refinement of DQI-A09**, not as a family.
6. **Prototype.** The liveness profiler of DQI-A09 on our cutoff-17 build: list the 889/890 windows by owner;
   if the owners outside multiply windows are inverse helpers, re-emit those calls with the 12-bit-norm
   inverse and measure.
7. **P = 0.35** (lite, for -0.05 to -0.35%); 20 h.

### F7. Cross-stage retention with measurement-based uncomputation

1. **Datapath.** Hold the output p = P n of a forward query MAC in a fresh word through convert and scan, so
   that the mirrored reverse MAC becomes d ^= p (CNOT) followed by an X-measurement of p and its phase debt
   Tr(m P n), a bilinear form of rank 24 (24 CCZ): 99 + 24 = 123 T instead of 2 x 99 per MAC pair. Same idea
   for pivot inverses (390 T per slot pair).
2. **Qubit peak.** +24 per held product (each word holds one product across the whole mirror).
3. **Toffolis.** -75 T per held product against +5.47% (~20k T-equivalent) per held word: break-even needs
   >= 267 products per word, and a word can hold one. Held inverses: -390 T for +24 Q, same verdict.
4. **Cost.** >= 3.64e11 (+5.4%) per held word.
5. **Red team.** NR-0184 (reusable products straddle the forward/reverse mirror), NR-0201 and NR-0209 (holding
   a 24-bit word through the middle stages): already measured. **Dead.** Reopen only for a held quantity of
   <= 2 qubits that saves >= 1,600 T of reverse recomputation.
6. **Prototype.** The DQI-A09 profiler can list reverse recomputations whose inputs are 1-2 bits.
7. **P = 0.01.**

### F8. Quartic-node CRT (GF(2^48) heads)

1. **Datapath.** Replaces the 35 quadratic nodes by 17 quartic nodes (GF(2^48) heads) plus one quadratic
   node; the Thiele chain then runs over 18 slots instead of 35 (the extension of DQI-B12's cubic idea).
2. **Qubit peak.** 840 core + 48-bit carry + >= 48 helpers for the GF(2^48) product and inverse windows:
   ~938 (+48).
3. **Toffolis.** Queries: sum_{i<=18} (i - 1) = 153 MACs per direction at mu(GF(2^48)) ~ 216 (GF(16)^12 on
   an elliptic curve over F16: 24 x mu_2(4) = 216; GF(64)^8: 15 x 15 = 225; Karatsuba over GF(2^24): 297):
   33.0k per direction against 58.9k, **-51.8k** for the pair. Per-slot overheads (inverses, scales, affine
   phases; ~36.5k today for 35 slots): a GF(2^48) inverse via its GF(2^24) norm costs ~600 T (3x today's
   195) and the scales grow with the head size: ~36.5k x (18/35) x 2.5 ~ 47k, **+10k**. Synthesis: same total
   degree, assumed unchanged. Net T ~324k. Failure classes would move (fewer, larger slots): unpriced.
4. **Cost.** (938/890)^2 x 324/365.9 = 1.111 x 0.886 = 0.984: ~3.40e11 (-1.6%) in the best case.
5. **Red team.** The cubic precedent was measured at about +48k T (NR-0280) with +24 Q (DQI-B12): the hidden
   costs (special charts, tail controller, synthesis with cubic pivots, larger inverses) exceeded the query
   saving. Quartic doubles both the saving and the width price. **Dead / marginal.**
6. **Prototype.** None.
7. **P = 0.03**; ~150 h (whole-decoder rewrite).

### Families considered and not written up

- Textbook decoders (in-place Berlekamp-Massey, Sugiyama, synchronized EEA, Gu et al.'s history-free BM):
  ~3,675 GF(2^12) products per direction for BM in coefficient form (154k T, above today's 101k forward key)
  and 2x width for the EEA; already a recorded null (Gu et al. 4.4e5 T at 1,837 Q; Khattar et al. 5.1e6 T).
- Fourier / extended-syndrome phase (phase = sum_k Tr(c_k S_k) over the 4,095 extended syndromes): 4,025
  LFSR steps of 35 products each, ~141k products: dead (the residue-phase null in `dqi.md` section 4).
- Interleaved even/odd syndromes (S_2i and S_2i+1 as two 35-syndrome codes with common locators): collaborative
  decoding reaches ~2/3 of 35 errors, below the weight law's support: dead (L6).
- Weight truncation to <= 34: 2.77% failed shots, R ~4.2: dead (L6).

## 4. Ranked recommendation

Ranking by expected value per hour (P x midpoint gain / effort of the next step):

| Rank | Family | Next step | Effort | P | Gain if it lands | Why this rank |
|---:|---|---|---:|---:|---|---|
| 1 | F3 convert-free Forney (DQI-B13) | MC screen of Tr(a/b) | 8 h | 0.03 | -11% (C_q = 20) to -20% (C_q = 12) | Only double-digit ceiling; the screen is cheap and settles it either way |
| 2 | F6-lite (888/889 Q via 12-bit-norm inverse helpers) | DQI-A09 profiler + re-emit | 20 h | 0.35 | -0.05 to -0.35% | Small but likely; reuses tooling the A-track needs anyway |
| — | F2(c) GF(2^24) rank <= 98 (= DQI-B08) | flip-graph / SAT search | 40 h | 0.10 | <= -1% | Lottery; only if a contributor has spare cycles |

No family deserves a scaled decoder prototype now: the only one whose paper price can fall below the frontier
(F3) depends on a single number, MC(Tr(a/b)), that a math screen measures directly. A decoder prototype comes
after the screen, and only if it passes.

### First experiment for F3 (DQI-B13): the MC screen

1. **Functions.** For n in {2, 3, 4, 6} (and 8 as a sanity point against the AES literature), the Boolean
   function f_n(a, b) = Tr(a * b^(2^n - 2)) on GF(2^n)^2 (2n variables; all Tr(c a/b) with c != 0 are affinely
   equivalent through a -> c a, so one function per n). Truth tables by NumPy over the field tables
   (`km.GF2Field(n)` gives the arithmetic).
2. **Exact MC for n <= 4** (4, 6, 8 variables): SAT encoding "there is an XOR-AND circuit with k AND gates
   whose output equals f_n", each AND input a GF(2)-linear combination of the variables and earlier AND
   outputs (the standard MC encoding of Calik, Sonmez Turan and Peralta), k increasing from the degree bound
   (n - 1) until SAT. Expected values to compare: tower bounds C(2) <= 3, C(3) <= 6, C(4) <= 12.
3. **Heuristic upper bounds for n = 6 and 8** (12 and 16 variables): XOR-AND resynthesis seeded with the tower
   circuit (tower bounds C(6) <= 24, C(8) <= 39); record the best found.
4. **Decision.** Ratio rho = best found / tower bound at n = 4 and 6. rho >= 0.75: extrapolated C(12) >= 50,
   F3 recorded as a null (NULLS row below). rho <= 0.45: one more day on n = 12 (seeded search, candidate
   checked exhaustively over all 2^24 inputs, bit-sliced, seconds). Build only on a verified C_q <= 25.
5. **Record** run time, peak memory and every MC value in `experiments/LEDGER.md` with the verdict.

### First experiment for F6-lite (as a step of DQI-A09)

1. Extend `tools/kmb.py` with an exact live-qubit profile (allocation intervals per instruction) and run it on
   our cutoff-17 build `cc0d7d25` (889 Q).
2. List every window at width 889 and its owners (generator call site by instruction pattern: GF(2^24)
   inverse helpers, out-of-place multiply windows, flags, tail controller).
3. Go if the 889th wire outside multiply windows is always an inverse helper: re-emit those inverse calls with
   the 12-bit-norm closing (DQI-A02 gadget) and measure static T, then the official verifier on 8 salts.
   Accept at 888 Q if the Toffoli increase is below ~860 (one qubit's worth at 889 Q).

### What would change the verdict on dqi families

- MC(Tr(a/b)) over GF(2^12) <= ~25 (F3 build, -11% to -20%).
- A bilinear rank of GF(2^24) multiplication <= 96 (B08, ~-1%).
- An in-place convert below ~1,000 products per direction (NR-0199's reopen condition), or any in-place map
  (sigma, Omega) -> (sigma, R) cheaper than today's half-Euclid frame.
- A convergent representation that stores one pair plus O(1) words and recovers the previous pair in O(k)
  products (F1, then queries and synthesis both shrink).
- Otherwise the in-place CRT/Thiele line looks close to the floor of its class (L1-L6): the leaders' -30% is
  more likely on chemistry, or on a dqi idea outside this class that none of the public material shows.

## 5. Proposed register updates

Not applied (this step writes only this file). The caller applies them with the changelog entry.

### `Analysis/IDEAS.md` (dqi, category B)

| ID | Idea | Est. gain | P | Effort (h) | Status | Notes |
|---|---|---|---:|---:|---|---|
| DQI-B13 | Convert-free Forney by a Boolean quotient phase: scan cell `[w = 0] * Tr(a/b)` on (E, O, Omega), convert + finish + unconvert deleted (-122k static T) | -20% (C_q = 12) to -11% (C_q = 20); break-even C_q = 30.4; +42% at our tower C_q = 69 | 0.03 | 8 (screen) + 120 (build) | OPEN, screen first | F3 of this file; lower bound C_q >= 11-12 (degree, bent), construction 69 |
| DQI-B14 | Explicit or blocked convergents in the extractor ("free queries") | -60k T only if width were free | 0.01 | — | NULL (paper) | F1: 70 + 2k + 4 words at slot k (+840 Q at k = 35); blocks of 12 need 624 Q |
| DQI-B15 | Re-based query MAC: monic pivots, shared-pivot CRT batching, inversion-free pivots | — | 0 | — | NULL (paper) | F2: +82k; break-even at 2.8 heads with +24 Q per head; +110k. Rank part stays DQI-B08 |
| DQI-B16 | Reciprocal convert (Horiguchi-Koetter with the synthesized penultimate denominator) + inverse-trace scan | — | 0.01 | — | NULL (paper) | F4: in-place multiply = in-place divide (L2); dominated by DQI-B13 |
| DQI-B17 | Decode-re-encode reverse | — | 0 | — | NULL (paper) | F5: in-place encode needs the Hankel LU (= key stage); out of place +840 Q |
| DQI-B18 | Subfield (12-qubit) carry bank | — | 0 | — | NULL (paper) | F6 full: +118k T for -12 to -24 Q. Lite part becomes a DQI-A09 sub-step |
| DQI-B19 | Hold forward products/inverses across the mirror, MBU on the reverse copy (99 -> 24) | — | 0.01 | — | NULL (paper) | F7: break-even >= 267 products per held 24-qubit word |
| DQI-B20 | Quartic-node CRT (GF(2^48) heads) | -1.6% best case | 0.03 | 150 | NULL (paper) | F8: queries -51.8k, width +48 Q, overheads +10k; cubic precedent NR-0280 |

Status changes for existing rows:

- **DQI-B05** (Bezout reuse): PARTLY NULL (paper) — Bezout reuse gives Omega^-1 or 1/e, and turning that into
  e in place costs what the convert costs (L2); only the per-position route of DQI-B13 remains.
- **DQI-B09** (explicit incremental synthesis): NULL (paper) unless a one-pair representation is found
  (exact count in F1).
- **DQI-B10** (two slots, shared pivot): NULL (paper) for m <= 4 slots (F2b arithmetic).
- **DQI-B11** (12-qubit carry): NULL (paper) for the carry (+118k); see DQI-B18.
- **DQI-B12** (cubic nodes): unchanged (NR-0280 already); quartic extension recorded as DQI-B20.
- **DQI-A09**: add the F6-lite step (re-emit inverse helpers at non-multiply peak windows with the 12-bit-norm
  closing).

### `experiments/NULLS.md` (dqi, "From our research")

| Ref | What was tried | Result | Source |
|---|---|---|---|
| DQI-B14 | Explicit/blocked convergents in the CRT/Thiele extractor (free public evaluations instead of chains) | width 70 + 2k + 4 words at slot k: +840 Q at slot 35 (cost x3.3); blocks: break-even b ~ 6, b = 12 needs 624 Q of scratch | `dqi_families_2026-11-19.md` F1 |
| DQI-B15 | Monic (subfield) pivots; CRT batching of shared-pivot products; inversion-free pivots | 2 x 84 T instead of 99 (+82k); (2m+1) x 42 vs 99m, break-even m = 2.8 heads at +24 Q each; +110k | F2 |
| DQI-B16 | Reciprocal (Horiguchi-Koetter) convert with inverse-trace scan | in-place multiply by O^2 mod sigma costs the same as the division it replaces; +4,025 x (>= 10) T per scan | F4 |
| DQI-B17 | Decode-re-encode instead of replaying the reverse key | in-place encode = Hankel LU (the key stage); out of place +840 Q and still needs a decode to erase | F5 |
| DQI-B18 | 12-qubit (subfield) carry bank | renormalizing the chain pair each step: +~118k T for -12 to -24 Q | F6 |
| DQI-B19 | Holding forward products across the mirror, MBU on the reverse copy | -75 T per product vs +20k T-equivalent per held 24-qubit word: needs >= 267 products per word | F7 |
| DQI-B20 | Quartic-node CRT (GF(2^48) heads) | queries -51.8k but +48 Q (+11.1%) and GF(2^48) inverses (+10k): -1.6% at best, cubic precedent +48k | F8 |
| (law) | Phase-oracle pricing | a per-position phase costs MC(f) Toffolis (ANDs computed, uncomputed by measurement); an in-place map and its inverse cost the same | section 1 |

If the DQI-B13 screen fails (rho >= 0.75), add: "DQI-B13 convert-free Forney by Boolean quotient phase |
MC(Tr(a/b)) at n = 4, 6 within 75% of the tower; extrapolated C(12) >= 50 > break-even 30.4 | screen".

### `Analysis/PLAN.md`

- S2.21: DONE on paper. 8 families, 7 structural nulls, 1 conditional survivor (DQI-B13, screen 8 h).
- M2.C go/no-go: not met unconditionally. Proposed: run the DQI-B13 screen as a short step (S2.21b) before
  M2.D; move the January family effort to chemistry first; dqi stays on the A-track (A01 source recovery,
  A07/A17 failure budget, A09 with the F6-lite step).

