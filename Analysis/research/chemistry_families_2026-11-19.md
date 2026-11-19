# Chemistry (dfthcsos) families on paper, with an integrated red team (S2.23, GEN-B11)

Date: 2026-11-19. Scope: plan step S2.23, the chemistry counterpart of S2.18 (`invsqrt_families_2026-11-19.md`)
and S2.21 (`dqi_families_2026-11-19.md`). Question: can any qualitatively new family lower the dfthcsos cost by
−20% to −35% (the scale the leaders' ≈ +12% score suggests), and which 1–2 families deserve the January P100 time?

Base (public frontier, pinned in `third_party/public-frontier-2026-11-19/PINS.md`): **294,642,169.5 =
λ_eff 46.9216 × T 14,706 × Q 427**, Lera's (R,B,C) = (8,143,44) coded table, maoxin123's 15 loaders at 4,022
Toffolis per pass, weights `6ef06081`. All prices below are against this base unless stated. Sources read: the
scorer (`verify.py`, `scoreinfra/costs_prepare.py`, `scoreinfra/circuits.py`, `SCORING_README.md`),
`Analysis/research/chemistry.md`, `experiments/NULLS.md` (chemistry rows, full text from
`third_party/manwithacat-vqcc-negative-results/ledger.csv`), `public_frontier_asterios_2026-11-19.md`.

## Summary

1. **No −20% chemistry family exists on paper under this scorer.** Q is pinned at 427 for every viable shape
   (the inner-QROAM block size is chosen by the cost model, not by us), λ_eff is 0.28% above a hard floor, and the
   Givens and state-preparation Toffolis (44% of T) are fixed functions of (R, B, C, brot). Only the loaders
   (55% of T) and the shape are movable, and both are boxed: see §0.
2. **New measurement (this study): the loader product count does not depend on the number of table rows.**
   The sum of per-plane low×high core ranks of the shipped (8,143,44) table is 8,443; dropping a whole group
   gives 8,444, dropping two gives 8,445, keeping only 512 rows gives 8,395 (script and numbers in §0.3). Only
   the per-circuit address cache (≈ one AND per 16-row block, per circuit) scales with rows. **A row is worth
   1.9 Toffolis (0.013%), not 6.6 (0.045%)** as `chemistry.md` §2.3 assumed.
3. Consequence: the shape moves everyone has been chasing are worth 3–4× less than believed: **R 8 → 7 at C 50 ≈
   −2.0%** (not −7.5%), **B 143 → 129 ≈ −1.6%** (not −5.2%), R 8 → 6 ≈ −3.5% (not −13%). brot 8 → 7 is worth
   ≈ −5 to −6% (not −10%). Only the C axis keeps its value (−0.53% per unit of C). With accuracy odds of 0.03–0.05
   (NR-0197, NR-0225, NR-0236, NR-0276), the row-axis moves leave the January list.
4. The single biggest prize in the cost model, forcing the inner QROAM block from 2^4 to 2^3 (Q 427 → 276,
   T +20.5%, **−22% cost**), is unreachable: it needs brot ≤ 3 or R·C·B < 2,432 (F2, structural null).
5. Seven families were priced and red-teamed (table in §1). Three nulls on paper (F2 structural, F4 killed by
   the one-body rows plus the NR-0249 tail numbers, F5 a net loss), one re-priced null (F1), three survivors, all
   small: non-bilinear loaders F3 (ceiling −15%, realistic −1 to −3%, P 0.1, CPU only), rank-budgeted codes on the
   C axis F6 (−1 to −3%, P 0.2, P100), helper/cache parking F7 (≤ −1%, P 0.2, CPU).
6. **Realistic chemistry ceiling for January: −2% to −5% task cost (+0.7% to +1.7% score).** Chemistry cannot
   supply the leaders' −30% total-cost jump; if they are far ahead on chemistry, they have a much better fit
   (continuous e2 well below 13 mHa at the frontier shape), not a cost-model trick visible from here.
7. Recommendation for January (§4): (1) **F3** first, on CPU, now: an exact small-case screen of whether
   non-bilinear AND circuits beat the low×high pairing on random multi-output tables (the NR-0231/NR-0268 reopen
   condition); (2) **F6** on the P100: one budget optimiser (CHEM-B05) that spends the 20 mHa on the C axis and on
   core rank at the exchange rates of §0.2, starting from a warm C = 42 fit. R 7, B 129 and brot 7 runs should not
   get P100 time.

## 0. The box: what the scorer leaves movable

### 0.1 Cost anatomy in the units the loaders are really priced in

`score = λ_eff × T × Q` (lower is better). Frontier values:

| Term | Formula (scorer) | Value | Share | Who moves it |
|---|---|---:|---:|---|
| Loader address caches | ≈ 2 passes × K circuits × H high blocks (H = ⌈rows/16⌉ = 77; 75–80 ANDs per circuit measured, NR-0231, NR-0128) | ≈ 2,310 | 15.7% of T | rows (shape), K = ⌈75/J⌉ |
| Loader products | ≈ 2 × 0.42 × Σ core rank (NR-0128 calibration, reproduces 4,168/pass; NR-0231: 195 products for core-rank sum 470) | ≈ 5,734 (Σ core rank ≈ 6,830) | 39.0% | word table (codes, gauges), loader architecture |
| Inner state prep | `RC·⌈(B+1)/2^k2⌉ + 2^k2·b2 + QR(RC,B+1) + (B+1) + N + b2 + 8⌈log2(B+1)⌉`, k2 = 4 | 4,095 | 27.8% | C, R (and B in steps of 16) |
| Givens rotations | `2(N−1)(2·brot − 1)` | 2,250 | 15.3% | brot only |
| Outer state prep | alias sampling over `X_o = N + RC` | 317 | 2.2% | RC (weakly) |
| **T** | | **14,706** | | |
| Inner temporary | `⌈log2(RC·⌈(B+1)/16⌉)⌉ + 15·b2 + 6` | 303 | 71% of Q | nothing viable (§F2) |
| Persistent (outer + inner) | `b1 + χ + uniform` / `b2 + χ + ⌈log2 R⌉ + 2 + uniform` | 84 | 20% | R (one qubit at R ≤ 4) |
| Word registers | `J·brot` | 40 | 9% | J, brot |
| **Q** | `84 + 40 + max(90, 303, ≤ 77)` | **427** | | |
| λ_eff | `sqrt(E_gap(2Λ − E_gap))` ≥ `sqrt(559.13 · 3.914898)` = 46.79 | 46.9216 | 0.28% above floor | weights only |

The products line uses the frontier's own split (8,044 − 2,310); its exact Σ core rank depends on the rank
definition (helper planes, codes), but every conclusion below uses only ratios.

### 0.2 Exchange rates (task cost, at the frontier point)

| Unit | Toffolis (both passes) | Task cost |
|---|---:|---:|
| 1 Toffoli | 1 | −0.0068% |
| 1 core-rank unit (one plane, one rank) | 0.84 | −0.0057% |
| 1 AND removed from every circuit's cache | 30 | −0.20% |
| 1 table row (1/16 of a cache AND per circuit) | 1.9 | −0.013% |
| 1 leaf (R·C − 1) | ≈ 9.7 | −0.066% (C − 1 at R 8 = −0.53%) |
| 1 qubit of Q | — | −0.234% |
| 1% of λ_eff | — | −1% |

**What −20% costs:** −2,941 Toffolis at fixed Q and λ_eff (37% of all loader Toffolis), or −85 qubits, or a
combination. §F2 shows Q cannot move by more than a few qubits; §F3 shows the loaders' information-theoretic floor
is ≈ −15%. So −20% is out of reach even before accuracy is considered.

### 0.3 Measured: core rank is invariant under row removal (new)

Script `core_rank.py` (scratchpad of this session; 40 lines, NumPy, 1 s, < 50 MB): for each of the 15 batches and
40 delta planes (`delta[k] = batch[k] XOR batch[k−1]`, J = 5, brot 8), build the H × 16 matrix (16-row high block ×
4 low address bits), Möbius-transform over the low bits, drop the constant coefficient and take the GF(2) rank
(≤ 15). Padded cells set to 0. Table: Syed Asad Ali's v12 (8,143,44) `angle_words` (1,220 × 75).

| Table | Rows | H | Σ core rank | Planes at rank 15 |
|---|---:|---:|---:|---:|
| Shipped | 1,220 | 77 | 8,443 | 545 / 600 |
| Random words, same shape (control) | 1,220 | 77 | 9,000 | 600 |
| Last group dropped (7 groups, re-packed) | 1,077 | 68 | 8,444 | 545 |
| Two groups dropped (6 groups) | 934 | 59 | 8,445 | 545 |
| 14 random rows per group dropped (B 129) | 1,108 | 70 | 9,000 | 600 |
| First 512 rows only | 512 | 32 | 8,395 | 545 |

Mean rank per bit (bit 0 … bit 7): 7.6, 15, 15, 15, 15, 15, 15, 15. The rank of an H × 15 matrix saturates at
15 as soon as H ≥ 15 blocks are generic, so **the product side does not shrink with the row count** (until
fewer than ≈ 250 rows); random row removal even destroys the block-aligned bit-0 code (8,443 → 9,000). Only the
cache side (≈ H ANDs per circuit) scales with rows. This contradicts the linear "3.34 Toffolis per row per pass"
assumption of `chemistry.md` §2.3 and the 0.934× stake quoted for (7,143,50) in NR-0236.

### 0.4 Re-priced shape and precision moves (cost model + core-rank loader model)

| Move | Cost-model ΔT | Loader ΔT | ΔQ | Task cost (new) | Old estimate | Accuracy evidence |
|---|---:|---:|---:|---:|---:|---|
| R 8 → 7, C 50 (RC 350) | −20 | −270 (H 77 → 68) | 0 | **−2.0%** | −7.5% | NR-0236: 19.1 mHa continuous (needs ≤ 13); NR-0276 |
| R 8 → 6, C 59 | +20 | −540 (H 59) | 0 | **−3.5%** | −13% | NR-0197: +6.5 mHa over the gate; NR-0078 |
| B 143 → 129 (R 8, C 44) | −28 (QR dag, B+1) | −210 (H 70) | 0 | **−1.6%** | −5.2% | NR-0225: 19.96 mHa continuous |
| C 44 → 43 / 42 / 40 | −78 / −156 / −312 | 0 | 0 | −0.53% / −1.06% / −2.1% | same | Asterios B18 (fixed words: 21.1 mHa); EY floors 8.4 / 10 / 12.5 mHa |
| brot 8 → 7 | Givens −300 | −250 to −480 (drops the cheap coded bit-0 planes) | −5 | **−4.9% to −6.4%** | −10% | NR-0124: δ7 ≥ 19.3 mHa, e2 ≈ 23–26 mHa |
| J 5 → 4 / 6 | 0 | +4 / −2 circuits × 77 × 2 = +616 / −308 | −8 / +8 | +2.2% (worse) / −0.3% at best | — | J 6 needs the ≈ 77-wire cache in 69 workspace qubits; NR-0128 measured +5.8% with borrowed wires |

## 1. Comparison table

| # | Family | Kind | Terms moved | Ceiling | Realistic cost | vs 294.64M | P | Effort | Red-team verdict |
|---|---|---|---|---:|---:|---:|---:|---|---|
| F1 | Row-axis shape moves (R 7, R 6, B 129) | factorization regime | caches (T), small cost-model T | −3.5% | −1.6 to −2.0% if accurate | 288.8M–290.0M | 0.03–0.05 | 60 h P100 | **NULL (re-priced)**: value cut 3–4× by §0.3; accuracy nulls stand |
| F2 | Force a smaller inner QROAM block (k2 = 3), J, shared max | cost-model terms | Q −151, T +3,016 | −22% | 0 | — | 0 | — | **STRUCTURAL NULL**: needs brot ≤ 3 or RCB < 2,432 |
| F3 | Non-bilinear (multi-level) loaders | loader architecture | products + caches | −15% (counting floor ≈ 195 ANDs per circuit vs 268) | −1 to −3% | 286M–292M | 0.10 | 80 h CPU | **SURVIVES** (only lever with no accuracy cost; long shot) |
| F4 | Block-localized dictionary (16-row blocks share an orbital region) | factorization regime + word structure | products via core rank (≤ active blocks per column), rounding | −3% (regions of 7 orbitals) to −12% (4 orbitals) | 0 | — | 0.01 | (1 h CPU screen, not recommended) | **NULL ON PAPER**: the 76 one-body rows keep ≥ 5 generic blocks in every plane; NR-0249 truncation at k ≤ 8 is 1.2–2.3 Ha |
| F5 | Rows for precision (dithered twin vectors) | precision encoding | rounding noise → C budget | +4% to −1% | net loss | — | 0.02 | — | **NULL ON PAPER**: twins cost more rows/B-steps than the C cut they buy |
| F6 | Rank-budgeted codes on the C axis (B01 + B05 + B08/B12 re-priced) | word codes + shape | core rank, leaves | −4% | −1 to −3% | 286M–292M | 0.20 | 50 h P100 + 10 h CPU | **SURVIVES** (best P100 use) |
| F7 | Cache amortization by parking more address functions in gauge-free word bits | cost-model term (per-circuit cache) | caches | −1% | −0.2 to −0.6% | 293M–294M | 0.20 | 15 h CPU | **SURVIVES, SMALL** (rule risk 757522) |

The families stack only partly: F3 and F7 act on loaders for any table; F6 and F4 change the table; F1 and F6
spend the same 20 mHa budget.

## 2. Families

### F1. Row-axis shape moves: R 8 → 7 or 6, B 143 → 129 (re-priced)

1. **Idea and terms moved.** Fewer two-body rows (R·B) at about the same pair-space rank R·min(B,C): drop a group
   and widen C (R 7, C 50; R 6, C 59), or drop 14 vectors per group (B 129, the last B that keeps b2 = 19 and
   k2 = 4). The cost model barely moves (RC within ±2; `⌈(B+1)/16⌉` stays 9; Q unchanged, `⌈log2 R⌉` stays 3), so
   the whole gain is in the loaders. §0.3 shows the product side is invariant under row removal (Σ core rank
   8,443 → 8,444 → 8,445); only the per-circuit cache shrinks (one AND per 16 rows per circuit).
2. **Accuracy.** All three are measured nulls at the continuous stage: (7,143,50) warm merge 19.1 mHa continuous
   against ≤ 13 needed (NR-0236); (7,143,60) 14.2 continuous but 23.6 rounded (NR-0276); (6,143,59) +6.5 mHa over
   the gate (NR-0197); (8,129,44) 19.96 continuous (NR-0225); even B 142 fails (NR-0247).
3. **Cost.** R 7: −2.0% (288.8M); R 6: −3.5% (284.3M); B 129: −1.6% (289.9M) — table §0.4. The linear row model
   had promised −7.5%, −13% and −5.2%.
4. **Red team.** Strongest objection: the move's value collapsed while its accuracy risk did not. At P ≈ 0.04 the
   expected gain of R 7 is −0.08%, below the cost of one P100 session. Hidden term checked: with fewer rows the
   `selection` register can drop to 10 bits (R 6: 934 rows), which frees one workspace qubit under the 128 cap but
   nothing in Q_total. **Verdict: NULL (re-priced).** The accuracy reopen conditions of NR-0236/NR-0276 stay
   valid, but they are no longer worth chasing for the row count alone.
5. **First experiment (done here).** The core-rank invariance check (§0.3, 1 s on CPU). Optional confirmation
   (CPU, 10 min): compile Syed Asad Ali's planner (`dfthcsos_source/compiler/`) on the 7-group truncated table and
   check that each loader's cache count falls by about 9 while its product count stays within ±5.
6. **P and effort.** P(accuracy) 0.03–0.05; effort 60 h of P100 per shape. Not recommended.

### F2. Move the cost-model terms: a smaller inner QROAM block, J, the shared max

1. **Idea and terms moved.** Q is 71% the inner-QROAM temporary register, `⌈log2(RC·⌈(B+1)/2^k2⌉)⌉ + (2^k2 − 1)·b2`.
   At k2 = 3 instead of 4: inner temp 303 → 152, **Q 427 → 276 (−35.4%)**, while `for_qroam` grows from
   352·9 + 16·19 = 3,472 to 352·18 + 8·19 = 6,488 (**T +3,016, +20.5%**). Net **−22.1%**: the largest single prize
   in the scorer, because the model minimises Toffolis under a qubit ceiling while the score prices T × Q.
   Related terms: J·brot (40 qubits), the shared `max(outer_temp, inner_temp, angles_temp)` (the ANGLE circuits could
   use 303 workspace qubits "for free"), the uncounted `selection` register.
2. **Feasibility (scorer rule, not accuracy).** `find_qroam_size_to_minimize_global_temp_ancilla` starts from
   `k2_opt = QR(R·C·B, b2)[0]` (= 6 here) and keeps the **largest** k whose temp is `< brot·(N − 1)`. k2 = 3 is chosen
   only if (a) temp(k = 4) ≥ 75·brot, i.e. `15·b2 + 12 ≥ 600` at brot 8, so b2 ≥ 40 (B > 2^28), or brot ≤ 3; or (b)
   `k2_opt ≤ 3`, i.e. `R·C·B < 128·b2 ≈ 2,432`, which with R·C ≥ 290 means B ≤ 8 (the small-B regime; leaf ranks
   ≥ 76 are needed, `chemistry.md` §2.5). Neither is physical. The ANGLE circuits cannot use the 303-qubit max
   because of the 128-qubit cap per circuit (angles_temp ≤ 128 − 11 − 40 = 77). J: §0.4 (+2.2% at J 4; J 6 does not
   fit the cache, NR-0128). `selection` is fixed at ⌈log2(76 + R·B)⌉ bits and must be unchanged after every
   circuit, so it cannot hold workspace or caches.
3. **Cost.** 0 (no legal point reaches the prize).
4. **Red team.** The objection is the rule itself; no hidden term opens it: `inner_state_prep` lowers k2 for
   B < 2^k2 on the Toffoli side only (B < 16, physics-dead), and `inner_state_prep_qubit_cost` uses the same k2 as
   the Toffoli side for B ≥ 2^k2. **Verdict: STRUCTURAL NULL.** The inner temporary register stays at 303 ± 2
   qubits for every shape with B ∈ [129, 256], brot ∈ [4, 8] and R·C·B ≥ 2,432, so Q moves only through J·brot and
   a few persistent bits.
5. **First experiment.** None needed (closed form). A 1-minute grid over `SpinFreeData_V1` for R ∈ [1, 16],
   B ∈ [65, 512], C ∈ [20, 120], brot ∈ [6, 8] can be attached to the register entry as proof.
6. **P and effort.** P = 0.

### F3. Non-bilinear (multi-level) loaders: escape the low×high pairing model

1. **Idea and terms moved.** Every public loader is a low×high bilinear construction: an address cache (about one
   AND per 16-row block) plus products of (low form ⊕ cache form) pairs, each product delivering two rank-one
   low×high terms (Syed Asad Ali's matroid-parity pairing). Its own counting floor is 190 paired products + 80
   caches = 270 per circuit (NR-0268), where the frontier already is (268). A general AND/XOR circuit (XAG) whose
   AND inputs are dense linear forms over all live wires is not bound by that floor. **Counting bound (ours):** a
   delta table of 35 non-helper planes × 1,220 rows carries ≈ 42,700 bits (≈ 38,000 after the codes); an AND in a
   circuit with ≤ 88 live wires (11 address + ≤ 77 workspace) is described by ≈ 2 × 88 + 35 (fan-out to planes) + 7
   (wire freed) ≈ 218 bits, minus log2(K!) for gate orderings, so K ≥ ≈ 185–200 ANDs per circuit. NR-0050
   (188–203 products, heuristic) and Asterios Panoras (fully mixed floor ≈ 187) agree. Moves T only (caches and
   products together).
2. **Accuracy.** None: the loader loads the same words (any table, any weights). Correctness is proven by the
   existing all-branch prover and the verifier.
3. **Cost.** Ceiling: 268 → ≈ 195 ANDs per circuit = −73 × 30 = −2,190 T = **−14.9%** (250.8M). Realistic, if a
   construction recovers 10–20% of the gap: −7 to −15 ANDs per circuit = **−1.4% to −3.1%** (285.5M–290.5M).
4. **Red team.** (a) Counting floors are not constructive: for random functions the best known explicit
   multi-output constructions are split-and-cache schemes like the frontier's, and NR-0268 found flip-graph walks
   stall at 277–280 products on both the real and a random tensor. (b) The table is random-like (Σ core rank within
   6% of random, §0.3), so no data structure helps. (c) The 128-qubit cap limits live wires to about 77, which
   shrinks the bits per AND and pushes the floor up. (d) Small-n evidence cuts the other way: every 6-variable
   Boolean function has multiplicative complexity ≤ 6 (Çalık, Sönmez Turan, Peralta), while a 3/3 bilinear split
   of a random 6-variable function costs about 11–15 ANDs, so the bilinear method is far from optimal at small n;
   whether that gap survives at n = 11, m = 35 is the open question. **Verdict: SURVIVES** as the only loader lever
   without accuracy cost, but as a long shot.
5. **First experiment (CPU, 2–4 h, no P100).** Exact minimum-AND synthesis (a SAT encoding of an XAG with K ANDs,
   or a published exact-MC tool if installable offline) on random multi-output functions with don't-cares at
   growing size: (n, m) = (5, 2), (6, 2), (6, 3), (7, 2), 20 instances each; compare with the bilinear low×high
   count at the same size from our core-rank model. **Go** if the exact optimum stays ≥ 15% below the bilinear
   count up to n = 7 with no sign of closing; then E2: a 3-level construction (low 3 / mid 4 / high 4 with cached
   mid×high partial products) priced by the same counting on one real loader. **Kill** if the gap closes below 5%
   by n = 7. This is the NR-0231 reopen condition ("exact/SAT search on 2–3 outputs that finds cross-plane
   sharing").
6. **P and effort.** P(≥ −1%) ≈ 0.10, P(≥ −3%) ≈ 0.03; 80 h CPU in total (E1 4 h, E2 20 h, generator 50 h).

### F4. Block-localized dictionary: make planes low-rank by construction

1. **Idea and terms moved.** A plane's core rank is at most the number of 16-row blocks in which it is not
   constant. Train a dictionary in which every 16-row block (rows are a free gauge inside a group) is supported on
   a shared orbital region S of s orbitals forming a subtree of the Givens tree: outside S the block's words are
   exactly 0 (no rotation) or 64 (π/2 on the path from the root), so they are constant in the block and add no
   rank. If regions spread evenly, each column is active in about s blocks, so the rank of its planes is ≤ s.
   Exact 0 and 64 words also carry no rounding error. Moves the product side of T (and possibly the rounding budget).
2. **Accuracy.** NR-0249 measured hard truncation of each row to its own top-k orbitals: e2 = 2,257 / 1,616 /
   1,246 / 882 / 656 / 399 mHa at k = 4 / 6 / 8 / 12 / 16 / 24, against 20 mHa; block-shared regions are a stricter
   constraint than per-row top-k. Even with the ≈ 20× recovery a long refit gave after a row drop (NR-0225: 378 →
   19.96 mHa), k ≤ 8 lands at ≥ 60 mHa.
3. **Cost.** The 76 one-body rows (high blocks 0–4) must be eigenvectors of h1', which are delocalized in the
   split-localized basis, so every plane keeps about 5 generic blocks: rank ≤ 5 + s', with s' the active two-body
   blocks per column. The frontier's average non-helper plane rank is about 13, so the family pays only for
   s' ≤ 7: s' = 7 → −8% products → **−3%**; s' = 4 → −31% products → **−12%** (259M). A separate sub-loader for the
   one-body rows would cost about +80 ANDs per circuit (+16%), so it cannot remove that floor.
4. **Red team.** Two independent kills: the one-body floor (regions must be ≤ 7 orbitals to pay at all) and the
   tail numbers (regions of ≤ 8 orbitals are 1.2–2.3 Ha off before refit). Also, a plane's rank drops only when
   fewer than 15 blocks are active on it (all or nothing), and important orbitals (Fe 3d) would be covered by many
   blocks. **Verdict: NULL ON PAPER.** Reopen only with a loader that handles a block-sparse table below the rank
   model, or a scoring change that loosens e2 (NR-0249's own reopen conditions).
5. **First experiment (optional, CPU, 1 h).** On the shipped table: per-row supports at 99% / 99.9% amplitude,
   greedy within-group permutation into blocks by support overlap, active blocks per column, and the post-hoc
   block-union truncation e2 through `evaluate_words`. Kill if the median active blocks per plane exceed 9 at any
   support level whose truncated e2 is below 1 Ha.
6. **P and effort.** P ≈ 0.01; 1 h screen, then 60 h P100 if it ever passed.

### F5. Rows for precision: dithered twin vectors

1. **Idea and terms moved.** The binding accuracy term is 8-bit rounding noise (≈ 12.5–14.6 mHa in quadrature at
   the pick, NR-0236/NR-0276; NR-0124's reopen condition asks for a 2× lower rounding term). Replace a sensitive
   dictionary vector u (weight w) by two grid vectors v1, v2 in the same group, each angle rounded down or up
   independently in each twin, with weights w/2 each: per angle the mean first-order error takes the values
   −f, 1 − f or ½ − f (f the fractional part in LSB), i.e. a half-LSB grid, so the twin pair behaves like a 9-bit
   vector to first order (the second-order εεᵀ term is small). Lower rounding noise frees continuous budget,
   spent on C (−0.53% per unit of C). Since §0.3 made rows cheap (1.9 T), the idea deserved a re-check.
2. **Accuracy.** Twins must sit in the same group (leaves only combine within a group). Inside B ≤ 143 a twin
   evicts a vector: NR-0247 shows even one vector per group costs more continuous fit than the refit recovers.
   Outside, B 143 → 159 gives 16 twins per group (128 in total).
3. **Cost.** B 159: inner `for_qroam` +352 (one more 16-block per leaf), QR-dagger and `(B+1)` terms +32, caches
   +8 per circuit = +240; Q unchanged (b2 19, k2 4, log term 12): **+624 T = +4.2%**. Benefit: if the 128 most
   sensitive vectors carry 50% of the rounding power, twinning halves their amplitude, rounding falls 12.5 → 9.9
   mHa and the continuous headroom rises by 1.8 mHa (√(400 − 98) − √(400 − 156)); at ≈ 0.5 mHa of continuous fit
   per unit of C (NR-0153: C 42 → 40 at R 9 cost 1.1 mHa), that buys C −3.6 = −1.9%. Net **+2.3% (a loss)**. Even
   if those 128 vectors carried all of the rounding power (rounding 6.25 mHa, headroom +3.4 mHa, C −6.8 = −3.6%), the
   net is +0.6%. Evicting 8 vectors per group instead (no B step) gives √(17.4² + 11.0²) = 20.6 > 19.8 mHa: worse.
4. **Red team.** The arithmetic is the objection; also NR-0277 (rounding damage tracks vectors per generator, and
   adding vectors to gain slack costs more loader Toffolis than it saves), and joint (Babai) rounding already
   exploits error cancellation between near-duplicate vectors, of which twins are a special case.
   **Verdict: NULL ON PAPER** for any concentration of the rounding damage.
5. **First experiment.** None needed. Optional (CPU, 30 min): the per-vector first-order rounding-sensitivity
   concentration curve of the shipped table, recorded with the null.
6. **P and effort.** P ≈ 0.02 (only if the mHa-per-C rate were several times better than measured).

### F6. Rank-budgeted codes on the C axis (CHEM-B01 + B05 + B08/B12, re-priced)

1. **Idea and terms moved.** One optimiser (the CHEM-B05 budget optimiser) that spends the 20 mHa budget at the
   §0.2 rates instead of filling the gates with λ_eff: (i) the C axis, the only shape axis whose value survives
   §0.3 (−0.53% per unit of C, −8 leaves: inner QROAM, its dagger, the outer prep); (ii) core rank, priced at
   0.0057% per unit, bought only on the planes with the lowest accuracy price per rank unit (sensitivity-weighted
   codes such as bit-0 RM(1,4) per block or Lera's XOR-of-five codes at bits 1–2, extended plane by plane), trained
   in rather than imposed after rounding. Moves T (inner prep, loader products); λ_eff stays at its floor.
2. **Accuracy.** C 43 at fixed words is 21.1 mHa (Asterios B18), so a retrain is required. EY floors: C 42 ≈ 10 mHa,
   C 40 ≈ 12.5 mHa at RC 336 / 320. The 9-group lineage fits C 42 at 12.86 mHa continuous (Lera) and C 40 at
   13.93 (NR-0153): R 9 has more per-group capacity, R 8 C 42 is untested warm. Cold starts at C 42 failed
   (25.8 mHa, NR-0052/NR-0146); warm starts from a coded C 44 table were never tried (CHEM-B08 note).
3. **Cost.** C 42: −1.06% (291.5M); C 40: −2.1% (288.4M); 200–400 core-rank units: −1.1% to −2.3%. Realistic
   combination under one budget: **−1% to −3%** (285.8M–291.7M).
4. **Red team.** (a) The budget is already spent: error_1 19.61 and error_2 19.9988 mHa on our v12, so every C unit
   or code is paid only by a better fit. (b) Codes imposed after rounding are infeasible (NR-0278: 3,288 → 3,277
   products per pass, infeasible at fixed weights; NR-0185) and plane-constrained QAT from points far off the grid
   fails (NR-0140); its reopen condition is a near-grid start, which a warm C 44 → 42 table is. (c) The
   address-permutation gauge (CHEM-B03) cannot help much: within-group row permutations carry ≈ 6,970 bits
   (8 × log2 143! + log2 8! + log2 76!), and lowering one plane's rank by one needs its 77 block coefficient
   vectors orthogonal to one nonzero low form (≈ 77 − 15 = 62 bits), so at most ≈ 110 rank units (−0.6%) before
   any conflict between the 600 planes or with the block-aligned codes (NR-0154 found the weight-based version
   worse). (d) λ_eff is 0.28% above its floor, so the budget cannot be bought back through λ.
   **Verdict: SURVIVES** — the best use of the P100.
5. **First experiment.** CPU (30 min, now): exact Σ core rank by (batch, word, bit) on the frontier table plus the
   linearized 8-bit rounding damage of forcing each plane to rank ≤ 4 (first-order e2 increase with weights
   refit), giving a "mHa² per rank unit" curve and the cheapest 50 planes. P100 (January, ≈ 8 h): warm C 44 → 42
   by merging one leaf pair per group (least-squares projection of the removed leaves, as NR-0236 did for groups),
   compact-group refit to 50k iterations; **go** if continuous e2 ≤ 13.5 mHa, then coded rounding with the cheapest
   planes from the CPU curve. **Kill** if continuous e2 > 14.5 mHa at 50k.
6. **P and effort.** P(≥ −1%) ≈ 0.2, P(≥ −3%) ≈ 0.07; 10 h CPU + 50 h P100.

### F7. Cache amortization: park more address functions in gauge-free word bits

1. **Idea and terms moved.** The per-circuit address cache (≈ 77 ANDs × 15 circuits × 2 passes = 2,310 T, 15.7%
   of T) is rebuilt in every circuit because workspace must return to |0⟩. The only storage that persists between
   circuits is the word registers: under the reachable-state reading (host ruling, topic 753092), a word bit that is
   equal in batch k − 1 and batch k of the same register can carry an address function across circuits. The
   frontier parks 5 such helpers (the top bits, via the π gauge) holding quadratic low monomials. Each further
   address function parked in a register-wide free bit and used by all 15 circuits saves one cache AND per circuit:
   **30 T (−0.20%)**. Moves the cache term only.
2. **Feasibility.** Not accuracy (the gauges are exact) but gauge capacity: a helper needs a bit that can be set
   freely in all 15 columns of that register under the parity conditions of the four word forms (θ, θ + π, −θ,
   π − θ). Bit 7 is taken. The π − θ form flips bit 6 for most words but rewrites bits 0–5 too, which breaks the
   coded bit-0 planes (+rank); NR-0243 found that dead-word and shear freedom does not reach the high-bit planes.
3. **Cost.** 1–3 more helpers: **−0.2% to −0.6%** (293.0M–294.0M); ceiling 5 more: −1%.
4. **Red team.** Rule risk: topic 757522 asks for a written ruling; a strict reading of the ANGLE_k contract would
   void every parked helper at once, the frontier's included. Code conflicts may leave zero usable bits.
   **Verdict: SURVIVES, SMALL** (this is CHEM-A10 with a price and a test).
5. **First experiment (CPU, 2 h).** On the frontier table, enumerate the legal word forms per (row, column) under
   the parity conditions and count, per register and bit position, whether a register-wide helper can be set
   without changing any coded plane. Output: the number of extra parkable helpers. Go if ≥ 2.
6. **P and effort.** P ≈ 0.2; 15 h CPU.

## 3. One-line closures (checked, no family)

| Idea | Why it is closed |
|---|---|
| Classical shortcut: copy the address into classical bits and load the words with classically controlled X (zero Toffolis) | kickmix has no Z-basis measurement into a bit: `CircuitBuilder` offers `hmr` (X basis, uniformly random outcome, phase kickback) and `reset` (dephases a nonzero qubit); `bit_store` writes constants only |
| Condition loader Toffolis on random measured bits so they count ½ on average | A nonlinear product must be present on every shot; Clifford gates under a random classical condition cannot create it and HMR only creates phases, so no gadget computes an AND with expected < 1 counted Toffoli |
| Cheaper reverse pass | The reverse traversal runs the same circuits; every Toffoli runs unconditionally in both passes |
| Loaders that are wrong on rare measurement branches | Excluded: incorrect circuits; the grader's secret salt draws fresh branches and one fault scores the whole submission 0 |
| λ_eff below 46.79 | `2Λ − E_gap ≥ 559.13 + 2 Σ positive eigenvalues of h1'` (triangle inequality on λ_h2) and `E_gap > 3.914898`; the frontier's 46.9216 leaves 0.28% |
| Codebook columns (CHEM-N11: load a 7-bit index, decode to the 8-bit word) | Under the core-rank model the cost is planes × rank, not the alphabet: the index saves one plane per word (≈ 5.5 products per word, 27 per circuit) while a 7 → 8-bit decoder per column costs ≈ 35–60 ANDs (bilinear, 128 entries) and the XOR-delta contract needs the decode of both batches: 350–600 ANDs per circuit, a loss of more than 10× |
| brot 7 | Re-priced to −4.9% to −6.4% (§0.4: Givens −300 T, Q −5, only the cheap coded bit-0 planes disappear) and accuracy-killed (NR-0124: δ7 ≥ 19.3 mHa) |

## 4. Ranked recommendation for January

The honest expectation is small: E[gain] ≈ 0.1 × 2% (F3) + 0.2 × 2% (F6) + 0.2 × 0.4% (F7) ≈ −0.7% chemistry
(+0.23% score). Chemistry cannot close the gap to the leaders; the conceivable maximum (F3 at its counting floor,
−14.9%, plus C 40, −2.1%) is ≈ −17%, and it needs a loader construction nobody knows.

1. **F3 — non-bilinear loaders (proposed CHEM-B19), CPU, can start before the break.** First experiment: SAT-exact
   minimum-AND synthesis on random multi-output functions with don't-cares at (n, m) = (5, 2), (6, 2), (6, 3),
   (7, 2), 20 instances each, against the bilinear low×high count from our core-rank model at the same sizes
   (2–4 h, 12 cores). Go: exact optimum ≥ 15% below bilinear up to n = 7 without closing. Kill: gap < 5% at n = 7.
   Why first: it is the only lever that needs no accuracy budget and no P100, it stacks with every table, and it
   is the one place with a −15% ceiling.
2. **F6 — rank-budgeted C-axis retrain (proposed CHEM-B21), P100 in January.** First experiment on CPU now (30 min):
   per-plane Σ core rank and the linearized rounding price of forcing each plane to rank ≤ 4 on the frontier table
   (the "mHa² per rank unit" curve). Then on the P100 (≈ 8 h): warm C 44 → 42 by merging one leaf pair per group,
   compact-group refit to 50k iterations; go if continuous e2 ≤ 13.5 mHa, kill above 14.5. Why second: the C axis
   is the only shape axis whose value survived the re-pricing, and NR-0140's reopen condition (near-grid start)
   is met by a warm start from the coded C 44 table.

Side item (cheap, CPU, 2 h): F7's helper count; implement only if ≥ 2 extra helpers fit.
Not recommended for P100 time: R 7, R 6, B 129 (F1), brot 7, F4, F5.

## 5. Proposed register updates

New `IDEAS.md` rows (chemistry, B = rewrite with a different technique):

| ID | Idea | Est. gain | P | Effort | Status |
|---|---|---|---:|---|---|
| CHEM-B19 | Non-bilinear (multi-level XAG) ANGLE_k loaders toward the generic counting floor (≈ 195 ANDs per circuit vs 268); supersedes the search framing of A12/B13 | −1.4% to −3.1% (ceiling −14.9%) | 0.10 | 80 h CPU | OPEN — E1 (exact small-case screen) next |
| CHEM-B20 | Block-localized dictionary: 16-row blocks share an orbital region so planes are low-rank by construction | −3% to −12% if feasible | 0.01 | — | NULL (paper, S2.23): one-body rows keep ≥ 5 generic blocks per plane; NR-0249 tails |
| CHEM-B21 | Rank-budgeted C-axis retrain: B05 optimiser with core-rank pricing (0.0057% per unit) and C 44 → 42/40 from a warm coded C 44 start; absorbs B01, B02, B05, B08, B12 as its arms | −1% to −3% | 0.20 | 10 h CPU + 50 h P100 | OPEN — CPU rank-price curve, then P100 in January |
| CHEM-B22 | Dithered twin vectors (rows for precision) to cut 8-bit rounding noise | net +0.6% to +2.3% | 0.02 | — | NULL (paper, S2.23) |

Status changes to existing rows:

- **CHEM-B07** (B 139/129) and **CHEM-B10** (R 7): value cut to −1.6% and −2.0% by the core-rank model (§0.3,
  §0.4); deprioritized below the January bar.
- **CHEM-B11** (brot 7): value −4.9% to −6.4%, not −10%; still accuracy-killed (NR-0124).
- **CHEM-B03** (address permutation): counting bound ≤ about −0.6%.
- **CHEM-B15 / B17** (J): J 4 +2.2% (worse); J 6 ≤ −0.3% even if the cache fitted (it does not, NR-0128).
- **CHEM-A10** (park more helpers): priced by F7 (−0.2% to −0.6%), first experiment defined.
- **CHEM-N11** (codebook columns): NULL (paper, S2.23), §3.

New rows for `experiments/NULLS.md` (chemistry, ours):

| Null | Number that killed it | Source |
|---|---|---|
| Row-axis shape moves at the linear row price | Σ core rank invariant under row removal (8,443 / 8,444 / 8,445 / 8,395 for 1,220 / 1,077 / 934 / 512 rows); only caches scale (1.9 T per row): R 7 −2.0%, R 6 −3.5%, B 129 −1.6% (not −7.5 / −13 / −5.2%) | this document §0.3–0.4 |
| Force the inner QROAM block to 2^3 (Q 427 → 276, T +20.5%, −22.1%) | The model keeps the largest k with temp < 75·brot: needs brot ≤ 3 or R·C·B < 2,432 | §F2 |
| Block-localized dictionary (CHEM-B20) | 76 one-body rows keep ≥ 5 generic blocks per plane (pays only for regions ≤ 7 orbitals); NR-0249 truncation 1.2–2.3 Ha at k ≤ 8 | §F4 |
| Dithered twin vectors (CHEM-B22) | B 159 costs +4.2% for 128 twins; the C budget they free is worth −1.9% (−3.6% at 100% damage concentration) | §F5 |
| Codebook-index loaders (CHEM-N11) | Index saves ≈ 27 products per circuit; per-column decoders for both batches cost ≈ 350–600 ANDs per circuit | §3 |
| Classical copy of the address | No Z-basis measurement in kickmix (only HMR and R) | §3 |
| Random-condition Toffoli halving | No gadget computes an AND with expected < 1 counted Toffoli | §3 |
| Address-permutation gauge alone (CHEM-B03) | ≤ about 110 rank units (−0.6%) by counting (6,970 bits of freedom, 62 bits per rank unit) | §F6 |

Correction for `Analysis/research/chemistry.md` §2.3 (to be applied by the step owner): replace the row-linear
shape rates by the core-rank model (loader T per pass ≈ 0.42 × Σ core rank + K × ⌈rows/16⌉; one row = 1.9 T =
0.013%) and the re-priced table of §0.4; the summary item 7 there ("R 7 (−7%), brot 7 (−10%)") becomes
"R 7 (−2%), brot 7 (−5% to −6%)".

Tooling note for the changelog entry: `core_rank.py` (session scratchpad; 40 lines, NumPy) runs in about 1 s with
< 50 MB on the Windows `.venv`; it should be copied to `tools/` if F6's CPU curve is built on it.
