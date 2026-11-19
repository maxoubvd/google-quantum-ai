# Public frontier survey, 2026-11-19 (all public work except Asterios Panoras's 28.23162 notebook)

Written 2026-11-19 by a research agent. Scope: the live Kaggle notebooks at the ~28.23 public plateau and those
published since 2026-11-17, excluding `asteriospanoras/qcc-28-23162-what-we-tried-what-failed-and-why` (covered by a
separate report). The forum snapshot in `Instructions/discussion/` was used only as background.

**Sources pulled (read only, nothing executed).** Leaderboard: `downloads/lb/2026-11-19/` (CSV
`...publicleaderboard-2026-11-19T08_27_12.csv`). Notebooks (`kaggle kernels pull -m`) in `downloads/kernels/`:
`maoxin123__qcc-28-23-4-022-toffoli-chemistry-loaders`, `boristown__qcc-gap-3-residue-repair-and-c43-native-compiler__2026-11-19`,
`alexchilton__28-23-lb-which-circuits-survive-private-salt`, `leotang09__vqcc-t7421-submission`,
`ozertuu__28-23163-sota-verified-compilation-top-portfolio`, `ozertuu__quantum-ai-fast-runner-starter-notebook`,
`ghazarosbarseghyan91__vqcc-28-23441-verified-portfolio`, `drpatrickchan__qcc-28-23-frontier-portfolio`,
`junmiaotiao__qcc-28-23-public-frontier-combo-y`, `mizeroluckygall__quantum-circuits-28-23-keep-the-and-frontier`,
`mizeroluckygall__quantum-circuits-28-21-hybrid-explained`, `manvenderraparia__vqcc-28-23163-lp-invsqrt-and-df-thc-refit`,
`ahoirg__qcc-invk6-k6c` (+ `output/result.json`), `atakanaldemir__qcc-shared-products-and-controlled-refits__2026-11-19`,
`varianceofx__a-25-cheaper-dqi-decoder-is-not-a-rewrite`, `aniruddhapawar05__the-df-thc-retrain-ceiling-is-1-32x-not-2x`,
`arizalfirdaus123__vqcc-df-thc-schur-nullspace-refit`, `adityaravi2__qcc-28-18965-lb-guarded-invsqrt-phase-refit`,
`balramsrinivasan__vqcc-cand-28199-submit`, `balramsrinivasan__vqcc-questions-31-5`,
`vumaniacom__qcc-verified-portfolio-and-where-the-gap-to-32-is`, `koushikrudra__qcc-v60-frontier-verified`,
`evgendvorkin__google-quantum-ai-deep-eda-with-analogies`, `issacizrealdavis__vqcc-literal-knob-bench-exact-if-then-search`,
`manwithacat__the-real-150k-was-the-friends-we-made`, `manwithacat__vqcc-hedge-wizardry-notebook`,
`manwithacat__vqcc-negative-results`. Dataset files: `downloads/datasets/balramsrinivasan__vqcc-public-work-1118/README.md`,
`downloads/datasets/manvenderraparia__vqcc-lp-refit-portfolio/{README.md, run_invsqrt_r3_align.py}`; file listings of
`ozertuu/quantum-sim-linalg-cache` and `balramsrinivasan/vqcc-public-work-1118`. `ozertuu/residual-sibling-gating-apex-better`
turned out to be a notebook for another competition (Name That LLM) and is ignored.

SHA-256 of the notebooks whose numbers are quoted most: `ahoirg/qcc-invk6-k6c.ipynb` `e136f05b…`,
`maoxin123/qcc-28-23-…ipynb` `c6b7bae8…`, `varianceofx/…ipynb` `1bc6ae1b…`, `aniruddhapawar05/…ipynb` `f2c3ec2f…`.

All costs below are official-verifier task costs at the default salt as reported by the authors (DQI: clean draw
unless stated). Ratios use the reference costs invsqrt 2,076,774.2337, dqi 24,439,463,583,672.93,
dfthcsos 3,355,599,588.3155.

---

## 1. Board snapshot (2026-11-19 08:27 UTC)

3,377 teams. Medal lines (Kaggle rule for ≥ 1,000 teams: gold = top 10 + 0.2 %, silver = top 5 %, bronze = top 10 %):

| Line | Rank | Score |
|---|---:|---:|
| Gold | 16 | 28.82439 |
| Silver | 168 | 28.23444 |
| Bronze | 337 | 28.23162 |

Top 10 (rank, team, score, submissions; rank on 2026-11-17 21:13 UTC in brackets):

| # | Team | Score | Subs | 17 Nov |
|---:|---|---:|---:|---|
| 1 | Z³C¹oherence | 31.74002 | 54 | 2nd, 31.42723 |
| 2 | ThemBones (was aadi_stonks) | 31.64282 | 36 | 1st, 31.44520 |
| 3 | Rithvik Redrouthu | 31.02835 | 5 | 4th, 30.13063 |
| 4 | Yicheng Guang & Taiseiu | 30.71253 | 104 | 5th, 30.00881 |
| 5 | Entangled Team | 30.28512 | 99 | 3rd, unchanged |
| 6 | Christoffer Thimsen | 30.18947 | 18 | 8th, 29.34887 |
| 7 | Schrödinger's cat lived | 29.97900 | 53 | unchanged |
| 8 | Slingshot Effect | 29.57281 | 43 | unchanged |
| 9 | mliston | 29.33122 | 31 | 23rd, 28.51675 |
| 10 | yu4u | 29.23215 | 9 | 16th, 28.71339 |

Only 12 teams are above 29.0 and 29 above 28.5. **Plateau:** the public best-per-task assembly. Exact-score clusters:
28.23162 (153 teams, ranks 281–433; d0227f3c + 15a4b740 + 4d63be3d), 28.23444 (53 teams; ecf908fe + e39b1396 + 6ef06081,
the current public best), 28.23441 (31; ecf908fe + 15a4b740 + 6ef06081), 28.23438 (11; maoxin123's v5 set),
28.18965 (62), 28.09486 (51), 28.12769 (30), 27.00209 (31; Syed Asad Ali's v12, also our submission E009).
448 teams are at ≥ 28.23, 542 at ≥ 28.20; on 17 Nov only 49 were at ≥ 28.23. The plateau moved
28.1807 → 28.2344 (+0.19 %) in 36 hours, all from public files.

Gap from the public plateau (28.23445) to the leader (31.74002): factor 1.1242 on the score, **1.421 on the
three-task cost product (−29.6 % total cost)**. The top 3 hold unpublished circuits.

---

## 2. Public frontier per task

| Task | Best public file (SHA-256) | Q | Cost | Ratio | Author chain (original builder first, as the NOTICEs require) | Generator public? |
|---|---|---:|---:|---:|---|---|
| invsqrt | `ecf908fec987b61a…` (178,411 B) | 44 | 74,383.70 (mean T 1,450.94, 33 rotations, worst err/tol 0.99946) | 27.9197 | Syed Asad Ali (generator line, v12 46 Q source; v13/v14 44 Q compiled) → Ali Haydar Özdağ (64-row table) → Ertuğrul Özer (angle refit) → Asterios Panoras, James Barlow (measured-CZ / dead-write) → Jeff Emerson, Leo Tang (AND merges) → maoxin123 (44-qubit rescheduled builds 0af21cb8 → 563d0583 → 864f81cd → d0227f3c → **ecf908fe**, 2026-11-19 03:47 UTC) | **No.** The 44-qubit file is embedded as compressed bytes; only the 46-qubit v12 generator (Syed) and the 46-qubit K=6 fork (ahoirg, §3 A2) have source. |
| dqi | `e39b1396fbdf9aae…` (473,123,147 B) | 890 | 345,255,477,425 (mean T 365,892.78; 0/0 errors default salt **and** public board) | 70.7866 | Syed Asad Ali (889 Q v12 source, 890 Q v13 `46446af1` compiled) → Hayato Hikawa (C8 product retention `4b350932`) → boristown (control merging, cached CSWAPs, measured flag uncomputations, product-cache chains; idea credited to busyaprime) → Balram Srinivasan (norm-preserving GF(2^24) closings `66b4f18b`; continued-fraction idea credited to James Day) → Grant Prater (merge, closing at all 49 sites, 4 retained pairs `575b043f`) → boristown (PW4 `5f901c55`) → VUMANIAcom + manwithacat (HW-30, 5 retention patches) → Balram Srinivasan (`40294f00`, converter splice) → Ertuğrul Özer (merge `15a4b740`) → boristown (one pair-window AND retention, **e39b1396**, shipped first by Leo Tang 2026-11-19 04:21 UTC) | **No.** Every step after v12 is a binary patch on a compiled file. Patch *tools* are public (Grant Prater's `vqcc_pair_window.py`, `vqcc_measure_out.py`; Balram's `m_patch.kmx` + `transplant.py`; boristown's `patch_dqi_ccc2ab84.py`), the v13 decoder source is not. |
| dfthcsos | npz `6ef06081cb9c5752…` + 15 loaders (`angle_k_0` `ed648708…`, `angle_k_14` `fe949998…`) | 427 | 294,642,169.51 (λ_eff 46.9215626, T 14,706, 4,022 per pass) | 11.3887 | Lera / marculera (`dfthc-c44-build`: (8,143,44) table, coded words, tree; ideas from ahmedberatozer, somejoint, gengsr, kozykappa, heliosli, Syed Asad Ali, …) → maoxin123 (15 sequential 4,022-Toffoli loaders, CC0; k2..k14 also credited to hayatohikawa) → boristown / Atakan Aldemir (`703843d3`) → Leo Tang (`4d63be3d`) → Manvender Singh (`6ef06081`) | **Partly.** Weight-refit scripts are public (Manvender `run_dfthcsos_v2_opt.py`; boristown `refit_dfthc_v2.py`); Lera's planner/loader builder and maoxin123's loader compiler are not (Lera's dataset ships only the npz; maoxin123 ships compiled loaders). |

Plateau score from these three: **28.234447**. Runner-up files: invsqrt `d0227f3c` 74,405.67 and Manvender's `f7841eb4`
74,405.66; dqi `15a4b740…` 345,256,821,434 (+3.9e-6), `e208d350` 345,252,882,081 at the default salt but **one failed shot on
the public draw** (scores 28.18863); chemistry `4d63be3d` 294,642,170.57, maoxin123's own weights 294,643,137.73.

**Our position for comparison** (E002–E004, v12 from source): invsqrt 80,162.73 (+7.77 % vs frontier), dqi 3.6073e11
(+4.48 %), chemistry 296,635,967.59 (+0.68 %); projected 27.078.

---

## 3. Table A — techniques that worked (since the ~28.18 plateau of 2026-11-17, plus the K=6 source release)

| # | Task | Author | What | Measured gain | Source public? | Licence | IDEAS IDs |
|---|---|---|---|---|---|---|---|
| A1 | invsqrt | maoxin123 | 44-qubit rescheduled builds on Syed's line with Özdağ's table; latest two: the division clean-up hands precomputed AND terms to the second table load, then one more AND removed from that hand-off | 864f81cd 1,453.375 T → d0227f3c 1,451.4375 → ecf908fe 1,450.9375; cost 74,490.95 → 74,383.70 (−0.14 %); since 0af21cb8 (16 Nov, 74,797.37) −0.55 % | No (compiled, embedded) | Apache-2.0 / CC0 by his text | INV-A03, INV-A04, INV-A17, INV-A23 → new **INV-N11** |
| A2 | invsqrt | Ali Haydar Özdağ (ahoirg), notebook `qcc-invk6-k6c`, run 2026-11-19 01:55 UTC | **K = 6 (64-row) table integrated into Syed's v12 generator from source**: `patch_k6.py` (20 kB: own loader certificates for A0..A12, S1..S4; slope word split over the R wire and in-place A11/A12 wires; 6 coefficient rows on an 8-bit window with deferred overflow; own decoder and table clean-up, no T1/T2 replay, no paid-R block), `hold_pass_k6.py`, `k6data.json`, endpoint angles refit by LP, fork of 234 pinned files; exact proof on all 3,133,452 x and all HMR transcripts | Q 46, raw 1,491.75 → final 1,480.75 Toffolis, 33 rotations, **78,765.64 vs 80,162.73 (−1.74 %)**, worst err/tol 0.9494, official PASS (200,540,928 shots), SHA `4560ac0c…` | **Yes, full source** (base64 in the notebook; read, not executed) | Kaggle public code (rules require an OSI licence); patch text says it adapts Syed's Apache-2.0 layer | **INV-A02** (exactly our paused S2.9/S2.10), INV-A12 |
| A3 | invsqrt | Ertuğrul Özer; Manvender Singh | Angle-only refits: 33-bit activation mask, cutting-plane HiGHS minimax LP, exact 128-bit integer refinement | margin 38,858,846 → 38,858,848 units; cost 74,405.67 → 74,405.66 (the 0.01 is the seed grinding of B4, not the LP) | Yes (`run_invsqrt_r3_align.py`, Apache-2.0) | Apache-2.0 | INV-A12 (angle part) |
| A4 | dqi | Balram Srinivasan (`40294f00`) | Syed's v13 decoder with a converter ("CSWAP converter") and phase-kick splice taken from an earlier variant of the same decoder line | 367,110 → 366,725 T (−385, −0.10 %); 3.46041e11 | No (author: "exact gadget set not independently confirmed") | Apache-2.0 + Syed NOTICE | DQI-A01 (partial) → new **DQI-N11** |
| A5 | dqi | Ertuğrul Özer (merge) on Balram, Grant Prater, boristown, VUMANIAcom/manwithacat | `15a4b740`: 32 more GF(2^24) closings (−672 T), 8 measured uncomputations (−8), PW4 retention (−44), HW-30 retention (−116) on `40294f00` | 366,725 → 365,894.20 T (−0.23 %); 3.452568e11, clean on default and board | Compiled; patch tools public (see §2) | Apache-2.0 | DQI-A02, DQI-A05, DQI-A06, DQI-A10 |
| A6 | dqi | boristown | Pair-window AND retention: a repeated `CCX(a,b,t)` with a, b unchanged in between is parked on a wire proven \|0⟩ (last touched by an unconditional HMR, untouched in the window) and erased by HMR + conditional CZ. Scanner checked all 281 repeated CCX pairs of the 57 M-op file: **1 site** (`CCX(q880,q14,q50)`, ops 39,436,586 / 39,462,171, free wire q889) | −1.42 mean T → 365,892.78; cost 345,255,477,425 (−3.9e-6); 0 errors default, board, and (Pablo Diego Pereira) salt private-1 at 25,600 shots | Method described; file compiled | Apache-2.0 | DQI-A06, GEN-A03 → refinement **DQI-N13** |
| A7 | dqi | Koushik Rudra (v60) | Claims 18 further pair-window sites on `e208d350` | not reported (no outputs in the public version) | No | — | DQI-A06 (unverified) |
| A8 | dqi | Atakan Aldemir | Commutes each paid control product into a CZ conditioned on the ending X measurement plus a CCZ conditioned on an earlier measurement (measurement count unchanged) | −6.74 / −6.85 mean T on his V11 DQI (two clean 102,400-shot runs, GM 28.18976 / 28.18975); **not in the frontier file** | Worker code embedded (zlib) | Apache-2.0 | DQI-A14 (close) → new **DQI-N12** |
| A9 | dfthcsos | maoxin123 | 15 sequential ANGLE loaders on XOR differences of consecutive word batches; k0 292 and k1 287 Toffolis; CNOT layers −64 % | 4,022 per pass (Syed v12: 4,074; T 14,810 → 14,706); frontier since 17 Nov | No (compiled loaders) | CC0 | CHEM-A09, CHEM-A03 |
| A10 | dfthcsos | Leo Tang, Manvender Singh, maoxin123, Arizal Muluk | Weight refits at fixed words: Schur-nullspace reduced-gradient L-BFGS with feasibility 5e-12 (2,063 steps, Manvender), corank-2 reduced gradient keeping both error gates and the gap within 2.65e-10 Ha (maoxin123), PCG projected Gauss-Newton + L1 epigraph LP (Arizal) | 703843d3 294,642,299.76 → 4d63be3d 294,642,170.57 → 6ef06081 294,642,169.51 (**−4.4e-7 relative**); λ_eff 46.92158 → 46.9215626 | Manvender's scripts yes | Apache-2.0 | CHEM-A13 (saturated) |
| A11 | packaging | Leo Tang, boristown, Patrick Chan, Manvender Singh | Choose among near-equal DQI files by their **public-board draw** (e208d350 fails one shot there, e39b1396 is clean) | public score 28.18863 vs 28.23444; no effect on the private expectation | — | — | GEN-A04, DQI-A11 |

---

## 4. Table B — nulls reported (with numbers)

| # | Author | Task | Null | Numbers | Scope | IDEAS / NULLS |
|---|---|---|---|---|---|---|
| B1 | varianceofx (*A 25 % cheaper DQI decoder is not a rewrite*, on `61c8348d`) | dqi | (a) R is not a lever | R ≥ 1.1913 at 0 failures (δ = 2^-128, n = 102,400); perfect reliability is worth at most ×1.06 on the task | general (verifier) | already in `dqi.md` §4 (R floor) |
| | | | (b) Cheaper GF(2^12) multiplier | hot block = truncated multiply-accumulate in **42 ANDs**; schoolbook 144, Karatsuba 81, tower 54 / 45; bilinear bound 23 unattained for n = 12 | Syed decoder line | consistent with `dqi.md` §4; weakens DQI-B08 only indirectly |
| | | | (c) Dead / affine Toffolis | exact ANF of the 42-, 42- and 173-CCX templates: degree 2 with 12 / 12 / 48 quadratic outputs; no affine block | family-specific | DQI-A15, GEN-B06 weakened |
| | | | (d) Global compute/uncompute cancellation | 362,016 CCX; 43,800 identical CCX on an unchanged target, **0** with both controls untouched; only **875 of 367,993** Toffolis (0.24 %) are skipped per shot | family-specific | DQI-A07, DQI-A17 (little conditional mass to save) |
| | | | (e) Predicate the 11-AND Chien zero test | `push_condition` takes a classical bit and there is no qubit→bit copy; error weight pinned (mean 31.05, p99 35, P(w ≤ 20) = 5.5e-5) | general (ISA + weight law) | DQI-B03, DQI-B09 (scan) must find another route |
| | | | Census | 100 % of CCX touch three distinct 12-qubit words; 1,392 whole-multiply triples hold 124,698 ANDs (34.4 %), the rest are partial products of 1–4 ANDs; to gain 1.34× on the score from DQI alone: T 367,118 → 274,434, i.e. **≤ 5,462 GF(2^12) multiplications instead of ~7,668** (−29 %) | — | quantitative target for DQI-B02/B05/B09 |
| B2 | Aniruddha Pawar (*The DF-THC retrain ceiling is 1.32x, not 2x*) | dfthcsos | Retrain to a cheaper shape | accuracy rank is **R·min(B, C)**, not R·C: (1,129,352) has rank 129, Eckart-Young floor 292 mHa vs 20 cap; rank 290 floor = 19.86 mHa; best feasible shape (3,129,118); whole retrain worth **≤ 1.32×** on chemistry (≈ 29.9 overall) with public-quality loaders and **0.98× (a loss)** with textbook QROAM loaders (495 vs 268 Toffolis per circuit; starter kit ≈ 1,219); Q is not a lever (k2 = 4 for every B ≥ 129, B < 129 → k2 = 5 and Q ≈ 2×); λ_eff converged (three refits within 4e-7) | general for the chemistry cost model | refines NR-0020; CHEM-B04/B07/B08/B10/B12 keep B ≫ C so rank = R·C still holds for them, but each must keep R·C ≥ ~290–352 and B ≥ 129 |
| B3 | Balram Srinivasan (dataset `vqcc-public-work-1118`, notebook `vqcc-questions-31-5`) | dfthcsos | brot 7 table (9,127,44) | continuous-fit floor error_2 15.21 mHa, after rounding **27.14–28.78 mHa (FAIL)**, error_1 12.1–12.3 PASS; λ_eff 48.4 (8-bit 46.92); dead ends: progressive rounding, single-word flips + refit (28.52 → 28.43), exact block-optimal weights (+0.07), copy-count sweeps C 44 → 127 (26.5), twin replacement (floors 21.6 / 25.4); loaders shrink only **4,022 → 3,876 per pass from 8 to 6 bits (−3.6 %)**; code includes `qat.py`, `adaround.py` | general for low-bit words at this table size | **CHEM-B11** NULL; CHEM-B01 value cut |
| B4 | Manvender Singh (`run_invsqrt_r3_align.py`) | invsqrt | "Skeleton seed alignment": flip low bits of one angle word (delta 0..N) to pick the file whose **default-salt** SHA-256 seed gives the fewest measured Toffolis | 1,451.4375 → 1,451.4359 mean T, cost 74,405.67 → 74,405.662; the gain is a property of the empty salt and **does not transfer** to the hidden salt | general (any HMR-dependent count) | GEN-A04: rank candidates on expected counts (our `tools/kmb.py` exact model) or multi-salt means, never on one draw |
| B5 | VUMANIAcom, alex chilton, Patrick Chan, Manvender Singh | dqi | Salt robustness of the patched 890 Q line | 2683225b: 1 failed shot at 1 of 5 salts (≈ 3/262,144 per shot, ≈ 20 % risk per unseen salt, −0.043 GM); f425ea20 fails 4 of 8 salts, 66b4f18b 1 of 8; e208d350 fails the public draw; each failure +0.46 % cost (−0.15 % GM) | general for this line | GEN-A04, DQI-A07, DQI-A11, DQI-A17 |
| B6 | boristown | dqi | Pair-window retention at fixed width is nearly exhausted on the frontier file | 281 repeated CCX pairs, **1** window with a provably free wire | family-specific (890 Q file) | DQI-A06: little left at 890 Q unless width changes |
| B7 | alex chilton (28.23 salt notebook) | all | Observation | "The last two weeks of public gains are refits and single-gate rewrites worth 0.01 to 0.1 % each. No public notebook … found a design change worth 2× on any task." | general | — |
| B8 | Issac Davis (*Literal Knob Bench*) | invsqrt | Precision-knob cube (d13, sq27, p27, r28) on Google's wide Newton-Raphson reference | 16 exhaustive rows; transfer to the 44 Q family explicitly unknown | reference family only | not applicable to our line |
| B9 | Atakan Aldemir | dqi | Previously captured rare decoder failures remain after his control-product rewrite | documented, no new numbers | family-specific | DQI-A07 |

---

## 5. Table C — new ideas not in our register

| Proposed ID | Task | Idea | Evidence | Expected gain | Notes |
|---|---|---|---|---|---|
| **INV-N11** | invsqrt | Cross-stage AND hand-off: keep the AND terms computed during the divider clean-up (DIV_CLEAN) and feed them to the second table load (LOAD2) instead of recomputing them; then minimise the hand-off | maoxin123 d0227f3c / ecf908fe: −2.44 mean T vs 864f81cd (−0.14 % cost), compiled only | 0.1–0.2 % task cost | Neighbour of INV-A23 (share LOAD1/LOAD2 logic) and INV-A17 (clean-up premium); only meaningful after INV-A02 (+A03) since the hand-off depends on the table address bits |
| **DQI-N11** | dqi | Stage-wise best-of splicing between decoder variants of the same line (e.g. v12 / v13 / earlier variants): take each stage (converter, phase kick, …) from the variant where it is cheapest, when stage boundaries carry the same registers | Balram Srinivasan `40294f00`: −385 T (−0.10 %) from a converter + phase-kick splice | 0.1–0.3 % | For us: our v12 generator can emit variants (cutoff, flags); a source-level version of this is part of DQI-A01 |
| **DQI-N12** | dqi | Commute paid control products into measurement-conditioned phase fix-ups: CZ conditioned on the ending X measurement plus a CCZ conditioned on an earlier measurement | Atakan Aldemir: −6.7 to −6.9 mean T on V11, two clean runs | ≈ 0.002 % per site group; cheap once a matcher exists | Close to DQI-A14 (phase-debt merging); needs GEN-A10 local proof |
| **DQI-N13** | dqi / cross-cutting | Free-wire certification rule for retention: a wire is free over a window if its last touch is an unconditional HMR and nothing uses it in between; scan all repeated CCX pairs for such windows | boristown: 1 site in 281 on e39b1396; Koushik Rudra claims 18 on e208d350 | small (≈ 1–20 T) | Refinement of DQI-A06 / GEN-A03; trivial to add to a `tools/` scanner over our own outputs |
| **CHEM-N11** | dfthcsos | Codebook columns: many of the 75 angle columns take only 65 or ~129 distinct 8-bit values (effective 6–7 bits); load a short index and decode it with a cheap reversible map | Aniruddha Pawar measurement on the public table | unknown, ≤ 3.6 % of loader T if it behaves like brot 6 (B3) | First check whether the bit-0/bit-1 affine codes already exploit exactly this structure (likely partly) |

Ideas seen but already in the register: gauge freedom of the word table (u_b sign, permutation of the B vectors in a
block, choice of `tree_edges`; Aniruddha Pawar) = CHEM-B03 + CHEM-B09 (+ gengsr's four word forms, already used);
board-draw file selection = GEN-A04 / DQI-A11 (no private value); daily public-circuit manifests (manwithacat
*Friends*, Leo Tang T-series, Patrick Chan) = GEN-A06 feeds.

---

## 6. Different circuits, or small techniques on the same few families?

**Answer: an accumulation of small techniques on three generator families. No public family change since about
2026-11-15/16 on any task.**

Evidence from the plateau:
- All 53 + 31 + 11 + 153 teams in the four top public clusters (28.2316–28.2344) ship byte-identical members from the
  same three lines: invsqrt from Syed Asad Ali's 44/46-qubit generator line, dqi from Syed's 890-qubit decoder line,
  chemistry from Lera's (8,143,44) coded table with maoxin123's loaders. Every portfolio notebook read (Leo Tang,
  boristown, Mizero Lucky Gall, Patrick Chan, junmiaotiao, Ghazaros Barseghyan, ozertuu, Manvender Singh, alex
  chilton, Koushik Rudra, Arizal Muluk) says "no new circuit design" or adds one patch.
- Per-task movement 16 → 19 Nov: invsqrt 74,797.37 → 74,383.70 (**−0.55 %**, Q 44 and 33 rotations unchanged, 8.4 fewer
  mean Toffolis); dqi 367,120 → 365,892.78 mean T (**−0.33 %**, Q 890 unchanged); chemistry 294,762,512.86 →
  294,642,169.51 (**−0.04 %**, Q 427 and T 14,706 unchanged since 17 Nov). The score moved 28.12 → 28.23.

Lineage of each frontier file back to its family, with the real family changes (dates approximate, from notebook run
times and leaderboard clusters; score milestones from Evgen Dvorkin's table as re-told by Mizero Lucky Gall):

**invsqrt** (`ecf908fe`)
1. Google reference: wide Newton-Raphson, cost 2,076,774.
2. Family change 1 — table + polynomial/Horner lines with signed digits and borrowed carries (Sumit Chongder,
   anazemcev 93 Q, Rudra Bhaskar 95 Q, Kirill Tymchishin's integration; costs not recorded here).
3. Family change 2 — compact QROM-seed + Newton lines at 63–65 Q (mrmeows, busyaprime), then 62 Q (yoyobar, LB 22.99,
   ~9 Nov); score milestones 22.70 → 25.45 at 62–63 Q.
4. **Family change 3 — Syed Asad Ali's v12 "measure then rebuild" 46-qubit generator** (public by 14 Nov at the latest,
   when ~70 teams sat at 27.31; the 25.45 → 27.08 milestone; ≈ −14 % task cost against the 62 Q line, estimated from
   62 × 1,509 Toffoli-eq; source public, our base).
5. **Width step 46 → 44** (Syed's v13/v14 `47b31838`, maoxin123's v27 `0af21cb8`, 15–16 Nov; −6.7 %; keeps root bits
   v10/v11 and exponent flags off the multiplier). Same family, compiled only.
6. Özdağ's 64-row table (−1.74 % alone at 46 Q), dead-write CZ (−1 T), AND merges (−2 T), angle LPs, the DIV_CLEAN →
   LOAD2 hand-off: −0.55 % in total since 16 Nov.

**dqi** (`e39b1396`)
1. Google reference, cost 2.44e13.
2. Family change 1 — textbook decoders (Berlekamp-Massey port of Gu et al., arXiv:2610.22367, 1,837 → 1,789 Q with
   kozykappa's gadgets; synchronized EEA of Khattar et al. at 1,938 Q).
3. **Family change 2 — continued-fraction / Schur (CRT-Thiele) in-place decoder** (jday96314 "QuantumAI" round 6,
   reverse synthesis at 883 Q): roughly halves Q, so ~4× on Q²; the later lines all descend from this idea.
4. SomeJoint's 892 Q decoder with flag carry; busyaprime's exact product caching (−0.65 % there).
5. **Syed Asad Ali's v12 generator, 889 Q, 383k T** (source public, our base) → **v13, 890 Q, ~368.7k T (−3.6 %)**,
   compiled only. Same family; the v13 delta is the last large step.
6. Patch era on the v13 file (retention, measured uncomputation, GF(2^24) closings, converter splice, pair-window
   retention): 368.7k → 365.9k T (−0.76 %).

**dfthcsos** (`6ef06081` + maoxin123's loaders)
1. Google reference: plain unary-iteration QROM loaders (~1,219 Toffolis per circuit, ~18,285 per pass), cost 3.36e9.
2. Family change 1 — sparse SOS factors and banked two-address XOR QROM (anazemcev, kozykappa).
3. **Family change 2 — Lera's compact (9,143,40/42) tables with a 128-qubit loader emitter** (~16k T; manwithacat's
   kernel-relation rebuild −1.8 %), five words per circuit, borrowed top bit (heliosli).
4. **Family change 3 — the (8,143,44) coded-word table** (somejoint's bit-0 affine codes, Lera's codes on the XOR of
   five words, maoxin123's re-coded tables v25–v27, ~13–14 Nov): T 14,840 → 14,712; Q 427.
5. Loader polishing 4,074 (Syed v12) → 4,038 → 4,033 → 4,030 → 4,028 → 4,022 per pass, then weight refits that now
   agree to 4e-7. The (R, B, C, brot) shape has not changed since step 4.

**What it means for the leaders.** The top 3 (31.03–31.74) need a 1.32–1.42× cheaper cost product than the plateau.
Aniruddha Pawar bounds a chemistry retrain at 1.32× even with public-quality loaders, Balram Srinivasan shows brot 7
failing at 27 mHa, varianceofx shows DQI needs ≤ 5,462 field multiplications (−29 %) for a 1.34× step, and
VUMANIAcom's split puts 32 out of reach of any single task's micro-levers. So the leaders almost surely run **new
families on at least two tasks** (or a new DQI key-equation algorithm plus a new chemistry loader architecture),
none of which is public. The public record has no evidence of what they are.

---

## 7. Proposed status per affected IDEAS ID

Status values as in `Analysis/IDEAS.md`; "for us" is the from-source status under decision S1.3 (track B).

| ID | Proposed status | Reason |
|---|---|---|
| INV-A01 | EXCLUDED (unchanged) | frontier file is now `ecf908fe` (74,383.70), not `563d0583` |
| **INV-A02** | **DONE (public, by Ali Haydar Özdağ, from source: 46 Q, 78,765.64, −1.74 %)**; for us STILL OPEN → unblocked | `ahoirg/qcc-invk6-k6c` publishes the whole K=6 patch layer on Syed's v12 generator (`patch_k6.py`, `hold_pass_k6.py`, `k6data.json`, LP angles, pipeline) with exact proof; it solves the blocker recorded in S2.9 (table-specific width patches). Reuse with credit (Syed Asad Ali first, then ahoirg) or use as the reference for our own patch; our S2.7 table/angles can be compared row by row |
| INV-A03 | STILL OPEN | 46 → 44 moves exist only in compiled files (Syed v13/v14, maoxin123) |
| INV-A04 | STILL OPEN | public micro-edits are compiled; total value since 44 Q ≈ −0.55 % |
| INV-A12 | angle-only part DONE (public, by Ertuğrul Özer and Manvender Singh: margin saturated at 38,858,848 units); joint table-rounding part STILL OPEN | |
| INV-A17, INV-A23 | STILL OPEN (public evidence of value: maoxin123's hand-off, see INV-N11) | |
| DQI-A01 | STILL OPEN | v13 source still not public; Balram's `40294f00` shows variant splices of the line pay (−385 T) |
| DQI-A02 | DONE (public, by Balram Srinivasan, Grant Prater, Ertuğrul Özer, compiled gadget transplant at 41 + 8 + 32 sites); for us DEFERRED (unchanged) | |
| DQI-A05 | DONE (public, compiled: boristown, Grant Prater, ozertuu MBUC); for us STILL OPEN | |
| DQI-A06 | DONE (public, compiled: Hikawa, boristown PW3/PW4/pair-window, Grant Prater, VUMANIAcom HW-30); nearly exhausted at 890 Q (1 free window in 281 pairs); for us STILL OPEN from source | |
| DQI-A07, DQI-A11, DQI-A17 | STILL OPEN | failure rate now quantified (~3/262,144 per shot; 1 failure = +0.46 % cost); only 0.24 % of Toffolis are conditionally skipped, so a guard cannot buy much T |
| DQI-A14 | STILL OPEN (positive public evidence: Atakan Aldemir −6.8 T, see DQI-N12) | |
| DQI-A15 | STILL OPEN, weakened (public: varianceofx — heavy templates are degree-2, no affine blocks) | re-synthesis to fewer ANDs is not excluded, deletion is |
| DQI-A13, DQI-B08 | STILL OPEN | varianceofx covers the GF(2^12) multiply (42 ANDs), not the 66-T inverse or the GF(2^24) rank |
| DQI-B02, DQI-B05, DQI-B09 | STILL OPEN; target set | to matter at the top: ≤ 5,462 GF(2^12) multiplications (≈ 274k T at 890 Q) |
| DQI-B03 | STILL OPEN | predicating the Chien zero test is not expressible (B1e); a fused finish must avoid a qubit-valued condition |
| CHEM-A01 | EXCLUDED (unchanged) | frontier weights now `6ef06081` |
| CHEM-A03 | STILL OPEN | Lera's planner and maoxin123's compiler are not public |
| CHEM-A09 | DONE (public, by maoxin123: k0 292, k1 287, compiled); for us STILL OPEN | |
| **CHEM-A13** | **NULL (public: saturated — Leo Tang, Manvender Singh, maoxin123, Arizal Muluk agree to 4e-7; last step −4.4e-7)** | |
| CHEM-B01 | STILL OPEN, value cut | loaders shrink only 3.6 % from 8 to 6 bits (Balram) |
| CHEM-B03, CHEM-B09 | STILL OPEN | gauge freedom confirmed exact (Aniruddha Pawar), untried publicly |
| CHEM-B04, B07, B08, B10, B12 | STILL OPEN, with the rank rule R·min(B,C) ≥ ~290 (352 shipped) and B ≥ 129 as hard filters; a full retrain is capped at 1.32× and needs frontier-quality loaders | Aniruddha Pawar |
| **CHEM-B11** | **NULL (public: Balram Srinivasan — brot 7 (9,127,44) rounds to error_2 27.1 mHa > 20 even from a 15.2 mHa continuous floor; QAT/AdaRound code tried)** | reopen only with a new rounding idea that keeps error_2 ≤ 20 |
| GEN-A04 | STILL OPEN (in progress); reinforced | seed grinding (B4) and board-draw selection (A11) show why ranking must use expected or multi-salt costs |
| GEN-A06 | STILL OPEN | ready-made feeds: manwithacat's *Friends* manifest, Leo Tang's lineage headers, alex chilton's salt table |
| NR-0020 (NULLS) | refine | rank is R·min(B, C), not R·C (matters only for shapes with B < C) |
