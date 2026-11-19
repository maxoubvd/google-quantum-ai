# Experiment ledger

One row per attempt on a circuit, a generator or a tool: what changed, the measured result, the verdict. It is
the short index of our work; the details (commands, numbers, analysis) are in the linked `changelog/` entry.
Rows are append-only: a correction is a new row that names the row it corrects.

Before starting an idea, search the known nulls: `experiments/NULLS.md` (community ledger and our research), for
example `grep -i "cutoff" experiments/NULLS.md`. A null is closed only in its design family.

**Verdicts:** `BASELINE` (reference measurement), `KEPT` (the change stays: measured gain or needed tool),
`NULL` (no gain, with the reason), `PARTIAL` (some gain, line closed short), `FAIL` (did not pass a gate),
`WAITING` (a run or a decision is pending).

**Costs** are the official task costs (invsqrt Q × Toffoli-eq; dqi R(ε) × Q² × Toffolis, at the mean over the
salts stated; chemistry λ_eff × T × Q). Score = geometric mean of reference / cost.

| ID | Date | Step | Idea | Task | What changed | Measured result | Gate | Verdict | Changelog |
|---|---|---|---|---|---|---|---|---|---|
| E001 | 2026-11-17 | S1.4 | — | all | Syed Asad Ali's v12 builders rebuilt from source (CPython 3.13.18, numpy 2.2.4) | 17 files byte-identical to the published SHA-256 | build checks | BASELINE | `2026-11-17-2040-s1.4-rebuild-v12-generators.md` |
| E002 | 2026-11-17 | S1.5 | — | invsqrt | v12 `321aff83` | Q 46, 1,511.57 Toffolis, 33 rotations, cost 80,162.73 (full run), worst error/tolerance 0.97358 | official full run + exact proof | BASELINE | `2026-11-17-2102-s1.5-v12-baseline-verification.md`, `2026-11-18-0847-…` |
| E003 | 2026-11-17 | S1.5 | — | dqi | v12 `7cbe8ecf` | Q 889, 383,156 Toffolis; errors on 8 salts 0,2,1,0,0,0,0,0; cost 3.6073e11 clean, mean 3.6132e11 | official, 8 salts | BASELINE | `2026-11-17-2102-s1.5-v12-baseline-verification.md` |
| E004 | 2026-11-17 | S1.5 | — | dfthcsos | v12 (15 loaders, 4,074 Toffolis/pass) | λ_eff 46.9073, T 14,810, Q 427, cost 296,635,967.59 | official | BASELINE | `2026-11-17-2102-s1.5-v12-baseline-verification.md` |
| E005 | 2026-11-17 | S1.6 | — | all | v12 packaged (`v12-baseline.zip`, SHA-256 `231a55f5…`) | projected score 27.078 (clean dqi), 27.049 (4-salt mean), 27.002 (worst salt) | package checks | BASELINE (submitted 2026-11-18, pending) | `2026-11-17-2114-s1.6-build-and-package-tools.md` |
| E006 | 2026-11-18 | S2.1 | GEN-A01 | all | kickmix-engine 0.6.0 as evaluator | identical to the official verifier on every field (dqi errors included); invsqrt 12×, dqi 6–7× (0.55 GB), chemistry 1× | cross-check | KEPT | `2026-11-18-0847-s2.1-s2.2-engine-and-kmb-reader.md` |
| E007 | 2026-11-18 | S2.2 | GEN-A02 | all | `tools/kmb.py` static reader | invsqrt mean 1,511.580 ± 0.008 vs 1,511.570; chemistry 4,074/pass exact; dqi static 384,003; < 1 s | vs verifier | KEPT | `2026-11-18-0847-s2.1-s2.2-engine-and-kmb-reader.md` |
| E008 | 2026-11-18 | S2.3 | GEN-A04, A07 | all | `tools/gate.py` (engine multi-salt, exact proof, cache) | baseline all PASS in 176 s cold / 1.5 s cached; planted invsqrt and dqi errors caught | planted defects | KEPT | `2026-11-18-0857-s2.3-gate-command.md` |
| E009 | 2026-11-18 | S1.7 | — | all | E005 scored on Kaggle (corrects E005's "pending") | public score 27.00209 = local worst-salt projection (2 dqi failed shots on the public split); rank ≈ 1,318 / 3,305 | Kaggle | BASELINE | `2026-11-18-0905-s1.7-s1.8-phase-1-review.md` |
| E010 | 2026-11-18 | S2.5 | DQI-A03 | dqi | Team layer `generators/dqi/team_build.sh` + `src/team/build_team.py` (run-time cutoff patch) at j_min 15 | byte-identical to v12 `7cbe8ecf` | SHA | KEPT (neutral) | `2026-11-18-0950-s2.5-dqi-cutoff.md` |
| E011 | 2026-11-18 | S2.5 | DQI-A03 | dqi | Cutoffs 16 / 17 / 18 built from source | static T −116 per step (383,887 / 383,771 / 383,655); j17 = Özer's public `cc0d7d25` byte for byte | SHA | — | `2026-11-18-0950-s2.5-dqi-cutoff.md` |
| E012 | 2026-11-18 | S2.5 | DQI-A03 | dqi | Targeted test (`tools/dqi_targeted.py`, 2,048 shots per weight 14–20, same inputs) | cutoff j fails all shots of weight < j and nothing else; v12 fails weight 14 | official sampler + simulator | — | `2026-11-18-0950-s2.5-dqi-cutoff.md` |
| E013 | 2026-11-18 | S2.5 | DQI-A03 | dqi | Cutoff choice on 16 salts + Poisson expected-cost model (base failures 0.35/draw pooled) | expected cost j16 −0.026%, **j17 −0.048%**, j18 −0.048% (sd 0.283% vs 0.264%); 16-salt means within noise | 16 salts + official PASS (3.605355e11) | KEPT (j17) | `2026-11-18-0950-s2.5-dqi-cutoff.md` |
| E014 | 2026-11-18 | S2.6 | DQI-A02 | dqi | Assessment of Balram Srinivasan's norm-preserving closing | public form is a compiled gadget (`m_patch.kmx`) transplanted into a finished circuit; generator re-implementation priced at 4–8 h for −1,000 static T (−0.26% task, +0.09% score) | — | DEFERRED | `2026-11-18-1000-s2.6-deferred.md` |
| E015 | 2026-11-18 | S2.7 | INV-A02 | invsqrt | `team/table_search.py`: exact integer model parameterised by the table; coefficient formula D(v) = A + Σ(1−2v_i)((S>>(L−i))+q_i) | reproduces the frozen K = 7 model: 4,096 coefficients, seeds, minimum margin 31,966,138 exactly (4.3 s) | regression | KEPT | `2026-11-18-1011-s2.7-k6-table.md` |
| E016 | 2026-11-18 | S2.7 | INV-A02 | invsqrt | K = 6 table search (64 rows, A 13 bits, S 7 bits), per-row coordinate search on 12 processes (164 s) + exhaustive check of the 12 tightest rows | base angles: 1 failing input (row 45, −3,703,417 units), other rows ≥ 68.9M | exact, all inputs | PARTIAL | `2026-11-18-1011-s2.7-k6-table.md` |
| E017 | 2026-11-18 | S2.7 | INV-A02, INV-A12 | invsqrt | Angle refit (active-set LP, 34 deltas) on the K = 6 table; method validated on K = 7 (124,397,056 vs Syed's 124,370,688) | **0 failures, minimum margin 115,688,300 units, worst error/tolerance 0.920** (6.4 s) | exact, all inputs | KEPT | `2026-11-18-1011-s2.7-k6-table.md` |
| E018 | 2026-11-18 | S2.8 | INV-A02 | invsqrt | `team/loader_synth.py`: bilinear-split XOR-AND synthesis of the 20 output bits (A 13, S 7) of the 64-row table | lower bound 19–20 products + 8–26 monomial ANDs; row-wise construction 44 ANDs (split {0} / {1..5}); rank-greedy finds nothing (rank-1 products rarely lie in the target span alone) | exact over 64 addresses (by construction) | PARTIAL | `2026-11-19-0800-s2.8-loader-synthesis-first-pass.md` |
| E019 | 2026-11-19 | S2.9 | INV-A02, INV-A03 | invsqrt | Map of the build: candidate10 base + 13 patches to reach Q 46 | the width patches rely on facts of the 128-row table (bit identities, measured full-domain facts, frozen loader/decoder certificates); a 64-row table invalidates them; Q 47 costs +2.2% > the table's −1.7% | analysis | PAUSED (Phase 3, own generator) | `2026-11-19-0830-invsqrt-paused.md` |
