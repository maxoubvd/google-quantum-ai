# Idea register

Created 2026-11-18 (plan step S1.9). This file merges the four research reports in `Analysis/research/`
(`invsqrt.md`, `dqi.md`, `chemistry.md`, `general.md`, written on 2026-11-17/18 by four research agents from
our own circuits, the papers, the forum, the live public notebooks and the web) into one prioritized list.
The research files hold the details of every idea (mechanism, risks, sources), the cost anatomy of each of our
circuits and the known nulls; this file holds the ranking and the status. Phases 2 to 4 of `Analysis/PLAN.md`
draw their milestones from it.

## 1. How to read it

- **Categories.** `A` = improve our code (incremental changes to the v12 generators, ports of public
  improvements into our source, tooling); `B` = rewrite with a different technique (new algorithm,
  architecture, factorization regime or pipeline). IDs: `INV-`, `DQI-`, `CHEM-` (task ideas), `GEN-`
  (cross-cutting methods and tooling).
- **Gain as stated** is the research report's own estimate: % of the **task** cost for INV/DQI/CHEM ideas, %
  of the **score** (or an "equivalent score %" for tooling, a stated judgement) for GEN ideas.
- **Score gain %** puts every idea on one scale: for a task idea, `100 × ((1/(1 − g/100))^(1/3) − 1)` ≈ g/3
  (the score is the geometric mean of three ratios).
- **P** is the estimated probability that the idea works as described. **Effort** is in hours of focused work,
  including the exact proof and a verifier run.
- **Priority = P × score gain % ÷ effort** (user decision, 2026-11-17: rank by expected value per hour, show P).
  Every table is sorted by this priority, highest first. The numbers are estimates made before any build;
  they are revised when an idea is tried (the status column and the changelog record the measured result).
- **Status:** `TODO`, `IN PROGRESS`, `DONE (measured gain)`, `NULL (measured, with reason)`, `EXCLUDED
  (decision)`, `MERGED (into another idea)`.
- Regenerate the tables of sections 3 to 7 after a research file changes with
  `python tools/ideas_rank.py --md <file>` (then paste them here, keeping the statuses).

## 2. What the research says, in short

- **Where the gap is.** Our base (Syed Asad Ali's v12 generators, score ≈ 27.06) is about 4% under the public
  plateau (28.18) and 15% under the leaders (31.4). The public frontier's advance over v12 lives only in
  compiled files (invsqrt 46 → 44 qubits, −7.0% task cost; dqi 889 → 890 qubits and −16k Toffolis, −4%;
  chemistry 4,074 → 4,022 Toffolis per pass, −0.7%). Rebuilding those changes **in our own source** is the
  safest gain available (INV-A02/A03/A04, DQI-A01/A02/A04, CHEM-A03).
- **Tooling multiplies everything else** and ranks highest: a cross-checked fast simulator (GEN-A01,
  kickmix-engine, 5–9× author-reported), a nulls ledger (GEN-A05), our own static circuit reader (GEN-A02),
  frontier watching (GEN-A06) and a multi-salt gate (GEN-A04).
- **Beyond the frontier,** each circuit is locally tight (0 removable Toffolis found by several scans), so the
  incremental ideas are worth 0.1–1% each. Large steps need a new algorithm family (INV-B*, DQI-B02/B04,
  CHEM-B01/B10/B11), found most plausibly by a divergent brainstorm plus an evolutionary agent loop over the
  generators (GEN-B11, GEN-B01, GEN-B02), as in Gu et al. (arXiv:2610.22367).
- **Leaders' gap:** 31.4 needs about a 28% cost cut on one task (e.g. invsqrt near 32 qubits, or dqi near
  264k Toffolis). Chemistry is boxed in by its cost model (λ_eff 0.25% above a hard floor; realistic total
  −2.5 to −3.5%).
- **Decision S1.3 (track B only)** excludes submitting compiled third-party files: INV-A01 and CHEM-A01 are kept
  for the record as `EXCLUDED`; their value is carried by the from-source ports.

## 3. Global ranking (top 30, all tasks)

| ID | Task | Idea | Gain as stated | Score gain % | P | Effort h | Priority | Status |
|---|---|---|---|---:|---:|---:|---:|---|
| GEN-A01 | cross-cutting | kickmix-engine as our fast evaluator | 7–9× faster verification (author-reported), dqi memory 9 GB → < 1 GB so 4+ workers fit in WSL; ≈ 0.5% eq. | 0.50 | 0.85 | 3 | 0.1417 | DONE (S2.1: identical results; dqi 6–7×, invsqrt 12×, chemistry 1×) |
| GEN-A05 | cross-cutting | Nulls ledger seeded from the community | Avoids repeating closed work; ≈ 0.3% eq. | 0.30 | 0.95 | 3 | 0.0950 | DONE (S2.4: `experiments/NULLS.md`, `experiments/LEDGER.md`) |
| GEN-A02 | cross-cutting | Static IR + exact static cost | Evaluator tier T0 in seconds instead of minutes; ≈ 0.4% eq. | 0.40 | 0.90 | 6 | 0.0600 | DONE (S2.2: `tools/kmb.py`, < 1 s per circuit) |
| GEN-A06 | cross-cutting | Daily frontier-ingestion agent | Keeps us on the moving plateau without manual reading; ≈ 0.6% eq. | 0.60 | 0.80 | 8 | 0.0600 | TODO |
| GEN-A04 | cross-cutting | Salt-robustness gate for every HMR-dependent edit | Avoids −0.15% per dqi failed shot and the 0-score case; ≈ 0.3% eq. (risk-weighted) | 0.30 | 0.90 | 5 | 0.0540 | IN PROGRESS (S2.3: multi-salt gate in `tools/gate.py`; targeted dqi inputs in M2.3) |
| GEN-A07 | cross-cutting | I/O and caching hygiene | 1.5–3× on load-heavy dqi runs, no repeated runs; ≈ 0.1% eq. | 0.10 | 0.90 | 2 | 0.0450 | DONE (S2.3: work on ~/work ext4, result cache keyed by SHA and salt) |
| GEN-B11 | cross-cutting | Divergent-family brainstorm + red team, weekly | ≈ 1% (one new family that beats the incumbent on one stage) | 1.00 | 0.25 | 8 | 0.0312 | TODO |
| INV-A02 | invsqrt | Port Ali Haydar Özdağ's 64-row coefficient table into our generator | 1.5–2.0 | 0.59 | 0.75 | 16 | 0.0277 | TODO (effort re-estimated at 24–40 h on 18 Nov: no table or loader search tools in the v12 source; steps S2.7–S2.10) |
| GEN-B01 | cross-cutting | Evolutionary agent loop over our generators | 2–5% over the remaining weeks (Gu et al. found 4–39× over **literature** baselines; ours are crowd-optimized, so far less) ≈ 3% | 3.00 | 0.35 | 40 | 0.0262 | TODO |
| INV-A03 | invsqrt | Port Syed's v13/v14 46 → 44 moves into the source | 4.3–4.8 | 1.56 | 0.50 | 30 | 0.0261 | TODO |
| GEN-B02 | cross-cutting | Principle sweep, stage by stage | ≈ 1.5% (a few stage-level wins) | 1.50 | 0.30 | 20 | 0.0225 | TODO |
| INV-A04 | invsqrt | Public micro-edit pack as generator passes | 0.5–0.9 | 0.23 | 0.70 | 8 | 0.0205 | TODO |
| GEN-A09 | cross-cutting | Kaggle batch-verification orchestrator | Frees the laptop; 48-salt dqi sweeps and full invsqrt runs in parallel; ≈ 0.15% eq. | 0.15 | 0.80 | 6 | 0.0200 | TODO |
| CHEM-A02 | chemistry | Daily "parts market" watch: pull new public chemistry loaders and weights, re-verify, adopt the cheapest compatible pair | 0.3 (0.1-0.7) | 0.10 | 0.70 | 4 | 0.0175 | TODO |
| DQI-A01 | dqi | Recover Syed's v13 delta and port it into the v12 generator | 3.0 to 3.6 | 1.12 | 0.55 | 40 | 0.0155 | TODO |
| GEN-B04 | cross-cutting | Qubit-cap ladder harness | invsqrt −1 qubit = 0.77% | 0.77 | 0.20 | 10 | 0.0154 | TODO |
| GEN-A03 | cross-cutting | Task-agnostic finishing pipeline | 0.01–0.1% per task per round; ≈ 0.25% | 0.25 | 0.60 | 12 | 0.0125 | TODO |
| GEN-A10 | cross-cutting | Local proof per rewrite | Catches rare-branch phase bugs early, fewer wasted full runs; ≈ 0.2% eq. | 0.20 | 0.60 | 10 | 0.0120 | TODO |
| GEN-A12 | cross-cutting | Design-space optimizer over generator parameters | ≈ 0.3% from better trade-offs | 0.30 | 0.40 | 10 | 0.0120 | TODO |
| DQI-A02 | dqi | Norm-preserving GF(2^24) closing in `gf24_nested_norm.emit_closing` | 0.21 to 0.26 | 0.08 | 0.85 | 6 | 0.0111 | TODO |
| DQI-A03 | dqi | Extractor cutoff j_min 15 -> 18 | 0.05 to 0.08 net | 0.02 | 0.95 | 2 | 0.0103 | DONE (S2.5: cutoff 17 kept, expected −0.048% task cost; measured, not the estimated −0.05 to −0.08%) |
| GEN-A08 | cross-cutting | Peak map + functional-dependency mining | ≈ 0.5% (invsqrt qubit 0.77%, dqi 0.075%) | 0.50 | 0.20 | 10 | 0.0100 | TODO |
| GEN-B12 | cross-cutting | One end-to-end error budget, fitted to the exact quantized computation | ≈ 0.4% (one table bit or a few rotations) | 0.40 | 0.25 | 12 | 0.0083 | TODO |
| INV-A05 | invsqrt | Forward normaliser with measurement-conditioned clean-up, like the inverse | 0.2–0.5 | 0.12 | 0.40 | 6 | 0.0078 | TODO |
| GEN-B03 | cross-cutting | Width-optimal rescheduling by spooky pebbling DP/SAT | ≈ 0.77% (one invsqrt qubit; open best case −0.37% cost) | 0.77 | 0.20 | 20 | 0.0077 | TODO |
| INV-A10 | invsqrt | Clean up rotation bookkeeping | 0.01–0.06 | 0.01 | 0.60 | 1 | 0.0070 | TODO |
| GEN-A13 | cross-cutting | Scaled-down proxy instances | Faster triage of new families; ≈ 0.2% eq. | 0.20 | 0.40 | 12 | 0.0067 | TODO |
| INV-A08 | invsqrt | LLM agent loop over the patch layer | 1–5 | 1.02 | 0.25 | 40 | 0.0064 | TODO |
| INV-A09 | invsqrt | Re-synthesise the PHASE endpoint ring | 0.3–1.2 | 0.25 | 0.30 | 12 | 0.0063 | TODO |
| INV-A06 | invsqrt | Domain-aware window re-synthesis scanner | 0.2–1.0 | 0.20 | 0.50 | 16 | 0.0063 | TODO |

## 4. invsqrt

Details (mechanism, risks, sources): `Analysis/research/invsqrt.md`.

### A. Improve our code / tooling (25 ideas)

| ID | Idea | Gain as stated | Score gain % | P | Effort h | Priority | Status |
|---|---|---|---:|---:|---:|---:|---|
| INV-A01 | Submit the public frontier file `563d0583` now, credited, while A02–A04 rebuild it from source | 7.0 | 2.45 | 0.95 | 1 | 2.3261 | EXCLUDED |
| INV-A02 | Port Ali Haydar Özdağ's 64-row coefficient table into our generator | 1.5–2.0 | 0.59 | 0.75 | 16 | 0.0277 | TODO (effort re-estimated at 24–40 h on 18 Nov: no table or loader search tools in the v12 source; steps S2.7–S2.10) |
| INV-A03 | Port Syed's v13/v14 46 → 44 moves into the source | 4.3–4.8 | 1.56 | 0.50 | 30 | 0.0261 | TODO |
| INV-A04 | Public micro-edit pack as generator passes | 0.5–0.9 | 0.23 | 0.70 | 8 | 0.0205 | TODO |
| INV-A05 | Forward normaliser with measurement-conditioned clean-up, like the inverse | 0.2–0.5 | 0.12 | 0.40 | 6 | 0.0078 | TODO |
| INV-A10 | Clean up rotation bookkeeping | 0.01–0.06 | 0.01 | 0.60 | 1 | 0.0070 | TODO |
| INV-A08 | LLM agent loop over the patch layer | 1–5 | 1.02 | 0.25 | 40 | 0.0064 | TODO |
| INV-A09 | Re-synthesise the PHASE endpoint ring | 0.3–1.2 | 0.25 | 0.30 | 12 | 0.0063 | TODO |
| INV-A06 | Domain-aware window re-synthesis scanner | 0.2–1.0 | 0.20 | 0.50 | 16 | 0.0063 | TODO |
| INV-A12 | Margin-maximising joint refit (MILP) | 0.15–0.6 | 0.13 | 0.50 | 10 | 0.0063 | TODO |
| INV-A07 | SAT-minimal rebuild programs with domain don't-cares | 0.1–0.5 | 0.10 | 0.50 | 8 | 0.0063 | TODO |
| INV-A11 | Remove the inverse-product premium | 0.4–1.5 | 0.32 | 0.30 | 16 | 0.0060 | TODO |
| INV-A13 | Trim the consuming root | 0.2–0.8 | 0.17 | 0.35 | 12 | 0.0049 | TODO |
| INV-A14 | Minimum-cost even-shift normaliser | 0.3–1.0 | 0.22 | 0.35 | 16 | 0.0048 | TODO |
| INV-A15 | Measurement-conditioned repair instead of reverse computation | 0.2–0.8 | 0.17 | 0.30 | 12 | 0.0042 | TODO |
| INV-A16 | Replace the last divider row by a rotation on the remainder sign | 0.3–1.0 | 0.22 | 0.15 | 8 | 0.0041 | TODO |
| INV-A17 | Reduce the divider clean-up premium | 0.3–1.8 | 0.35 | 0.25 | 24 | 0.0037 | TODO |
| INV-A18 | Loader-aware table rounding | 0.2–0.7 | 0.15 | 0.35 | 16 | 0.0033 | TODO |
| INV-A19 | Table segmentation sweep | 0.3–1.2 | 0.25 | 0.25 | 20 | 0.0031 | TODO |
| INV-A20 | Our own rescheduler pass | 0–0.3 | 0.05 | 0.60 | 10 | 0.0030 | TODO |
| INV-A21 | Hide more root bits across the multiplier | 0.5–2.3 | 0.47 | 0.15 | 24 | 0.0029 | TODO |
| INV-A23 | Share address logic between the two table loads | 0.1–0.5 | 0.10 | 0.25 | 10 | 0.0025 | TODO |
| INV-A22 | Cubic / basis-change parking census | 0–2.3 | 0.39 | 0.10 | 16 | 0.0024 | TODO |
| INV-A24 | 44 → 43 by removing peak carries | 1.0–2.3 net | 0.56 | 0.12 | 40 | 0.0017 | TODO |
| INV-A25 | Split the coefficient load | −2 to +4 | 0.34 | 0.12 | 40 | 0.0010 | TODO |

### B. Rewrite with a different technique / new pipeline (9 ideas)

| ID | Idea | Gain as stated | Score gain % | P | Effort h | Priority | Status |
|---|---|---|---:|---:|---:|---:|---|
| INV-B01 | Evolutionary LLM search over algorithm families | 0–25 | 4.55 | 0.08 | 80 | 0.0046 | TODO |
| INV-B02 | Fused reciprocal-square-root digit recurrence | 5–20 | 4.55 | 0.05 | 80 | 0.0028 | TODO |
| INV-B05 | Tolerance-approximate carries | 0.5–3 | 0.59 | 0.10 | 30 | 0.0020 | TODO |
| INV-B03 | Compressed body encoding | 2.3–4.6 | 1.18 | 0.10 | 60 | 0.0020 | TODO |
| INV-B04 | Odd-shift normalisation with parity-selected table | 0–3 | 0.51 | 0.15 | 40 | 0.0019 | TODO |
| INV-B08 | Narrow Horner hybrid | −5 to +10 | 0.85 | 0.07 | 60 | 0.0010 | TODO |
| INV-B07 | Division-free Goldschmidt with streamed phase | −5 to +8 | 0.51 | 0.08 | 50 | 0.0008 | TODO |
| INV-B06 | History-free reverse divider (hybrid) | 0–3 | 0.51 | 0.05 | 40 | 0.0006 | TODO |
| INV-B09 | Coarser normalisation (shift by 4^2k) | −1 to +2 | 0.17 | 0.10 | 30 | 0.0006 | TODO |

## 5. dqi

Details (mechanism, risks, sources): `Analysis/research/dqi.md`.

### A. Improve our code / tooling (19 ideas)

| ID | Idea | Gain as stated | Score gain % | P | Effort h | Priority | Status |
|---|---|---|---:|---:|---:|---:|---|
| DQI-A01 | Recover Syed's v13 delta and port it into the v12 generator | 3.0 to 3.6 | 1.12 | 0.55 | 40 | 0.0155 | TODO |
| DQI-A02 | Norm-preserving GF(2^24) closing in `gf24_nested_norm.emit_closing` | 0.21 to 0.26 | 0.08 | 0.85 | 6 | 0.0111 | TODO |
| DQI-A03 | Extractor cutoff j_min 15 -> 18 | 0.05 to 0.08 net | 0.02 | 0.95 | 2 | 0.0103 | DONE (S2.5: cutoff 17 kept, expected −0.048% task cost; measured, not the estimated −0.05 to −0.08%) |
| DQI-A04 | boristown's native patches inside the generator | 0.10 to 0.25 | 0.06 | 0.70 | 12 | 0.0034 | TODO |
| DQI-A08 | Width / Toffoli Pareto pass | 0 to 1.0 | 0.17 | 0.30 | 16 | 0.0031 | TODO |
| DQI-A05 | Generator-level measurement-based uncomputation audit | 0.05 to 0.30 | 0.06 | 0.45 | 10 | 0.0026 | TODO |
| DQI-A07 | Failure-budget audit of rare-branch handlers | 0 to 0.6 | 0.10 | 0.35 | 14 | 0.0025 | TODO |
| DQI-A06 | Product retention / caching with exact lifetimes from the source | 0.10 to 0.40 | 0.08 | 0.60 | 20 | 0.0025 | TODO |
| DQI-A10 | Closing 15-AND product: measured uncomputation and a second helper | 0.013 to 0.18 | 0.03 | 0.50 | 8 | 0.0020 | TODO |
| DQI-A11 | Salt-robust selection and failure attribution | 0 to 0.05 (risk control) | 0.01 | 0.90 | 4 | 0.0019 | TODO |
| DQI-A12 | Atakan Aldemir's four-Toffoli phase component | 0.01 to 0.03 | 0.01 | 0.80 | 3 | 0.0018 | TODO |
| DQI-A13 | Cheaper GF(2^12) inverse (66 T today) | 0.10 to 0.40 | 0.08 | 0.30 | 16 | 0.0016 | TODO |
| DQI-A09 | Liveness profiler and peak shave of the 889 windows | 0 to 0.45 | 0.08 | 0.25 | 14 | 0.0013 | TODO |
| DQI-A14 | Phase-debt merging audit | 0.05 to 0.30 | 0.06 | 0.30 | 14 | 0.0013 | TODO |
| DQI-A15 | Exhaustive re-synthesis of small recurring kernels | 0.10 to 0.50 | 0.10 | 0.20 | 20 | 0.0010 | TODO |
| DQI-A16 | Better middle-product ranks in the Toeplitz solver | 0 to 0.10 | 0.02 | 0.30 | 6 | 0.0008 | TODO |
| DQI-A17 | Cheap guard for the dominant failure class | 0.05 to 0.15 | 0.03 | 0.35 | 16 | 0.0007 | TODO |
| DQI-A19 | Operand-sharing block decompositions | 0 to 0.50 | 0.08 | 0.10 | 20 | 0.0004 | TODO |
| DQI-A18 | Larger adaptive mask gauge | 0.005 to 0.015 | 0.00 | 0.60 | 6 | 0.0003 | TODO |

### B. Rewrite with a different technique / new pipeline (12 ideas)

| ID | Idea | Gain as stated | Score gain % | P | Effort h | Priority | Status |
|---|---|---|---:|---:|---:|---:|---|
| DQI-B01 | Agentic evolutionary loop over the generator's kernels | 0.5 to 3 | 0.59 | 0.30 | 40 | 0.0044 | TODO |
| DQI-B03 | Fuse the projective finish into the synthesis output basis | 0.5 to 2.0 | 0.42 | 0.15 | 30 | 0.0021 | TODO |
| DQI-B02 | Make the reverse key as cheap as the forward key | 1.0 to 3.8 | 0.81 | 0.15 | 60 | 0.0020 | TODO |
| DQI-B04 | Quadratic-node CRT (K24) form for the convert stage | 1 to 6 | 1.19 | 0.10 | 100 | 0.0012 | TODO |
| DQI-B05 | Bezout reuse from the extractor's convergents | 1 to 5 | 1.02 | 0.07 | 60 | 0.0012 | TODO |
| DQI-B07 | Recover branch decisions from the state instead of storing flags | 0.2 to 1.0 | 0.20 | 0.15 | 30 | 0.0010 | TODO |
| DQI-B06 | Exact-synthesis search for one fused extractor slot | 0.5 to 3 | 0.59 | 0.08 | 50 | 0.0009 | TODO |
| DQI-B08 | GF(2^24) bilinear rank below 99 | 0.3 to 1.5 | 0.30 | 0.10 | 40 | 0.0008 | TODO |
| DQI-B09 | Explicit incremental synthesis so that queries become free evaluations | 0 to 10 | 1.72 | 0.03 | 80 | 0.0006 | TODO |
| DQI-B10 | Two slots per chain with a shared pivot operand | 0 to 3 | 0.51 | 0.05 | 40 | 0.0006 | TODO |
| DQI-B11 | 12-qubit carry bank | -1 to 2.7 | 0.28 | 0.05 | 60 | 0.0002 | TODO |
| DQI-B12 | Cubic-node CRT (GF(2^36) heads) | -2 to 4 | 0.34 | 0.08 | 150 | 0.0002 | TODO |

## 6. Chemistry (dfthcsos)

Details (mechanism, risks, sources): `Analysis/research/chemistry.md`.

### A. Improve our code / tooling (13 ideas)

| ID | Idea | Gain as stated | Score gain % | P | Effort h | Priority | Status |
|---|---|---|---:|---:|---:|---:|---|
| CHEM-A01 | Adopt the public frontier parts: Lera's coded table + maoxin123's 4,022 loaders + Atakan Aldemir's weights `703843d3` | 0.67 (0.6-0.7) | 0.22 | 0.95 | 2 | 0.1066 | EXCLUDED |
| CHEM-A02 | Daily "parts market" watch: pull new public chemistry loaders and weights, re-verify, adopt the cheapest compatible pair | 0.3 (0.1-0.7) | 0.10 | 0.70 | 4 | 0.0175 | TODO |
| CHEM-A03 | Rebuild Lera's table loaders with Syed Asad Ali's compiler (our source) | 0.6 (0.5-0.7) | 0.20 | 0.50 | 16 | 0.0063 | TODO |
| CHEM-A04 | Code the one-body rows harder | 0.3 (0.1-0.6) | 0.10 | 0.35 | 8 | 0.0044 | TODO |
| CHEM-A05 | Port the codes (bit-0 affine, XOR-of-five at bits 1-2) onto our own table | 0.75 (0.5-1.0) | 0.25 | 0.40 | 24 | 0.0042 | TODO |
| CHEM-A06 | Use the spare qubits in loaders k1..k14 (126 of 128) | 0.15 (0.05-0.3) | 0.05 | 0.50 | 6 | 0.0042 | TODO |
| CHEM-A07 | Per-slot best-of between maoxin123's and Lera's loaders on the same words | 0.05 (0-0.1) | 0.02 | 0.50 | 2 | 0.0042 | TODO |
| CHEM-A08 | Search the low/high address split and low-block partition per loader | 0.3 (0.1-0.5) | 0.10 | 0.40 | 12 | 0.0033 | TODO |
| CHEM-A09 | Bring loaders k0..k3 down to the k4..k14 level | 0.65 (0.3-1.2) | 0.22 | 0.30 | 20 | 0.0033 | TODO |
| CHEM-A10 | Park more helper monomials in free word bits | 0.25 (0.1-0.4) | 0.08 | 0.30 | 10 | 0.0025 | TODO |
| CHEM-A11 | Peephole / kernel-relation search on compiled loaders | 0.2 (0.05-0.4) | 0.07 | 0.30 | 10 | 0.0020 | TODO |
| CHEM-A12 | Exact re-planning toward the counting floor | 1.5 (0.5-5) | 0.51 | 0.15 | 40 | 0.0019 | TODO |
| CHEM-A13 | Global (multi-start) weight refit at fixed words | 0.03 (0.01-0.1) | 0.01 | 0.30 | 6 | 0.0005 | TODO |

### B. Rewrite with a different technique / new pipeline (18 ideas)

| ID | Idea | Gain as stated | Score gain % | P | Effort h | Priority | Status |
|---|---|---|---:|---:|---:|---:|---|
| CHEM-B01 | Mixed-precision coded words, trained with a loader-cost surrogate | 2.5 (1-4) | 0.85 | 0.25 | 50 | 0.0042 | TODO |
| CHEM-B02 | Row-sensitivity-weighted coding | 0.8 (0.3-1.5) | 0.27 | 0.25 | 20 | 0.0034 | TODO |
| CHEM-B03 | Address-permutation gauge | 1.0 (0.5-2) | 0.34 | 0.20 | 30 | 0.0022 | TODO |
| CHEM-B05 | Exchange-rate-aware error budget | 1.0 (0.5-2) | 0.34 | 0.25 | 40 | 0.0021 | TODO |
| CHEM-B04 | New table at C = 43 | 0.47 (0.4-0.53) | 0.16 | 0.40 | 30 | 0.0021 | TODO |
| CHEM-B06 | Better initializers for lower-C tables | 0.85 (0.6-1.1) | 0.28 | 0.20 | 30 | 0.0019 | TODO |
| CHEM-B07 | Fewer basis vectors per group: B = 139 (R 8, C 44) | 1.5 (1.0-1.9) | 0.51 | 0.15 | 40 | 0.0019 | TODO |
| CHEM-B08 | New table at C = 42 | 0.85 (0.6-1.1) | 0.28 | 0.25 | 40 | 0.0018 | TODO |
| CHEM-B09 | Co-design the Givens tree, edge order and gauge with the codes | 1.0 (0.5-2) | 0.34 | 0.20 | 40 | 0.0017 | TODO |
| CHEM-B10 | R = 7, B = 143, C = 48-50 | 7 (6-7.5) | 2.45 | 0.04 | 60 | 0.0016 | TODO |
| CHEM-B11 | brot 7 with codes and QAT | 9 (8-10) | 3.19 | 0.03 | 60 | 0.0016 | TODO |
| CHEM-B12 | New table at C = 40 | 1.8 (1.5-2.1) | 0.61 | 0.10 | 50 | 0.0012 | TODO |
| CHEM-B13 | Loader architecture from XAG multiplicative-complexity synthesis | 2 (1-5) | 0.68 | 0.10 | 60 | 0.0011 | TODO |
| CHEM-B14 | Use the word registers as dirty workspace | 0.5 (0-2) | 0.17 | 0.10 | 30 | 0.0006 | TODO |
| CHEM-B15 | J = 4 with larger caches | 0.4 (-1 to +1) | 0.13 | 0.10 | 30 | 0.0004 | TODO |
| CHEM-B16 | Direction-symmetric word features | 1 (0.5-3) | 0.34 | 0.05 | 40 | 0.0004 | TODO |
| CHEM-B17 | J = 6 or 7 | 0.2 (-0.5 to +0.5) | 0.07 | 0.15 | 30 | 0.0003 | TODO |
| CHEM-B18 | Multi-basin lambda refit on a new table | 0.1 (0.05-0.25) | 0.03 | 0.20 | 30 | 0.0002 | TODO |

## 7. Cross-cutting methods and tooling

Details (mechanism, risks, sources): `Analysis/research/general.md`.

### A. Improve our code / tooling (15 ideas)

| ID | Idea | Gain as stated | Score gain % | P | Effort h | Priority | Status |
|---|---|---|---:|---:|---:|---:|---|
| GEN-A01 | kickmix-engine as our fast evaluator | 7–9× faster verification (author-reported), dqi memory 9 GB → < 1 GB so 4+ workers fit in WSL; ≈ 0.5% eq. | 0.50 | 0.85 | 3 | 0.1417 | DONE (S2.1: identical results; dqi 6–7×, invsqrt 12×, chemistry 1×) |
| GEN-A05 | Nulls ledger seeded from the community | Avoids repeating closed work; ≈ 0.3% eq. | 0.30 | 0.95 | 3 | 0.0950 | DONE (S2.4: `experiments/NULLS.md`, `experiments/LEDGER.md`) |
| GEN-A02 | Static IR + exact static cost | Evaluator tier T0 in seconds instead of minutes; ≈ 0.4% eq. | 0.40 | 0.90 | 6 | 0.0600 | DONE (S2.2: `tools/kmb.py`, < 1 s per circuit) |
| GEN-A06 | Daily frontier-ingestion agent | Keeps us on the moving plateau without manual reading; ≈ 0.6% eq. | 0.60 | 0.80 | 8 | 0.0600 | TODO |
| GEN-A04 | Salt-robustness gate for every HMR-dependent edit | Avoids −0.15% per dqi failed shot and the 0-score case; ≈ 0.3% eq. (risk-weighted) | 0.30 | 0.90 | 5 | 0.0540 | IN PROGRESS (S2.3: multi-salt gate in `tools/gate.py`; targeted dqi inputs in M2.3) |
| GEN-A07 | I/O and caching hygiene | 1.5–3× on load-heavy dqi runs, no repeated runs; ≈ 0.1% eq. | 0.10 | 0.90 | 2 | 0.0450 | DONE (S2.3: work on ~/work ext4, result cache keyed by SHA and salt) |
| GEN-A09 | Kaggle batch-verification orchestrator | Frees the laptop; 48-salt dqi sweeps and full invsqrt runs in parallel; ≈ 0.15% eq. | 0.15 | 0.80 | 6 | 0.0200 | TODO |
| GEN-A03 | Task-agnostic finishing pipeline | 0.01–0.1% per task per round; ≈ 0.25% | 0.25 | 0.60 | 12 | 0.0125 | TODO |
| GEN-A10 | Local proof per rewrite | Catches rare-branch phase bugs early, fewer wasted full runs; ≈ 0.2% eq. | 0.20 | 0.60 | 10 | 0.0120 | TODO |
| GEN-A12 | Design-space optimizer over generator parameters | ≈ 0.3% from better trade-offs | 0.30 | 0.40 | 10 | 0.0120 | TODO |
| GEN-A08 | Peak map + functional-dependency mining | ≈ 0.5% (invsqrt qubit 0.77%, dqi 0.075%) | 0.50 | 0.20 | 10 | 0.0100 | TODO |
| GEN-A13 | Scaled-down proxy instances | Faster triage of new families; ≈ 0.2% eq. | 0.20 | 0.40 | 12 | 0.0067 | TODO |
| GEN-A11 | Gadget library + matcher | 0.05–0.2% | 0.12 | 0.45 | 10 | 0.0056 | TODO |
| GEN-A15 | Branch-cost balancing | ≈ 0.08% | 0.08 | 0.30 | 6 | 0.0040 | TODO |
| GEN-A14 | CUDA bit-sliced simulator (RTX 4060) | ≤ 2–3× over 12 cores; ≈ 0.05% eq. | 0.05 | 0.30 | 12 | 0.0013 | TODO |

### B. Rewrite with a different technique / new pipeline (17 ideas)

| ID | Idea | Gain as stated | Score gain % | P | Effort h | Priority | Status |
|---|---|---|---:|---:|---:|---:|---|
| GEN-B11 | Divergent-family brainstorm + red team, weekly | ≈ 1% (one new family that beats the incumbent on one stage) | 1.00 | 0.25 | 8 | 0.0312 | TODO |
| GEN-B01 | Evolutionary agent loop over our generators | 2–5% over the remaining weeks (Gu et al. found 4–39× over **literature** baselines; ours are crowd-optimized, so far less) ≈ 3% | 3.00 | 0.35 | 40 | 0.0262 | TODO |
| GEN-B02 | Principle sweep, stage by stage | ≈ 1.5% (a few stage-level wins) | 1.50 | 0.30 | 20 | 0.0225 | TODO |
| GEN-B04 | Qubit-cap ladder harness | invsqrt −1 qubit = 0.77% | 0.77 | 0.20 | 10 | 0.0154 | TODO |
| GEN-B12 | One end-to-end error budget, fitted to the exact quantized computation | ≈ 0.4% (one table bit or a few rotations) | 0.40 | 0.25 | 12 | 0.0083 | TODO |
| GEN-B03 | Width-optimal rescheduling by spooky pebbling DP/SAT | ≈ 0.77% (one invsqrt qubit; open best case −0.37% cost) | 0.77 | 0.20 | 20 | 0.0077 | TODO |
| GEN-B10 | Conditional-routing detector | ≈ 0.5% where it hits | 0.50 | 0.15 | 12 | 0.0062 | TODO |
| GEN-B08 | Conditionally clean ancillae everywhere | one invsqrt qubit = 0.77% | 0.77 | 0.12 | 16 | 0.0058 | TODO |
| GEN-B13 | Shared library of optimal primitives | ≈ 0.3% | 0.30 | 0.30 | 16 | 0.0056 | TODO |
| GEN-B06 | XAG multiplicative-complexity resynthesis with reachable-state don't-cares | ≈ 0.3% | 0.30 | 0.30 | 20 | 0.0045 | TODO |
| GEN-B07 | SAT superoptimizer for small windows | ≈ 0.2% | 0.20 | 0.30 | 16 | 0.0037 | TODO |
| GEN-B09 | Evolve the search heuristics, not the circuits | ≈ 0.3% | 0.30 | 0.20 | 16 | 0.0037 | TODO |
| GEN-B05 | Minimum-CCZ phase polynomials | ≈ 0.15% | 0.15 | 0.25 | 12 | 0.0031 | TODO |
| GEN-B15 | ZX-calculus (PyZX) Toffoli reduction | ≈ 0.1% | 0.10 | 0.03 | 6 | 0.0005 | TODO |
| GEN-B14 | RL / MCTS over rewrite sequences | ≈ 0.3% | 0.30 | 0.05 | 40 | 0.0004 | TODO |
| GEN-B17 | Lean-verified rewrites (DeepMind style) | ≈ 0.05% | 0.05 | 0.10 | 40 | 0.0001 | TODO |
| GEN-B16 | Train AlphaTensor-Quantum from scratch | ≈ 0.2% | 0.20 | 0.03 | 60 | 0.0001 | TODO |

## 8. Combination rules

The score is a product of the three task ratios, so the tasks are independent: the best version of each task
is always combined, with no further test. Combinations matter **within a task** (user decision,
2026-11-17): ideas that pass are stacked one at a time, re-verifying after each addition, and two ideas are
tested jointly only when they touch the same part of the circuit or share a budget. The research reports
(section 6 of each) give the details; the binding constraints are:

**invsqrt**
- The frontier path is one ordered chain: INV-A02 (64-row table) → INV-A03 (46 → 44 moves) → INV-A04
  (micro-edits). Re-run the scanners after each structural step instead of replaying public line patches.
- Width ideas (A21, A22, A24, A25, with A20's rescheduler first) only pay together: 43 qubits needs a cover of
  every 44-wide window within one budget of ~39 Toffolis.
- Rotation margin is shared by A02, A12, A16, A18, A19 and B05: run A12 (margin refit) last, and A16/B05 only
  after A12.
- Table families are mutually exclusive (A02 uniform, A19 non-uniform, B04 parity-selected, B09 wider range):
  pick one, then apply A18 and A12 on it.
- Divider ideas (A09, A16, A17, B06) edit the same rows; A03 first. Normaliser: A14 supersedes A05.

**dqi**
- DQI-A01 (port Syed's v13 delta) comes first; the public patches (A02, A04, A06) are re-derived on top of it.
- Every change to rare-branch handling (A03 cutoff, A07, A17) is judged on the multi-salt expected cost and
  on Atakan Aldemir's adversarial low-weight inputs, never on one draw.
- B ideas that restructure the key or convert stage (B02, B04, B05, B09) reset the site-specific patches.

**Chemistry**
- One table line: either Lera's (8,143,44) coded table rebuilt from our source (CHEM-A03) or our own table with
  ported codes (A05/A09). Not both.
- Everything that changes words or weights spends the same 20 mHa budget (A04, A05, A09, B01–B03, and all shape
  moves B04, B07, B08, B10, B12): build the budget optimizer (B05) before stacking more than two of them.
- Shape changes reset the loaders, words and weights (full rebuild, Kaggle P100). Loader-only ideas (A06, A07,
  A08, A10, A11, A12, B13) stack with any table.
- λ_eff ideas (A13, B18) are capped at 0.25% in total.

**Cross-cutting**
- GEN-A01, A02, A05 and A07 are prerequisites of the search ideas; GEN-B01 is the vehicle in which GEN-B02,
  B04, B10 and B12 run as mutation operators, seeded by GEN-B11.

## 9. Known nulls

Do not start an idea without checking the nulls of its task: `Analysis/research/invsqrt.md` §4 (28 nulls),
`dqi.md` §4, `chemistry.md` §4, `general.md` §3 (23 nulls). A null is closed only in its design family; a new
family can reopen it, but must say why it escapes the recorded bound.
