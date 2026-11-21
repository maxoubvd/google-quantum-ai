# Action plan: Google Quantum AI — Verifiable Quantum Circuit Compilation

**Current plan, version of 2026-11-21 11:51 UTC** (re-plan of 21 Nov on the ThemBones source base). This file holds
only the plan in force; earlier versions are in `Analysis/plan_history/` (index in its README), the work done in
`changelog/`, the ideas in `Analysis/IDEAS.md`, the measured attempts in `experiments/LEDGER.md` and the known
negative results in `experiments/NULLS.md`.

**The calendar is indicative.** Phases end on their goal ("done when"), not on their date: if a goal is reached
early, the next phase starts early; if it slips, we cut scope rather than skip verification. Only the Kaggle
deadlines are fixed: team merger and entry **2027-01-13**, final submission **2027-01-20 23:59 UTC**, and the rule of
submitting at least 24 h before any deadline.

**Targets are relative, never fixed scores.** The board rises fast (28.23 → 33.10 on 20–21 Nov), so every goal is
stated against the board of the day: the public frontier, a medal line, the gap to the leader. Section 1 is a dated
snapshot, refreshed at the start of each phase.

| Phase | Dates (indicative) | Target, measured on the board at that time | Status |
|---|---|---|---|
| 1 | 17 → 20 Nov | A first reproducible submission | DONE (27.00209) |
| 2 | 21 → 30 Nov | At the public frontier (credited public files), then the ThemBones source base rebuilt in our repo and pushed beyond the public files | IN PROGRESS |
| — | 1 Dec → 3 Jan | Break | — |
| 3 | 4 → 16 Jan | Back at the frontier after the December rise, then beyond it from our own source | — |
| 4 | 17 → 20 Jan | Best private-board expectation: robustness on the private salt first, then score | — |

## 1. Where the competition stands (snapshot 2026-11-21)

On 20 Nov team ThemBones (Danish Ali) published the circuits and generators of its 31.68 level (Apache-2.0). Within
a day, other authors improved them and the public frontier rose from 28.234 to **33.10186**. The board top is a tie
at about 33.10–33.13; the leader is at 33.348. Analysis: `Analysis/research/public_frontier_2026-11-21.md`.

| Task | Best public file (21 Nov) | Qubits | Toffolis | Cost | Generator source |
|---|---|---:|---:|---:|---|
| invsqrt | `d2148a7c` (Rithvik Redrouthu, on Syed Asad Ali's line via ThemBones and Almond Labs) | 37 | 1,406.3 | 60,955 | none public for 37 Q; ThemBones' 39 Q generator is public |
| dqi | `4b572141` (Tony Li, on ThemBones' designs) | 900 | 316,131 | 3.050e11 | ThemBones kit (their file A never rebuilt); tonylica from recovered KMX assembly |
| dfthcsos | table `fc4b4c01` (Ertugrul Ozer) + 15 loaders chained through dead words (goodpjw2008), (R,B,C) = (9,127,46), brot 7 | 405 | 13,290 | 2.525e8 | ThemBones chemistry tools and table |

Our submissions: 27.00209 (18 Nov, v12 base), 28.23444 (19 Nov, public set), and the 33.10 public set packaged on
21 Nov (S2.25). Exchange rates at the new frontier:

| Lever | Task cost change | Score change |
|---|---|---|
| invsqrt −1 qubit (37 → 36) | −2.7% | +0.90% |
| invsqrt −14 Toffoli-equivalents (1% of ~1,647) | −1% | +0.33% |
| dqi −1 qubit (Q² term, 900 Q) | −0.22% | +0.074% |
| dqi one extra failed shot on 102,400 | +0.3 to +0.5% | −0.1 to −0.15% |
| chemistry −1 qubit of Q_total = 405 | −0.25% | +0.082% |
| chemistry −1% Toffolis or λ_eff | −1% | +0.33% |

A 10% cost cut on one task is +3.57% on the score (geometric mean of three ratios).

## 2. Constraints

- 5 submissions per day; scoring takes about 45 min, and queues of 8–9 hours have been reported. Never rely on a
  same-day submission before a deadline.
- 2 final submissions are selected for the private board (the other half of the test set, another secret salt).
  dqi failures move with the salt; invsqrt and chemistry do not when they are proven on every input.
- One failed task scores the whole submission 0.
- Prize eligibility needs a generator that reproduces the submitted circuits from Apache-2.0 source we hold
  (host, topic 752010; rules B.5, B.7). Whether recovered assembly (KMX) counts as source is to be asked (M3.6).
- Laptop: 12 threads, 15.7 GB RAM (WSL capped at 8 GB), RTX 4060 Laptop 8 GB. Kaggle: CPU sessions (30 GB RAM, 9 h),
  P100 / T4 ×2 (30 h of GPU per week).

## 3. Licensing, provenance and where to run what

- `LICENSE` (Apache-2.0, our code), `NOTICE` (our copyright, then every upstream work with author and licence),
  `third_party/` (upstream material with its LICENSE/NOTICE verbatim; `PROVENANCE.md` records author, URL, licence
  and SHA-256 of every reused file). Syed Asad Ali's NOTICE requires naming him first in every lineage line.
- Not ours to license: `Scientific Papers/`, `Instructions/`, `Analysis/Gemini Analysis.pdf`.
- The GitHub repository stays **private** until the end (rules C.5.d).

| Workload | Place | Why |
|---|---|---|
| Verifiers, generators, search loops | Laptop, WSL (`~/venvs/gqai`, Python 3.14; `~/venvs/gqai313` for the public generators; `~/venvs/sat` for SAT) | 12 local threads, no upload round-trip |
| Full dqi verification at 4 workers, long sweeps | Kaggle CPU sessions | WSL is capped at 8 GB |
| Chemistry fits in float64 | Kaggle P100 | float64 at 1/2 rate (1/32–1/64 on the RTX 4060 and the T4) |
| LLM proposer for a search loop (if M3.4a passes) | API (Haiku, Sonnet escalation, spend cap); a local model only as a fallback | `Analysis/research/agent_loop_local_llm_2026-11-21.md` |

## 4. How we work

- **Milestones** (`M<phase>.<n>`) implement and test ideas from `Analysis/IDEAS.md` (or deliver tooling or a
  submission), with a time box, an expected gain and a **go / no-go** criterion. **Steps** (`S<phase>.<n>`) are half
  a day to a day of work, each with a goal, a "done when" test and a status (`TODO`, `IN PROGRESS`, `DONE`,
  `DROPPED`, `PAUSED`, `WAITING`). Steps of the next milestones are written out; later ones get theirs at the review
  that opens them.
- **Start of a step:** state its goal and "done when" in the chat. **End of a step:** a `changelog/` entry, the
  status updated here and in `IDEAS.md`, a `LEDGER.md` row, a short review with the user and a proposed commit
  message. The user can reorder, add or drop steps at every review.
- **Order** by the priority of `IDEAS.md` (P × score gain ÷ effort), adjusted by what we measure.
- **Combinations:** the tasks are independent. Within a task, ideas that pass are stacked one at a time with a
  re-verification after each; two ideas are tested jointly only when they touch the same part of the circuit or
  share a budget (rotation margin, 20 mHa error budget, a qubit cap).
- **Speed is a deliverable, never at the cost of exactness** (rule in `CLAUDE.md`): iterate on static costs, reduced
  runs, the generators' exact proofs and kickmix-engine; the official verifiers run at gates.
- **Gates before a submission:** byte-identical rebuild from source when the file is ours, exact proof where it exists
  (invsqrt, chemistry loaders), `tools/gate.py` (several salts + official verifiers), `tools/package.py`; then the
  user submits.
- **Nulls state their domain.** A negative result is recorded with the conditions under which it holds, and a
  ceiling is checked cheaply (scorer grid, exact small case) before it is quoted. Four 19 Nov nulls failed this
  test (`experiments/NULLS.md`, corrections of 21 Nov).
- **Exploring a new family:** paper (datapath, qubit peak, Toffoli count, checked against the nulls) → scaled
  prototype verified exhaustively → full scale. Time-boxed, go/no-go at each gate.
- **Parallel work:** independent milestones can run in subagents when the user agrees; the main session reviews.
- **Phase review:** retrospective, fresh board snapshot, next milestones re-ranked from `IDEAS.md`.

## 5. Phase 1 — done (17 → 20 Nov)

Licensing in place; Syed Asad Ali's v12 generators rebuilt byte for byte and verified; `build/` and
`tools/package.py`; first submission **27.00209** (18 Nov); research and idea register (`Analysis/research/`,
`Analysis/IDEAS.md`, 128 ideas at the time). Details: `changelog/` 17–18 Nov and `plan_history/`.

## 6. Phase 2 — public frontier, then the ThemBones source base (21 → 30 Nov)

**Done before this plan (18–20 Nov, details in `changelog/`):** evaluation toolkit (kickmix-engine cross-checked,
`tools/kmb.py`, `tools/gate.py`, ledger and nulls index); dqi cutoff 17 from source; our own K = 6 invsqrt table;
second submission **28.23444** (credited public set, 19 Nov); paper studies of new families for the three tasks
(19 Nov; their ceilings were partly wrong, see §4); parity divider tail E1 GO (INV-B10); DQI-B13 SAT screen,
parked (`tools/mc_sat.py`). The v12-based steps still open then (invsqrt K = 6 integration, chemistry port of Lera's
table, dqi guard on v12) are dropped with the v12 base.

**Decisions (user, 21 Nov):** third submission at the 33.10 public level; new base for all three tasks = the
ThemBones source line (Apache-2.0, byte-for-byte rebuilds), with the public edits rebuilt as generator passes and
then pushed beyond (invsqrt 37 Q from source, dqi 890 Q); our v12 generators are retired as bases, their tools carry
over. The LLM search loop was priced on paper and gets a one-day go/no-go in January (M3.4a); fine-tuning is
excluded before the deadline.

| Milestone | Content | Expected gain | Time box | Go / no-go |
|---|---|---|---|---|
| M2.E Public level 33.10 | S2.25 | 28.23 → about 33.10 | 0.5 day | All gates PASS; the user submits |
| M2.F ThemBones base from source | S2.26, S2.27, S2.29, S2.31 | prize-eligible source at the public level | 4–5 days | Each task rebuilt from our repo at or below its public cost |
| M2.G Beyond the public set | S2.28, S2.30 | +0.5 to +2% | rest of Phase 2 | A task below its public file under the full gates |
| M2.H Freeze before the break | S2.32 | — | 0.5 day | A fresh session rebuilds and resubmits alone |

| Step | Milestone | Goal | Done when | Status |
|---|---|---|---|---|
| S2.25 | M2.E | Adopt the 33.10 public set (goodpjw2008's output, 18 members), pins and credits (`third_party/public-frontier-2026-11-21/`), gate (8 salts + official verifiers), zip `artifacts/submissions/public-33102.zip` | all gates PASS | **DONE, READY** (21 Nov): all PASS; dqi fails 1–2 shots on 3 of 8 salts (projected 33.078 at the salt mean); the user submits |
| S2.26 | M2.F | Reproduce the ThemBones line in `~/venvs/gqai313` (add numba): invsqrt `c639a045`, dqi file A `2c043a9a`, chemistry loaders and npz `e05fb18f`; record run time and memory | byte-for-byte SHA match per task (dqi: first independent rebuild) | TODO (next) |
| S2.27 | M2.F | invsqrt from source to the public 39 Q level: LEND29, T-fusion, then tonylica's divider/carry ideas as passes (INV-A26) | ≤ 63,016 (then ≤ 61,696), exact proof on all inputs | TODO |
| S2.28 | M2.G | invsqrt 39 → 38 → 37 Q from source (INV-A27): liveness diff of `d2148a7c` vs `c639a045`, cap-38/37 `game` DP; INV-B10 re-screen | 38 Q at C ≤ 1,604 or 37 Q at C ≤ 1,647 | TODO |
| S2.29 | M2.F | dqi from source: Almond Labs' five edits, then tonylica's D16 formulas (DQI-A20) | ≤ 3.3248e11, then ≤ 3.05e11 on 16 salts, no more failures than the parent | TODO |
| S2.30 | M2.G | dqi width 900 → 890 (DQI-A08/A09) | ≥ 4 qubits by 30 Nov (−0.22% each) | TODO |
| S2.31 | M2.F | chemistry from source: dead-word criterion, our chain DP, table exchange; CHEM-A10 on the new loaders | 13,290 at λ ≤ 46.9182 from our repo | TODO |
| S2.32 | M2.H | Freeze: `STATUS.md` (score, files, open ideas, dead ends), build and package commands per task, Phase 2 review | a fresh session rebuilds and resubmits | TODO |

## 7. Break (around 1 Dec → 3 Jan)

No work required. Optional: forum questions that would unblock January (host clarifications, e.g. KMX as source),
keep the Kaggle token valid, and, if the user wants it, a scheduled frontier watch (GEN-A06) so that January starts
with a ready list.

## 8. Phase 3 — catch-up, then beyond the public set from source (around 4 → 16 Jan)

Re-ranked at M3.1 from `IDEAS.md`, `STATUS.md` and the January board.

| Milestone | Content | Expected gain | Time box | Go / no-go |
|---|---|---|---|---|
| M3.1 Snapshot and catch-up | Read `changelog/` and `STATUS.md`; board and public notebooks (GEN-A06); re-run S2.25 on whatever moved | back to the frontier of the day | 1–2 days | Score at the frontier of the day |
| M3.2 invsqrt ≤ 36 Q | Search over per-block programs of the ThemBones generator under a wire cap, judged by the 5–8 s exact all-input check (INV-B01, GEN-B01); INV-B10 if it survives S2.28 | −2% to −5% on invsqrt | 3–4 days | Exact proof PASS below the board file |
| M3.3 dqi kernels and width | Multi-output shared-operand bilinear formulas (flip graphs) at every field site (DQI-A20), 890 Q (DQI-A08/A09), failure guard (DQI-A07/A17); DQI-B13 only with a structured resynthesis | −3% to −8% on dqi | 3–4 days | 16+ salts below the board file, no worse failure rate |
| M3.4a LLM loop go/no-go | GEN-B18 (Haiku/Sonnet API proposer, spend cap ≤ $25), GEN-B23 (exemplar bank), GEN-B20 (loader golf on a fixed truth table with an exact ms check, on a ThemBones chemistry or invsqrt block), raced against a non-LLM search of equal time | enabler | 1 day | ≥ 2 ANDs better than the non-LLM baseline within 600 proposals |
| M3.4b LLM loop production | GEN-B01 with GEN-B02 operators on the blocks of M3.2/M3.3/M3.5, ≤ $100, only on a go; API diff proposals, Claude Code subagents only for ideation and gates | +0.3 to +1% | alongside M3.2–M3.5 | Frozen after 3 rounds without gain |
| M3.5 Chemistry fit and loaders | CHEM-B23 scorer grid over (R, B, C, brot, J) first (1 min); rounding-aware, loader-aware fit at brot 7 / C 45 and a brot 6 attempt (CHEM-B24, GPU); CHEM-B19 non-bilinear loaders with dead-word don't-cares (CPU) | −1% to −5% on chemistry | 2–3 days | e1, e2 < 20 mHa with λ ≤ 46.95, cost below the board file |
| M3.6 Team decision (**before 13 Jan, hard deadline**) | Includes asking the hosts early in January whether recovered assembly (KMX) counts as source for the prize | — | — | Decision recorded |

## 9. Phase 4 — final submission (17 → 20 Jan)

Only the method is fixed here; the score aimed for is whatever the January board makes realistic.

| Milestone | Date | Goal | Done when |
|---|---|---|---|
| M4.1 Code freeze | Sun 17 | No more circuit changes | Freeze recorded |
| M4.2 Full verifications | Sun 17 | invsqrt exhaustive proof and full verifier, dqi on 48 salts plus adversarial inputs, dfthcsos | All PASS, failure spread known |
| M4.3 Choose the two finals | Mon 18 | (A) safest, fewest dqi failures across salts; (B) best expected score | Both already scored on Kaggle |
| M4.4 Select on Kaggle | Tue 19 | Selection, with a buffer for a failed scoring or a queue delay | Selection confirmed |
| M4.5 Winner package | 17–20 | Generators, `LICENSE`, `NOTICE`, `PROVENANCE.md`, a methodology write-up per task, one command that rebuilds `submission.zip` with the same SHA-256 (rules B.5 and B.7) | Rebuild from a clean clone matches |
| — | Wed 20 | No new submission unless scored and verified before. Deadline 23:59 UTC | — |

## 10. Risks

| Risk | Mitigation |
|---|---|
| One task fails on the private salt → the whole submission scores 0 | Exhaustive invsqrt proof; dqi and chemistry checked on many salts; final A is the safest portfolio |
| dqi failed shots on the private salt (−0.1 to −0.15% each; the public dqi line fails on 3 of 8 salts) | Rank on 16–48 salts; failure guard (DQI-A07/A17); prefer robustness over the last 0.1% |
| Part not rebuildable from source → prize ineligible | ThemBones source base (M2.F); provenance table; replace any file without a generator; ask the hosts about KMX |
| The public edits cannot all be rebuilt as generator passes | Time boxes; keep the passes that work; the credited public file stays on the board meanwhile |
| A paper null or ceiling is quoted outside its domain (four such errors on 19 Nov) | Every null states its domain; run the cheap check before calling a ceiling |
| Research estimates are optimistic | Every gain is measured before it counts; `IDEAS.md` records the measured value |
| LLM loop: API spend runaway, evaluator reward hacking, local LLM starving the evaluator | Provider-side spend cap; read-only SHA-pinned verifiers, final check on the untouched official verifier; LLM on Windows, evaluator in WSL |
| License breach on reused code | `third_party/` keeps upstream LICENSE/NOTICE; changed files marked; Syed Asad Ali named first; repo private |
| Kaggle queue delay at a deadline | Submit at least 24 h before every deadline |
| Laptop memory (8 GB WSL; the gate peaked at 7.5 GB on 21 Nov) | kickmix-engine for iteration; dqi with `--max-workers 2`; heavy runs on Kaggle CPU notebooks |
| Usage limit stops a session or its agents | `autoContinueAtUsageLimit`; agents save progressively; `changelog/` and this plan let a new session resume |
| Frontier moves away during December | M3.1 catch-up is planned, not optional; optional frontier watch during the break |
