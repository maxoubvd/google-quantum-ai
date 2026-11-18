# Action plan: Google Quantum AI — Verifiable Quantum Circuit Compilation

Written 2026-11-17, updated 2026-11-17 20:00 UTC with the live leaderboard (Kaggle API, 3,274 teams)
and the top public notebooks, and **restructured 2026-11-18 into milestones** drawn from the idea register
`Analysis/IDEAS.md` (step S1.10). Other sources: `Instructions/` (overview, data, rules, forum snapshot of
2026-11-15), `Analysis/Gemini Analysis.pdf`, the verifiers and `CLAUDE.md`.

**The calendar is indicative.** Phases move with progress: if a phase's goal is reached early, the
next one starts early; if it slips, we cut scope (fewer ideas per track) rather than skip verification.
Only the Kaggle deadlines are fixed: team merger and entry 2027-01-13, final submission 2027-01-20
23:59 UTC, and the rule of submitting at least 24 h before any deadline. Each phase ends on its goal
("done when"), not on its date.

**Targets are relative, never fixed scores.** The board rises every day (+11% between 12 and 17 Nov),
so every goal is stated against the board at the time: the public plateau, a medal line, or the gap
to the leader. Section 1 is a dated snapshot, refreshed at the start of each phase (download the
leaderboard and list the top notebooks with the Kaggle CLI); the numbers in it are examples, not goals.

| Phase | Target, measured on the board at that time |
|---|---|
| 1 | At the public plateau (the score of the best public portfolio) |
| 2 | Back to the plateau from our own source, then strictly above it, aiming for the gold line (top 10 + 0.2%) |
| 3 | Back above the plateau after the December rise, then as close to gold as possible |
| 4 | Best private-board expectation: robustness on the private salt first, then score |

## 1. Where the competition stands (2026-11-17)

**Snapshot 2026-11-18 09:00 UTC (Phase 1 review, 3,305 teams):** leader 31.617 (Z³C¹oherence), 6 teams ≥ 30,
gold line (rank 16) 28.741, silver (rank 165) 28.190, bronze (rank 330) 28.181, median 25.148; 537 teams at
28.18 or more. **Our score 27.00209 (rank ≈ 1,318, tied with 36 teams on the same public files).** The tables
below are the 17 Nov snapshot.

### Leaderboard (public, 3,274 teams)

| | Rank | Score |
|---|---:|---:|
| Leader (aadi_stonks) | 1 | **31.445** |
| 5 teams ≥ 30 | 5 | 30.131 |
| Gold line (top 10 + 0.2%) | 16 | **28.713** |
| Silver line (top 5%) | 163 | **28.1807** |
| Bronze line (top 10%) | 327 | **28.151** |
| Median | 1,637 | 25.148 |

- **Public plateau: ≈ 28.18.** 442 teams are at 28.15 or more; the silver line and rank 100 are both at
  28.1807. A medal on today's board needs the plateau at least; gold needs +1.9% above it.
- The top two jumped to 31.4 on 17 Nov, 11.6% above the plateau: someone has a structural gain that is
  not public (equivalent to a 28% cost cut on one task).
- **Ties are broken by the earliest submission** (rules C.7.b). Joining the plateau late ranks us behind
  everyone already on it.
- The public frontier moved +3.3% between 14 and 17 Nov (27.30 → 28.18), after +7.2% on 12–14 Nov.

### Best public file per task

| | invsqrt | dqi | dfthcsos (chemistry) |
|---|---|---|---|
| Reference cost | 260 Q × 7,988 = 2,076,880 | ≈ 2.44 × 10^13 | 3,355,599,588 |
| Best public file | `563d0583` (maoxin123, 17 Nov; on Syed Asad Ali's generator with Ali Haydar Özdağ's 64-row table) | `575b043f` (GrantPrater; merge of boristown's patches and Balram Srinivasan's GF(2^24) closings on Syed Asad Ali's v13 decoder, plus four product-retention pairs) | maoxin123's 15 loaders at 4,022 Toffolis per pass (CC0) with Atakan Aldemir's weights `703843d3`, on Lera's (8,143,44) table |
| Qubits / cost | **44 Q** / 74,560.91 (1,454 mean Toffolis) | **890 Q** / 3.4641 × 10^11 (367,116 Toffolis, clean draw) | 427 Q / 2.94642 × 10^8 (T 14,706, λ_eff 46.92) |
| Ratio | 27.85 | 70.6 | 11.39 |
| Change since 14 Nov | −5.3% (46 → 44 qubits) | −4.0% (889 → 890 Q, −15k Toffolis) | −0.17% |

Sources: alex chilton's *28.15 LB | Which Circuits Survive Private Salt?*, maoxin123's *QCC 28.17 |
4,022-Toffoli chemistry loaders*, Asterios Panoras's *QCC 28.18084 | Frontier + 3 DQI Products*
(all public notebooks, 17 Nov).

### Exchange rates (what one unit is worth on the final score)

The score is a geometric mean, so a 10% cost cut on any one task is +3.57% on the score. Reaching the
gold line from the plateau (+1.9%) needs a 5.5% cost cut on one task; matching the leader (+11.6%) needs 28%.

| Lever | Task cost change | Score change |
|---|---|---|
| invsqrt −1 qubit (44 → 43) | −2.3% | **+0.77%** |
| invsqrt −17 Toffoli-equivalents (1% of ~1,695) | −1% | +0.33% |
| dqi −1 qubit (Q² term) | −0.22% | +0.075% |
| dqi one extra failed shot on 102,400 | +0.46% | −0.15% |
| chemistry −1 qubit of Q_total = 427 | −0.23% | +0.078% |
| chemistry −1% Toffolis or λ_eff | −1% | +0.33% |

invsqrt still has by far the best exchange rate per qubit and no hard width floor (260 → 63 → 46 → 44
so far). dqi is boxed in by 840 incompressible syndrome qubits and the R(ε) floor of 1.19126.
Chemistry is boxed in by λ_eff ≥ ~46.85, the 20 mHa gates and the analytic cost model.

### What the forum says is left (17 Nov)

- alex chilton: "the last two weeks of public gains are refits and single-gate rewrites worth 0.01 to
  0.1% each. No public notebook found a design change worth 2× on any task." The 31.4 teams show that
  bigger changes exist, but privately.
- **The private salt is the real risk.** A circuit that passes at the default salt can fail at another;
  per alex chilton, **one failed task scores the whole submission 0**. Atakan Aldemir published
  weight-16/17 dqi inputs that a "full pass" decoder fails. The best public dqi files draw a failed shot
  on 1 to 4 of 8 salts tested (−0.15% each).

## 2. Assessment of the Gemini analysis

Its numbers are those of the 14 Nov forum snapshot and are now outdated (see section 1). Points to
correct or weigh:

- **Milestone Prize (20 Nov): out of reach.** The leader is at 31.4. Our 20 Nov goal is a solid,
  eligible base, not the prize.
- **Single-task prizes:** realistic only if one of our tasks beats the private leaders, which we cannot
  see. Treat them as a bonus, not a target.
- **Prize eligibility is stricter than "use public files".** The host (Tanuj Khattar, topic 752010) says
  a generator that reproduces the exact circuits is required for all prizes. Every file we submit must be
  rebuilt byte-for-byte from open-source code we hold, with author, license, NOTICE and SHA-256 recorded.
- Its tooling claims (LLM brands, "6×10^10 ops/s") are anecdotal and do not change the plan.
- Its useful core, which this plan adopts: (1) gains come from changing generators, not editing their
  output; (2) measure-then-rebuild, conditional routing and failure budgets are the three recurring
  levers; (3) evaluate dqi on many salts and rank by expected cost; (4) a fast local loop matters more
  than any single idea.

## 3. Constraints that shape the schedule

- 5 submissions per day; scoring takes ~45 min, and queues of 8–9 hours have been reported (topic
  757551). **Never rely on a same-day submission before a deadline.**
- 2 final submissions selected for the private leaderboard (the other 50% of the test set, a different
  secret salt). dqi failures move with the salt; invsqrt and chemistry do not if the circuits are proven
  on every input and measurement transcript.
- Laptop: 12 cores, WSL capped at 8 GB, RTX 4060 Laptop (8 GB). Full invsqrt verification ≈ 1 h; dqi
  needs `--max-workers 2`. Kaggle CPU notebooks (9 h, 30 GB RAM) are a second verification machine.
- Entry and team-merger deadline: **2027-01-13**. Final submission deadline: **2027-01-20 23:59 UTC**.

## 4. Prerequisites

1. ~~Kaggle API token~~ Done 2026-11-17: `~/.kaggle/access_token` on Windows and WSL; the CLI
   (`kaggle==2.2.5`, plus `six`) is installed in both environments and pinned in `requirements.txt`.
2. ~~Rules accepted~~ Done: the account is entered in the competition (first submission 18 Nov).
3. `Notebooks/` holds our own notebooks (versioned), usually run on Kaggle.
4. **Licensing** (section 4.1), before any third-party code enters the repo.

### 4.1 Licensing

The rules (B.5, C.6.b) license our submission and the code that generates it under **Apache-2.0**, and
the upstream generators we build on are Apache-2.0 (Syed Asad Ali, Lera, prvsiyan, …) or CC0 (maoxin123).
Apache-2.0 §4 requires keeping their license and NOTICE files and marking files we change. So:

- `LICENSE` at the repo root: the Apache-2.0 text, stating that it covers our own code (generators,
  tools, notebooks), not the third-party material listed below.
- `NOTICE` at the root: our copyright line, then every upstream work we include, with its author and
  license.
- `third_party/<author>-<work>/`: upstream code kept with its original `LICENSE`/`NOTICE` verbatim;
  changed files get a header line saying so. `third_party/PROVENANCE.md` records author, URL, version,
  license and SHA-256 of every reused file.
- Not ours to license, and not to be redistributed if the repo becomes public: `Scientific Papers/`
  (publisher copyright), `Instructions/` (Kaggle pages and forum posts), `Analysis/Gemini Analysis.pdf`.
  The competition data keeps its own Apache-2.0 license.
- Keep the GitHub repository **private** until the end of the competition (it is, as of 17 Nov):
  sharing code outside the team is only allowed if made available to all participants on the forums
  (rules C.5.d).

### 4.2 Where to run what

| Workload | Best place | Why |
|---|---|---|
| Verifiers, generators, search loops (kickmix) | **Laptop, WSL** | kickmix is CPU-only; 12 local cores against 4 in a Kaggle CPU session; no upload round-trip |
| Full dqi verification at 4 workers (~9 GB) | Kaggle CPU (30 GB RAM) or WSL with a higher memory cap | WSL is capped at 8 GB |
| Long jobs in parallel (e.g. 48-salt dqi sweeps) | Kaggle CPU sessions, several at once (9 h each) | Frees the laptop |
| Chemistry refits in float64 (JAX/PyTorch) | **Kaggle P100** | The RTX 4060 Laptop and the T4 run float64 at 1/32–1/64 of float32; the P100 runs it at 1/2 |
| TPU | Not useful | No float64, and kickmix does not run on it |

## 5. How we work: milestones and steps

Phases 2 to 4 are organised in **milestones**; each milestone is worked as a sequence of small **steps**,
like the stories of an agile sprint.

- **A milestone** (`M<phase>.<n>`) implements and tests one idea or a group of ideas from `Analysis/IDEAS.md`
  (or delivers a piece of tooling or a submission gate). It states its ideas, its time box, its expected gain
  and a **go / no-go** criterion: the measured result that makes us keep the change, or drop it and log a null.
- **A step** (`S<phase>.<n>`, numbered in order within the phase) is half a day to a day of work: usually one
  idea, one tool or one combination test. Each step has a goal, a "done when" test and a status: `TODO`,
  `IN PROGRESS`, `DONE`, `DROPPED` (with the reason) or `WAITING` (on a decision or a run). Steps of the next
  one or two milestones are written out; later milestones list their ideas and get their steps at the review
  that opens them.
- **Start of a step:** state its goal and its "done when" in the chat; split it further if it looks bigger than
  a day. **End of a step:** a `changelog/` entry (what, why, results), the status updated here and in
  `IDEAS.md` (`DONE (measured gain)` or `NULL (reason)`), and a short review with the user, who can reorder,
  add or drop steps and milestones at every review.
- **Order:** by the priority of `IDEAS.md` (P × score gain ÷ effort), adjusted at each review by what we
  measured. A milestone that delivers more than expected can get more steps; one that stalls is cut at its
  time box.
- **Combinations** (`IDEAS.md` section 8): the tasks are independent (the best version of each is always
  combined). Within a task, ideas that pass are stacked one at a time with a re-verification after each
  addition, and two ideas are tested jointly only when they touch the same part of the circuit or share a
  budget (rotation margin, 20 mHa error budget, a 43-qubit cover). Joint tests are steps of their own.
- **Speed is a deliverable.** Every tool is written for minimal run time (algorithmic complexity first, then
  vectorization and the 12 cores; rule in `CLAUDE.md`), and its measured run time is part of its "done when",
  but never at the cost of exactness: a faster tool must reproduce the reference results before it replaces
  it. The full official verifiers run only at gates; iteration uses static costs, reduced runs and the
  generators' exact proofs.
- **Gates before a submission:** generator build from source (byte-identical rebuild), the generator's exact
  proof where it exists (invsqrt, chemistry loaders), the official verifiers (invsqrt full run, dqi on several
  salts, chemistry), `tools/package.py`, then the user submits.
- **Parallel work:** independent milestones of different tasks can run at the same time in subagents
  (separate git worktrees), when the user agrees; the main session reviews and merges.
- **Experiments are time-boxed.** An idea gets its step; if it shows no measurable gain by the end, it is
  logged as a null and dropped, not stretched.
- **Phase review:** at the end of each phase, a short retrospective (what worked, what to change), a fresh
  board snapshot (section 1), and the next phase's milestones re-ranked from `IDEAS.md`.

## 6. Phase 1 — our own reproducible base, submitted (17 → 20 Nov)

Goal: the public generators rebuilt in our repo, byte for byte, verified, and submitted as our first
score. This is the base all later work forks. **Decision S1.3 (user, 17 Nov): track B only**, no
submission of compiled third-party files. Consequence: our first score (≈ 27.06) starts about 4%
under the public plateau (≈ 28.18); Phase 2 closes that gap by porting the public improvements
into our generators, then goes beyond.

| Step | Goal | Done when | Status |
|---|---|---|---|
| S1.1 | Licensing: `LICENSE`, `NOTICE`, `third_party/PROVENANCE.md` | Files in place, scope stated | DONE (17 Nov) |
| S1.2 | Source inventory: who published what, with or without a generator | `PROVENANCE.md` lists the 5 key sources and their generator status | DONE (17 Nov): only Syed Asad Ali's v12 builders (≈ 27.08 level) ship full source; the 28.18 parts are compiled files or patches of compiled files |
| S1.3 | Decision: track A (credited public frontier) and/or track B (own base from the v12 builders) | User decision recorded here | DONE (17 Nov): **track B only** |
| S1.4 | Unpack Syed Asad Ali's v12 builders into `third_party/`, read them, run their `build.sh` in WSL | Outputs byte-identical to the published SHA-256 | DONE (17 Nov): all 17 files byte-identical (CPython 3.13.18, numpy 2.2.4, WSL; 3 min in all) |
| S1.5 | Verify the v12 outputs with the official verifiers (dfthcsos, full invsqrt, dqi on 8 salts); record our baseline costs | All PASS; baseline table in the changelog | DONE (17 Nov): all PASS; expected score 27.06 (27.08 clean dqi draw, 27.00 worst of 8 salts) |
| S1.6 | Our build entry point `build/` (one command rebuilds all three tasks) and packaging tool `tools/package.py` (builds `submission.zip`, runs the metric-only checks of topic 751951) | The tool rejects a planted bad zip and accepts the good one | DONE (17 Nov): `build/build_all.sh`, `tools/package.py` (self-test rejects 4 planted defects; packaging 24 s, 49 MB zip) |
| S1.7 | First submission (user) and comparison with the local estimate | Kaggle score ≈ 27.0, matching the local estimate | DONE (18 Nov): public score **27.00209**, exactly the local worst-salt projection (2 dqi failed shots on the public split); rank ≈ 1,318 / 3,305 |
| S1.8 | Phase 1 review: retrospective, board snapshot | Plan updated | DONE (18 Nov): see the snapshot below section 1's header and `changelog/2026-11-18-0905-…` |
| S1.9 | Research and idea register: our code (cost anatomy of each circuit) and the state of the art (papers, public notebooks, open-source tools), by four parallel research agents (invsqrt, dqi, chemistry, cross-cutting); merged into `Analysis/IDEAS.md`, as long a list as possible, by task and by category (A: improve our code; B: rewrite with a different technique), ranked by priority = P(success) × gain ÷ effort | `Analysis/research/*.md` and `Analysis/IDEAS.md` written | DONE (18 Nov): 128 ideas (INV 34, DQI 31, CHEM 31, GEN 32) in `Analysis/IDEAS.md`, details in `Analysis/research/` |
| S1.10 | Plan refactor: phases 2–4 rewritten as milestones drawn from `IDEAS.md`, keeping every existing step | Plan refactored and reviewed with the user | DONE (18 Nov), pending the user's review |

## 7. Phase 2 — back to the frontier from our own source, then beyond (21 → 30 Nov)

Goal: **strictly above the public plateau, aiming for the gold line of the day** (on 17 Nov: 28.71, +1.9%
over the plateau), with every gain coming from our own generators. The research (S1.9) shows that the
public frontier's advance over our v12 base exists only as compiled files: rebuilding it in our source is
the safest gain (it brings us from ≈ 27.06 to about the plateau, ≈ +4%), and everything beyond it starts from
there. Milestones are ordered by priority; the expected gains are research estimates, revised as we measure.
The 10 days will not cover every milestone: what is left moves to January at the Phase 2 review.

| Milestone | Ideas (`IDEAS.md`) | Expected score gain | Time box | Go / no-go |
|---|---|---|---|---|
| M2.1 Evaluation toolkit | GEN-A01, A02, A04, A05, A07 (and the former S2.1–S2.4) | enabler (loop 5–10× faster) | 1.5 days | Each tool reproduces the official verifier's costs and verdicts on the v12 baseline before it is used |
| M2.2 invsqrt back to the frontier | INV-A02 → A03 → A04, A10 (chain, in this order) | +2.3% (task −7.0%) | 4 days | Each link: exact proof PASS and cost below the previous link. Target 44 qubits at ≤ 74,560.91, from source |
| M2.3 dqi back to the frontier | DQI-A01, then A02, A03, A04, A12 | +1.3% (task ≈ −4%) | 3 days | Multi-salt expected cost ≤ 3.4641e11 with no more failures than v12 |
| M2.4 Chemistry back to the frontier | CHEM-A03, A07, A06; A02 as a watch (ports only) | +0.3% (task ≈ −0.8%) | 1.5 days | Loader all-branch proof PASS, verifier PASS, cost ≤ 2.94642e8 from source |
| M2.5 Beyond the frontier, round 1 | invsqrt A05, A09, A06, A07, A12, A11; dqi A05, A07, A08, A06; chemistry A08, A04 | +0.5 to +1.5% (several 0.1–0.5% wins) | rest of the sprint | One idea per step, kept if the measured cost drops; stacked one at a time; nulls logged |
| M2.6 Divergent families (seed for January) | GEN-B11 per task (≥ 5 qualitatively different designs for the dominant stage, red-teamed against the nulls) | none yet (prepares B ideas) | 1 day | 1–2 surviving families per task, priced on paper, written into `IDEAS.md` |
| M2.7 Submissions (runs alongside) | — | — | continuous | After each milestone that lifts the score: full gates, then the user submits (about 3 a day at most, 2 kept for safety) |
| M2.8 Freeze before the break | — | — | 0.5 day | `build/` regenerates the best version byte for byte; `STATUS.md` (score, files, open ideas, dead ends) lets a fresh session rebuild and resubmit alone; Phase 2 review done |

### Steps of M2.1 and M2.2 (written out; later milestones get theirs at their opening review)

| Step | Milestone | Goal | Done when | Status |
|---|---|---|---|---|
| S2.1 | M2.1 | GEN-A01: install manwithacat's `kickmix-engine` in WSL and cross-check it against the official verifiers on our three v12 circuits | Same costs and verdicts on all three tasks; measured speed-up and memory recorded (or dropped as a null if it disagrees) | DONE (18 Nov): identical on all fields; dqi 6–7× faster (16–18 s, 0.55 GB vs 110–119 s, 4.7 GB), invsqrt full run 12× (15 s vs 183 s), chemistry no gain (44 vs 42 s) |
| S2.2 | M2.1 | GEN-A02: our static `.kmb` reader (struct-of-arrays, modelled on GrantPrater's `vqcc_ir.py`, credited): qubits, exact expected Toffolis, liveness per stage, in < 1 s | Agrees with the verifier's Toffoli mean on the three baselines | DONE (18 Nov): `tools/kmb.py`; invsqrt 1,511.580 ± 0.008 vs verifier 1,511.570 (0.24 s); chemistry 4,074/pass exact (0.1 s); dqi static 384,003 and 889 Q (0.85 s; its mean stays with the engine) |
| S2.3 | M2.1 | GEN-A04 + former S2.2/S2.3: one gate command per task: invsqrt exact proof (the v12 tree's `prove_any`), dqi on N salts (16 screening, 48 finalists) plus Atakan Aldemir's adversarial inputs, chemistry all-branch proof + verifier; GEN-A07 work dirs on ext4 and a result cache keyed by (SHA, salt) | One command per task gives PASS/FAIL, cost and failure spread on the baseline; a planted error is caught | DONE (18 Nov): `tools/gate.py` (all three tasks in 176 s cold, 1.5 s from cache; planted invsqrt and dqi errors caught). Atakan Aldemir's targeted dqi inputs moved to M2.3 (they matter when the cutoff changes) |
| S2.4 | M2.1 | GEN-A05 + former S2.4: experiment ledger `experiments/LEDGER.md` (one line per attempt, linked to its changelog entry), seeded with the research nulls and manwithacat's negative-results dataset | Ledger exists with the known nulls | DONE (18 Nov): `experiments/LEDGER.md` (8 rows so far), `experiments/NULLS.md` (297 community rows + 95 from our research, by task), generated by `tools/nulls_index.py`. **M2.1 complete** |
| S2.5 | M2.2 | INV-A02: re-synthesise Ali Haydar Özdağ's 64-row table loaders (27 + 29 ANDs) in our generator, S4 measured out, angle refit | 46 qubits at ≈ 78,766 (−1.7%), exact proof PASS | TODO |
| S2.6 | M2.2 | INV-A03 (part 1): stage map and diff of v12 vs Syed's compiled 44-qubit file `47b31838` (which wires leave the multiplier, how v10/v11 and e2/e3 are rebuilt) | Written derivation of the 46 → 44 moves | TODO |
| S2.7 | M2.2 | INV-A03 (part 2): implement the moves as generator patches, qubit cap 44 | 44 qubits at ≈ 75,171 (−4.6%), exact proof PASS | TODO |
| S2.8 | M2.2 | INV-A04 + A10: public micro-edits as passes (dead-write fix, AND merges, CSWAP rings as CCX, rescheduling, angle LP re-solve), rotation bookkeeping | ≤ 74,560.91 from source; gates passed; submission candidate | TODO |

The former Phase 2 steps are all kept, as follows: S2.1 fast evaluator → S2.1–S2.2; S2.2 invsqrt exhaustive
checker and S2.3 dqi robustness harness → S2.3; S2.4 ledger → S2.4; S2.5 port 46 → 44 → M2.2; S2.6 width
profile → done by the research (invsqrt.md §2.2), deepened by GEN-A08 in M3.3; S2.7 error budget → INV-A12
(M2.5) and GEN-B12 (M3.4); S2.8 narrowing to 43 → M3.3; S2.9 wild ideas → M2.6 and M3.4; S2.10 dqi port → M2.3;
S2.11 cutoff sweep → DQI-A03 (M2.3) and DQI-A07 (M2.5); S2.12 GF multiplication cuts → DQI-A02/A13/B08; S2.13
chemistry cost-model study → done by the research (chemistry.md §2); S2.14 loader Toffolis → M2.4/M2.5; S2.15
final weight refit → M3.5; S2.16 submissions → M2.7; S2.17 freeze and S2.18 review → M2.8.

## 8. Break (around 1 Dec → 3 Jan)

No work required. Optional: pass on any question in the forum that would unblock January (host
clarifications), keep the Kaggle token valid, and, if the user wants it, let a scheduled routine run the
frontier watch (GEN-A06: list new public notebooks, re-verify, record what to port) so that January starts
with a ready list.

## 9. Phase 3 — catch-up, then new designs (around 4 → 16 Jan)

Sketched, not planned: its milestones are re-ranked at M3.1 from `IDEAS.md`, `STATUS.md` and the January board.

| Milestone | Ideas (`IDEAS.md`) | Expected score gain | Time box | Go / no-go |
|---|---|---|---|---|
| M3.1 Snapshot and catch-up | GEN-A06; former S3.1–S3.3 | back to the plateau of the day | 2 days | Read `changelog/` and `STATUS.md`, download the board, list the top notebooks, read the December forum; port every public advance we can rebuild from source; resubmit; rewrite sections 1 and 9 |
| M3.2 Evolutionary agent loop | GEN-B01 with GEN-B02 operators, families seeded by M2.6; INV-A08/B01, DQI-B01 | +1 to +3% | 3–4 days | A family beats its incumbent on a task under the full gates; otherwise frozen after 3 rounds without improvement |
| M3.3 invsqrt 43-qubit campaign (joint test) | INV-A20 → A21 + A22 + A24 + A25, GEN-B04, GEN-B03, GEN-A08 | +0.77% per qubit | 3 days | A cover of every 44-wide window within ~39 Toffoli-equivalents; otherwise null |
| M3.4 New designs per task | best surviving B ideas, e.g. INV-B03/B04/B02, DQI-B02/B03/B04, GEN-B12, GEN-B10 | +1 to +5% if one lands | 3 days | Priced prototype beats the incumbent stage; otherwise null with the reason |
| M3.5 Chemistry error-budget campaign (joint test) | CHEM-B05 (budget optimizer) then B01/B02, B03, B04 (C = 43), final weight refit (former S2.15), on a Kaggle P100 | +0.5 to +1% | 2–3 days | All gates (20 mHa, Λ > 49, E_gap) passed with a lower cost |
| M3.6 Team decision (**before 13 Jan, hard deadline**) | — | — | — | Decision recorded (merge only if it brings a task where we are weak) |

## 10. Phase 4 — final submission (17 → 20 Jan)

Only the method is fixed here; the score aimed for is whatever the January board makes realistic.

| Milestone | Date | Goal | Done when |
|---|---|---|---|
| M4.1 Code freeze | Sun 17 | No more circuit changes | Freeze recorded |
| M4.2 Full verifications | Sun 17 | invsqrt exhaustive proof and full verifier, dqi on 48 salts plus adversarial inputs, dfthcsos | All PASS, failure spread known |
| M4.3 Choose the two finals | Mon 18 | (A) safest, fewest dqi failures across salts; (B) best expected score | Both already scored on Kaggle |
| M4.4 Select on Kaggle | Tue 19 | Selection, with a buffer for a failed scoring or a queue delay | Selection confirmed |
| M4.5 Winner package | 17–20 | Generators, `LICENSE`, `NOTICE`, `PROVENANCE.md`, a methodology write-up per task, one command that rebuilds `submission.zip` with the same SHA-256 (rules B.5 and B.7) | Rebuild from a clean clone matches |
| — | Wed 20 | No new submission unless scored and verified before. Deadline 23:59 UTC | — |

## 11. Backlog

The idea backlog is `Analysis/IDEAS.md` (128 ranked ideas). Items from the earlier backlog: the agent loop
(GEN-B01), the dqi working qubits (DQI-A08/A09) and the chemistry tree/QAT co-design (CHEM-B09) are in it;
exploiting impossible x values is a recorded null (invsqrt.md N1: the verifier checks every x). Not yet in it:

- invsqrt: direct synthesis of the phase bits as low-degree Boolean functions of x (to be priced by M2.6).

## 12. Risks

| Risk | Mitigation |
|---|---|
| One task fails on the private salt → the whole submission scores 0 | Exhaustive invsqrt proof; dqi and chemistry checked on many salts; final A is the safest portfolio |
| Part not rebuildable from source → prize ineligible | Provenance table from day 1; replace any file without a generator |
| The compiled-only frontier moves cannot be re-derived in our source (INV-A03 P = 0.5, DQI-A01 P = 0.55) | Time boxes; keep the links that work (e.g. INV-A02 + A04 without A03); move the rest to M3.1 |
| Research estimates are optimistic | Every gain is measured before it counts; `IDEAS.md` statuses record the measured value |
| dqi failed shots on the private salt (−0.15% each) | Rank on 48 salts; prefer robustness over the last 0.1% |
| License breach on reused code | `third_party/` keeps upstream LICENSE/NOTICE; changed files marked; repo private |
| Kaggle queue delay at a deadline | Submit at least 24 h before every deadline |
| Laptop memory (8 GB WSL) | dqi with `--max-workers 2` (or kickmix-engine, < 1 GB); heavy runs on Kaggle CPU notebooks |
| Usage limit stops a session or its agents | `autoContinueAtUsageLimit`; agents save progressively; `changelog/` and this plan let a new session resume |
| Frontier moves away during December | M3.1 catch-up is planned, not optional; optional frontier watch during the break |
