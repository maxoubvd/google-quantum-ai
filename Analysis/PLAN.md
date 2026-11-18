# Action plan: Google Quantum AI — Verifiable Quantum Circuit Compilation

Written 2026-11-17, **updated 2026-11-17 20:00 UTC** with the live leaderboard (Kaggle API, 3,274 teams)
and the top public notebooks. Other sources: `Instructions/` (overview, data, rules, forum snapshot of
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
| 2 | Strictly above the plateau, aiming for the gold line (top 10 + 0.2%) |
| 3 | Back above the plateau after the December rise, then as close to gold as possible |
| 4 | Best private-board expectation: robustness on the private salt first, then score |

## 1. Where the competition stands (2026-11-17)

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
2. ~~Rules accepted~~ Done: the account is entered in the competition (no submission yet).
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

## 5. How we work: one step at a time

The work is split into small **steps** (half a day to a day each), worked one at a time, in order,
like the stories of an agile sprint.

- **Each step has** an ID (`S<phase>.<n>`), a goal, a "done when" test and a status:
  `TODO`, `IN PROGRESS`, `DONE`, `DROPPED` (with the reason) or `WAITING` (on a decision or a run).
- **Start of a step:** state its goal and its "done when" in the chat; split it further if it looks
  bigger than a day.
- **End of a step:** a `changelog/` entry (what, why, results), the status updated here, and a short
  review with the user: what we got, what we learned, and which step comes next. The user can reorder,
  add or drop steps at every review.
- **Speed is a deliverable.** Every tool is written for minimal run time (algorithmic complexity
  first, then vectorization and the 12 cores; rule in `CLAUDE.md`), and its measured run time is
  part of its "done when", but never at the cost of exactness: a faster tool must reproduce the
  reference results before it replaces it. The full official verifiers run only at gates; iteration uses static
  costs, reduced runs and the generators' exact proofs.
- **Experiments are time-boxed.** An idea gets one step; if it shows no measurable gain by the end, it is
  logged as a null and dropped, not stretched.
- **Backlog:** ideas that come up mid-step go to the backlog (section 11) instead of interrupting the
  step in progress.
- **Phase review:** at the end of each phase, a short retrospective (what worked, what to change), a
  fresh board snapshot (section 1), and the next phase's steps rewritten.

## 6. Phase 1 — our own reproducible base, submitted (17 → 20 Nov)

Goal: the public generators rebuilt in our repo, byte for byte, verified, and submitted as our first
score. This is the base all later work forks. **Decision S1.3 (user, 17 Nov): track B only**, no
submission of compiled third-party files. Consequence: our first score (≈ 27.08) starts about 1.1
points under the public plateau (≈ 28.18); Phase 2 closes that gap by porting the public improvements
into our generators, then goes beyond.

| Step | Goal | Done when | Status |
|---|---|---|---|
| S1.1 | Licensing: `LICENSE`, `NOTICE`, `third_party/PROVENANCE.md` | Files in place, scope stated | DONE (17 Nov) |
| S1.2 | Source inventory: who published what, with or without a generator | `PROVENANCE.md` lists the 5 key sources and their generator status | DONE (17 Nov): only Syed Asad Ali's v12 builders (≈ 27.08 level) ship full source; the 28.18 parts are compiled files or patches of compiled files |
| S1.3 | Decision: track A (credited public frontier) and/or track B (own base from the v12 builders) | User decision recorded here | DONE (17 Nov): **track B only** |
| S1.4 | Unpack Syed Asad Ali's v12 builders into `third_party/`, read them, run their `build.sh` in WSL | Outputs byte-identical to the published SHA-256 | DONE (17 Nov): all 17 files byte-identical (CPython 3.13.18, numpy 2.2.4, WSL; 3 min in all) |
| S1.5 | Verify the v12 outputs with the official verifiers (dfthcsos, full invsqrt, dqi on 8 salts); record our baseline costs | All PASS; baseline table in the changelog | DONE (17 Nov): all PASS; expected score 27.06 (27.08 clean dqi draw, 27.00 worst of 8 salts) |
| S1.6 | Our build entry point `build/` (one command rebuilds all three tasks) and packaging tool `tools/package.py` (builds `submission.zip`, runs the metric-only checks of topic 751951) | The tool rejects a planted bad zip and accepts the good one | DONE (17 Nov): `build/build_all.sh`, `tools/package.py` (self-test rejects 4 planted defects; packaging 24 s, 49 MB zip) |
| S1.7 | First submission (user) and comparison with the local estimate | Kaggle score ≈ 27.0, matching the local estimate | TODO |
| S1.8 | Phase 1 review: retrospective, board snapshot | Plan updated | TODO |
| S1.9 | Research and idea register: our code (cost anatomy of each circuit) and the state of the art (papers, public notebooks, open-source tools), by four parallel research agents (invsqrt, dqi, chemistry, cross-cutting); merged into `Analysis/IDEAS.md`, as long a list as possible, by task and by category (A: improve our code; B: rewrite with a different technique), ranked by priority = P(success) × gain ÷ effort | `Analysis/research/*.md` and `Analysis/IDEAS.md` written | DONE (18 Nov): 128 ideas (INV 34, DQI 31, CHEM 31, GEN 32) in `Analysis/IDEAS.md`, details in `Analysis/research/` |
| S1.10 | Plan refactor: phases 2–4 rewritten as milestones, each milestone = implementing and testing one idea or a group of ideas from `IDEAS.md` (with go/no-go criteria), keeping every existing step (tooling, submissions, freeze, January catch-up, final choice) as or inside milestones; combinations tested within a task by stacking the ideas that pass, plus joint tests where two ideas touch the same part of the circuit | Plan refactored and reviewed with the user | TODO |

Submit by 19 Nov evening at the latest; never miss the 24 h buffer before a deadline.

## 7. Phase 2 — improvement sprint (21 → 30 Nov)

Goal: **strictly above the public plateau, aiming for the gold line of the day** (on 17 Nov: 28.71,
+1.9% over the plateau), with every gain coming from our own generator changes. Effort split by
exchange rate, revised as the board moves: **invsqrt ~50%, dqi ~25%, chemistry ~25%** to start.
If a new public jump lands during the sprint, rebase on it first (a new step), then keep improving.

### Tooling steps (first)

| Step | Goal | Done when |
|---|---|---|
| S2.1 | Fast evaluator per task (`tools/eval_*.py`): static cost from `operation_statistics()` and qubit count in seconds; reduced verifier run in about a minute | One command per task prints cost, qubits and Toffolis; agrees with the full verifier on the baseline |
| S2.2 | invsqrt exhaustive checker: every x and every measurement branch | Proves the baseline; catches a planted one-input error |
| S2.3 | dqi robustness harness: N salts (16 screening, 48 finalists) plus Atakan Aldemir's adversarial low-weight inputs | Expected cost and failure spread for the baseline |
| S2.4 | Experiment ledger `experiments/LEDGER.md` (one line per attempt, linking its changelog entry), seeded with the forum nulls (topic 753185, manwithacat's negative-results dataset, alex chilton's notebook) | Ledger exists; known nulls listed |

### invsqrt steps (main track)

| Step | Goal | Done when |
|---|---|---|
| S2.5 | Port the public 46 → 44 qubit improvements (64-row table, rescheduling, angle refit) into our generator | Our build reaches 44 qubits at ≤ the public cost, from source |
| S2.6 | Width profile: map live qubits per stage, list every point where 44 are live | Profile and a ranked list of narrowing candidates |
| S2.7 | Error budget: re-split the 2^-24/√x tolerance (table, divider truncation, Z_POW precision) | Measured cost of each split; best kept or logged as null |
| S2.8 | Narrowing attempts (measure-then-rebuild on remaining peaks, table re-synthesis): one candidate per step | 43 qubits within ~39 extra Toffoli-equivalents, or null |
| S2.9 | One wild idea per step (low-degree phase-bit synthesis, impossible x values): time-boxed to one day | Gain measured or null logged |

### dqi steps

| Step | Goal | Done when |
|---|---|---|
| S2.10 | Port the public 889 → 890 improvements (cutoff, retained products, GF(2^24) closings) into our pipeline from source | Our decoder at ≤ the public expected cost |
| S2.11 | Failure budget: sweep the extractor cutoff against multi-salt expected cost | Best cutoff by expected cost and by robustness |
| S2.12 | Toffoli cuts on GF multiplications and retained products (topic 756476) | Gain measured or null |

### Chemistry steps

| Step | Goal | Done when |
|---|---|---|
| S2.13 | Cost-model study: what makes up T_total = 14,706 and Q_total = 427, which knob moves which term | A written breakdown; a list of terms nobody has attacked |
| S2.14 | Loader Toffolis (4,022/pass): gauges, don't-cares, ANF banking, tree order | Gain measured or null |
| S2.15 | Final weight refit at the 20 mHa edge (P100 on Kaggle) | Refit applied to our best loaders |

### Closing steps

| Step | Goal | Done when |
|---|---|---|
| S2.16 | Submissions: about 3 a day for verified candidates, 2 kept for safety (runs alongside the steps above) | Every candidate fully verified before it is sent |
| S2.17 | Freeze before the break: `build/` regenerates the best version byte-for-byte; `STATUS.md` (score, files, open ideas, dead ends) | A fresh session can rebuild and resubmit from `STATUS.md` alone |
| S2.18 | Phase 2 review: retrospective and backlog triage for January | Plan updated |

The order inside each track is a default. At each review we pick the next step by expected gain on the
score, so a promising track can take more steps than its share.

## 8. Break (around 1 Dec → 3 Jan)

No work required. Optional: pass on any question in the forum that would unblock January (host
clarifications), and keep the Kaggle token valid.

## 9. Phase 3 — catch-up and new ideas (around 4 → 16 Jan)

This phase is sketched, not planned: its steps depend on where the board and the public frontier are
in January. S3.1 rewrites them.

| Step | Goal | Done when |
|---|---|---|
| S3.1 | Snapshot: read `changelog/` and `STATUS.md`, download the board, list the top notebooks, read the December forum; rewrite sections 1 and 9 of this plan | Plan updated with January numbers and steps |
| S3.2 | Catch-up: pull every public part better than ours; rebase our generators on the new frontier | Our build at the new public frontier |
| S3.3 | Resubmit to recover our position | Score at or above the plateau of the day |
| S3.4… | New-idea steps, one per backlog item, highest expected gain first | Gain measured or null |
| S3.x | **Before 13 Jan (hard deadline):** team decision (merge only if it brings a task where we are weak) | Decision recorded |

## 10. Phase 4 — final submission (17 → 20 Jan)

Only the method is fixed here; the score aimed for is whatever the January board makes realistic.

| Step | Date | Goal | Done when |
|---|---|---|---|
| S4.1 | Sun 17 | Code freeze | No more circuit changes |
| S4.2 | Sun 17 | Full verifications: invsqrt exhaustive, dqi on 48 salts, dfthcsos | All PASS, failure spread known |
| S4.3 | Mon 18 | Choose the two finals: (A) safest, fewest dqi failures across salts; (B) best expected score | Both already scored on Kaggle |
| S4.4 | Tue 19 | Select them on Kaggle; buffer for a failed scoring or a queue delay | Selection confirmed |
| S4.5 | 17–20 | Winner package: generators, `LICENSE`, `NOTICE`, `PROVENANCE.md`, a methodology write-up per task, one command that rebuilds `submission.zip` with the same SHA-256 (rules B.5 and B.7) | Rebuild from a clean clone matches |
| — | Wed 20 | No new submission unless scored and verified before. Deadline 23:59 UTC | — |

## 11. Backlog

Ideas not yet scheduled as steps. Each review may promote one.

- invsqrt: direct synthesis of phase bits as low-degree Boolean functions of x.
- invsqrt: exploit x values that cannot occur in the algorithm (topic 753451; check the host answer).
- dqi: reduce the 50 working qubits above the 840 syndrome qubits (one qubit ≈ 830 Toffolis).
- chemistry: co-design the Givens tree order and the QAT training so the XOR word differences compress.
- Infrastructure: an agent loop that evolves the generators against the evaluator, several families in
  parallel (Gu et al., arXiv:2610.22367).

## 12. Risks

| Risk | Mitigation |
|---|---|
| One task fails on the private salt → the whole submission scores 0 | Exhaustive invsqrt proof; dqi and chemistry checked on many salts; final A is the safest portfolio |
| Part not rebuildable from source → prize ineligible | Provenance table from day 1; replace any file without a generator |
| dqi failed shots on the private salt (−0.15% each) | Rank on 48 salts; prefer robustness over the last 0.1% |
| License breach on reused code | `third_party/` keeps upstream LICENSE/NOTICE; changed files marked; repo private |
| Kaggle queue delay at a deadline | Submit at least 24 h before every deadline |
| Laptop memory (8 GB WSL) | dqi with `--max-workers 2`; heavy runs on Kaggle CPU notebooks |
| Frontier moves away during December | Catch-up days in January are planned, not optional |
