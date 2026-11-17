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

## 5. Phase 1 — a first good, eligible submission (17 → 20 Nov)

Goal: a submission **at the public plateau** (≈ 28.18 on 17 Nov, bronze-to-silver zone) whose three parts
are rebuilt from source by our own scripts, verified locally at several salts, and documented for
eligibility. This is the base everything else forks.

| Day | Work | Done when |
|---|---|---|
| Tue 17 | `LICENSE`, `NOTICE`, `third_party/` skeleton (section 4.1). Pull the sources behind the three best files: Syed Asad Ali's dataset `syed3000/vqcc-own-circuit-portfolio` (invsqrt and dqi builders), maoxin123's notebooks (44-qubit invsqrt, 4,022-Toffoli loaders), GrantPrater's `grantprater/vqcc-dqi-61c8348d-merged-closings` and notebook, Balram Srinivasan's and boristown's dqi patches, Atakan Aldemir's `atakanaldemir/qcc-shared-4022-loaders-and-own703-weights`, Lera's `marculera/dfthc-c44-build`, prvsiyan's refit code. Asterios Panoras's QCC Circuit Atlas and manwithacat's circuits dataset serve as SHA-256 cross-checks. | Sources in `third_party/` with `PROVENANCE.md`; for each best file we know whether a generator exists |
| Wed 18 | Run each generator in WSL and check the output SHA-256 against the published one. Wrap them in our own `build/` scripts. dfthcsos verifier (1 min); full invsqrt verification (1 h). | Byte-identical rebuilds; dfthcsos PASS at the expected cost |
| Thu 19 | dqi on the default salt plus 7 extra salts (`--max-workers 2`); pick the dqi file with the fewest failed shots across salts, not the cheapest on one draw. `tools/package.py` builds `submission.zip` and runs the metric-only checks of topic 751951 (one copy of each task, `angle_k_*` without gaps, named registers, little-endian `.npz` arrays, `max(qubits, bits) − 1 < 699,050`). **Submit (you) by Thursday evening.** | Kaggle score at the plateau of the day, matching the local estimate |
| Fri 20 | Buffer for queue delays or a failed submission. If time remains, swap in any newer public part that beats ours and is rebuildable from source. | Phase 1 base frozen |

Decision rule: if a part cannot be rebuilt from source by 19 Nov, submit the published file anyway to get
a score, and list it in `PROVENANCE.md` as "not yet eligible" to fix in Phase 2. dqi lines built by
patching a compiled file need special care: the patch script plus the pinned input must reproduce the
file (topic 752010 asks the hosts exactly this; no answer yet).

## 6. Phase 2 — improvement sprint (21 → 30 Nov)

Goal: **strictly above the public plateau, aiming for the gold line of the day** (on 17 Nov: 28.71,
+1.9% over the plateau), with every gain coming from our own generator changes. Effort split by
exchange rate, revised as the board moves: **invsqrt ~50%, dqi ~25%, chemistry ~25%** to start.
If a new public jump lands during the sprint, rebase on it first, then keep improving on top of it.

### 6.1 Infrastructure first (21–22 Nov)

- **Fast evaluator per task** (`tools/eval_*.py`): static cost from `operation_statistics()` and qubit
  count in seconds; a reduced verifier run in about a minute (`--shots-per-x 1 --batch-size 12`,
  `--num-shots`); the full verifier only before a submission.
- **invsqrt exhaustive checker**: every x and every measurement branch, as Syed Asad Ali does, so that
  a local PASS is a proof and the private salt cannot fail us.
- **dqi robustness harness**: N salts (16 for screening, 48 for finalists), ranked by mean cost; plus the
  adversarial low-weight inputs published by Atakan Aldemir.
- **Experiment ledger** (`experiments/LEDGER.md`): one line per attempt with the idea, the measured
  cost and the verdict, including nulls. Check the forum null lists first (topic 753185, manwithacat's
  negative-results dataset, alex chilton's gold notebook) to avoid repeating them.
- **Agent loop**: Claude Code works on the generators against the evaluator, several families in
  parallel (Gu et al., arXiv:2610.22367). One task per session, with the ledger as shared memory.

### 6.2 invsqrt (main track)

Starting point: the 44-qubit line (in-place square root, table plus correction, 24-row divider, phase per
quotient bit, measure-then-rebuild at width peaks, 33 fine rotations, worst error/tolerance 0.99956).

1. Map the width profile per stage and list every point where 44 qubits are live. A 43-qubit build is
   worth +0.77% and can afford up to ~39 extra Toffoli-equivalents.
2. Error budget: the margin is nearly used up (0.99956), so any saving must come with a re-split of
   the 2^-24/√x tolerance between the table, the divider truncation and the Z_POW synthesis precision.
   Fewer divider rows or a coarser table with a better correction may pay.
3. Measure-then-rebuild on the remaining peaks; table re-synthesis (ANF, don't-cares) on the loaders.
4. Wilder ideas, one day at most each: direct synthesis of the phase bits as low-degree Boolean functions
   of x; exploiting x values that cannot occur (topic 753451, check the host answer first). The 31.4
   teams suggest a structural change exists; invsqrt is the likeliest place.

### 6.3 dqi

1. Robustness first: a decoder that fails on fewer salts is worth more on the private board than one
   that is 0.1% cheaper.
2. Failure budget: sweep the extractor cutoff and other "rare case" branches against the multi-salt
   expected cost.
3. Toffoli cuts on the GF(2^12)/GF(2^24) multiplications and retained products (topic 756476 for exact
   multiplication complexities).
4. Width: 50 working qubits above the 840 syndrome qubits; one qubit is worth ~830 Toffolis.

### 6.4 Chemistry

1. Read `scoreinfra/costs_prepare.py` and `SCORING_README.md` until we can say exactly what makes up
   T_total = 14,706 and Q_total = 427, and which term each knob (R, B, C, brot, J) moves. Look for a term
   that the forum has not attacked ("the verifier's own cost model fixes most of the qubits").
2. Loader Toffolis (4,022 per pass): gauge choices (±θ, π − θ), dead addresses as don't-cares, ANF
   banking, joint choice of tree order and word coding.
3. Weight refits are already at the 20 mHa edge (gates at 0.0199999985); only as final polish.

### 6.5 Submissions and wrap-up

- Use about 3 submissions a day for real candidates and keep 2 for safety. Submit a candidate only after
  a full local verification.
- **Wrap-up before the break (around 29–30 Nov).** Freeze the best version, make sure `build/`
  regenerates it byte-for-byte, write `STATUS.md` (score, files, open ideas, dead ends) so January
  starts in minutes.

## 7. Break (around 1 Dec → 3 Jan)

No work required. Optional: pass on any question in the forum that would unblock January (host
clarifications), and keep the Kaggle token valid.

## 8. Phase 3 — catch-up and new ideas (around 4 → 16 Jan)

This phase is sketched, not planned: its content depends on where the board and the public frontier
are in January. Its first step is to rewrite sections 1, 6 and 8 of this plan from a fresh snapshot.

1. **Catch-up (2 days).** Read the forum and the December briefings, list the public notebooks by score
   (`kaggle kernels list --competition … --sort-by scoreDescending`), download the leaderboard, and pull
   every public part better than ours. Rebase our generators on the new frontier. Resubmit to recover our
   position.
2. **New ideas.** Pick the 2–3 highest-value ideas from the `STATUS.md` backlog and the January forum.
   Keep the same infrastructure and ledger.
3. **Before 13 Jan (hard deadline): team decision.** Decide whether to merge with another team (only
   worth it if it brings a task where we are weak). After 13 Jan no merger and no new entrant.
4. Last improvement round, still on full verification.

## 9. Phase 4 — final submission (17 → 20 Jan)

Only the method is fixed here (freeze, verify, choose two); the score aimed for is whatever the January
board makes realistic.

| Date | Work |
|---|---|
| Sun 17 | Code freeze. Full verifications: invsqrt exhaustive, dqi on 48 salts, dfthcsos. |
| Mon 18 | Choose the **two final submissions**: (A) safest, fewest dqi failures across salts; (B) best expected score. Both already scored on Kaggle. |
| Tue 19 | Select them on Kaggle. Last-chance buffer for a failed scoring or a queue delay. |
| Wed 20 | No new submission unless one was scored and verified before. Deadline 23:59 UTC. |

In parallel: the winner package. Apache-2.0 repository with the generators, `LICENSE`, `NOTICE`,
`third_party/PROVENANCE.md`, a methodology write-up per task and one command that rebuilds
`submission.zip` with the same SHA-256 (rules B.5 and B.7).

## 10. Risks

| Risk | Mitigation |
|---|---|
| One task fails on the private salt → the whole submission scores 0 | Exhaustive invsqrt proof; dqi and chemistry checked on many salts; final A is the safest portfolio |
| Part not rebuildable from source → prize ineligible | Provenance table from day 1; replace any file without a generator |
| dqi failed shots on the private salt (−0.15% each) | Rank on 48 salts; prefer robustness over the last 0.1% |
| License breach on reused code | `third_party/` keeps upstream LICENSE/NOTICE; changed files marked; repo private |
| Kaggle queue delay at a deadline | Submit at least 24 h before every deadline |
| Laptop memory (8 GB WSL) | dqi with `--max-workers 2`; heavy runs on Kaggle CPU notebooks |
| Frontier moves away during December | Catch-up days in January are planned, not optional |
