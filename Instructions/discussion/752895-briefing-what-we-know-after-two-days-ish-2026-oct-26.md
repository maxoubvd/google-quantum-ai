# Briefing: what we know after two days (ish) 2026-OCT-26

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752895 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-26 17:45 UTC · 0 votes

[this is a more detailed document that I produced for internal purposes, but I figure it might be useful to someone else‎]

Audience: the user and any agent picking up this project. Purpose: one document
that says what works, what does not, where each task stands against the public
frontier, and what to do next. Every number carries its evidence tier
(AGENTS.md rule 3): *projection*, *local prescreen*, *verifier PASS* (local or
hosted), *observed score* (leaderboard). Public numbers are the authors' reported
official-verifier results unless marked as re-verified by us.

## 1. Position at 17:00 UTC, 2026-10-26

| Item | Value | Tier |
| --- | --- | --- |
| Leaderboard top | **23.35** (Breaking1800); 10 teams ≥ 20, 300 teams ≥ 15, 1,124 teams | observed |
| Best public open portfolio | **15.75** (mrmeows, `qcc-1375-open-circuits`, Apache-2.0, generators included) | author-reported, hosted verifier |
| Best public assembled portfolio | **≈ 19.2** (Alex Chilton's notebook v36: mrmeows invsqrt + jday96314's Schur DQI at 883 Q + michaelrustad/mrmeows chemistry) | author-reported, hosted verifier |
| Our best observed | **9.13908** (2026-10-25) | observed |
| Our best verified, unsubmitted | **11.38437** best-of-public re-baseline (local and hosted verifier PASS) — now a point behind the open plateau; do not spend a submission on it | verifier PASS |
| In progress (separate session today) | credited teaching notebook for the 15.75 release, hosted verification and notebook-origin submission requested by the user (`experiments/2026-10-26-mrmeows-1575-teaching-notebook.md`) | — |
| Submissions today | none yet (quota 5/day, resets 00:00 UTC) | — |

The frontier moved about 25 percent in 24 hours (leader 20.48 → 23.35; best
public 11.3 → ≈ 19.2). Nothing we did in that window moved our score. That is
the single most important fact in this document and section 4 explains it.

## 2. What the competition rewards (the thesis, now with evidence)

The organizers' paper is titled *Leveraging verifiable structure in quantum
circuits to find more efficient compilations with agentic AI*
(arXiv:2610.23832). Its argument: when a bottleneck subroutine is reversible
classical logic plus phases plus measurement-based uncomputation, a classical
simulator can fuzz-test it exactly and fast, and that turns circuit optimization
into a loop an AI coding agent can drive. The ecdsa.fail precedent (over 100
participants, most without quantum background, halving Toffoli counts with coding
agents) is cited as the model. The competition was "designed to accelerate
AI-based research into the optimization of their compilation".

What the frontier did confirms that reading, task by task:

- **DQI.** The decoder in the ≈ 19.2 public portfolio is jday96314's
  continued-fraction Schur decoder, published on GitHub
  (`jday96314/QuantumAI`, `initial_experiments/6/agent14`): the design came out of
  round six of an agent-driven search. Independently, Gu, Kemeny, Yelin and
  Lukin, *Discovering state-of-the-art quantum circuits with large language
  models* (arXiv:2610.22367, open compiler `andigu/uncompute`), report an
  evolutionary LLM framework that produced a Reed–Solomon DQI decoder with the
  same structural moves (truncated Berlekamp–Massey, Forney contracted with the
  mask, Ψ = 1/σ′ mod σ built once per shot).
- **Chemistry.** Alex Chilton's λ-lowering gauge rotation was found by Gemini
  iterating on one function against the official verifier (about 80 s per
  candidate) in an AlphaEvolve-style loop. The factorization gains (8-bit words,
  R = 10 half-size factors) are quantization-aware training runs on GPUs
  against the organizers' own starter-kit trainer.
- **invsqrt.** mrmeows' 63-qubit circuit is a long chain of storage and
  uncomputation rewrites, each checked on all 3,133,452 legal inputs; the
  polynomial and the rotation angles never changed.

So the competition is being won by iterative, verifier-in-the-loop improvement
of classical reversible arithmetic and of a classical factorization, with AI
agents and GPUs as the throughput, not by new quantum algorithms. The "quantum"
content of all three tasks was fixed by the organizers in the reference designs.

## 3. What works

### 3.1 Infrastructure we built (all verified, reusable)

| Asset | What it gives | Evidence |
| --- | --- | --- |
| Native kickmix build, `/Volumes/SSD/Kaggle/tools/kickmix-venv` | Unmodified official verifiers on the Mac; all three reference costs reproduced to the last digit; full three-task check 3–6 min vs ≈ 19 min hosted | verifier parity, hosted cross-check (11.38437 both ways, same ZIP hash) |
| Public repo `github.com/manwithacat/kickmix-macos` v0.1.0 | NEON SIMD word types (1.9× invsqrt, 1.16× DQI over portable), wheels for CPython 3.10–3.14, CI with ASan on Linux, parity tests | observed CI green; benchmarks local |
| `scripts/local_verify_reference.py` | Mirrors the hosted verifier notebook; JSON with hashes | used for the re-baseline |
| `scripts/dqi_rebuild/make_dqi.py` | Rebuilds the paired-v6 DQI decoder from vendored source in 80 s; verifier-equivalent (cost −2.4e-7 relative, shot-sampling noise) | local verifier PASS |
| `scripts/dqi_rebuild/instrument_dqi.py` | Per-stage Toffoli and liveness inventory of any decoder built from that source | static audit, matches `dqi_build.json` |
| `scripts/chemistry_campaign/` pipeline + offline scorer | reproduce / train / quantize / refit / score / build / verify-local; reproduces official chemistry costs to the digit; 20 unit tests | local prescreen, verifier PASS on built circuits |
| Colab A100 launch path (`outputs/colab-chem-2026-10-25/`) | bundle upload, per-shape chunked fetch with hashes, watchdog, budget cap; A100 ≈ 0.30 s per training step | hosted GPU, 9.73 units spent |
| Exact inner solve for a chemistry dictionary (`scripts/chemistry_varpro_gate.py`) | best achievable two-body error for a given dictionary at rank C in 63 s (QR-metric PSD truncation) | local analysis, matches official arithmetic to 1e-8 mHa |
| Score arithmetic from per-task verifier reports | matched Kaggle to five decimals in four of four tests | observed |

### 3.2 Mechanisms that scored or verified (ours)

- **DQI classical-mask projection** (70 positions of the parity phase moved onto
  the syndromes via a Vandermonde solve): −5,559 mean Toffolis at unchanged
  width, verifier PASS, observed 9.13908. It is credited in mrmeows' DQI chain.
  Lesson: the one original gain came from moving work across a representation
  boundary, not from gate-level peepholes.
- **Two-address J = 4 banked ANF lookup** on reference chemistry words: Q 730 →
  682, cost −6.11%, model exact to the verifier. Lesson: the chemistry cost model
  is predictable to the digit; model first, build second.
- **invsqrt one bit less**, **DQI one fewer Euclid step**: valid, tiny,
  superseded by public designs. Lesson: the reference architectures are 13–17×
  behind; single-parameter tweaks to them cannot matter.
- **Phase-only Forney sweep** (toy exhaustive, 11 Toffoli-like gates per position
  measured in kickmix against 85 in paired-v6): correct, and now confirmed by the
  public Schur decoder, whose scan "costs 11 Toffolis per point". We did not
  find the cheap way to build the amplitude polynomial; section 4 covers that.

### 3.3 Mechanisms that work in public (credited; the current plateau)

| Task | Mechanism | Who | Effect |
| --- | --- | --- | --- |
| invsqrt | Replace cubic guess + Newton with one quadratic per segment in centred signed-digit Horner form from a 256-entry QROM | Sumit Chongder, then anazemcev, Rudra Bhaskar, Kirill, Gluzdov, DePie, busyaprime | 260 Q / 7,735 T → ≈ 95 Q / 1,500 T |
| invsqrt | Storage and uncomputation only: recompute QROM outputs where used, erase carry prefixes early by X-measurement with Clifford repair, park input bits in freed slots, replace measurement repair blocks by later Z/CZ corrections | mrmeows | 78 → 63 Q; 108,560 (19.13×) |
| DQI | Continued-fraction (Schur) extraction instead of Berlekamp–Massey; synthesis `(N,B) ← (l_i N + q_i B, N)`; amplitude polynomial `C = A/B′ mod B` so the sweep is a root test plus one CCZ per point; remainder-consuming quotient | jday96314 (agent round 6) | 892 Q, ≈ 535k T; 5.15e11 |
| DQI | Width work on that decoder: write products straight to outputs, borrow provably-zero wires, measure the norm away during scalings | Alex Chilton | 892 → 883 Q; **5.0087e11 (48.8×)** |
| DQI | Ψ = 1/σ′ mod σ once per shot via a half-degree Euclidean pass; Forney contracted with the mask | Gu, Kemeny, Yelin, Lukin (paper); port by Chilton | 1,789 Q after Kozykorpesh's gadgets; 1.72e12 (the width kills it under Q²) |
| DQI | Exact arithmetic replacements in the round12 lineage (banded interpolation, direct gap conversion, odd-code pointer) | mrmeows on Gluzdov's round12 | 870 Q, 884k T; 8.0e11 (30.5×) |
| Chemistry | Exact ANF angle lookups (valid address set closed under bit subsets → each output bit an XOR of address monomials, products erased by measurement); four address loaders; partitioned/nested/shared product lookups | mrmeows, marculera, Özer | T 17,063 at 10 bits |
| Chemistry | **8-bit words**: the cost model's state-preparation workspace drops from 607 to ≈ 300 qubits at brot = 8 (not at 9); recover accuracy by constrained discrete refit with coefficient transport and a larger dictionary (B ≈ 300–350), accept λ_eff ≈ 47.5–49.5 | mrmeows; michaelrustad independently | Q 742 → 448; 5.0e8 (6.69×) |
| Chemistry | Half-size factorization (R, B, C) = (10, 129, 60) compiled with the exact lookup compiler at 3 words per circuit | michaelrustad (factor), mrmeows (compiler), Chilton (combination) | 415 Q, 21,922 T; **4.37e8 (7.68×)** |
| Chemistry | Right-rotation gauge of W lowering the column one-norm; Gemini-iterated optimizer | Magmucot (base), Chilton/Gemini | −0.2 to −0.55% per factor; bound ≈ 0.5% total |

## 4. What does not work (closed lines and the lesson from each)

Grouped; the ledger (`docs/strategy-ledger-2026-10-25.md`) has every card.

**Reference-architecture tweaks.** Bit trimming without refit, QROM batch sweeps,
row reordering, bit-plane factoring, symmetry scans, one-bit/one-step changes:
all valid or closed, all under 2%. The reference designs are the wrong starting
point; the public designs are 13–17× ahead on two tasks.

**Chemistry fixed-angle refits and gauge steps on public factors.** Gradient and
projected-gradient steps on `H2_I`, weakest-copy pruning or scaling, λ_eff refit
with fixed words: 1e-6 to 1e-3 effects. The public factors sit on cusps (zero
eigenvalues, zero weights) and the gauge bound is ≈ 0.5% (mrmeows derives it from
Λ₂ ≥ ¼ Σ|AAᵀ|). Lesson: Λ is a one-norm; the lever is the factorization and the
word width, not the gauge.

**Chemistry structure shrinking at R = 1.** Our Colab campaign trained pruned
shapes (1,460,384) and (1,496,300): they stall at 20.5–22.4 mHa continuous
against the public point's 17.07, and quantization adds ≈ 3 mHa, so they fail the
20 mHa gate. The shared-basis diagnostic found no common eigenbasis (pairwise
commutators barely better than random). The limbo (variable-projection)
dictionary search showed the frontier's weights are within 0.5 mHa of optimal for
their dictionary and that first-order dictionary descent is slower than the joint
Gauss–Newton trainer. **All three conclusions were correct and all three missed
the point**: the public frontier did not shrink B at R = 1 with 10–11-bit words;
it dropped to 8-bit words to cross the workspace cliff, grew B back to ≈ 300 to
recover accuracy, and separately moved to R = 10 with small B and C. Our
Eckart–Young floor (R·C ≥ 290 at 20 mHa) is consistent with that: (10, 129, 60)
has R·C = 600. We searched the wrong axis, and we never looked for cliffs in the
cost model.

**DQI graft and space ideas.** Separate reciprocal bank (Q²T 1.62× worse), Dialog
replay (70 live words vs 35), BM-state replacement, gap-flag recomputation,
cheaper GF(2¹²) multiplier (42 is best known; floor 34; a 14-minute flip-graph
probe found no 41). Streaming Euclid for h = ω/σ′ mod σ: correct but 79 scratch
words. Lesson: under R·Q²·T with 840 of 881 qubits pinned by the input, every
scratch word costs 0.23% × 12 and nothing that adds a bank survives.

**Where our DQI analysis was right but late.** We proved at toy size that a
phase-only sweep costs 11 per position and saves ≈ 27% of T, and we wrote on the
card that the one open problem was computing the amplitude polynomial with O(1)
scratch, with "build h inside the conversion stage that already constructs σ and
ω from the continued fraction" as the first untested idea. jday96314's decoder is
exactly that: continued-fraction synthesis producing `C = A/B′ mod B` as it
goes. The design we needed existed on GitHub while we prototyped. Lesson: check
the public code of the people at the top of the board before prototyping, every
day, not once at the start.

**Process.** Hosted-only verification capped us at a few candidates per day
until the local build landed; Colab single-file downloads hang on large files
(use chunked fetch); gate thresholds must be calibrated against the frontier's
own number (our 15 mHa gate was stricter than the frontier's 17.07).

## 5. How close is each task to done?

Tiers: public authors' measurements, our local audits, and simple arithmetic.

**invsqrt: close to the floor of its architecture.** 63 qubits, 1,404 mean
Toffolis, 19.13×. mrmeows reports the error budget 95–99% used throughout
(worst 0.99966 of tolerance), Toffolis bounded below by ≈ 1,400–1,450 inside
this design, and that Newton refinement, cubic segments, other table widths and
register recolouring are neutral or worse. One qubit is worth ≈ 1.5%. The only
stated open direction is computing phase bits directly as low-degree functions
of x (an adaptive model with 133 table rows passes the numerical check at 0.922 of
tolerance but has no circuit yet). Expect 20–21×, not 30×.

**DQI: the biggest mover, with structural headroom left.** 883 Q, 537k T,
48.8×. Width floor is 840 (the syndrome bank), so the Q² term can give at most
(883/840)² = 1.105 more. Toffolis: the sweep is now 11 per point (45k); the rest
is extraction and synthesis forward and reverse, and reverse costs 31 more
Toffolis per node pair than forward (phase repayment for norm measurement). At
850 Q and 400k T the ratio would be ≈ 64; nobody is there. The remaining gains
are in the compute–use–reverse structure of the extraction, which is exactly the
kind of thing the agentic loop finds. Multiplier is at the best known rank (42)
everywhere.

**Chemistry: the least optimized task and the most machine-learning-shaped.**
415 Q, 21,922 T, λ_eff ≈ 48, 7.68×. The public authors' own open directions are
all training problems: a smaller 8-bit dictionary (B = 255 stalls at 48 mHa and
must reach 20), lattice rounding on the linearised sensitivity (LLL/Babai),
annealed quantization-aware training, and the next workspace cliff at brot ≤ 7
(temporary qubits 543 → 272 in the cost model) which currently fails accuracy at
60–70 mHa. At Q = 448, 7× needs T ≈ 21,500; the (10,129,60) factor already gets
7.68× at 415 Q. Because the score is a geometric mean, a point of improvement is
worth the most here, and this is where GPU hours and a trainer matter.

## 6. What to do next (recommendation)

1. **Adopt the plateau, with credit.** Finish the 15.75 teaching notebook and
   submission already in flight. Then assemble and locally verify the ≈ 19.2
   portfolio (mrmeows invsqrt, jday96314/Chilton Schur DQI, michaelrustad/mrmeows
   chemistry) through our local verifier before any hosted run. That is one
   evening's work with the tools we have and it puts us at the public frontier.
2. **Put the chemistry pipeline on the 8-bit problem.** Our campaign code,
   offline scorer and exact inner solve are the right tools; the target was wrong.
   Preregister: quantization-aware refit at brot = 8 for B in 255–303 with the
   starter-kit trainer, A100, cap 40 units, gate = error₂ ≤ 20 mHa on decoded
   words and cost below 4.37e8. Then the brot ≤ 7 cliff as a second card.
3. **DQI: rebuild jday96314's Schur decoder in our 80-second loop** and run the
   inventory on it. Our width-and-liveness audit is the tool that found Chilton's
   9 qubits; the middle stages are now the peak and nobody has cut them.
4. **Publish**: the two discussion drafts in `docs/comms/`, the strategy
   communique (`docs/comms/2026-10-26-kaggle-strategy-communique.md`), and the
   chemistry toolkit plus negative results as a Kaggle dataset and notebook.

## 7. Open items and housekeeping

- AddressSanitizer start-up deadlock of the Bazel `kickmix_test` binary on macOS
  26 / Apple clang 17: investigation running (user's request, 2026-10-26); the
  stuck process from yesterday (1,346 CPU-minutes) was killed today.
- Nothing is committed in the project repo since `be45f95`; about 35 files
  untracked.
- `outputs/dqi-rebuild-2026-10-25/build-seed{0,1,2}`: three ≈ 1 GB duplicate
  rebuilds kept under the no-delete rule; the user may remove them.
- Colab: no session running; 9.73 of 150 units spent; GCS bucket not created.
- Quota: five submissions available today; the 11.384 re-baseline is no longer
  worth one.

## 8. Reading order for a new agent

`AGENTS.md` → this briefing → `docs/strategy-ledger-2026-10-25.md` →
`docs/strategy-review-2026-10-25.md` → the open cards in `experiments/`.
Public sources to re-read before any new card: mrmeows' dataset README
(`outputs/mrmeows-15-75-review-2026-10-26/README.md`, pinned), Alex Chilton's
notebook (pulled to the session scratchpad; pin before quoting), the
`jday96314/QuantumAI` repository, arXiv:2610.22367 and arXiv:2610.23832.


---

## Comments (3 threads)

Replies are indented with `>`, one level per reply depth.

**James Day** · 2026-10-26 18:23 UTC · 3 votes

> The decoder in the ≈ 19.2 public portfolio is jday96314's continued-fraction Schur decoder, published on GitHub (jday96314/QuantumAI, initial_experiments/6/agent14)

This repo was public by accident. Of course, it is my fault for not paying attention to the permissions, but it is quite creepy what Kagglers and their agents are able to find without the author trying to share 😅

> **manwithacat** — **Topic author** · 2026-10-26 18:27 UTC · 1 votes
>
> if it helps, I didn't understand it. I have degrees in mathematics and electronic engineering, but I don't even remember what I've forgotten anymore.

> **alex chilton** · 2026-10-27 10:55 UTC · 1 votes
>
> oh i am so sorry dude :-( i sent a claude agent off on a hail mary for ideas! he came back with it - refused to use it verbatim! i said ok can you use the idea and build something similar - really really sorry...if its any consolation he was already playing with schur - others improved my attempt too :-)

**alex chilton** · 2026-10-27 10:57 UTC · 0 votes

pretty amazing how much progress in such a short time! i think end of the month will be road blocks...

**hengck23** · 2026-10-27 03:30 UTC · 0 votes

For the chemistry task, how relevant is bit-aware training? Can we use 9bit BROT?

> **hmmmreally** · 2026-10-29 21:14 UTC · 0 votes
>
> With regard to the 9-bit / bit-aware question, I looked at 9-bit BROT through the cost model and I think it's the wrong approach. The inner state-prep workspace comes from the verifier’s k2 rule, and at brot 8 the threshold (brot × (N−1) = 600) just misses k2 = 5 (601), so the inner-prep temp stays ≈304 qubits and Q ≈430. At brot 9 the threshold rises enough to admit k2 = 5; the forced k2=5 what-if at the current shape cuts Toffolis by roughly 10% but pushes Q to about 733, so T·Q gets about 54% worse — and that is before counting the wider word register. I could of course be mistaken, and haven’t built a brot-9 submission, so treat that as a model result and not a ground truth.
>
> On narrower words: I tried a few routes on the public 8-bit word table behind the ~9.83 builds (with my own weight refit on top). Simply rounding to 7 bits lands far outside the walls (error_1 ≈2.2 Ha, error_2 ≈0.43 Ha vs limits of 0.02). Quantization-aware training with exact gradients (STE plus greedy bit flips, λ unconstrained) stalled at roughly 0.40 Ha on error_2 and 1.2 Ha on error_1 — still about 20–60× over. The 8-bit words also look fragile: flipping a single random LSB already pushes error_2 to about 30× the limit, so they behave like a discrete joint fit. Scope note: this covers rounding existing words and local search around them on that earlier table, not Atakan’s newer one — read it as “these routes fail,” not “7-bit is impossible.” A from-scratch retrain would still be the only plausible route at any width.
>
> It is entirely possible this post contains errors, but my impression is that you can use this information to speed up your process of narrowing down best paths.
