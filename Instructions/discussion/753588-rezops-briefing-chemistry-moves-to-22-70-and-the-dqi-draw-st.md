# RezOps Briefing: chemistry moves to 22.70, and the DQI draw sticks to the circuit — 2026-10-28

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753588 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-28 17:16 UTC · 1 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps runs the background work of
our research programme: it scans new public notebooks, re-verifies their circuits, maintains our public circuits
dataset and the notebook that blends the best of it, and drafts this briefing. A human (me) sets the direction,
reviews and approves every post. Numbers are labelled: *observed* (leaderboard), *author-reported* (the author's
official-verifier run), or *re-verified by us* (the unchanged official verifier, run locally at its default salt).
Window: 27 Oct 17:00 UTC to 28 Oct 17:00 UTC. Corrections welcome.

## 1. Leaderboard (28 Oct, 17:13 UTC, observed)

| | 27 Oct | 28 Oct |
| --- | --- | --- |
| Top score | 24.53 | **26.14** (Breaking1800) |
| Teams ≥ 25 | 0 | 2 |
| Teams ≥ 23 | 7 | **25** |
| Teams ≥ 22.5 | — | 370 |
| Teams ≥ 20 | 481 | 734 |
| Teams | 1,275 | 1,437 |

The top ten are 23.72–26.14 (Breaking1800, Jack, Z³uperposition, Thinking..., Schrödinger's cat lived, Bingran You,
Distortion!!, Kevin Bogner, Danial Motlagh, VAGO Superposition). The median team is at 20.47.

**Approximate medal lines today** (Kaggle's rules for 1,000+ teams: gold top 10 + 0.2%, silver top 5%, bronze top
10%, percentages rounded down; see Kaggle's published
[medal rules](https://www.kaggle.com/progression/competitions); on 1,437 teams):

| | Rank needed | Score at that rank now |
| --- | --- | --- |
| Gold | ≤ 12 | 23.67 |
| Silver | ≤ 71 | 22.68142 |
| Bronze | ≤ 143 | 22.67253 |

Two public blends dominate the medal zone. Today's best (22.70030 on the board) is already held by 10 teams, inside
silver. The earlier one (22.67253) is held by 58 teams and sits exactly on the bronze line. The competition rules
break ties by whichever submission was entered first
([rules](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/rules)). Both board scores
match our local figures for those circuits, as section 4 would predict. These are snapshots: the field is still growing, and final medals come
from the private leaderboard.

The public plateau moved from about 21.5 to about 22.7 in a day, and 370 teams now sit at or above 22.5. The top two
moved further, to 25.5 and 26.1, with nothing public near them.

## 2. Public frontier by task (re-verified by us unless stated)

| Task | Best public | Qubits | Cost | Ratio to reference | Source |
| --- | --- | --- | --- | --- | --- |
| invsqrt | mrmeows's 63-qubit circuit with busyaprime's re-synthesised lookup (1,334.5 Toffolis plus rotations = 1,629.6 Toffoli-equivalents per shot) | 63 | 102,664 | 20.23 | [Akmal Xodarev (busyaprime), *LB 22.33981 · 22x Cheaper Than Google's Circuits*](https://www.kaggle.com/code/busyaprime/lb-22-33981-22x-cheaper-than-google-s-circuits); core design [mrmeows, *21.0 circuits, let's all reach the plateau!*](https://www.kaggle.com/code/mrmeows/21-0-circuits-let-s-all-reach-the-plateau) |
| DQI | heliosli's 889-qubit narrowing of SomeJoint's 897-qubit fused decoder | 889 | 4.1478 × 10¹¹ | 58.92 | [heliosli, *The dataset is now public now. Enjoy!*](https://www.kaggle.com/code/heliosli/the-dataset-is-now-public-now-enjoy); core design [Agniv Banerjee (somejoint), *22.34 Local GM: 897-Qubit DQI*](https://www.kaggle.com/code/somejoint/quantum-circuits-22-34-local-gm-897-qubit-dqi) |
| Chemistry | somejoint's R8/B143/C60 factorization (8-bit words, full-XOR five-word loaders) with ten of fifteen angle circuits rebuilt by heliosli | 430 | 3.4191 × 10⁸ (16,306 Toffolis, λ_eff 48.764) | 9.81 | [heliosli, *The dataset is now public now. Enjoy!* (version 12)](https://www.kaggle.com/code/heliosli/the-dataset-is-now-public-now-enjoy); core design [Agniv Banerjee (somejoint), *22.67 Local GM: Full-XOR Chem*](https://www.kaggle.com/code/somejoint/quantum-circuits-22-67-local-gm-full-xor-chem) |

Best-per-task geometric mean: **22.70031** (re-verified by us; heliosli's notebook reports 22.7003). All three
notebooks are Apache-2.0.

## 3. What's new

- **Chemistry at rank 60 with full-XOR loaders.** somejoint's team pruned their R8 factorization to B = 143 vectors
  at rank C = 60 with 8-bit words, then refitted weights and offsets for the lowest λ_eff inside the error and gap
  limits. Every loader is now a full-XOR five-word lookup, and each low product is computed once. Within hours there were
  four variants on the same factorization, all re-verified by us:
  - [Atakan Aldemir, *C60 Reachable Cache*](https://www.kaggle.com/code/atakanaldemir/quantum-circuits-c60-reachable-cache), 3.4348 × 10⁸;
  - [Kirill Tymchishin, *Exact Gauge + 5-Word QROM*](https://www.kaggle.com/code/kirilltym/quantum-circuits-exact-gauge-5-word-qrom), 3.4481 × 10⁸;
  - [avazaku, *C60 Gauge Descent + Full-XOR QROM*](https://www.kaggle.com/code/avazaku/c60-gauge-descent-full-xor-qrom), 3.4496 × 10⁸;
  - [AlekseiProvorov, *Full-XOR Chemistry: No Borrowed Bits*](https://www.kaggle.com/code/leoprovorov/full-xor-chemistry-no-borrowed-bits), 3.4512 × 10⁸.
- **And then 22.70 (heliosli, version 12, late afternoon).** heliosli kept somejoint's factorization file byte for
  byte and five of its angle circuits, and rebuilt the other ten: 60 fewer Toffolis (16,366 → 16,306), chemistry
  cost −0.37%, best blend 22.70031 (re-verified by us).
- **889-qubit DQI from heliosli**, 1.8% cheaper than somejoint's 897-qubit decoder (re-verified by us). Comparing
  the two files, the change looks like a rearranged schedule for releasing and recomputing scratch registers. It removes
  8 qubits for only 30 more Toffolis, far below the roughly 1,100 Toffolis a qubit is worth at this size (our
  reading of the diff, not the author's description).
- **busyaprime's invsqrt** is mrmeows's 63-qubit circuit with one lookup re-synthesised from 186 to 172 Toffolis.
  The rest of the file is byte-identical to mrmeows's (we checked by diff), and the cost falls 0.85%.
- **Portfolios.** [alex chilton's *Verified Three-Task Portfolio*](https://www.kaggle.com/code/alexchilton/22-57-lb-verified-three-task-portfolio)
  scored 22.57358 on the board with the previous chemistry (observed). Its current version carries
  the 22.67 components. Its version table is the clearest record of who contributed what.
- **Audits and calculators.**
  - [ENES CAN (povaj1), *QCC Score Audit*](https://www.kaggle.com/code/povaj1/qcc-score-audit-local-vs-real-lb-scores)
    pairs six submissions' local reports with the scores they got (section 4 uses it).
  - [Harshit Kumar, *Marginal Value Calculator*](https://www.kaggle.com/code/harshitsama/quantum-circuits-marginal-value-calculator)
    ([753527](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753527)):
    removing one qubit is worth 1.6% of task cost on invsqrt and 0.23% on DQI and chemistry.
  - [Ak, *Reproducibility Audit*](https://www.kaggle.com/code/abisheksrivastav/quantum-circuit-portfolio-reproducibility-audit).
- **From source.** [Dmitrii Gluzdov, *Quantum Circuit Compilation from Source*](https://www.kaggle.com/code/dmitriigluzdov/quantum-circuit-compilation-from-source)
  builds its circuits from generators. That matters given the hosts' answer on prize eligibility in section 4.

## 4. From the forum

- **kickmix is on PyPI** ([753307](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753307),
  Tanuj Khattar, host): `pip install kickmix` now has prebuilt wheels for Windows, Linux and macOS (Apple Silicon
  and Intel), Python 3.10–3.14. This is the easiest way to verify locally.
- **Generators are required for prizes** ([752010](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752010)).
  Tanuj Khattar: "a circuit generator that produces the exact circuits must be provided", for all benchmarks and
  all prizes. Kaggle staff are confirming the details for reused upstream artifacts. The related question about
  reusing Apache-2.0 circuits ([753316](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753316))
  is still open.
- **Provenance of the frontier DQI.** Agniv Banerjee (somejoint) answered a question on
  [their notebook](https://www.kaggle.com/code/somejoint/quantum-circuits-22-67-local-gm-full-xor-chem/comments):
  all the DQI code in it was written by their team and is released under Apache 2.0, and earlier decoders are
  credited as design lineage only, implemented from a written description rather than from code. With
  heliosli's 889-qubit modification also Apache-2.0, today's public frontier is open-licensed throughout.
- **More measured nulls** ([753185](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753185),
  alex chilton). The chemistry Givens tree: re-encoding the same C46 vectors on 18 other trees and rounding to 8 bits,
  the starter tree rounds best. It reaches 219 mHa against 236 for a binary heap; relabelled orbitals are 17–29% worse,
  and chain or star trees are 80% worse. A staged round-and-refit on the best rivals is still running. He also points
  to the Eckart–Young floor we posted, R·C ≥ 290 for the 20 mHa gate, against 460–496 in public shapes (480 for
  today's best chemistry): room on paper, if anything can train there. On DQI, his list says peak-shaving the
  877-qubit build finds nothing under the exchange rate. heliosli's −8 qubits on the 897-qubit decoder suggests that
  build had slack the 877 didn't.
- **"What is the ceiling?"** ([753404](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753404)).
  Guesses range from 25 to 35. One reply notes that the pace has slowed since the first days, another that small
  randomness could decide near-ties. On that point:
- **The DQI score is a fixed draw per circuit, not per submission.** The verifier seeds every shot from
  SHA-256(circuit hash + salt) (`verify_submission_dqi.py`). The leaderboard numbers behave as if its salt is fixed:
  - Every pair of local report and board score we could find fits a whole number of hosted DQI failures to within
    0.00001: povaj1's six submissions, plus alex chilton's and busyaprime's shared 22.33981. busyaprime's notebook
    makes the same inference for its own entry.
  - The same circuit gets the same count every time, including across teams. The 897-qubit decoder drew 2 failures
    (−0.28%) for both teams that submitted it; locally it draws 0.
  - The 889-qubit decoder used by today's best portfolios draws 0 on the board, the best possible R.
  - A check on our own account: we predicted 22.70031 for today's best blend before submitting it, and the board
    returned 22.70030 (it shows five decimals).
  - So resubmitting the same circuit won't change its score. Any byte change re-rolls it, and final scoring uses a
    private salt, so the public draw tells you little about the final one
    ([753424](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753424)).

## 4b. Position risk: the final DQI draw, circuit by circuit

Final standings are rescored with a private salt. For DQI that is one fresh draw of failures for each circuit file,
and it can move a team's score by about 0.15% per failure. Two things follow from the verifier's design:

- **Identical files share a fate.** The seed is the file's hash plus the salt, so every team that submitted the same
  DQI bytes gets the same private draw. Copies of one portfolio rise or fall together; the risk is relative to teams
  with different DQI files.
- **A lucky public draw can only get worse.** A circuit showing 0 failures on the board can't improve at the final
  rescoring; one showing 2 probably will.

Per circuit, from our salt check (9 salts × 102,400 shots each, re-run daily; the leaderboard count is inferred
from board scores as in section 4). The failure rate behind each row is uncertain: for the frontier 889-qubit file
the 95% range gives P(worse) from 0.03 to 0.55.

| DQI file (SHA-256 prefix) | Qubits | Board failures | P(final worse) | P(final better) | Expected change |
| --- | --- | --- | --- | --- | --- |
| `4ee7a17b4cbc` (today's frontier) | 889 | 0 | 0.20 | 0 | −0.03% |
| `ef334fee63f3` | 883 | 0 | 0.28 | 0 | −0.05% |
| `7beea291a5de` | 889 | 0 | 0.63 | 0 | −0.14% |
| `f3d335c3f693` | 877 | 0 | 0.67 | 0 | −0.16% |
| `b3354f60e459` | 897 | 2 | ≈ 0 | 0.96 | +0.23% |

**What one failure is worth in rank.** On the 17:13 board, 97 teams sit within 0.153% (one extra DQI failure) below
the silver line, and 30 below the bronze line; the gold line has only 2 within reach. So around silver and bronze,
the final DQI draw can move a team by dozens of places. Teams whose DQI file differs from their neighbours' carry
that risk independently.

Our estimate, local re-verification; invsqrt and chemistry scores barely move with the salt. Tools: the
[DQI Salt Check](https://www.kaggle.com/code/manwithacat/dqi-salt-check-public-frontier) notebook measures every DQI
file in our dataset daily, so you can look up your own file by hash.

## 5. Tools and notes from us

- **[VQCC Apache-2.0 public circuits](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits)**
  (dataset, Apache-2.0). It holds 60 DQI, invsqrt and chemistry circuits from public Apache-2.0 notebooks, with
  chemistry factorizations, de-duplicated by hash. Each entry has a credit and licence check and our local verifier
  score. Credit has two levels: the notebook that published the exact file, and, where we've compared files, who
  designed the circuit's core and what the publishing notebook changed. We add new circuits as they appear.
- **[The Real $150K Was the Friends We Made](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made).**
  It picks the cheapest verified circuit per task from that dataset, checks hashes, verifies, and writes a
  submission zip. It reruns whenever the dataset updates.
- **[DQI Salt Check](https://www.kaggle.com/code/manwithacat/dqi-salt-check-public-frontier)** re-measures the
  public DQI decoders under new salts each day.
- **ANGLE_k for every measurement outcome.** A formal check is
  [here](https://www.kaggle.com/code/manwithacat/prove-angle-k-for-every-measurement-outcome), with a technical
  note at [doi:10.5281/zenodo.23721373](https://doi.org/10.5281/zenodo.23721373).
- **Our kickmix fork: work in progress, more coming.**
  ([GitHub](https://github.com/manwithacat/kickmix-macos), Apache-2.0;
  [753230](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753230)).
  Faster and leaner local verification has been our side project this week. It isn't finished, and new features
  are on the way. Where it stands:
  - **Wins.** Results are identical to the official build on all three official verifiers. Large DQI circuits
    verify in a fraction of the memory. Today's release adds a batched phase read-out, which on a Mac makes a full
    invsqrt verification about 3.4× faster than the official build in our tests (the driver that uses it isn't
    published yet). Releases now mirror to Kaggle automatically as the
    [kickmix-fork-wheels](https://www.kaggle.com/datasets/manwithacat/kickmix-fork-wheels) dataset, with a
    [kickmix_tools](https://www.kaggle.com/code/manwithacat/kickmix-tools) script to install and run it offline.
  - **Problems.** Speed doesn't carry over cleanly to Kaggle's Intel CPUs: in our invsqrt timings there, one circuit
    verified 1.6× faster and another 19% slower. The memory saving is the dependable part. With the official
    wheels now on PyPI, use those unless memory is your problem.
  - **Upstream.** One change is with the Qualtran maintainers: streaming `reaction_depth` instead of materialising
    every op ([PR #1994](https://github.com/quantumlib/Qualtran/pull/1994), for
    [issue #1993](https://github.com/quantumlib/Qualtran/issues/1993)). We'll offer more upstream if the maintainers
    want them, within their policy on AI-assisted contributions. All of this code was written with AI assistance.
- **Chemistry cold start on a T4 GPU: negative**, added to the null list in
  [753185](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753185).

- **How this gets made.** Every circuit in today's frontier table was pulled, hashed, licence-checked and re-verified
  by RezOps between 14:30 and 16:30 UTC. When a new one beat the old, the dataset was updated and the blend notebook
  re-ran on it. We're automating more of that loop. If you spot a credit or number that's wrong, say so
  in the comments: corrections go straight into the dataset.

## 6. What we think it means (opinion)

Today's public gain came almost entirely from chemistry, and on one shared factorization. DQI and invsqrt barely
moved. The top two teams are 13–15% above the public 22.70, a gap a single task can't easily close. With
generators now required for prizes, we'd expect more of the public work to become source rather than circuits. And
because the board's DQI draw is fixed per circuit, small public differences between near-identical portfolios
deserve less weight than they get.
