# RezOps Briefing: every task moved — lattice-rounded C = 39 chemistry, an 891-qubit DQI, a refit tool — 2026-11-03

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755495 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-11-03 19:38 UTC · 0 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps runs the background work of
our research programme: it scans new public notebooks, re-verifies their circuits, maintains our public circuits
dataset and the notebook that blends the best of it, and drafts this briefing. A human (me) sets the direction,
reviews and approves every post. Numbers are labelled: *observed* (leaderboard), *author-reported* (the author's
official-verifier run), or *re-verified by us* (the unchanged official verifier, run locally at its default salt).
DQI circuits are ranked by their mean cost over 48 fixed draws, not by one draw. Window: 02 Nov 17:00 UTC to 03 Nov 19:30 UTC (leaderboard tables from the 16:25 UTC snapshot). Corrections welcome.

## 1. Leaderboard (03 Nov, 16:25 UTC, observed)

| | 02 Nov 16:24 UTC | 03 Nov, 16:25 UTC |
| --- | --- | --- |
| Top score | 27.99617 | **28.82853** (James & Geremie & Tony) |
| Teams ≥ 26 | 10 | 13 |
| Teams ≥ 25 | 24 | 29 |
| Teams ≥ 24 | 169 | 599 |
| Teams ≥ 23 | 978 | 1,131 |
| Teams ≥ 22.5 | 1,241 | 1,363 |
| Teams ≥ 22 | 1,279 | 1,397 |
| Teams ≥ 20 | 1,379 | 1,497 |
| Teams ≥ 15 | 1,445 | 1,567 |
| Median | 22.99251 | 23.17204 |
| Teams | 2,043 | 2,168 |

The top ten are 26.35–28.83 (James & Geremie & Tony, Schrödinger’s cat lived, Breaking1800, Jack, Thinking..., Bingran You, 田村さんを囲む会, trust the transpiler, taiseiu, Adisak Sukul). The median team is at 23.17.

**Approximate medal lines today** (Kaggle's rules for 1,000+ teams: gold top 10 + 0.2%, silver top 5%, bronze top
10%, percentages rounded down; see Kaggle's published [medal rules](https://www.kaggle.com/progression/competitions);
on 2,168 teams):

| | Rank needed | Score at that rank now | Previous | Teams within one DQI failure (0.153%) below |
| --- | --- | --- | --- | --- |
| Gold | ≤ 14 | 25.82849 | 25.71721 | 1 |
| Silver | ≤ 108 | 24.43253 | 24.14832 | 59 |
| Bronze | ≤ 216 | 24.27055 | 23.94701 | 22 |

Ties at the frontier (observed): 24.45812 held by 9 teams (ranks 73-81); 24.43253 held by 32 teams (ranks 101-132); 24.43209 held by 53 teams (ranks 135-187); 24.18056 held by 6 teams (ranks 299-304); 24.18054 held by 4 teams (ranks 305-308); 24.18048 held by 6 teams (ranks 309-314); 24.17647 held by 29 teams (ranks 326-354); 24.17645 held by 21 teams (ranks 355-375).
Ties are broken by the earlier submission ([rules](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/rules)).

James & Geremie & Tony took the top last night (28.83, first seen at #1 between 20:26 and 21:25 UTC on 02 Nov) and
lead Schrödinger's cat lived by 2.7%. The public plateau jumped about 2% today, because every task's best public
circuit moved (section 2). The most-held lines near it match public notebooks: 24.43209 (53 teams) is Syed Asad Ali's
*VQCC 24.43 Portfolio* exactly; 24.43253 (32 teams) and 24.45812 (9 teams) sit just above it, probably the same
portfolio with one of the day's later chemistry or invsqrt improvements (we have not matched them to a file). Teams at 24 or more went from 169 to 599 in a day.

## 2. Public frontier by task (re-verified by us unless stated)

| Task | Best public | Qubits | Cost | Ratio to reference | Source | Change since last briefing |
| --- | --- | --- | --- | --- | --- | --- |
| invsqrt | `6b929acdcb55` | 62 | 98,410 | 21.10 | [Syed Asad Ali, *VQCC 24.43 Portfolio*](https://www.kaggle.com/code/syed3000/vqcc-24-43-portfolio) (version 3), after Kirill Tymchishin's two-product windows | was `bce2948ff24e` at 100,568 (−2.15%) |
| dqi | `fc3ba81614e1`: flags carried by measurement, measurement-based erasure, product caching on idle wires | 891 | 3.7570 × 10^11 (48-draw mean) | 65.05 | [Syed Asad Ali, *VQCC 24.43 Portfolio*](https://www.kaggle.com/code/syed3000/vqcc-24-43-portfolio) (version 1) | was `e3928244a1c1` at 3.7763 × 10^11 (−0.51%) |
| chemistry | `19c045a0cf04`: weights-only refit of Agniv Banerjee's loaders for Syed Asad Ali's C = 39 table | 428 | 3.1367 × 10^8 | 10.70 | [Isabe11aM00re, *QCC v58x submission*](https://www.kaggle.com/code/isabe11am00re/qcc-v58x-submission) (version 4) | was `5b4482b2e74c` at 3.2679 × 10^8 (−4.02%) |

Best-per-task geometric mean: **24.48873** (re-verified by us, 18:20 UTC); at the last briefing 23.94157.
All three frontier sources are Apache-2.0 (notebook footer or dataset licence, checked today). Our re-verified costs
agree with the authors' own figures: Syed Asad Ali reports 98,843.84, 891 qubits and 314,740,179 for his three
circuits, and his notebook's observed public score of 24.43209 matches our zero-failure local estimate (24.43210), so
its leaderboard DQI draw had no failed shot.

## 3. What's new

**For a few hours this afternoon Syed Asad Ali held all three task records at once**, with his first public notebook.
Since then the frontier has stacked several improvements on top of his work: Kirill Tymchishin's invsqrt rewrite,
Heliosli's weight refit on his chemistry table, and Agniv Banerjee's new loaders for that table. DQI is still his.
In order (credited, re-verified by us):

- **Syed Asad Ali, [*VQCC 24.43 Portfolio*](https://www.kaggle.com/code/syed3000/vqcc-24-43-portfolio)** (new
  competitor, v1) moved all three tasks at once. Chemistry: a new (R, B, C) = (9, 143, 39) table, one column fewer than
  every earlier public table, rounded with Babai's nearest-plane lattice method and refit with an augmented
  Lagrangian; loaders built by exact linear-matroid-parity product selection (296 ANDs per coded batch), −0.84%.
  DQI: pole flags carried by measurement (an idea he credits to Ali Haydar Özdağ, via SomeJoint's 892-qubit line),
  measurement-based erasure and exact product caching on idle wires, 902 → 891 qubits, −0.51% by the 48-draw mean.
  The notebook says one rare high-weight error class is handled approximately; in our 48 draws all pass and 8 have one
  value error, priced in the mean. invsqrt: adder carries erased by measurement plus a conditional CZ instead of
  reverse ripples, −1.7%. The notebook also lists what did not work (7-bit angles, other loader splits, the DQI
  key-equation route).
- **Kirill Tymchishin, [*Quantum Circuits | Four Loaders, Eight Toffolis*](https://www.kaggle.com/code/kirilltym/quantum-circuits-four-loaders-eight-toffolis)**
  (v2) rewrote five three-product windows of that invsqrt as two-product windows with an XOR factoring that reuses a
  product already computed for a later target: 98,534 (−0.31%), with an exhaustive check of the identity.
- **Heliosli, [*The dataset is public now. Enjoy!*](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy)**
  refit weights on each new table within hours, as on the previous days: on Syed's table, −0.005%.
- **Agniv Banerjee (SomeJoint), [*QCC ladder parts*](https://www.kaggle.com/datasets/somejoint/qcc-ladder-parts)**
  (v4) rebuilt the loaders of Syed's table (44 fewer Toffolis) with new weights: 313,678,954, −0.33%, today's chemistry
  record. The same dataset carries a 98,743 invsqrt, also from the notebook
  [*Quantum Circuits 24.47 local GM, split carries*](https://www.kaggle.com/code/somejoint/quantum-circuits-24-47-local-gm-split-carries).
- **This evening** two more records landed: Syed Asad Ali's notebook v3 brought invsqrt to **98,410** (−0.13% on
  Kirill's), and [Isabe11aM00re](https://www.kaggle.com/code/isabe11am00re/qcc-v58x-submission) (v4) refit the weights
  of Agniv's circuit to 313,669,863 (−0.003%). Kirill Tymchishin's v3 (98,433) and Syed's v3 chemistry update were
  also added.

All 40 circuits added since the last briefing, each re-verified and credited, are in the
[circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits) (`manifest.csv`) and in the
[blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made)'s new-circuits table.

## 4. From the forum

Three questions to the hosts today matter to everyone, because they ask what the final evaluation will check beyond
the published verifier. No host answer yet in the window.

- **What does Task III's loader have to do off the verified path?**
  [755427](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755427)
  (Pizzaboi): the README describes ANGLE_k as an XOR oracle on any target register, but the verifier only tests the word
  states its traversal reaches. A minimal reproducer adds three CX gates to one loader: it passes the full verifier at
  the reference cost but differs from arbitrary-target XOR on 7,952 of 15,904 extra inputs. Are invariants of the
  reachable states allowed?
- **Is R (reset) allowed in place of uncomputation?**
  [755302](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755302)
  (Jyoti Negi): the verifiers run basis-state inputs only, so clearing input-dependent scratch with R passes, yet on
  superposition inputs it would measure the scratch and break the phase oracle.
- **What must winners deliver for parts they did not write?**
  [752010](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752010):
  the host (Tanuj Khattar) wrote a week ago that prize eligibility needs a generator that produces the exact circuits,
  for all benchmarks and prizes; Houssem Alayet asked today how that applies to an Apache-2.0 compiled input with no
  generator.

Also in the window: alex chilton's list of about 30 measured nulls across all three tasks
([753185](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753185)) drew
replies; snuffkin posted a free introductory guide to FTQC circuit optimisation for this competition
([754970](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754970)); and
there are new threads on steering agents
([755303](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755303)) and on
tools newcomers use
([755404](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755404)).

## 4b. Position risk: the final DQI draw, circuit by circuit

| DQI file (SHA-256 prefix) | Qubits | Board failures | P(final worse) | P(final better) | Expected change | One more failure |
| --- | --- | --- | --- | --- | --- | --- |
| `4ee7a17b4cbc` | 889 | 0 | 0.20 | 0.00 | -0.033% | -0.153% |
| `b3354f60e459` | 897 | 2 | 0.00 | 0.96 | +0.233% | -0.116% |
| `7beea291a5de` | 889 | 0 | 0.63 | 0.00 | -0.142% | -0.153% |
| `f3d335c3f693` | 877 | 0 | 0.67 | 0.00 | -0.158% | -0.153% |
| `ef334fee63f3` | 883 | 0 | 0.28 | 0.00 | -0.050% | -0.153% |
| `0161fd2a6e50` | 889 | 1 | 0.31 | 0.33 | -0.005% | -0.129% |

Our estimate, from the hosted salt-check notebook's draws. A DQI failed shot raises that task's cost through the
repetition penalty, which moves a portfolio score by about 0.15% (the "one more failure" column, and the
within-one-failure counts in section 1). Today's frontier file `fc3ba81614e1` is not in the table yet: its
leaderboard draw had no failure (the observed 24.43209 matches our zero-failure estimate), and in our 48 local draws
8 have one value error, so its final draw is more likely to be worse than its public score than better.

## 4c. Prize watch (leaderboard-timed prizes)

**Prize watch** (from 79 leaderboard snapshots; latest 03 Nov 16:25 UTC; our estimate)

- **Leader:** James & Geremie & Tony at 28.82853 (best submission 03 Nov 01:42 UTC); next is Schrödinger’s cat lived at 28.07608, 2.68% behind.
- **Mid-Competition Milestone Prize** (highest public score at 20 Nov 23:59 UTC, notebook published): 17.3 days left.
- **Longest Leader Prize** (most time at #1 by 4 Jan 23:59 UTC): James & Geremie & Tony has led since between 02 Nov 20:26 and 02 Nov 21:25 UTC, i.e. 0.8–0.9 days. If it stays #1, no other team can catch its total after **10 Dec 08:26–10 Dec 09:25 UTC** (start of its stint plus half the 75-day window; worst case, eligibility not checked).

## 5. Tools and notes from us

- The [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits) added 40 circuits
  since the last briefing (versions 79–87), each re-verified and credited.
- **Our kickmix-engine is now on Kaggle as wheels**
  ([kickmix-engine-wheels](https://www.kaggle.com/datasets/manwithacat/kickmix-engine-wheels), version 0.6.0), with
  the [`kickmix_tools`](https://www.kaggle.com/code/manwithacat/kickmix-tools) utility script able to verify all three
  tasks with it. It runs the competition's own verifier code and reproduces the official verifier's results field for
  field (checked on passing and failing controls). In the
  [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made), verifying all three
  current picks on Kaggle took 90 s for invsqrt and 33 s for DQI (the unchanged verifier takes roughly 14 and 4
  minutes there), and the whole notebook 3.4 minutes; `FORCE_FULL_VERIFY = True` runs the unchanged official verifier
  instead. Anyone can attach the same datasets.
- The blend notebook again shows, for each task, the top five public circuits under the headline, as well as a table
  of every circuit added on the dataset's latest day.
- **New today: a chemistry weight-refit tool.** [`dfthc_refit`](https://www.kaggle.com/code/manwithacat/dfthc-refit) is a
  Kaggle utility script (Apache-2.0) that takes any chemistry submission folder (the 15 angle loaders and
  `submission.npz`), keeps the words and loaders fixed (so qubits and Toffolis do not change), and lowers λ_eff within
  the error_1, error_2 and E_gap limits. It alternates a constrained refit with the exact column-gauge step and stops
  on its own; every number comes from the official verifier code, and the output is never worse than the input. The
  [demo notebook](https://www.kaggle.com/code/manwithacat/df-thc-weight-refit-demo) runs it on the current best
  public chemistry family and verifies before and after with kickmix-engine; on Kaggle CPUs it took about 34 minutes
  and lowered λ_eff on `4a9adf368a10` by 4.3e-5 relative, with the after-check agreeing with the official `verify()`. Two things we learned building it,
  useful to anyone refitting: published weights carry thousands of tiny deliberate entries (zeroing them costs about
  1e-5 in λ_eff), and refitted points often sit where two eigenvalues of h1′ are almost exactly zero, a kink a plain
  refit cannot cross; the tool holds that block and moves along it.
- **If you work on chemistry, refit before you publish or submit a new circuit.** A new table or new loaders usually
  arrive with weights that are not yet at their best for those loaders, and today every new chemistry record was
  followed within hours by someone else's refit of it. Running the refit first captures that gain in your own entry.
  Example, re-verified by us: on Agniv Banerjee's `4a9adf368a10` the tool reached 313,664,330 (−4.7e-5) in about 9
  minutes locally, released as our Hedge Wizardry HW-15 (hosted entry gate PASS); a separate public refit of the same
  circuit reached −2.9e-5.

## 6. What we think it means (opinion)

**The mathematics has changed.** In late October the public chemistry work was numerical linear algebra: double
factorisation and tensor hypercontraction, eigen-decompositions, gradient refits and Jacobi sweeps (plane rotations
that zero one pair at a time; an exact two-column Jacobi gauge was the clearest public example). DQI started from
textbook algebraic decoding over GF(2¹²), with decoders written as Jacobi continued fractions, so the name turned up
in two tasks for different reasons; invsqrt from polynomial approximation and fixed-point arithmetic. The middle
period was reversible-circuit engineering: measurement-based uncomputation, conditionally clean ancillas, spooky
pebbling, and GF(2) linear algebra for the lookups. Today's work brought methods new to this competition's public
record (not new to mathematics): rounding as a closest-vector problem in a lattice (Babai's nearest-plane method),
augmented-Lagrangian refits that spend the error budget on λ, linear matroid parity to choose shared products, and
measurement-carried flags. Our reading: the skill that pays now is treating each design choice as a discrete
optimisation problem with an exact solver (lattices for rounding, matroids for sharing, pebbling for qubits), and
the fast exact evaluators that make searching those spaces practical.
