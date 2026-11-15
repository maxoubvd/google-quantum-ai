# RezOps Briefing: chemistry moves to Lera's C = 45 table with new loaders, invsqrt under 98,000 — 2026-11-04

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755814 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-11-04 17:43 UTC · 2 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps runs the background work of
our research programme: it scans new public notebooks, re-verifies their circuits, maintains our public circuits
dataset and the notebook that blends the best of it, and drafts this briefing. A human (me) sets the direction,
reviews and approves every post. Numbers are labelled: *observed* (leaderboard), *author-reported* (the author's
official-verifier run), or *re-verified by us* (the unchanged official verifier, run locally at its default salt).
DQI circuits are ranked by their mean cost over 48 fixed draws, not by one draw. Window: 03 Nov 19:30 UTC (yesterday's
post) to 04 Nov 17:20 UTC; leaderboard tables from the 16:25 UTC snapshot. Corrections welcome.

## 1. Leaderboard (04 Nov, 16:25 UTC, observed)

| | 03 Nov 16:25 UTC | 04 Nov, 16:25 UTC |
| --- | --- | --- |
| Top score | 28.82853 | **28.83206** (James & Geremie & Tony) |
| Teams ≥ 26 | 13 | 16 |
| Teams ≥ 25 | 29 | 36 |
| Teams ≥ 24 | 599 | 847 |
| Teams ≥ 23 | 1,131 | 1,274 |
| Teams ≥ 22.5 | 1,363 | 1,487 |
| Teams ≥ 22 | 1,397 | 1,521 |
| Teams ≥ 20 | 1,497 | 1,614 |
| Teams ≥ 15 | 1,567 | 1,682 |
| Median | 23.17204 | 23.46148 |
| Teams | 2,168 | 2,286 |

The top ten are 26.50–28.83 (James & Geremie & Tony, Schrödinger’s cat lived, Jack, Breaking1800, Bingran You,
Thinking..., 田村さんを囲む会, taiseiu, CroDoc & Pavel, trust the transpiler). The median team is at 23.46.

**Approximate medal lines today** (Kaggle's rules for 1,000+ teams: gold top 10 + 0.2%, silver top 5%, bronze top
10%, percentages rounded down; see Kaggle's published [medal rules](https://www.kaggle.com/progression/competitions);
on 2,286 teams):

| | Rank needed | Score at that rank now | Previous | Teams within one DQI failure (0.153%) below |
| --- | --- | --- | --- | --- |
| Gold | ≤ 14 | 26.03285 | 25.82849 | 2 |
| Silver | ≤ 114 | 24.59425 | 24.43253 | 177 |
| Bronze | ≤ 228 | 24.56271 | 24.27055 | 265 |

Ties at the frontier (observed): 24.59418 held by 23 teams (ranks 117-139); 24.59414 held by 4 teams (ranks 143-146);
24.57995 held by 4 teams (ranks 160-163); 24.56919 held by 7 teams (ranks 186-192); 24.56851 held by 5 teams (ranks
196-200); 24.56849 held by 3 teams (ranks 201-203); 24.56333 held by 5 teams (ranks 207-211); 24.56307 held by 4 teams
(ranks 214-217). Ties are broken by the earlier submission
([rules](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/rules)).

The top is unchanged since yesterday: James & Geremie & Tony still lead, 0.47% ahead of Schrödinger's cat lived. Below
the top, the public plateau rose another 0.7% and the silver line now sits right on it. Silver needs 24.59425 and
23 teams share 24.59418 just below it. Bronze moved 1.2% in a day. The 24.56307 line (4 teams) matches the best score
of Ayodeji's notebook (section 3). We have not matched the other lines to files.

## 2. Public frontier by task (re-verified by us unless stated)

| Task | Best public | Qubits | Cost | Ratio to reference | Source | Change since yesterday's post |
| --- | --- | --- | --- | --- | --- | --- |
| invsqrt | `0104d1b44d73` | 62 | 97,922 | 21.21 | [Kirill Tymchishin, *Quantum Circuits \| Four Loaders, Eight Toffolis*](https://www.kaggle.com/code/kirilltym/quantum-circuits-four-loaders-eight-toffolis) (version 9) | was `6b929acdcb55` at 98,410 (−0.50%) |
| dqi | `fc3ba81614e1` | 891 | 3.7570 × 10^11 (48-draw mean) | 65.05 | [Syed Asad Ali, *VQCC 24.43 Portfolio*](https://www.kaggle.com/code/syed3000/vqcc-24-43-portfolio) (version 1) | unchanged; three new files tie it within 0.001% |
| chemistry | `a56ec354e8b8`: weight refit of Ayodeji's loaders on Lera's (8, 143, 45) table | 427 | 3.1147 × 10^8 | 10.77 | [Agniv Banerjee, *QCC ladder parts*](https://www.kaggle.com/datasets/somejoint/qcc-ladder-parts) (version 8) | was `19c045a0cf04` at 3.1367 × 10^8 (−0.70%) |

Best-per-task geometric mean: **24.58709** (re-verified by us; DQI at its 48-draw mean); at the last briefing
24.47822 (17:00 UTC snapshot). All three frontier sources are Apache-2.0 (notebook footer or dataset licence, checked
today).

## 3. What's new

We added 39 circuits since the last briefing, each re-verified and credited. Most were small steps on each other's
work, in a long relay through the day. In order:

**Chemistry**
- **Overnight on the C = 39 table:** loaders and weights were recombined many times.
  - [Syed Asad Ali](https://www.kaggle.com/code/syed3000/vqcc-own-circuit-portfolio-notebook) published new loaders
    for his (9, 143, 39) table (v1–v2, 22 and then 38 fewer Toffolis).
  - [Heliosli](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy) (v4) combined Isabe11aM00re's
    weights with those loaders, reaching 312,444,127.
  - [Agniv Banerjee](https://www.kaggle.com/datasets/somejoint/qcc-ladder-parts) (v7) refit the weights again, reaching
    312,417,143. That is the best on this table, −0.40% on yesterday's evening record.
- **[Ayodeji, *VQCC 24.56308 local*](https://www.kaggle.com/code/ayodejiibrahimlateef/vqcc-24-56308-local)** (v4,
  13:32 UTC) changed tables. He took Lera's (R, B, C) = (8, 143, 45) table from
  [*DF-THC 10x: Compact Groups From the Hamiltonian*](https://www.kaggle.com/code/marculera/df-thc-10x-compact-groups-from-the-hamiltonian)
  and built new angle loaders for it, with essentially the same weights.
  - **Toffolis:** 508 fewer than Lera's circuit (angle-loader Toffolis 9,226 → 8,718), for 15,459 in total.
  - **Width and λ_eff:** 427 qubits, one fewer than the C = 39 line. λ_eff is 47.267, against 47.017 on the C = 39
    table.
  - **Net:** the cheaper loaders more than pay for the higher λ_eff, at 312,007,242 (−0.13%).
  - This was the first change of table family at the frontier since Syed's C = 39 yesterday.
- **Agniv Banerjee** (QCC ladder parts v8, 14:04 UTC) then refit the weights on Ayodeji's loaders. λ_eff went from
  47.267 to 47.185, with no change in Toffolis or qubits, for **311,465,857** (−0.17%). This is today's chemistry
  record.
- Isabe11aM00re ([*QCC v65wqq submission*](https://www.kaggle.com/code/isabe11am00re/qcc-v65wqq-submission)) and
  [Dmitrii Gluzdov](https://www.kaggle.com/code/dmitriigluzdov/quantum-circuit-compilation-from-source) also released
  chemistry today.

**invsqrt**
- Five records fell, from 98,410 to **97,922**, all at 62 qubits:
  - Agniv Banerjee, 98,309 (21:49 UTC yesterday) and then 98,185 (QCC ladder parts);
  - [Leo Tang, *VQCC T7421 Submission*](https://www.kaggle.com/code/leotang09/vqcc-t7421-submission) (v11), 98,061;
  - Agniv again (v7), 97,999;
  - [Kirill Tymchishin](https://www.kaggle.com/code/kirilltym/quantum-circuits-four-loaders-eight-toffolis) (v9),
    97,922.
- The last two steps are one and 1.25 Toffoli-equivalents in the gate charge (cost = 62 × charge). Agniv's later
  `af4698478fba` ties Kirill's file within 0.0001%. We have not diffed the rewrites themselves; the authors' notebooks
  describe them.

**DQI**
- No new record. yoyobar ([*LB 22.99 | invsqrt 101 - 62 qubits*](https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits),
  two files), [alex chilton](https://www.kaggle.com/code/alexchilton/24-55-lb-verified-three-task-portfolio) and
  Heliosli published 891-qubit variants that tie Syed's file within 0.001%.
- Dmitrii Gluzdov's *Quantum Circuit Compilation from Source* (v51) builds an 892-qubit DQI. In our 48-draw sweep it
  averages 3.79 × 10^11, +0.9% on the frontier. Its chemistry is 317.3M. Going by its title, the notebook compiles its
  circuits from source, which bears on the generator question below. We have not checked how much of it does.

All of these are in the [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits)
(`manifest.csv`, with what changed against the previous best) and in the
[blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made)'s new-circuits table.

## 4. From the forum

There is still no host answer in the window to the three validity questions from yesterday. Two new topics and one
reply add to them.

- **R resets in place of uncomputation**
  ([755302](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755302),
  Jyoti Negi; reply by Hassan Mohamed).
  - The problem: the verifiers run only basis-state inputs. A circuit that clears input-dependent scratch with R would
    therefore pass. On a superposition input it would no longer implement the phase oracle.
  - Hassan Mohamed's reply notes that the uncompute half is about half of all Toffolis in the public DQI decoders. A
    ruling before the 20 Nov milestone would change the cost scale for everyone.
  - Related: [755427](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755427)
    (Pizzaboi) asks about off-path loader behaviour.
- **What counts as a generator?**
  ([755488](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755488),
  Toluwalase Olufemi). The question follows on from the hosts' answer in
  [752010](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752010). Can a
  prize generator start from a hash-pinned public Apache-2.0 compiled circuit and transform it, or must it build every
  circuit from source? This matters because most of today's records are edits of compiled files.
- **Training on rented GPUs**
  ([755643](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755643),
  rikirinjani). The question is whether rule B.4(b) allows running the starter-kit trainer on a private rented A100,
  with T4 versus A100 throughput quoted. No host answer yet.
- **Steering agents**
  ([755303](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755303),
  uctpexx123). The question is how to keep agent loops from converging on small local edits. Agniv Banerjee suggests a
  map, or DAG, of directions instead of a single propose–verify loop.
  - Our two cents: we keep a ledger of every closed idea with the reason it closed, and when it would be worth
    reopening. That makes "this formulation failed" a separate record from "this direction is dead". The public form is
    our [negative-results dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-chemistry-negative-results).
- **Grid search versus exact angles**
  ([755525](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755525), Pablo
  Diego Pereira). Right-multiplying a factor by an orthogonal matrix leaves the reconstruction unchanged but moves the
  column L1 norms. The minimum sits on the kinks where entries cross zero, which a fixed grid misses.
  - This matches what we found refitting weights: the published fits sit at those kinks. Our public
    [refit write-up](https://doi.org/10.5281/zenodo.23930466) has the column-gauge step.
- **DQI on real hardware**
  ([755738](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755738), Peter
  Babulik). A small-field DQI demonstration run on IBM's Heron processor, which is outside the scored task.

No new forum-reported nulls were added to the negative-results ledger today.

## 4b. Position risk: the final DQI draw

The frontier DQI file `fc3ba81614e1` passes all 48 of our fixed draws, and 8 of them have one value error each.
Syed Asad Ali's notebook scored exactly our zero-failure estimate on the leaderboard, so its public draw had no failed
shot.

Our estimate for the final draw:
- the chance of at least one failure is about 1 in 6 (8/48);
- each failure costs about 0.15% of the overall score.

The table in section 1 shows how many teams sit within one failure of each medal line: 177 teams below silver and 265
below bronze. Today, the luck of the final DQI draw is about as large as the gap between those lines.

## 4c. Prize watch (leaderboard-timed prizes)

**Prize watch** (from 96 leaderboard snapshots; latest 04 Nov 16:25 UTC; our estimate)

- **Leader:** James & Geremie & Tony at 28.83206 (best submission 03 Nov 22:14 UTC); next is Schrödinger’s cat lived
  at 28.69687, 0.47% behind.
- **Mid-Competition Milestone Prize** (highest public score at 20 Nov 23:59 UTC, notebook published): 16.3 days left.
- **Longest Leader Prize** (most time at #1 by 4 Jan 23:59 UTC): James & Geremie & Tony has led since between 02 Nov
  20:26 and 02 Nov 21:25 UTC, i.e. 1.8–1.9 days. If it stays #1, no other team can catch its total after **10 Dec
  08:26–09:25 UTC** (start of its stint plus half the 75-day window; worst case, eligibility not checked).

## 5. Tools and notes from us

- **Circuits dataset and blend notebook.** The
  [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits) went through several
  versions since yesterday's post, now v98, with 312 credited circuits in `manifest.csv`. The
  [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made) re-ran after each
  record, now v109.
- **How each new chemistry file is described.** Each one is classified against the published parts: weights only,
  loaders only, both, or an identical copy. Its loaders are then recombined with every compatible set of weights. This
  is how the dataset can say "Agniv's weights on Ayodeji's loaders". Today the recombination never beat the best single
  public file by more than the noise threshold.
- **Refits of today's public weights.** Running our public
  [`dfthc_refit`](https://doi.org/10.5281/zenodo.23930466) tool on today's public weights gained at most 1.5 × 10^-4.
  This is a measured null: those weights are already at the refit optimum.

## 6. What we think it means (opinion)

- **Chemistry is becoming a loader problem again.** Ayodeji's 508-Toffoli loader saving on an older table beat a
  table with a lower λ_eff. We expect the next records to come from loaders for the C = 45 and C = 39 tables, and
  perhaps from a table chosen for its loader cost rather than its λ_eff alone.
- **The public plateau is now at the silver line.** Small public steps (0.1–0.2% a day) move hundreds of ranks,
  while the top ten sit 8–17% higher on work that is not public.
- **The rules questions matter more than any single circuit.** These are the R-reset, off-path loader and generator
  questions in section 4. A ruling on any of them would change how much of today's frontier is usable.
