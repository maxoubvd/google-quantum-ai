# RezOps Briefing: chemistry becomes a parts market, hengck23's 889-qubit DQI takes the lead, plateau at 22.87 — 2026-10-29

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753923 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-29 19:56 UTC · 2 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps runs the background work of
our research programme: it scans new public notebooks, re-verifies their circuits, maintains our public circuits
dataset and the notebook that blends the best of it, and drafts this briefing. A human (me) sets the direction,
reviews and approves every post. Numbers are labelled: *observed* (leaderboard), *author-reported* (the author's
official-verifier run), or *re-verified by us* (the unchanged official verifier, run locally at its default salt).
Window: 28 Oct 17:16 UTC to 29 Oct 19:45 UTC. Corrections welcome.

## 1. Leaderboard (29 Oct, 19:45 UTC, observed)

| | 28 Oct 17:13 UTC | 29 Oct 19:45 UTC |
| --- | --- | --- |
| Top score | 26.14113 | **26.14113** (Breaking1800) |
| Teams ≥ 25 | 2 | 5 |
| Teams ≥ 24 | 6 | 15 |
| Teams ≥ 23 | 25 | 67 |
| Teams ≥ 22.5 | 370 | 698 |
| Teams ≥ 20 | 734 | 916 |
| Median | 20.51 | 21.88 |
| Teams | 1,437 | 1,588 |

The top ten are 24.44–26.14 (Breaking1800, Schrödinger's cat lived, Jack, Bingran You, Thinking..., Z³uperposition,
Distortion!!, Kevin Bogner, VAGO Superposition, Adisak Sukul).

**Approximate medal lines** (Kaggle's rules for 1,000+ teams: gold top 10 + 0.2%, silver top 5%, bronze top 10%,
rounded down; [medal rules](https://www.kaggle.com/progression/competitions); on 1,588 teams):

| | Rank needed | Score at that rank now | Yesterday |
| --- | --- | --- | --- |
| Gold | ≤ 13 | 24.25946 | 23.66860 |
| Silver | ≤ 79 | 22.97381 | 22.68142 |
| Bronze | ≤ 158 | 22.86842 | 22.67253 |

The public plateau rose from 22.70 to about 22.87 and the field filled in behind it: 698 teams are now at 22.5 or
better, almost double yesterday. The bronze line sits exactly on a public blend: 22.86842 is held by 7 teams, and
22.86836 by 9 more. Yesterday the best public blend reached silver; today silver (22.974) is about 0.5% above the best
public blend, and everything above it is private work. Ties are broken by the earlier submission
([rules](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/rules)). The top score
hasn't moved, but the group above 24 grew from 6 to 15 teams.

## 2. Public frontier by task (re-verified by us)

| Task | Best public | Qubits | Cost | Ratio to reference | Since yesterday |
| --- | --- | --- | --- | --- | --- |
| invsqrt | busyaprime's second re-synthesis on mrmeows's circuit ([Akmal Xodarev (busyaprime), *LB 22.72431*](https://www.kaggle.com/code/busyaprime/lb-22-72431-22x-cheaper-than-google-s-circuits), v10) | 63 | 102,349 | 20.29 | −0.31% |
| DQI | hengck23's 889-qubit decoder, iteration 014 ([hengck23, *DQI: kmb at 59.56854*](https://www.kaggle.com/code/hengck23/dqi-kmb-at-59-56854); file published in [Kirill Tymchishin, *Exact Gauge + 5-Word QROM*](https://www.kaggle.com/code/kirilltym/quantum-circuits-exact-gauge-5-word-qrom), v25) | 889 | 4.1027 × 10¹¹ | 59.57 | −1.09% |
| Chemistry | A recombination of public parts: Atakan Aldemir's word table, loaders by Atakan Aldemir, heliosli and Dmitrii Gluzdov, weights by busyaprime ([alex chilton, *Verified Three-Task Portfolio*](https://www.kaggle.com/code/alexchilton/22-57-lb-verified-three-task-portfolio), v62) | 430 | 3.3906 × 10⁸ (16,266 Toffolis, λ_eff 48.476) | 9.90 | −0.83% |

Best-per-task geometric mean: **22.87030** (re-verified by us); yesterday 22.70031. All the notebooks above are
Apache-2.0. Every task moved, chemistry in 17 small steps, DQI in 5, invsqrt in 1.

## 3. What's new

We added 45 new public circuits to the dataset in the window and re-verified each one. Of the new chemistry builds,
one (`6fbdd8519c33`) failed the official verifier's gap floor (E_gap 3.91481 against the strict 3.914898) and was left
out.

**DQI: five steps, the last one the biggest.**
- [Atakan Aldemir, *Iso-Gap + DQI h5*](https://www.kaggle.com/code/atakanaldemir/quantum-circuits-iso-gap-dqi-h5)
  (889 qubits, −0.19%), then
  [Agniv Banerjee (somejoint), *22.75 Local GM: 897-Qubit DQI*](https://www.kaggle.com/code/somejoint/quantum-circuits-22-75-local-gm-897-qubit-dqi)
  (−0.26%: rewritten convert/unconvert stages at 897 qubits),
  [Mudit Chaturvedi, *895Q Streamed Decoder*](https://www.kaggle.com/code/muchat007/895q-streamed-decoder-measured-node-choices)
  (−0.07%) and
  [Atakan Aldemir, *888Q Repair and Native Chemistry Co-Design*](https://www.kaggle.com/code/atakanaldemir/qcc-888q-repair-and-native-chemistry-co-design)
  (888 qubits, −0.07%).
- Then **hengck23's iteration 014** (−0.51%). It first appeared in a verifier log in
  [*DQI: kmb at 59.56854*](https://www.kaggle.com/code/hengck23/dqi-kmb-at-59-56854) at 12:56 UTC, with the binary in
  hengck23's dataset [hengck23-dqi-59-57](https://www.kaggle.com/datasets/hengck23/hengck23-dqi-59-57). Kirill Tymchishin
  then published the file at 14:08 as the "889Q hengck23" entry of his decoder comparison, which he integrated and
  audited. It is described there as projective extraction and reverse, rolling continuants and tail/head folds. Compared
  with Mudit's 895-qubit file it removes 6 qubits and about 11,500 static Toffolis. Extract and unextract no longer
  use controlled swaps, and unsynthesis is about 13,000 Toffolis cheaper. Credit note: our dataset first credited this
  file to Kirill, the notebook it came from. It now names hengck23 as the designer and Kirill as integrator.
- hengck23 explains the workflow in
  [753607](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753607). The
  decoder is broken into steps, each step profiled for qubits and Toffolis, and an agent refines the steps. A
  "circuit surgery" pass then cuts and pastes better modules from public `.kmb` files, keeping each one only if an A/B
  test says it helps. The generator code comes afterwards, by reverse engineering.

**Chemistry: loaders and weights now improve separately.** A chemistry submission is a set of parts:
- a word table: the angle words the circuits load;
- fifteen loader circuits, one per angle slot, each switching the word registers from one block of words to the next;
- the weights file `submission.npz`, which sets λ_eff.

Today authors improved one kind of part and borrowed the others, so the frontier moved in two separate threads:
- **Fewer Toffolis (loaders).** heliosli rebuilt loaders four times: 16,306 → 16,296 → 16,286 → 16,276 → 16,272
  Toffolis per shot, in [*The dataset is public now*](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy)
  and its predecessor. [Dmitrii Gluzdov](https://www.kaggle.com/code/dmitriigluzdov/quantum-circuit-compilation-from-source)
  then rebuilt ten slots for 16,266, the fewest so far. Gluzdov builds his circuits from source.
- **Lower λ_eff (weights).** Gluzdov 48.632 → 48.614,
  [Ahmed Berat Özer](https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression) 48.590, Gluzdov
  48.573, Atakan Aldemir 48.560, Ahmed Berat Özer **48.476** (−0.17%, the largest single chemistry step today), and
  [busyaprime](https://www.kaggle.com/code/busyaprime/lb-22-86841-22x-cheaper-than-google-s-circuits) 48.476 (−0.001%).
  sigmaborov, Baidalin Adilzhan, Nilanjan_Saha123 and AlekseiProvorov also published weights.
- **A new word table.** [Atakan Aldemir](https://www.kaggle.com/code/atakanaldemir/qcc-888q-repair-and-native-chemistry-co-design)
  shipped the day's first new word table, with all fifteen loaders and the weights redesigned together (16,270
  Toffolis). Six of its fifteen word blocks are new. Its loaders are now part of the frontier.
- **Mixing parts.** Because loaders and weights come from different authors, the cheapest chemistry is often a
  combination that no one has published yet. alex chilton published two such combinations, each the best single
  chemistry release at the time:
  - [`06c6d23ddcc4`](https://www.kaggle.com/code/alexchilton/22-57-lb-verified-three-task-portfolio) at 02:23 UTC;
  - `bcf29775f82b` at 19:10 UTC, today's frontier.

  Both match the assembly our blend notebook builds from the dataset, byte for byte. The dataset credits each part
  to its first public author.
- **A different idea.** [Lera, *Merge Groups Before Dropping Them*](https://www.kaggle.com/code/marculera/df-thc-merge-groups-before-dropping-them)
  reaches λ_eff 48.17, the lowest we've seen, but at 690 qubits and 25,434 Toffolis, so its cost is 8.45 × 10⁸. It's
  not competitive as it stands. It is the only public build today that moves λ_eff by more than a fraction of a
  percent.

**invsqrt.** busyaprime re-synthesised a second region of mrmeows's 63-qubit circuit (window S1, 25 AND gates instead
of 30). Static Toffolis fall from 1,402 to 1,397 and cost by 0.31%. That was the only invsqrt change.

**Portfolios and assemblies.**
- [Nilanjan_Saha123, *Auto-Assemble the Best QCC Portfolio*](https://www.kaggle.com/code/nilanjansaha123/auto-assemble-the-best-qcc-portfolio-22-70)
  assembles a portfolio automatically.
- [Vyas Senthilkumar, *Verified Assembly of Public QCC Circuits*](https://www.kaggle.com/code/vyassenthilkumar/verified-assembly-of-public-qcc-circuits)
  set yesterday evening's chemistry step (−0.15%).
- [sigmaborov, *22.85569 LB: Exact Gauge Search*](https://www.kaggle.com/code/sigmaborov/22-85569-lb-exact-gauge-search)
  publishes its board score in the title.

## 4. From the forum

- **A floor for λ_eff** ([753633](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753633),
  Lynn Langit). Using only the official verifier's definitions, the gap floor forces λ_eff ≥ about 46.8 for any
  factorization. All the remaining slack is the triangle inequality in the two-body term. Today's public λ_eff is 48.476,
  about 3.6% above that floor, so closing the whole gap is worth at most about 1.2% on score. Lynn's conclusion,
  which today's numbers bear out: chemistry gains have to come from T × Q, not λ_eff.
- **Where to look next in chemistry** ([753738](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753738)).
  theredbluepill (rank 11) suggests the DFT word width, 7-bit or even 6-bit words instead of 8. Given the λ_eff
  floor above, narrower words, which cut T, look like the more promising direction. No public build has tried it yet.
- **hengck23's workflow** ([753607](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753607)),
  summarised in section 3. hengck23 also replied on
  [753403](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753403) that our
  invsqrt viewer is really a data-capture tool. The idea: record the register traces, then search them for qubits with
  the same pattern that could share a wire. It's a good suggestion.
- **Longest Leader Prize timing** ([753654](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753654),
  Geremie Yeo). The publication deadline for that prize is 4 January, 16 days before the end. The suggestion is to
  extend it to the competition end. No host reply yet.
- **Search loops for kickmix** ([753741](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753741),
  Raitul Hasan Priyo) asks whether bigger gains come from early ancilla uncomputation or algebraic rewrites. Today's
  DQI evidence: the large steps came from restructuring whole stages (convert/unconvert, projective extraction), not
  from local peepholes.

## 4b. The DQI draw for today's frontier file

The board still behaves as a fixed draw per circuit file. The 22.86842 held by 7 teams matches our local GM at the
default salt to five decimals, including the one value failure hengck23's 889-qubit file draws there. So that file
carries one failure on the public board, worth about 0.13–0.15% of score. Final standings use a private salt, and a
file that draws one failure publicly can move either way. Our daily salt check doesn't include this file yet, so we
can't put odds on it today. Our 28 Oct table for the older files still applies to them.

## 5. Tools and notes from us

- **[VQCC Apache-2.0 public circuits](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits)**:
  - **105 circuits** (60 yesterday), each re-verified, hash-checked and credited.
  - **New `parts/` folder** for chemistry, holding:
    - every loader by slot and every weights file, with who published each first;
    - the word tables;
    - a family check that flags when a release changes the word table or ships weights that only fit its own words;
    - `assembly_best.json`, the cheapest verified combination of public parts.

  The file-level descriptions on the dataset page are lagging behind; `manifest.csv` and `parts/` are current.
- **[The Real $150K Was the Friends We Made](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made)**
  re-ran about a dozen times today, once per dataset update.
  - Each run picks the cheapest verified circuit per task, or the parts assembly for chemistry, and writes a submission
    zip.
  - It now prints a **receipt**: hash, official verifier and match-to-record gates.
  - It draws **frontier charts** showing each task's public frontier over time and where each gain came from (qubits or
    Toffolis).
  - A **chemistry-as-parts chart** shows the best loader per slot over time.
- **[kickmix-engine](https://github.com/manwithacat/kickmix-engine)** v0.5.0
  ([753642](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753642)) now
  verifies chemistry too. On a Kaggle CPU notebook it's about 4.9× faster than the official verifier when warm, at 4
  threads, with 21% less memory. It matches the official verifier field for field, and `--cross-check` runs both and
  diffs them. Wheels are on GitHub and in the
  [kickmix-engine-wheels](https://www.kaggle.com/datasets/manwithacat/kickmix-engine-wheels) dataset.
- **How this gets made.** RezOps scanned public notebooks every few hours through the day. Each new circuit was pulled,
  hashed, licence-checked and re-verified, then the dataset was updated and the blend notebook re-ran. If a credit is
  wrong, say so; today's DQI credit was corrected that way.

## 6. What we think it means (opinion)

The public frontier now improves in small steps from many hands. Chemistry in particular behaves like a market in
parts: one author's loaders and another's weights combine within hours. It is fast, collaborative and honest about
credit, but the steps are small: chemistry gained 0.83% over the day in 17 releases. The biggest single gain came from
a redesigned DQI stage, not from tuning. On Lynn's bound, chemistry's λ_eff has at most about 3.6% left. That leaves
T × Q, and narrower words, as the open direction. The top of the board hasn't moved, and silver now sits above the best
public blend. The plateau is good enough for bronze today; the teams above it are working beyond the public frontier.


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**manwithacat** — **Topic author** · 2026-10-29 20:27 UTC · 1 votes

**Full disclosure:** we have something that extends the public frontier. With the official verifier run locally, our portfolio comes to **22.96**, against **22.87** for the best public blend. It hasn't been scored on the leaderboard yet.

I'm not making it public yet, because I don't understand it. Literally: it looks like voodoo. So there's a short delay while Claude explains linear algebra to me. I hope to publish the tooling I used to find it fairly soon.

In the meantime, I'd encourage you to look at [kickmix-engine](https://github.com/manwithacat/kickmix-engine) ([753642](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753642); wheels in [kickmix-engine-wheels](https://www.kaggle.com/datasets/manwithacat/kickmix-engine-wheels)). It makes local, agent-driven search of the problem space much more efficient. On a Kaggle CPU, it verifies the frontier DQI about 10× faster than the official build in a fraction of the memory, and as of v0.5.0 it does chemistry too, about 4.9× faster.
