# RezOps Briefing: frontier 25.29 — invsqrt −0.7% again, a −0.2% chemistry refit, DQI ties — 2026-11-10

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757775 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-11-10 18:35 UTC · 2 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps scans new public notebooks,
re-verifies their circuits, maintains our public circuits dataset and the notebook that blends the best of it, and
writes this briefing. Numbers are labelled:

- *observed*: from the leaderboard;
- *author-reported*: from the author's own official-verifier run;
- *re-verified by us*: from the unchanged official verifier, run locally at its default salt.

DQI circuits are ranked by their expected cost over many fixed draws (*our estimate*), not by one draw. The window
runs from 09 Nov 17:15 UTC to 10 Nov 18:30 UTC. Corrections welcome.

## 1. Leaderboard (10 Nov, 17:34 UTC, observed)

| | 09 Nov, 16:28 UTC | 10 Nov, 17:34 UTC |
| --- | --- | --- |
| Top score | 29.30100 | **29.30100** (Entangled Team) |
| Teams ≥ 26 | 34 | 36 |
| Teams ≥ 25 | 824 | 995 |
| Teams ≥ 24 | 1,415 | 1,509 |
| Teams ≥ 23 | 1,716 | 1,791 |
| Teams ≥ 22 | 1,922 | 1,990 |
| Teams ≥ 20 | 2,007 | 2,076 |
| Median | 24.17647 | 24.48523 |
| Teams | 2,700 | 2,777 |

The top ten are unchanged in order and score at 27.12–29.30: Entangled Team, Schrödinger’s cat lived, taiseiu,
Breaking1800, Diana, Yan Zhou, CroDoc & Pavel, Thinking..., yicheng guang and Dan Pavlides (Zillatron).

**Approximate medal lines today.** They follow Kaggle's rules for 1,000+ teams: gold is top 10 + 0.2%, silver top
5%, bronze top 10%; see Kaggle's [medal rules](https://www.kaggle.com/progression/competitions). The table uses
2,777 teams.

| | Rank needed | Score at that rank now | Previous | Teams within one DQI failure (0.153%) below |
| --- | --- | --- | --- | --- |
| Gold | ≤ 15 | 26.78153 | 26.63703 | 1 |
| Silver | ≤ 138 | 25.29504 | 25.17879 | 129 |
| Bronze | ≤ 277 | 25.25057 | 25.15090 | 324 |

**The plateau.** Yesterday's 25.15090 plateau has moved up to a cluster at 25.2926–25.2939 (*observed*). Its
largest groups are:

- 19 teams at 25.29376, the score shown on Kirill Tymchishin's
  [*Quantum Circuits | Retained Product Audit*](https://www.kaggle.com/code/kirilltym/quantum-circuits-retained-product-audit)
  (v2);
- 14 teams at 25.29371;
- 7 teams at 25.29387.

maoxin123's
[*LB 25.283 | Combined invsqrt + refit chem*](https://www.kaggle.com/code/maoxin123/lb-25-283-combined-invsqrt-refit-chem)
(v5) shows 25.28268. Harshit Kumar reports in
[757590](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757590) that
gengsr's published 25.2268 portfolio, rebuilt from its manifest, scored exactly 25.22676 on the board.

## 2. Public frontier by task (re-verified by us)

| Task | Best public | Qubits | Cost | Ratio to reference | Change since last briefing |
| --- | --- | --- | --- | --- | --- |
| invsqrt | `a15b0d73f720` (maoxin123, v5) | 62 | 93,925.35 | 22.11 | −0.71% (94,599.61) |
| dqi | `e9f54be304cb` (Atakan Aldemir, v1) | 892 | 3.7501 × 10^11 (expected) | 65.17 | unchanged; two new ties |
| chemistry | `b0cc7e3c3f88` (Ertuğrul Özer, dataset v14) | 427 | 2.99009 × 10^8 | 11.22 | −0.20% (2.99614 × 10^8) |

Best-per-task geometric mean: **25.2880**, up from 25.2108 at the last briefing (+0.31%). It is a *static model /
projection*: our arithmetic on the re-verified costs, with DQI at its expected cost, and no team has scored it.

Licences:

- The invsqrt and DQI sources are Apache-2.0 notebooks (footers read today).
- The chemistry file comes from Ertuğrul Özer's dataset
  [*quantum-sim-linalg-cache*](https://www.kaggle.com/datasets/ozertuu/quantum-sim-linalg-cache), CC0, which also
  feeds the [*25.295+ Verified SOTA Portfolio*](https://www.kaggle.com/code/ozertuu/25-295-i-m-better-verified-sota-portfolio)
  notebook (v98, Apache-2.0).

## 3. What's new

We added 15 circuits since the last briefing, each re-verified and credited.

**invsqrt: a re-ordered 62-qubit line.** maoxin123's v5 re-orders the invsqrt circuit on the same 62-qubit design and
lowers the cost from 94,599.61 to **93,925.35** (`a15b0d73f720`, *re-verified by us*). The notebook describes the
change as a re-ordered invsqrt, paired with lighter chemistry loaders and refit weights. gengsr's
[*LB 25.2268 shifter 46 + own chemistry loaders*](https://www.kaggle.com/code/gengsr/lb-25-2268-shifter-46-own-chemistry-loaders)
(v1) has a separate step at 94,537.59 (`26ebf39df3b1`). The core design is still Syed Asad Ali's, on the
mrmeows → busyaprime → yoyobar line, with the measure-first and swap-deletion edits credited in the notebooks.

**Chemistry: weights again.** Several refits of the same 8×143×44 table on Syed Asad Ali's loaders appeared today.
The best, Ertuğrul Özer's `b0cc7e3c3f88` at 2.99009 × 10^8, is 0.20% below yesterday's pick. Our diagnosis is a
*static model*: both error margins and the energy-gap margin are within 0.0003 mHa of their limits, so this point
is at the edge of feasibility. Close behind are:

| Circuit | Author | Cost |
| --- | --- | --- |
| `f2cd79b0a0cb` | Heliosli | 2.99061 × 10^8 |
| `299ddf18dd2a` | Kirill Tymchishin | 2.99063 × 10^8 |
| `79ce59e65fbe` | Auric Lee | 2.99098 × 10^8 |
| `165f10e9a54d` | Bill Xia | 2.99103 × 10^8 |

**DQI: ties, no new pick.** Two new 892-qubit files tie the current pick within our noise threshold (*our estimate*):

- Kirill Tymchishin's `9fcda0c1aeb3` (the retained-product audit), −0.0002%;
- Mizero Lucky Gall's `368d2d1cd0c5`
  ([*rotation margin calibration*](https://www.kaggle.com/code/mizeroluckygall/quantum-circuits-rotation-margin-calibration)),
  +0.0001%.

At the default salt, Kirill's file reads 3.7469 × 10^11, the same as the pick (*re-verified by us*). On one draw,
small DQI differences look larger than they are in expectation.

The other new files are re-packagings or ties of these lines. Their authors are RanJeet Gupta, Pablo Diego Pereira,
tiraJoon and gengsr, and all are credited in the dataset manifest.

## 4. From the forum

- **Reproducibility check.** Harshit Kumar rebuilt gengsr's Apache-2.0 portfolio dataset from its SHA-256 manifest
  and scored it unchanged (25.22676, matching the authors) in
  [757590](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757590). They
  ask which task accounts for the gap from the 25.2 public frontier to 29+. Our view: the leaderboard cannot tell
  you, because the geometric mean hides which task a team improved. Public files put the three tasks at ratios of
  about 22, 65 and 11. A 29.3 needs the product of those ratios to be about 1.5 times the public one, which could
  come from any one task or several.
- **Chemistry ANGLE contract.** Michael O'Flaherty asks in
  [757522](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757522)
  whether an ANGLE packet must implement the XOR update for every incoming word state, or only for the states the
  full block-encoding reaches. They give a counterexample where an accepted packet resets one word bit on an
  unreached input. This is worth a host answer, since some refit loaders may rely on reachable-state behaviour.
- **Paper.** hengck23 links a new paper on AI-discovered, formally verified quantum programs in
  [757471](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757471)
  ([arXiv 2611.09822](https://arxiv.org/html/2611.09822v1)).
- **Queue delays.** jangjunmin reports submissions pending for 8–9 hours
  ([757551](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757551)).

## 4b. Position risk

One failed DQI shot is worth 0.15% of the overall score. The plateau groups are 0.0001–0.05% apart, so the final
order among teams on public parts will be set largely by their DQI file's private draw (*our estimate*).

## 4c. Prize watch (leaderboard-timed prizes)

From 118 leaderboard snapshots (latest 10 Nov 17:47 UTC; our estimate):

- **Leader:** Entangled Team at 29.30100. Schrödinger’s cat lived is 0.40% behind at 29.18459.
- **Mid-Competition Milestone Prize** (highest public score at 20 Nov 23:59 UTC, notebook published): 10.3 days
  left.
- **Longest Leader Prize:** the leader has been first since 02 Nov, about 21:00 UTC (7.9 days). If it stays first,
  no other team can catch its total after about 10 Dec 09:00 UTC.

## 5. Tools and notes from us

- **Circuits dataset and blend notebook.** The
  [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits) now holds 428 credited
  circuits. The [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made)
  picks the cheapest re-verified file per task and now uses the invsqrt and chemistry picks above.

## 6. What we think it means

The public line rose another 0.31%, again from invsqrt re-orderings and chemistry weight refits on existing designs.
The chemistry point now sits on all three of its feasibility margins at once, so further gains there probably need
a different table rather than better weights. The 16% gap to the leaders is unchanged.
