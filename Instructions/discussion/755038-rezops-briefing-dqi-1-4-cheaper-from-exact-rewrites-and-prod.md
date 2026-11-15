# RezOps Briefing: DQI 1.4% cheaper from exact rewrites and product caching, chemistry weight refits — 2026-11-02

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755038 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-11-02 17:18 UTC · 2 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps scans new public notebooks,
re-verifies their circuits, maintains our public circuits dataset and the notebook that blends the best of it, and
drafts this briefing. I set the direction and review and approve every post. Numbers are labelled: *observed*
(leaderboard), *author-reported* (the author's own run), or *re-verified by us* (the unchanged official verifier, run
locally). DQI circuits are ranked by their mean cost over 48 fixed draws, not by one draw. Window: 01 Nov 17:00 UTC
to 02 Nov 17:10 UTC. Corrections welcome.

## 1. Leaderboard (02 Nov, 16:24 UTC, observed)

| | 01 Nov 16:24 UTC | 02 Nov 16:24 UTC |
| --- | --- | --- |
| Top score | 27.15971 (Breaking1800) | **27.99617** (Breaking1800) |
| Teams ≥ 26 | 8 | 10 |
| Teams ≥ 25 | 18 | 24 |
| Teams ≥ 24 | 91 | 169 |
| Teams ≥ 23 | 800 | 978 |
| Median | 22.80578 | 22.99251 |
| Teams | 1,931 | 2,043 |

The top ten are 26.12–28.00: Breaking1800, Jack, Schrödinger’s cat lived, Simulation, James & Tony, Thinking...,
Bingran You, Adisak Sukul, trust the transpiler, 田村さんを囲む会. Breaking1800 moved from 27.16 to 28.00 last night and
leads Jack by 0.74%.

The public plateau moved up about 0.47% today. The most-held scores near it (observed) are 23.94675 (33 teams),
23.94695 (24), 23.94667 (19) and 23.94682 (15). They match today's public portfolios on a DQI draw with no failed shot:
yoyobar's 892-qubit DQI (§3) with one of the chemistry weight refits, which differ from each other by less than
0.001%.

Approximate medal lines on 2,043 teams (Kaggle's [medal rules](https://www.kaggle.com/progression/competitions); our
estimate):

| | Rank needed | Score at that rank now | Yesterday |
| --- | --- | --- | --- |
| Gold | ≤ 14 | 25.71721 | 25.47565 |
| Silver | ≤ 102 | 24.14832 | 23.98765 |
| Bronze | ≤ 204 | 23.94701 | 23.83433 |

The bronze line now sits on the public plateau, as it did yesterday. Silver is about 0.86% above the best public
portfolio's expected score.

## 2. Public frontier by task (re-verified by us)

| Task | Best public | Qubits | Cost | Ratio to reference | Change since 01 Nov |
| --- | --- | --- | --- | --- | --- |
| Inverse square root | `bce2948ff24e` by [Ehan Mallick](https://www.kaggle.com/code/ehanmallick/invsqrt-20-65-two-3-to-2-and-merges-on-y6225) | 62 | 100,568 | 20.65 | unchanged |
| DQI | `e3928244a1c1` by [yoyobar](https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits), on busyaprime's rewrites of SomeJoint's decoder | 892 | 3.7763 × 10¹¹ (mean of 48 draws) | 64.72 | −1.36% |
| Chemistry | `5b4482b2e74c` by [Alex Chilton](https://www.kaggle.com/code/alexchilton/23-87-lb-verified-three-task-portfolio): weights on Gluzdov's loaders and Lera's (9, 143, 40) table | 428 | 3.2679 × 10⁸ | 10.27 | −0.01% |

Best-per-task geometric mean: **23.94157** (re-verified by us, DQI at its 48-draw mean), up from 23.83176 yesterday.
All source notebooks are Apache-2.0. Our [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made)
picks these three.

## 3. What's new

**DQI: four steps at 892 qubits.**
- [SomeJoint](https://www.kaggle.com/code/somejoint/quantum-circuits-23-85-local-gm-892-qubit-dqi) (Agniv Banerjee's
  team; yesterday's late arrival) `b7ca9bce749d`: about 1,883 fewer Toffolis per shot than their previous 892-qubit
  decoder (author-reported); −0.21% against yoyobar's `6d58bd6d0232` on all 48 of our draws, same failing shots.
- [busyaprime](https://www.kaggle.com/code/busyaprime/23-906-local-gm-exact-product-caching-892q-dqi) (Akmal Xodarev)
  **exact product caching** on that decoder (`3779a3431e9c`): a symbolic pass finds Toffolis whose two controls hold
  the same product, keeps one copy on a wire proven idle, turns later uses into CNOTs, and clears the copy by
  measurement-based uncomputation; an exact linear program chooses which chains fit in the free wires. 2,811 fewer
  static Toffolis, −0.65% on every draw. busyaprime credits Ahmed Berat Özer for the idea.
- busyaprime's [**exact DQI rewrites**](https://www.kaggle.com/code/busyaprime/23-910-local-gm-exact-dqi-rewrites-892-qubits)
  (`92750bd64525`, from the dataset [33-qcc-v8-members](https://www.kaggle.com/datasets/busyaprime/33-qcc-v8-members)):
  24 flag-clearing blocks are replaced by measurement-based uncomputation (12 of them with our HW-07 ladder), then two
  caching passes. −0.044% on our 48 draws. Emmanuel Mireku Mensah's
  [portfolio](https://www.kaggle.com/code/emmanuelmirekumensah/q-busya-23-909-portfolio) ships a near-identical file
  (`54d512bfcdac`, +0.0015%).
- [yoyobar](https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits) (with GPT-6 Astra) added
  **lead-product reuse** and more caching on top (`e3928244a1c1`): mean Toffolis per shot 399,970 → 398,139, −0.46% on
  our 48 draws with the same seven failing shots. Its later version (`148bc9bcc77d`) ties it. The author notes the
  decoder keeps approximate behaviour on rare branches and that finite tests do not prove equivalence on every input
  (author-reported); our draws found no new failures.

Dmitrii Gluzdov and Ahmed Berat Özer also published new 892-qubit DQI files; none beats the line above.

**Chemistry: weight refits on one table.** Every new chemistry record keeps the same words and the same 16,233
Toffolis and refits only the table weights (λ_eff 47.0393 → 47.0356), each step a few thousandths of a percent:
[Isabe11aM00re](https://www.kaggle.com/code/isabe11am00re/qcc-v49-submission) twice,
[Heliosli](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy) three times,
[busyaprime](https://www.kaggle.com/code/busyaprime/23-910-local-gm-exact-dqi-rewrites-892-qubits),
[Ahmed Berat Özer](https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression) twice, then
[Alex Chilton](https://www.kaggle.com/code/alexchilton/23-87-lb-verified-three-task-portfolio) (the record in §2).
Later refits by [yoyobar](https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits) and
[sigmaborov](https://www.kaggle.com/code/sigmaborov/22-85569-lb-exact-gauge-search) tie it within 1e-5, so the earlier
file keeps the pick. The whole day's chemistry gain is about 0.01%.

**Inverse square root:** no new record. Ehan Mallick confirmed Apache-2.0 for his two patches and pointed to the
notebook that rebuilds `bce2948ff24e` byte for byte (§4).

## 4. From the forum

No new host posts in the window.

- **Licences and generators for final submissions** ([754677](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754677)).
  Ehan Mallick released his invsqrt patches under Apache-2.0, linked the generator that rebuilds `bce2948ff24e` and
  passes the official verifier on every input, and noted that the base circuit needs yoyobar's own confirmation;
  yoyobar said a formal licence will be added in the next notebook version. A follow-up asks the hosts whether a new
  Apache-2.0 circuit whose design follows unlicensed earlier work can win, whether a generator may start from an
  already-compiled Apache-2.0 circuit, and whether a GPU fit must be rerun. Rui Codex asked the same about trained
  chemistry tables in [752010](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752010),
  where the host's earlier answer stands: a generator that produces the exact circuits is required.
- **LB anatomy** ([754956](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754956), Harshit Kumar):
  a snapshot with open data (949 of 2,018 teams at 23+, 153 at 24+) and a list of 276 public notebooks with the scores
  their titles claim. pjmathematician replied on submission timing at the top.
- **Scoring time** ([754805](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754805), Tony Li):
  several top teams report scoring takes 20–60 minutes per submission (author-reported).
- **New tools and guides:** yoch's [SymQ](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754815)
  checks reversible rewrites in context, including dirty ancillas and measurement-based corrections; snuffkin's
  [introductory guide](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754970)
  to fault-tolerant circuit optimisation (Japanese PDF, CC BY 4.0).
- **Agents** ([754652](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754652),
  [754914](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754914)):
  hengck23 suggests testing yourself by reproducing a high-scoring notebook's result without showing it to the agent;
  Agniv Banerjee recommends "clean" agents with only the problem description and kit. Ehan Mallick asks what people
  score per task ([754965](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754965));
  the §2 ratios are the public answer.

## 5. Position risk: the final DQI draw

The final ranking uses a different DQI salt, so each DQI file gets one fresh draw of failing shots. For the frontier
DQI file (our estimate from our 48 draws): about 0.15 failing shots per draw, so about a 14% chance that the final
draw has at least one failure, each costing 0.153% of score; the expected change is −0.02%. Teams within one DQI
failure (0.153%) below a medal line now (observed): 1 for gold, 17 for silver, 149 for bronze.

## 6. Prize watch (our estimate, from our leaderboard snapshots)

- **Leader:** Breaking1800 at 27.99617 (best submission 01 Nov 21:05 UTC); Jack is 0.74% behind at 27.79189.
- **Mid-Competition Milestone Prize** (highest public score at 20 Nov 23:59 UTC, notebook published): 18.3 days left.
- **Longest Leader Prize:** Breaking1800 has led since 01 Nov 21:20–22:37 UTC. If it stays first without a break,
  nobody can catch its total from about 9 Dec (worst case; eligibility not checked).

## 7. Tools from us

- **[Hedge Wizardry](https://www.kaggle.com/code/manwithacat/vqcc-hedge-wizardry-notebook) HW-09.** On busyaprime's
  `92750bd64525`, 24 AND chains borrow in-use wires as conditionally clean ancillae (Khattar and Gidney, 2024): a
  borrowed wire is used as if clean where the circuit guarantees it is clean whenever it matters, and our value
  tracker proves every such condition first. −294 Toffolis per shot, −0.074% DQI cost on all 48 draws with the same
  failing shots; the official verifier passed it on Kaggle (notebook version 10). yoyobar's `e3928244a1c1` has since
  moved the frontier further, so HW-09 is a small, citable step rather than the current best. HW-08 is also in the
  dataset as a reference record.
- The [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made) now shows its
  score as the best public portfolio's expected value, with error bounds for the DQI draw (no failure, one failure,
  two failures) instead of a single number.
- The [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits): 237 public
  Apache-2.0 circuits (28 new today) with author, version, local verifier report and 48-draw DQI costs.
- [kickmix-engine](https://github.com/manwithacat/kickmix-engine) (Apache-2.0), our Rust re-implementation of
  kickmix's simulator with a `verify` command; it ran today's 48-draw rankings.
- Our own recent nulls (DQI width reductions at 892 qubits and chemistry word-table compression) go into
  [VQCC Negative Results](https://www.kaggle.com/datasets/manwithacat/vqcc-chemistry-negative-results) with its next
  version.

## 8. What we think it means (opinion)

DQI is now won by exact bookkeeping at fixed width: find repeated products, prove which wires are free, and clear
temporaries by measurement. Three authors built on each other's files within a day, and the published notebooks
explain their methods. Chemistry has settled into weight refits worth thousandths of a percent each; the next real
step there likely needs a new table or new loaders. The leaders are about 17% ahead of the public frontier, up from 14%
yesterday.


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**(deleted user)** · 2026-11-03 06:30 UTC · 0 votes

*(deleted)*
