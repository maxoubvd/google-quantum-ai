# RezOps Briefing: DQI at 892 qubits three times over, a parts-assembled chemistry record, and a chemistry primer — 2026-11-01

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754793 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-11-01 17:28 UTC · 2 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps scans new public notebooks,
re-verifies their circuits, maintains our public circuits dataset and the notebook that blends the best of it, and
drafts this briefing. I set the direction and review and approve every post. Numbers are labelled: *observed*
(leaderboard), *author-reported* (the author's official-verifier run), or *re-verified by us* (the unchanged official
verifier, run locally). DQI circuits are ranked by their mean cost over 48 fixed draws, not by one draw. Window:
31 Oct 18:25 UTC to 01 Nov 17:03 UTC. Corrections welcome.

## 1. Leaderboard (01 Nov, 17:03 UTC, observed)

| | 31 Oct 16:34 UTC | 01 Nov 17:03 UTC |
| --- | --- | --- |
| Top score | 26.76409 (Breaking1800) | **27.15971** (Breaking1800) |
| Teams ≥ 26 | 5 | 8 |
| Teams ≥ 25 | 12 | 18 |
| Teams ≥ 24 | 36 | 94 |
| Teams ≥ 23 | 597 | 803 |
| Median | 22.70030 | 22.80578 |
| Teams | 1,810 | 1,934 |

The top ten are 25.90–27.16: Breaking1800, Schrödinger’s cat lived, Jack, Bingran You, Thinking..., trust the
transpiler, Simulation, 田村さんを囲む会, Adisak Sukul, Z³uperposition. Jack led for a few hours overnight; first and
second are now 0.35% apart.

The public plateau moved several times today (observed ties at 17:03): 23.83433 (20 teams), 23.83512 (25),
23.83707 (5) and 23.83790 (4, ranks 137–140, including our Hedge Wizardry notebook's portfolio, §7). They match public portfolios on a DQI draw with no failed shots: 23.81418 is bce2 + bdee + Heliosli's 0d75,
23.83433 swaps in yoyobar's 6d58, and 23.83512 also takes Heliosli's ffcc. The §2 parts assembly scores
23.83707 on the same draw, and with our HW-06 DQI file (§7) 23.83790.

Approximate medal lines on 1,934 teams (Kaggle's [medal rules](https://www.kaggle.com/progression/competitions); our
estimate):

| | Rank needed | Score at that rank now | Yesterday |
| --- | --- | --- | --- |
| Gold | ≤ 13 | 25.47565 | 24.96177 |
| Silver | ≤ 96 | 23.98765 | 23.52905 |
| Bronze | ≤ 193 | 23.83433 | 23.36706 |

The bronze line now sits on the newest public plateau. Silver is about 0.6% above the best public portfolio.

## 2. Public frontier by task (re-verified by us)

| Task | Best public | Qubits | Cost | Ratio to reference | Change since 31 Oct |
| --- | --- | --- | --- | --- | --- |
| Inverse square root | `bce2948ff24e` by [Ehan Mallick](https://www.kaggle.com/code/ehanmallick/invsqrt-20-65-two-3-to-2-and-merges-on-y6225) | 62 | 100,568 | 20.65 | unchanged |
| DQI | `6d58bd6d0232` by [yoyobar](https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits), on SomeJoint's flag carry | 892 | 3.8283 × 10¹¹ (mean of 48 draws) | 63.84 | −2.9% |
| Chemistry | public parts: [Dmitrii Gluzdov's](https://www.kaggle.com/code/dmitriigluzdov/quantum-circuit-compilation-from-source) loaders with [Heliosli's](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy) weights, on [Lera's](https://www.kaggle.com/code/marculera/df-thc-10x-compact-groups-from-the-hamiltonian) (9, 143, 40) table | 428 | 3.2682 × 10⁸ | 10.27 | −1.3% |

Best-per-task geometric mean: **23.83176** (re-verified by us, DQI at its 48-draw mean), up from 23.49622 yesterday.
All source notebooks are Apache-2.0. The chemistry row is an assembly: Gluzdov's cheaper angle loaders plus
Heliosli's weights, verified together by us. Neither published circuit alone reaches it (Heliosli's best single
family, `ffcc1ee812ac`, is 0.025% dearer). Our [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made)
picks these three.

## 3. What's new

**DQI: 900 → 892 qubits, then cheaper at 892.**
- [Ahmed Berat Özer](https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression) applied his
  history caches to SomeJoint's 900-qubit circuit: 576 Toffolis whose product was computed earlier become CNOTs from a
  copy kept in an idle wire, cleared by measurement-based uncomputation (`909410f50a54`, −0.14% on all 48 of our draws).
- [Ali Haydar Özdağ](https://www.kaggle.com/code/ahoirg/score-23-64254-reproducible-circuit-portfolio) rebuilt
  SomeJoint's decoder sources with a homogeneous route for groups 6–16 and a carry bank for the pole flags
  (`fb2da4f299c4`, 893 qubits, −1.68%).
- SomeJoint's [**flag carry**](https://www.kaggle.com/code/somejoint/quantum-circuits-23-77-local-gm-892-qubit-dqi)
  (`bdee18489fac`, 892 qubits): pole flags are measured on the way forward and rebuilt on the way back from the
  previous flag and the group's zero test, with the measurement repaid by a CZ. −0.92% on every draw. The author
  credits Özdağ's flag idea and says Claude wrote the rewrite.
- [yoyobar](https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits) (with GPT-6 Astra) added an
  "infinity-lead product reuse" to SomeJoint's circuit (`6d58bd6d0232`): 1,026 fewer static Toffolis at the same
  width, −0.25% on all 48 draws with the same failing shots.

**Chemistry: a (9, 143, 40) table, then weights and loaders.**
[Lera](https://www.kaggle.com/code/marculera/df-thc-10x-compact-groups-from-the-hamiltonian) took the compact-groups
factorization to C = 40 (`b163954eb5ac`): 428 qubits, 16,237 Toffolis, −1.11% against the best circuit on the C42
table. Weights-only records followed on the same words and loaders from
[kawh1111](https://www.kaggle.com/code/kawh1111/qcomp-harvest-v32), then
[Heliosli](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy) three times (latest `ffcc1ee812ac`,
λ_eff 47.0401); [sigmaborov's](https://www.kaggle.com/code/sigmaborov/22-85569-lb-exact-gauge-search) exact gauge
search tied one of them within noise. [Dmitrii Gluzdov](https://www.kaggle.com/code/dmitriigluzdov/quantum-circuit-compilation-from-source)
recompiled angle loaders 0–3 and 14 on the same words (−4 Toffolis); with Heliosli's weights that is the record in §2.

**Inverse square root:** no new record.

**Late arrivals (after 16:58 UTC, being re-verified):** [SomeJoint](https://www.kaggle.com/code/somejoint/quantum-circuits-23-85-local-gm-892-qubit-dqi)
published a new 892-qubit DQI, `b7ca9bce749d`. Its local official-verifier run passes at the default salt, and on its
zero-failure cost it is 0.28% below yoyobar's `6d58bd6d0232`. Our 48-draw comparison is running; tomorrow's briefing
will rank it. The same notebook ships a chemistry family equal in cost to the §2 parts assembly.

## 4. From the forum

No host posts in the window.

- **Licences for final submissions** ([754677](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754677), Bhupendra Singh). The rules require a winning submission
  and the code that generated it to be released under Apache-2.0 and to be reproducible. The post asks the authors of
  two widely reused parts to confirm their licence and point to the generator: the 62-qubit inverse-square-root
  circuits (yoyobar's and Ehan Mallick's), and the 893–898-qubit DQI circuits, whose plan files are public but whose
  wiring code is not. No replies yet. It matters to everyone building on public circuits, us included: our dataset
  lists each circuit's notebook, version and licence.
- **How do people plan to win?** ([754557](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754557), started by yoyobar). pjmathematician runs bulk AI experiments,
  then tunes by hand. Agniv Banerjee tried an agentic design loop and hit cost and scope problems even with a faster
  verifier (author-reported). Zero found LLM-driven research giving diminishing returns and evolutionary search
  impractical with a slow validator, and after trying our Rust
  [kickmix-engine](https://github.com/manwithacat/kickmix-engine) said a faster verifier could help a lot. hengck23 asked whether, without a better mathematical formulation, the rest of
  the competition is wiring optimisation. Our view: not yet. The two biggest moves this week were formulation
  changes, Lera's compact-groups factorization and the DQI flag reclaiming in §3.
- **What to learn from this competition** ([754652](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754652), hengck23): with most entrants using the same AI
  agents, what do people take away? Suggestions include steering agents well, distilling forum and notebook
  knowledge, and agent swarms. Geremie Yeo replied that the roadmap is being discovered as we go.
- **Replies to our last briefing** ([754605](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754605)): hengck23 agrees that AlphaTensor-Quantum would need
  modifying for the angle loaders and should not be expected to give large gains. Thanks to Syed Asad Ali for the
  kind words.
- Maaz Rahim Khan posted a beginner walkthrough of a load → optimise → verify pipeline ([754696](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754696)).

Three of these went into our negative-results ledger (NR-0164 to NR-0166, author-reported): Zero's validator-speed
limit on search, Agniv Banerjee's agentic-loop costs, and hengck23's AlphaTensor estimate.

## 5. Position risk: the final DQI draw

The final ranking uses a different DQI salt, so each DQI file gets one fresh draw of failing shots. For the frontier
DQI file (our estimate from our 48 draws): about 0.15 failing shots per draw, so about a 14% chance the final draw
has at least one failure, each costing 0.153% of score; the expected change is −0.02%. Teams within one DQI failure
(0.153%) below a medal line now (observed): 0 for gold, 11 for silver, 157 for bronze.

## 6. Prize watch (our estimate, from our leaderboard snapshots)

- **Leader:** Breaking1800 at 27.15971 (best submission 01 Nov 02:49 UTC); Schrödinger’s cat lived is 0.35% behind
  at 27.06470.
- **Mid-Competition Milestone Prize** (highest public score at 20 Nov 23:59 UTC, notebook published): 19.3 days left.
- **Longest Leader Prize:** Breaking1800 regained first between 03:18 and 07:51 UTC today. If it stays first without a
  break, nobody can catch its total from about 8 Dec (worst case; eligibility not checked).

## 7. Tools from us

- **New: a [chemistry primer](https://www.kaggle.com/code/manwithacat/vqcc-chemistry-primer-df-thc-task-explained).**
  One notebook that explains the chemistry task: what you submit, how the cost λ_eff × Toffolis × qubits is built and
  what a 1% saving in each factor is worth, the DF-THC-SOS factorization and its gates, angle loaders as QROM lookup
  tables, word tables and blocks, and which levers move which factor. Its numbers come live from our dataset.
  Corrections welcome.
- **New: [Hedge Wizardry](https://www.kaggle.com/code/manwithacat/vqcc-hedge-wizardry-notebook) HW-06.** On yoyobar's
  892-qubit DQI, 12 flag ANDs that were uncomputed with 4 Toffolis each (they borrow an in-use wire) are cleared by
  measurement-based uncomputation instead: measure the flag in the X basis, then a CCZ on its three inputs when the
  measured bit is 1. −42 Toffolis per shot, −0.0104% DQI cost over 48 draws with the same failing shots; the official
  verifier passed it on Kaggle. A value tracker that proves which wires are zero found the 12 sites. Small, but
  everyone building on this circuit can take it.
- **New: HW-07**, on top of HW-06. The 12 uncomputes of the six-input normalize flag, which also borrowed an in-use
  wire and cost 16 Toffolis every shot, are replaced by measurement-based uncomputation with a two-wire ladder: measure
  the flag; only when the outcome is 1, rebuild the AND on two clean wires to apply the phase fix, then clear them.
  About 2.25 Toffolis per site on average. −165 Toffolis per shot against HW-06; over 48 draws −0.051% DQI cost
  against yoyobar's `6d58bd6d0232`, cheaper on every draw with the same failing shots (local engine; official verifier
  PASS locally). One caveat: on its own default-salt draw this file has one failing shot, so a single run can look
  worse than HW-06 even though it is cheaper on average.
- The [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made) is slimmer: an
  at-a-glance table (each task's best, its change over 24 hours and 7 days, who holds it), the chemistry chart for
  the current word table and how the frontier moved between tables.
- [kickmix-engine](https://github.com/manwithacat/kickmix-engine) (Apache-2.0): a Rust re-implementation of kickmix's
  simulator, checked shot for shot against kickmix, with a `verify` command that runs the official verifier logic for
  all three tasks (and `--cross-check` against the official build). On Kaggle CPU it runs chemistry verification
  4.9× faster warm with 21% less memory (our benchmark, in its README). Wheels are also on Kaggle as
  [kickmix-engine-wheels](https://www.kaggle.com/datasets/manwithacat/kickmix-engine-wheels).
- The [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits): 209 public
  Apache-2.0 circuits with author, version, local verifier report and DQI 48-draw costs, plus verified parts
  assemblies for chemistry.
- New negative results from today go into [VQCC Negative Results](https://www.kaggle.com/datasets/manwithacat/vqcc-chemistry-negative-results)
  with its next version.

## 8. What we think it means (opinion)

DQI now moves by reclaiming flag qubits and then trimming Toffolis at fixed width: 8 qubits in a day, worth about 1.8%
on that task because cost grows with the square of the width. The next qubit at 892 looks hard: in our own attempt,
the normalize step of each epoch blocks 891. Chemistry has shifted from new tables to recombining parts: the record
is loaders from one author with weights from another. The leaders are about 14% ahead of the public frontier, much as
yesterday.
