# RezOps Briefing: a new chemistry family from compact groups, a 900-qubit DQI, and the public blend at 23.505 — 2026-10-31

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754605 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-31 18:28 UTC · 0 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps scans new public notebooks,
re-verifies their circuits, maintains our public circuits dataset and the notebook that blends the best of it, and
drafts this briefing. I set the direction and review and approve every post. Numbers are labelled: *observed*
(leaderboard), *author-reported* (the author's official-verifier run), or *re-verified by us* (the unchanged official
verifier, run locally). DQI circuits are ranked by their mean cost over 48 fixed draws, not by one draw. Window:
30 Oct 17:09 UTC to 31 Oct 18:25 UTC. Corrections welcome; one from today is in §5.

## 1. Leaderboard (31 Oct, 18:25 UTC, observed)

| | 30 Oct 17:09 UTC | 31 Oct 18:25 UTC |
| --- | --- | --- |
| Top score | 26.76409 (Breaking1800) | **26.76409** (Breaking1800) |
| Teams ≥ 26 | 3 | 5 |
| Teams ≥ 25 | 6 | 12 |
| Teams ≥ 24 | 21 | 36 |
| Teams ≥ 23 | 236 | 614 |
| Teams ≥ 22.5 | 842 | 1,008 |
| Median | 22.48595 | 22.70119 |
| Teams | 1,696 | 1,818 |

The top ten are 25.28–26.76: Breaking1800, Schrödinger’s cat lived, Jack, Thinking..., Bingran You, trust the
transpiler, Earthwalker, Cache Me If You Can, 田村さんを囲む会, VAGO Superposition. Schrödinger’s cat lived is 0.44%
behind first.

The public plateau climbed all day; the ties show each step (observed): 23.36477 (22 teams), 23.40252 (7),
23.46148 (5), 23.48785 (6) and **23.50521** (4 teams, ranks 98–101). The last is the public blend described in §2.

Approximate medal lines on 1,818 teams (Kaggle's [medal rules](https://www.kaggle.com/progression/competitions); our
estimate):

| | Rank needed | Score at that rank now | Yesterday |
| --- | --- | --- | --- |
| Gold | ≤ 13 | 24.96177 | 24.74027 |
| Silver | ≤ 90 | 23.52905 | 23.13695 |
| Bronze | ≤ 181 | 23.36706 | 23.03386 |

The best public portfolio now sits about 0.1% below the silver line.

## 2. Public frontier by task (re-verified by us)

| Task | Best public | Qubits | Cost | Ratio to reference | Change since 30 Oct |
| --- | --- | --- | --- | --- | --- |
| Inverse square root | `bce2948ff24e` by [Ehan Mallick](https://www.kaggle.com/code/ehanmallick/invsqrt-20-65-two-3-to-2-and-merges-on-y6225), on yoyobar's 62-qubit line | 62 | 100,568 | 20.65 | −0.37% |
| DQI | `38f0ac173b0f` by [SomeJoint](https://www.kaggle.com/code/somejoint/quantum-circuits-23-46-local-gm-900-qubit-dqi) (Agniv Banerjee) | 900 | 3.9440 × 10¹¹ (mean of 48 draws) | 61.97 | −2.2% |
| Chemistry | `02f8f7d9d75a`, our [Hedge Wizardry](https://www.kaggle.com/code/manwithacat/vqcc-hedge-wizardry-notebook) HW-04 refit of [Lera's](https://www.kaggle.com/code/marculera/df-thc-10x-compact-groups-from-the-hamiltonian) (9, 143, 42) table | 428 | 3.3103 × 10⁸ | 10.14 | −2.0% |

Best-per-task geometric mean: **23.49622** (re-verified by us, DQI at its 48-draw mean), up from 23.13587 yesterday.
Our [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made) v59 picks these three
and scored **23.50521** on the leaderboard (observed), the same as its local run. All source notebooks are Apache-2.0.

## 3. What's new

**Chemistry: a new family, built from the Hamiltonian.** [Lera](https://www.kaggle.com/code/marculera/df-thc-10x-compact-groups-from-the-hamiltonian)
seeds the factorization from pivoted-Cholesky leaves clustered by orbital into *compact groups* (each group touches few
orbitals and overlaps the others little), then fits directly at the scoring shape with no merge or shrink. At
(R, B, C) = (9, 143, 43) that gave 428 qubits and 10.07× the reference, −1.22% against the best circuit on the previous
word table (re-verified by us). Lera then published (9, 143, 42): λ_eff 47.0809, 16,464 Toffolis, −0.48% more. The C43
table was refined the same day by [Atakan Aldemir](https://www.kaggle.com/code/atakanaldemir/qcc-gap-3-residue-repair-and-c43-native-compiler)
(offset integration, with Kirill Tymchishin's exact-pair polish), [Ahmed Berat Özer](https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression)
(column and Jacobi gauge) and [sigmaborov](https://www.kaggle.com/code/sigmaborov/22-85569-lb-exact-gauge-search). Our
HW-04 is a constrained refit of Lera's C42 table with the words and angle loaders unchanged: λ_eff 46.9767, −0.22%.

Earlier in the window, weights-only records on Atakan Aldemir's 8-bit table came from
[isabe11am00re](https://www.kaggle.com/code/isabe11am00re/qcc-v47-submission) (a per-rank orthogonal column gauge,
which leaves error₁, error₂ and the gap unchanged), our HW-03, sigmaborov and [Heliosli](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy).

**DQI: three steps on the 888/900-qubit lines.**
- [yoyobar](https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits): first 889 → 888 by omitting the gi = 8
  pole branch (`25db915d86ab`; over our 48 draws it adds no failures), then "F3", the F flag spilled (`e20e6abea7e2`,
  −56 Toffolis).
- [Ahmed Berat Özer](https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression): 130 exact history
  caches (`fe98677c61ef`). A Toffoli whose product was computed earlier becomes a CNOT from a copy kept in a wire that
  sits idle in between; the copy is cleared by measurement-based uncomputation. −0.031%, same failing shots.
- [SomeJoint](https://www.kaggle.com/code/somejoint/quantum-circuits-23-46-local-gm-900-qubit-dqi): a "narrow return"
  for the reverse of extraction (`38f0ac173b0f`). Each node keeps its terminal value while undoing its chain, retires
  the old numerator in two halves and reuses the freed wires; two inversions that existed only to repay signs are gone.
  11 more qubits, about 15,900 fewer Toffolis per shot: −1.10% against Özer's circuit on all 48 of our draws, with the
  same rare failures. The author credits GPT-6 Astra Pro for the idea and Claude for the build.

**DQI failures diagnosed.** [Atakan Aldemir](https://www.kaggle.com/code/atakanaldemir/qcc-gap-3-residue-repair-and-c43-native-compiler)
traced a captured failure to a degree-5 → 2 Euclid gap in residue conversion and repaired it (guarded Bézout
completion; 891 qubits, +797 Toffolis). We checked it against the 12 failing shots of the current best over 48 draws:
none is that case, which is real but rare. Eight of the twelve are chance 24-bit zeros on ordinary nodes in groups
4–7.

**Inverse square root: 62 qubits by two routes.** [Ehan Mallick](https://www.kaggle.com/code/ehanmallick/invsqrt-20-65-two-3-to-2-and-merges-on-y6225)
rewrote two windows where three Toffolis only XOR into an accumulator as two (a GF(2) quadratic-form rank argument),
−2 Toffolis. [Leo Tang's](https://www.kaggle.com/code/leotang09/vqcc-t7421-submission) 62-qubit circuit is the earlier
67482db6c538 with 2,838 Clifford gates cancelled; our [Tangfastic](https://www.kaggle.com/code/manwithacat/tangfastic)
notebook explains it. A sweep we ran over every accumulator window of `bce2948ff24e` found no further merge that fits
in 62 qubits.

## 4. Method notes and what didn't work

- **Ranking DQI.** The hosted verifier draws with its own salt, so one local run is one sample. Our dataset and blend
  rank DQI by mean cost over 48 fixed draws; within 1e-5 relative, the earliest public circuit keeps the pick.
- **Chemistry levers checked against the official verifier.** Symmetry shifts (multiples of N̂ and N̂²) are
  overwritten by the verifier's baseline shift; orbital rotations fail error₁ and error₂ because the reference
  integrals are fixed in their basis. Both are logged as observations so nobody spends time on them.
- **Angle loaders (LUTs).** hengck23 asked ([discussion](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754466))
  whether AlphaTensor-Quantum is worth trying on the chemistry lookup circuits. Our numbers: the angle-loader cache is
  already at its minimum; a counting estimate puts the ceiling near −5% for a typical table of this structure, but
  about 130 million local moves found nothing cheaper than the public loaders, and AlphaTensor-Quantum minimises
  T-count by trading Toffolis for T and CS gates, which these loaders cannot use. We would expect 1% or less.
- **Terminology.** The fifteen ANGLE_k circuits are QROM lookup tables (LUTs); we call them angle loaders. Word tables
  now carry stable labels in our dataset (W01–W20, by first appearance) so their history can be followed.
- New negative results and observations from today are in our ledger (NR-0131 to NR-0152) and go into the
  [VQCC Negative Results](https://www.kaggle.com/datasets/manwithacat/vqcc-chemistry-negative-results) dataset with its
  next version.

## 5. Corrections

Our circuits dataset missed Lera's C43 circuits for about three hours (11:39 to about 14:45 UTC): they sat in an output
folder named `dfthcsos_3`, and our scanner only recognised `dfthcsos/`. Fixed; anyone who blended from our dataset in
that window had the older chemistry.

## 6. Prize watch (our estimate, from 35 leaderboard snapshots)

- **Leader:** Breaking1800 at 26.76409 since 30 Oct 04:23 UTC; Schrödinger’s cat lived is 0.44% behind at 26.64635.
- **Mid-Competition Milestone Prize** (highest public score at 20 Nov 23:59 UTC, notebook published): 20 days left.
- **Longest Leader Prize:** if Breaking1800 stays first without a break, nobody can catch its total from about
  6 Dec (worst case; eligibility not checked).

## 7. Tools from us

- The [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits): 191 public
  Apache-2.0 circuits with author, version, local verifier report and DQI 48-draw costs; new today: stable word-table
  labels with lineage (which angle loaders each table rebuilt) and a leaderboard summary file (no team names).
- The [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made) now shows which
  word tables have held the chemistry frontier and how the public frontier compares with the leaderboard.
- [Hedge Wizardry](https://www.kaggle.com/code/manwithacat/vqcc-hedge-wizardry-notebook): small public improvements
  on top of others' circuits, released at once; each release leads with the technique it demonstrates.

## 8. What we think it means (opinion)

Chemistry moved by changing the shape of the factorization rather than polishing weights on a fixed one, and the
gauge-type gains on the old table had shrunk to parts per million. DQI moved by redesigning one stage, the reverse of
extraction, and paying width for Toffolis. The leaders are still about 14% ahead on score, and the public frontier is
now within 0.1% of the silver line.


---

## Comments (2 threads)

Replies are indented with `>`, one level per reply depth.

**Syed Asad Ali** · 2026-10-31 19:39 UTC · 0 votes

Thank you for these daily doses! Helps a lot in analysis and understanding!

**hengck23** · 2026-10-31 18:46 UTC · 0 votes

AlphaTensor will need to be modified, but don't expect large gains
