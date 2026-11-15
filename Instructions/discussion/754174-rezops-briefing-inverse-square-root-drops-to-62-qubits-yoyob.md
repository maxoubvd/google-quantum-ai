# RezOps Briefing: inverse square root drops to 62 qubits (yoyobar), the DQI failures get fixed, and the public blend reaches 23.135 — 2026-10-30

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754174 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-30 17:12 UTC · 0 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps scans new public notebooks, re-verifies their circuits, maintains our public circuits dataset and the notebook that blends the best of it, and drafts this briefing. I set the direction and review and approve every post. Numbers are labelled: *observed* (leaderboard), *author-reported* (the author's official-verifier run), or *re-verified by us* (the unchanged official verifier, run locally at its default salt). Window: 29 Oct 17:00 UTC to 30 Oct 17:09 UTC. Corrections welcome. One correction from today is already below (§4).

## 1. Leaderboard (30 Oct, 17:09 UTC, observed)

| | 28 Oct 21:44 UTC | 30 Oct 17:09 UTC |
| --- | --- | --- |
| Top score | 26.14113 | **26.76409** (Breaking1800) |
| Teams ≥ 26 | 1 | 3 |
| Teams ≥ 25 | 2 | 6 |
| Teams ≥ 24 | 7 | 21 |
| Teams ≥ 23 | 29 | 236 |
| Teams ≥ 22.5 | 448 | 842 |
| Median | 21.10349 | 22.48595 |
| Teams | 1,461 | 1,696 |

The top ten are 24.74–26.76: Breaking1800, Jack, Schrödinger’s cat lived, Thinking..., Bingran You, trust the transpiler, Kevin Bogner, Adisak Sukul, Cache Me If You Can, Z³uperposition. Jack held first briefly early this morning, and Breaking1800 retook it at 04:23 UTC with 26.76409.

The public plateau rose twice today, in steps you can read from the ties (observed):
- **22.99252**: 52 teams. yoyobar's 62-qubit inverse square root with the public DQI and chemistry.
- **23.00037**: 15 teams. The same, plus newer chemistry weights.
- **23.13505**: 17 teams by 17:09 UTC (ranks 95–111). The public DQI failure fix and chemistry refit (§3), which our blend notebook picks up. The silver line has moved just above it (23.13695 at rank 84).

Approximate medal lines at 17:09 UTC, on 1,696 teams (Kaggle's [medal rules](https://www.kaggle.com/progression/competitions); our estimate):

| | Rank needed | Score at that rank now | Two days ago |
| --- | --- | --- | --- |
| Gold | ≤ 13 | 24.74027 | 23.71909 |
| Silver | ≤ 84 | 23.13695 | 22.71082 |
| Bronze | ≤ 169 | 23.03386 | 22.70030 |

Ties are broken by the earlier submission ([rules](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/rules)).

## 2. Public frontier by task (re-verified by us)

| Task | Best public | Qubits | Cost | Ratio to reference | Change since 28 Oct |
| --- | --- | --- | --- | --- | --- |
| Inverse square root | `32d3e7d26c49` by [Dmitrii Gluzdov](https://www.kaggle.com/code/dmitriigluzdov/quantum-circuit-compilation-from-source), a hair below Ahmed Berat Özer's `6dbb98b34c3b`, both on [yoyobar's](https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits) 62-qubit `a568bc793481` | 62 | 100,941 | 20.57 | −1.38% (was 63 qubits, 102,349) |
| DQI | `2ed7ccf189f1`, [VQCC Quantum Voodoo](https://www.kaggle.com/code/manwithacat/vqcc-quantum-voodoo), on the SomeJoint / heliosli lineage | 888 | 4.0334 × 10¹¹ | 60.59 | −1.69% (was hengck23's 889-qubit `7a378f`) |
| Chemistry | `4eabe756bc36` by [Ahmed Berat Özer](https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression): new weights on Atakan Aldemir's word table and the existing loaders | 430 | 3.3780 × 10⁸ | 9.93 | −0.40% |

Best-per-task geometric mean: **23.13587** (re-verified by us), up from 22.86836 at the last briefing. That is just short of the silver line at 17:09 UTC (23.13695). All three source notebooks are Apache-2.0.

## 3. What's new

**Inverse square root: 63 → 62 qubits (yoyobar).** [yoyobar's *LB 22.99 | invsqrt 101 - 62 qubits*](https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits) reworks the encoded lookup, normalization and tail of busyaprime's circuit (on mrmeows / team selenium's design). The author says it was developed with GPT-6 Astra. Measured by us, it keeps the same angles and first 239 operations and makes every 63-wide stage one wire narrower, for +4.5 Toffoli-equivalents against a qubit worth about 26. [Ahmed Berat Özer](https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression) then replaced one controlled swap by an affine identity (a free wire, two CX, a measurement and a CZ fix-up) for −1 Toffoli. [heliosli](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy) republished both with careful re-verification and credit. Details and our exact check of 6dbb are in topic 754152.

**DQI: fewer failures by design.** A DQI file's cost includes a repetition factor R set by how many shots fail on a fixed random draw. Each extra failure costs about 0.46% of DQI. The draw is seeded from the file's own hash, so any byte change re-rolls it. We traced most failures of the 888-qubit line to one mechanism: a chance 24-bit zero entering a field inverse. Rebuilding four extraction groups as homogeneous chains removes it and also saves 1,460 Toffolis. Over 48 draws, failing draws go from 13 to 6 (re-verified by us; local engine). The files and the explanation are public in [VQCC Quantum Voodoo](https://www.kaggle.com/code/manwithacat/vqcc-quantum-voodoo), which also covers the 889-qubit splice and the 888-qubit release.

**Chemistry: a weights market.** Today's chemistry gains are almost all new weights on the same word table and loaders, each 0.001–0.2%. The authors were heliosli (four steps, including `7c6c875a7516` at −0.204%), Dmitrii Gluzdov, sigmaborov, Aamod06, Ahmed Berat Özer, Pablo Diego Pereira and Dmitry Abulkhanov (re-verified by us). Two new structures appeared: [Den Pugovkin's *Twelve bits fails*](https://www.kaggle.com/code/denpugovkin/twelve-bits-fails) (423 qubits, a different trade) and [Leo Tang's](https://www.kaggle.com/code/leotang09/vqcc-t7421-submission) new word table (742 qubits). Neither is at the frontier yet.

**A new independent DQI line (Leo Tang).** [Leo Tang's *VQCC T7421 Submission*](https://www.kaggle.com/code/leotang09/vqcc-t7421-submission) ships a 1,797-qubit DQI decoder that is not from any public lineage. It passes the official verifier (re-verified by us) but costs about 6× the frontier (ratio about 10 vs 60), and the file is 1.46 GB, so we list and credit it but don't ship the file in the dataset. An independent design is worth watching: it may turn into something interesting as it narrows.

All 42 new circuits are in the dataset (one by reference only), with credits and local verifier reports.

## 4. From the forum

- **754152 (ours), with a correction.** I first credited the 62-qubit circuit to heliosli. yoyobar pointed out that it is theirs, and heliosli's own notebook says so. The post, the dataset credits and our notebook are corrected. Apologies again to yoyobar. Our pipeline credited whichever notebook ran most recently first; it now flags any new circuit found in several notebooks for a human origin check.

## 4c. Prize watch (our estimate, from 17 leaderboard snapshots)

- **Leader:** Breaking1800 at 26.76409; Jack is 1.91% behind at 26.26331.
- **Mid-Competition Milestone Prize** (highest public score at 20 Nov 23:59 UTC, notebook published): 21 days left.
- **Longest Leader Prize** (time at #1 is combined across a team's stints, per the overview). Jack's 1.2 hours at #1 this morning (03:11–04:23 UTC) ended Breaking1800's first stint, which puts a number on it.
  - **Breaking1800's banked time:** about 3.5–5.0 days from that first stint (it took #1 between 25 Oct 02:24 and 26 Oct 14:25 UTC by our snapshots), plus 0.5 days since 04:23 today. That totals **about 4.0–5.5 days**.
  - **Jack's banked time:** Jack also led before that (top of our 24 Oct 22:25 export), so it holds roughly 0.2–1.7 days plus this morning.
  - **If Breaking1800 now holds #1 without a break,** nobody can catch its total from about **30 Nov – 2 Dec**. Our estimate; it ignores anyone who led before our first export, and eligibility needs public notebooks, which we have not checked.

## 5. Tools and notes from us

- The [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits): 142 public Apache-2.0 circuits, each with its author, version and local verifier report. Credits corrected today (§4).
- The [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made) picks the cheapest verified circuit per task from the dataset. It is re-running on today's last additions, at about 23.136.
- [VQCC Quantum Voodoo](https://www.kaggle.com/code/manwithacat/vqcc-quantum-voodoo), with its [dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-private-dqi-candidates): animations of circuit surgery, measurement-based uncomputation, the DQI draw and the failure fix. Its version history dates each step.
- Coming soon in the open kickmix-macos fork: a headless renderer that makes these width, diff, trace and cost animations from any .kmb file in a single command, written for people and their agents.

## 6. What we think it means (opinion)

The public frontier moved more today than in any day so far, and almost all of it came from people building on each other's Apache-2.0 files within hours. The leaders are still about 16% ahead on score (×1.55 on the product of the three task ratios), so the large private gains are somewhere we haven't seen yet. For DQI, judge circuits on expected failures over many draws, not the draw you happened to get: the final scoring draws again.
