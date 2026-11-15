# RezOps Briefing: the 24.94 plateau — chemistry −1.0% on Syed Asad Ali's re-planned loaders; silver and bronze 0.04% apart

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756558 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-11-06 17:56 UTC · 1 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps runs the background work of
our research programme: it scans new public notebooks, re-verifies their circuits, maintains our public circuits
dataset and the notebook that blends the best of it, and drafts this briefing. A human (me) sets the direction,
reviews and approves every post. Numbers are labelled: *observed* (leaderboard), *author-reported* (the author's
official-verifier run), or *re-verified by us* (the unchanged official verifier, run locally at its default salt).
DQI circuits are ranked by their mean cost over 48 fixed draws (*our estimate*), not by one draw. Window: 05 Nov 17:00
UTC to 06 Nov 17:06 UTC. Corrections welcome.

## 1. Leaderboard (06 Nov, 17:06 UTC, observed)

| | 05 Nov 16:35 UTC | 06 Nov, 17:06 UTC |
| --- | --- | --- |
| Top score | 29.05916 | **29.30100** (James & Geremie & Tony) |
| Teams ≥ 26 | 24 | 25 |
| Teams ≥ 25 | 59 | 96 |
| Teams ≥ 24 | 1,029 | 1,157 |
| Teams ≥ 23 | 1,399 | 1,499 |
| Teams ≥ 22.5 | 1,594 | 1,685 |
| Teams ≥ 22 | 1,626 | 1,718 |
| Teams ≥ 20 | 1,720 | 1,807 |
| Teams ≥ 15 | 1,786 | 1,870 |
| Median | 23.81128 | 23.90494 |
| Teams | 2,384 | 2,481 |

The top ten are 26.87–29.30 (James & Geremie & Tony, Schrödinger’s cat lived, Jack, Breaking1800, taiseiu,
Thinking..., Diana, Yan Zhou, Bingran You, CroDoc & Pavel). The median team is at 23.90.

**Approximate medal lines today** (Kaggle's rules for 1,000+ teams: gold top 10 + 0.2%, silver top 5%, bronze top
10%, percentages rounded down; see Kaggle's published [medal rules](https://www.kaggle.com/progression/competitions);
on 2,481 teams):

| | Rank needed | Score at that rank now | Previous | Teams within one DQI failure (0.153%) below |
| --- | --- | --- | --- | --- |
| Gold | ≤ 14 | 26.49520 | 26.35447 | 0 |
| Silver | ≤ 124 | 24.95166 | 24.81410 | 390 |
| Bronze | ≤ 248 | 24.94093 | 24.79796 | 260 |

Ties at the frontier (observed): 24.94111 held by 8 teams (ranks 152-159); 24.94109 held by 14 teams (ranks 160-173);
24.94103 held by 3 teams (ranks 175-177); 24.94101 held by 4 teams (ranks 178-181); 24.94099 held by 27 teams (ranks
182-208); 24.94096 held by 13 teams (ranks 211-223); 24.94093 held by 32 teams (ranks 227-258); 24.94090 held by 83
teams (ranks 259-341). Ties are broken by the earlier submission
([rules](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/rules)).

The public plateau rose about 0.55% in a day, and both medal lines rose with it. Silver (24.95166) and bronze
(24.94093) are now only 0.043% apart. The bronze cut falls inside a tie group of 32 teams at 24.94093. Just below it,
83 teams share 24.94090: 0.0001% below the line, and ranked only by submission time. That value is where today's
public parts land. Public notebooks list 24.94 in their titles: maoxin123's *LB 24.94* and MAIDANG's *LB 24.94096*
(*author-reported*). Our blend notebook computes 24.94088 locally (default salt). The top also moved:
James & Geremie & Tony rose 0.83% (best submission 06 Nov 00:25 UTC) and now lead Schrödinger's cat lived by 2.11%.

## 2. Public frontier by task (re-verified by us unless stated)

| Task | Best public | Qubits | Cost | Ratio to reference | Source | Change since last briefing |
| --- | --- | --- | --- | --- | --- | --- |
| invsqrt | `c27e83acb576`: the public 62-qubit lineage with one fewer Toffoli (by our diff against `19a0386ae02f`: −1 static Toffoli, 688 fewer other operations) | 62 | 97,747 | 21.25 | [Andrew Oscar, *QCC 24.86 fork*](https://www.kaggle.com/code/andrewoscar/qcc-24-86-fork) (version 1) | was `19a0386ae02f` at 97,809 (−0.06%) |
| dqi | `7b6cbb92a797`: Atakan Aldemir's 17-product reuse on the 892-qubit decoder of Syed Asad Ali's line (by our diff, 17 fewer static Toffolis than `1ab3ac53d442` in one stage) | 892 | 3.7508 × 10^11 (48-draw mean) | 65.16 | [Atakan Aldemir, *QCC: 17-Product DQI Reuse at 892 Qubits*](https://www.kaggle.com/code/atakanaldemir/qcc-17-product-dqi-reuse-at-892-qubits) (version 3) | was `1c43b7ff2fcf` at 3.7566 × 10^11 (−0.15%) |
| chemistry | `da4d23796b13`: Syed Asad Ali's (8, 143, 44) table with loaders 4–14 re-planned from 202 to 195 products each; 427 qubits | 427 | 2.9967 × 10^8 | 11.20 | [Syed Asad Ali, *VQCC Own-Circuit Portfolio Notebook*](https://www.kaggle.com/code/syed3000/vqcc-own-circuit-portfolio-notebook) (version 9) | was `eff6a4863f08` at 3.0276 × 10^8 (−1.02%) |

Best-per-task geometric mean: **24.93402** (*static model / projection*: our arithmetic on the re-verified costs, DQI
at its 48-draw mean; no team has scored it); at the last briefing 24.83106 (+0.41%). All three frontier sources are
Apache-2.0 (notebook footers read today).

Three notes on credit:
- **invsqrt:** the same file also appears in Isabe11aM00re's dataset
  [*qcc-v67a11-circuits*](https://www.kaggle.com/datasets/isabe11am00re/qcc-v67a11-circuits) (version 1, 05 Nov 22:05
  UTC) and in Syed Asad Ali's portfolio dataset. We have not yet settled who released it first; corrections welcome.
  maoxin123's `c17971aff984` ties it.
- **chemistry:** Heliosli's `1f0eb655dd13`
  ([*The dataset is public now. Enjoy!*](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy), v14)
  uses the same 15 loader files with refit weights and is 2,052 cheaper (299,671,051). That is within our 10^-5 tie
  rule, so the earlier file keeps the slot.
- **DQI:** the 892-qubit decoder was one qubit wider than yesterday's 891 but cheaper on the 48-draw mean. We first
  saw it in [Andrew Oscar's fork](https://www.kaggle.com/code/andrewoscar/qcc-24-86-fork) as `1ab3ac53d442`.

## 3. What's new

We added 21 circuits since the last briefing, each re-verified and credited.

**DQI: a 892-qubit line, then product reuse on it.**
- `1ab3ac53d442` (Syed Asad Ali's line, seen first in Andrew Oscar's *QCC 24.86 fork* and in Syed's portfolio
  dataset) took the slot from `1c43b7ff2fcf` at −0.15%.
- Atakan Aldemir's [*17-Product DQI Reuse*](https://www.kaggle.com/code/atakanaldemir/qcc-17-product-dqi-reuse-at-892-qubits)
  (v3) removed 17 more Toffolis by reusing products, for −0.004%.
- MAIDANG's [*LB 24.94096 C44 and Exact DQI Product Reuse*](https://www.kaggle.com/code/daifanhao/lb-24-94096-c44-and-exact-dqi-product-reuse)
  (`b4bbf86d504e`) and yoyobar (`768495b4ce0f`) published ties.
- milos's [*QCC 24.80 | Best per Task + Where the Cost Goes*](https://www.kaggle.com/code/milosg/qcc-24-80-best-per-task-where-the-cost-goes)
  adds an 891-qubit variant and a cost breakdown per task that is worth reading.

**Chemistry: re-planned loaders, then a weight refit, then many ties.**
- Syed Asad Ali (v9) kept the (8, 143, 44) table and re-planned loaders 4–14: 195 products each (from 202), with a
  partial low cache, at 270–274 Toffolis per loader. **299,673,103** (−1.02%).
- Heliosli refit the weights on those loaders: 299,671,051 (a tie under our rule).
- Ties within 0.001% followed from maoxin123, Atakan Aldemir, MAIDANG, Isabe11aM00re and Pablo Diego Pereira
  (whose [*Google Quantum Compilation Starter Notebook*](https://www.kaggle.com/code/pablodiegoo/google-quantum-compilation-starter-notebook)
  ships a byte-for-byte copy of Syed's family).
- Atakan Aldemir's gap-3 notebook and Andrew Oscar's fork carry 302.20M variants; Or4k2!'s
  [*Quantum circuit compilation pipeline*](https://www.kaggle.com/code/orecord/quantum-circuit-compilation-pipeline)
  carries two byte variants of yesterday's 302.76M family.

## 4. From the forum

- **Does clearing scratch with R pass?**
  ([755302](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755302)).
  Jyoti Negi asked whether resetting input-dependent scratch with R, instead of uncomputing it, is allowed; Satvik
  Verma, Hassan Mohamed and David Pilný asked for a host ruling before 20 Nov.
  - Syed Asad Ali posted a test: R in kickmix randomises the phase when the qubit is |1>, and the verifier runs 64
    seeded shots per x. A passing invsqrt circuit with `CX q0 q61; R q61` appended fails ("x=5: phase error 5.000e-01
    turns"); the same circuit with `CX q0 q61; CX q0 q61` appended passes (*author-reported*).
  - So in practice the published verifier already rejects R on input-dependent values. A host ruling is still
    outstanding.
- **Exact finite-field multiplication counts**
  ([756476](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756476)).
  田村 posted computer-assisted proofs that μ₄(5) = μ₄^sym(5) = 11 and μ₄^sym(6) = 14, with code in a public
  notebook. pjmathematician replied with an independent proof of the asymmetric case, μ₄(6) = 14: no 13-term
  GF(4)-bilinear multiplier for GF(4096) exists at all (*author-reported*; we have not checked the proofs). This
  bounds what any GF(4)-bilinear multiplier for GF(2¹²), the DQI decoder's field, can save.

## 4b. Position risk: the final DQI draw, circuit by circuit

| DQI file (SHA-256 prefix) | Qubits | Board failures | P(final worse) | P(final better) | Expected change | One more failure |
| --- | --- | --- | --- | --- | --- | --- |
| `4ee7a17b4cbc` | 889 | 0 | 0.20 | 0.00 | -0.033% | -0.153% |
| `b3354f60e459` | 897 | 2 | 0.00 | 0.96 | +0.233% | -0.116% |
| `7beea291a5de` | 889 | 0 | 0.63 | 0.00 | -0.142% | -0.153% |
| `f3d335c3f693` | 877 | 0 | 0.67 | 0.00 | -0.158% | -0.153% |
| `ef334fee63f3` | 883 | 0 | 0.28 | 0.00 | -0.050% | -0.153% |
| `0161fd2a6e50` | 889 | 1 | 0.31 | 0.33 | -0.005% | -0.129% |

The table covers older files only. The new frontier file `7b6cbb92a797` is not in it, because we do not know how
many failures it drew on the public leaderboard.
- **What we know (our estimate).** It passes all 48 of our fixed draws, and 9 of them have one failed shot each (8
  value errors, 1 phase error): the same picture as yesterday's frontier file. So the chance of at least one failure
  on a new draw is about 1 in 5, and each failure costs about 0.15% of the overall score.
- **What this means for the medal lines.** Silver and bronze are now 0.043% apart, so one failed DQI shot in the
  final draw is worth more than three times the gap between them.

## 4c. Prize watch (leaderboard-timed prizes)

**Prize watch** (from 108 leaderboard snapshots; latest 06 Nov 17:06 UTC; our estimate)

- **Leader:** James & Geremie & Tony at 29.30100 (best submission 06 Nov 00:25 UTC); next is Schrödinger’s cat lived
  at 28.69687, 2.11% behind.
- **Mid-Competition Milestone Prize** (highest public score at 20 Nov 23:59 UTC, notebook published): 14.3 days left.
- **Longest Leader Prize** (most time at #1 by 4 Jan 23:59 UTC): James & Geremie & Tony has led since between 02 Nov
  20:26 and 02 Nov 21:25 UTC, i.e. 3.8–3.9 days. If it stays #1, no other team can catch its total after **10 Dec
  08:26–10 Dec 09:25 UTC** (start of its stint plus half the 75-day window; worst case, eligibility not checked).

## 5. Tools and notes from us

- **Circuits dataset and blend notebook.** The
  [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits) is now v107, with 366
  credited circuits in `manifest.csv`. The
  [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made) is now v118. Its
  chemistry chart was rebuilt: with 16 switches between word tables the old side panel had become unreadable, so the
  recent switches now get their own panel, and the list is a table under the chart.
- **Our scanner was missing large notebooks.** Notebooks whose outputs exceed 1 GB, or that have hundreds of files,
  had been skipped on every scan. The scanner now downloads only their circuit files. That found 10 circuits we
  had missed, from Or4k2!, Ali Haydar Özdağ, MAIDANG and Atakan Aldemir; they are now in the dataset with credit.
- **A refit at the constraint edge is not byte-reproducible, and that is expected.** Our
  [dfthc-refit reproduction notebook](https://www.kaggle.com/code/manwithacat/dfthc-refit-reproduce-hw-21-and-hw-23)
  (v2) now reports a finding.
  - A weight refit that holds both residuals right at the 20 mHa caps stops on threshold tests.
  - Last-bit floating-point differences (CPU, library versions, thread count) flip those tests. That changes how
    many evaluations the optimiser gets, and so where it ends.
  - On one Mac, changing only the thread count gave four different weight files; with Kaggle's run, the end points
    spread by 2.1 × 10^-5 relative, all inside the caps.
  - The second Kaggle run reproduced the first byte for byte, so each setup is deterministic; the split is between
    setups.
  - In practice: such refits reproduce to the verifier's pass, not to the byte; ship the weights, and a few thread
    counts give a free multi-start.
- **Hedge Wizardry HW-24** (hosted official-verifier PASS, released this morning) re-planned loaders 4–14 on Syed
  Asad Ali's earlier family for −0.29%. Syed's own re-plan above superseded it within hours.

## 6. What we think it means (opinion)

- **The bronze line is now a submission-time race.** 83 teams sit 0.0001% below it, all on public parts. A team at
  the plateau can only clear the line by being earlier, or by finding something the public parts do not have.
- **Chemistry loaders are near their floor.** Syed Asad Ali's re-plan put every coded loader below the best our own
  planner reaches. Weight refits on the new loaders now gain about 10^-7 relative. The next real step probably needs
  a new table, as on 05 Nov, not more refitting.
- **DQI gains are now single products.** The day's DQI records were −0.15% for a new decoder and −0.004% for 17
  reused products. The new field-multiplication bounds in 756476 suggest the multipliers are close to optimal too.


---

## Comments (3 threads)

Replies are indented with `>`, one level per reply depth.

**Navneet** · 2026-11-07 04:17 UTC · -1 votes

Thank you for the RezOps Briefing @manwithacat

**Zero** · 2026-11-07 03:53 UTC · 0 votes

Hi

My agent and I found a potential issue in the third-track FeMoCo circuits.

In the ANGLE_1 circuit from the yoyobar v30 dataset, our independent test found that 761 out of 1,220 addresses failed the documented XOR lookup behavior when starting with a zero payload. However, the circuit still passed the forward and reverse sequences tested by the official verifier.

You should check whether the circuits in your collection satisfy the XOR lookup contract for arbitrary initial payloads, rather than just the sequences covered by the official verifier.

It's also worth checking whether kickmix-engine independently validates this behavior, since relying on the official verifier alone could reproduce the same coverage gap.

> **manwithacat** — **Topic author** · 2026-11-07 05:44 UTC · 0 votes
>
> Thanks Zero, we reproduce this: 761/1,220 on ANGLE_1 for the current top chemistry tables too, with workspace clean and no phase. It's expected, though. These loaders use the fact that the word registers already hold the previous batch, and the hosts ruled in discussion 753092 that ANGLE_k only has to be correct on states reachable in the forward and reverse traversals, which is what the verifier checks. A plain QROM family passes the zero-payload test at every k, so the test itself is a good way to tell strict loaders from reachable-only ones. For coverage beyond the verifier's 64 sampled shots, --contract reachable mode in our ANGLE_k checker proves the reachable contract for all measurement outcomes.

> > **Zero** · 2026-11-07 06:39 UTC · 0 votes
> >
> > Thanks for the explanation!

**(deleted user)** · 2026-11-07 11:37 UTC · 0 votes

*(deleted)*
