# RezOps Briefing: chemistry falls 2.8% via Atakan Aldemir's pi gauge and Syed Asad Ali's C = 44 table; invsqrt under 97,810

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756210 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-11-05 17:19 UTC · 2 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps runs the background work of
our research programme: it scans new public notebooks, re-verifies their circuits, maintains our public circuits
dataset and the notebook that blends the best of it, and drafts this briefing. A human (me) sets the direction,
reviews and approves every post. Numbers are labelled: *observed* (leaderboard), *author-reported* (the author's
official-verifier run), or *re-verified by us* (the unchanged official verifier, run locally at its default salt).
DQI circuits are ranked by their mean cost over 48 fixed draws (*our estimate*), not by one draw. Window: 04 Nov 17:00
UTC to 05 Nov 17:00 UTC; leaderboard tables from the 16:35 UTC snapshot. Corrections welcome.

## 1. Leaderboard (05 Nov, 16:35 UTC, observed)

| | 04 Nov 16:25 UTC | 05 Nov, 16:35 UTC |
| --- | --- | --- |
| Top score | 28.83206 | **29.05916** (James & Geremie & Tony) |
| Teams ≥ 26 | 16 | 24 |
| Teams ≥ 25 | 36 | 59 |
| Teams ≥ 24 | 847 | 1,029 |
| Teams ≥ 23 | 1,274 | 1,399 |
| Teams ≥ 22.5 | 1,487 | 1,594 |
| Teams ≥ 22 | 1,521 | 1,626 |
| Teams ≥ 20 | 1,614 | 1,720 |
| Teams ≥ 15 | 1,682 | 1,786 |
| Median | 23.46148 | 23.81128 |
| Teams | 2,286 | 2,384 |

The top ten are 26.71–29.06 (James & Geremie & Tony, Schrödinger’s cat lived, Jack, Breaking1800, Thinking..., Yan
Zhou, Bingran You, taiseiu, 田村さんを囲む会, Diana). The median team is at 23.81.

**Approximate medal lines today** (Kaggle's rules for 1,000+ teams: gold top 10 + 0.2%, silver top 5%, bronze top
10%, percentages rounded down; see Kaggle's published [medal rules](https://www.kaggle.com/progression/competitions);
on 2,384 teams):

| | Rank needed | Score at that rank now | Previous | Teams within one DQI failure (0.153%) below |
| --- | --- | --- | --- | --- |
| Gold | ≤ 14 | 26.35447 | 26.03285 | 1 |
| Silver | ≤ 119 | 24.81410 | 24.59425 | 299 |
| Bronze | ≤ 238 | 24.79796 | 24.56271 | 216 |

Ties at the frontier (observed): 24.80731 held by 7 teams (ranks 125-131); 24.80295 held by 6 teams (ranks 144-149);
24.80243 held by 6 teams (ranks 157-162); 24.80233 held by 14 teams (ranks 163-176); 24.80186 held by 11 teams (ranks
181-191); 24.80176 held by 9 teams (ranks 193-201); 24.80146 held by 5 teams (ranks 209-213); 24.80144 held by 3 teams
(ranks 214-216). Ties are broken by the earlier submission
([rules](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/rules)).

The top moved more than on any recent day: James & Geremie & Tony rose 0.79% (best submission 05 Nov 01:24 UTC) and
now lead Schrödinger's cat lived by 1.26%. The gold line rose 1.2% and eight more teams are above 26. Below that, the
public plateau rose about 0.9% in a day, and it carried both medal lines with it. Silver (24.81410) and bronze
(24.79796) are now only 0.065% apart, with about 120 teams packed between them, mostly in the tie groups above. The bronze
line sits inside the cluster of public notebooks that list 24.79–24.80 in their titles or text (malachi w moody's
24.79280, maoxin123's 24.798, Kaliden Aiarys' 24.80179, alex chilton's 24.80; all *author-reported*). The silver line
was above every public notebook score we had seen before 16:25 UTC. We have not matched the tie lines to files.

## 2. Public frontier by task (re-verified by us unless stated)

| Task | Best public | Qubits | Cost | Ratio to reference | Source | Change since last briefing |
| --- | --- | --- | --- | --- | --- | --- |
| invsqrt | `19a0386ae02f`: Syed Asad Ali's rewrites on the public invsqrt lineage (mrmeows → busyaprime → yoyobar) | 62 | 97,809 | 21.23 | [Syed Asad Ali, *VQCC Own-Circuit Portfolio Notebook*](https://www.kaggle.com/code/syed3000/vqcc-own-circuit-portfolio-notebook) (version 4) | was `0104d1b44d73` at 97,922 (−0.11%) |
| dqi | `1c43b7ff2fcf`: Alex Chilton's 04 Nov version of the 891-qubit decoder (Syed Asad Ali's line), released by Leo Tang | 891 | 3.7566 × 10^11 (48-draw mean) | 65.06 | [Leo Tang, *VQCC T7421 Submission*](https://www.kaggle.com/code/leotang09/vqcc-t7421-submission) (version 17) | was `fc3ba81614e1` at 3.7570 × 10^11 (−0.01%) |
| chemistry | `eff6a4863f08`: Syed Asad Ali's own (8, 143, 44) table, loaders and weights; λ_eff 46.8937, 15,120 Toffolis | 427 | 3.0276 × 10^8 | 11.08 | [Syed Asad Ali, *VQCC Own-Circuit Portfolio Notebook*](https://www.kaggle.com/code/syed3000/vqcc-own-circuit-portfolio-notebook) (version 4) | was `a56ec354e8b8` at 3.1147 × 10^8 (−2.80%) |

Best-per-task geometric mean: **24.83106** (*static model / projection*: our arithmetic on the re-verified costs, DQI
at its 48-draw mean; no team has scored it); at the last briefing 24.58709 (+0.99%). All three frontier sources are
Apache-2.0 (notebook footers read today). Syed Asad Ali's notebook reports a local geometric mean of 24.8373
(*author-reported*).

## 3. What's new

We added 33 circuits since the last briefing, each re-verified and credited. Chemistry moved most, in three steps:
a sign gauge with new loaders, a run of weight refits on it, then a new table. In order:

**Chemistry**
- **Evening of 04 Nov, on Ayodeji's loaders:** Isabe11aM00re
  ([*QCC v66s33w submission*](https://www.kaggle.com/code/isabe11am00re/qcc-v66s33w-submission)) rebuilt loader 10
  (311.43M); then qwertykeypadapple's
  [*QCC Kernel Search 24.65*](https://www.kaggle.com/code/qwertykeypadapple/qcc-kernel-search-24-65) cut 110 Toffolis
  from the loaders, for 309.21M. [Agniv Banerjee](https://www.kaggle.com/datasets/somejoint/qcc-ladder-parts),
  [Kaliden Aiarys](https://www.kaggle.com/code/kalidenaiarys/quantum-circuit-compilation-public-lb-24-59425) and
  [Kirill Tymchishin](https://www.kaggle.com/code/kirilltym/quantum-compilation-19) published ties of these.
- **Syed Asad Ali's own (8, 143, 45) table** (v3, 21:28 UTC): λ_eff 47.020 at 15,321 Toffolis, 307.61M (−0.52%).
  Kirill Tymchishin (*Quantum Compilation 19*, 00:12 UTC) took it to 15,277 Toffolis.
- **Soft words** (01:09 UTC): *i want to run gemma 9b* published a new word table of the same shape, with new
  loaders for all 15 slots and new weights: 15,175 Toffolis, λ_eff 47.171, **305.65M** (−0.64%). Their write-up,
  [*VQCC DF-THC loaders: what worked and what didn't*](https://www.kaggle.com/code/iwanttorungemma9b/vqcc-df-thc-loaders-what-worked-and-what-didn-t),
  goes with the [dataset](https://www.kaggle.com/datasets/iwanttorungemma9b/vqcc-dfthc-soft-word-loaders).
  Akmal Xodarev ([busyaprime](https://www.kaggle.com/code/busyaprime/lb-24-73387-top-30-qcc-notebooks-re-verified),
  v18) refit Syed's weights to λ_eff 46.944.
- **Atakan Aldemir's pi gauge** ([*QCC: Gap-3 Residue Repair and C43 Native Compiler*](https://www.kaggle.com/code/atakanaldemir/qcc-gap-3-residue-repair-and-c43-native-compiler),
  run at 04:10 UTC). Atakan re-parameterised Syed's table with a new Givens tree in which some angles carry an extra
  π (section 5 explains why that is free), and built 15 loaders for it: 15,181 Toffolis, down from 15,277. Kirill
  Tymchishin integrated it in *Quantum Compilation 19*. A relay of weight refits on these loaders followed, all at 427
  qubits and 15,181 Toffolis:
  - [maoxin123](https://www.kaggle.com/code/maoxin123/lb-24-786-pi-gauge-loaders-transferred-weights) (v2, 05:48 UTC)
    moved busyaprime's weights across, 304.31M;
  - [alex chilton](https://www.kaggle.com/code/alexchilton/24-79-lb-verified-three-task-portfolio) (v109), 304.30M;
  - [Heliosli](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy) (v10, 08:35 UTC) refined
    busyaprime's refit to λ_eff 46.876, **303.87M** (−0.58% on the soft words);
  - [malachi w moody](https://www.kaggle.com/code/malachiwmoody/vqcc-lb-24-7928-public-parts-our-refits) (v9) refit
    on Syed's table and carried the weights across the gauge, 304.18M;
  - [maoxin123](https://www.kaggle.com/code/maoxin123/lb-24-798-pi-gauge-loaders-refit-chem-weights) (v4, 11:48 UTC),
    303.86M, a tie with Heliosli's within 0.002%.
- **Syed Asad Ali's (8, 143, 44) table** (v4, 16:25 UTC). This is one weight column fewer, from their own continuous
  fit and lattice-aware rounding, with their own loader compiler at 4,229 Toffolis per pass. The notebook says the
  compiler also uses the π-gauge freedom of each address (bit 7) to carry helper products from one batch to the next.
  λ_eff 46.894, 15,120 Toffolis, 427 qubits: **302,756,737** (−0.37% on Heliosli's). This is today's chemistry
  record, −2.80% on yesterday's. Syed credits Heliosli with first showing publicly that weights can be refit right up
  to the error gate.
- [Dmitrii Gluzdov](https://www.kaggle.com/code/dmitriigluzdov/quantum-circuit-compilation-from-source) (v52, v54),
  [Pablo Diego Pereira](https://www.kaggle.com/code/pablodiegoo/google-quantum-compilation-starter-notebook),
  [qiuqiu](https://www.kaggle.com/code/qiuqiuh/qcc-refresh-v70z-k5-ay3-bound-dqi) and
  [matousakura](https://www.kaggle.com/code/matousakura/best-of-public-portfolio-vqcc-24-0x-class) also released
  chemistry in the window.

**invsqrt**
- [malachi w moody](https://www.kaggle.com/code/malachiwmoody/vqcc-lb-24-7928-public-parts-our-refits) (v9, 09:15 UTC)
  made two exact edits to Kirill Tymchishin's `0104d1b44d73`. One replaces a CCZ that fires on half the shots with
  one- and two-qubit phase gates. The other replaces a 19-step fix-up ladder that matters only for x = 1, 2 and 3 with
  a smaller block. Together they save 0.75 Toffoli-equivalents per shot (*author-reported*): 97,875.13 (−0.048%).
- [yoyobar](https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits) (v31) ported a conditional phase
  rewrite onto the current circuit, 97,921.23 (−0.0005%).
- [Syed Asad Ali](https://www.kaggle.com/code/syed3000/vqcc-own-circuit-portfolio-notebook) (v4, 16:25 UTC): **97,809.27**
  (−0.067% on malachi's), at 62 qubits. The notebook lists its rewrites: exact window resynthesis, a measurement
  clean-up, and several further proven rewrites, each checked on all 3,133,452 inputs and all measurement outcomes
  (*author-reported*). Their earlier v3 (98,208) is also in the dataset.

**DQI**
- Akmal Xodarev ([busyaprime, *33 QCC V8 members*](https://www.kaggle.com/datasets/busyaprime/33-qcc-v8-members), 23:29
  UTC) published an exact rewrite of Syed Asad Ali's `fc3ba81614e1` (−0.004% on the 48-draw mean). yoyobar's
  product-cache rewrite `e87b2de64313` ties it.
- [Leo Tang](https://www.kaggle.com/code/leotang09/vqcc-t7421-submission) (v17, 04:45 UTC) released Alex Chilton's
  04 Nov version of the 891-qubit decoder, `1c43b7ff2fcf`, crediting Alex Chilton for it: 3.7566 × 10^11 on the
  48-draw mean (−0.007% on busyaprime's). This is the new DQI frontier.
- Dmitrii Gluzdov's v54 builds an 892-qubit decoder from source, at 3.79 × 10^11 on the 48-draw mean (+1.0%).

**Worth reading.** malachi w moody's notebook has a "what did not work" table with a written kill line for each idea.
It also has a clear account of why a sign gauge cannot change the cost terms, and of how a DQI file's draw is tied to
its own hash. Syed Asad Ali's notebook also lists what did not work (7-bit angles, other loader address splits).

All of these are in the [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits)
(`manifest.csv`, with what changed against the previous best) and in the
[blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made)'s new-circuits table.

## 4. From the forum

A quieter day on the forum. We found no host answer in the window to any of the open rules questions.

- **R resets in place of uncomputation**
  ([755302](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755302)).
  Satvik Verma asked the hosts again for a ruling: is resetting input-dependent scratch with R allowed if the
  verifier passes? This, the generator question
  ([755488](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755488)) and
  the off-path loader question
  ([755427](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755427)) are
  still open.
- **Field notes from a crowd-built entry**
  ([755972](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755972),
  Lynn Langit). This is a visual summary of two weeks of re-verifying and blending: what each task charges, eight
  lessons, and ten closed levers.
  - The most useful lesson for everyone: gate deletions in invsqrt passed every input at the default seed, then failed
    the grader with a half-turn phase error. The deleted gates were phase fixes for rare measurement outcomes.
  - Lynn's advice is to check any edit that touches measurement-conditioned gates on several salts. We agree; it is
    the same reason we rank DQI files over 48 draws.
- **Prize logistics**
  ([755996](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755996),
  radidani). This asks whether a winner in Japan can complete verification and payment remotely and in writing. No
  host answer yet.
- **How do you work?**
  ([755953](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755953)). An
  entrant who does not code reports reaching 24.7 by asking an agent to do its best, and asks how others work.
- Pablo Diego Pereira updated the grid-search versus exact-angles post with better numbers
  ([755525](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755525)).

We added four author-reported nulls to our negative-results ledger: Lynn Langit's invsqrt gate deletion, and three
from malachi w moody's notebook (a (4, 143, 109) chemistry table, a pruned (9, 124, 39) table that is parked, and four
small DQI savings). They will appear in the next version of the
[negative-results dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-chemistry-negative-results).

## 4b. Position risk: the final DQI draw, circuit by circuit

| DQI file (SHA-256 prefix) | Qubits | Board failures | P(final worse) | P(final better) | Expected change | One more failure |
| --- | --- | --- | --- | --- | --- | --- |
| `4ee7a17b4cbc` | 889 | 0 | 0.20 | 0.00 | -0.033% | -0.153% |
| `b3354f60e459` | 897 | 2 | 0.00 | 0.96 | +0.233% | -0.116% |
| `7beea291a5de` | 889 | 0 | 0.63 | 0.00 | -0.142% | -0.153% |
| `f3d335c3f693` | 877 | 0 | 0.67 | 0.00 | -0.158% | -0.153% |
| `ef334fee63f3` | 883 | 0 | 0.28 | 0.00 | -0.050% | -0.153% |
| `0161fd2a6e50` | 889 | 1 | 0.31 | 0.33 | -0.005% | -0.129% |

The table above covers older files only. The new frontier file `1c43b7ff2fcf` is not in it yet, because we do not
know how many failures it drew on the public leaderboard.
- **What we know (our estimate).** It passes all 48 of our fixed draws, and 9 of them have one failed shot each (8
  value errors, 1 phase error). So the chance of at least one failure on a new draw is about 1 in 5, and each failure
  costs about 0.15% of the overall score.
- **Another entrant's check.** malachi w moody reports the same picture: one value error on one of four local salts
  (*author-reported*).
- **What this means for the medal lines.** The silver and bronze lines are now 0.065% apart, so one failed DQI shot
  in the final draw is worth more than twice the gap between them. The 299 and 216 teams in the section 1 table sit
  within that distance of a line.

## 4c. Prize watch (leaderboard-timed prizes)

**Prize watch** (from 105 leaderboard snapshots; latest 05 Nov 16:35 UTC; our estimate)

- **Leader:** James & Geremie & Tony at 29.05916 (best submission 05 Nov 01:24 UTC); next is Schrödinger’s cat lived
  at 28.69687, 1.26% behind.
- **Mid-Competition Milestone Prize** (highest public score at 20 Nov 23:59 UTC, notebook published): 15.3 days left.
- **Longest Leader Prize** (most time at #1 by 4 Jan 23:59 UTC): James & Geremie & Tony has led since between 02 Nov
  20:26 and 02 Nov 21:25 UTC, i.e. 2.8–2.9 days. If it stays #1, no other team can catch its total after **10 Dec
  08:26–10 Dec 09:25 UTC** (start of its stint plus half the 75-day window; worst case, eligibility not checked).

## 5. Tools and notes from us

- **Circuits dataset and blend notebook.** The
  [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits) is now v105, with 345
  credited circuits in `manifest.csv`. The
  [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made) re-ran after each
  record and is now v116.
- **Hedge Wizardry** ([notebook](https://www.kaggle.com/code/manwithacat/vqcc-hedge-wizardry-notebook) v24 and v25,
  dataset v23 and v24). Three small, credited entries, each run through the unmodified official verifier on Kaggle
  (*official verifier PASS, hosted*):
  - **HW-21 (chemistry).** A weight refit of Heliosli's `aba448aeadb5` at fixed words and loaders: 303,854,940.78
    (−0.004%).
  - **HW-22 (DQI).** We took Akmal Xodarev's `fd0ef09efd45` and Alex Chilton's `1c43b7ff2fcf` (as released by Leo
    Tang). Where the two files are byte-identical, we joined the unsynthesis section of the first to the second. The
    result has 14 fewer static Toffolis at 891 qubits, and a 48-draw mean of 375,643,651,078.69 (−0.004%, our
    estimate). It passed on Kaggle with 0 value and 0 phase errors.
  - **HW-23 (chemistry).** HW-21's weights with four loaders re-planned at 276 Toffolis (from 277; 8 fewer
    Toffolis in all) and emitted in Atakan Aldemir's pi gauge: 303,694,816.97 (−0.05%). It was superseded within hours by Syed Asad Ali's new table.
- **Why the pi gauge is free.** Each vector u in the chemistry table is stored as a chain of Givens angles, and the
  cost model only ever uses u through u uᵀ. Adding π to selected angles flips signs inside the chain without changing
  u uᵀ, so it is an exact symmetry: λ_eff, both errors and the gap stay the same to rounding (we measured a largest
  projector change of 6.7 × 10^-16).
  - In the gauged table that Atakan Aldemir's loaders use, 20.6% of the angle words sit off the plain signed convention. The rows fall into 8 flip
    patterns.
  - A loader can therefore let a word's top bit carry the pattern and choose the cheaper form of each word. That is
    where Atakan's Toffoli saving comes from, and why weights refit on Syed's plain table carry across exactly, as
    malachi w moody showed.
  - Syed Asad Ali's new compiler now uses the same freedom (bit 7) for its own purpose.
- **The fast verifier in use.** Syed Asad Ali's notebook reports verifying with our Rust
  [kickmix-engine](https://github.com/manwithacat/kickmix-engine), and Lynn Langit's field notes mention it too. We are
  glad it is useful. As always, the unmodified official verifier on Kaggle is the authority.

## 6. What we think it means (opinion)

- **The public plateau now sets the medal lines.** Silver and bronze moved about 0.9% in a day on public parts alone,
  and are now 0.065% apart. The best public parts sit just above silver by our arithmetic (24.831 against 24.814).
  From here the final DQI draw matters as much as the next public record.
- **Chemistry is now a table-and-loader problem, not a refit problem.**
  - Refits at fixed words are close to exhausted. malachi w moody's restarts gained 0.0017% on the score, and our
    HW-21 gained 0.004%.
  - The two real moves today came from structure: a gauge that the loaders can exploit (Atakan Aldemir), and a smaller
    table built together with its loaders (Syed Asad Ali).
  - We expect the next records to come from tables chosen for their loader cost.
- **The top is pulling away again.** The leader gained 0.8%, and the top ten sit 8–17% above the public frontier. The
  unanswered rules questions in section 4 still matter more for prize entries than any single public circuit.
