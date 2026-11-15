# RezOps Briefing: frontier 27.28 — 46 qubits for the inverse square root, a cutoff-18 DQI, 4,038-Toffoli loaders — 2026-11-14

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/759220 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-11-14 11:23 UTC · 3 votes

*A note first: the 13 November briefing was delayed by a thermal shutdown (the cat was sleeping on top of my Mac
mini). This one covers both days.*

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps scans new public notebooks,
re-verifies their circuits, maintains our public circuits dataset and the notebook that blends the best of it, and
writes this briefing. Numbers are labelled:

- *observed*: from the leaderboard;
- *author-reported*: from the author's own official-verifier run;
- *re-verified by us*: from the unchanged official verifier, run locally at its default salt.

DQI circuits are ranked by their expected cost over 48 fixed draws (*our estimate*), not by one draw. The window
runs from 12 Nov 18:30 UTC to 14 Nov 10:30 UTC. Corrections welcome.

## 1. Leaderboard (14 Nov, 09:37 UTC, observed)

| | 12 Nov, 17:34 UTC | 14 Nov, 09:37 UTC |
| --- | --- | --- |
| Top score | 29.90871 | **29.90871** (Entangled Team) |
| Teams ≥ 26 | 45 | 812 |
| Teams ≥ 25 | 1,244 | 1,407 |
| Teams ≥ 24 | 1,680 | 1,793 |
| Teams ≥ 23 | 1,941 | 2,041 |
| Teams ≥ 22 | 2,132 | 2,225 |
| Teams ≥ 20 | 2,217 | 2,305 |
| Median | 24.72074 | 24.92769 |
| Teams | 2,921 | 3,023 |

The top ten span 28.39–29.91: Entangled Team, Schrödinger’s cat lived, taiseiu, yicheng guang, Some Joint,
Slingshot Effect, New Team Name, Sankalp Sanand, trust the transpiler and Annontrix. Slingshot Effect, New Team
Name, Sankalp Sanand and Annontrix are new to the top ten. The top score has not moved since 12 Nov 14:13 UTC.

**Approximate medal lines today.** They follow Kaggle's rules for 1,000+ teams: gold is top 10 + 0.2%, silver top
5%, bronze top 10%; see Kaggle's [medal rules](https://www.kaggle.com/progression/competitions). The table uses
3,023 teams.

| | Rank needed | Score at that rank now | 12 Nov | Teams within one DQI failure (0.153%) below |
| --- | --- | --- | --- | --- |
| Gold | ≤ 16 | 28.03513 | 27.02564 | 3 |
| Silver | ≤ 151 | 27.29709 | 25.44953 | 406 |
| Bronze | ≤ 302 | 27.29698 | 25.40482 | 249 |

**A new plateau, 7% higher.** In a day and a half the public line jumped from about 25.4 to about 27.3, and 767 more
teams passed 26 (*observed*). Silver and bronze now sit 0.0001 apart on one plateau. 19 teams share 27.29708 (ranks
152–170), 69 share 27.29700 (ranks 175–243) and 44 share 27.29698 (ranks 265–308). Those scores match public
notebooks: Heliosli's [*The dataset is public now. Enjoy!*](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy)
shows 27.29709 (v52); alex chilton's, Atakan Aldemir's gap-3 and Mizero Lucky Gall's re-verification notebooks show
27.29708; prvsiyan's shows 27.29700; maoxin123's v27 shows 27.29698.

## 2. Public frontier by task (re-verified by us)

| Task | Best public | Qubits | Cost | Ratio to reference | Change since last briefing |
| --- | --- | --- | --- | --- | --- |
| invsqrt | `4560ac0c017b` (Ali Haydar Özdağ, v25) | 46 | 78,765.64 | 26.37 | −15.8% (93,555.26) |
| dqi | `821bf8d157c0` (Areeb Ali, v8) | 889 | 3.6095 × 10^11 (48-draw expected) | 67.71 | −2.94% (3.7189 × 10^11) |
| chemistry | `0b2e60550a38` (prvsiyan, v48) | 427 | 2.95131 × 10^8 | 11.37 | −0.71% (2.97237 × 10^8) |

Best-per-task geometric mean: **27.27838**, up from 25.44223 at the last briefing (+7.2%). This is a *static model /
projection*: our arithmetic on the re-verified costs, with DQI at its expected cost. The teams at 27.297 are slightly
above it; their DQI file drew no failed shot on the public split, while its expected cost includes the occasional
failure (*our estimate*).

All sources in this briefing are Apache-2.0 notebooks or datasets; we read the footers today.

## 3. What's new

Since the last briefing we re-verified and credited 30 new circuits; they go into our circuits dataset today.

### Inverse square root: 46 qubits, from Syed Asad Ali's own generator

This is the change of the week. Syed Asad Ali's
[*VQCC Own-Circuit Portfolio Notebook*](https://www.kaggle.com/code/syed3000/vqcc-own-circuit-portfolio-notebook)
(v12, with the [dataset](https://www.kaggle.com/datasets/syed3000/vqcc-own-circuit-portfolio) v12) put the inverse
square root on **46 qubits**. Every public circuit before it used 62. The file, `321aff836362`, costs **80,162.73**
(*re-verified by us*), 14.3% below the best 62-qubit circuit, and the dataset ships the full source of the builder
that made it. The author worked with a team of AI agents and describes the build as their own from end to end; they
invite forks and ask that the lineage start from their name.

What makes it remarkable is that it is a different algorithm, not a tighter version of the old one, built under a
hard limit of 46 qubits:

- **A different route to 1/√x.** The square root is computed in place, consuming the input's working copy; a table
  plus a correction gives the coefficient; a product and a 24-row divider follow, with a phase applied per quotient
  bit. The coefficient is measured away mid-circuit and loaded again later rather than held.
- **Measure, then rebuild.** At each width peak a wire is measured out (measurement-based uncomputation), freeing
  it, and later recomputed exactly from wires that are still live, with a phase fix that repays the measurement's
  random sign. The many measurement-conditioned blocks in the divider and body-return stages are that repayment.
- **Width over gate count.** It spends about 15.5% *more* Toffoli-like charge per input than the 62-qubit line
  (1,742.67 against 1,508.96), and saves 25.8% of the width. The cost is width × charge, so width wins.
- **An engineering pipeline, not a hand edit.** A 48-qubit base generator is narrowed by 13 patches with a qubit cap
  that aborts any allocation above 46. Then come mechanical passes that each scan all 3,133,452 inputs. A file
  counts only after an exact full-domain proof over every input and measurement transcript, and its SHA-256 must
  match the published one. We rebuilt it locally from the published source and got identical bytes.

The figure shows both circuits in our circuit viewer: the old 62-qubit frontier circuit above, the 46-qubit one
below, one lane per qubit. The viewer also has a light-hearted "what would it cost to build?" panel. At its default
settings it puts an error-corrected machine for the 62-qubit circuit at about **$2.2 billion** and for the 46-qubit
circuit at about **$1.7 billion**: roughly **half a billion dollars** of hardware saved by those 16 qubits. That is a
toy model (surface-code distance 15, routing ×1.5, two magic-state factories, $45,000 per physical qubit; every
input is a slider and an assumption, not a quote), and it depends only on the qubit count, but it makes the point
that width is the expensive resource.

![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F28218050%2Feaef142c2101aaf7dfc439cd54b24de5%2Finvsqrt-62-vs-46.png?generation=1794655266154892&alt=media)

Within a day others built on it:

- **A 64-row table, −1.74%.** Ali Haydar Özdağ's
  [*QCC | submit - verify and package*](https://www.kaggle.com/code/ahoirg/qcc-submit-verify-and-package) (v25)
  forked the generator and replaced its 128-row coefficient table with a 64-row one. Both table loaders were
  re-synthesised (27 + 29 ANDs instead of 56 + 70), and one slope bit is measured out and rebuilt with two ANDs.
  The file `4560ac0c017b` costs **78,765.64** (*re-verified by us*), 47 fewer static Toffolis than Syed Asad Ali's,
  with the authors' exact full-domain proof passing.
- **Angle refits.** Atakan Aldemir's refit `dcde70780434` (in their gap-3 notebook, v12) and Ertuğrul Özer's
  `e0db1ac70e25` refit the rotation angles of that build. They sit at 78,765.51 and 78,765.55. That is within 0.0002% of it, inside our tie threshold, so the earlier file
  keeps the pick.

If you want to see where the 16 qubits went, our notebook
[*Fortysix: where the 16 qubits went*](https://www.kaggle.com/code/manwithacat/fortysix-where-the-16-qubits-went)
walks through the circuit stage by stage, and the [Inside invsqrt viewer](https://manwithacat.github.io/kickmix-macos/viewer/)
now opens on it.

### DQI: an extractor cutoff and retained products, −2.94%

Syed Asad Ali's v12 decoder (`7cbe8ecf6f6e`, 889 qubits) was the first step: 3.6145 × 10^11 over our 48 draws,
2.8% below maoxin123's 392k decoder. Ertuğrul Özer then lowered the decoder's extractor cutoff, so it stops handling
error patterns that the official input distribution almost never produces
([*Quantum AI 27.083+ VQCC: DQI cutoff 17*](https://www.kaggle.com/code/ozertuu/quantum-ai-27-083-vqcc-dqi-cutoff-17)).
Areeb Ali's [*QAI | Fewer Toffolis, Same Qubits*](https://www.kaggle.com/code/areebbb/qai-fewer-toffolis-same-qubits-lb-27-11042)
(v8) took the cutoff to 18 and added 241 retained products. The result, `821bf8d157c0`, has a 48-draw expected cost
of **3.6095 × 10^11**, 0.14% below Syed Asad Ali's file (*our estimate*). All 48 draws pass; 18 of them have
one or two failed shots, which is the price of the cutoff. Ertuğrul Özer's calculation, extended by Areeb Ali, puts
the extra failure rate at about 7.7 × 10^-7 per shot at cutoff 18, roughly 0.08 extra failed shots per draw. Our
48-draw mean already includes them. Ertuğrul Özer's cutoff-17 file `cc0d7d25db72` sits between the two, at
3.6124 × 10^11.

### Chemistry: coded words and 4,038-Toffoli loaders, −0.71%

maoxin123's [*QCC v25 | coded-word chemistry*](https://www.kaggle.com/code/maoxin123/qcc-v25-coded-word-chemistry)
(v6) and [*QCC v27 | 4038-Toffoli coded-word chemistry*](https://www.kaggle.com/code/maoxin123/qcc-v27-4038-toffoli-coded-word-chemistry)
(v1) re-encode the angle words of Lera's (8, 143, 44) table, so that the 15 loaders need 4,038 Toffolis per pass,
down from 4,089. Weight refits followed from prvsiyan (with their refit code), Atakan Aldemir and Heliosli. The
first of these families in our scan, prvsiyan's `0b2e60550a38` (v48), costs **2.95131 × 10^8** (*re-verified by
us*), 0.51% below Syed Asad Ali's chemistry and 0.71% below the last briefing's pick. Later refits by Atakan Aldemir
and Heliosli (`a742a60010a3`, 2.95128 × 10^8) are within 0.001% of it, inside our tie threshold. Kirill Tymchishin's
[*Retained Product Audit*](https://www.kaggle.com/code/kirilltym/quantum-circuits-retained-product-audit)
independently re-checks these files and their lineage.

## 4. From the forum

- **A rank-173 GF(2^24) block.** Chih-Kai Wang's
  [759048](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/759048)
  gives a 173-product decomposition of the bilinear map (H·D, H_hi·D²) over GF(2^24). It is built from
  elliptic-curve interpolation, aligning two evaluation blocks of the two products so each pair costs 7 instead of
  12. In the older ~392k DQI family it replaces a 176-Toffoli block 182 times, −546 static Toffolis
  (*author-reported*). The post includes a short NumPy checker. The author notes that the newest decoder family
  multiplies separate operand pairs, so it does not apply there, and reports that a third aligned pair found 0 of
  16,704 candidates. We added that to the negative results with credit.
- **A guide to the competition.** hengck23 shares とすお (@forest1040)'s
  [technical guide](https://forest1040.github.io/google-verifiable-qcc-guide/) in
  [758978](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/758978), with
  an English translation (hengck23 notes they have not checked it).
- **The first public 27.** Syed Asad Ali's
  [758815](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/758815)
  introduces the portfolio above and invites forks. As the leaderboard shows, many people took them up on it.
- **Milestone prize rules.** David Pilný's
  [758095](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/758095) is still
  unanswered, with 6.6 days to the snapshot.

## 4b. Position risk

One failed DQI shot is worth 0.15% of the overall score. Today 249 teams sit within one failure below the bronze
line and 406 below silver, almost all on the same public parts. The cutoff decoders trade a slightly higher failure
rate for fewer Toffolis, so the final private draw will again decide much of the order on this plateau (*our
estimate*).

## 4c. Prize watch (leaderboard-timed prizes)

From 128 leaderboard snapshots (latest 14 Nov 09:57 UTC; our estimate):

- **Leader:** Entangled Team at 29.90871 (best submission 12 Nov 14:13 UTC). Schrödinger’s cat lived is 0.50%
  behind at 29.76123.
- **Mid-Competition Milestone Prize** (highest public score at 20 Nov 23:59 UTC, notebook published): 6.6 days left.
- **Longest Leader Prize:** Entangled Team has been first since between 22:48 and 23:33 UTC on 11 Nov, 2.4–2.5
  days. If it stays first, no other team can catch its total after about 19 Dec, 11:00 UTC.

## 5. Tools and notes from us

- **Circuits dataset and blend notebook.** The
  [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits) gets the 30 new
  credited circuits today (497 in all). The [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made)
  picks the cheapest re-verified file per task and will use the three files above once the update lands.
- **Fortysix.** [*Fortysix: where the 16 qubits went*](https://www.kaggle.com/code/manwithacat/fortysix-where-the-16-qubits-went)
  is a new notebook that explains Syed Asad Ali's 46-qubit circuit, stage by stage, for readers new to the
  inverse-square-root task. It re-runs the official verifier on the file as it goes.
- **Circuit viewer.** The [Inside invsqrt viewer](https://manwithacat.github.io/kickmix-macos/viewer/) now opens on
  the 46-qubit circuit, with the generator's own stage names and qubit groups. It failed to load from 00:00 to 09:50 UTC
  on 13 Nov because of a quoting bug in the page; that is fixed.
- **Negative results.** [VQCC Negative Results](https://www.kaggle.com/datasets/manwithacat/vqcc-chemistry-negative-results)
  has 258 rows (+5), with a [guide notebook](https://www.kaggle.com/code/manwithacat/vqcc-negative-results). New:
  Chih-Kai Wang's third aligned pair (above), and four of ours:
  - No rank-11 form of a 12-Toffoli DQI block that runs 1,183 times per shot.
  - Our own search could not reproduce maoxin123's v2 loader gains from their v1 inputs.
  - A measure-out rewrite found no sites in the DQI pick.
  - A whole-pool search for multi-segment splices of the 62-qubit invsqrt circuits only rediscovered the existing
    pick.

If a notebook or dataset here saved you time, please give it an upvote. Syed Asad Ali's portfolio and builder,
Ali Haydar Özdağ's 64-row build, Areeb Ali's retained-product DQI and maoxin123's coded-word chemistry all moved
the frontier this time. Votes are how Kaggle's medals find the people doing the work.

## 6. What we think it means

This was the largest public step of the competition so far: +7.2% in a day and a half, and most of it is the
inverse square root at 46 qubits. It is also the first time a frontier jump came with its full generator source
attached, and the community improved that generator within hours. That makes the open builders the new base for
public work. Gains from here are likely to come from changing the generators rather than editing their output.
The leaders are about 9.6% ahead of the public line, down from 17.5%.


---

## Comments (5 threads)

Replies are indented with `>`, one level per reply depth.

**Syed Asad Ali** · 2026-11-14 17:21 UTC · 1 votes

@manwithacat why does your table start with Teams ≥ 26 and below? The leaderboard has moved a lot now since my 27+ stack release, so from now onwards please include Teams ≥ 27, Teams ≥ 28, Teams ≥ 29 and also Teams ≥30 (Zero at the moment, but will soon cross that).

> **manwithacat** — **Topic author** · 2026-11-14 18:26 UTC · 0 votes
>
> yes, good point. i hard-coded the score boundaries and forgot about them. maybe i should switch it to a bar chart display.

**D Chakraborty** · 2026-11-14 14:06 UTC · 1 votes

>Slingshot Effect, New Team Name, Sankalp Sanand and Annontrix are new to the top ten.

Not new lol this is from 26th of October

![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F5819143%2Fad42b1329723341b97b76dd02cfc1328%2FScreenshot%202026-11-14%20at%205.04.43PM.png?generation=1794665118234410&alt=media)

> **manwithacat** — **Topic author** · 2026-11-14 16:11 UTC · 0 votes
>
> that one is on me. I only started collecting granular data about a week into the competition.

**pjmathematician** · 2026-11-14 11:27 UTC · 1 votes

Cat sitting on top of the mac mini is the only acceptable reason for the delay😼😼😼

**manwithacat** — **Topic author** · 2026-11-14 19:32 UTC · 0 votes

Thanks for the notes, all fair. Here's what we've changed:

1. **Score bands.** The "Teams ≥" table now follows the board instead of fixed cut-offs. It shows six bands from just above the top score, so ≥ 30 down to ≥ 25 today.
2. **"New to the top ten."** This is now only worked out against the previous briefing's snapshot. If that snapshot is missing, we make no claim. Apologies to Slingshot Effect, New Team Name, Sankalp Sanand and Annontrix, who've been there since 26 Oct.
3. **Credits, builder first.** Each circuit's credit now starts with whoever wrote the builder or decoder, then each change in order. For the 46-qubit inverse square root and the current DQI decoder that's Syed Asad Ali first, then Ali Haydar Özdağ / Ertuğrul Özer / Areeb Ali, then later refits. Syed's Apache-2.0 NOTICE is now carried verbatim in our public circuits dataset and in the Hedge Wizardry dataset. The Friends notebook and Hedge Wizardry text have been updated to match.

The credit chains are getting quite complex, with files forked, refit and republished by several people within hours. We'll make a best effort to keep them coherent, but we will get some wrong. If you spot an incorrect or missing credit, please say so in the discussion and we'll fix it.

(Prepared with AI assistance.)

> **Syed Asad Ali** · 2026-11-14 19:54 UTC · 0 votes
>
> Thank you, this reads correctly now: builder first, then each change in order, with the NOTICE carried along. I appreciate you taking the time to rework the bands and the top-ten wording too. Credit chains will keep shifting as files get forked and refit; if I spot anything off I’ll flag it here, and the NOTICE in our dataset stays the reference for the original builders.

**Syed Asad Ali** · 2026-11-14 12:44 UTC · 0 votes

@manwithacat Thanks for the careful re-verification. A lineage correction for your dataset and next briefing.

Your own table shows teams at 26 or above going from 45 to 812 between 12 Nov 17:34 and 14 Nov 09:37 UTC, the 1.5 days after our portfolio release. As your briefing notes, most of that +7.2 % step came from our 46-qubit inverse square root, and many people forked the portfolio.

The same lineage holds for DQI: 821bf8d157c0 is our v12 decoder (7cbe8ecf6f6e) rebuilt at extractor cutoff 18 (Ertuğrul Özer’s idea) with Areeb Ali’s 241 retained products, so the decoder itself, and every file derived from it, starts with our builder. Likewise 4560ac0c and its refits are built on our 46-qubit generator.

Both builders are Apache-2.0, which asks that derived copies keep the original attribution and NOTICE. We’d be grateful if notebooks shipping these files credit the builder first, then the improvements. Thanks to Ertuğrul, Areeb and Ali Haydar for theirs.
