# RezOps Briefing: the plateau moves to 25.09 — measure-first rewrites cut invsqrt 1.8%, DQI loses 14 Toffolis, refit chemistry weights

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756824 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-11-07 19:07 UTC · 0 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps runs the background work of
our research programme: it scans new public notebooks, re-verifies their circuits, maintains our public circuits
dataset and the notebook that blends the best of it, and drafts this briefing. A human (me) sets the direction,
reviews and approves every post. Numbers are labelled: *observed* (leaderboard), *author-reported* (the author's
official-verifier run), or *re-verified by us* (the unchanged official verifier, run locally at its default salt).
DQI circuits are ranked by their mean cost over 48 fixed draws (*our estimate*), not by one draw. Window: 06 Nov 17:06
UTC to 07 Nov 18:51 UTC. Corrections welcome.

## 1. Leaderboard (07 Nov, 18:51 UTC, observed)

| | 06 Nov, 17:06 UTC | 07 Nov, 18:51 UTC |
| --- | --- | --- |
| Top score | 29.30100 | **29.30100** (James & Geremie & Tony) |
| Teams ≥ 26 | 25 | 30 |
| Teams ≥ 25 | 96 | 353 |
| Teams ≥ 24 | 1,157 | 1,282 |
| Teams ≥ 23 | 1,499 | 1,600 |
| Teams ≥ 22.5 | 1,685 | 1,775 |
| Teams ≥ 22 | 1,718 | 1,808 |
| Teams ≥ 20 | 1,807 | 1,896 |
| Teams ≥ 15 | 1,870 | 1,955 |
| Median | 23.90494 | 23.97140 |
| Teams | 2,481 | 2,570 |

The top ten are 27.03–29.30 (James & Geremie & Tony, Schrödinger’s cat lived, taiseiu, Jack, Breaking1800, Diana,
CroDoc & Pavel, Thinking..., Yan Zhou, Bingran You). The median team is at 23.97.

**Approximate medal lines today** (Kaggle's rules for 1,000+ teams: gold top 10 + 0.2%, silver top 5%, bronze top
10%, percentages rounded down; see Kaggle's published [medal rules](https://www.kaggle.com/progression/competitions);
on 2,570 teams):

| | Rank needed | Score at that rank now | Previous | Teams within one DQI failure (0.153%) below |
| --- | --- | --- | --- | --- |
| Gold | ≤ 15 | 26.61949 | 26.49520 | 1 |
| Silver | ≤ 128 | 25.09052 | 24.95166 | 95 |
| Bronze | ≤ 257 | 25.03102 | 24.94093 | 103 |

Ties at the frontier (observed): 25.09052 held by 33 teams (ranks 114-146); 25.09026 held by 3 teams (ranks 147-149);
25.09022 held by 31 teams (ranks 152-182); 25.08978 held by 3 teams (ranks 189-191); 25.06895 held by 5 teams (ranks
208-212); 25.06885 held by 3 teams (ranks 213-215); 25.06819 held by 3 teams (ranks 224-226); 25.06570 held by 7 teams
(ranks 230-236). Ties are broken by the earlier submission
([rules](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/rules)).

The public plateau moved from 24.94 to 25.09 (+0.6%) in a day: 353 teams are now at 25 or above, against 96 yesterday.
The silver line sits on the largest tie group, 33 teams at 25.09052. That is the public score of Kirill Tymchishin's
[*Pi-Gauge & Weight Refinement*](https://www.kaggle.com/code/kirilltym/quantum-circuits-pi-gauge-weight-refinement)
(v21). Next come maoxin123's [*LB 25.09*](https://www.kaggle.com/code/maoxin123/lb-25-09-syed-circuits-refit-chem-weights)
at 25.09022 (v3, 31 teams), alex chilton's
[*24.94 LB: Verified Three-Task Portfolio*](https://www.kaggle.com/code/alexchilton/24-94-lb-verified-three-task-portfolio)
at 25.06895 (v122), and Syed Asad Ali's
[*VQCC Own-Circuit Portfolio Notebook*](https://www.kaggle.com/code/syed3000/vqcc-own-circuit-portfolio-notebook) at
25.06570 (v11); all are public scores shown on the notebooks (*observed*). The top ten did not change score.

## 2. Public frontier by task (re-verified by us unless stated)

| Task | Best public | Qubits | Cost | Ratio to reference | Source | Change since last briefing |
| --- | --- | --- | --- | --- | --- | --- |
| invsqrt | `95621d7ee3e1`: Syed Asad Ali's v11 circuit (`829b2b88d3b1`) with hosen42's measure-first rewrite applied to two write-only runs (−4 Toffolis on every shot) | 62 | 96,020 | 21.63 | [maoxin123, *LB 25.09 \| Syed circuits + refit chem weights*](https://www.kaggle.com/code/maoxin123/lb-25-09-syed-circuits-refit-chem-weights) (version 3) | was `c27e83acb576` at 97,747 (−1.77%) |
| dqi | `2c94a1d57125`: a measurement peephole (14 Toffolis whose target is X-measured become a measurement plus a conditional CZ) and caching of products that are equal as Boolean functions | 892 | 3.7504 × 10^11 (48-draw mean) | 65.17 | [Syed Asad Ali, *VQCC Own-Circuit Portfolio*](https://www.kaggle.com/datasets/syed3000/vqcc-own-circuit-portfolio) (dataset version 10) | was `7b6cbb92a797` at 3.7508 × 10^11 (−0.01%) |
| chemistry | `dcd5500ea868`: Syed Asad Ali's table and 15 loaders with maoxin123's refit weights (λ_eff 46.89199) | 427 | 2.99662 × 10^8 | 11.20 | [maoxin123, *LB 25.09 \| Syed circuits + refit chem weights*](https://www.kaggle.com/code/maoxin123/lb-25-09-syed-circuits-refit-chem-weights) (version 3) | was `da4d23796b13` at 2.99673 × 10^8 (−0.004%) |

Best-per-task geometric mean: **25.08389** (*static model / projection*: our arithmetic on the re-verified costs, DQI
at its 48-draw mean; no team has scored it); at the last briefing 24.93402 (+0.60%). All frontier sources are
Apache-2.0 (notebook footers and the dataset page read today).

Notes on credit and on the board:
- **invsqrt** went down in three steps today, all on Syed Asad Ali's line:
  - 97,747 → 96,958 (`af846e54422c`, Syed's dataset v10): shorter clean-up and shift networks.
  - → 96,268 (`829b2b88d3b1`, v11): more proven rewrites.
  - → 96,020 (`95621d7ee3e1`): maoxin123 applied the measure-first rewrite that Md. Hamid Hosen (hosen42) first used in
    [*VQCC Our Submission*](https://www.kaggle.com/code/hosen42/vqcc-our-submission) (v5) to Syed's v11 file.
  - Kirill Tymchishin's notebook ships the same invsqrt file.
- **chemistry:** alex chilton's `d3650d7ecc04` reaches the same cost within 10^-10 relative, a tie under our rule.
- **DQI on the public board.** The cheapest DQI file in expectation is not the one that scores best on the public
  leaderboard today.
  - The public board uses one fixed draw per file. maoxin123 reports that a portfolio with `2c94a1d57125` scored about
    two failed shots there, so they used Ahmed Berat Özer's `77bace9de69c` instead, which had zero failures on four
    salts they ran (*author-reported*).
  - A submission of the public parts with `2c94a1d57125` scored 24.93934 (*observed*). That is exactly the
    two-failure value of the official formula.
  - Kirill Tymchishin's 25.09052 uses `337cb0497b4d`, the same 14-site measurement identity applied to Özer's
    `77bace9de69c`, with a written proof (see section 3).
  - Over our 48 fixed draws the five 892-qubit files are within 0.035% of each other. `2c94a1d57125` is still the
    cheapest in expectation, which is what matters for the private leaderboard.

## 3. What's new

We added 22 circuits since the last briefing, each re-verified and credited.

**A measurement identity for Toffolis whose target is measured out.**
- If a CCX(a, b, t) is followed by X-basis measurements of a and t (with t not read in between), the CCX can be
  replaced: measure t first, apply CZ(a, b) conditioned on that outcome, then measure a. The phase is the same for
  every pair of outcomes, so the circuits agree on every input, including superpositions, at the same width.
- Kirill Tymchishin wrote out the proof in
  [*Pi-Gauge & Weight Refinement*](https://www.kaggle.com/code/kirilltym/quantum-circuits-pi-gauge-weight-refinement)
  (v21). He applied it at 14 sites of Özer's DQI `77bace9de69c` (−14 Toffolis), checked it with an MPS simulator
  over all 32 local basis/outcome combinations and 2,048 decoder inputs, and published the composition as
  `337cb0497b4d` (*author-reported*; the file passes our re-verification).
- Syed Asad Ali's `2c94a1d57125` (dataset v10, released earlier the same day) contains the same peephole at 14 sites,
  together with caching of products that are equal as Boolean functions.
- hosen42's measure-first rewrite on invsqrt is the same idea for write-only runs. A wire that is only ever written
  and then measured can be measured first, with its writes turned into phases. That removes 4 Toffolis on every shot.

**invsqrt: Syed Asad Ali's v10 and v11.**
- Syed's notebook lists the new rewrites with a proof note for each:
  - the k = 9 and k = 10 clean-up cascades removed by holding class bits in idle qubits;
  - the normalise/unnormalise shift networks shortened from 52 to 47 and 52 to 48 controlled swaps;
  - a 5-product coefficient decode after the hi pass re-synthesised with 4 products.
- Each is reported exact on all 3,133,452 inputs and all measurement outcomes (*author-reported*).
- Mean Toffolis 1,281.48 → 1,268.75 in v10. The v11 file is 96,268 (re-verified by us).

**Chemistry: refit weights and an estimate of what is left.**
- maoxin123 refit W and H₂ᴵ on Syed's unchanged words and loaders: λ_eff 46.89276 → 46.89199. They credit a chain
  that starts with Isabe11aM00re's weights, Heliosli's first public refit to the gate, and our `dfthc-refit`.
- They also report a bound: at fixed words, about 0.09% of chemistry cost is left (λ_eff ≥ 46.8497 from a convex
  relaxation), worth about +0.008 on the score (*author-reported*).
- The same notebook's "What we measured" section lists two more nulls:
  - invsqrt's 41 rotation angles are at the exact Chebyshev optimum;
  - holding a coherent 24-bit DQI denominator through the middle stages needs about 25 free wires that the 892-qubit
    file does not have.

**Other new files.**
- Clare Laurent's [*Observatory - Three Audited Circuit Studies*](https://www.kaggle.com/code/clarelaurent/observatory-three-audited-circuit-studies)
  (v6) adds a chemistry variant (`8f46845528ad`).
- KOUSHIK KUMAR DINDA's [notebook](https://www.kaggle.com/code/koushikkumardinda/verifiable-quantum-circuit-compilation)
  (v1) ships the reference invsqrt.
- samanyu1808's [*VQCC Credited Stack v4*](https://www.kaggle.com/code/samanyu1808/vqcc-credited-stack-v4) is a
  credited redistribution of the public parts.

## 4. From the forum

- **Do the top chemistry loaders satisfy the XOR-lookup contract?** (asked under our last briefing,
  [756558](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756558)).
  - Zero found that yoyobar's ANGLE_1 fails a zero-payload XOR test on 761 of 1,220 addresses while passing the
    official verifier.
  - We reproduce the 761/1,220 on the current top tables. It is expected: these loaders rely on the word registers
    already holding the previous batch.
  - The hosts ruled in [753092](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753092)
    that ANGLE_k only needs to be correct on states reachable in the forward and reverse traversals, which is what
    the verifier checks.
  - A plain QROM family passes the zero-payload test at every k, so the test is a good way to tell strict loaders
    from reachable-only ones.
- **Formal proofs of the GF(4)-bilinear multiplication counts**
  ([756476](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756476)).
  - Following 田村's computer-assisted proofs and pjmathematician's asymmetric sketch (last briefing), Dyst-0 posted
    Lean 4 + Mathlib proofs of μ₄(5) = 11 and μ₄(6) = 14, the asymmetric degree-six case included, with the source on
    GitHub.
  - Three finite searches use `native_decide`; everything else is kernel-checked (*author-reported*).
  - So three independent groups now agree that no 13-term GF(4)-bilinear multiplier for GF(4096) exists.
- **Prize eligibility for portfolios built from public Apache-2.0 notebooks**
  ([756745](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756745)).
  - Mondler asks whether a byte-for-byte assembly of public outputs is eligible for medals.
  - They also ask what Rules B.5/B.7 require when an upstream lineage has no stated licence.
  - No host answer yet. Given how many teams now sit on public parts, the answer matters to a large part of the board.

## 4b. Position risk: the final DQI draw, circuit by circuit

| DQI file (SHA-256 prefix) | Qubits | Board failures | P(final worse) | P(final better) | Expected change | One more failure |
| --- | --- | --- | --- | --- | --- | --- |
| `4ee7a17b4cbc` | 889 | 0 | 0.20 | 0.00 | -0.033% | -0.153% |
| `b3354f60e459` | 897 | 2 | 0.00 | 0.96 | +0.233% | -0.116% |
| `7beea291a5de` | 889 | 0 | 0.63 | 0.00 | -0.142% | -0.153% |
| `f3d335c3f693` | 877 | 0 | 0.67 | 0.00 | -0.158% | -0.153% |
| `ef334fee63f3` | 883 | 0 | 0.28 | 0.00 | -0.050% | -0.153% |
| `0161fd2a6e50` | 889 | 1 | 0.31 | 0.33 | -0.005% | -0.129% |

The table covers older files. For today's 892-qubit files (*our estimate*):
- **`2c94a1d57125`** drew 2 failed shots on the public board. Over our 48 fixed draws it had 0 failures in 39 and 1
  in 9, so a fresh draw has about a 17% chance of at least one failure and about 1.6% of two or more. Its final
  score is therefore likely to be better than its public one.
- **Files with 0 public failures** (e.g. `77bace9de69c`, `337cb0497b4d`) have about the same per-shot failure rate, so
  their final draw is more likely to be worse than their public one.
- **What this means for the medal lines.** One failed shot is worth 0.15% of the overall score. The plateau groups
  are 0.001–0.03% apart, so the final ordering among teams on public parts will be set largely by the DQI draw.

## 4c. Prize watch (leaderboard-timed prizes)

**Prize watch** (from 112 leaderboard snapshots; latest 07 Nov 18:51 UTC; our estimate)

- **Leader:** James & Geremie & Tony at 29.30100 (best submission 06 Nov 00:25 UTC); next is Schrödinger’s cat lived
  at 28.69687, 2.11% behind.
- **Mid-Competition Milestone Prize** (highest public score at 20 Nov 23:59 UTC, notebook published): 13.2 days left.
- **Longest Leader Prize** (most time at #1 by 4 Jan 23:59 UTC): James & Geremie & Tony has led since between 02 Nov
  20:26 and 02 Nov 21:25 UTC, i.e. 4.9 days. If it stays #1, no other team can catch its total after **10 Dec
  08:26–10 Dec 09:25 UTC** (start of its stint plus half the 75-day window; worst case, eligibility not checked).

## 5. Tools and notes from us

- **Circuits dataset and blend notebook.** The
  [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits) is now v110, with 388
  credited circuits in `manifest.csv`. The
  [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made) is now v121. Its
  "what each pick is built from" table no longer shows "?" for picks whose lineage we have not reviewed yet; it says
  so in words.
- **Kaggle's DQI verifier fails exactly the same shots as a local run, seed by seed.**
  - We replayed 8 seeds of `2c94a1d57125` with the unmodified verifier and the organisers' kickmix build in a Kaggle
    CPU notebook, then ran the same seeds locally on both the organisers' build (compiled on a Mac) and our
    kickmix-engine.
  - All three agreed on every seed, failures and costs: 2 seeds with one value error, 6 with none.
  - So local salt sweeps predict hosted draws exactly, and a file's public score fixes its public draw: resubmitting
    the same bytes cannot re-roll it.
- **Chemistry loader census.** We ran our `dfthc-refit` analysis tools on all 964 distinct public ANGLE_k loaders;
  843 are in the 8-bit-word families the tools support.
  - Where several public loaders implement the same words, the cheaper one never used fewer caches. In 85 of 86 cases
    with a gap of 5 or more Toffolis, it exploited more linear relations among the table's rows.
  - Every loader on today's best table sits at the bound for its own structure, with no spare low caches.
  - The run also found a bug in our tool: structure recovery assumed 5 words per loader, so loaders with 3 or 4 words
    were misread. It is fixed in our repository.

## 6. What we think it means (opinion)

- **Measurement identities are the theme of the day.** Measure-first on invsqrt and the measured-target Toffoli on DQI
  are both local, provable and free in width. They are easy to search for: any CCX whose target is next X-measured
  without being read is a candidate. We expect more of them in every task, including chemistry clean-up.
- **The public board now ranks DQI draws as much as circuits.** The five best 892-qubit DQI files are within 0.035%
  in expectation, while one failed shot costs 0.15%. Choosing a file by its public draw helps on the public board but
  not on the private one, where every file gets a fresh draw.
- **Chemistry at fixed words looks nearly exhausted.** maoxin123's bound (about 0.09% left) and our loader census
  point the same way. As on 05 Nov, the next real chemistry step probably needs a new table.
