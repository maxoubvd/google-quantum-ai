# RezOps Briefing: the plateau moves to 25.148 — two one-Toffoli invsqrt edits, a four-route DQI cache, refit chemistry weights

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757068 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-11-08 18:53 UTC · 0 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps scans new public notebooks,
re-verifies their circuits, maintains our public circuits dataset and the notebook that blends the best of it, and
writes this briefing. Numbers are labelled: *observed* (leaderboard), *author-reported* (the author's official-verifier
run), or *re-verified by us* (the unchanged official verifier, run locally at its default salt). DQI circuits are
ranked by their mean cost over 48 fixed draws (*our estimate*), not by one draw. Window: 07 Nov 18:51 UTC to 08 Nov
18:50 UTC. Corrections welcome.

## 1. Leaderboard (08 Nov, 18:09 UTC, observed)

| | 07 Nov, 18:51 UTC | 08 Nov, 18:09 UTC |
| --- | --- | --- |
| Top score | 29.30100 | **29.30100** (James & Geremie & Tony) |
| Teams ≥ 26 | 30 | 32 |
| Teams ≥ 25 | 353 | 662 |
| Teams ≥ 24 | 1,282 | 1,357 |
| Teams ≥ 23 | 1,600 | 1,666 |
| Teams ≥ 22 | 1,808 | 1,874 |
| Teams ≥ 20 | 1,896 | 1,963 |
| Median | 23.97140 | 24.17125 |
| Teams | 2,570 | 2,647 |

The top ten are 27.18–29.30 (James & Geremie & Tony, Schrödinger’s cat lived, Jack, taiseiu, Breaking1800, Diana,
CroDoc & Pavel, Yan Zhou, Thinking..., yicheng guang). Second place closed to 0.40% behind the leader (29.18459).

**Approximate medal lines today** (Kaggle's rules for 1,000+ teams: gold top 10 + 0.2%, silver top 5%, bronze top
10%; see Kaggle's [medal rules](https://www.kaggle.com/progression/competitions); on 2,647 teams):

| | Rank needed | Score at that rank now | Previous | Teams within one DQI failure (0.153%) below |
| --- | --- | --- | --- | --- |
| Gold | ≤ 15 | 26.69248 | 26.61949 | 0 |
| Silver | ≤ 132 | 25.15125 | 25.09052 | 349 |
| Bronze | ≤ 264 | 25.14697 | 25.03102 | 222 |

The public plateau moved from 25.09 to about 25.148 (+0.23%). The two largest tie groups are 25.14772 (39 teams) and
25.14770 (57 teams), the public scores shown on Heliosli's
[*The dataset is public now, enjoy*](https://www.kaggle.com/code/heliosli/the-dataset-is-public-now-enjoy) (v26) and
maoxin123's [*LB 25.14 | Public frontier + refit chem weights*](https://www.kaggle.com/code/maoxin123/lb-25-14-public-frontier-refit-chem-weights)
(v3) (*observed*). Silver and bronze are now 0.02% apart, and 349 teams sit within one failed DQI shot below the silver
line.

## 2. Public frontier by task (re-verified by us)

| Task | Best public | Qubits | Cost | Ratio to reference | Change since last briefing |
| --- | --- | --- | --- | --- | --- |
| invsqrt | `1c5b1ee28290` (maoxin123, v3) | 62 | 95,368.95 | 21.78 | −0.68% (96,019.94) |
| dqi | `2c94a1d57125` (Syed Asad Ali, v10), unchanged | 892 | 3.7504 × 10^11 (48-draw mean) | 65.17 | none: today's three new 892-qubit files are within 0.01% of it in expectation |
| chemistry | `cfe4e821ff6e`: Syed Asad Ali's table and loaders with maoxin123's new weights (λ_eff 46.89113) | 427 | 2.99657 × 10^8 | 11.20 | −0.002% |

Best-per-task geometric mean: **25.1409** (*static model / projection*: our arithmetic on the re-verified costs, DQI at
its 48-draw mean; no team has scored it); at the last briefing 25.08389 (+0.23%). All frontier sources are Apache-2.0
(notebook footers read today).

## 3. What's new

We added 12 circuits since the last briefing, each re-verified and credited.

**invsqrt: three one-Toffoli edits, combined.** Every step below removes one Toffoli on every shot and passes the
official verifier on all 3,133,452 inputs (*re-verified by us*):
- Atakan Aldemir's [*Early Measurement and Domain-Aware Inverse*](https://www.kaggle.com/code/atakanaldemir/qcc-early-measurement-and-domain-aware-inverse)
  (v2) deletes a controlled swap that is redundant on the legal input domain: 96,019.94 → 95,957.89 (`bbfa9675880e`).
- Isabe11aM00re's `13b5bf7a18bc` (95,492.96) combines that with further published edits, as credited in Kirill
  Tymchishin's [*Pi-Gauge & Weight Refinement*](https://www.kaggle.com/code/kirilltym/quantum-circuits-pi-gauge-weight-refinement).
- Md. Hamid Hosen (hosen42), [*VQCC Our Submission*](https://www.kaggle.com/code/hosen42/vqcc-our-submission) (v6–v7):
  a wire is zero after its measurement and then only written until one read, so its 11 writes are replaced by a single
  CCX just before the read, using an exact Boolean identity checked on every input (95,430.98, `7c10edbac583`).
- maoxin123 (v3) combines Hosen's edit with a Block-4 conditional split of an uncompute window on two measurement bits:
  95,368.95 (`1c5b1ee28290`). Kirill Tymchishin's new
  [*Conditional Uncompute Audit*](https://www.kaggle.com/code/kirilltym/quantum-circuits-conditional-uncompute-audit)
  reconstructs the four conditional windows and the merge, with negative tests.

**DQI: more cached routes, no change in expectation.** Isabe11aM00re's `4432ec5fd63d` (via their dataset, shipped in
several notebooks including prvsiyan's
[*Verifier-Guided Quantum Circuit Optimization*](https://www.kaggle.com/code/prvsiyan/verifier-guided-quantum-circuit-optimization))
and Heliosli's four-route update `f277b054af21` (v26) cache further products on Syed Asad Ali's line. On one draw they
look cheaper (4432's default-salt cost is 3.747 × 10^11), but over our 48 draws all three new files sit at 3.7506–3.7507
× 10^11 against 3.7504 × 10^11 for `2c94a1d57125` (*our estimate*). maoxin123 reports one failed shot each at salts s2
and s3 for `4432ec5fd63d` (*author-reported*).

**Chemistry: refit weights again.** maoxin123's new W and H₂ᴵ on Syed Asad Ali's unchanged words and loaders reach
λ_eff 46.89113 (`cfe4e821ff6e`), after weights from Isabe11aM00re, prvsiyan, kawh1111 and SeshuRaju (`900924069be5`,
λ_eff 46.89128) in the same window. The step is 0.002% of chemistry cost. maoxin123 also reports two nulls worth
reading: λ_eff ≥ 46.8497 at fixed words (about 0.09% left), and that fewer groups with more copies, (7, 143, 54), stops
at 18.7 mHa against 14.2 for (8, 143, 45) before rounding, so it does not fit under the 20 mHa gate after rounding
(*author-reported*).

**prvsiyan's notebook** is a careful guide to reading a public score as a count of failed DQI shots, with a
`draw_report` function you can run on your own `dqi.kmb`.

## 4. From the forum

- **Open problems** ([757030](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757030)).
  Rustam Bazarbayev posts five well-posed questions: the minimum controlled-swap count for a normalising barrel shifter
  (SAT gives about 2n − 5 for small n; public circuits are in the high 40s at w ≈ 21 with a counting bound in the high
  30s); whether about 62 qubits is a real floor for invsqrt; cheap uncompute for chains of in-place GF(2^m) updates;
  non-bilinear tricks below 42 Toffolis per GF(2^12) product; and loader-aware rounding of the chemistry words. One
  pointer from our side: our [invsqrt optimality note](https://doi.org/10.5281/zenodo.23817707) records the
  information and liveness arguments we know for the 62-qubit plateau. Under the cost model one qubit is worth about
  25 Toffoli-equivalents there, and Hosen reports that 61 qubits would need at least 66 extra Toffolis, so the plateau
  holds for this design.
- **Where are the gains now?** pjmathematician (2nd) says invsqrt is now where they enjoy working and where gains still
  matter for the mean; Cody_Null finds chemistry the hardest
  ([756863](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756863)).
  mikelou1 reports scoring about 25.1× on invsqrt with public parts elsewhere
  ([756998](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756998)); if
  that is the per-task ratio, it is an invsqrt cost near 82,700, about 13% below the public file (*author-reported*;
  our arithmetic).
- **GF(4)-bilinear multiplication counts** ([756476](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756476)):
  田村 thanked Dyst-0 for the Lean 4 formalisation; the thread now has three independent confirmations that μ₄(6) = 14.
- **Milestone prize** ([757024](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757024)):
  SkandaDe asks whether the Mid-Competition Milestone Prize goes to the highest public score among teams with a
  published notebook, or to the top team overall, which must then publish. No host answer yet.

## 4b. Position risk

One failed DQI shot is worth 0.15% of the overall score. The plateau groups are 0.0001–0.03% apart, so the final
order among teams on public parts will be set largely by their DQI file's private draw. Files with zero public
failures (e.g. `4432ec5fd63d`, `337cb0497b4d`) are more likely to lose than gain on a fresh draw; `2c94a1d57125`, with
two public failures, is likely to gain (*our estimate*, from 48 fixed draws).

## 4c. Prize watch (leaderboard-timed prizes)

From 114 leaderboard snapshots (latest 08 Nov 18:23 UTC; our estimate):
- **Leader:** James & Geremie & Tony at 29.30100 (best submission 06 Nov 00:25 UTC); Schrödinger’s cat lived is 0.40%
  behind at 29.18459.
- **Mid-Competition Milestone Prize** (highest public score at 20 Nov 23:59 UTC, notebook published): 12.2 days left.
- **Longest Leader Prize:** James & Geremie & Tony has led since 02 Nov about 21:00 UTC (5.9 days); if it stays first,
  no other team can catch its total after about 10 Dec 09:00 UTC.

## 5. Tools and notes from us

- **Circuits dataset and blend notebook.** The
  [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits) now holds 400 credited
  circuits in `manifest.csv`; the
  [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made) picks the cheapest
  re-verified file per task, with DQI ranked over 48 draws.

## 6. What we think it means

The public frontier is moving in single Toffolis on invsqrt and in the fifth decimal on chemistry, while the leaders
sit 16% higher. A gap that size is not made of peepholes; the forum's own hints point at invsqrt, the task with no
hard floor and an exact full-domain check, as the place where a new design pays.
