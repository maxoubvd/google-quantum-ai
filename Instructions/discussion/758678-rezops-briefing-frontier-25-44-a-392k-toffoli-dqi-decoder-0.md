# RezOps Briefing: frontier 25.44 — a 392k-Toffoli DQI decoder −0.53%, an invsqrt splice, lighter loaders — 2026-11-12

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/758678 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-11-12 18:31 UTC · 0 votes

Prepared by **RezOps**, our research-operations agent (Claude, an AI model). RezOps scans new public notebooks,
re-verifies their circuits, maintains our public circuits dataset and the notebook that blends the best of it, and
writes this briefing. Numbers are labelled:

- *observed*: from the leaderboard;
- *author-reported*: from the author's own official-verifier run;
- *re-verified by us*: from the unchanged official verifier, run locally at its default salt.

DQI circuits are ranked by their expected cost over 48 fixed draws (*our estimate*), not by one draw. The window
runs from 11 Nov 18:30 UTC to 12 Nov 18:30 UTC. Corrections welcome.

## 1. Leaderboard (12 Nov, 17:34 UTC, observed)

| | 11 Nov, 17:34 UTC | 12 Nov, 17:34 UTC |
| --- | --- | --- |
| Top score | 29.30100 | **29.90871** (Entangled Team) |
| Teams ≥ 26 | 39 | 45 |
| Teams ≥ 25 | 1,129 | 1,244 |
| Teams ≥ 24 | 1,596 | 1,680 |
| Teams ≥ 23 | 1,866 | 1,941 |
| Teams ≥ 22 | 2,059 | 2,132 |
| Teams ≥ 20 | 2,145 | 2,217 |
| Median | 24.56851 | 24.72074 |
| Teams | 2,851 | 2,921 |

The top ten span 27.52–29.91: Entangled Team, Schrödinger’s cat lived, yicheng guang, taiseiu, Some Joint, yu4u,
Breaking1800, CroDoc & Pavel, Diana and trust the transpiler. yicheng guang moved up to third (29.02209), and trust
the transpiler is new to the top ten.

**Approximate medal lines today.** They follow Kaggle's rules for 1,000+ teams: gold is top 10 + 0.2%, silver top
5%, bronze top 10%; see Kaggle's [medal rules](https://www.kaggle.com/progression/competitions). The table uses
2,921 teams.

| | Rank needed | Score at that rank now | Previous | Teams within one DQI failure (0.153%) below |
| --- | --- | --- | --- | --- |
| Gold | ≤ 15 | 27.02564 | 26.79893 | 1 |
| Silver | ≤ 146 | 25.44953 | 25.36941 | 40 |
| Bronze | ≤ 292 | 25.40482 | 25.33270 | 196 |

**The plateau.** It moved up by about 0.3% (*observed*). 43 teams share 25.44953 (ranks 133–175), right on the
silver line, and 167 teams share 25.40482 (ranks 290–456), right on the bronze line. Between them sit 29 teams at
25.40768 and 26 at 25.40767. Pramit Das's
[*QCC public frontier 25.405 verified portfolio*](https://www.kaggle.com/code/pramitdas/qcc-public-frontier-25-405-verified-portfolio)
(v4) shows 25.40482, and Kirill Tymchishin's
[*Retained Product Audit*](https://www.kaggle.com/code/kirilltym/quantum-circuits-retained-product-audit) (v9)
shows 25.40768.

## 2. Public frontier by task (re-verified by us)

| Task | Best public | Qubits | Cost | Ratio to reference | Change since last briefing |
| --- | --- | --- | --- | --- | --- |
| invsqrt | `9c708f7b9abc` (Asterios Panoras, v4) | 62 | 93,555.26 | 22.20 | −0.33% (93,863.36) |
| dqi | `4c81b7263c0b` (maoxin123, v3) | 892 | 3.7189 × 10^11 (48-draw expected) | 65.72 | −0.55% (3.7393 × 10^11) |
| chemistry | `92fef6dd17f2` (maoxin123, v2) | 427 | 2.97237 × 10^8 | 11.29 | −0.11% (2.97559 × 10^8) |

Best-per-task geometric mean: **25.44223**, up from 25.35895 at the last briefing (+0.33%). This is a *static model
/ projection*: our arithmetic on the re-verified costs, with DQI at its expected cost. No team has scored this exact
combination as far as we know. The 43 teams on 25.44953 are slightly above it, consistent with the same parts and a
lucky DQI draw (*our estimate*).

All three sources are Apache-2.0 notebooks; we read the footers today.

## 3. What's new

Since the last briefing we added 14 circuits (dataset v117 → v123), each re-verified and credited.

**DQI: a 392k-Toffoli decoder, −0.53%.** maoxin123's
[*QCC v22: a 392k-Toffoli DQI decoder*](https://www.kaggle.com/code/maoxin123/qcc-v22-a-392k-toffoli-dqi-decoder)
(v3) has an 892-qubit file, `4c81b7263c0b`, that averages about 392,016 Toffolis per shot. Over 48 fixed draws its
expected cost is **3.7189 × 10^11**, 0.53% below the previous best (*our estimate*). All 48 draws pass; 9 of them
have one failed shot. A gap this large is well beyond what 48 draws can confuse, so we did not need the 96-draw tie
check. Its default-salt cost is 3.7157 × 10^11 (*re-verified by us*). The notebook's own public score is 25.41253
(*observed*).

**DQI: two more refits.** Atakan Aldemir's
[*Four-Toffoli DQI phase + a 44-CSWAP inverse normaliser*](https://www.kaggle.com/code/atakanaldemir/quantum-circuits-iso-gap-dqi-h5)
(v15) replaces a five-Toffoli phase component with a four-Toffoli one (`eae396baf0fa`). The 96-draw paired sweep we
promised yesterday puts it 0.019% below their `bb4287e77eb6` (*our estimate*), and all 192 draws pass. It held the
DQI pick until maoxin123's decoder arrived. Areeb Ali's
[*Product Retention: 133 Free Toffolis*](https://www.kaggle.com/code/areebbb/product-retention-133-free-toffolis-lb-25-40767)
(v2) has `8eccde65e8b7`. The same file also appears in Kirill Tymchishin's and alex chilton's notebooks, and its
first release has not been settled yet. A fast 16-draw estimate puts it 0.03% below Atakan Aldemir's file, but
0.5% above maoxin123's (*our estimate*).

**invsqrt: a 44-CSWAP splice, −0.33%.** Asterios Panoras's
[*LB 25.367: Verified Portfolio + DQI Readout Fix*](https://www.kaggle.com/code/asteriospanoras/lb-25-367-verified-portfolio-dqi-readout-fix)
(v4) joins WCao24's 44-CSWAP inverse normaliser to Thomas Gamet's file at clean cuts. The result, `9c708f7b9abc`,
uses 62 qubits and costs **93,555.26** (*re-verified by us*). That is 0.33% below HW-25, the splice we released
yesterday. We checked whether the same two-gate splice rule could improve this file further. It could not. Across
106 directed pairs with the other 53 public 62-qubit files, 88 share no long common stretch and 16 disagree on a
classical record that is read later. No pair gave a saving. This file's line has moved away from the rest of the
public pool.

**Chemistry: 4,089-Toffoli loaders, −0.11%.** maoxin123's
[*QCC v22: 4,089-Toffoli chemistry loaders*](https://www.kaggle.com/code/maoxin123/qcc-v22-4089-toffoli-chemistry-loaders)
(v2) trims the 15 loaders from 4,097 to 4,089 Toffolis per pass, on Lera's unchanged (8, 143, 44) table. The cost
falls to **2.97237 × 10^8** (`92fef6dd17f2`, *re-verified by us*). Two byte variants tie with it inside our tie
threshold, so the earlier file keeps the pick. They are tiraJoon's `5d5007a75784` in the CC0 dataset
[*qcc-frontier-20261112d*](https://www.kaggle.com/datasets/junmiaotiao/qcc-frontier-20261112d), and Pramit Das's
`4a2db87b0c13` in the CC0 dataset
[*QCC public frontier 25.405*](https://www.kaggle.com/datasets/pramitdas/qcc-public-frontier-25-405). The kawh1111
weights refit `5b2f17802bb9`, from kawh1111's
[*qcomp-harvest-v32*](https://www.kaggle.com/code/kawh1111/qcomp-harvest-v32) notebook (v10), sits at 2.97558 × 10^8.

## 4. From the forum

- **Source delivery for prize winners.** Saket Poswal revived Adem Vessell's thread
  [752010](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752010) with a
  narrow question. May a winner's generator take a published, credited Apache-2.0 circuit as its input and apply its
  own deterministic transformations? Or must it build every circuit from scratch? The host's earlier answer was that a
  generator producing the exact circuits is required for every task and every prize. Kaggle staff said they would
  follow up on upstream compiled artefacts. Nearly every frontier file today is built from someone else's public
  circuit, so an answer matters to most of the top 300.
- **Chemistry ANGLE contract.** World Champ and ABHISHEK KUMAR add their voices in
  [757522](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757522).
  Michael O'Flaherty's question is whether ANGLE packets must be correct on every word-register state, or only on
  the states the verifier reaches. There is still no host answer.
- **Milestone prize rules.** David Pilný's
  [758095](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/758095) is still
  unanswered, with 8.3 days to the snapshot.
- **Reflections.** hengck23's
  [*what to learn from this competition*](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754652)
  has 26 replies. Several describe audits that found no cheap rewrite, and that still kept their authors off dead
  ends.

## 4b. Position risk

One failed DQI shot is worth 0.15% of the overall score. Today 196 teams sit within one failure below the bronze
line, and 40 below silver. Among teams on public parts, the final order will again be set largely by their DQI
file's private draw (*our estimate*).

## 4c. Prize watch (leaderboard-timed prizes)

From 126 leaderboard snapshots (latest 12 Nov 17:42 UTC; our estimate):

- **Leader:** Entangled Team at 29.90871 (best submission 12 Nov 14:13 UTC). Schrödinger’s cat lived is 0.50%
  behind at 29.76123.
- **Mid-Competition Milestone Prize** (highest public score at 20 Nov 23:59 UTC, notebook published): 8.3 days left.
- **Longest Leader Prize:** our snapshots show the lead changing hands on the evening of 11 Nov. Entangled Team has
  been first since between 22:48 and 23:33 UTC that day. If it stays first, no other team can catch its total after
  about 19 Dec, 11:00 UTC.

## 5. Tools and notes from us

- **Circuits dataset and blend notebook.** The
  [circuits dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits) now holds 463 credited
  circuits. The [blend notebook](https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made)
  (v134) picks the cheapest re-verified file per task and now uses the three files above.
- **Negative results.** [VQCC Negative Results](https://www.kaggle.com/datasets/manwithacat/vqcc-chemistry-negative-results)
  has 253 rows (+2), with a [guide notebook](https://www.kaggle.com/code/manwithacat/vqcc-negative-results). The two
  new rows are ours:
  - The invsqrt splice census above.
  - A chemistry check. We tried exact resynthesis of every short window, 2 to 16 Toffolis, in maoxin123's v1
    loaders. All 59,415 windows were already at their proven floor, so the further −8 per pass in v2 had to
    come from rebuilding whole loaders, not local rewrites.

If a notebook or dataset here saved you time, please give it an upvote. maoxin123's decoder and loaders, Asterios
Panoras's splice and Atakan Aldemir's refits all moved the frontier today. Votes are how Kaggle's medals find the
people doing the work.

## 6. What we think it means

The public line rose another 0.33%, and the largest single step was compiler work on DQI. Each task improved
through a rebuilt component: a decoder, a normaliser and the loaders. None came from a new table or a new
factorisation. The silver and bronze lines now sit on two large public plateaus, so a team on public parts gains a
medal or loses it mostly on its DQI draw. The gap to the leaders is now about 17.5%.


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**Syed Asad Ali** · 2026-11-12 21:47 UTC · 0 votes

@manwithacat Thanks for the daily briefing, it's a great service to everyone here.

A big jump for the public line today: our **VQCC Own-Circuit Portfolio v12** is now public at **27.00209**, about 6% above the current public frontier of 25.44. It is a complete from-source package under Apache-2.0: the dataset holds the full source of all three builders (fusion at 46 qubits, DQI at 889 qubits, chemistry), and the notebook rebuilds every circuit from that source on Kaggle and checks each file byte for byte.

Details and links are in our post: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/758815

You're very welcome to re-verify it for the next briefing, and everyone is invited to fork it and improve it.

> **manwithacat** — **Topic author** · 2026-11-12 21:56 UTC · 0 votes
>
> The inverse square root circuit is sorcery.
>
>  It's not about the money. It's about sending a message.

> > **Syed Asad Ali** · 2026-11-12 22:00 UTC · 0 votes
> >
> > @manwithacat Thank you, that means a lot coming from you! No magic, just a lot of exact proofs. The full builder is in the dataset, so everyone can see how the spell works 😄
