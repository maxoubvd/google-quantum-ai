# Re-plan of 19 Nov: public-frontier analysis, public level first, then exploration of new families

- **Date:** 2026-11-19 08:41 UTC
- **Plan step:** re-plan (Phase 2 rewritten from M2.A on; Phase 3 rewritten); decision S1.3 revised
- **Files:** `Analysis/research/public_frontier_asterios_2026-11-19.md` and
  `Analysis/research/public_frontier_others_2026-11-19.md` (written by two research agents; downloads under
  `downloads/kernels/`, `downloads/lb/2026-11-19/`, `downloads/datasets/{balramsrinivasan__vqcc-public-work-1118,
  manvenderraparia__vqcc-lp-refit-portfolio}`, read only); `Analysis/IDEAS.md` (22 statuses revised, new section 10
  with 8 ideas); `Analysis/PLAN.md` (snapshot of 19 Nov, Phase 2 re-plan with milestones M2.A–M2.D, Phase 3 table,
  targets, S1.3 row)

## What was done

- Two research agents read the public work at the ~28.2 frontier: Asterios Panoras's notebook *QCC 28.23162 | What
  we tried, what failed and why* (77 cells) and ~25 other recent or top notebooks, plus the live board.
- Their findings were mapped onto `IDEAS.md`: ideas already in the public frontier files are marked IN PUBLIC
  FRONTIER, measured public failures NULL (public) or PARTLY NULL, and eight new ideas were added (INV-N01, DQI-N01,
  DQI-N02, DQI-N11–N13, CHEM-N11, GEN-N01).
- Re-plan with the user's decisions: (1) decision S1.3 revised: credited public Apache-2.0/CC0 files may be
  submitted for the leaderboard; (2) a second submission at the public level first (M2.A); (3) phases of
  exploration of other families, invsqrt first, then dqi, chemistry in January (M2.B, M2.C, M3.2–M3.5), with a
  paper → scaled prototype → full-scale gate per family.

## Why

The user asked to drop from the plan what other teams already tested, to rebuild the plan from the new public
base, to submit at the public level, and to explore other families in order to stand out.

## Results and analysis

- Board, 2026-11-19 08:27 UTC (3,377 teams): leader 31.74; gold 28.824 (rank 16); silver 28.2344; bronze 28.2316;
  153 teams at 28.23162. Public frontier 28.23445: invsqrt `ecf908fe` (44 Q, 74,383.70), dqi `e39b1396` (890 Q,
  3.4526e11), chemistry `6ef06081` + 4,022-T loaders (294,642,169.5). Against our v12 base: −7.2%, −4.3%, −0.7%.
- **Accumulation, not new families:** every top public portfolio is a stack of patches on Syed Asad Ali's invsqrt
  and dqi generators and Lera's coded chemistry table; no public family change since 15–16 Nov; the frontier files
  are exhausted by every public scan. Past big jumps came from family changes (invsqrt 62 → 46 qubits; successive
  dqi decoders; coded chemistry words). The leaders need about −30% total cost, i.e. new families.
- New public nulls folded into `IDEAS.md`: angle-only refits saturated (INV-A12 part), rotation merges (INV-A10),
  loader Toffoli removal (CHEM-A11), weight refits (CHEM-A13), brot 7 (CHEM-B11, 27.1 mHa), dropping rare dqi
  slots (DQI-A07 variant).
- Ali Haydar Özdağ published the source of his 64-row table on Syed's v12 generator (`ahoirg/qcc-invk6-k6c`, 19 Nov),
  which would unblock S2.9; with the revised strategy it moves to the backlog (prize-eligible own line only).

## Next

M2.A / S2.17: assemble, pin, gate and package the public frontier set (28.23445) for the user's second submission.
