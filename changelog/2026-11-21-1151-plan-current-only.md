# PLAN.md cleaned to the current plan only; plan history archived; new rule in CLAUDE.md

- **Date:** 2026-11-21 11:51 UTC
- **Plan step:** project rules (no step change; S2.26 stays next)
- **Files:** `Analysis/PLAN.md` (rewritten), `Analysis/plan_history/` (new: README and four archived versions),
  `CLAUDE.md` (new rule, repository description)

## What was done

- Archived the plan versions, unchanged: `PLAN_2026-11-17_initial.md` (git `1d3c03b`),
  `PLAN_2026-11-18_milestones.md` (git `dcbb4e6`), `PLAN_2026-11-19_replan.md` (git `821b50e`), all exported with
  `git show`, and `PLAN_2026-11-21_thembones-replan.md` (the composite file of 11:45 that held every revision
  inline). `Analysis/plan_history/README.md` indexes them.
- Rewrote `Analysis/PLAN.md` with the plan in force only: phases with their status, the 21 Nov snapshot and
  exchange rates at the new frontier (invsqrt 37 Q: −1 qubit = +0.90% score), constraints, licensing and where to run
  what, the method (including the rule that nulls state their domain), Phase 1 and the done part of Phase 2 in one
  paragraph each, the current Phase 2 milestones M2.E–M2.H and steps S2.25–S2.32, the break, Phase 3 (M3.1–M3.6),
  Phase 4 and the current risks. Removed: older snapshots, the Gemini assessment, the 17/18/19 Nov milestones and
  the S2.1–S2.24 step tables (all kept in the archive and in `changelog/`).
- `CLAUDE.md`: new rule — `PLAN.md` holds only the current plan; on each revision the replaced version is copied
  unchanged into `Analysis/plan_history/` with a README row, then `PLAN.md` is rewritten.

## Why

User request (21 Nov): the plan file must show only the current plan at any time, with old versions kept as
versioned documents.

## Results and analysis

No measurement. `PLAN.md` went from 436 lines to about 200.
