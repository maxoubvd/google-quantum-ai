# Sprint 2a reverted: back to the sequential plan

- **Date:** 2026-11-18 09:10 UTC (just after the launch entry `2026-11-18-0908-sprint-2a-launch.md`)
- **Plan step:** none (corrects the Sprint 2a launch)
- **Files:** `Analysis/PLAN.md` (Sprint 2a subsection removed: restored to the committed version);
  `generators/` (removed); WSL `~/work/inv`, `~/work/dqi`, `~/work/chem` (removed)

## What was done

The three engineering agents (invsqrt, dqi, chemistry) were stopped a few minutes after launch, at the user's
request, and the sprint was undone. They had only copied Syed Asad Ali's v12 source trees into
`generators/<task>/` (identical to `third_party/`, checked with `diff -rq`) and added two small build helper
scripts; nothing was built, gated or changed in the sources. All of it was deleted. The plan is back to the
committed version: Phase 2 runs milestone by milestone, one step at a time, with a review after each step.

## Why

Asked whether the sprint was a good idea, the assessment was: parallel work across independent tasks is
sound, but the 19 Nov deadline brings almost nothing (the milestone prize is out of reach at 27.0 vs 31.6; the
public board today does not decide the final ranking), the time pressure invites rushed choices (e.g. a dqi
cutoff accepted too fast raises the private-salt risk), it bypasses the user's step-by-step reviews, and three
heavy agents risk exhausting the usage quota (four research agents hit the 5-hour limit the night before). The
user chose to stop the agents and continue as originally planned.

## Results and analysis

No circuit or score change. Lesson recorded: a deadline must bring value of its own before it justifies
changing the way of working; and parallel agents are a choice to make with the user, step by step.

## Next

M2.2, step S2.5 (INV-A02, the 64-row table), done in the main session as originally planned.
