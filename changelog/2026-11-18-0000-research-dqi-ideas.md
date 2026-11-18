# Research: prioritized dqi ideas

- **Date:** 2026-11-18 00:00 UTC
- **Plan step:** none (research for the Phase 2 idea register)
- **Files:** `Analysis/research/dqi.md` (new)

## What was done

Wrote `Analysis/research/dqi.md`: cost anatomy of our v12 dqi decoder (per stage, from the generator's own
source ledgers and our default-salt verifier run), state of the art with URLs, known nulls, and 31 ideas
(19 incremental, 12 rewrites) scored by P x gain / effort, with combination notes.

## Why

To choose the dqi work of Phase 2 from a single ranked list instead of ad hoc ports.

## Results and analysis

- v12 static Toffolis: key stage 216,928 (56.5%; forward 101,087, reverse 115,841), of which 117,810 are 1,190
  GF(2^24) multiply-accumulates at 99 T; convert 114,423 (29.8%); scan 44,286 (11.5%); finish 7,856.
- Peak 889 = 840 syndromes + 24-bit carry + up to 25 helpers (special inverse returns).
- Top items: recover Syed's v13 delta (-3.0 to -3.6%), norm-preserving closings (-0.21 to -0.26%), cutoff 18
  (net -0.05 to -0.08%), boristown/GrantPrater patches in the generator.
- No circuit was built or verified in this step; all gains are estimates.

## Next

Pick the first dqi step from the ranked list (likely the `.kmb` reader/aligner and the v12/v13 diff).
