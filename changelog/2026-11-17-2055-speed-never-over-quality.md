# Speed never over correctness or result quality

- **Date:** 2026-11-17 20:55 UTC
- **Plan step:** none (working method), during S1.5
- **Files:** `CLAUDE.md` (rule "Optimize run time" reworded), `Analysis/PLAN.md` (section 5)

## What was done

Added the limit to the run-time rule: an optimized tool must reproduce the reference method's
results (same circuit bytes, costs and verdicts) on known cases before it replaces it; fast proxies
(static cost, reduced runs, sampling) only screen candidates and never replace the full check at a
gate; no search is cut short or approximated in a way that could miss a better circuit unless the
trade-off is stated and agreed with the user.

## Why

The user asked that run-time optimizations never degrade the quality of the results.

## Results and analysis

Applied already in S1.5: the full invsqrt verifier run was replaced by a reduced run only because the
generator's exact proof over all 3,133,452 inputs and all HMR transcripts already establishes
correctness; the full run stays mandatory before a submission is chosen as final.

## Next

Finish S1.5.
