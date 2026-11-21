# Re-plan of 21 Nov: ThemBones source base, null corrections, LLM loop priced

- **Date:** 2026-11-21 11:45 UTC
- **Plan step:** re-plan (milestones M2.E–M2.H, steps S2.25–S2.32; Phase 3 milestones M3.1–M3.6 rewritten)
- **Files:** `Analysis/research/public_frontier_2026-11-21.md` and `Analysis/research/agent_loop_local_llm_2026-11-21.md`
  (written by two research agents); `Analysis/PLAN.md` (§1 snapshot, §7 re-plan of 21 Nov, §9 Phase 3, §12 risks);
  `Analysis/IDEAS.md` (59 status changes, 12 new rows, 21 Nov summary); `Analysis/research/invsqrt.md`, `dqi.md`,
  `chemistry.md` (§4: 4 corrections and 17 public nulls); `experiments/NULLS.md` (regenerated);
  `experiments/LEDGER.md` (E028, E029)

## What was done

- **Public frontier audit** (agent; notebooks and the 29 MB ThemBones generator dataset downloaded read only to
  `downloads/kernels-2026-11-21/`, nothing executed): methods of the 33.10 set per task, generator sources and
  licences, why our 19 Nov studies missed them, remaining headroom, a proposed plan.
- **LLM search loop study** (agent; web sources, nothing installed): cost of an evolve-style loop, the best
  open-weight model for this laptop (RTX 4060 Laptop 8 GB, i5-12450HX, 15.7 GB RAM), fine-tuning feasibility.
- Applied both reports' register proposals and rewrote the plan.

## Why

User request (21 Nov): resubmit at the public level, refactor the plan around the best public methods, price the
idea of a model-driven search loop, find the best local model and assess fine-tuning.

## Results and analysis

- **Every public gain is in-family:**
  - invsqrt 44 → 39 Q by measuring the divider digits out and re-making them in the undo, plus lent wires
    (LEND29); 39 → 37 Q for +27 Toffolis;
  - chemistry (8,143,44) → (9,127,46) at constant R·B (B ≤ 128 shrinks the inner temporary 303 → 287) and brot 7
    by staged rounding with refit (Q 427 → 405), loaders co-designed with the table (dead words as don't-cares);
  - dqi multi-output shared-operand field kernels (−6% ThemBones, −8% tonylica's D16 sites).
- **Four of our 19 Nov nulls were true only on a restricted domain:** the invsqrt width floor (stored digits),
  "Q pinned at 427" (B ≥ 129; the one-minute scorer grid that would have shown it was proposed but never run), brot 7
  (roundings of 8-bit fits), dqi operand sharing (single products). Recorded as corrections with their domains; the
  old rows stay. New rule in the plan's risks: a null states its domain, and a ceiling is checked cheaply before it
  is quoted.
- **Sources:** ThemBones has the only complete Apache-2.0 generator line with byte-for-byte rebuilds (dqi file A
  never rebuilt by them). Rithvik Redrouthu's 37 Q invsqrt has no public generator; tonylica builds from recovered
  KMX assembly (whether that counts as source for the prize is a question for the hosts, M3.6).
- **LLM loop:** one proposal costs about $0.004 (Haiku 5.5), $0.08 (Sonnet 5.5), $0.16 (Opus 5.5) via the API,
  billed outside the subscription; our evaluators screen 60–200 whole-task candidates per hour. Best local model:
  Qwen3.6-35B-A3B (UD-IQ3_XXS, 13.2 GB, est. 15–30 tok/s with MoE offload), no cheaper than Haiku per proposal and
  5–10× slower, so a fallback only. Fine-tuning EXCLUDED before the deadline (memory, no bf16 on P100/T4, too few
  accepted edits). Expected gain of the loop +0.3% to +1% score, P ≈ 0.25; one-day go/no-go (M3.4a) against a
  non-LLM search of equal time. Its proposed first target (K = 6 invsqrt loaders) is retired with the v12 base and
  must be re-chosen on the ThemBones base.
- **New plan (rest of Phase 2):** S2.25 submission (done), S2.26 reproduce ThemBones byte for byte, S2.27/S2.29/S2.31
  public edits as generator passes, S2.28 invsqrt 37 Q from source, S2.30 dqi 890 Q, S2.32 freeze.

## Next

S2.26: add numba to `~/venvs/gqai313` and rebuild the ThemBones invsqrt, dqi and chemistry files byte for byte.
