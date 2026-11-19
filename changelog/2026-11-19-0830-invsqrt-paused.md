# invsqrt paused for Phase 2 (S2.9–S2.13); Phase 2 re-oriented to chemistry and dqi; own invsqrt generator in Phase 3

- **Date:** 2026-11-19 08:30 UTC
- **Plan step:** S2.9 (stopped at the analysis), S2.10–S2.13 PAUSED; new steps S2.14–S2.16; new milestone M3.4b
- **Files:** `Analysis/PLAN.md`, `Analysis/IDEAS.md` (INV-A02, INV-A03 PAUSED), `experiments/LEDGER.md` (E019)

## What was done

Started S2.9 (integrating our certified 64-row table into Syed Asad Ali's v12 invsqrt generator) by mapping the
build: `tools/build_variant.py` replays candidate10's generator and applies 13 patch modules
(`src/patches/PATCH_ORDER.txt`) that bring the circuit to 46 qubits.

## Why we stop here (option A rejected, user decision 2026-11-19)

The 46-qubit width is not a property of the algorithm alone; it is built from facts that are **specific to the
128-row table** and were measured or searched for it:

- **Bit identities of the table words.** `patch_coef2_q46` uses "H4 = D12 ⊕ H3 ⊕ D11·H3" (from the pinned
  `slope_h3_h4_retirement_v3.json`) to keep D12 off the wires in the second coefficient pass.
- **Full-domain measured facts.** `patch_h2ret` (retire the second-table metadata bit H2 and rebuild it with 3
  Toffolis), `patch_seedinv` (recover raw root bits from the seed rows), `patch_w47_body` (early X-measurement of
  a6) and `patch_hret_sat_s15k7a` (SAT-minimal restore programs) were each validated on all inputs with the
  128-row coefficients.
- **Frozen certificates.** Both table loaders (`first_cap20_56b8.json`, `second_cap19_70nb1.json`, work caps 20
  and 19 under the 46 peak) and the 21-Toffoli root-value decoder are search results for those 128 rows, pinned
  by the model's SHA-256; the search tools that produced them are not in the v12 source.
- **The exchange rate makes a partial port worthless.** Without these tricks the circuit no longer fits in 46
  qubits; one more qubit (47) costs +2.2% of the task, more than the whole gain of the 64-row table (−1.7%,
  Ali Haydar Özdağ's public result).

Re-deriving each trick for a new table is several days of reverse engineering with a high risk of not reaching
46 qubits again, for at most +0.6% of score. INV-A03 (44 qubits) sits on the same stack. The S2.7 work is kept:
the certified 64-row table (refit margin 115.7M units, 0 failures) and the tools `team/table_search.py` (exact
integer model, table search, angle refit) and `team/loader_synth.py`.

## New direction (option B now, option C for January)

- **Phase 2:** tasks whose generators we can change safely. S2.14–S2.15: CHEM-A03 (loaders for Lera's
  (8,143,44) coded table compiled by our copy of Syed Asad Ali's chemistry compiler, ≈ −0.7% task cost), then
  CHEM-A06/A08. S2.16: DQI-A07 (the ≈ 0.35 failed shots per draw of the v12 decoder line: a guard would cut both
  the expected cost and the private-draw risk).
- **Phase 3, milestone M3.4b:** an invsqrt generator of our own, written so that we understand and control every
  width trick, starting from the S2.7 tools and the best B ideas; the route to the large invsqrt gains.

## Results and analysis

No circuit change. Lesson for the plan: before porting an idea into a third-party generator, check whether the
target structure (here the table) is an input the generator is built to vary, or a constant its optimisations
were tuned on.

## Next

S2.14: CHEM-A03 part 1.
