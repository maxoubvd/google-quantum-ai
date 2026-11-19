# Team changes to the invsqrt generator

This folder is our working copy of Syed Asad Ali's v12 invsqrt generator (Apache-2.0; `LICENSE` and `NOTICE` kept
verbatim; Syed Asad Ali is the original builder of this generator line and is credited first). The pristine source
is `third_party/syed3000-vqcc-own-circuit-portfolio-v12/invsqrt/invsqrt_source/`. No file of the original tree is
edited; our additions live in `team/`.

| File | What it does | Added |
|---|---|---|
| `team/table_search.py` | Coefficient-table search and certification (INV-A02): exact re-statement of the generator's integer model (`whole_model.py`) parameterised by the table; regression check against the frozen K = 7 model (`--check-k7`), K = 6 row search (`--search`), angle-refit feasibility (`--refit`) | 2026-11-18 (S2.7) |
| `team/models/model_k6_v1.json` | 64-row table (13-bit intercept A, 7-bit slope S, 6 low bits), with D, C, U0 for all 4,096 roots; SHA-256 `eac13e49…5981` | 2026-11-18 (S2.7) |
| `team/models/refit_k6_v1.json` | Angle-refit feasibility on that table: 34 angle deltas, exact minimum margin 115,688,300 units, 0 failures | 2026-11-18 (S2.7) |
| `team/models/k6_refine.json` | Exhaustive local search on the 12 tightest rows (confirms the coordinate search optimum) | 2026-11-18 (S2.7) |
