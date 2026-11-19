# Public frontier as of 18 November: Asterios Panoras's "what we tried, what failed and why"

Written 2026-11-19. Source: the public notebook `asteriospanoras/qcc-28-23162-what-we-tried-what-failed-and-why`
(Asterios Panoras, board score 28.23162, updated 2026-11-18; 77 cells; pulled into
`downloads/kernels/asteriospanoras__qcc-28-23162-what-we-tried-what-failed-and-why/`, read as data, not executed).
Its earlier version *QCC 28.12769 | What Works, What Fails and Why* (same folder family) is a strict subset of it
and was already ingested into `Analysis/research/*.md` on 17-18 November; this note adds what is new (its sections 8
and "Open problems", the 18 November frontier files) and maps everything onto `Analysis/IDEAS.md` and
`experiments/NULLS.md`. Companion notebooks read for exact numbers: `maoxin123/qcc-28-23-4-022-toffoli-chemistry-loaders`
(version 5), `ozertuu/28-23163-sota-verified-compilation-top-portfolio`, and the HW-30 entry of
`manwithacat/vqcc-hedge-wizardry-notebook`. Numbers are the authors' own (official verifier, default salt, unless
marked). Our base for "gain beyond v12" is our rebuild of Syed Asad Ali's v12: invsqrt 80,162.73, dqi 3.6073e11
(default salt; 3.6132e11 mean of 8 salts), chemistry 296,635,967.59.

## 1. Summary

1. **The 28.23162 set (18 Nov)** = invsqrt `d0227f3c548ae7f1a409159a143cf4a401de56a2414a5230a3efc4b38e0a8160`
   (maoxin123, 44 Q, 1,451.44 mean Toffolis, **74,405.67**; Asterios prints 74,405.73), dqi
   `15a4b7403356092e0ef6cc61c9f18644864964daae9630575d77c29314f75c8b` (473,123,120 bytes; Ertuğrul Özer's merge, 890 Q,
   **365,894.20** mean Toffolis, **3.45257e11**, 0 failures at the default salt, clean on the board), chemistry = 15
   loaders at 4,022 Toffolis per pass (`angle_k_0` SHA prefix `ed648708`, from junmiaotiao's dataset
   `qcc-frontier-20261118w`, maoxin123's loaders on the `906fa220` coded-word table) + weights `4d63be3d` (leotang09,
   λ_eff 46.92156), T 14,706, Q 427, **294,642,170.57**.
2. Against our v12 base: invsqrt **−7.18%**, dqi **−4.29%** (−4.45% against our 8-salt mean), chemistry **−0.67%**:
   the whole 27.06 → 28.23 gap. Since then (19 Nov) maoxin123's invsqrt `ecf908fe` reaches 74,383.70 (−0.03%), and
   boristown's dqi `e39b1396` adds one more retained AND.
3. **Who built what:** invsqrt = Syed Asad Ali's v13 44-qubit move + Ali Haydar Özdağ's 64-row table + maoxin123's
   rescheduling chain (0af21cb8 → 563d0583 → 864f81cd → d0227f3c → ecf908fe) + small edits (Asterios's dead write
   and angle LP, jeffemerson's AND merge, boristown's f2cdab0b). dqi = Syed's v13 decoder + Balram Srinivasan's
   closings and converter (`40294f00`) + GrantPrater's merge and closings + boristown's PW4 retention + VUMANIAcom's
   retention (HW-30, merged by manwithacat) + Özer's final merge. Chemistry = Lera's coded-word line + maoxin123's
   4,022 loaders + successive weight refits (boristown, Atakan Aldemir, leotang09).
4. **Asterios's own measured gains are tiny** (dead write −0.058%, angle LP −0.0007%, 3 retained dqi products,
   a 62-qubit splice −0.20% later made obsolete); all the large steps are other people's compiled files.
5. **Main lessons:** (a) every circuit on the frontier is locally exhausted by the scans that exist (0 candidates for
   affine/constant/dead/cancelling Toffolis, angle LP at its bound, no retainable dqi pair left, no removable loader
   Toffoli); (b) the next 1% per task needs a new design, not a scan; (c) one failed dqi shot costs 0.13-0.15 points
   and the board draw depends on the file hash, so a changed dqi file is a gamble until drawn; (d) the input model
   is part of every claim (31 "dead" loader Toffolis were live because loader k starts from batch k−1).
6. **For us (track B, own source):** the most useful items are the open problems that *need the source build* and
   that nobody public can do on compiled files: the dqi one-qubit peak shave (DQI-A09, estimated +0.010 to +0.021
   points), and the exact per-adder pebbling DP (GEN-B03). Two ideas are new to our register: a cheaper
   GF(2^24) → GF(2^12) norm (DQI-N01) and the terminal readout-window rewrite as a generator pass (GEN-N01).

Exchange rates printed by the notebook (at the 18 Nov operating point): one invsqrt qubit ≈ 39 Toffoli-equivalents;
one dqi qubit ≈ 830-860 dqi Toffolis; 1% of one task's cost ≈ +0.1 points ≈ 18 invsqrt mean Toffolis ≈ 79 loader
Toffolis (per pass) ≈ 4,200 dqi Toffolis; one dqi Toffoli ≈ +0.00002 to +0.00003 points; one loader Toffoli
(−2 in T) ≈ +0.0013 points; one failed dqi shot ≈ +0.46% dqi cost ≈ −0.13 to −0.15 points.

## 2. Table A: techniques that worked publicly (beyond our v12 base)

"Gain" is the change of that task's cost against the file it was applied to (author's official numbers).
"Code" says whether a generator / source is public or only the compiled circuit.

| # | Technique | Task | Author | Measured gain (task cost) | Code public? | IDEAS ID |
|---|---|---|---|---|---|---|
| A1 | 64-row coefficient table (address 7 → 6 bits, loaders 56+70 → 27+29 ANDs, slope bit measured out) | invsqrt | Ali Haydar Özdağ (`4560ac0c`) | −1.74% on the 46 Q line (our measurement 80,162.73 → 78,765.64) | compiled; loader certificates from a private search | INV-A02 |
| A2 | Measure-first dead write: `CCX a b t` whose target is only measured next → `HMR t` first + `CZ a b if m` | invsqrt | Asterios Panoras (14 Nov, on `dcde7078`) | −1 Toffoli, 78,765.51 → 78,719.53 (−0.058%); already inside all 44 Q files | rule documented; detector in `qccscope` (Apache-2.0); circuit compiled | INV-A04 (N25: 0 sites in v12) |
| A3 | Syed's v13 46 → 44 qubits (root bits v10/v11 and flags e2/e3 off the wires) | invsqrt | Syed Asad Ali (`47b31838`) | 78,719.53 → 75,170.62 (−4.5%); Q −2, C 1,711.3 → 1,708.4 | compiled only (v13 source not public) | INV-A03 |
| A4 | Exact LP on the rotation angles (phase is linear in the angles; HiGHS, constraint generation, exact integer recheck) | invsqrt | Asterios Panoras (15 Nov, `09000960`); Ertuğrul Özer's angle refit in maoxin123's line | margin +0.16%, 75,170.62 → 75,170.13 (−0.0007%) | method described; LP code not published (only a margin recompute) | INV-A04, INV-A12 |
| A5 | Rescheduling / micro-edit chain on the 44 Q circuit, 34 → 33 rotations, AND merge cut | invsqrt | boristown `f2cdab0b` (75,006.38, −0.22%), maoxin123 `0af21cb8` (74,797.37, −0.28%), `563d0583` (74,560.91, −0.32%), `864f81cd` (74,490.95, −0.094%), `d0227f3c` (74,405.67, −0.11%); jeffemerson AND-merge (−2 T) | 75,170.13 → 74,405.67 = −1.02% in all | compiled, embedded in notebooks (maoxin123: Apache-2.0 / CC0); one-line mechanism notes only | INV-A04 |
| A6 | Division clean-up hands precomputed AND terms to the second table load (one more AND removed) | invsqrt | maoxin123 `ecf908fe` (19 Nov) | 1,451.4375 → 1,450.9375 mean T, 74,405.67 → 74,383.70 (−0.03%); margin 38,858,847 | compiled (embedded, Apache-2.0 / CC0) | INV-A23 / **INV-N01** |
| A7 | 62 Q splice of two public circuits at clean cuts (shortest path over sync points) | invsqrt | Asterios Panoras (`9c708f7b`, 12 Nov) | 1,216.83 → 1,213.87 T (−0.20%), then obsoleted by the 46 Q line | splice-search cell in the notebook (Apache-2.0) | none (NULLS N24, NR-0253/0254) |
| A8 | Extractor cutoff 18 + 241 retained products on Syed's v12 decoder | dqi | Ertuğrul Özer (cutoff), Areeb Ali (retention), `821bf8d1` | 3.602e11 vs our v12-at-15 3.6073e11 (−0.15%) | compiled; notebooks describe | DQI-A03 (DONE at 17), DQI-A06 |
| A9 | Terminal readout window: `CSWAP c a b` + measure both outputs → `CX a b; HMR a→m1; CZ c b if m1; HMR b→m2; BIT_INVERT m1 if m2` (0 Toffoli) | dqi | Atakan Aldemir (`e9f54be3`); transplant by Asterios onto `bb4287e7` | −24 T per window, −48.28 T/shot (−0.012%) on that older line | rule documented; circuits compiled | DQI-A05 / **GEN-N01** |
| A10 | Syed's v13 890 Q decoder (+1 qubit buys −13,900 T) | dqi | Syed Asad Ali (`46446af1`) | 3.602e11 → 3.479e11 (−3.4%) | compiled only | DQI-A01 |
| A11 | Norm-preserving GF(2^24) closing (12-bit norm kept, inverted via GF(2^6)), 41 sites | dqi | Balram Srinivasan (`66b4f18b`) | 3.479e11 → 3.471e11 (−0.23%, −820 T) | compiled; gadget described in notebook | DQI-A02 |
| A12 | Control merging, cached CSWAPs, measured flag uncomputations, product-cache chains | dqi | boristown (`ccc2ab84`, `00ac37d4`) | ≈ −0.2% | patch notebook (binary patches) | DQI-A04, DQI-A05 |
| A13 | Merge of patch sets, closing carried to 8 more sites (−152 T), 49 measured uncomputations of a 15-AND product (−49 T), 4 window-proven retained products | dqi | GrantPrater (`575b043f`) | board 28.15096 file | tools public in his dataset (`vqcc_*.py`) | DQI-A04, DQI-A10, DQI-A06 |
| A14 | Retained products "second kind": same three wires at both Toffolis, free wire over ~46k ops, `CX(t,b)` after the first and before the second (6 sites incl. 3 in condition blocks) | dqi | boristown PW4 (`5f901c55`); Asterios found the 3 condition-block ones independently (`327b0751`, −3 T) | −6 static T, 3.4640e11 (+0.00014 points) | patch notebook; Asterios's scan described | DQI-A06 |
| A15 | Sibling merge: align two files op by op, keep the cheaper side of every patch; 5 retention patches from VUMANIAcom's `2683225b` (HW-30) | dqi | VUMANIAcom; merged by manwithacat | −104.93 mean T, 48-draw expected cost −0.0286% | dataset `manwithacat/vqcc-hedge-wizardry` (Apache-2.0 + NOTICE) | DQI-A06 |
| A16 | Super-merge on Balram's `40294f00` (v13 + converter / phase-kick splice from an earlier variant, "CSWAP converter"): 32 GF(2^24) closings (−672 T), 8 MBUC (−8 T), PW4/retention (−44 T), HW-30 (−116 T) | dqi | Ertuğrul Özer (`15a4b740`), on Balram Srinivasan's `40294f00` | 3.46404e11 → 3.45257e11 (−0.33%); 365,894.20 T | compiled (dataset `ozertuu/quantum-sim-linalg-cache`) | DQI-A01, DQI-A02, **DQI-A10 (probable)**, DQI-A05, DQI-A06, **DQI-N02** |
| A17 | Coded-word table, 15 loaders (14,712 T total) + refit of her weights (λ_eff 46.92176) | chemistry | Lera (marculera) + boristown | 294.88M → 294.76M (−0.04%) | dataset `marculera/dfthc-c44-build` (Apache-2.0, NOTICE) | CHEM-A01 (EXCLUDED), CHEM-A03 |
| A18 | 4,022-Toffoli loaders on the XOR of consecutive batches (k0 292, k1 287, CNOT layers −64%) | chemistry | maoxin123 | T 14,712 → 14,706 (−0.04%) | compiled, embedded (CC0) | CHEM-A09, CHEM-A03 |
| A19 | Fixed-word weight refits (GRG / corank-2 reduced gradient, gates within ~1e-9 Ha) | chemistry | Atakan Aldemir `703843d3` (λ 46.92158), leotang09 `4d63be3d` (λ 46.92156, also rewrote two loaders), maoxin123 `56ef2562` | 703843d3 → 4d63be3d: −0.00004% | npz files public; Asterios's refit toolkit in the Atlas dataset | CHEM-A13 |

Note on A16: −672 T = 48 × 14, which matches GrantPrater's open item "14 more Toffolis per closing, waiting for a
second idle helper wire" (DQI-A10). This is an inference from the arithmetic, not stated by the authors; a
`kickmix-engine diff` of `5f901c55` against `15a4b740` would settle it. The "HW-30 (−116 T)" figure is a retention
patch, not a cutoff move (checked in the HW-30 entry), even though one cutoff step also removes 116 static Toffolis.

## 3. Table B: techniques tried and failed (nulls)

"Our row" proposes how the IDEAS entry should read. A null closes only its design family (IDEAS §9).

| # | Technique | Task | Author | Number that killed it | IDEAS ID | Our row |
|---|---|---|---|---|---|---|
| B1 | 45 qubits on the 46 Q line (rescheduler, in-place adders, rebuilds) | invsqrt | Asterios | peak at 46 in 411 runs (485 episodes before rescheduling, 261 after, 257 opened by an alternating-carry adder); break-even +38 T, cheapest exact cover +196 T, old attempt +257; hiding the rare shift flag (63 inputs) in affine parities: 0 candidates | INV-A24, GEN-B04 | stays open: the 44 Q peak has not been mapped by anyone (Asterios says the 46 Q analysis does not transfer); NR-0275 kills one 43 Q route |
| B2 | Drop / add / re-tune rotations (46 Q and 44 Q) | invsqrt | Asterios | 46 Q: 32 of 33 drops infeasible, the feasible one costs 79,369.55 vs 78,765.26; adding a 34th reaches margin 167M vs 332.6M needed. 563d0583: all 33 drops infeasible (break-even 10.0M units, best LP bound −55.7M); LP gains 2 units of 31,789,811 | INV-A12 (angle part), INV-A16 | angle-only LP: NULL (public); INV-A12 open only for the joint table-word MILP; INV-A16 needs margin that is not there (0af21cb8 uses 99.96% of its tolerance) |
| B3 | Merge two rotations / drop or flip the T gate | invsqrt | Asterios (on 563d0583) | no pair of free angles has equal or complementary firing bits (closest differs on 938 inputs); T gate fires on 54 inputs, dropping it LP bound −1.0e18, flipping to T† −2.3e18 | INV-A10 | NULL (public) on the frontier layout; recheck once on our own build (cheap) |
| B4 | Micro-edit scans (constant/equal/complementary controls, affine products, CCX → HMR, cancelling pairs) | invsqrt | Asterios (46 Q, 44 Q, 563d0583: 1,548 Toffolis, 16,000 sampled inputs × 2 seeds) | 0 candidates; rarest Toffolis fire on 1, 3, 7 inputs | INV-A04 (scan part), INV-A06, INV-A07, GEN-A03 | the frontier file is exhausted for these patterns; on our source they remain the porting checklist |
| B5 | Skip the recomputations of the backward half / pebbling by retention at 46 wires | invsqrt | Asterios | uncompute 791 vs 689 forward T (+102: 64 replays of a 13-wire divider register, ~69 conditional carry fix-ups); 490 keepable products each cross ≥ 100 moments with no free wire (median 846); exact packing keeps 0 at 46 Q and 2 at 47 | INV-A17, INV-A15 | open but weakened (each excess Toffoli is the price of a width saving; NR-0261 agrees) |
| B6 | 47 qubits with a third carry slot | invsqrt | Asterios | 141 of 142 conditional CCX compute values that never existed on a wire; ceiling 43.5 saved T vs 37.4 needed (hard ceiling 66.5); exact DP over the 31 adder groups not done | GEN-B03, GEN-A12 | open (already NULLS N8; the DP is the open part) |
| B7 | Splice of 8 public 46 Q variants | invsqrt | Asterios | best path = base file (1,716.259 T-eq) | none | NULL (public); irrelevant to track B |
| B8 | Choose between equivalent files by the hash-seeded Toffoli mean | invsqrt | Asterios (declined) | seed noise ≈ ±0.03 in cost | none | excluded on principle (gaming the verifier's randomness) |
| B9 | Cheaper field multipliers: the "rank-126 region → rank-99 GF(2^24) product" | dqi | Asterios | the region is three GF(2^12) rank-42 products (a·x, b·y, (a+b)(x+y)) with 36 outputs; none of the 12 runs of 126 CCX has output rank 24; 0 Toffolis | DQI-B08 (unaffected), DQI-A19 | the 126 → 99 sub-idea is NULL (public); DQI-A19 inherits the open three-product question (Table C) |
| B10 | Keep products across the return path instead of recomputing | dqi | Asterios (889 Q) | 620 GF(2^24) products (61,380 T) belong to an in-place three-term recurrence (bijection, no garbage); keeping one needs 24 qubits where no wire is free; the 82 repeated 12-bit norms (3,444 T) need +12 qubits (≈ −10k T in value) | DQI-B02, DQI-A06 | family-specific null (NR-0184, NR-0201 agree); B02 stays open as a redesign |
| B11 | Skip rare routing gates that look like no-ops | dqi | Asterios | 9,724 of 12,240 CSWAPs never swap on 1,024 shots, but on 65,536 official-distribution shots 948 of 1,020 routing slots fire (≈ 1 per 4,096 frames); dropping one fails ≈ 2.5e-4 of shots; 16,369 "never firing" gates on the 890 Q file are sample artifacts | DQI-A07 | the "drop rare gates" variant is NULL (public), as NR-0170/NR-0187; A07 stays open only as exact pricing of whole handlers |
| B12 | Extractor cutoff 19 and beyond | dqi | Asterios | −116 static T per step vs +3.46e-6 failures per shot (0.35 per draw, ≈ −0.041 points per failure) | DQI-A03 | already a null (dqi.md); our choice of 17 stands |
| B13 | Repair the rare decoder failures cheaply | dqi | Asterios (on `7cbe8ecf`, cutoff 15) | 6 failures in 25 salts × 102,400 shots (2.3e-6/shot), 128 of 175 flagged shots fail deterministically; no cheap exact fix found; expected loss ≈ 0.019 points | DQI-A17, DQI-A07 | open (exact handling is listed as an open problem) |
| B14 | More retained products on the frontier file (symbol pairs, condition blocks, dirty helpers, cross-block pairs, "second kind") | dqi | Asterios (on `575b043f`, 57.5M ops) | 3,978 unconditioned pairs: none has a free wire over its span (median span 53,714 ops); 206 in-block pairs → 3 usable; dirty helpers 0/4,184 and 0/606; second kind: 279,577 consecutive identical-wire pairs → 24 hold one product → the same 3 as PW4 | DQI-A06 | the frontier file is exhausted under these definitions; for our source build A06 remains the port of all public retention (and new pairs may appear in a different schedule) |
| B15 | 0 cancelling Toffoli/CCZ/CSWAP pairs among 58.4M ops of the 890 Q file | dqi | Asterios | 0 | GEN-A03 | expectation ≈ 0 on frontier-class files |
| B16 | Re-plan the loaders (block matcher, greedy, local and flip-graph search) | chemistry | Asterios | public loaders 212 products (193 for loader 4) vs ≈ 210 counting estimate; block matcher 229, greedy 240-250; removing any product or pair leaves the target uncovered | CHEM-A12, CHEM-A11 | NULL (public) for greedy/local families (NR-0268 agrees); CHEM-A12 open only as exact ILP/SAT with fully mixed products |
| B17 | Delete 31 "dead" loader Toffolis (dead when the angle registers start at zero) | chemistry | Asterios | rejected by the official verifier ("uncorrected phase"): loader k starts from the words of batch k−1 | CHEM-A11, GEN-A10 | NULL (public) for the dead-gate idea; a reminder that any loader edit must model z = previous batch |
| B18 | Smaller shape at fixed words: C 44 → 43; R 8 → 7 | chemistry | Asterios | C = 44 two-body floor 19.797 mHa (table 0.20 mHa above it); C = 43 floor 21.10 mHa > 20; R = 7 reaches 256 mHa | CHEM-B04, CHEM-B10 | at fixed words NULL (public), as NR-0210; B04 (retrained words) and B10 (R 7 with C 48-50) stay open |
| B19 | Continue the weight refit (GRG from boristown's and Atakan's weights; augmented Lagrangian) | chemistry | Asterios | λ_eff unchanged (46.92176; 46.92158); AL worse (46.9107 vs 46.9082, 46.9108 vs 46.9082, 46.9146 vs 46.9092); continuing from 703843d3 changes λ by ≈ 1e-9; 0bd9c68f (46.92158430) worse than 703843d3 (46.92158337) | CHEM-A13, CHEM-B18 | NULL (public) for local refits at fixed words; multi-start untested but capped at ≈ 1e-6 relative |
| B20 | Remove loader Toffolis / shrink the cache / lower Q | chemistry | Asterios (15 loaders, 4,022 per pass) | no removable Toffoli; cache at its minimum; Q fixed by the cost model | CHEM-A10, CHEM-A11, CHEM-A06 | A11 NULL (public); A10 and A06 (Q) unaffected only if they change the construction or the cost-model point |
| B21 | Re-submit a changed dqi file for the board draw | dqi | Asterios (board table) | of 5 files drawn, one lost a shot (`f425ea20`: 28.0884); `e208d350` (365,890 T, cheaper than 15a4b740) scored 28.18862 after a lost shot | GEN-A04, DQI-A11 | evidence for the multi-salt gate; the private draw is independent of the public one |

## 4. Table C: open problems they list, mapped to our register

| Open problem (their words, condensed) | Their estimate | IDEAS ID | Comment for us |
|---|---|---|---|
| One qubit less on the 44 Q invsqrt; peak of the 44 Q circuits not mapped | 1 qubit ≈ 39 T-eq ≈ 2.3% | INV-A24, GEN-B04, GEN-B08 | We have the generator; a liveness map of our own 44 Q build (once INV-A03 is done) is the prerequisite |
| Exact per-adder pebbling DP over the 31 adder groups of the 46 Q design | best case −0.37% | GEN-B03 | Ours to do at source level; small, exact, never done publicly |
| A better invsqrt approximation with more phase margin (table or interpolation design) | each margin doubling ≈ 4.6 T-eq (33 rotations); one rotation fewer ≈ 7.3 T-eq | INV-B01..B09, GEN-B12, INV-A12, INV-A18, INV-A19 | "angle tuning is exhausted": all margin must come from the arithmetic |
| dqi peak: 509 blocking episodes on the 889 Q file in seven families, one-wire reallocation each; "needs the decoder's source build" | +0.010 to +0.021 points | DQI-A09 | The strongest track-B opportunity in the notebook: public teams cannot do it on compiled files |
| Cheaper three-product map (a·x, b·y, (a+b)(x+y)), i.e. multiplicative complexity below 3 × 42 | 126-T region, small | DQI-A19 | Karatsuba-shaped triple; rank additivity likely, low P |
| GF(2^24) → GF(2^12) norm computed with a rank-42 product; 42 → 36 | ≈ 1.2k T, ≈ 0.3% (estimate) | **DQI-N01 (new)** | Not in our register; see below |
| An exact dqi decoder (abnormal Euclid steps handled exactly) | worth ≈ 0.019 points in expectation | DQI-A17, DQI-A07 | Pays only if cheaper than ≈ 0.019 points (≈ 1,500-2,000 T) |
| Chemistry optimum beyond the local one (different start, different table, global method) | small | CHEM-A13, CHEM-B18 | Capped by B19; low priority |
| Rebuild at C = 40 (Lera's (9, 143, C) builds suggest it) | −1.8% cost, ≈ +0.16 points | CHEM-B12 | Needs new factorization, words and loaders |
| Fully mixed loader products: counting floor ≈ 187 vs 212 per loader | ≈ −25 products per loader | CHEM-A12, CHEM-B13 | Nobody has a construction; exact ILP/SAT on one loader is the cheap test |

New ideas not in our register:

| Proposed ID | Idea | Evidence | Estimated gain | Notes |
|---|---|---|---|---|
| INV-N01 | Hand precomputed AND terms from the division clean-up to the second table load (and drop the duplicated AND), as in maoxin123's `ecf908fe` | maoxin123 v5, 19 Nov: −0.5 mean T, −0.03% | 0.03% (more if other stage boundaries allow it) | Close to INV-A23 but across the divider/LOAD2 boundary; implement in the generator once A02/A03 are in |
| DQI-N01 | Cheaper GF(2^24) → GF(2^12) norm. In the tower basis z = u + v·t, N(z) = u² + u·v + α·v²: the squares are free CX, so the norm is a *quadratic* map whose only cost is u·v (42 ANDs today). Search for a lower AND count of the 12 quadratic output forms (quadratic, not bilinear, complexity; products of linear forms in u and v together) | Asterios's open problem (42 → 36 estimate) | ≈ 0.3% task cost if 36 is reachable at every norm site (≈ 1.2k T) | Our dqi generator computes norms in `gf24_nested_norm` and in the K24 path; check the number of norm calls in v12 first; a 576-pair exhaustive check proves a candidate in seconds |
| DQI-N02 | Port Balram Srinivasan's `40294f00` converter / phase-kick splice ("CSWAP converter") into the v12 generator | base of the 18 Nov frontier; mechanism not described in the notebooks read here | part of the −0.33% of `15a4b740` vs `5f901c55` | First read `balramsrinivasan/vqcc-cand-28199-submit` (not downloaded); may overlap DQI-A01 |
| GEN-N01 | Terminal readout-window rewrite as a generator pass: any `CSWAP c a b` whose two outputs are only measured afterwards becomes `CX a b; HMR a→m1; CZ c b if m1; HMR b→m2; BIT_INVERT m1 if m2` (0 Toffoli, same record and phase) | Atakan Aldemir's windows, −24 T per window, transplant −48.28 T/shot (Asterios) | small per site; free where it applies | Add to the GEN-A11 matcher and to the finishing pipeline (GEN-A03); applies to any task with CSWAP + measurement |

## 5. Their method and tools

**Procedure** (their section 4, "profile first, bound second, search last"): (1) exchange rates in points for one
qubit, Toffoli, rotation and bit of phase margin; (2) census and attribution of every Toffoli to a gadget, descend only
into the dominant one; (3) width profile: peak episodes, which wires open them, which sit idle; (4) distance to a
bound for each gadget (bilinear rank, n−1 for an adder, LP-optimal angles, counting floor) — a gadget at its bound is
closed; (5) Amdahl ceiling of each lever; (6) motif matching and value fingerprinting (repeated, constant, affine,
dead Toffolis); (7) check against recorded negatives; (8) cheapest falsification first, count nothing until the
official verifier passes. Two practical traps they document: a peephole scan that skipped input-register wires
missed the dead write; and a splice that looked 0.228 T-eq cheaper on 512 shots was only 0.109 on 4,098 (conditional
blocks fire on random outcomes: re-measure every candidate on a large fresh sample).

**Pattern catalog** (motif → best known → status on the files): temporary AND (1 T, erase by HMR + CZ; 0 left
anywhere), dead write (0 T; fixed), constant/equal controls (0 in exhaustive invsqrt scans), repeated product
(retention; 479 in the 46 Q file, 0 fit), n-bit adder (n−1 T; sets the invsqrt peak), GF(2^k) product (rank 42 for
k = 12; at the bound in dqi), QROM loader (counting floor of the decomposition; at the floor of its method),
rotations (LP-optimal; Hamming-weight phasing for equal angles).

**Tools.**

| Tool | What it is | Public? Licence | Reusable for us? |
|---|---|---|---|
| `kmsim` | ~150-line numpy simulator of the kickmix subset (bool arrays, one column per shot, phase in 2^-64 turn, hook after any op) plus `check()` with the invsqrt verifier rules | embedded in the notebook (cell 9), Apache-2.0 per the notebook | Nothing new for us: `kickmix-engine`, Syed's bit-sliced `kmsim.py` / `prove_any.py` and our `tools/kmb.py` cover it and are faster |
| `qccscope` v0.1 | ~1,650 lines pure Python/numpy, 9 modules (`kmx` streaming parser to a numpy op table, `cost` exchange rates in points, `liveness` width profile / peak episodes / idle-gap cover, `structure` static motif detectors, `values` bit-sliced sampled simulation and product fingerprints, `motifs` catalog with bilinear-rank bounds, `report` markdown); exhaustive invsqrt value scan in seconds; chunked for the 53M-op dqi | embedded as a base64 zip in cell 30 (28,232 bytes, SHA-256 `29df0ad66e7992cfbcf3e49af28032d1d6212f721733508535d014801114225b`), Apache-2.0 per the notebook | Useful as a **design reference** for GEN-A08 (peak map), GEN-A11 (motif matcher) and the dead-write / repeated-product detectors. Repo rule: downloaded code is untrusted, read as documentation and re-implement in `tools/` (it also hard-codes a stale `GM_PUBLIC_BEST = 27.29`) |
| QCC Circuit Atlas (`asteriospanoras/qcc-circuit-atlas`) | dataset: index of public circuits (`circuits.csv`), `.kmx` copies, `best/` files, `tools/refit` (GRG2 weight-refit model after prvsiyan's method, `atlas_io`) | public Kaggle dataset; licence of the dataset itself not checked here | The refit toolkit duplicates CHEM-A13 work that is now shown to be at a local optimum; low value |
| Splice search | shortest path over clean cuts (empty condition stack, no classical bit crossing) of public circuits | appendix cell, Apache-2.0 | Not for track B (we do not submit third-party files) |

## 6. Proposed status per IDEAS ID

Proposals only; `Analysis/IDEAS.md` is not edited here.

| ID | Proposed status | Justification |
|---|---|---|
| INV-A01 | EXCLUDED (unchanged) | Track B decision; the public file is now `d0227f3c` (74,405.67) or `ecf908fe` (74,383.70) |
| INV-A02 | STILL OPEN (PAUSED) | −1.74% confirmed publicly; all frontier files carry it; compiled only |
| INV-A03 | STILL OPEN (PAUSED) | −4.5% confirmed (75,170.62); v13 source not public |
| INV-A04 | STILL OPEN | The public chain is now 75,170.13 → 74,383.70 (−1.05%); scans on the frontier file find 0 more, so the value of A04 is the port, not new edits |
| INV-A06, INV-A07 | STILL OPEN (lower P) | 0 candidates in all public scans of 46 Q / 44 Q files |
| INV-A10 | NULL (public: no mergeable rotation pair, T gate cannot be dropped or flipped, on 563d0583) | Recheck once on our own build, 1 h |
| INV-A12 | STILL OPEN (angle part NULL (public): LP bound equals the file, +2 units) | Only the joint table-word + angle MILP remains |
| INV-A16 | STILL OPEN (lower P) | No margin left (99.96% of tolerance used; every rotation drop infeasible) |
| INV-A17, INV-A15 | STILL OPEN (lower P) | 46 Q: every excess uncompute Toffoli pays for width; retention keeps 0 products at 46 Q |
| INV-A23 | STILL OPEN | Related public evidence in `ecf908fe` (see INV-N01) |
| INV-A24 | STILL OPEN | The 44 Q peak is unmapped publicly; the 46 Q analogue costs +196 T vs 38 |
| INV-N01 | NEW, TODO | Clean-up → LOAD2 AND hand-off, −0.03% public |
| GEN-B03 | STILL OPEN | Exact per-adder pebbling DP named as open; best case −0.37% |
| GEN-B04, GEN-B08 | STILL OPEN | One invsqrt qubit ≈ 39 T-eq; no public cover at 44 → 43 |
| GEN-B12, INV-B* | STILL OPEN | "A better invsqrt needs a different approximation or arithmetic design" |
| DQI-A01 | STILL OPEN | −3.4% confirmed; compiled only |
| DQI-A02 | STILL OPEN (DEFERRED) | Public value grew: −820 (Balram) −152 (GrantPrater) and probably −672 (A16) T on closings |
| DQI-A03 | DONE (ours, cutoff 17) | Public confirms ≥ 19 is negative; the 889 Q line was at 18 (the cutoff of the v13-based `15a4b740` is not stated) |
| DQI-A04, DQI-A05 | STILL OPEN | Ports of boristown's patches, GrantPrater's 49 + Özer's 8 measured uncomputations, readout windows |
| DQI-A06 | STILL OPEN (port only) | The frontier file is exhausted by three scan definitions; the public total to port is Areeb 241 + GrantPrater 4 + PW4 6 + HW-30 5 (+ boristown's 1 of 19 Nov) |
| DQI-A07 | STILL OPEN (drop-rare-gate variant NULL (public)) | 948 of 1,020 rare routing slots fire; dropping one fails ≈ 2.5e-4 of shots |
| DQI-A09 | STILL OPEN (raise P) | Asterios's 7-family plan (+0.010 to +0.021 points) needs exactly what we have: the source build |
| DQI-A10 | DONE (public, probably by Ertuğrul Özer in `15a4b740`: −672 T = 48 × 14) | Verify by a diff before recording; still to port into our generator |
| DQI-A11 | STILL OPEN | Board table: 1 of 5 files lost a shot; `e208d350` lost 0.04 points against `15a4b740` |
| DQI-A17 | STILL OPEN | Exact decoder worth ≈ 0.019 points; no cheap fix found publicly |
| DQI-A19 | STILL OPEN (target refined) | The 126-T region is the triple (a·x, b·y, (a+b)(x+y)); its complexity below 126 is open |
| DQI-B02 | STILL OPEN | Return-path retention is null at 889 Q (no free wire), so only a redesign can help |
| DQI-B08 | STILL OPEN | Unaffected: the 126 → 99 sub-case does not exist, the real rank-99 sites remain |
| DQI-N01 | NEW, TODO | Norm map 42 → 36 ANDs, ≈ 0.3% if reachable |
| DQI-N02 | NEW, TODO (read first) | Balram's `40294f00` converter / phase-kick splice |
| CHEM-A01 | EXCLUDED (unchanged) | Frontier parts are now the 4,022 loaders + `4d63be3d` |
| CHEM-A02 | STILL OPEN | Weights moved three times in four days, each by < 1e-6 relative |
| CHEM-A03, CHEM-A09 | STILL OPEN | maoxin123's 4,022 loaders (k0 292, k1 287) show the target per loader |
| CHEM-A11 | NULL (public: no removable Toffoli in the 15 loaders; "dead" Toffolis are live under the batch k−1 input model) | Reopen only for loaders we synthesize ourselves |
| CHEM-A12 | STILL OPEN (greedy / local / flip-graph variants NULL (public)) | Floor ≈ 187 vs 212; only an exact or fully mixed construction is untested |
| CHEM-A13 | NULL (public: local optimum from three starts, AL worse, leotang09's gain 4.4e-7 relative) | Multi-start is untested but capped far below one loader Toffoli |
| CHEM-B04 | STILL OPEN (fixed-word variant NULL (public): floor 21.10 mHa > 20) | Only a retrained table at C = 43 remains |
| CHEM-B10 | STILL OPEN | R = 7 at C = 44 fails (256 mHa); the idea is R = 7 with C = 48-50 |
| CHEM-B12 | STILL OPEN | Their estimate −1.8% (+0.16 points), untried |
| CHEM-B13 | STILL OPEN | Matches their "fully mixed loader products" open problem |
| GEN-A03 | STILL OPEN (lower expectation) | 0 cancelling pairs in 58.4M dqi ops; 0 invsqrt micro-edits |
| GEN-A04 | IN PROGRESS (supported) | Board-draw table: a changed dqi file can lose 0.13-0.15 points |
| GEN-A08, GEN-A11 | STILL OPEN | `qccscope` (liveness, motifs) is a good reference design; re-implement, do not execute |
| GEN-A10 | STILL OPEN | The 31 "dead" loader Toffolis show why every rewrite needs a proof under the true input model |
| GEN-N01 | NEW, TODO | Terminal readout-window rewrite as a generator pass |
