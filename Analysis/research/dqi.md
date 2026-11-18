# dqi: research notes and prioritized ideas

Written 2026-11-17. Scope: the `dqi` task (reversible RS(4095, 4025) syndrome decoder over GF(2^12) for DQI).
Base: our rebuild of Syed Asad Ali's v12 generator, `artifacts/v12/dqi.kmb` (`7cbe8ecf`), 889 qubits,
383,156.3 mean Toffolis (384,003 static), cost 3.6073e11 at 0 failures, 3.6132e11 mean over 8 salts.
Public frontier: GrantPrater `575b043f` (890 Q, ~367,116 T, 3.4641e11 on a clean draw), i.e. 3.97% below us.

Exchange rates used throughout (cost = R(eps+) x Q^2 x T, baseline T = 383,156):

| Lever | Task cost | In Toffolis |
|---|---|---|
| -1,000 mean Toffolis | -0.261% | 1,000 |
| -1 qubit (889 -> 888) | -0.225% | ~862 |
| +1 expected failed shot per 102,400-shot draw | +0.46% (first), +0.39% (second) | ~1,760 |
| failure probability p per shot | +0.046% per 1e-6 | p x 1.8e8 |

A 10% task-cost cut is +3.57% on the score; the whole task is 1/3 of the log-score.

## 1. Summary

1. Our v12 decoder costs 3.6073e11 (889 Q, 383,156 T); the public frontier is 3.97% cheaper (890 Q, ~367,116 T),
   built on Syed Asad Ali's v13 file, whose source is not public.
2. Anatomy (v12 ledger): key stage 56.5% (forward 101k, reverse 116k), of which 30.7% is 1,190 GF(2^24)
   multiply-accumulates at 99 T in the Thiele/CRT extractor; convert (Forney half-Euclid) 29.8%; scan 11.5%
   (~11 T per position, near its floor); projective finish 2%.
3. Width: 840 + 24-bit carry + up to 25 helpers; the 889 peak is in the extractor's special inverse returns.
4. Failures: ~0.37 per draw (+0.17% expected), mostly weight-30-35 "zero word at a step" classes; cutoff 18 is
   the optimum of the cutoff trade.
5. Best odds: recover the v13 delta into our generator (A01, -3.0 to -3.6%), then port the public closings,
   cutoff, MBU and retention patches (A02-A06, A10, A12): about -4 to -4.5% in total, from source.
6. Beyond the frontier: failure-budget audit (A07), width/Toffoli Pareto and peak shaving (A08, A09).
7. Structural bets (B): reverse-key redesign (B02), CRT form of the convert (B04), explicit-synthesis queries
   (B09, blocked by width so far), GF(2^24) rank below 99 (B08), explored with an agentic loop (B01).
8. 31 ideas (19 in A, 12 in B), each sorted by P x gain / effort; nulls in section 4 to avoid repeats.
9. Shared prerequisites: our own `.kmb` reader/aligner and liveness profiler, and kickmix-engine.
10. Matching the 31.4 leaders needs ~28% on one task; on dqi only category B could supply that.

## 2. Cost anatomy of the v12 decoder

### How it was obtained

- **Stage ledger from the generator itself.** The v12 sources keep an exact static-Toffoli ledger per component:
  `src/controller_work/after_v3/integration_v6/SOURCE_MODEL_LEDGER.json` (v6, 383,992 static T) plus the v7 deltas
  hard-coded in `integration_v7/v7_pipeline.py::source_ledger` (v7 = 384,003 static, +11 for the adaptive gauge).
  The sub-split of the extractor comes from `night_c27/arithmetic/raw_query_extractor.py::source_ledger`
  (formula `99 n^2 + 951 n - 198` for the arithmetic, n = 35 slots), which predates the last v4-v7 trims, so its
  shares are approximate (model 157k against 151.5k shipped).
- **Pipeline order** from `night_c27/synthesis/fixed_pipeline.py::emit` (the v7 builder rebinds the providers but
  keeps this order): mask gauge open -> public CRT map -> extractor prefix -> tail gauge -> scalar seed ->
  synthesis prefix (34 nodes) -> outer/convert prefix -> projective finish prefix -> scan -> finish closing ->
  outer/convert closing -> key closing (synthesis returns, seed, gauge, extractor return) -> public map inverse ->
  mask gauge close.
- **Executed counts** from our default-salt verifier run (`artifacts/v12-verify/dqi_default.out`).
- **Width** from `integration_v5b/v6/v7 WIDTH_COMPATIBILITY.md` and the extractor ledger. No full liveness
  profile of the 889-qubit tape exists yet (idea DQI-A08 builds one).

### What the algorithm is (in plain words)

1. **Mask gauge (classical, free).** The phase `(-1)^(sum_j <e_j, M_j>)` is unchanged if the masks `M` are shifted by
   any vector of the dual (RS) code, provided the linear phase `<s, y>` on the syndromes is applied (CZ/Z only).
   The 70 GF(2^12) degrees of freedom zero the effective mask at 70 fixed positions (1..70); an adaptive classical
   selector over 64 candidate gauges (322 classical bits, 0 qubits) zeroes a few more on average (-21.65 T mean).
   Classical-bit logic costs no Toffolis.
2. **Public CRT map (CNOT only).** The 70 syndromes are mapped linearly to 35 residues modulo 35 public quadratic
   node polynomials; each residue is a GF(2^24) element ("head", 24 qubits). 35 x 24 = 840.
3. **Extractor ("key", forward).** A Thiele-type continued-fraction / rational interpolation over the 35 GF(2^24)
   nodes, in place: for slot i, the head is updated through all earlier pivots by `query_step`
   `(n, d) -> ((d + P_l(alpha_i) n) / q_l(alpha_i), n)` = one GF(2^24) multiply-accumulate (99 Toffolis, bilinear
   rank 99) plus public scalings (free); then a deferred GF(2^24) inverse and a normalizer product into the
   24-bit carry bank. A rolling "tail" controller handles error weights below 35 from slot `j_min` = 15 on
   (cutoff 17/18 in the public forks); weights below `j_min` are not decoded.
4. **Synthesis** (34 nodes): rebuilds the numerator/denominator (evaluator and locator) polynomials of degree
   35/36 over GF(2^12) from the continued-fraction coefficients ("fused CRT synthesis").
5. **Convert / preimage (Forney).** Locator `sigma = E^2 + x O^2` (char 2: `sigma' = O^2`). A half-size extended
   Euclid between O and E (`raw_gap2_bezout`, degree gaps up to 2 covered) gives the polynomial
   `R = Omega / O^2 mod sigma`, whose value at each root is the error value. 
6. **Projective finish** and **scan**: an additive FFT (public linear maps, free) evaluates sigma and R at all
   field points; at each of the ~4,024 positions not zeroed by the gauge a 12-bit zero test of sigma(X) and a
   phase `Tr(M'_j R(X))` cost 11 Toffoli-likes; the zero tests are uncomputed by measurement (free).
7. Everything is uncomputed in reverse order (closings), mostly with measurement-based uncomputation (HMR plus
   Clifford or cheap quadratic phase repayments, "masks" kept in classical bits).

### Toffoli budget per stage (static Toffoli-likes)

| Stage | Forward | Reverse | Pair | Share of 384,003 | Notes |
|---|---:|---:|---:|---:|---|
| Key = extractor + synthesis | 101,087 | 115,841 | 216,928 | 56.5% | reverse is 14,754 T (+14.6%) dearer than forward |
| - extractor (rolling gauge, padded seed) | | | 151,531 | 39.5% | v6 151,852 minus 321 in v7 |
| -- query MACs, 1,190 x GF(2^24) at 99 T | 58,905 | 58,905 | 117,810 | 30.7% | the single largest item; O(n^2) in the 35 slots |
| -- consuming scales (99 x (4n-2)) | | | ~13,662 | 3.6% | model |
| -- deferred K24 inverses (390 per slot pair) | | | ~13,650 | 3.6% | model; matched 191/195-T inverses on the first 3 charts |
| -- affine quotient phases (264 per slot) | | | ~9,240 | 2.4% | model |
| -- rolling tail controller (138 per slot above j_min, -45) | | | ~2,715 | 0.7% | 138 T per cutoff step |
| -- seed, gauge, special charts, trims | | | rest | | |
| - synthesis (34 nodes) | | | 65,076 | 16.9% | |
| Convert / Forney preimage (outer, odd/even half-Euclid) | 57,034 | ~57,389 | 114,423 | 29.8% | v6 115,025 minus 602 in v7 |
| Scan (additive FFT + zero tests) | 44,286 | 0 | 44,286 | 11.5% | ~4,024 positions x 11 T; one pass only |
| Projective finish | | | 7,856 | 2.0% | |
| Seam delta | | | -92 | | |
| **Total static** | | | **384,003** | 100% | mean executed 383,156 (847 T skipped by measured conditions) |

Executed per shot (default salt): CCX 363,794, CSWAP 12,384, CCZ 6,978 (= 383,156 Toffoli-likes), HMR 70,959,
CX 31.1 M, SWAP 266k, `num_ops` 52.9 M, reaction depth 41,048, classical bits 53,911.

### Width (Q = 889)

- 840 syndrome/head qubits (incompressible: they are the input and are restored).
- 24-bit carry bank (running pivot reciprocal of the extractor), so 864 are live through most of the key stage.
- Peak 889 = 864 + return predicate + 24 inverse helpers at the special (first three charts) returns; the
  quotient-62 phase uses 22 helpers; the extractor model gives 890 with the tail controller, 888 without.
- The public frontier line (Syed v13 and its forks) is at 890: it traded +1 qubit (+0.225%) for -14.5k T.

### Failures

Our v12 file: 3 value errors in 8 x 102,400 shots (3.7e-6 per shot, about 0.37 failed shots per draw, about
+0.17% expected cost). All public 889/890 failures are deterministic in the input (manwithacat 753424,
GrantPrater). Classes: (a) error weight below the cutoff (weight <= 16 at cutoff 17: 1.6e-7 of shots; weight <= 17
at cutoff 18: 7.7e-7); (b) "half-Euclid" classes at weights 30-35 (a register word at zero when an extractor step
begins, or a convert-stage degree gap > 2), about 3e-6 per shot, which dominates.

DQI weight law (verifier, `w_k^2`): P(k <= 16) = 1.6e-7, P(17) = 6.1e-7, P(18) = 2.7e-6, P(19) = 1.1e-5,
P(20) = 4.0e-5; mean weight 31.05; P(35) = 2.77%.

## 3. State of the art

### Competition and DQI papers

- **Khattar, Shutty, Gidney, Zalcman, Yosri, Maslov, Babbush, Jordan, "Verifiable quantum advantage via optimized
  DQI circuits", arXiv:2510.10967** (local: `Scientific Papers/2510.10967v1.pdf`). The reference design: reversible
  extended Euclid (synchronized EEA, explicit vs implicit Bezout access, Proos-Zalka style registers), Chien search,
  Forney; GF(2^b) Karatsuba multipliers written as "parity-control Toffolis"; 5.72M Toffolis, 1,893 qubits for this
  instance. Its in-place space argument is the starting point of every public decoder.
  https://arxiv.org/abs/2510.10967
- **Jordan et al., "Optimization by Decoded Quantum Interferometry", Nature 2025 / arXiv:2408.08292** (local:
  `Scientific Papers/Nature 2025.pdf`). Defines the weight law `w_k` that the verifier samples (mean weight 31.05,
  P(weight <= 17) = 7.7e-7) and why the decoder only needs to be right on typical errors.
  https://arxiv.org/abs/2408.08292
- **Khattar, Rubin, Gidney, Broughton, Low, Babbush, "Leveraging verifiable structure in quantum circuits to find
  more efficient compilations with agentic AI", arXiv:2610.23832** (the competition paper,
  `google-verifiable-quantum-circuit-compilation/verifiable_circuit_compilation.pdf`): kickmix gate set, fuzz
  testing by Fiat-Shamir seeding, the dqi contract and R(eps) pricing.
- **Gu, Kemeny, Yelin, Lukin, "Discovering state-of-the-art quantum circuits with large language models",
  arXiv:2610.22367.** LLM-evolved circuits. For this exact OPI instance: Berlekamp-Massey locator solve at about
  140k Toffolis (vs 1.9M for the synchronized EEA), with (a) both BM branches merged into one shared update because
  the locator is only needed up to a nonzero scalar, (b) the old discrepancy recovered from the locator's constant
  coefficient instead of stored, (c) a failure budget: one history entry (19 qubits instead of 665), the second
  needed with probability 1e-4. Whole decoder (T, Q) = (4.4e5, 1,837). Also a reversible divsteps (Bernstein-Yang)
  modular inversion that recovers branch decisions from the state (trailing zeros), with a stored-history cutoff
  and a scheduled step count below the worst case (failure 1e-5). Both ideas (branch recovery from state, bounded
  history with a priced failure budget) transfer to the dqi key stage. https://arxiv.org/abs/2610.22367

### Finite-field multiplication (the 99-T and 42-T kernels)

- **Barbulescu, Detrey, Estibals, Zimmermann, "Finding optimal formulae for bilinear maps", ePrint 2012/110**: exact
  small bilinear ranks and search methods. https://eprint.iacr.org/2012/110
- **Cenk, Ozbudak, "On multiplication in finite fields", J. Complexity 26(2), 2010**: best known bilinear ranks of
  F_{q^n}, q in {2,3,4}, n <= 18 (so mu_2(12) and mu_4(12) are tabulated; n = 24 is not).
- **Randriambololona, "Bilinear complexity of algebras and the Chudnovsky-Chudnovsky interpolation method",
  arXiv:1107.0336** (asymmetric interpolation with different divisors), and **Ballet, Bonnecaze, Tukumuli,
  arXiv:1303.7082** (elliptic Chudnovsky algorithms). https://arxiv.org/abs/1107.0336 , https://arxiv.org/abs/1303.7082
- **Forum 756476** (Tamura, pjmathematician, Dyst-0): mu_4(5) = 11 and mu_4(6) = 14 exactly, with Lean proofs
  (https://github.com/Dyst-0/mu4-lean-proofs). So 3 x 14 = 42 is optimal **via GF(4)**; mu_2(12) itself is not
  settled below 42 by this.
- **Forum 759048** (Chih-Kai Wang): rank-173 decomposition of the paired map (H D, H_hi D^2) over GF(2^24) by
  elliptic-curve interpolation over F8 (14 rational points x 6 + one degree-2 place at 15 = 99 for H D), with
  aligned evaluation pairs costing 7 instead of 12. Shows that 99 for a single GF(2^24) product is the plain
  F8-curve value, and that products sharing operands can be cheaper as a block.
- **van Hoof, "Space-efficient quantum multiplication of polynomials for binary finite fields with sub-quadratic
  Toffoli gate count", arXiv:1910.02849**: in-place Karatsuba with no ancilla, O(n^1.58) Toffolis. Relevant wherever
  the decoder multiplies two quantum polynomials (not scalar x polynomial). https://arxiv.org/abs/1910.02849
- **Banegas, Bernstein, van Hoof, Lange, "Concrete quantum cryptanalysis of binary elliptic curves", TCHES 2021 /
  ePrint 2020/1296**: reversible GF(2^n) arithmetic at minimal width, Bernstein-Yang and Fermat/Itoh-Tsujii
  inversion, squaring as free CNOT maps. https://eprint.iacr.org/2020/1296
- **Itoh-Tsujii inversion** (1988) and tower fields GF((2^6)^2), GF((2^12)^2): x^-1 = conj(x) N(x)^-1 with the norm in
  the subfield; this is what Balram Srinivasan's cheaper GF(2^24) closing uses.

### Euclid, BM and RS decoder structure

- **Bernstein, Yang, "Fast constant-time gcd computation and modular inversion", ePrint 2019/266** (divsteps,
  jumpdivsteps / recursive half-gcd with transition matrices). https://eprint.iacr.org/2019/266
- **Sarwate, Shanbhag, "High-speed architectures for Reed-Solomon decoders", IEEE TVLSI 2001**: inversionless
  reformulated BM (RiBM), a key-equation solver on a fixed 3t+1 register array (hardware analogue of an in-place
  decoder).
- **Horiguchi (1989) / Koetter (1996) error-value formula**: e_j = c / (B(X_j) sigma'(X_j)) with B the auxiliary
  (previous) locator, avoiding the evaluator Omega. It yields 1/e_j as a polynomial value (no inversion), the
  opposite of what the trace phase needs (see DQI-B06).
- **Gao, Mateer, "Additive fast Fourier transforms over finite fields", IEEE TIT 2010** and **Lin, Chung, Han,
  arXiv:1404.3458**: the additive FFT the v12 scan already uses (all basis changes public, CNOT only).

### Uncomputation and ancilla techniques

- **Gidney, "Halving the cost of quantum addition", arXiv:1709.06648**: measurement-based uncomputation of AND
  (HMR + CZ), the basis of every "free uncompute" here. https://arxiv.org/abs/1709.06648
- **Khattar, Gidney, "Rise of conditionally clean ancillae for optimizing quantum circuits", arXiv:2407.17966**:
  borrowing in-use wires that are provably clean when it matters (manwithacat's HW-09 used it on DQI, -294 T).
  https://arxiv.org/abs/2407.17966
- **Luongo, Miti, Narasimhachar, Sireesh, "Measurement-based uncomputation of quantum circuits for modular
  arithmetic", arXiv:2407.20167**: MBU lemma; uncomputing a 1-bit output with a correction needed half the time.
  https://arxiv.org/abs/2407.20167
- ECC/Shor resource papers for reversible field arithmetic and pebbling: Roetteler et al. arXiv:1706.06752;
  Haner et al. arXiv:2001.09580; Litinski arXiv:2306.08585; Gidney arXiv:2505.15917; Chevignard, Fouque,
  Schrottenloher ePrint 2024/222 (residue representations that trade qubits for recomputation).

### Public competition work on this decoder line (Apache-2.0 unless noted)

| Who | What | Gain | Where |
|---|---|---|---|
| Syed Asad Ali | v12 generator (889 Q, source public); v13 compiled file `46446af1` (890 Q, ~368,700 T, source **not** public) | v13 is -3.56% vs v12 | dataset `syed3000/vqcc-own-circuit-portfolio` v12/v13 |
| Ertugrul Ozer, Areeb Ali | extractor cutoff 17 / 18; 241 retained products (`821bf8d1`) | -0.14% (48-draw) | notebooks `ozertuu/quantum-ai-27-083-vqcc-dqi-cutoff-17`, `areebbb/qai-fewer-toffolis-same-qubits-lb-28-06015` |
| Balram Srinivasan | GF(2^24) closing keeps the 12-bit norm, inverts it via GF(2^6): 340 instead of 360 T per closing, 41 sites (`66b4f18b`) | -820 T | notebook `balramsrinivasan/vqcc-28-10-gf24-closings-in-syed-v13-dqi` |
| boristown | control merging, cached CSWAPs, measured flag uncomputations (HW-27 CCX -> HMR), 37 then 144 product-cache chains (`d54d3387`, `bb9bf415`) | ~-0.2% | notebook `boristown/qcc-gap-3-residue-repair-and-c43-native-compiler` |
| GrantPrater | merge of both patch sets, closing carried to the 8 unmatched sites (-152 T), 49 measured uncomputations of a 15-AND product (-49 T), 4 window-proven retained products (`575b043f`) | frontier | dataset `grantprater/vqcc-dqi-61c8348d-merged-closings` (`vqcc_stream_diff.py`, `vqcc_measure_out.py`, `vqcc_pair_window.py`) |
| Asterios Panoras | 3 more retained products in condition blocks (`327b0751`, -3 T) | tiny | notebook `asteriospanoras/qcc-28-18084-could-be-top-1-public-try-it` |
| Atakan Aldemir | four-Toffoli DQI phase component (from five), shared products; weight-16/17 targeted failing inputs | -0.02% | notebooks `atakanaldemir/quantum-circuits-iso-gap-dqi-h5`, `.../qcc-shared-products-and-controlled-refits` |
| busyaprime, Ozer | exact product caching (symbolic pass + LP over idle wires), on the older 892 Q line | -0.65% there | notebook `busyaprime/23-906-local-gm-exact-product-caching-892q-dqi` |
| manwithacat | HW-06/07/09: MBU of flag ANDs, two-wire ladder, conditionally clean ancillae; Rust `kickmix-engine` (about 10x faster dqi verification, 0.4 GB) | | https://github.com/manwithacat/kickmix-engine |

GrantPrater's open items (dataset README): 14 more Toffolis per closing (686 in total) wait for a second idle helper
wire; a cheap guard for the lost-shot class at the forward extractor; a finder for retainable pairs inside
condition blocks.

## 4. Known nulls (do not repeat without a new angle)

| Null | Number that killed it | Source |
|---|---|---|
| GF(2^12) multiplier below 42 Toffolis | none found; mu_4(6) = 14 proved, so 42 is optimal via GF(4); only lower bounds for mu_2(12) | 753185, 756476, 757030 #4, 755972 |
| Rank-11 form of a 12-T block run 1,183 times per shot | none exists | manwithacat, 759220 |
| Third aligned pair in the rank-173 J0 block | 0 of 16,704 candidates | 759048 |
| Per-position Forney (inverse + 12 CCZ per scan position) instead of the convert stage | needs an inverse at <= 16 T; current 66 T; the AES 8-bit inverse needs 32 ANDs | Balram Srinivasan notebook |
| Cheaper inverse for the reverse extraction | reverse 223 T per pair vs 192 forward; every route found costs the same | alex chilton, 753185 |
| Narrower reverse division (one qubit less) | +1,518 T (+0.32%) > one qubit (0.225%) | 753185 |
| Peak shaving by reorder/recompute/borrow (877 Q build) | 2,046 peak windows in 684 kernel calls; cheapest narrowing 23-42 T per window | 753185 |
| Static ancilla remap (892 Q line) | all ancillas co-live from extract to unextract | rikirinjani / dcvxzyk, 753185 |
| Caching the reverse inversions | costs more than it saves for 12-240 cache qubits | dcvxzyk, 753185 |
| Dropping the jet zero flag / the BM gap flag | +2 value errors (+0.85%); 24 colliding post-step states | 753185 |
| Freeing the constant-1 border word in reverse folds | 30 value errors, R +7.1% | 753185 |
| Newton/Hensel phase-only sweep; Fermat/Itoh-Tsujii everywhere | does not converge; about 1.16M T one way | manwithacat, 753185 |
| Reusing the half-degree "Dialog" division in the reciprocal | divides modulo E where O^2 mod sigma is needed | manwithacat, 753185 |
| Porting BM-decoder gadgets to the CRT/Schur decoder | structures absent | 753185 |
| Measure-out pass on the current pick | no new unconditioned sites (the rest sit inside condition blocks) | 759220, GrantPrater |
| Sampled screening faster than the verifier; GF(2) tensor-core simulation | no gain | 753185 |
| R floor | R(0 failures) = 1.19126 is set by 102,400 shots; nothing goes below it | rikirinjani, 753185 |
| Clearing input-dependent scratch with `R` | caught: kickmix `R` randomizes the phase of a 1 | Syed Asad Ali, 755302 |
| Raising the cutoff above 18 | weight 18 has P = 2.7e-6 (0.28 failures per draw, +0.13%) against 138 T (0.036%) saved | our arithmetic from the weight law |
| Rewriting on a textbook decoder (BM of Gu et al., synchronized EEA of Khattar et al.) | 4.4e5 T at 1,837 Q and 5.1e6 T at 1,938 Q: an order of magnitude worse in Q^2 T than the CRT/Thiele line | arXiv:2610.22367, 2510.10967 |
| Phase-only Forney by residues (sum over roots of F/sigma' = top coefficient of F mod sigma) | F would include the interpolated mask polynomial of degree ~4,094: reducing it mod a quantum sigma costs millions of Toffolis | our analysis |
| Faster scan zero test | 11 Toffoli-likes per position is the AND count of a 12-input zero test with the phase folded in; the gauge can zero only 70 positions (840 GF(2) degrees of freedom, 12 constraints each) | our analysis; v12 ledger |
| "Re-rolling" the circuit hash for a lucky public draw | changes the public draw only; the private salt is independent, expected cost unchanged | 753424 |
| Z-basis measurement or `R` on data-dependent scratch | physically collapses the superposition (rule risk, and `R` is already caught by the verifier); excluded on purpose | 755302 |

## 5. Ideas

Conventions. Gain = expected task-cost reduction in % (positive is better), against our v12 base unless the row says
otherwise. Priority score = P(success) x midpoint of the gain range / effort in hours (units: % per hour). Effort
includes building, the official verifier on at least 8 salts (WSL, `--max-workers 2`, about 2 min per salt) and
attribution of any new failed shot. Every row's test protocol ends the same way: (1) the generator rebuilds
deterministically, (2) component checker from the v12 tree where one exists, (3) official verifier on default +
`private-1..7` salts, (4) compare mean cost and failures with the v12 8-salt baseline (3.6132e11, 3 failures).

Downloaded tools (GrantPrater's `vqcc_*.py`, boristown's patchers) are untrusted: we read them as documentation
and re-implement, we do not execute them (repo rule). Our own `.kmb` reader is therefore a shared prerequisite
(about 6 h, counted once, in DQI-A01).

### A. Improve our code (v12 generator, ports of public patches, failure budget)

| ID | Idea | Mechanism | Est. task-cost gain % | P(success) | Effort (h) | Priority score | Risks / interactions | Sources |
|---|---|---|---|---:|---:|---:|---|---|
| DQI-A01 | Recover Syed's v13 delta and port it into the v12 generator | Write our own `.kmb` reader and an instruction-stream aligner (same idea as GrantPrater's `stream_diff`: 64-bit op keys, resync on 48-op common runs), diff `7cbe8ecf` (v12) against `46446af1` (v13, Apache-2.0), attribute each differing region to a v12 source function (stage boundaries from `fixed_pipeline.emit`), then re-implement the changes in the generator. v13 is 890 Q and ~14.5k T lower | 3.0 to 3.6 | 0.55 | 40 | 0.0454 | v13 source is not public: some regions may be hard to read back into source (fallback: transplant compiled regions, which is not prize-eligible without a generator). +1 qubit is part of the trade. v13's failure rate is ~0.33/draw, similar to v12. Everything else in category A is measured on top of this | `downloads/datasets/syed3000__vqcc-own-circuit-portfolio` (v13 hash in SHA256SUMS); GrantPrater README; Balram notebook ("v12 has 50 calls, v13 41 exact") |
| DQI-A02 | Norm-preserving GF(2^24) closing in `gf24_nested_norm.emit_closing` | Keep the 12-bit norm N alive and invert via the GF(2^6) subfield: N^-1 = conj(N) C^-1, C = N_6(N), conjugation free in the tower basis. 340 instead of 360 static T per closing, same 24 helpers. v12 calls it 50 times (41 + 9 variants with the 67-instruction zero-test block) | 0.21 to 0.26 | 0.85 | 6 | 0.0333 | Must keep the rare-branch behaviour (norm 1, norm in GF(64), zero source); Balram stress-tested each class. Additive with A01 only if v13 did not already change these closings (v13 still uses the old gadget) | Balram Srinivasan notebook; GrantPrater README (8 extra sites) |
| DQI-A03 | Extractor cutoff j_min 15 -> 18 | Rolling tail controller costs 138 T per slot pair above j_min; 3 fewer slots = -414 T (-0.108%). Weights <= 17 are no longer decoded: +7.7e-7 per shot = +0.079 failures per draw (+0.036%) | 0.05 to 0.08 net | 0.95 | 2 | 0.0309 | Generator asserts `j_min in (15,16,17)` in `fixed_pipeline.emit_key_prefix`; 18 needs the same check widened (Areeb Ali did it). Raises failure variance on the single private draw (weight-17 inputs are known: Atakan Aldemir's targeted set). Cutoff 19 is a null (row in section 4) | Ozer cutoff-17 and Areeb Ali cutoff-18 notebooks; 759220 |
| DQI-A04 | boristown's native patches inside the generator | Control merging (merge consecutive gates sharing a classical condition), cached controlled swaps, HW-27 `CCX -> HMR + CZ if m` for flag uncomputations, written at the emitting functions instead of by binary patching | 0.10 to 0.25 | 0.70 | 12 | 0.0102 | Measured-flag uncomputations are "equal on samples only" in the public files (an identity of the field arithmetic): prove each with an exact symbolic check before trusting it. Check 8 salts for new failure classes | boristown notebook; GrantPrater README "What is proven" |
| DQI-A05 | Generator-level measurement-based uncomputation audit | Instrument the builder: every CCX whose role is erasing a value (flags, products, AND ladders) becomes `HMR` + Clifford fix-up; every repeated AND rebuild in the reverse pass that only serves erasure becomes a two-wire ladder executed only when the outcome is 1 (manwithacat HW-07: ~2.25 T per site on average). Today only 847 of 384,003 static T are measurement-conditioned | 0.05 to 0.30 | 0.45 | 10 | 0.0079 | Phase fix-ups must be exact on every classical branch (Khattar's 50%-failure signature, 752344). Shares sites with A04 (count once) | 752344 (host answer); manwithacat HW-06/07; Gidney arXiv:1709.06648 |
| DQI-A06 | Product retention / caching with exact lifetimes from the source | The generator knows which wires are idle and what each holds, so retained products (Areeb 241, busyaprime LP, GrantPrater window pairs, Asterios 3) can be found without a value tracker: a product a&b computed twice with an idle clean wire between is kept (`CCX(a,b,w) CX(w,t1) ... CX(w,t2) HMR(w) CZ(a,b) if m`); include pairs inside condition blocks (GrantPrater open item 3) | 0.10 to 0.40 | 0.60 | 20 | 0.0075 | A retained wire raises the live count where it is held: forbidden across the 889 peak windows (needs A09's profiler). Collisions with A04 sites | busyaprime notebook; Areeb Ali notebook; GrantPrater `vqcc_pair_window` README; Asterios notebook |
| DQI-A07 | Failure-budget audit of rare-branch handlers | List every component that exists only for rare inputs (gap-2 coverage in `raw_gap2_bezout`, zero-source admissions, first-three-chart special inverses 6 x 191 + 6 x 195 + 620 + 654 T, adaptive count cap, tail-controller corner cases). For each: branch probability p (classical simulation of the decoder on 1e7 sampled inputs, numpy, no kickmix) and its Toffoli cost X. Drop it iff X > 1.8e8 x p; conversely add a handler for any class with 1.8e8 x p > its cost (Gu et al.'s bounded-history pricing) | 0 to 0.6 | 0.35 | 14 | 0.0075 | Zero pivots occur in ~0.85% of shots (35 slots x 1/4096) and must stay handled. Every drop adds private-draw variance; prefer drops with p < 1e-7 | Gu et al. arXiv:2610.22367; 753424; 752344 |
| DQI-A08 | Width / Toffoli Pareto pass | Price each component at k and k+1 helpers (e.g. quotient phase 62 T with 22 helpers vs 63 T with 18; narrower reverse division +1,518 T per qubit saved, 753185). Choose Q globally: one extra qubit pays when it buys > 862 T anywhere; one qubit less pays when it costs < 862 T at every peak window | 0 to 1.0 | 0.30 | 16 | 0.0075 | Needs A09's liveness profile. Overlaps A01 (v13 is exactly such a trade) | v12 WIDTH_COMPATIBILITY docs; 753185 |
| DQI-A09 | Liveness profiler and peak shave of the 889 windows | Build an exact live-qubit profile of the emitted tape (allocation intervals from our `.kmb` reader); list every window at 889/888 and its owner (special returns with 24 inverse helpers + predicate, tail controller, quotient phase). If the 889 peak lives only in the 6 special returns, re-schedule them (inverse with 12-bit retained norm needs fewer simultaneous helpers, or a measured predicate) for 888 | 0 to 0.45 | 0.25 | 14 | 0.0040 | The 877 Q line had 2,046 peak windows (null in section 4); v12 may be similar. Profiler is reused by A06, A08 | 753185 peak-shaving null; v12 WIDTH docs |
| DQI-A10 | Closing 15-AND product: measured uncomputation and a second helper | GrantPrater: each closing computes one 15-AND product twice, 1,659 instructions apart; one product held on an idle wire saves 1 T per closing; 14 more per closing (686 total) need a second helper idle in the span. In the generator we can allocate it, if the span is below the peak | 0.013 to 0.18 | 0.50 | 8 | 0.0063 | Needs A02 first (the new closing). Helper must not raise Q: check with A09 | GrantPrater README "What is open" item 1 |
| DQI-A11 | Salt-robust selection and failure attribution | Install `kickmix-engine` (Apache-2.0, Rust, ~10x faster, 0.4 GB) as a second simulator; rank candidates on 48 salts; replay every failed input on our build and on the reference decoder to attribute it (GrantPrater section 5 method) | 0 to 0.05 (risk control) | 0.90 | 4 | 0.0056 | Must cross-check kickmix-engine against the official verifier before trusting it (its `--cross-check` flag). Speeds up every other row | 753642; 753424; GrantPrater notebook |
| DQI-A12 | Atakan Aldemir's four-Toffoli phase component | Replace the five-Toffoli phase component by the four-Toffoli one (shared products, CCZ) at every site of the generator | 0.01 to 0.03 | 0.80 | 3 | 0.0053 | Small; verify the derivation's 4,096 input/mask checks | `atakanaldemir/quantum-circuits-iso-gap-dqi-h5`; 758678 |
| DQI-A13 | Cheaper GF(2^12) inverse (66 T today) | Tower GF((2^6)^2) with GF(2^6) inverse via GF((2^3)^2), squarings as free CNOT maps, outputs uncomputed by MBU; or Itoh-Tsujii on the GF(64) norm only. ~100 calls (convert, K24 norms, quotient phases) | 0.10 to 0.40 | 0.30 | 16 | 0.0047 | Helpers per inverse must not exceed current; zero-input flags must be preserved | Banegas et al. ePrint 2020/1296; Balram notebook (66 T figure) |
| DQI-A14 | Phase-debt merging audit | Several HMR debts over a common operand can be repaid by one affine phase (as `phase_affine_quotient` already does for the extractor return). Search the generator for debts repaid separately over the same divisor or operand, in the convert and synthesis returns | 0.05 to 0.30 | 0.30 | 14 | 0.0037 | Product debts over GF(2^12) have full rank 12 for a nonzero mask, so savings only come from merging, not from rank | `raw_query_extractor.emit_closing`; `paired_converter_scale` |
| DQI-A15 | Exhaustive re-synthesis of small recurring kernels | SAT / exhaustive search (fixed gate count, all inputs and all measurement branches) on kernels that repeat thousands of times: 12-bit zero test + phase (11 T, ~4,024 calls), tail controller step (138 T), 15-AND product, special-chart prefix/return (620/654 T) | 0.10 to 0.50 | 0.20 | 20 | 0.0030 | Many kernels are already at their bound (manwithacat: no rank-11 form of a 12-T block); pick kernels whose count lacks a lower-bound argument | 759220; 757030 #1 (SAT encodings) |
| DQI-A16 | Better middle-product ranks in the Toeplitz solver | `unit_toeplitz.cross` uses evaluation at a+b-1 points of GF(2^12) (42 T each) or schoolbook. Check every (a, b, d) shape actually used: Karatsuba-type or point-at-infinity formulas, or GF(4)-structured interpolation, may beat 42(a+b-1) for the small shapes that fall back to schoolbook | 0 to 0.10 | 0.30 | 6 | 0.0025 | Small; must keep zero workspace | `unit_toeplitz.py`; Barbulescu et al. ePrint 2012/110 |
| DQI-A17 | Cheap guard for the dominant failure class | Dominant losses (~3e-6 per shot) have a register word at zero when a forward extractor step begins (half-Euclid classes, weights 30-35). Test it where a step already holds a subfield norm (GrantPrater item 2) and route it to the existing special-chart handler | 0.05 to 0.15 | 0.35 | 16 | 0.0022 | Pays only if the guard costs < ~500 T (0.37 expected failures x 1,760 T). Also reduces private-draw variance, which the score average hides | GrantPrater README items 2 and 4; Balram notebook "ideas wanted" 2 |
| DQI-A18 | Larger adaptive mask gauge | Classical logic is free: extend the selector from 64 candidates / cap 7 to thousands of candidate gauges and pick the one zeroing most scan positions (each zero saves 11 T). Today -21.65 T mean | 0.005 to 0.015 | 0.60 | 6 | 0.0010 | Grows `num_ops`, classical bits (cap 699,050) and file size; Poisson tail makes gains tiny | `integration_v7/README.md`; `adaptive_gauge64.py` |
| DQI-A19 | Operand-sharing block decompositions | Find pairs of products sharing an operand in v12 (extractor return: current x carry and current x next word; synthesis) and search joint decompositions below 2 x 99 (Chih-Kai Wang's aligned-pair method gives 7 instead of 12 per aligned F8 pair, but only for twisted pairs) | 0 to 0.50 | 0.10 | 20 | 0.0010 | Generic (a b, a c) tensors are probably rank-additive; only structured pairs gain | 759048 |

### B. Rewrite with a different technique

| ID | Idea | Mechanism | Est. task-cost gain % | P(success) | Effort (h) | Priority score | Risks / interactions | Sources |
|---|---|---|---|---:|---:|---:|---|---|
| DQI-B01 | Agentic evolutionary loop over the generator's kernels | Gu et al.'s recipe: an LLM proposes rewrites of one emitting function at a time (query step, closing, synthesis node, convert frame, scan cell); an evaluator rebuilds, runs the v12 component checkers, then kickmix-engine on reduced shots, then the official verifier at gates. Fitness = R x Q^2 x T with failures priced. Several families in parallel (PLAN backlog) | 0.5 to 3 | 0.30 | 40 | 0.0113 | Agent loops found little on chemistry gauges (753185), but Gu et al. report 12x on this decoder family vs EEA. Token cost; overfitting to one salt (rank on 8+ salts) | arXiv:2610.22367; 754557; PLAN section 11 |
| DQI-B02 | Make the reverse key as cheap as the forward key | Reverse key 115,841 T vs forward 101,087 (+14,754). The closing recomputes the pivot inverse, two products and an affine quotient phase per slot (223 vs 192 T per pair). Redesign the forward pass to leave exactly the classical masks the return needs (deferred-mask style, as the extractor already does for inverses): e.g. keep the normalizer and next-word products as measured debts repaid by one affine phase instead of recomputing them | 1.0 to 3.8 | 0.15 | 60 | 0.0050 | alex chilton: "every route we found costs the same" for a cheaper inverse alone. Rustam Bazarbayev open problem 3. Interacts with A01 (v13 may already touch the closings) | 753185; 757030 #3; `raw_query_extractor.emit_closing` |
| DQI-B03 | Fuse the projective finish into the synthesis output basis | Synthesis emits sigma and R in a coefficient basis; the projective finish (7,856 T) and basis changes then feed the additive FFT. Emitting synthesis directly in the FFT's triangular basis, or in a projective form the scan accepts, may remove most of the finish | 0.5 to 2.0 | 0.15 | 30 | 0.0050 | The finish also normalises the locator; the zero test is scale-invariant but R's scale is not (trace phase), so one scalar division remains | `projective_finish.py`; `additive_fft_phase_scan.py` |
| DQI-B04 | Quadratic-node CRT (K24) form for the convert stage | The extractor's gain is algorithmic: an O(n^2) chain over n/2 packed GF(2^24) nodes costs (n/2)^2 x 99 instead of n^2 x 42 (-41%). The convert's half-Euclid (O vs E, degree ~17) is still in GF(2^12) coefficient form (`raw_gap2_bezout`, scalar x polynomial updates). A Thiele/CRT formulation of that half-Euclid and of R = x Omega E^-2 mod sigma could cut its quadratic part by ~40% | 1 to 6 | 0.10 | 100 | 0.0035 | The convert already uses evaluation-based middle products (Toom, 42 T per point) in its Toeplitz solver, which are quasi-linear: the gain applies only to the remaining quadratic scalar updates. Keep or price the gap-2 coverage (A07) | `raw_gap2_bezout.py`, `unit_toeplitz.py`; Khattar et al. 2510.10967 section II.B |
| DQI-B05 | Bezout reuse from the extractor's convergents | The continuant identity A_k B_(k-1) - A_(k-1) B_k = c prod(q) gives the inverse of the previous convergent modulo sigma (Horiguchi-Koetter). If E^-1 mod sigma (or O^-2) can be obtained from quantities the key stage still holds before its closing, the half-Euclid shrinks or disappears | 1 to 5 | 0.07 | 60 | 0.0035 | Our analysis: the identity yields 1/e_j as a polynomial while the trace phase needs e_j, so one inversion remains; useful only if a half-size inverse replaces the full one. Needs the previous convergent alive (width) | Horiguchi 1989, Koetter 1996; 753185 (Dialog null) |
| DQI-B06 | Exact-synthesis search for one fused extractor slot | Instead of composing query steps one product at a time, search (SAT/ILP on small fields, then lift) a fused kernel for the last 2-3 query steps of a slot plus its normalizer, exploiting that the public node constants make the chain affine in each pivot | 0.5 to 3 | 0.08 | 50 | 0.0028 | Products of distinct quantum operands do not compress (rank additivity); gains can only come from cancelling public scales and swaps | `raw_query_extractor.query_step`; 759048 |
| DQI-B07 | Recover branch decisions from the state instead of storing flags | Gu et al.'s divsteps trick: a stored branch bit (tail-controller predicates, zero-pivot flags, the third chart's 37 pending masks, return predicates) is re-derived on the way back from registers that already encode it (a zero test of a word that the step leaves zero exactly on that branch). Frees qubits at the peak or removes predicate upkeep | 0.2 to 1.0 | 0.15 | 30 | 0.0025 | Must be exact on every branch or priced as failures (A07). The qubit part only helps at peak windows (A09). Gap and jet flags were shown to carry information (null) | arXiv:2610.22367; 753185 |
| DQI-B08 | GF(2^24) bilinear rank below 99 | 99 = 14 x 6 + 15 (elliptic curve over F8) = 11 x 9 (GF(16)^6) = 3 x 33 (GF(4)^12 if mu_4(12) <= 33). Search rank 96-98: genus-2 curves over F8 (up to 18 points) with asymmetric divisors, mu_4(12) from the Cenk-Ozbudak tables, flip-graph search seeded with the 99-term scheme. Each -1 rank saves ~1,190 T in the queries alone | 0.3 to 1.5 | 0.10 | 40 | 0.0023 | No optimality proof exists either way; mu_8(8) >= 16 and the point counts make 96 hard. A 576-pair bilinear check tests a scheme in 1 s | 759048; 756476; Cenk-Ozbudak 2010; arXiv:1107.0336; ePrint 2012/110 |
| DQI-B09 | Explicit incremental synthesis so that queries become free evaluations | With the convergent polynomials held explicitly during extraction, a slot's query is a public evaluation (CNOT only) plus ~2 products instead of i products: the 117,810-T query bill could fall by up to ~100k. Obstacle: the forward recurrence needs two consecutive convergents (numerator and denominator), ~48k bits at slot k while the consumed heads free only ~24k bits. Look for a hybrid (explicit only for late slots) or a representation storing one convergent and recovering the other cheaply | 0 to 10 | 0.03 | 80 | 0.0019 | Probably blocked by width, which is why the line uses implicit Thiele access; listed because its payoff is the largest of any row | our analysis; Khattar et al. (explicit vs implicit Bezout access); `fused_crt_synthesis.py` |
| DQI-B10 | Two slots per chain with a shared pivot operand | Run the query chains of slots i and i+1 together so each earlier pivot P_l multiplies two currents; search a joint decomposition of (P b, P c) below 2 x 99 | 0 to 3 | 0.05 | 40 | 0.0019 | Generic (ab, ac) is very probably rank-additive; slot i+1 also needs P_i, which exists only after slot i ends | 759048 |
| DQI-B11 | 12-qubit carry bank | Keep the running reciprocal in a form that needs 12 qubits (pivots normalised into the GF(2^12) subfield, or the carry folded into the current head): -12 qubits = -2.7% if Toffolis do not grow | -1 to 2.7 | 0.05 | 60 | 0.0008 | Each normalisation into the subfield costs a norm (~42 T) per slot; the 24 inverse helpers still set the peak unless reworked too (A09) | `raw_query_extractor.py`; WIDTH docs |
| DQI-B12 | Cubic-node CRT (GF(2^36) heads) | 70 = 23 x 3 + 1 nodes: queries 552 x ~153 T (mu_16(9) x mu_2(4)) = 84.5k vs 117.8k (-33k), but a 36-bit carry and 36-bit inverse helpers add ~24 qubits (+5.4%, ~20.7k T equivalent) and every inverse grows | -2 to 4 | 0.08 | 150 | 0.0005 | Whole-decoder rewrite (extractor, synthesis, convert interfaces), mixed node degrees | our arithmetic from the extractor ledger |

## 6. Combination notes

- **Order of work.** A11 (fast simulator, failure attribution) and the `.kmb` reader of A01 come first: they serve
  almost every row. Then A01 (v13 delta, the only multi-percent item with good odds), then the cheap ports in one
  build on top of it: A02 + A10 (closings), A03 (cutoff 18), A04 + A05 (MBU), A12 (phase component), and A06
  (product retention) last, because it consumes idle wires that other edits may want. Expected stack: about
  -4.0 to -4.5% task cost, at or slightly below the public frontier (3.4641e11), but from source.
- **Do not double count.** The frontier's -3.97% against us is v13 (-3.56%) + Balram closings (-0.22%) +
  boristown, GrantPrater and retention patches (about -0.2%). A01 + A02 + A04 + A06 + A10 reproduce that, not more.
  Beating the frontier needs A07, A08, A09 (failure budget and width) and category B.
- **Width edits conflict with retention.** Retained products (A06, A10) hold extra wires; decide any qubit removal
  (A09, B07, B11) first, then let retention fill only the windows below the new peak.
- **The failure budget is shared.** Cutoff 18 (A03) adds ~0.08 failures per draw; each handler dropped in A07 adds
  its own. Keep a running total of expected failures per draw (target <= 0.5 including the line's ~0.33) and
  attribute every new failed shot by replay before accepting a build. A17 (guard) can buy budget back for A07.
- **Re-measure the anatomy after A01.** Section 2 describes v12; v13's -14.5k T changes the stage shares, so B02
  (reverse key) and B04 (convert) must be re-priced on the v13-equivalent build.
- **Where a large gain can come from.** Matching the 31.4 leaders needs about -28% on one task. Nothing in category
  A comes close; on dqi only B02, B04, B09 (or an idea not listed here) could. B01 is the engine to explore B02-B06
  cheaply.
- **Rules.** Every public patch is re-implemented in our generator (prize eligibility needs the generating source)
  and credited: Syed Asad Ali first (his NOTICE, Apache-2.0 section 4(d)), then each patch author.

