# invsqrt: research notes and prioritized ideas

Written 2026-11-17. Sources: the local forum snapshot (`Instructions/discussion/`, 2026-11-15), the live list of
public Kaggle notebooks (pulled 2026-11-17 into `downloads/kernels/`, read as data only), the official verifier, the
v12 generator in `third_party/syed3000-vqcc-own-circuit-portfolio-v12/invsqrt/invsqrt_source/`, the competition paper
and web sources. Our own measurements are labelled *ours*; numbers quoted from other teams are labelled with their
author.

## 1. Summary

1. **Where the cost is (our v12, 46 Q, C = 1,742.67):** the divider and its clean-up take 34% of the charge
   (598 Toffolis), coefficient/seed/product arithmetic 22% (389), rotations 13% (231 for 33 rotations, ~7 each),
   the square root 11% (192), the PHASE ring 6% (105), the table loads 8% (133), normalisation 5% (94). Half of all
   Toffolis are uncomputation. The 46-wide peak is spread over 12 stages (8,137 ops), not one plateau.
2. **What the public frontier did (74,560.91, −7.0% vs us):** Özdağ's 64-row table (−1.74%), Syed's own v13 move to
   44 Q (root bits v10, v11 and flags e2, e3 kept off the wires; −4.6%), then micro-edits and rescheduling (−0.8%).
   None of it has public source after v12.
3. **Most promising:** rebuild that chain in our generator (A02 → A03 → A04, about 7%, ~55 h), with the compiled
   frontier file as a leaderboard stop-gap (A01). After that, the circuit is locally tight (0 affine Toffolis,
   0 dead writes, parking fails): expect 1–3% from many small, proven passes (A05–A20), and only a new arithmetic
   family (B01–B04) or a 43-qubit cover (A21–A25) could give a large step. The leaders' private 28%-class gain, if it
   is on invsqrt, must be such a family.

## 2. Cost anatomy of our circuit

### 2.1 How the cost is computed (from the verifier)

`cost = Q × (Toffolis + T/4)`, where Toffolis is the mean number of executed `CCX`/`CCZ`/`CSWAP` per shot (a gate
inside a `PUSH_CONDITION` on a measured bit counts only on the shots where it fires), and every executed `Z_POW`
whose angle is not a multiple of 1/8 turn costs `0.56·log2(1/eps) + 5.3` T, with
`eps = π · min_margin / 2^64 / n_rot`. `min_margin` is the smallest `tol(x) − |error(x)|` over **all** 3,133,452
inputs, in absolute units of 2^-64 turns, and `n_rot` is the mean number of executed rotations per shot.
Consequences used below:

- One rotation costs about 7.0 Toffoli-equivalents at our margin (v12: 230.85 / 33 = 6.995, i.e. 28 T each).
- Doubling `min_margin` saves `0.56/4 × n_rot` = 4.6 Toffoli-equivalents at 33 rotations; perfect rotations (zero
  arithmetic error) would save at most about 17 (jeffemerson's estimate).
- Syed's Q48 base generator (candidate10 README) had a margin of 124,370,688 units, and the v12 rotation charge
  (230.85) matches that margin exactly; Syed's v13 44-qubit file 22,067,537,
  re-solved by Asterios Panoras to 22,102,971; the frontier `563d0583` runs at worst error/tolerance 0.99956. A
  smaller table pushed the margin down by about 5×, costing about +11 Toffoli-equivalents of rotation charge, which
  the table saving more than paid for.
- One qubit is worth `C/(Q−1)` Toffoli-equivalents: 38.7 at 46 Q (our circuit), 39.4 at 44 Q (frontier).
- The verifier checks every integer x in [1, 3,133,452], so x values that cannot occur in the physics (non-sums of
  three squares) cannot be exploited (topic 753451).

### 2.2 Our circuit `321aff83` (v12, 46 qubits) stage by stage

How obtained (*ours*): converted `artifacts/v12/invsqrt.kmb` to `.kmx` with kickmix, then ran a small bit-sliced
basis-state simulator (no phase tracking) over 20,000 inputs (x = 1..64 plus random x) with random HMR outcomes.
It reproduces the exact mean Toffoli count (1,511.55 sampled vs 1,511.57 exact), restores x and cleans all scratch.
Stage boundaries are the generator's own cut points as published by manwithacat (*Fortysix*, section 4.9; the four
divider boundaries are approximate to about ±30 ops). "Width" is the number of qubits that are non-zero on at least
one sampled input (a lower bound on live width, equal to it at the peaks). Script:
`<scratchpad>/anat.py` + `stages.py` (not kept in the repo).

| Stage | Ops | What it does | Mean Toffolis | Static Toffoli-like | HMR | Width (min–max) | Ops at 46 |
|---|---|---|---:|---:|---:|---|---:|
| NORM | 0–609 | shift x left by 2k, write k (0..10) in 4 flag wires | 51.00 | 51 | 19 | 22–36 | 0 |
| ROOT | 609–1,638 | consuming square root: 22 wires → root v (12) + residual h (10) | 96.00 | 96 | 84 | 26–39 | 0 |
| LOAD1 | 1,638–5,248 | first 128-row table load (address = top 7 bits of v): D13 + H3 | 56.00 | 56 | 37 | 27–46 | 2,468 |
| COEF1 | 5,248–5,563 | coefficient refined by low root bits | 44.00 | 44 | 36 | 42–46 | 97 |
| META1 | 5,563–5,592 | metadata measured away | 1.00 | 1 | 6 | 38–43 | 0 |
| SEED | 5,592–5,704 | three-row seed | 10.00 | 10 | 7 | 38–46 | 8 |
| PROD | 5,704–6,707 | product D × h → signed 15-bit P | 103.00 | 103 | 94 | 43–46 | 361 |
| CRET | 6,707–6,720 | 13 coefficient wires HMR'd | 0.00 | 0 | 13 | 31–43 | 0 |
| DIV_FWD | 6,720–9,877 | 24-row non-restoring divider by s = (4096+v)/4096, keeps quotient history | 270.00 | 270 | 239 | 32–46 | 704 |
| PHASE | 9,877–10,668 | 33 rotations + 1 T; 63 CSWAPs rotate the quotient by k | 105.00 | 105 | 40 | 40–46 | 88 |
| DIV_EXIT | 10,668–10,808 | divider exit | 5.00 | 5 | 4 | 40–45 | 0 |
| DIV_CLEAN | 10,808–17,306 | divider rows backwards; history signs retired and rebuilt | 323.00 | 364 | 334 | 31–46 | 485 |
| LOAD2 | 17,306–21,361 | second table load (for the inverse product) | 71.00 | 71 | 48 | 31–46 | 3,264 |
| H3 / COEF2 / META2 | 21,361–21,766 | second coefficient | 44.00 | 44 | 38 | 43–46 | 97 |
| SEED_CLR | 21,766–24,203 | inverse product, seed cleared | 129.99 | 147 | 139 | 37–46 | 194 |
| CSRC / A_REST | 24,203–24,734 | coefficient inverse | 39.00 | 39 | 36 | 39–46 | 48 |
| R_PAID | 24,734–26,132 | raw root bits restored | 18.00 | 18 | 3 | 40–46 | 323 |
| BODY_RET | 26,132–33,310 | table clean-up and phase repayment (thousands of CZ/Z) | 6.31 | 30 | 26 | 26–45 | 0 |
| ROOT_INV | 33,310–34,523 | square root undone | 96.00 | 96 | 84 | 26–39 | 0 |
| NORM_INV | 34,523–34,778 | normalisation undone | 43.25 | 49 | 10 | 22–30 | 0 |
| **Total** | 34,789 ops | | **1,511.55** | 1,599 | 1,297 | 22–46 | 8,137 |

Grouped (Toffolis, share of the 1,742.67 charge):

| Group | Toffoli-eq | Share | Peak width |
|---|---:|---:|---:|
| Normalise fwd + inv | 94.25 | 5.4% | 36 |
| Square root fwd + inv | 192.00 | 11.0% | 39 |
| Table loads (LOAD1 + LOAD2) + BODY_RET | 133.31 | 7.6% | 46 |
| Coefficient / seed / product arithmetic, fwd + inv | 389.00 | 22.3% | 46 |
| Divider fwd + exit + clean-up | 598.00 | 34.3% | 46 |
| Phase endpoint CSWAP ring | 105.00 | 6.0% | 46 |
| Rotations (33) + 1 T | 231.10 | 13.3% | — |
| **Charge C** | **1,742.67** | 100% | **Q = 46** |

Observations (*ours* unless credited):

1. **The divider is a third of the charge.** Its clean-up costs 53 more than its forward pass (323 vs 270), the
   price of retiring history signs to fit in 46 wires (jeffemerson: inverse-adder segments about +22, sign rebuilds
   about +33).
2. **The inverse product costs 27 more than the forward product** (SEED_CLR 130 vs PROD 103), for the same reason.
3. **The circuit is an almost exact mirror.** Root 96 + 96, normalisation 51 + 43, loads 56 + 71, coefficient 44 + 43.
   Everything in the second half is uncomputation (775.6 of the 1,511.6 Toffolis, after the last rotation at op
   10,594).
4. **The peak is spread out.** 8,137 ops (23%) sit at width 46 across LOAD1, COEF1, PROD, DIV_FWD, PHASE, DIV_CLEAN, LOAD2,
   COEF2, SEED_CLR, A_REST, R_PAID. jeffemerson's X-ray: 450 short windows at width 46, each the top of a carry
   staircase (one carry on top of a 45-wide body); 89 windows (74 in the multiplier) have no parkable wire.
5. **Only 88 Toffolis are measurement-conditioned** (static 1,599 vs mean 1,511.6); they sit in DIV_CLEAN,
   SEED_CLR, BODY_RET and NORM_INV.
6. **Gate mix:** 1,461 CCX, 129 CSWAP, 9 CCZ, 38 Z_POW (33 generic rotations, 1 T, 4 Clifford), 1,297 HMR, 895
   conditioned blocks.

### 2.3 The public frontier `563d0583` (maoxin123, 44 qubits) compared

How obtained (*ours*): decoded the base64/zlib blob embedded in maoxin123's notebook
*QCC 28.17 | 4,022-Toffoli chemistry loaders* (SHA-256 `563d0583f736afc2d8111f376bc1bd5ce6aef686525b0dd281fc1f2b4eb09321`,
Apache-2.0 for the invsqrt per the notebook footer; Syed Asad Ali's NOTICE applies), same simulator.

| | ours `321aff83` | frontier `563d0583` | Syed v13/v14 `47b31838` |
|---|---:|---:|---:|
| Qubits | 46 | 44 | 44 |
| Mean Toffolis | 1,511.57 | 1,454.0 (sim 1,454.09) | 1,458.25 |
| Charge C | 1,742.67 | 1,694.6 | ~1,708.4 |
| Cost | 80,162.67 | 74,560.91 | 75,170.62 |
| Ops | 34,789 | 20,659 | — |
| Static Toffoli-like | 1,599 (1,461 CCX + 129 CSWAP + 9 CCZ) | 1,548 (1,518 CCX + 30 CCZ, **0 CSWAP**) | — |
| HMR | 1,297 | 1,277 | — |
| Ops at full width (sampled) | 8,137 | 4,299 | 435 windows (jeffemerson) |
| Source available | yes (v12) | no (compiled only) | no (compiled only) |

What the frontier contains (from its credits and from jeffemerson's X-ray of `47b31838`):
- Syed's own v13/v14 move to 44: root bits v10 **and** v11 are absent across the whole multiplier and rebuilt
  afterwards (v10 from the top 4 bits of D plus v9), and exponent flags e2, e3 are kept off the wires while h is
  explicit.
- Ali Haydar Özdağ's 64-row table (6 address bits, 13 intercept + 7 slope bits): both loaders re-synthesised as
  27 + 29 ANDs instead of 56 + 70; slope bit S4 measured out before the inverse product and rebuilt with 2 ANDs;
  angle refit. Alone on the 46-qubit build it gave −47 static Toffolis (−1.74%, `4560ac0c`).
- Asterios Panoras's dead-write fix (−1), jeffemerson's AND-merge cut (−2), Ertuğrul Özer's / Panoras's angle LP refit,
  manwithacat's measure-out edits, and maoxin123's rescheduling (34,789 → 20,659 ops, all CSWAPs written as CCX
  networks).

So the gap from us to the public frontier is **−7.0% task cost** (−4.3% from the two qubits, −2.8% from the charge),
and none of it is in our source tree yet.

Two quick scans (*ours*): (a) a Toffoli whose target is next touched only by an `HMR` (Asterios Panoras's
"dead-write" pattern) occurs **0 times** in v12 and 0 times in `563d0583`; 15 more sites in v12's DIV_CLEAN have a
`CZ`/`Z` read of the target in between, so they are not free. (b) Neither circuit has a classically conditioned
`Z_POW`, so all 33 rotations are paid on every shot.

## 3. State of the art

### 3.1 Papers and techniques

| Reference | What is relevant to invsqrt |
|---|---|
| Competition paper, `google-verifiable-quantum-circuit-compilation/verifiable_circuit_compilation.pdf` (arXiv:2610.23832), §A | Reference = normalise, cubic seed from QROM, one Newton step, phase kickback (260 Q, cost 2.08M). The authors explicitly invite "a simpler polynomial over a larger range" and "a good polynomial approximation over the input domain". |
| Rubin et al., *Quantum computation of stopping power for inertial fusion target design*, PNAS 2024 (arXiv:2308.12352), in `Scientific Papers/` | Origin of the task: x = a²+b²+c², 24-bit relative precision. |
| Su, Berry, Wiebe, Rubin, Babbush, *Fault-tolerant quantum simulations of chemistry in first quantization*, PRX Quantum 2021, https://arxiv.org/abs/2105.12767 | The QROM-seed + Newton construction the reference follows. |
| Häner, Roetteler, Svore, *Optimizing quantum circuits for arithmetic*, https://arxiv.org/abs/1805.12445 | Reversible fixed-point inverse square root via Newton / piecewise polynomials; pebbling-based qubit/gate trade-offs; dirty-ancilla tricks. Their designs are far wider than 44 Q. |
| Häner, Soeken, Roetteler, Svore, *Quantum circuits for floating-point arithmetic*, https://arxiv.org/abs/1807.02023 | Floating point saves qubits when range or relative precision is large; our normalise-by-4^k is exactly a float with an exponent of about 3.5 bits. |
| Sanders et al., *Compilation of fault-tolerant quantum heuristics for combinatorial optimization*, https://arxiv.org/abs/2007.07391 | QROM-indexed piecewise-linear function evaluation with uneven segments; phase by addition into a phase-gradient register. |
| Sun et al., *Low depth phase oracle using a parallel piecewise circuit*, https://arxiv.org/abs/2409.04587 | QROM-loaded interpolation then phase-gradient kickback; cost tables (note the published erratum). |
| Gidney, *Halving the cost of quantum addition*, https://arxiv.org/abs/1709.06648 | Temporary logical-AND with free measurement-based uncompute: the basis of every adder here. |
| Gidney, *Windowed quantum arithmetic*, https://arxiv.org/abs/1905.07682 | Table lookups replacing chunks of a multiplication; relevant to the coefficient product. |
| Gidney, *Approximate encoded permutations and piecewise quantum adders*, https://arxiv.org/abs/1905.08488 | Oblivious carry runways and tolerated approximation in arithmetic, matching our tolerance-driven design. |
| Gidney, *Spooky pebble games and irreversible uncomputation*, https://algassert.com/post/1905 (arXiv:1904.07356) | "Measure then rebuild", the core of the 46-qubit design, is a spooky pebble game; the paper gives optimal strategies. |
| Khattar & Gidney, *Rise of conditionally clean ancillae for optimizing quantum circuits*, https://arxiv.org/abs/2407.17966 | Comparators, incrementers and adders using dirty / conditionally clean ancillae: the tool to remove a carry wire at a peak without a clean qubit. |
| Takahashi, Tani, Kunihiro, https://arxiv.org/abs/0910.2530 ; Cuccaro et al., https://arxiv.org/abs/quant-ph/0410184 | Ancilla-free / one-ancilla ripple adders: trade Toffolis for the carry wire. |
| Muñoz-Coreas & Thapliyal, non-restoring square root https://arxiv.org/abs/1712.08254 ; integer division https://arxiv.org/abs/1809.09732 | Garbage-free non-restoring root (2n+1 qubits) and dividers; baselines for our ROOT (96) and divider (24 rows). |
| Ercegovac & Lang, *Digital Arithmetic* (Morgan Kaufmann 2004), ch. 5–7 | Digit-recurrence division, square root and reciprocal square root (radix 2/4, SRT, redundant remainders, on-the-fly conversion); the classical menu for divider redesigns. |
| Gu, Kemeny, Yelin, Lukin, *Discovering state-of-the-art quantum circuits with large language models*, https://arxiv.org/abs/2610.22367 | LLM evolutionary loop with a programmatic evaluator; reported 4×–9,100× space-time gains on primitives. The method to run on our generator's components. |
| DeepMind, *Discovering formally verified quantum programs with AI*, https://arxiv.org/html/2611.09822v1 (forum 757471) | Agent-driven discovery with formal verification; same loop idea. |
| Qualtran arithmetic bloqs, https://github.com/quantumlib/Qualtran (`qualtran/bloqs/arithmetic`, `qualtran/bloqs/rotations/phase_gradient.py`) | Reference adders, comparators, QROM/QROAM and phase-gradient rotations; the reference circuit was built with it. |

### 3.2 Public notebooks and forum material (live list pulled 2026-11-17)

| Source | Content |
|---|---|
| Syed Asad Ali, *VQCC Own-Circuit Portfolio* (notebook + dataset v12–v14), https://www.kaggle.com/code/syed3000/vqcc-own-circuit-portfolio-notebook | v12 = our base (source, Apache-2.0; NOTICE requires naming him first). v13/v14 ship a 44-qubit `47b31838` (1,458.25 Toffolis, 75,170.62) **compiled only**: the dataset file list (checked 2026-11-17) has no new source. |
| maoxin123, *QCC v27 / a 44-qubit invsqrt on the 28.09 frontier*, https://www.kaggle.com/code/maoxin123/qcc-v27-a-44-qubit-invsqrt-on-the-28-09-frontier and *QCC 28.17 / 4,022-Toffoli chemistry loaders* | `0af21cb8` (74,797.37) then `563d0583` (74,560.91, 1,454.0 Toffolis): Syed's line + Özdağ's 64-row table + Özer's angle refit + rescheduling. Compiled only, embedded as base64. |
| Ali Haydar Özdağ, *QCC / submit - verify and package* (v25), https://www.kaggle.com/code/ahoirg/qcc-submit-verify-and-package | 64-row table, loaders 27 + 29 ANDs, S4 rebuilt with 2 ANDs, robust hold pass; built in a private kernel (`ahoirg/qcc-invk6-k6c`) from an open fork of Syed's generator. |
| manwithacat, *Fortysix: where the 16 qubits went*, https://www.kaggle.com/code/manwithacat/fortysix-where-the-16-qubits-went | Stage map of v12 (used in section 2), the 48 → 46 patch ladder, leave-one-out replays (without patch 8: +69 cost; without patch 13: +92). |
| jeffemerson, *invsqrt X-ray: 46-qubit dead ends and ideas*, https://www.kaggle.com/code/jeffemerson/invsqrt-x-ray-46-qubit-dead-ends-and-ideas and *28.059 public stack + 44-qubit invsqrt X-ray* | Whole-domain simulator; parking census (46 → 45 needs ≥ 58 Toffolis and leaves 89 windows uncovered; 44 → 43 needs ≥ 49 and leaves 88); 14 priced nulls; 5 open directions. |
| jeffemerson, *28.071: 44-qubit stack + our AND-merge cut*, https://www.kaggle.com/code/jeffemerson/28-071-44-qubit-stack-our-and-merge-cut | Two AND merges on wire q28: −2 Toffolis. |
| Asterios Panoras, *QCC 27.30 / free CZ fix in invsqrt* and *QCC 28.12769 / what works, what fails and why* | Dead-write fix (−1 Toffoli); LP re-solve of the 34 rotation angles (the phase is exactly linear in the angles); 45 Q break-even +38, a 47 Q third carry slot loses. |
| Atakan Aldemir, *Early measurement and domain-aware inverse* | Deletes a CSWAP that is redundant on the legal domain (62-qubit line). |
| maoxin123, CSWAP-network re-synthesis on the 62-qubit line (topic 757363) | 226 CSWAPs → 205 at +9 CCX. |
| hosen42, *VQCC Our Submission* | Measure-first for write-only runs; 11 writes replaced by one CCX before the read. |
| Kirill Tymchishin, *Conditional uncompute audit*, *Pi-gauge & weight refinement* | Measured-target Toffoli identity (CCX followed by X-measurements of target and a control → measure first + conditioned CZ). |
| manwithacat, viewer https://manwithacat.github.io/kickmix-macos/viewer/ and `render.py` in https://github.com/manwithacat/kickmix-macos | Per-op traces, width and cost per stage for any `.kmb`. |
| Forum 757030 (Rustam Bazarbayev, open problems) | Minimum-CSWAP barrel shifter (w ≈ 21: public high 40s, counting bound high 30s; SAT gives about 2n−5 for small n). |
| Forum 759333 | A team reports a "directly evaluated piecewise quadratic" invsqrt (154 Q, not competitive). |

## 4. Known nulls (do not repeat without a new angle)

| # | Idea | Measured result | Source |
|---|---|---|---|
| N1 | Exploit x values that cannot occur (non-sums of three squares) | The verifier checks every integer x; the only effect would be ~0.08% rotation slack if the verifier changed. | 753451 (manwithacat) |
| N2 | Drop or retune rotations | Dropping 1/2/3 rotations → worst error/tol 1.354/2.214/3.937 (63 Q line); at 44–46 Q dropping one or adding a 34th is infeasible or loses; angles already at the LP optimum (LP re-solve: +0.16% margin, −0.0007% cost). | 753185, 756824, Panoras 28.12769 |
| N3 | Trade rotation precision regionally | Precision-saturated everywhere; no region with slack. | 753185 (manwithacat) |
| N4 | Perfect rotations (zero arithmetic error) | Saves at most ~17 Toffoli-eq. | jeffemerson #7 |
| N5 | Smaller table (32 or 16 rows) | 32 rows: 1.137× tolerance (1.049× with an extra table bit); 16 fail. 64 is the optimum for this arithmetic. | jeffemerson #2 |
| N6 | Narrower table words (12-bit A, slope −1 bit) | 1.338× tolerance. | jeffemerson #3 |
| N7 | Drop one of 32 Horner table bits (63 Q line) | Best variant 1.228× tolerance. | 753185 |
| N8 | Let the divider use 47–48 wires | Cap 48 saves 30 Toffolis but needs 75; 47 needs 37. A 47 Q third carry slot: 43.5 saved vs 37.4 needed. | jeffemerson #4, Panoras |
| N9 | Other dividers: radix-4, reciprocal by multiply, 2Q/(Q²+N), skipping the clean-up | All above the 706-Toffoli incumbent (init 13 + rows 257 + endpoint 115 + clean-up 320). | jeffemerson #5 |
| N10 | MBUC for the square root or the divider undo | The root is an in-place bijection with no garbage; the undo borrows carry information. | jeffemerson #6 |
| N11 | Hide the 4 exponent flags in the zero low bits of n | Only possible while h is explicit; the divider consumes h. (Syed's v13 did keep e2, e3 off the wires while h is explicit.) | jeffemerson #8, 28.059 X-ray |
| N12 | Bit-serial quotient rotation, split endpoint, absolute-frame division, phase-gradient rotations, q×d reverse | 8 redesigns, each priced above the incumbent. | jeffemerson #9 |
| N13 | Reorder the lazy lookup clean-up | 6.28 vs best alternative 7.00; ~0.5 Toffoli slack. | jeffemerson #10 |
| N14 | Cap the build at 45 and let the divider pay | Multiplier registers already total 46 (v 12 + flags 4 + C 13 + acc/quotient 16 + carry 1). | jeffemerson #11 |
| N15 | Re-load the coefficient C for every digit | +1,100–1,400 Toffolis to free 13 wires worth 484. | jeffemerson #12 |
| N16 | Drop the exponent flags and recompute them | Information-bearing (573,648 conflicting states for q22); ~3.5 bits must live somewhere. | jeffemerson #13, Fortysix §5 |
| N17 | Table-free methods (digit recurrence, Newton) | Hold three ~26-bit numbers: ≥ 74 Q before any gate. | jeffemerson #14 |
| N18 | Park one quadratic-rebuildable wire per peak (46 → 45, 44 → 43) | ≥ 58 Toffolis + 89 uncovered windows (46); ≥ 49 + 88 (44). One verified relation (q9 = q20 ⊕ q30 ⊕ q20·q32) clears the forward division only. | jeffemerson |
| N19 | Delete Toffolis whose product is already affine in live wires | 0 of 1,599 on v12; 0 of 1,550 on HW-26; 0 of 1,501 on the 62 Q file. | jeffemerson, 754152 |
| N20 | Cache a recomputed product | No product is recomputed more than 5 times (about 38 needed to pay a qubit). | jeffemerson |
| N21 | Whole-circuit T optimisers (PyZX `full_reduce`, TODD) | T-count −0–22% but 0 Toffolis removable. | jeffemerson #1 |
| N22 | Dirty-wire AND for a row-opening gate (63 Q line) | Borrowed wires read up to 8,474 times before the replay closes: invalid. | 753185 |
| N23 | Bipartite/multipartite tables, higher-degree polynomials with fewer segments, more MBUC on multi-bit sums | Partial products grow faster than the table shrinks; recompute eats the saving (62 Q line). | 757030 |
| N24 | Multi-segment splices between public files | Only rediscovers the existing pick. | 758678, 759220 |
| N25 | Dead-write measure-first sites left in v12 or the frontier | 0 strict sites in either. | *ours*, section 2.3 |
| N26 | GPU/tensor-core simulation of linear blocks; threaded verifier | No invsqrt linear blocks long enough; threads ~1.2× only. | 753185 |

Caveat: most nulls are scoped to a design family (jeffemerson states the family for each). N8, N9, N12, N14, N15
and N17 are the ones that close the obvious width moves; any B-category idea below must explain why it escapes them.

Two more nulls reasoned here (*ours*, not built):

| # | Idea | Why it fails | Source |
|---|---|---|---|
| N27 | Make rotations cheaper by putting them under a classical condition (the verifier counts executed `Z_POW`, so a conditioned one counts on about half the shots) | The phase θ·q depends on a data qubit; it cannot be moved behind a measurement outcome without a controlled rotation, which costs more. | verifier semantics |
| N28 | Replace Bennett uncomputation of the whole body by "HMR everything after the phase, repay the debt" | The debt is a random parity of the intermediate values; repaying it needs the values again (a forward pass plus its own clean-up), so it costs at least as much. The in-place stages (normalise, root, divider remainder) cannot be measured at all. | spooky pebbling argument |


Nulls on paper from the family studies of 19 Nov (S2.18, S2.21):

| # | Idea | Measured result | Source |
|---|---|---|---|
| F2 | Divider-free piecewise polynomial (Horner, any segmentation) | ≥ 50 Q (multiplier, multiplicand and accumulator coexist), ≥ +17% | `invsqrt_families_2026-11-19.md` |
| F3 | Absolute-frame streaming quotient without the CSWAP ring | ≥ 46 Q, +10 divider rows, ≥ +19%; history cannot be dropped in feed rows (lemma 5.2) | same |
| F4 | Table-free reciprocal-root recurrences (CORDIC, Newton, digit recurrence) | carry n, y and a residual at once: ≥ 74 Q, ≥ +200% | same |
| F5 | Coarse exponent-indexed phase table + fine part | 45–46 Q and +13 rotations: +9% to +13% | same |
| F6 | Logarithmic number system | two hard function evaluations instead of one: ≥ 55 Q, ≥ +120% | same |
| F7 | Direct Boolean synthesis of the phase bits | ≥ +490% | same |
| L5.2 | History-free rewrite of divider rows that take a dividend bit | impossible (no invariant interval when a bit is fed): feed rows must store their digits | same, lemma 5.2 |

## 5. Ideas

Units and conventions:
- **Gain** is the expected cut of the invsqrt task cost, in %, as a range. One Toffoli-equivalent is 0.059% of task
  cost at the 44-qubit frontier (0.057% on our 46-qubit build); one qubit is 2.3% (2.2% at 46).
- **Priority = P × (mid-point of gain) / effort**, effort in hours of focused work including the exact full-domain
  proof (`tools/prove_any.py` in the v12 tree) and a WSL verifier run.
- Gains are against the step before it in the same chain where stated (A02 → A03 → A04 form the path to the
  public frontier); other gains are on top of the frontier-equivalent build unless noted.
- **Test protocol for every idea:** build with `build.sh` (WSL, CPython 3.13 as the build script requires), run the
  generator's exact proof (`PASS_EXACT_NATIVE_TEXT_ALL_REQUESTED_INPUTS_ALL_HMR_TRANSCRIPTS`), then the official
  verifier at the default salt and at least two salts (the cost does not depend on the salt, but correctness on rare
  HMR branches does), then compare cost. The sampled anatomy script of section 2.2 gives the per-stage Toffoli delta
  in seconds and is the quick screen.

### A. Improve our code (v12 generator, ports of public improvements)

| ID | Idea | Mechanism (how it saves cost) | Est. gain % | P | Effort (h) | Priority | Risks / interactions | Sources |
|---|---|---|---|---:|---:|---:|---|---|
| INV-A01 | Submit the public frontier file `563d0583` now, credited, while A02–A04 rebuild it from source | Q 46 → 44 and 1,511.6 → 1,454.0 Toffolis: 80,162.67 → 74,560.91 | 7.0 | 0.95 | 1 | **6.65** | Compiled only, no generator: prize eligibility needs a generator that produces the exact circuit (host, topic 752010; topic 755488 on pinned compiled inputs still open). Leaderboard-only until A02–A04 land. Record SHA-256 `563d0583…9321`, Apache-2.0, NOTICE naming Syed Asad Ali first, then Özdağ, Özer, Panoras, manwithacat, jeffemerson, maoxin123. | maoxin123 *QCC 28.17*; Mizero Lucky Gall walkthrough |
| INV-A02 | Port Ali Haydar Özdağ's 64-row coefficient table into our generator | Table address 7 → 6 bits; both loaders re-synthesised (56 + 70 → 27 + 29 ANDs, LOAD1/LOAD2 in section 2.2); slope bit S4 measured out before the inverse product and rebuilt with 2 ANDs; angle LP refit. Measured on the 46-qubit build: −47 static Toffolis, 80,162.73 → 78,765.64 | 1.5–2.0 | 0.75 | 16 | **0.082** | Loader certificates were found by a private search (`ahoirg/qcc-invk6-k6c`); we must re-synthesise the 64-row loaders ourselves (SAT/ANF search under a 20/19-wire cap). Margin falls about 5×, raising the rotation charge by ~10 (already netted). Prerequisite of A03. Test: rebuild, compare LOAD1/LOAD2 stage counts and worst error/tolerance. | ahoirg *submit - verify and package* v25; RezOps 759220 |
| INV-A03 | Port Syed's v13/v14 46 → 44 moves into the source | Keep root bits v10 **and** v11 off the wires across the whole multiplier and rebuild them after (v10 from the top 4 bits of D plus v9); keep exponent flags e2, e3 off the wires while h is explicit; re-cap the loaders' working space for 44. Public effect: 78,765.64 → 75,170.62 | 4.3–4.8 | 0.5 | 30 | **0.077** | Only the compiled `47b31838` exists; the moves must be re-derived from jeffemerson's trace and a diff of the two files (stage map of section 2.2 helps). Every other 46-wide window must also drop by 2: 8,137 sampled ops sit at 46 today. Depends on A02 (the "slope table"). | jeffemerson *28.059 X-ray* §3; Syed v13 notebook |
| INV-A04 | Public micro-edit pack as generator passes | Measure-first dead-write fix (Panoras, −1), two AND merges on q28 (jeffemerson, −2), measure-out edits (manwithacat), CSWAP networks rewritten as CCX (frontier has 0 CSWAP), maoxin123's rescheduling, angle LP re-solve (Panoras/Özer). Public effect from `47b31838` to `563d0583`: 75,170.62 → 74,560.91 | 0.5–0.9 | 0.7 | 8 | **0.061** | Each edit is site-specific: re-run the scanners after A02/A03 rather than replay line patches. Must stay inside `pipeline.sh` so the build stays deterministic. | Panoras 27.30; jeffemerson 28.071; maoxin123 v27 |
| INV-A05 | Forward normaliser with measurement-conditioned clean-up, like the inverse | NORM costs 51 for 51 static Toffolis; NORM_INV costs 43.25 because 49 static gates are partly conditioned. Rewrite the forward prefix-AND erasures so their repair Toffolis run only on m = 1 | 0.2–0.5 | 0.4 | 6 | **0.023** | Patch 1 (`patch_fwd_narrow`) already measures prefixes early; the gain may be smaller. Touches the same stage as A14. Test: NORM stage mean in the anatomy script. | section 2.2 |
| INV-A06 | Domain-aware window re-synthesis scanner | For every window of ≤ 8–12 Toffoli-like gates, compute the exact reachable input/output relation over all 3,133,452 inputs and HMR branches, then SAT-minimise the AND count with unreachable states as don't-cares (this is how Aldemir's domain-aware CSWAP deletion and the q28 merge were found) | 0.2–1.0 | 0.5 | 16 | **0.019** | Chemistry loaders were all at their floor under the same test (758678), so expect small. Conditioned blocks need care (754152: apparent hits inside fix-up blocks were refuted on other branches). | 754152; Aldemir; jeffemerson |
| INV-A07 | SAT-minimal rebuild programs with domain don't-cares | Re-derive every "measure then rebuild" program (a6, v9, v10, v11, H2 3-AND, D12, S4 2-AND, history signs 13–15) as the minimum AND program over the reachable states only, and try rebuilding from different live wires | 0.1–0.5 | 0.5 | 8 | **0.019** | H2 already proved minimal at 3 (1 and 2 fail); signs 13–15 are SAT-minimal. Gains come from v10/v11 (new in A03) and from choosing other source wires. | Fortysix §6; HISTORY12 notes in `src/work/c10` |
| INV-A08 | LLM agent loop over the patch layer | Agents propose new patch modules (function replacements, as the 13 existing patches do); the evaluator is build + exact proof + cost; keep the Pareto set (Q, Toffolis). Directly reuses the generator's guard/cap machinery | 1–5 | 0.25 | 40 | **0.019** | Infrastructure cost is front-loaded; also serves A06/A07/A11/A17. Watch the token budget (forum 752858, 757915). | Gu et al. arXiv:2610.22367; PLAN.md §11 |
| INV-A09 | Re-synthesise the PHASE endpoint ring | 63 CSWAPs (part of PHASE's 105 Toffolis) rotate the quotient word cyclically by k. Search the minimum controlled-permutation network for the 11 legal k values, using the positions that are provably zero or constant for each k and sharing controls with decoded flags | 0.3–1.2 | 0.3 | 12 | **0.019** | Frontier already rewrote CSWAPs as CCX; measure its PHASE stage first. Interacts with A03 (e2, e3 off-wire changes which flag predicates exist). | Fortysix §4.6; 757030 problem 1 |
| INV-A10 | Clean up rotation bookkeeping | Fold the residual T gate and Clifford `Z_POW`s into merged features, check that no rotation pair acts on the same Boolean feature (frontier has 4 T-like gates, v12 has 1) | 0.01–0.06 | 0.6 | 1 | **0.018** | Tiny; do it as part of A04. | verifier |
| INV-A11 | Remove the inverse-product premium | SEED_CLR costs 130 against 103 for the forward product (+27): the inverse product regenerates conditionally retired data. Try holding the retired rounding carries in idle flag wires off-peak, or a different inverse multiplier (e.g. subtract partial products in reverse order with measurement-conditioned regeneration) | 0.4–1.5 | 0.3 | 16 | **0.018** | The premium buys width (patches 5, 8, 9); any saving must keep the 44 cap. Shares wires with A03's hidden root bits. | section 2.2; Fortysix §6 |
| INV-A12 | Margin-maximising joint refit (MILP) | Choose table word roundings (floor/ceil of each intercept/slope word, 64 × 2 words) together with the 33 angles to maximise `min_margin`; each doubling of the margin saves 4.6 Toffoli-eq of rotation charge. The frontier sits at worst error/tolerance 0.99956, so there is room | 0.15–0.6 | 0.5 | 10 | **0.018** | Table words feed the loaders: rounding changes loader cost (do jointly with A18). The phase is exactly linear in the angles (Panoras) but not in the words, so use MILP or local search with an exact LP inner loop. | Panoras 28.12769; manwithacat 753185 |
| INV-A13 | Trim the consuming root | ROOT costs 96 each way. Use exact range analysis on reachable states to drop always-zero carry bits in early rows, exploit the known top bit (n ≥ 2^20), and try conditionally clean carries for the top rows | 0.2–0.8 | 0.35 | 12 | **0.015** | Syed already went 99 → 96 per direction (RN3 "source-aware root"). Root inverse must mirror it exactly. | `research/root_normalization/` in v12 tree |
| INV-A14 | Minimum-cost even-shift normaliser | Treat forward (51) and inverse (43.25) normalisation as a controlled-permutation synthesis problem for 11 shift amounts on 22 bits with legal-domain don't-cares (top bits of x are never both 1); SAT with symmetry breaking, layered by shift bit | 0.3–1.0 | 0.35 | 16 | **0.014** | Open problem 1 of 757030: an 11-hour CaDiCaL run on a similar instance was inconclusive. Syed's RN2/RN2-B plans (leading-one route, −1) are unexecuted ideas in the tree worth reading first. | 757030; `RN2_LEADING_ONE_PLAN.md` |
| INV-A15 | Measurement-conditioned repair instead of reverse computation | Where a multi-Toffoli temporary is erased by recomputation, HMR it and run the repair as a phase oracle only when m = 1: about half the Toffolis on average. Candidates: CSRC (3), R_PAID (18), the parts of DIV_CLEAN that rebuild values only to erase them | 0.2–0.8 | 0.3 | 12 | **0.013** | Only applies to values not needed again; in-place stages are excluded (N10). Increases HMR count (no cost). | Gidney 1709.06648; Fortysix §4.4 |
| INV-A16 | Replace the last divider row by a rotation on the remainder sign | The final rows' quotient digits are remainder signs; skipping the last forward/reverse row (~11 + ~13 Toffolis) and rotating the final remainder's top bit costs one rotation (~7) if the error budget allows | 0.3–1.0 | 0.15 | 8 | **0.012** | jeffemerson: "late-row truncation is already at its edge"; the margin is near 0.9996, so probably infeasible unless A12 frees margin. | jeffemerson #7 |
| INV-A17 | Reduce the divider clean-up premium | DIV_CLEAN costs 323 + 5 against 270 forward. Reorder sign retirements, rebuild signs from cheaper predicates, borrow idle (dirty) wires for reverse-row carries via conditionally clean ancillae | 0.3–1.8 | 0.25 | 24 | **0.011** | jeffemerson found ≤ 5 Toffolis in the current layout. The premium is what pays for the width; must not reopen a 45-wide peak. | jeffemerson direction 4; Khattar & Gidney 2407.17966 |
| INV-A18 | Loader-aware table rounding | Pick floor/ceil (and the XOR frame) of each 64-row word to minimise the AND count of the two loaders (27 + 29), under the exact error gate | 0.2–0.7 | 0.35 | 16 | **0.010** | Same idea as chemistry problem 5 of 757030. Interacts with A12 (margin) and A02. | 757030 |
| INV-A19 | Table segmentation sweep | Non-uniform rows over v (Sanders-style uneven segments), different intercept/slope bit split, or a table indexed by v and one flag; target the same margin with fewer loader ANDs | 0.3–1.2 | 0.25 | 20 | **0.009** | 32 rows fail (N5); 64 uniform is the known optimum for this arithmetic, so only non-uniform splits are new. | Sanders et al. 2007.07391; jeffemerson #2 |
| INV-A20 | Our own rescheduler pass | Move each HMR to just after the last read and each first write to just before the first use (maoxin123 cut ops 34,789 → 20,659 and full-width ops 8,137 → 4,299) | 0–0.3 | 0.6 | 10 | **0.009** | Small direct gain; its value is enabling A21/A22/A24 by shrinking the set of peak windows. | maoxin123; jeffemerson direction 3 |
| INV-A21 | Hide more root bits across the multiplier | Extend Syed's v10/v11 trick to v9/v8 (rebuild from D, v-low and h with short AND programs) to empty the multiplier windows that block 43 | 0.5–2.3 | 0.15 | 24 | **0.009** | Needs a rebuild of ≤ ~15 Toffolis per bit; the bits get less redundant as you go down. Pairs with A22 and A24. | jeffemerson 28.059 X-ray |
| INV-A22 | Cubic / basis-change parking census | Extend the X-ray parking test to degree-3 formulas and affine basis changes of the live wires; find a cover of the 88 uncovered 44-wide windows (74 in the multiplier) | 0–2.3 | 0.1 | 16 | **0.008** | Analysis first; only pays if a cover under ~39 Toffolis exists. jeffemerson rates the odds at ~5%. | jeffemerson open question 1 |
| INV-A23 | Share address logic between the two table loads | LOAD1 and LOAD2 use the same address; keep part of the unary-iteration AND tree (or its measured parities) between them where width allows, or reload from a cheaper frame | 0.1–0.5 | 0.25 | 10 | **0.008** | The two loads are 10,000 ops apart at full width; holding anything costs width. The held-products pass may already do the feasible part. | section 2.2 |
| INV-A24 | 44 → 43 by removing peak carries | Every 44-wide window is the top of a carry staircase (one carry on a 43-wide body). Replace that carry by a conditionally clean or dirty borrowed ancilla, or compute the top sum bit without a carry-out; needs < ~39 Toffolis in total | 1.0–2.3 net | 0.12 | 40 | **0.005** | Panoras: break-even +38, an older attempt cost +257. Must be combined with A20/A21/A22 to cover all windows; the loader peaks need their own fix (A25 or a narrower loader). | jeffemerson; Khattar & Gidney; Takahashi et al. |
| INV-A25 | Split the coefficient load | Load C as two 64-row half-words (7 + 6 bits; loader cost barely depends on word width), run the product in two passes and HMR each half after use: about −6 wires in the multiplier window for ~+60–120 Toffolis | −2 to +4 | 0.12 | 40 | **0.005** | Pays only if the divider and loaders also fit at 43 or 42 (A24). Differs from N15 (reload per digit, +1,100). | jeffemerson #12, #11 |

### B. Rewrite with a different technique

| ID | Idea | Mechanism (how it saves cost) | Est. gain % | P | Effort (h) | Priority | Risks / interactions | Sources |
|---|---|---|---|---:|---:|---:|---|---|
| INV-B01 | Evolutionary LLM search over algorithm families | Run a Gu-et-al.-style loop whose evaluator is our exact whole-domain simulator plus the verifier's cost formula, seeded with the three known families (Newton reference, 62 Q table + Horner, 46/44 Q root + divider) and asked for Q ≤ 40 designs | 0–25 | 0.08 | 80 | **0.008** | Expensive and open-ended; the leaders' private 28%-class gain (PLAN.md §1) suggests a new family exists somewhere. Shares infrastructure with A08. | arXiv:2610.22367; 2611.09822 |
| INV-B02 | Fused reciprocal-square-root digit recurrence | Replace root + table + product + divider by one digit recurrence for 1/√n (Ercegovac–Lang) that consumes n in place and emits phase digits, removing the 192-Toffoli root and the two table loads | 5–20 | 0.05 | 80 | **0.0075** | N17: naive recurrences hold three ~26-bit numbers (≥ 74 Q). Only worth it if a two-register formulation exists (jeffemerson direction 5). | Ercegovac & Lang ch. 6–7 |
| INV-B03 | Compressed body encoding | (v, h, k) uses 26 wires for 22 bits of information. For k ≥ 1 (x < 2^20) the root/residual have spare range; encode the flags inside it and decode only where the ring needs k. Generalises Syed's "e2, e3 off the wires" | 2.3–4.6 | 0.1 | 60 | **0.006** | N11/N16 show the naive versions fail; the divider consumes h. Must keep every stage's arithmetic cheap on the encoded form. | Fortysix §5; jeffemerson #8, #13 |
| INV-B04 | Odd-shift normalisation with parity-selected table | Shift by any j (not only even), so n ∈ [2^21, 2^22) with an implicit top bit; fold the √2 factor for odd j into the table (one extra address bit) and the dividend constant. The range of s halves, so fewer rows or more margin; the implicit top bit frees a body wire, the parity flag costs one | 0–3 | 0.15 | 40 | **0.006** | Ring weights must stay powers of two (the √2 must live in the dividend, not the angles). Syed's RN2 plans already handle the leading one; check them first. | Häner et al. 1807.02023; RN2 plans |
| INV-B05 | Tolerance-approximate carries | Truncate carry chains (Gidney-style runways / approximate adders) where the induced phase error provably stays inside the budget on every legal x, removing carry wires at peaks | 0.5–3 | 0.1 | 30 | **0.005** | Budget is nearly used (0.99956): needs A12 first. Exhaustive verifier, so "rare" errors are not allowed. | Gidney 1905.08488 |
| INV-B06 | History-free reverse divider (hybrid) | With an odd divisor, a *restoring* row has q_i = r_i mod 2, so its reverse needs no history. Use restoring rows only where the history wires form the peak (the last rows), non-restoring elsewhere | 0–3 | 0.05 | 40 | **0.004** | *Ours, priced on paper:* a restoring row costs about 3n Toffolis forward and 2n reverse (comparator + controlled subtract, since the ± mask trick of non-restoring rows is lost) against about n + n·1.2 today, so a full switch adds ~+900 Toffolis (loses). Only a few rows can switch. | Thapliyal 1809.09732 |
| INV-B07 | Division-free Goldschmidt with streamed phase | 1/s = T(a)·(1 − e + e² − …) with e = s·T(a) − 1 ≈ 2^-6: kick each partial product's bits into the phase as they are formed instead of holding a 24-bit quotient | −5 to +8 | 0.08 | 50 | **0.003** | N9 priced "reciprocal by multiply" above the incumbent; streaming to the phase raises the rotation count (7 Toffoli-eq each). | jeffemerson #5 |
| INV-B08 | Narrow Horner hybrid | Return to the 62 Q line's table + Horner arithmetic (656 Toffolis, close to the divider's 593) but with the 46 Q discipline: load, use, HMR and reload the table words instead of holding 36 table wires | −5 to +10 | 0.07 | 60 | **0.003** | The Horner chain carries a 17-bit slope + 17-bit constant + 14 digit wires; may not fit below 44 even with reloads. | Fortysix §5.2; Tangfastic |
| INV-B09 | Coarser normalisation (shift by 4^2k) | Use 3 flag wires (k' = 0..5) and let the table/root cover n ∈ [2^18, 2^22): −1 flag wire for a 2× wider s range (more rows or one more divider row) | −1 to +2 | 0.1 | 30 | **0.002** | Probably a wash: +~40 Toffolis vs one qubit (~39). Only interesting combined with B04. | ours |

## 6. Combination notes

- **The frontier path is one chain: A02 → A03 → A04.** A03 needs A02's slope table (v10 is rebuilt from the top
  bits of D) and A04's micro-edits are site-specific, so re-run the scanners after each structural step instead of
  replaying the public line patches. A01 is the stop-gap for the leaderboard while the chain is rebuilt.
- **Width ideas must be tested together.** A21, A22, A24 and A25 each remove only part of the 44-wide windows
  (loaders, multiplier, divider rows, ring). A 43-qubit build needs a cover of every window, so prototype them on one
  branch with A20's rescheduler first (it halves the windows). Reaching 43 is worth ~39 Toffolis; each move's
  Toffoli price must be summed against that single budget.
- **Margin is a shared resource.** A02 (smaller table), A12 (refit), A16 (row truncation), A18 (loader-aware
  rounding), A19 (segmentation) and B05 (approximate carries) all spend or create rotation margin. Run A12 last in
  any combination, and evaluate A16/B05 only after A12 has maximised the margin.
- **Table design ideas are mutually exclusive variants:** A02 (uniform 64 rows), A19 (non-uniform), B04
  (parity-selected), B09 (wider range). Pick one table family, then apply A18 and A12 on it.
- **Divider ideas touch the same rows:** A09 (ring), A16 (last row), A17 (clean-up premium) and B06 (restoring tail)
  all edit DIV_FWD/PHASE/DIV_CLEAN; A03 also changes which flags are live there. Do A03 first.
- **Normaliser ideas:** A05 and A14 both rewrite NORM; A14 supersedes A05 if it succeeds. A13 (root) is independent
  of both but shares the "known top bit" fact with B04.
- **Product ideas:** A11 (inverse premium), A25 (split load) and A21 (hidden root bits) change the multiplier's
  register set together; do not combine A11 with A25 until each is proved alone.
- **Search infrastructure is shared:** the exact whole-domain simulator (section 2.2 script, extended with phase
  tracking) feeds A06, A07, A22 and A08/B01. Building it once in WSL (CPython 3.13, as `build.sh` requires) is the
  first engineering task.
- **Provenance:** anything derived from `563d0583`, `47b31838`, `4560ac0c` or the v12 tree must carry Syed Asad Ali's
  NOTICE and name him first; for prize eligibility each must be regenerated by our own generator (host answer in
  752010).

