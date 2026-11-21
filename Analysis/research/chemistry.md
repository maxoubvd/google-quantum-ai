# Chemistry (dfthcsos): cost anatomy and improvement ideas

Date: 2026-11-17. Scope: the `dfthcsos` task (DF-THC SOS block encoding of FeMoCo76o113e).
Base: our rebuild of Syed Asad Ali's v12 chemistry (`artifacts/v12/dfthcsos/`), cost 296,635,967.59.
Public frontier: about 2.94642e8 (maoxin123's 4,022-Toffoli loaders + Atakan Aldemir's weights `703843d3`).

## 1. Summary

1. Our cost 2.9664e8 = lambda_eff 46.907 x T 14,810 x Q 427. T is 55% loaders (8,148), 28% inner state
   prep (4,095), 15% Givens (2,250), 2% outer prep. Q is 71% the inner-QROAM temp (303), set by the shape.
2. lambda_eff is 0.25% above a hard floor (46.79); E_gap and both errors already sit on their gates.
3. The cost model has a cliff: (8,143,44) at brot 8 keeps the inner QROAM at k2 = 4 by one qubit
   (600 is not < 600). B <= 128 or brot 9 jump Q by +64%. B in 129..143 is the safe band.
4. The cost model prizes small B (0.4x-0.7x our cost), but we measured that FeMoCo's leading pair-space
   factors are full rank: B < 76 is off by orders of magnitude. That door is closed.
5. Fastest win: adopt the public frontier parts (CHEM-A01, -0.67%, 2 h), then rebuild them from source.
6. Main controllable lever: loader Toffolis via word codes (bit-level Reed-Muller constraints, sensitivity-
   weighted and mixed precision), worth 1-4% if the 20 mHa budget allows.
7. Shape moves: C 43 (-0.5%), C 42 (-1.1%), B 139 (-0.4%), R 7 (-2%), brot 7 (-5% to -6%) (re-priced by S2.23, see §2.3); all are capacity-
   or accuracy-limited and need a full rebuild on a P100.
8. A free gauge nobody has used publicly: permuting rows within groups (and groups, and one-body rows).
9. Realistic combined target: -2.5 to -3.5% task cost (+0.8 to +1.2% score). 31 ideas below.

## 2. Cost anatomy

All numbers below come from the scorer code (`verify.py`, `scoreinfra/costs_prepare.py`,
`scoreinfra/dfthc_math.py`) evaluated on our `artifacts/v12/dfthcsos/submission.npz`
(scripts run locally with the Windows `.venv`; `evaluate_words` on the stored `angle_words` at brot 8).

### 2.1 The formula

```
score      = lambda_eff * T_total * Q_total
lambda_eff = sqrt(E_gap * (2 Lambda - E_gap))
T_total    = 2 * ANGLE_per_pass                      (measured mean Toffolis, forward + reverse)
           + outer_state_prep(R,B,C,chi=10)          (SpinFreeData_V1, includes the dagger)
           + inner_state_prep(R,B,C,chi=10,beta)     (SpinFreeData_V1, includes the dagger)
           + 2 * (N-1) * (2*brot - 1)                (Givens rotations, both passes)
Q_total    = outer_persistent + inner_persistent     (SpinFreeData_V1)
           + J * brot                                (word registers)
           + max(outer_temp, inner_temp, angles_temp)
angles_temp = max_k num_qubits(ANGLE_k) - |selection| - J*brot
```

Note: the `selection` register (11 qubits at 1,220 addresses) is **not** counted in `Q_total`; it is
treated as the state-preparation's address register.

### 2.2 Current values (our v12 build, (R,B,C) = (8,143,44), brot 8, J 5, 15 circuits)

| Term | Formula detail | Value | Share |
|---|---|---:|---:|
| ANGLE Toffolis | 2 x 4,074 per pass (15 loaders, 263-293 each) | **8,148** | 55.0% of T |
| inner_state_prep | for_qroam `RC*ceil((B+1)/2^k2) + 2^k2*b2` = 352*9 + 16*19 = 3,472; for_qroam_dag `QR(RC,B+1) + (B+1) + N + b2` = 320+144+76+19 = 559; uniform 2*4*ceil(log2(B+1)) = 64 | **4,095** | 27.7% of T |
| givens_rotations | 2*75*(2*8-1) | **2,250** | 15.2% of T |
| outer_state_prep | X_o = N+RC = 428; 2*4*9 + QR(428,19)=164 +19 + QI(428)=42 +20 | **317** | 2.1% of T |
| **T_total** | | **14,810** | |
| inner temporary | `ceil(log2(RC*ceil((B+1)/2^k2))) + (2^k2-1)*b2` = 12 + 15*19 = 297, + uniform temp 6 | **303** | 71.0% of Q |
| inner persistent | b2 + chi + ceil(log2 R) + 2 = 34, + uniform 10 | 44 | 10.3% of Q |
| outer persistent | b1 + chi = 29, + uniform 11 | 40 | 9.4% of Q |
| word registers | J*brot = 5*8 | 40 | 9.4% of Q |
| outer temporary | (not binding) | 90 | - |
| angles temporary | <= 128 - 11 - 40 = 77 (not binding) | <= 77 | - |
| **Q_total** | 84 + 40 + max(303, 90, <=77) | **427** | |
| Lambda | lambda_h1 + lambda_h2 | 282.9736 | |
| E_gap | floor is 3.914898 (strict) | 3.9148982 | at the floor (+2e-7) |
| error_1 / error_2 | gate 20 mHa | 19.610 / 19.9988 mHa | both at the edge |
| **lambda_eff** | | **46.9073** | floor ~46.79 |
| **score** | | **296,635,967.6** | |

Here `b2 = ceil(log2 B) + chi + 1 = 19`, `chi = KEEP_BITS = 10` (fixed), and `k2` is the inner QROAM
block exponent chosen by the model (see 2.4).

### 2.3 Exchange rates (local derivatives at our point)

| Change | Task cost | Score (geo-mean) |
|---|---:|---:|
| -1 ANGLE Toffoli per pass | -0.0135% | -0.0045% |
| -52 Toffolis per pass (our 4,074 -> public 4,022) | -0.70% | +0.23% |
| -1 qubit of Q_total | -0.234% | +0.078% |
| -1% lambda_eff | -1% | +0.33% |
| C 44 -> 43 (cost model only, same loaders/row) | -0.53% | +0.18% |
| C 44 -> 42 | -1.06% | +0.35% |
| C 44 -> 40 | -2.12% | +0.71% |
| B 143 -> 139 (R=8, C=44; fewer rows, same k2) | -1.5% | +0.5% |
| B 143 -> 135 | -3.0% | +1.0% |
| B 143 -> 129 (lowest B that keeps k2=4) | -5.2% | +1.8% |
| R 8 -> 7 with B=143, C=48 | -7.5% | +2.6% |
| J 5 -> 4 at equal ANGLE Toffolis | -1.87% | +0.63% |
| J 5 -> 3 at equal ANGLE Toffolis | -3.75% | +1.3% |
| brot 8 -> 7 (cost model + 7/8 of loader bits) | -10.0% | +3.5% |
| B 143 -> 144 (ceil((B+1)/16) steps from 9 to 10) | +2.75% | -0.9% |
| B 143 -> 128 (k2 jumps 4 -> 5) | +41.7% | catastrophic |
| brot 8 -> 9 (k2 jumps 4 -> 5) | +74.6% | catastrophic |

The shape rows assume the loaders keep the current 3.34 Toffolis per row per pass (4,074/1,220); only
the ANGLE share scales with the row count `76 + R*B`.

**Correction (S2.23, 19 Nov, `chemistry_families_2026-11-19.md` §0.3–0.4):** that assumption is wrong. The sum of
per-plane core ranks of the loaded bit-planes (`tools/core_rank.py`) is invariant under row removal (8,443 / 8,444 /
8,445 / 8,395 for 1,220 / 1,077 / 934 / 512 rows); only the per-circuit address caches scale with rows. Loader T per
pass ≈ 0.42 × Σ core rank + K × ⌈rows/16⌉, so one row ≈ 1.9 T (−0.013%), not 3.34–6.6. Re-priced: B 143 → 129
−1.6% (not −5.2%), R 8 → 7 −2.0% (not −7.5%), R 8 → 6 −3.5%, brot 8 → 7 −4.9% to −6.4% (not −10%). The C rows
and the cost-model-only rows are unchanged.

### 2.4 The k2 cliff: why (8,143,44) at brot 8 is a sweet spot of the cost model

`find_qroam_size_to_minimize_global_temp_ancilla` starts from the Toffoli-optimal QROAM exponent
`k2_opt = QR(R*C*B, b2)[0]` (= 6 here) and takes the largest `k <= k2_opt` such that

```
ceil(log2(R*C*ceil(B/2^k))) + (2^k - 1)*b2  <  beta*(N-1) = 75*brot
```

At brot 8 the right side is 600. With `b2 = 19` (B in 129..256): `k=5` gives `ceil(log2(352*5)) + 31*19 =
11 + 589 = 600`, which is **not** `< 600`, so `k2 = 4` and the inner temp is 303. Any change that lets
`k = 5` pass makes the inner temp ~574-606 and Q_total ~700 (+64%):

- `b2` dropping to 18 (B <= 128): `31*18 + 11 = 569 < 600` -> k2 = 5.
- brot 9 (right side 675): k2 = 5. brot 7 (525): still k2 = 4 (fine).
- `R*C*ceil(B/32) <= 1024` (log term 10): k2 = 5.

Conversely the Toffoli side uses `ceil((B+1)/2^k2)`, so `B+1 = 144 = 9*16` is exactly tight at k2 = 4:
B = 144 costs +407 Toffolis. So, inside the B >= 129 regime, B in [129, 143] with B+1 close to a
multiple of 16 (B = 143, 127 is excluded) is optimal, and the current shape already sits on it. **The
cost model punishes every "nearby" shape change except dropping C, dropping R, or dropping B toward 129.**

Second inconsistency worth knowing (not exploitable as far as we can see): `inner_state_prep` lowers
`k2` when `B < 2^k2`, but `inner_state_prep_qubit_cost` does not, so for B < 2^k2_opt the qubits are
priced for a larger QROAM than the Toffolis.

### 2.5 The small-B regime the cost model loves, and why the physics forbids it

A grid over `RB in [900,1500]`, `RC in [290,520]` finds the cost model's minimum at tiny B (e.g.
(290,4,1) at 0.40x our cost, (143,8,3) at 0.64x, (72,16,5) at 0.68x, (59,20,6) at 0.72x): small B gives
small b2, a forced small k2, a tiny inner temp (~100-250 qubits) and a tiny inner QROAM. Between
B = 24 and B = 128 the cost model is 1.2x-1.5x worse (k2 = 5 forced).

We checked the physics with a quick experiment (pair-space eigendecomposition of the shifted target
`tei + (S x I + I x S)/2`, then rank-r truncation of each leaf, i.e. a DF-like (L, r, 1) shape):

| Check | Result |
|---|---|
| Eckart-Young floor of error_2, rank L in pair space | L=290: 19.9 mHa; 300: 16.7; 330: 10.0; 352: 7.4; 400: 4.1; 450: 2.3 |
| DF truncation, rank r per leaf (any L >= 290) | r=4: 10.8 Ha; r=16: 9.0 Ha; r=32: 7.1 Ha; r=64: 2.5 Ha; r=76: = EY floor |

The leading pair-space eigenvectors are full-rank 76x76 matrices, so any shape whose X^(rho c) have
rank < ~76 is off by orders of magnitude. DF-THC needs B >= ~76 per group (the current B = 143 is an
overcomplete THC-like dictionary shared by C = 44 weightings). Optimized fits beat naive truncation,
but not by 100x. **Conclusion: the B < 76 regime is closed; the realistic shape moves are C down, B down
to 129..142, or R 8 -> 7, all of which are capacity-limited.** The EY numbers also say RC = 352 leaves
~7 mHa of floor, and RC = 336 (C = 42) ~10 mHa, before the dictionary restriction and 8-bit rounding.

### 2.6 lambda_eff is pinned

From `dfthc_math.py`: `E_gap = const + 2*sum|neg eig(h1')| + 0.5*sum W^2` with `const = -559.13 Ha`, and
`Lambda >= sum pos eig + sum|neg eig| + 0.25*sum W^2` (triangle inequality on lambda_h2). Hence
`lambda_eff^2 >= E_gap * (E_gap + 559.13 - E_gap) = 559.13 * E_gap >= 559.13 * 3.914898`, i.e.
**lambda_eff >= 46.79 for every factorization** (Lynn Langit, topic 753633; re-derived here). We are at
46.907, 0.25% above. Our E_gap is already at the floor (3.9148982). The only lambda_eff slack left is the
triangle slack in lambda_h2: 60.5% of our (rho,c) columns have mixed-sign w_b. Closing all of it is worth
at most 0.25% of task cost.

### 2.7 What actually moves the cost (summary of the anatomy)

1. **ANGLE loaders: 55% of T.** Every Toffoli per pass is 0.0135%. This is the only term fully under our
   control without retraining. The public frontier is 52 per pass (0.70%) below us.
2. **Shape (R,B,C):** the inner QROAM (27.7% of T, 71% of Q) is fixed by the shape. C -1 = -0.53%;
   B -1 (in 129..143) ~ -0.37%; R -1 ~ -7%. Needs new fits at the 20 mHa edge.
3. **brot:** 7 bits would be -10% but every public attempt failed accuracy (error_2 110-350 mHa).
4. **J:** -1.87% per word removed, paid back by more loader overhead (break-even +138 Toffolis per pass
   for J=4).
5. **lambda_eff:** at most -0.25%.

## 3. State of the art

### 3.1 Papers (factorizations, block encodings, loaders)

| Reference | Why it matters here |
|---|---|
| Low, King, Berry, Han, DePrince, White, Babbush, Somma, Rubin, *Fast quantum simulation of electronic structure by spectrum amplification*, PRX 15, 041016 (2025), [arXiv:2502.15882](https://arxiv.org/abs/2502.15882); local `Scientific Papers/PRX 2025.pdf` | Defines the task: DF-THC in a spin-free SOS form, BLISS-type shifts, cost ~ sqrt(lambda * E_gap) (our `lambda_eff`), and the SpinFreeData cost model the verifier copies. |
| Competition paper, `google-verifiable-quantum-circuit-compilation/verifiable_circuit_compilation.pdf` (Sec. III.C, App. A.4, Eq. 72) | ANGLE_k XOR contract, the forward/reverse chain, why J*brot is persistent. |
| Lee, Berry, Gidney, Huggins, McClean, Wiebe, Babbush, *Even more efficient quantum computations of chemistry through tensor hypercontraction*, PRX Quantum 2, 030305 (2021), [arXiv:2011.03494](https://arxiv.org/abs/2011.03494) | THC block encoding, the QROAM/inner-prep costing style, FeMoCo (Reiher and Li Hamiltonians) baselines. |
| von Burg, Low, Haener, Steiger, Reiher, Roetteler, Troyer, *Quantum computing enhanced computational catalysis*, PRR 3, 033055 (2021), [arXiv:2007.14460](https://arxiv.org/abs/2007.14460) | Double factorization (DF) with Givens-rotation basis changes loaded by QROM: the ancestor of our ANGLE loaders. |
| Berry, Gidney, Motta, McClean, Babbush, *Qubitization of arbitrary basis quantum chemistry leveraging sparsity and low rank factorization*, Quantum 3, 208 (2019), [arXiv:1902.02134](https://arxiv.org/abs/1902.02134) | QROAM costs `ceil(L/k) + M(k-1)`, the exact form of `QR`/`QI` in `costs_prepare.py`; rotation-angle loading. |
| Babbush, Gidney, Berry, Wiebe, McClean, Paler, Fowler, Neven, *Encoding electronic spectra in quantum circuits with linear T complexity*, PRX 8, 041015 (2018), [arXiv:1805.03662](https://arxiv.org/abs/1805.03662) | Unary-iteration QROM (the reference loaders), measurement-based uncomputation. |
| Low, Kliuchnikov, Schaeffer, *Trading T-gates for dirty qubits in state preparation and unitary synthesis*, [arXiv:1812.00954](https://arxiv.org/abs/1812.00954) | SelectSwap/QROAM with clean or dirty qubits: background for using the word registers as borrowed workspace (B14). |
| Motlagh, Pocrnic, *Halving the cost of QROM*, [arXiv:2605.20334](https://arxiv.org/abs/2605.20334) (2026) | SelectCopy, about 2x cheaper than SelectSwap in the dirty-qubit regime. Our loaders are no longer plain QROMs, but the dirty-borrowing idea may transfer to the cache (B14). |
| *Database Reordering for Compact Grover Oracles*, [arXiv:2604.06578](https://arxiv.org/abs/2604.06578) (2026) | Permuting the address assignment of a lookup table before ESOP synthesis. Our (rho, b) row order is a free gauge (B03). |
| Gidney, *Halving the cost of quantum addition*, Quantum 2, 74 (2018), [arXiv:1709.06648](https://arxiv.org/abs/1709.06648) | Temporary AND + HMR/CZ uncomputation, used by every public loader (0 Toffolis to uncompute). |
| Meuli, Soeken, De Micheli, *Xor-And-Inverter graphs for quantum compilation*, npj QI 8, 7 (2022), [doi:10.1038/s41534-021-00514-y](https://doi.org/10.1038/s41534-021-00514-y); Meuli et al., *The role of multiplicative complexity in compiling low T-count oracle circuits*, [arXiv:1908.01609](https://arxiv.org/abs/1908.01609) | Loader cost = multiplicative complexity of a 40-output Boolean function; XAG rewriting, exact MC synthesis (B13, A12). |
| Lovasz, *Matroid matching and some applications*, JCTB 28 (1980) | Rank formula behind Syed Asad Ali's exact product-pairing (matroid parity) planner. |
| Loaiza, Izmaylov, *Block-invariant symmetry shift (BLISS)*, JCTC 2023, [arXiv:2304.13772](https://arxiv.org/abs/2304.13772) | The S, beta_1 shift; the verifier fixes it to the baseline's, so this lever is closed for us. |
| Cohn, Motta, Parrish, *Quantum filter diagonalization with compressed double-factorized Hamiltonians*, PRX Quantum 2, 040352 (2021), [arXiv:2104.08957](https://arxiv.org/abs/2104.08957) | Compressed DF (shared leaf bases fitted by optimization), close to the (R, B, C) structure. |
| Oumarou, Scheurer, Parrish, Hohenstein, Gogolin, *Regularized compressed double factorization*, Quantum 8, 1371 (2024), [arXiv:2212.07957](https://arxiv.org/abs/2212.07957) | L1-regularized CDF fits with lambda below THC: candidate initializer for lower-C tables (B06). |
| Rocca et al., *Symmetry-compressed double factorization*, JCTC 2024, [arXiv:2403.03502](https://arxiv.org/abs/2403.03502); Deka, Zak, [arXiv:2412.01338](https://arxiv.org/abs/2412.01338) | CDF + symmetry shift; a conference abstract reports about 50% runtime reduction against THC for FeMoCo. Initializer ideas only (our shift is fixed). |
| Reiher, Wiebe, Svore, Wecker, Troyer, PNAS 114, 7555 (2017), [arXiv:1605.03590](https://arxiv.org/abs/1605.03590) | FeMoCo as benchmark. |
| Li, Li, Dattani, Umrigar, Chan, *The electronic complexity of the ground-state of the FeMo cofactor of nitrogenase*, JCP 150, 024302 (2019), [arXiv:1809.10307](https://arxiv.org/abs/1809.10307); data [zhendongli2008/Active-space-model-for-FeMoco](https://github.com/zhendongli2008/Active-space-model-for-FeMoco) | Source of the 76o/113e FCIDUMP (byte-identical, forum 753185). Orbitals are split-localized and DMRG-ordered. |
| Nagel et al., *AdaRound*, [arXiv:2004.10568](https://arxiv.org/abs/2004.10568); Frantar et al., *GPTQ*, [arXiv:2210.17323](https://arxiv.org/abs/2210.17323); Babai nearest-plane rounding | Quantization-aware rounding. GPTQ was measured as a null here (753185). Syed Asad Ali used Babai rounding + augmented-Lagrangian refits. |

### 3.2 Public competition work (chemistry lineage, as of 2026-11-17)

| Item | Author | License | What it brought | Where |
|---|---|---|---|---|
| (8,143,44) table with codes, 4,025-Toffoli loaders, cost 2.948e8 | Lera (marculera) | Apache-2.0 | Codes on the XOR of the five words (bits 1, 2), bit-0 affine codes, sign-negation gauge lifts the pi-gauge limit | [dataset](https://www.kaggle.com/datasets/marculera/dfthc-c44-build), [notebook](https://www.kaggle.com/code/marculera/df-thc-11-38x-codes-on-the-xor-of-five-words), local `downloads/datasets/marculera__dfthc-c44-build/NOTICE` |
| 4,022-Toffoli loaders on Lera's table (frontier) | maoxin123 | CC0 (compiled only) | First two loaders + a free padding row | [notebook](https://www.kaggle.com/code/maoxin123/qcc-28-17-4-022-toffoli-chemistry-loaders) |
| Weights `703843d3` (lambda_eff about 46.92) | Atakan Aldemir | Apache-2.0 | Fixed-word refit to the gates | local `downloads/datasets/atakanaldemir__qcc-shared-4022-loaders-and-own703-weights/` |
| Coded-word chemistry v25/v27 (4,044 -> 4,038) | maoxin123 | CC0 | Reed-Muller codes on every batch from 3; full 75-column top-bit gauge parks 5 quadratic low monomials | [v25](https://www.kaggle.com/code/maoxin123/qcc-v25-coded-word-chemistry) |
| Own (8,143,44) table, compiler, 4,074 loaders (our base) | Syed Asad Ali | Apache-2.0 | Babai rounding, matroid-parity pairing, pi-gauge helper carry, exact all-branch proof | `third_party/syed3000-vqcc-own-circuit-portfolio-v12/dfthcsos/dfthcsos_source/` |
| Compact groups from the Hamiltonian, (9,143,43) / (8,143,45) | Lera | Apache-2.0 | Pivoted-Cholesky leaves clustered by orbital as the seed; the shape family everyone uses | [notebook](https://www.kaggle.com/code/marculera/df-thc-10x-compact-groups-from-the-hamiltonian) |
| Pi gauge (theta + pi is free under a sign-parity condition) | Atakan Aldemir | Apache-2.0 | New tree parametrization, -96 Toffolis on its day | [notebook](https://www.kaggle.com/code/atakanaldemir/qcc-gap-3-residue-repair-and-c43-native-compiler) |
| Four word forms (theta, theta+pi, -theta, pi-theta), product pattern for five dirty words in 128 qubits | gengsr | Apache-2.0 | | per Lera's NOTICE |
| Staged progressive rounding, sign-mixing measure D(w), exact offset gauge | Ahmed Berat Ozer | Apache-2.0 | | [notebook](https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression) |
| Constrained refit, column gauge, `dfthc_refit` | James Barlow (manwithacat) | Apache-2.0 | | [Zenodo write-up](https://doi.org/10.5281/zenodo.23930466) |
| Exact piecewise-trigonometric column gauge | Pablo Diego Pereira | | Lambda_2 polish at 1e-6 relative | forum 755525 |
| Refit code `refit_grg.py`, `chemlib.py` | prvsiyan | Apache-2.0 | GRG refit; shows the weights sit at a local optimum | [dataset](https://www.kaggle.com/datasets/prvsiyan/vqcc-verified-best-of-circuits-1106a) |
| Cost-model scans of B and (R, C) | Asterios Panoras | | "Fully mixed products": counting floor about 25 products per loader below the current method | [notebook](https://www.kaggle.com/code/asteriospanoras/qcc-28-12769-what-works-what-fails-and-why) |
| Formal checker of ANGLE_k on all measurement branches | manwithacat | | Reachable and strict readings | [notebook](https://www.kaggle.com/code/manwithacat/prove-angle-k-for-every-measurement-outcome) |
| Negative-results dataset | manwithacat | | NR ledger | [dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-chemistry-negative-results) |

Rule status: the host (Tanuj Khattar, topic 753092) ruled that ANGLE_k only has to be correct on the
reachable states (z = batch k-1 forward, z = batch k in reverse). Topic 757522 (9-15 Nov) asks again
for a written ruling for prize eligibility; no new host answer in our snapshot. Every frontier loader,
ours included, relies on the reachable reading (helper bits carried in word top bits).

## 4. Known nulls (do not repeat without a new angle)

| Null | Number that killed it | Source |
|---|---|---|
| lambda_eff below about 46.79 | Hard bound from the fixed shift and the gap floor; we are at 46.907 (0.25% above) | 753633 (Lynn Langit); re-derived in 2.6 |
| Further fixed-word weight refits | Local optimum on the active surface: 1e-6 to 1e-3% per step; an augmented-Lagrangian refit is worse | 753185 (manwithacat), Asterios Panoras notebook, Ahmed Berat Ozer v23 |
| brot 7 / 6 / 4 | error_2 349.7 / 794.5 mHa naive; QAT on older shapes plateaus at 30-40 mHa; 7-bit local search plateaus at 110 mHa; Syed Asad Ali: "7-bit angles did not work" | 753185 (alex chilton), 753687 (Wan Xiaolong), 755495 |
| Re-rounding 10-bit factorizations to 8 bits | error_2 162-201 mHa | 753185 |
| GPTQ-style rounding | 47.04 vs 46.86 mHa for a plain +-1 polish; the schedule matters, not the rule | 753185 |
| Linear change of the 11-bit address basis before the lookup | Needs 169 qubits (cap 128); streamed versions cost 350-974 Toffolis per circuit | 753185 |
| Reordering the Givens tree so batch XORs get simpler | Mixed ANF matrices have full rank, like random ones; 18 other trees round worse than the starter tree (219 vs 236+ mHa) | 753185 (alex chilton) |
| Orbital relabelling | Starter numbering worth about 20% in rounding error; spectral orders only tie the DMRG file order | 753185 |
| R 10 -> 9 merges, (9,143,42), (8,159,42) cold starts | 22.1, 25.8, 29.0 mHa continuous | 753185 (marculera, manwithacat) |
| Dropping one copy of (9,143,40) -> (9,143,39) | 14.2 mHa continuous, rounding not finished (inconclusive) | 753185 (rikirinjani) |
| R = 1 big dictionaries (1,460,384), (1,496,300) | 20.5-22.4 mHa continuous | 753185 |
| Eckart-Young floor | RC >= 290 for 20 mHa (we measured L = 290 -> 19.9 mHa, 352 -> 7.4 mHa) | 753185 + section 2.5 |
| Shared eigenbasis / joint diagonalization of the factors | Barely better than random | 753185 |
| Symmetry shifts (multiples of N, N^2) and orbital rotations | Overwritten by the verifier's fixed shift; reference integrals fixed in their basis | 754605 |
| Small-B regime (B < 76) that the cost model prices at 0.4x-0.7x | Leading pair-space eigenvectors are full rank: rank-r leaf truncation errors 10.8 Ha (r = 4) to 2.5 Ha (r = 64) | section 2.5 (ours) |
| B <= 128 or brot >= 9 | QROAM k2 jumps to 5: Q +270 (+64% cost) | section 2.4 (ours) |
| Other J / lowbits on (10,129,60) words | J1, J2, J4, J5 worse than J3 there (Toffolis); two-phase cache at J5 costs 508 vs 401 break-even | 753185 |
| AlphaTensor-Quantum on loaders | Trades Toffolis for T/CS gates, which the loaders cannot use; 130M local moves found nothing | 754605 |
| Narrower single loader | 126/127-wire loaders do not change Q_total (the temp is set by the inner prep) | Asterios Panoras notebook |
| LLM-evolved gauge optimizer | 0 improvements in 15 generations | 753185 |


Nulls on paper from the family study of 19 Nov (S2.23):

| # | Idea | Measured result | Source |
|---|---|---|---|
| S2.23-1 | Row-axis shape moves at the linear row price | Σ core rank invariant under row removal; only caches scale (1.9 T per row): R 7 −2.0%, R 6 −3.5%, B 129 −1.6% (not −7.5 / −13 / −5.2%) | `chemistry_families_2026-11-19.md` §0.3–0.4 |
| S2.23-2 | Force the inner QROAM block to 2^3 (Q 427 → 276, −22%) | the model keeps the largest k with temp < 75·brot: needs brot ≤ 3 or R·C·B < 2,432 | same, F2 |
| S2.23-3 | CHEM-B20 block-localized dictionary | 76 one-body rows keep ≥ 5 generic blocks per plane; NR-0249 truncation 1.2–2.3 Ha at k ≤ 8 | same, F4 |
| S2.23-4 | CHEM-B22 dithered twin vectors | B 159 costs +4.2% for 128 twins; the freed C budget is worth −1.9% (−3.6% at best) | same, F5 |
| S2.23-5 | CHEM-N11 codebook-index loaders | index saves ≈ 27 products per circuit; per-column decoders cost ≈ 350–600 ANDs per circuit | same, §3 |
| S2.23-6 | Classical copy of the address | no Z-basis measurement in kickmix (only HMR and R) | same, §3 |
| S2.23-7 | Random-condition Toffoli halving | no gadget computes an AND with expected < 1 counted Toffoli | same, §3 |
| S2.23-8 | CHEM-B03 address-permutation gauge alone | ≤ about 110 rank units (−0.6%) by counting | same, F6 |


Corrections and public nulls of 21 Nov (the corrected rows above stay, read with this domain):

| # | Idea | Measured result | Source |
|---|---|---|---|
| CORR-3 | Correction of "Q pinned at 427" (S2.23 F2) and of the k2 cliff "B ≤ 128 → Q +64%" | true only for B in [129, 256] at brot 8: at (9, 127, 46), brot 7, the inner temporary is 287 and Q 405 | `public_frontier_2026-11-21.md` §3 |
| CORR-4 | Correction of NR-0124 / CHEM-B11 (brot 7 accuracy null, 27.1 mHa) | true for roundings of 8-bit-optimal fits; a staged rounding with refit at (9, 127, 46) passes at e1 19.9959 / e2 19.9992 mHa | `public_frontier_2026-11-21.md` §3 |
| P21-12 | C 46 → 45 at 7 bits, fixed words | e2 21.6 mHa, λ > 48.2; discrete plateau 20.351 mHa | goodpjw2008, Almond Labs |
| P21-13 | Refit of (W, H2_I) on the 13,290 words | no measurable gain (three active limits) | goodpjw2008 |
| P21-14 | Loader exchange across 41 public families | −38 T is the limit (13,328 → 13,290) | goodpjw2008 |
| P21-15 | Loader Toffolis writing only dead words; Toffolis replaceable by CNOTs | none of 1,856 | goodpjw2008 |
| P21-16 | Cheaper public loaders at positions 6 and 14 that change live words | e2 58–1,448 mHa | goodpjw2008 |
| P21-17 | 257-T k0 by a fake-carrier dead-cell rewrite | invalid on the reverse (involution) pass | Almond Labs |

## 5. Ideas

Conventions. Gains are task-cost reductions against our base (296,635,968). Multiply by about 1/3 for the
score. `Priority = P x (mid gain %) / effort hours`. Effort counts our hours. GPU notes: float64 fits and
refits run at about 11 steps/s on a T4 vs about 400 on an A100 (753185); use a Kaggle **P100** (float64
about 4.7 TFLOPS) for anything that trains or refits a table. Every row that changes words or weights
must end at error_1, error_2 <= 20 mHa, Lambda > 49 and E_gap > 3.914898 (strict). Every row that
changes a loader must pass the all-branch proof (`verification/prove_loader_branches.py`) and several salts.

### A. Improve our code (loaders, weights, ports of public improvements)

| ID | Idea | Mechanism | Est. gain % | P | Effort h | Priority | Risks / interactions | Sources |
|---|---|---|---|---:|---:|---:|---|---|
| CHEM-A01 | Adopt the public frontier parts: Lera's coded table + maoxin123's 4,022 loaders + Atakan Aldemir's weights `703843d3` | -52 Toffolis/pass (T 14,810 -> 14,706), lambda_eff ~46.92 (+0.03%) | 0.67 (0.6-0.7) | 0.95 | 2 | **0.318** | Loaders are CC0 compiled-only, no generator: prize-eligibility risk (PLAN 12). Record SHA-256s, authors, licences. Supersedes A04/A09 on our table | 759220; maoxin123 notebook; `downloads/datasets/atakanaldemir__...`, `marculera__dfthc-c44-build/NOTICE` |
| CHEM-A02 | Daily "parts market" watch: pull new public chemistry loaders and weights, re-verify, adopt the cheapest compatible pair | Weights and loaders on the same word table are interchangeable; the frontier moved 0.1-1% per day in November | 0.3 (0.1-0.7) | 0.7 | 4 | **0.053** | Only same-word pairs combine. Each new part needs licence + SHA record. Pair with A03 for eligibility | 753923, 756210, 759220 (RezOps), manwithacat circuits dataset |
| CHEM-A03 | Rebuild Lera's table loaders with Syed Asad Ali's compiler (our source) | Same table as the frontier, but a from-source build; adapt the planner to the XOR-of-five codes and the 4-form gauge | 0.6 (0.5-0.7) | 0.5 | 16 | **0.019** | Syed's emitter asserts his own helper convention (`(w>>7)^(w>>6)` table equal across batches): needs rework. Gives a prize-eligible version of A01 | v12 `compiler/emit_plan.py`; Lera notebook section 5 |
| CHEM-A04 | Code the one-body rows harder | error_1 weights each one-body row by its eigenvalue: rows with small abs(eigenvalue) tolerate coarse or coded words. Add bit-0/bit-1 codes or constant words on those of the 76 rows, refit | 0.3 (0.1-0.6) | 0.35 | 8 | **0.013** | error_1 is at 19.61 (frontier tables 19.9975): little slack; the one-body rows also enter E_gap via h1' eigenvalues | `dfthc_math.reconstruct_h1_thc`; 753185 (mrmeows one-body refresh) |
| CHEM-A05 | Port the codes (bit-0 affine, XOR-of-five at bits 1-2) onto our own table | Re-round and refit our words under the code constraints, as maoxin123 did on Lera's | 0.75 (0.5-1.0) | 0.4 | 24 (P100) | **0.013** | Redundant if A01/A03 switch tables. Spends error budget at the 20 mHa edge | maoxin123 v25 notes; Lera notebook section 3 |
| CHEM-A06 | Use the spare qubits in loaders k1..k14 (126 of 128) | Two more cache wires: cache one more monomial / keep one product live, cutting 1-3 products per loader | 0.15 (0.05-0.3) | 0.5 | 6 | **0.013** | Free in Q_total (angles temp 75-77 < inner temp 303). Must keep the five-dirty-words layout | `ops.py` count: k0 128 qubits, others 126 |
| CHEM-A07 | Per-slot best-of between maoxin123's and Lera's loaders on the same words | Both families load Lera's words; take the cheaper loader per slot | 0.05 (0-0.1) | 0.5 | 2 | **0.013** | Helper-bit conventions must match across adjacent slots (reachable states!); verify end to end | maoxin123 4,022 vs Lera 4,025 |
| CHEM-A08 | Search the low/high address split and low-block partition per loader | Rerun Syed's planner (lowbits 4 or 5, partial low cache, five-block plans) per slot with more seeds | 0.3 (0.1-0.5) | 0.4 | 12 | **0.010** | Already tuned once (v9 202 -> 195 products) | v12 `plans/`, `research/scripts/planner.py`, 756558 |
| CHEM-A09 | Bring loaders k0..k3 down to the k4..k14 level | k0-k3 cost 287-293 vs 263-268 (about 92 Toffolis/pass extra): extend the codes to batches 0-2 and give k0 its own low-cache window | 0.65 (0.3-1.2) | 0.3 | 20 | **0.010** | Batches 0-2 hold the root columns of the tree (most sensitive); may cost error budget. Can reorder edges (any topological order) to move sensitive columns | maoxin123 (k0 292, k1-3 286-287); v12 `k0_windows.json` |
| CHEM-A10 | Park more helper monomials in free word bits | The top-bit gauge parks 5 monomials per batch; Lera's negation gauge frees every column whose bits 1-2 are free. Use them for 1-3 more cached monomials per loader | 0.25 (0.1-0.4) | 0.3 | 10 | **0.0075** | Relies on the reachable-state ruling (753092 vs open 757522) | Lera notebook finding 1; maoxin123 v25 |
| CHEM-A11 | Peephole / kernel-relation search on compiled loaders | Find products that are XOR-combinations of other products already computed (linear relations over GF(2) in the product span) and drop them | 0.2 (0.05-0.4) | 0.3 | 10 | **0.006** | Needs our own formal check after each edit | qwertykeypadapple (-110 Toffolis on an older family, 756210) |
| CHEM-A12 | Exact re-planning toward the counting floor | ILP/SAT on the product cover per loader, allowing mixed products (both factors dense combinations), not only low x high pairs | 1.5 (0.5-5) | 0.15 | 40 | **0.0056** | Asterios estimates the floor at about 25 products per loader below today; nobody has reached it | Asterios Panoras notebook (sections 5.3, 8) |
| CHEM-A13 | Global (multi-start) weight refit at fixed words | Basin hopping across sign patterns of w to close the triangle slack in lambda_h2 | 0.03 (0.01-0.1) | 0.3 | 6 (P100) | 0.0015 | Hard cap 0.25% (section 2.6); all gates active | 753633, Asterios Panoras ("local optimum") |

### B. Rewrite with a different technique (shape, loader architecture, regime, retraining)

| ID | Idea | Mechanism | Est. gain % | P | Effort h | Priority | Risks / interactions | Sources |
|---|---|---|---|---:|---:|---:|---|---|
| CHEM-B01 | Mixed-precision coded words, trained with a loader-cost surrogate | Give each (row, column, bit) a code depth according to its sensitivity (row weight L1 ranges 0.01-0.29; columns near the root matter more); train with a differentiable proxy of the loader product count (GF(2) rank proxy, Spearman 0.83 per Wan Xiaolong) | 2.5 (1-4) | 0.25 | 50 (P100) | **0.0125** | Spends the 20 mHa budget; must be co-tuned with lambda_eff. Biggest controllable lever (loaders are 55% of T) | 753687; maoxin123/Lera codes; our row-weight stats |
| CHEM-B02 | Row-sensitivity-weighted coding | Same as B01 restricted to rows: rows with small weights get extra coded bits | 0.8 (0.3-1.5) | 0.25 | 20 (P100) | **0.010** | Cheaper first step of B01 | our W statistics (no all-zero rows; 25% of rows have L1 < 0.19) |
| CHEM-B03 | Address-permutation gauge | Permuting b inside each rho (with W rows), permuting rho (with W, H2_I) and permuting one-body rows leaves every metric unchanged but changes the loader truth tables. Search permutations jointly with the codes, before final rounding | 1.0 (0.5-2) | 0.2 | 30 | **0.0067** | Breaks existing codes if applied after rounding (codes are per block of 16 addresses); not found in any public notebook or forum post | arXiv:2604.06578; scorer `split_u_vectors` / assignment of one-body rows |
| CHEM-B04 | New table at C = 43 | Cost model -0.53% at equal loaders; EY floor at RC = 344 about 8.4 mHa | 0.47 (0.4-0.53) | 0.4 | 30 (P100) | **0.0063** | Must keep the codes; loaders rebuilt (A03 pipeline). Syed went 45 -> 44 successfully | section 2.3; 756210 |
| CHEM-B05 | Exchange-rate-aware error budget | One optimizer that trades error budget between C, codes and lambda_eff at the exact rates (1 Toffoli/pass = 0.0135%, C-1 = 0.53%, lambda 1% = 1%) instead of filling the gates with lambda_eff | 1.0 (0.5-2) | 0.25 | 40 (P100) | **0.0063** | Infrastructure for B01, B04, B08 | section 2.3 |
| CHEM-B06 | Better initializers for lower-C tables | Seed from regularized compressed DF (RC-DF) / SCDF-style fits, or compact groups with orbital-locality priors (two-cubane structure), then fit at the scoring shape | 0.85 (0.6-1.1) | 0.2 | 30 (P100) | **0.0057** | Enabler for B08/B12; cold starts failed at C = 42 (753185) | arXiv:2212.07957, 2403.03502; Lera compact groups |
| CHEM-B07 | Fewer basis vectors per group: B = 139 (R 8, C 44) | -32 rows; k2 stays 4 down to B = 129; each B-1 is about -0.37% | 1.5 (1.0-1.9) | 0.15 | 40 (P100) | **0.0056** | Capacity: rows are the dictionary; B = 135 would be -3% (P ~0.07). B must stay >= 129 and B+1 should avoid crossing multiples of 16 | section 2.3, 2.4 |
| CHEM-B08 | New table at C = 42 | -1.06% from the cost model | 0.85 (0.6-1.1) | 0.25 | 40 (P100) | **0.0053** | Cold starts at (9,143,42), (8,159,42) failed (25.8-29 mHa) but warm starts from a coded C44 table were not tried | 753185 |
| CHEM-B09 | Co-design the Givens tree, edge order and gauge with the codes | Tree and edge order fix which columns go to which batch and which pi/negation gauges are free; choose them to maximize code-friendly columns per batch | 1.0 (0.5-2) | 0.2 | 40 | **0.005** | alex chilton found the starter tree best for rounding (19 trees); the trade is now rounding vs loader cost | PLAN 11 backlog; 753185; Lera finding 1 |
| CHEM-B10 | R = 7, B = 143, C = 48-50 | -143 rows: -6.6 to -7.5% | 7 (6-7.5) | 0.04 | 60 (P100) | **0.0047** | Capacity: RB 1001; Lera's compact-group seed would need a new grouping | section 2.3 |
| CHEM-B11 | brot 7 with codes and QAT | Cost model + loaders -10% | 9 (8-10) | 0.03 | 60 (P100) | **0.0045** | Every public attempt failed (110-350 mHa). The angle grid theta = pi*w/2^(brot-1) is fixed, so the only way in is a much better fit, not a smarter grid | 753185, 753687, 755495 |
| CHEM-B12 | New table at C = 40 | -2.1% | 1.8 (1.5-2.1) | 0.1 | 50 (P100) | **0.0036** | EY floor about 12.5 mHa before dictionary and rounding losses | section 2.5; Asterios estimate -1.8% |
| CHEM-B13 | Loader architecture from XAG multiplicative-complexity synthesis | Synthesize the 40-output delta function (with 828 don't-care addresses) as an XAG with minimum AND count, then map ANDs to CCX with HMR uncompute | 2 (1-5) | 0.1 | 60 | **0.0033** | The 128-qubit cap; most public effort already went into this space | arXiv:1908.01609; npj QI 8, 7 |
| CHEM-B14 | Use the word registers as dirty workspace | z holds known words; dirty-borrowing (toggle twice) gives extra cache width without clean qubits | 0.5 (0-2) | 0.1 | 30 | 0.0017 | Doubles the work on borrowed lines; reachable-state rule | arXiv:1812.00954, 2605.20334 |
| CHEM-B15 | J = 4 with larger caches | Q -8 (-1.87%); 19 loaders; +8 workspace qubits per loader. Break-even +138 Toffolis/pass; our simple model says +300 | 0.4 (-1 to +1) | 0.1 | 30 | 0.0013 | alex chilton: other J worse on (10,129,60) | 753185; section 2.3 |
| CHEM-B16 | Direction-symmetric word features | Use word bits that are equal in batch k-1 and batch k (beyond the top-bit helper) as free features of the address | 1 (0.5-3) | 0.05 | 40 | 0.0013 | Few such bits exist without spending accuracy | 753092 ruling; v12 helper table |
| CHEM-B17 | J = 6 or 7 | Fewer cache rebuilds (13 or 11 loaders) vs +8/+16 qubits; workspace 69/61 under the cap | 0.2 (-0.5 to +0.5) | 0.15 | 30 | 0.0010 | Cache of about 75 wires does not fit at J = 7 | section 2.3 |
| CHEM-B18 | Multi-basin lambda refit on a new table | Search sign patterns of w on a new table | 0.1 (0.05-0.25) | 0.2 | 30 (P100) | 0.0007 | Hard cap 0.25% | 753633 |

## 6. Combination notes

- **Do A01 first, as a hedge.** It is 2 hours for -0.67% and puts us on the public frontier. Keep our
  from-source v12 as the prize-eligible fallback until A03 (from-source loaders on Lera's table) exists.
- **A01/A03 and A05/A09 are alternatives**, not additive: the first switches to Lera's coded table, the
  second codes our own table. Pick one table line. Lera's line is ahead and has the codes already.
- **Everything that changes words or weights spends the same 20 mHa budget**: B01, B02, B03 (if applied
  after rounding), A04, A05, A09, and the shape moves B04, B07, B08, B10, B12. The gains do not add
  naively; build B05 (the budget optimizer) before stacking more than two of them.
- **Shape changes reset the loaders.** B04/B07/B08/B10/B12 need new words, new codes, new loaders
  (A03 pipeline) and new weights. Each is a full rebuild; do them only with the pipeline automated.
- **Loader-only ideas stack** with any table: A06, A07, A08, A10, A11, A12, B13. They are cheap to test
  with the all-branch prover and the verifier (1 min in WSL).
- **lambda_eff ideas (A13, B18) are capped at 0.25%** in total; ignore until everything else is done.
- **Combined realistic target:** A01 (-0.67%) + A06/A08/A11 (-0.3%) + B01/B02 codes (-1 to -2%) + C 43
  (-0.5%) gives about -2.5 to -3.5% task cost (+0.8 to +1.2% score). Getting -10% or more on chemistry
  requires R 7, brot 7 or a new loader architecture (B10, B11, B13), all low probability.
- **Risk:** chemistry is checked on 64 sampled measurement records per address; a rare-branch fault
  would score the whole submission 0. Every new loader must pass the exact all-branch prover and at
  least 4 salts. Keep the reachable-state assumption visible: if the hosts ever rule for the strict
  reading (757522), every frontier loader with parked helpers becomes invalid at once.

