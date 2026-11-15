# What didn't work: 30-odd measured nulls across invsqrt, dqi and the chemistry

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753185 — copied 2026-11-15 (raw Markdown, unedited).

**alex chilton** · 2026-10-27 12:39 UTC · 5 votes

Most of what gets posted here is what worked.
Here is the other half: everything we tried that did not work, with the number that killed it.
If you were about to try one of these, this might save you a day.

All numbers are from the official verifiers or from our own exact re-implementations of their formulas, run locally.
The full table, with the source write-up for each row, is in the "What we tried that did not work" section of our notebook: https://www.kaggle.com/code/alexchilton/21-39-lb-verified-three-task-portfolio

## Chemistry (dfthcsos)

- **Re-rounding public 10-bit factorizations to 8 bits** (Özer's (32,60,30), (26,60,30), (24,60,30), marculera's (10,120,54), Gluzdov's word files): best error_2 162 to 201 mHa against the 20 mHa gate. The gate is set by how the words were trained, not by the weight count.
- **Fewer angle bits** (brot 7: error_2 349.7 mHa; brot 6: 794.5 mHa; brot 4: three orders over) and quantization-aware retraining on (40,60,30), (50,60,30), (55,60,30): every arm stuck at 30 to 40 mHa.
- **Lower Lambda with a different fit**: every fit anyone has measured sits at Lambda 290 to 308. A 2x would need 75.8.
- **Other J / lowbits on Rustad's (10,129,60) words**: J1, J2, J4 and J5 at lowbits 4 or 5 are all worse than J3. Our two-phase cache extended to J5 costs 508 Toffolis per circuit against a break-even of 401.
- **A linear change of the 11-bit address basis before the lookup**: the best layout reaches 266 products per circuit but needs 169 qubits (cap 128); every map whose cache fits is above 313. Streaming the high monomials over 2 to 6 measurement phases fits the cap but costs 350 to 974 Toffolis per circuit.
- **Rounding for the circuit, or reordering the Givens tree so consecutive word batches XOR to simpler tables**: on Özer's v9, the mixed (low x high) ANF matrix of every one of the 525 output bits has full GF(2) rank 15, exactly like a random matrix. His 321 Toffolis per circuit sit on the 262-product rank bound plus the caches. Saving one Toffoli would need an output's 85 columns to collapse into a 14-dimensional subspace. Column-pair XORs are 0.95 to 1.0x as complex as single columns, so reordering buys nothing either.
- **GPTQ-style rounding** (the LLM-quantization trick: round one angle at a time and correct the rest of the row in closed form from the Hessian): it does lower the error of each rounding step, but the refit after each stage absorbs it. From the same continuous start and schedule it ended at 47.04 mHa against 46.86 for a plain +-1 polish. What decides the result is the schedule (how finely you freeze and refit), not the per-row rounding rule.
- **One group fewer (R 10 -> 9) on the current best chemistry**: Özer's v9 is at error_2 19.998 mHa and E_gap 3.9151 (floor 3.915), and marculera measured nine merged groups at 22.1 mHa after rounding.
- **LLM-evolving the gauge optimiser** with ShinkaEvolve and a local Qwen 30B: 0 improvements in 15 generations. (Gemini in a hand-driven loop improved on its first try.)
- **Our gauge on Özer's v9 chemistry**: -0.0004 %. His lambda loop leaves nothing to rotate away.

## DQI

- **888 qubits from three source edits**: passes, but 4.7590e11 against 4.4860e11 at 889.
- **Dropping the jet zero flag**: 2 extra value errors on identical inputs, R 1.19126 -> 1.20139 (+0.85 %), against a qubit worth 0.225 %.
- **Freeing the constant-1 border word in the reverse folds**: 30 value errors, R +7.1 %.
- **Narrower reverse division**: +1,518 Toffolis (+0.32 %), more than the qubit is worth.
- **Porting the BM-decoder gadgets (Kozykorpesh, mrmeows q870) to the Schur decoder**: the structures they rely on do not exist there.
- **A cheaper inverse for the reverse extraction**: reverse costs 223 Toffolis per pair against 192 forward, and every route we found costs the same.
- **Peak shaving on Özer's 877-qubit build** (reorder, recompute, dirty or clean borrowing): with cost ~ Q^2 T, one qubit is worth about 1,100 Toffolis. The 877 peak is hit in 2,046 windows across about 684 kernel calls, there are 0 known-value qubits to borrow at the 1,734 jet windows, and the cheapest narrowing costs at least 23 to 42 Toffolis per window. So: nothing under the exchange rate.

## Inverse square root

- **Borrowing zero or affine-redundant wires, dropping degenerate Toffolis** (on our 66-qubit build and on the public 63): 0 degenerate controls in about 1,500 Toffolis, no borrowable wire in the rows.
- **Retuning or dropping the 44 rotation angles**: the LP gives the identical margin; dropping 1, 2 or 3 rotations gives worst error/tolerance 1.354, 2.214, 3.937.
- **Dropping one of the 32 Horner table bits and refitting**: every variant fails, the best at 1.228x the tolerance.
- **A dirty-wire AND for a row-opening gate**: the cheap cover (+23 Toffolis) is not valid, because the borrowed wires are read up to 8,474 times before the replay closes.
- **63 -> 62 qubits at all**: manwithacat's "Is the 63-qubit inverse square root optimal?" closes it far better than we did. Worth reading.

## Literature that did not give anything here

QROM halving (the chemistry QROAM cost is the verifier's formula, not our circuit), a GF(2^12) multiplier below rank 42 (only lower bounds exist), the reference DQI design in arXiv 2510.10967 (5.72e6 Toffolis), ECC repos and the MBUC gadget catalogue.

Credit to everyone whose public work made these tests possible: mrmeows, Agniv Banerjee (somejoint), Ahmed Berat Özer, Michael Rustad, marculera, heliosli, manwithacat, Kirill Tymchishin, Dmitrii Gluzdov, Kozykorpesh, jday96314 and the rest.
If you have a null of your own, add it below. A list like this is most useful when it is long.


---

## Comments (4 threads)

Replies are indented with `>`, one level per reply depth.

**rikirinjani** · 2026-11-04 05:56 UTC · 1 votes

Thanks for the list. A few measurements of ours to add, all from the official verifiers or our own exact scripts on the box, done with AI assistance. Some are not nulls, just numbers that may save someone a run.

Chemistry (dfthcsos)

Dropping the least harmful copy of marculera's published (9,143,40) Adam table and refitting: L-BFGS 20,000 iterations gives error_2 14.464 mHa at (9,143,39), then 200,000 Adam steps (lr 1e-4 to 1e-5, kappa 5e-5) give 14.246 mHa, continuous, before 8-bit rounding. We stopped before the rounding finished, so we cannot say whether it passes the gate.
Cost of that run: on a Kaggle T4 the float64 Adam stage runs at about 11 steps/s (200k steps = 5.25 h), against roughly 400 steps/s reported on an A100. A T4's float64 rate is the limit; the second GPU of a T4x2 session is not used by this code.
DQI

The score's repetition factor has a floor: with 0 errors in the default 102,400 shots, the Clopper-Pearson bound at 2^-128 gives eps about 8.66e-4 and R = 1.1913. Any R below that needs more shots than the verifier draws, so only errors can move it, upward.

Salt check on one public 891-qubit decoder (the one in several 24.50 portfolios): default, s1 and s2 give 0 errors (cost 3.7538e11); salt s3 gives 1 value error, R 1.1913 to 1.1967, cost +0.46%. The version with one CCZ fewer (24.53 class) was clean on all four salts. Worth a few extra --salt runs before trusting a dqi file.

Static phase profile of that 891-qubit decoder (397,194 Toffoli-like gates): unextract 97.5k (24.5%), unconvert 55.8k, convert 54.3k, extract 53.4k, scan 44.3k (only 23 ancillas), unsynthesis 38.5k, synthesis 31.5k. All 51 ancillas are active over 94 to 95% of the circuit. That fits the peak-shaving null above.

Pointers to two closed DQI levers that dcvxzyk published as notebooks (credit to them): a sound static ancilla remap finds all 52 ancillas co-live from extract through unextract (best static Q' = 892), and caching the reverse inversions costs more than it saves for 12 to 240 cache qubits.

Data

The FeMoCo FCIDUMP in the zhendongli2008/Active-space-model-for-FeMoco repository (FCIDUMP.zip) is byte-identical to the competition's FeMoCo76o113e.fcidump (SHA-256 e5a48af7...e1bd, 128,128,348 bytes).
Credit to everyone whose public work these runs built on.

> **hengck23** · 2026-11-04 06:09 UTC · 0 votes
>
> Go to the data repo and find the paper. Maybe the chemical structure can help in factorisation and compression

> > **rikirinjani** · 2026-11-04 07:52 UTC · 0 votes
> >
> > Thanks hengck23, we followed this up. The repo's FCIDUMP is the competition file (byte-identical), and its README says the 76 orbitals are split-localized MOs ordered by the genetic algorithm in BLOCK, so the file order is already a structure-aware (DMRG) order.
> >
> > We tested whether re-ordering the orbitals could do better, on the two-body integrals (fraction of sum of squares with both pairs (p,q) and (r,s) within |index difference| <= 8): random permutations 0.299 +- 0.016 (300 draws), file order 0.612, weighted Fiedler order of |K_pq| 0.612 to 0.613, reverse Cuthill-McKee on the top 10% of |K_pq| edges 0.620 (lower than the file order at band 16). So spectral orderings only tie the file order; nothing near a 3% gain. The file also has ORBSYM = 1 for every orbital, so there is no exact point-group symmetry to use.
> >
> > Not tested: block seeding of the factors from the two-cubane locality, or the atom-class labels in the paper (the per-orbital class assignment is not in the files we could read). If anyone has the MO coefficients with atom labels, that would be the next thing to try.

**yoyobar** · 2026-10-30 04:03 UTC · 2 votes

yoo~~ gpt-6 Astra has just developed a Inverse square root circuit that only needs 62 qubits, with cost = 101,002.54.

https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits

Amazing!

![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F8293010%2Fffe12eaab1da914fb5c13528a6620f0d%2F692be7ae-6e81-47a5-9a04-1dd4b6563b7b.png?generation=1793333004714390&alt=media)

> **(deleted user)** · 2026-10-30 05:33 UTC · 0 votes
>
> *(deleted)*

> **alex chilton** — **Topic author** · 2026-10-30 07:06 UTC · 1 votes
>
> yay girl! :-) excellent :-)

> **(deleted user)** · 2026-11-03 12:43 UTC · -1 votes
>
> *(deleted)*

**Mahmoud saad Elgharib** · 2026-11-04 14:51 UTC · 0 votes

Getting the cost down to around 100k with 62 qubits is a really smart approach to the circuit compile.
😎 🆒️

**manwithacat** · 2026-10-27 19:56 UTC · 0 votes

Great list, thank you. A few more from our side, in no particular order. Numbers are from the local official verifiers or exact local analysis unless stated, and the work was done with AI assistance.

**Chemistry** (details and the tools to check them in our [chemistry negative results notebook](https://www.kaggle.com/code/manwithacat/chemistry-negative-results-what-fails-and-why))
- Rounding a fixed 12-bit factorization uniformly to 11-bit angles: verifier FAIL at 76.9 / 29.7 mHa. Precision has to be co-fitted with the weights.
- Gauge and offset descent on public factors: H2_I gradient steps worsen cost in both directions (cusps at zero eigenvalues and zero weights); active-set projection 1.2e-5; pruning the weakest copy +0.29%; a fixed-word λ_eff refit −0.08%, decaying to about −0.13%. Consistent with mrmeows' ~0.5% bound.
- A smaller dictionary at R = 1 with 10–11-bit words, (1,460,384) and (1,496,300): 20.5–22.4 mHa continuous against the ~17 needed, with Gauss–Newton gaining ~1e-4 mHa per iteration after 3,764 iterations. Capacity, not the optimiser.
- Eckart–Young floors from the pair-space spectrum: R·C ≥ 290 for 20 mHa (307 for 15, 330 for 10). No shape below that passes however it's trained.
- A shared eigenbasis for the factors: pairwise commutators barely better than random (0.148 vs 0.161), joint diagonalisation moves off-diagonal mass only 17.6% → 16.4%, and a greedy eigen-dictionary sits at 1,279 mHa at B = 496.
- Dictionary-only search with an exact inner solve for the weights: public weights are within 0.5 mHa of optimal for their dictionary (19.466 vs 19.990), and dictionary-only L-BFGS loses to joint Gauss–Newton at equal wall clock (17.585 vs 17.072 mHa).

**DQI**
- Dropping the gap flag from the Berlekamp–Massey state: even restricted to syndromes of correctable errors, 24 post-step states collide between the two flag branches, so the flag carries information.
- Reusing the half-degree "Dialog" division in the reciprocal step: it divides modulo the wrong polynomial (E, where the map needs O² modulo σ = E² + zO²), so it can't be dropped in.
- A phase-only sweep with Newton/Hensel inversion: doesn't converge (classical test). The Fermat/Itoh–Tsujii replacement costs about 1.16 M Toffolis one way (formula, not built).
- Hunting cheap GF multiplications in a public 870-qubit BM decoder: 72.3% of its Toffolis are 16,998 GF multiplies (14,750 full 12-bit at 42 each), and there are no squares or constant/classical-factor products left: the designers already route those through CNOT-only operations.

**Inverse square root**
- Trading rotation precision where the error has slack: mapping error/tolerance over the whole input domain (with the verifier's own chunking and seeding) found no region with enough slack. The circuit is precision-saturated everywhere, so the rotation price can't move.

**Tooling, for anyone building local loops**
- Sampled "screening" doesn't beat the official verifiers on broken circuits. They already abort at the first failing chunk (a broken invsqrt candidate is refused in ~1.3 s), and DQI's ~1 s of load and hash is a floor.
- Simulating linear (X/CX/SWAP) blocks as GF(2) matrix products on tensor cores: only ~5% of DQI's ops, and none of invsqrt's or chemistry's, sit in blocks long enough to pay off.
- A straightforward bit-sliced GPU kernel (Metal, exact against kickmix): about the throughput of all 10 CPU cores of an M2 Pro (~5e11 op·shots/s), so no win on a laptop.
- Running verifier batches on threads with the GIL released: 2.7× for DQI but only ~1.2× for invsqrt, where per-shot Python work dominates. Separate processes are better.

> **manwithacat** · 2026-10-28 08:51 UTC · 0 votes
>
> **Chemistry follow-up: block-aligned shapes at lower C.** We tried to make a (10,129,46) shape from a public portfolio cheaper by moving to the block-aligned (9,143,42) and (8,159,42), which the cost model prices about 5–6% lower. We trained both from random initialisation on a free Kaggle T4: float32 Adam for up to 295,000 steps, then an hour of float64 L-BFGS on the official error₂ residual. Best error₂ was 25.8 mHa for (9,143,42) and 29.0 mHa for (8,159,42), before the 8-bit word penalty (roughly another 4 mHa by our calibration), against the 20 mHa limit. Short warm starts from the public factor had already reached about 22.7 and 29.4 mHa, so the long cold start didn't help. At these shapes, capacity runs out below C = 46 before the cost saving pays. One seed and one optimiser recipe per shape, so this is a negative result for this recipe rather than a proof. Local analysis only; nothing went through the hosted verifier.

> > **alex chilton** — **Topic author** · 2026-10-28 11:37 UTC · 0 votes
> >
> > thanks, that's a really useful set. one more from us: the givens tree. we re-encoded the same continuous C46 vectors exactly on 18 other trees and rounded to 8 bits - the starter tree rounds best of all of them (219 mHa vs 236 for a binary heap, relabelled orbitals +17-29%, chain or star +80%). the starter's own orbital numbering is worth
> >  about 20%. a staged round-and-refit on the best rivals is running now, will post if it changes anything.
> >
> > the eckart-young floor is the interesting one for us - R·C >= 290 against 460-496 in the public shapes means there is room on paper, the question is just whether anything can train down there.
