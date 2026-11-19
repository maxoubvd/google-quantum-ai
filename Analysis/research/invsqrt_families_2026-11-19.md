# invsqrt: divergent family brainstorm, priced on paper (S2.18, GEN-B11)

Written 2026-11-19 for plan step S2.18 (milestone M2.B, idea GEN-B11 applied to invsqrt). Inputs: the official
verifier (`verifiers/invsqrt/verify_submission_invsqrt.py`, the specification), `Analysis/research/invsqrt.md`
(anatomy §2, nulls §4, ideas §5), `experiments/NULLS.md` (invsqrt, 57 community rows + 28 of ours),
`Analysis/research/public_frontier_asterios_2026-11-19.md`, `Analysis/research/public_frontier_others_2026-11-19.md`,
the team's exact integer model `generators/invsqrt/team/table_search.py`, and Syed Asad Ali's divider notes in the
v12 tree (`src/work/c10/narrow/arithmetic_round2/research/divider/`, AR004 result). Numbers measured by others are
credited; numbers marked *estimate* are ours, on paper. The two new claims (sections 5.1 and 5.2) were checked
exhaustively by a short script, reproduced in Appendix A.

## 1. Summary

1. **Budget.** The public frontier `ecf908fe` costs 44 Q × 1,690.5 = 74,383.7 (1,450.9 mean Toffolis + 33 rotations
   at 7.23 Toffoli-equivalents each + T gates). A family must reach `C < 74,384 / Q`: 1,771 at 42 Q, 1,617 at 46 Q,
   1,488 at 50 Q. Two floors bind every family: the 22 x qubits, and about **33 generic rotations** (the 25
   significant bits of the result slide over absolute weights 2^-4 … 2^-36 as the exponent k runs over 0..10;
   ≈ 239 Toffoli-eq at the frontier margin, ≥ 221 even with zero arithmetic error).
2. **Seven families were priced** (section 3). **Six are nulls on paper**, each above the frontier by +9% to more
   than +100%: divider-free polynomials (F2, width ≥ 50 Q), absolute-frame streaming quotients (F3), table-free
   reciprocal-root recurrences (F4, ≥ 74 Q), a phase-additive split with an exponent-indexed phase table (F5,
   +13 rotations), a logarithmic number system (F6) and direct Boolean synthesis of the phase bits (F7). They are
   worth recording because each null has a structural reason, not a tuning one.
3. **One family survives: F1, the history-free ("parity") divider.** A new division row — subtract the odd divisor E
   when the remainder is non-negative, add E + 1 otherwise — costs exactly one masked adder like a non-restoring row,
   but its quotient digit equals the least significant bit of the next remainder. Rows that take no dividend bit
   therefore need **no stored history**: they reverse with zero extra Toffolis, and their digits can be measured out
   (HMR) right after the phase stage, the phase debt being repaid for free (a classically conditioned Z on the
   remainder's low bit) during the reverse pass. Exhaustively checked for every odd divisor 3..199 (Appendix A).
   In the incumbent, 12 of the 24 divider rows (rows 24–35) take no dividend bit; freeing their 12 history wires
   during the divider clean-up should remove most of its width premium (clean-up 308–323 vs forward 270).
   *Estimate:* −35 to −60 Toffolis, **C ≈ 1,640, cost ≈ 72,200 (−2.1% to −3.6%)** at the same 44 Q.
4. **A companion impossibility result** (section 5.2): with a one-bit dividend feed per row, no division row whose
   two offsets have different parities can keep its remainder in a fixed interval. So the dividend's information
   (about 10–15 bits of h) must be stored as quotient history; only feed-free rows can be history-free. This
   explains why the jeffemerson/NR-0261 nulls hold for the feed rows, and bounds what F1 can gain.
5. **What this says about the leaders.** No family priced here reaches the −30% class on invsqrt. The best
   plausible invsqrt target from our own source is F1 on a 44-qubit base plus the in-family width moves
   (INV-A21–A24): −3% to −9%. Unless a family outside this list exists (INV-B01's search could find one), the
   leaders' 1.42× cost-product advantage is more likely on dqi and chemistry.
6. **Recommendation for S2.19–S2.20:** prototype **F1 only** (section 6), in two steps: (a) put the parity tail into
   the exact integer model and measure the phase margin on all 3,133,452 inputs (2–4 h); (b) a scaled kickmix
   divider with an 8-bit divisor, verified exhaustively, measuring the clean-up Toffolis at a width cap with and
   without the parity tail (6–10 h). The second recommendation is not a new family: build the scaled-proxy harness
   (GEN-A13) on which F1 and our own 44-qubit generator (M3.2) both depend. No second family deserves prototype
   time on this evidence.

## 2. The budget every family must beat

### 2.1 Cost model and exchange rates (from the verifier)

`cost = Q × (Toffolis + T/4)`. Each generic `Z_POW` costs `0.56·log2(1/eps) + 5.3` T with
`eps = π·M/2^64/n_rot`, `M` the smallest `tol(x) − |error(x)|` over all x (units of 2^-64 turn), `n_rot` the mean
number of executed generic rotations. At the frontier (M = 38,858,848 ≈ 2^25.21, n_rot = 33):
log2(1/eps) = 64 − 1.65 − 25.21 + 5.04 = 42.18, so 28.92 T = **7.23 Toffoli-eq per rotation**, 238.6 for 33.

| Lever | Value at the frontier (44 Q, C = 1,690.5) |
|---|---|
| 1 Toffoli | 44 cost units = 0.059% of task cost |
| 1 qubit | 1,690.5 cost units = 2.27% ≈ 39.3 Toffolis (at 44 → 43) |
| 1 generic rotation more or fewer | ≈ 7.2 Toffoli-eq (plus 0.56/4·log2 of the n_rot change on all of them) |
| Doubling the margin M | −0.14 Toffoli-eq per rotation, −4.6 for 33 |
| Zero arithmetic error (M = min tol ≈ 2^29.2) | 6.70 per rotation, 221 for 33 (−17.5): the ceiling of any precision gain |

Break-even charge `C* = 74,383.7 / Q`:

| Q | 40 | 42 | 43 | 44 | 46 | 48 | 50 | 53 | 62 | 74 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| C* | 1,859.6 | 1,771.0 | 1,729.8 | 1,690.5 | 1,617.0 | 1,549.7 | 1,487.7 | 1,403.5 | 1,199.7 | 1,005.2 |

### 2.2 The incumbent, stage by stage

Measured on our v12 build (46 Q, `invsqrt.md` §2.2) and approximated for the 44 Q frontier (*estimate*: v12 stages
scaled with the public deltas: 64-row table −47, 46 → 44 moves, micro-edits −0.8%).

| Stage group | v12 46 Q (measured) | Frontier 44 Q (*estimate*) | Live registers at the peak |
|---|---:|---:|---|
| Normalise fwd + inv (shift by 2k, 4 flags) | 94.3 | ≈ 94 | x/n 22, flags 4, prefix ANDs |
| Consuming square root fwd + inv | 192.0 | ≈ 192 | n → v 12 + h 10, flags 4, carries |
| Table loads + body return | 133.3 | ≈ 75 | address a (6) + loader work |
| Coefficient, seed, product, fwd + inv | 389.0 | ≈ 380 | v 12 + flags 4 + C 13 + acc 16 + carry 1 = 46 → 44 with v10, v11, e2, e3 hidden |
| Divider fwd + exit + clean-up | 598.0 | ≈ 605 | endpoint: 43 owners (v, flags, remainder, quotient/history, after retirements) + carries |
| Phase ring (rotate the quotient by k) | 105.0 | ≈ 105 | quotient 27 + the endpoint registers |
| **Toffolis** | **1,511.6** | **1,450.9** (measured) | |
| Rotations (33) + T | 231.1 | 239.6 | |
| **C** | **1,742.7** | **1,690.5** | **Q = 46 / 44** |

Two facts from Syed Asad Ali's divider notes (v12 tree, AR004, Q48 build) matter below: the divider forward pass
costs **270** and its clean-up **308** (the +38 is the price of width), and retiring a single history wire after the
phase endpoint, so that the inverse rows get one more clean carry, saved **9.5** expected Toffolis.

### 2.3 Three floors that every family meets

1. **x is 22 qubits** and must be restored; all work is extra.
2. **Rotation floor ≈ 33.** The phase is `2^k·y` with y the normalised result (25 significant bits) and k ∈ 0..10.
   Every absolute weight from 2^-4 (small x) to 2^-36 (x near 3.1M, tolerance 2^-34.8) must be carried by some
   wire with a fixed angle, i.e. ≈ 33 generic rotations, ≈ 239 Toffoli-eq (≥ 221 with a perfect margin). NR-0174
   (the shift range spreads the significant bits over ~34 fixed-angle positions) is the same observation. Any
   family that phases two separately computed terms pays a second set of positions.
3. **Width at the phase endpoint.** The rotations need the whole 27-digit quotient at once (after the ring). A
   divider design then holds, raw, v 12 + flags 4 + remainder 14 + quotient 27 = 57 owners; the incumbent gets to
   43 by hiding v10/v11/e2/e3 and retiring history. A polynomial design instead holds its accumulator plus the
   operands still needed (section 4, F2). There is no information-theoretic reason for 44 (x has 21.6 bits and the
   25-bit result alone determines x), but every known route from x to y and back keeps two encodings alive.

## 3. Comparison

Gap is against the frontier 74,383.7. P is the probability that the family, built from our own source at full
scale, beats the frontier. Effort is focused hours to a full-scale build with exact proof (prototype hours in
brackets). All Q and C are *estimates* from section 4 unless marked.

| Family | Q | C | Q × C | Gap | P | Effort (h) |
|---|---:|---:|---:|---:|---:|---|
| **F1 History-free ("parity") divider** on the root + table + product datapath | 44 | ≈ 1,640 (1,630–1,655) | ≈ 72,200 (71,700–72,800) | **−2.9%** (−3.6% to −2.1%) | 0.30 (0.6 that it saves ≥ 2% on any divider build of ours) | 40–60 inside M3.2 (8–14) |
| F2 Divider-free piecewise polynomial (Horner/Estrin, on n or on the root coordinates) | ≥ 50 | ≥ 1,740 | ≥ 87,000 | ≥ +17% | 0.02 | 80 (20) |
| F3 Absolute-frame streaming quotient (ring-free, fixed-angle digits) | ≥ 46 (likely 50+) | ≈ 1,930 | ≥ 88,800 | ≥ +19% | 0.01 | 60 (15) |
| F4 Table-free reciprocal-root recurrences (Ercegovac–Lang, multiplicative normalisation / CORDIC, Newton, Goldschmidt) | ≥ 74 | ≥ 3,000 | ≥ 220,000 | ≥ +200% | < 0.01 | — |
| F5 Phase-additive split: exponent-indexed coarse phase table + shorter divider | 45–46 | ≈ 1,800 | ≈ 81,000–83,000 | +9% to +11% | 0.03 | 50 (12) |
| F6 Logarithmic number system (log2 n, halve, exponentiate) | ≥ 55 | ≥ 3,000 | ≥ 165,000 | ≥ +120% | < 0.01 | — |
| F7 Direct Boolean synthesis of the phase bits (ANF / Walsh, low degree) | ≥ 44 | ≥ 10,000 | ≥ 440,000 | ≥ +490% | ≈ 0 | — |

## 4. The families

### F1. History-free ("parity") divider

**1. Datapath.** Keep the incumbent's mathematics (normalise by 4^k, consuming root n → (v, h), coefficient
C(v) from a 64-row table, product P = U0(v) − C·h/2, divide by the odd divisor E = 2Q + 1, ring by k, rotate,
uncompute), but change the division recurrence in every row that takes no dividend bit:

- non-restoring row (incumbent): `σ = sign(r)`, `r' = 2r + b − σ·E` (b = fed dividend bit); the reverse needs σ,
  so σ is stored (quotient = history);
- **parity row** (new): `p = [r ≥ 0]`, `r' = 2r − E` if p else `2r + E + 1`, on the interval [−E−1, E]. Because E is
  odd and E + 1 even, **p = r' mod 2**. Forward: one adder `r' = 2r + (E ⊕ mask(p)) + 1` (masked operand,
  carry-in 1), exactly the non-restoring cost. Reverse: copy the low bit of r' into a control wire (CX), add
  `E ⊕ mask(¬p)` with carry-in 0, halve (relabel), erase the control from the sign of r (CX + X). Zero history.

Stage list (registers live):

| Stage | Change from the incumbent | Live registers |
|---|---|---|
| NORM, ROOT, LOAD1, COEF, SEED, PROD | none | as incumbent (peak 44) |
| DIV_FWD feed rows (≈ 12–15 rows, where product bits enter) | none (non-restoring, history kept) | v, flags, [r \| q \| P-rest] |
| DIV_FWD feed-free rows (rows 24–35 in the model: 12 rows) | parity rows; each digit copied (CX) from the remainder's low bit into its quotient wire | same; quotient complete at the endpoint (43 owners) |
| PHASE | ring by k, 33 rotations, inverse ring (unchanged) | endpoint registers |
| **Retire tail digits** | **HMR the 12 parity-row digits right after the inverse ring** (outcomes m_j) | 43 − 12 = 31 owners |
| DIV_CLEAN, parity rows | history-free reverse; at row j apply `Z(low bit of r_{j+1})` conditioned on m_j to repay the HMR debt (0 Toffolis) | 31 owners + up to 13 free wires for clean carries |
| DIV_CLEAN, feed rows | unchanged arithmetic, now with ~12 more clean carries available | ≤ 43 |
| LOAD2 … NORM_INV | none | as incumbent |

Optional extension F1b: restructure the dividend so that the large v-dependent part is preloaded into the initial
remainder and only the small h-correction is fed (relative weights 2^-12 … 2^-27). Then the top ≈ 10 rows are also
feed-free, their digits also history-free; this generalises Syed's `sigma12 = v11 ⊕ v10 ⊕ (v10 ∧ sigma13)`
identity (AR004) from one digit to all top digits.

**2. Qubit peak.** Unchanged at **44**: the peaks are the multiplier (v 12 + flags 4 + C 13 + acc 16 + carry,
minus 4 hidden bits) and the divider endpoint (43 owners + carry), and F1 changes neither. What changes is the
divider clean-up, which today runs with 42–43 data owners and buys width with Toffolis; with the 12 tail digits
retired it runs with ≈ 31 owners, i.e. with up to 13 spare wires at the same cap. F1 alone does not reach 43 Q:
jeffemerson's census puts 74 of the 88 uncovered 44-wide windows in the multiplier.

**3. Toffolis per stage.** *Estimate*, against the frontier:

| Stage | Incumbent | F1 | Reason |
|---|---:|---:|---|
| DIV_FWD | ≈ 270 | 270 | parity rows cost one masked adder, as non-restoring rows |
| PHASE ring | ≈ 105 | 105 | unchanged (the ring needs all 27 digits) |
| Tail retirement + debt repayment | — | 0 | HMR + classically conditioned Z |
| DIV_CLEAN | ≈ 320–330 (premium 38 at Q48, 53 at Q46 in v12, at least as large at 44) | ≈ 270–285 | the premium is width-bound (AR004: one retired history wire → −9.5); 12 retired wires remove most of it |
| Everything else | ≈ 746 | 746 | unchanged |
| **Total** | **1,450.9** | **≈ 1,400 (1,390–1,415)** | −35 to −60 |

Rotations: 33 generic + T, unchanged in count. Error budget: the parity tail's digits are binary digits of
`(r + E + 1)/(2E + 1)` instead of `r/E`; written in the tail's own scale, `r/E = (2 + 1/E)·T − 1 − 1/E` (T the
tail value), so the tail digits need weights proportional to `(2 + 1/E)`. With E ∈ [2^13, 2^14) the
E-dependence of the weight is ≤ 2^-14 of a term that is itself ≤ 2^-13 of y: a systematic error ≤ 2^-27 of y
before refit, ≤ ≈ 2^-29 (≈ 3% of the tightest tolerance) after the LP absorbs its mean. The frontier margin is
≈ 6% of the tightest tolerance, so the tail may cost up to half of it: at worst +0.56/4 per rotation ×33 ≈
+4.6 Toffoli-eq (margin halved). Netted in the range above.

**4. Cost.** C ≈ 1,640 (1,630–1,655), **Q × C ≈ 72,200 (71,700–72,800), −2.9% (−3.6% to −2.1%)**. On our own
46-qubit v12 build it would give ≈ 78,100 (−2.6% against our 80,163), which is still above the frontier: F1 beats
74,384 only on a 44-qubit base built from source (INV-A03, or our own generator in M3.2).

**5. Nulls it must escape.**
- *N9 / jeffemerson #5* (other dividers above the 706-Toffoli incumbent) and *INV-B06* (restoring tail priced
  +3n forward, +2n reverse per row): escaped, because a parity row costs the same single masked adder as a
  non-restoring row in both directions; the selection is the sign of r, as in non-restoring division.
- *N10* (no measurement-based uncomputation for the divider undo: "the undo borrows carry information") and
  *NR-0261* (clean-up is width-bound; reopen with "a different divider algorithm"): F1 is that different
  algorithm for the feed-free rows. For feed rows, section 5.2 proves N10 is right: they need history.
- *N8* (give the divider 47–48 wires: loses) and *Asterios B5* (retention at 46 wires impossible): F1 frees wires
  inside the existing cap instead of adding qubits.
- *INV-A16 / Asterios B2* (late-row truncation at its edge, margin used to 99.95%): the risk. The incumbent's late
  rows also use truncated divisors (`Aa` in `phase_words`) and an odd-remainder representation; the parity rows
  must be re-derived with truncation (E_j forced odd per row, interval nesting E_{j+1} ≥ E_j). If the tail's
  systematic error does not fit, keep 1–2 more non-restoring rows at the end of the feed section.

**6. Risks and the scaled prototype.** Risks: (a) the margin (above); (b) how much of the clean-up premium is really
due to history width rather than to the carry debts of the product bits (AR004 says history width is binding, but
the returns diminish with each extra clean carry); (c) F1 needs a 44-qubit source base to matter. Prototype:
section 6.

**7. P and effort.** P(primitive correct) 0.95 (proved exhaustively for small divisors); P(tail error fits the
margin) 0.7; P(≥ 35 Toffolis saved at the cap) 0.5; P(a 44-qubit base from source by January) 0.6 (M3.2) →
**P ≈ 0.3** to beat the frontier, ≈ 0.6 to save ≥ 2% on whatever divider build we have. Effort: prototype 8–14 h
(S2.19–S2.20); full scale 40–60 h, best done inside our own generator (M3.2) rather than as a patch on v12.

### F2. Divider-free piecewise polynomial

**1. Datapath.** Replace "root + division" by a polynomial evaluated with table coefficients. Two variants:

- *F2a, on n directly:* normalise; address a = top 7 bits of n (128 centred segments), t = the other 15 bits;
  y = c0(a) + t·(c1(a) + t·(c2(a) + t·c3(a))). With |t/A| ≤ 2^-8 the quartic term is 2^-34 (dropped) and the cubic
  term (5/16)·2^-24 is needed.
- *F2b, on the root coordinates (the best case):* keep the root, so y = R(v) − P with P = h·C(v) as today and
  R(v) = 1/Q evaluated by a centred cubic in the 6 low bits l of v with four 64-row words T0..T3 on the 6 high bits:
  R = T0 − l·(T1 − l·(T2 − l·T3)). This replaces only the division by a polynomial.

Stages (F2b): NORM, ROOT, LOAD/COEF/PROD (as today, P in a 16-bit accumulator), load T3; Horner w ← T2 − l·w,
w ← T1 − l·w (w grows to 19 bits); Y ← T0 − l·w − P (27 bits); unwind w; ring by k; rotate; mirror.

**2. Qubit peak.** At the last Horner product: Y 27 + v 12 + flags 4 + w 19 + carry 1 = 63 raw. Hiding v10, v11,
e2, e3 (−4) and building w in two halves (−9) gives **≥ 50**. The reason is structural: a Horner step needs its
multiplier (l, or t in F2a) and its multiplicand at the same time as the accumulator, and nothing is consumed in
place. A divider instead converts its dividend into quotient digits inside one register ([r | q | P-rest]), so its
operands never coexist with a separate accumulator. F2a is worse: t (15 bits) must stay live through every step,
so n (22) cannot be consumed, ≥ 60 Q.

**3. Toffolis.** F2b forward: NORM 51 + ROOT 96 + P path 185 (load 28, coef 44, seed 10, product 103) + R(v) ≈ 375
(four loads ≈ 112, products 6×7 + 6×13 + 6×19 ≈ 234, final add ≈ 27) + combine and ring ≈ 75 = ≈ 780; reverse
≈ 780 if no width premium is paid → ≈ 1,560 Toffolis. Rotations unchanged (33). The division it replaces costs
270 forward; the polynomial costs ≈ 375, because its 19-bit middle coefficient is multiplied by a 6-bit operand
while the divider adds a 13-bit divisor once per digit. **C ≈ 1,800 (≥ 1,740 if the loads are shared).** F2a:
products 15×6 + 15×13 + 15×19 truncated ≈ 485, four 128-row loads ≈ 250, C ≥ 1,900.

**4. Cost.** ≥ 50 × 1,740 = **87,000 (+17%)**; likely 53 × 1,800 = 95,400 (+28%). A null on both axes.

**5. Nulls.** Balram Srinivasan's divider-free polynomial (41–42 Q for the mantissa, ≥ 48–51 Q once the rotations
need the result at k-shifted absolute weights; reported in the S2.18 brief, notebook not in our snapshot) agrees
with this width estimate. N23 (bipartite tables and higher degree: partial products grow faster than the table
shrinks), NR-0198 (alternatives re-priced at 62 Q all ≥ 11% worse), NR-0221/0222/0224 (the quadratic Horner family
sits at its accuracy edge; rounding noise needs a bit of width everywhere). F2 does not escape them. The only new
element (the in-place cyclic ring on the accumulator, which removes the extra ring frame) does not touch the
binding constraint, the Horner temporaries.

**6. Risks / prototype.** Not recommended. A cheap falsification, if wanted: a register-allocation model of F2b on
the scaled instance (n_p = 5) to confirm the ≥ 50 bound.

**7. P and effort.** P ≈ 0.02; 80 h (prototype 20 h).

### F3. Absolute-frame streaming quotient (ring-free)

**1. Datapath.** Remove the ring and, ideally, the stored quotient: make division row t always produce the digit
of absolute weight 2^-t, rotate it at once with a fixed angle, and drop it. The k-shift must then be applied to the
dividend instead of the quotient: (a) delay the dividend by s = 10 − k rows (barrel-shift P inside a 25-wire frame
before the division); or (b) divide by an un-normalised √x carried to 25 significant bits (a 36-bit fixed-point
divisor), which is the width disaster that normalisation exists to avoid.

**2. Qubit peak.** (a): 37 rows instead of 27, and **streaming is impossible**: because the feed start depends on
k, every row is a feed row for some x, so by section 5.2 every digit needs history. The endpoint holds 37 history
digits (+10 against today): raw v 12 + flags 4 + r 14 + q 37 = 67, ≈ 53 after the incumbent's kinds of retirement,
**≥ 46 optimistic, ≥ 50 likely**. (b): ≥ 80.

**3. Toffolis.** +10 rows × (11.25 forward + 12.8 clean-up) ≈ +240; a barrel shift of the dividend in a growing
frame (15 → 25 bits, ≈ 15 + 16 + 18 + 22 = 71 CSWAPs) costs as much as the quotient ring it replaces (63 CSWAPs),
so ±0; rotations unchanged (the absolute positions 2^-4 … 2^-36 are the same 33). **C ≈ 1,930.**

**4. Cost.** ≥ 46 × 1,930 = **88,800 (+19%)**, likely ≥ 96,500 (+30%).

**5. Nulls.** N12 (jeffemerson: bit-serial quotient rotation, split endpoint, absolute-frame division, phase-gradient
rotations; 8 redesigns all above the incumbent) and NR-0174 (the shift range spreads the significant bits over
~34 fixed-angle positions). Not escaped; section 5.2 adds the reason the streaming variant cannot shed its
history.

**6. Risks / prototype.** Not recommended.

**7. P and effort.** P ≈ 0.01; 60 h.

### F4. Table-free reciprocal-root recurrences

**1. Datapath.** One recurrence from n to y with no table and no division:
(a) Ercegovac–Lang digit recurrence for 1/√n (digits of y, residual w_j = 2^j(1 − n·y_j²), running product n·y_j);
(b) multiplicative normalisation, the hyperbolic-CORDIC form: choose s_j ∈ {0, 1} so that z = n·Π(1 + s_j 2^-j)² → 1,
in place z ← z + s_j(2^{1−j} z + 2^{−2j} z), and accumulate y ← y + s_j 2^-j y;
(c) a table seed with one Newton step y = y0(3 − n·y0²)/2 (the reference family), or Goldschmidt.

**2. Qubit peak.** (a), (b): n or z (26) + y (26) + residual or running product (26) + the s_j history (25, since a
step z ← z(1 + 2^-j)² is ambiguous to reverse exactly as a non-restoring row) → **≥ 74** (N17). (c): public lines
at 62–99 Q (NR-0235: Newton 3.02× at 99 Q).

**3. Toffolis.** (b): per digit two shifted adds on z and one on y, ≈ 78; × 25 digits ≈ 1,950 forward, ≈ 2,000
reverse, **C ≥ 4,000**. (c): n·y0² is a 22 × 24 product (≈ 530) plus y0² and y0·e (≈ 230): ≈ 760 forward,
C ≈ 1,850 at ≥ 62 Q.

**4. Cost.** (a)/(b) ≥ 74 × 3,000 = **≥ 220,000**; (c) ≈ 62 × 1,850 = 115,000 (+54%).

**5. Nulls.** N17, NR-0235, NR-0029/0030 (reference precision lattice), INV-B02, INV-B07. Not escaped. The common
reason: none of these recurrences consumes n in place the way the consuming root does (the root turns n into
(v, h) with no garbage for 96 Toffolis), so they carry n, y and a residual at once.

**6. Risks / prototype.** Not recommended. INV-B02 stays open only if a two-register formulation is found
(jeffemerson's direction 5); nothing here suggests one.

**7. P and effort.** P < 0.01.

### F5. Phase-additive split with an exponent-indexed coarse phase

**1. Datapath.** Use the linearity of the phase: `φ = 2^k·Y_c(a) + 2^k·(y − Y_c(a))`. The coarse term depends only
on (a, k) (a = 6 high bits of v), so it can be applied early and on its own: load Y_c(a) (7 significant bits) from a
64-row table, barrel-shift it by k into a 17-wire absolute frame, rotate the 17 wires, unshift, unload. The fine
term needs only ≈ 19 significant quotient digits: divide `N − Y_c(a)·E` by E (dividend correction by a 7 × 6
product plus a table word), with 19 rows instead of 24–27, a ring on 19 digits, and rotations on 19 + 10 = 29
positions.

**2. Qubit peak.** The coarse stage runs after the root: v 12 + h 10 + flags 4 + 17-wire frame + ≈ 3 shift and loader
wires = 46; the rest is the incumbent's 44. **Q = 45–46** (the frame could be shortened by rotating in two passes,
at the cost of more CSWAPs).

**3. Toffolis (change against the incumbent).** Coarse load, shift, unshift, unload ≈ +130 (two 28-AND loads, two
7-bit shifts into a 17-bit frame ≈ 35 CSWAPs each); dividend correction ≈ +90 (both directions); divider −8 rows
≈ −190; shorter ring ≈ −25; rotations 17 + 29 = 46 instead of 33: +13 × 7.23 + 46 × 0.067 ≈ +97. Net ≈ +100 to
+130, **C ≈ 1,800**.

**4. Cost.** 45 × 1,800 = 81,000 to 46 × 1,820 = 83,700: **+9% to +13%**.

**5. Nulls.** Floor 2 of section 2.3: two separately phased terms pay two sets of absolute positions (the 13 extra
rotations alone cost more than the 8 divider rows saved); N2/N3 (rotations saturated) and NR-0242 (folding coarse
bits into the fine word needs both live); N9 (the Y_c·E correction is a "reciprocal by multiply" step). Not
escaped. Sharing the rotations would require adding the coarse word into the quotient before the ring, which is
the incumbent's seed.

**6. Risks / prototype.** Not recommended.

**7. P and effort.** P ≈ 0.03; 50 h.

### F6. Logarithmic number system

**1. Datapath.** L = log2(n) in fixed point with 26 fractional bits (the integer part comes from k and the
normalisation), u = −L/2, y = 2^{frac(u)} by a second table + polynomial, then the k-ring and the rotations. The idea
is that the exponent becomes an addition and the halving a relabelling.

**2. Qubit peak.** log2(n) to 2^-26 cannot be computed in place: n (22) stays until L (26) is complete, plus the
polynomial temporaries (≈ 19) and flags: **≥ 55** (≥ 60 without heroic retirements).

**3. Toffolis.** Each of log2 and 2^u to 26 bits is an F2a-sized evaluation (64–128-row tables, cubic in a
15-bit offset, ≈ 800 forward), so ≥ 1,600 forward, ≥ 3,200 with the mirror, **C ≥ 3,000** even with shared loads.

**4. Cost.** ≥ 55 × 3,000 = **≥ 165,000 (+120%)**.

**5. Nulls.** NR-0235 (every non-Horner numeric family priced at ≥ 3× without an x-dependent branch). Not escaped:
LNS replaces one smooth evaluation with two of the same difficulty; it pays only when many multiplications share
one conversion, which is not the case here.

**6. Risks / prototype.** Not recommended.

**7. P and effort.** P < 0.01.

### F7. Direct Boolean synthesis of the phase bits

**1. Datapath.** This is the open backlog item of `PLAN.md` §11 ("direct synthesis of the phase bits as
low-degree Boolean functions of x"). Write the phase as a sum of rotations on Boolean features of x, either per
output bit (bit p of 2^k·y at absolute weight 2^-p, p = 4..36, each a Boolean function of 22 inputs) or as a
Walsh/ANF phase polynomial `Σ_S θ_S·Π_{i∈S} x_i`; compute the needed monomials with ANDs, rotate them, uncompute
by measurement.

**2. Qubit peak.** x 22 + monomial scratch; ≥ 44 for any useful amount of sharing.

**3. Toffolis and rotations.** The bits of weight 2^-12 and below are carries of an irrational expansion: each
depends on all 22 input bits with full degree and behaves like a pseudo-random function, so its multiplicative
complexity grows exponentially in the number of inputs. In the Walsh form every monomial with a distinct angle is
its own rotation (≈ 7.2 Toffoli-eq each), and the significant coefficients number in the thousands. NR-0198 priced
the Walsh family at **12.2×** the 62-qubit pick. **C ≥ 10,000.**

**4. Cost.** ≥ 44 × 10,000 = **≥ 440,000 (≥ +490%)**.

**5. Nulls.** NR-0198 (Walsh 12.2×), N19/NR-0085 (no affine or low-degree structure in the frontier's products).
Not escaped. Only the top 3–8 phase bits are low-degree functions of the high bits of x, and those are the cheap
Clifford/T positions anyway. This closes the backlog item as a null.

**6. Risks / prototype.** Not recommended.

**7. P and effort.** P ≈ 0.

## 5. Cross-cutting findings

### 5.1 The parity division row (new primitive)

Let E be odd and I = [−E−1, E]. Row: `p = [r ≥ 0]`; `r' = 2r − E` if p = 1, `r' = 2r + E + 1` if p = 0.

- *Invariance:* r ∈ [0, E] gives 2r − E ∈ [−E, E]; r ∈ [−E−1, −1] gives 2r + E + 1 ∈ [−E−1, E−1]. Both lie in I.
- *Bijection:* the p = 1 images are the E + 1 odd numbers of [−E, E], the p = 0 images the E + 1 even numbers of
  [−E−1, E−1]; together they are exactly I.
- *Digit readable from the next remainder:* `p = r' mod 2`.
- *Value:* after J rows, `Σ_j p_j·2^(J−1−j) = (2^J·r_0 + (E+1)(2^J − 1) − r_J) / (2E + 1)`: the digits are the binary
  quotient of `(r_0 + E + 1)/(2E + 1)`, with an error below one unit in the last place.
- *Cost:* forward `r' = 2r + (E ⊕ mask(p)) + 1` (since −E = ¬E + 1): one adder with a masked operand and carry-in 1,
  the same as a non-restoring row (which uses carry-in p). Reverse: `c ← r' mod 2` (CX), `r ← r' + (E ⊕ mask(¬c))`
  with carry-in 0 (adds E if c = 1, ¬E = −E−1 if c = 0), halve by relabelling, then clear c from the sign of r (CX, X).
  No stored digit.
- *Measurement-based retirement:* once a digit wire q_j is no longer needed (after the phase), `HMR q_j → m_j`
  leaves the debt (−1)^{m_j·p_j}; in the reverse pass p_j sits on the low wire of r_{j+1} when row j is undone, so
  `Z(low wire)` conditioned on m_j repays it. Zero Toffolis.
- The reverse row is a Montgomery (REDC) step (add a multiple of the odd modulus chosen by parity, then halve): the
  parity divider is Montgomery reduction run backwards, and a reversible circuit can run either direction at the
  same price.

Checked exhaustively for every odd E in 3..199 (all r_0, 12 rows), and on 20,000 random hybrids of 15
non-restoring feed rows followed by 12 parity rows with E in [2049, 4095], the tail reversed without stored
digits (Appendix A).

### 5.2 Feed rows cannot be history-free (impossibility lemma)

*Claim.* Consider a row `r' = 2r + b + A` (if r < t) or `r' = 2r + b − B` (if r ≥ t), with a fed dividend bit
b ∈ {0, 1}, a threshold t on r, and offsets of different parity (A ≢ B mod 2, the condition for reading the digit
from r' mod 2 once b is known). No interval [L, U] is invariant.

*Proof.* Both images must stay in [L, U] for both values of b. Left set: 2L + A ≥ L, so L ≥ −A, and
2(t − 1) + 1 + A ≤ U. Right set: 2t − B ≥ L, and 2U + 1 − B ≤ U, so U ≤ B − 1. Hence 2t ≤ U − A + 1 ≤ B − A and
2t ≥ L + B ≥ B − A, so 2t = B − A, which needs B − A even: a contradiction. ∎ (Non-restoring rows have A = B = E,
B − A = 0, t = 0. A selection that looks at b as well needs a comparison before the add: the restoring row of
INV-B06, about 3n instead of n Toffolis.) An exhaustive search over intervals and thresholds for E ≤ 61 finds none
(Appendix A).

*Consequence.* A fed bit adds one bit of information to a bounded remainder, so one bit per feed row must leave as
stored history; only feed-free rows can be history-free. In the incumbent the dividend (P, 15 bits, ≈ 10 bits of
information about h) enters in the first 12–15 rows, which therefore keep their history whatever the recurrence.
This confirms N10 and NR-0261 for those rows, and caps F1 at the 12 feed-free rows (≈ 22 with F1b).

### 5.3 Why the root + divider family is hard to beat

- The **consuming root** turns n into (v, h) in place, with no garbage, for 96 Toffolis per direction, and leaves a
  correction ε = h/Q² ≤ 2^-11 that a single 13 × 10 product covers. Without it, a table-based start needs 128–4,096
  rows or a cubic (F2, F5).
- The **divider** converts its dividend into quotient digits inside one register; its operands never coexist with a
  separate accumulator. A polynomial keeps multiplier, multiplicand and accumulator at once (F2: ≥ 50 Q), and its
  middle product (19 × 6 bits) costs more than a 13-bit add per digit.
- The **ring** is cheaper than any alternative way of applying the k-shift (dividend shift in F3: same CSWAPs and
  +10 rows; separate coarse phase in F5: +13 rotations).
- The **mirror** is forced: normalisation, root, product (in place on h) and division are in-place bijections, so
  they cannot be uncomputed by measurement (N10, N28), except in feed-free division rows (5.1).

So the families that differ in *how the phase is computed* all lose on paper, and the remaining gains are in
*how the same computation is uncomputed and laid out* (F1, the in-family width moves INV-A21–A24, F1b).

### 5.4 What it implies for the leaders

The best invsqrt we can price from our own source is F1 + the in-family width moves on a 44-qubit base:
−3% to −9% of task cost (≈ +1% to +3% on the score). That is far from the 1.42× cost-product advantage of the top
three. Unless a family outside this list exists (INV-B01's evolutionary search is the tool to look for one), the
leaders' gain is more likely on dqi and chemistry than on invsqrt.

## 6. Recommendation for S2.19–S2.20

### Rank 1: F1, the parity divider (prototype)

**E1 — exact-model check of the error budget (S2.19, 2–4 h, runs in minutes).** A new script beside
`generators/invsqrt/team/table_search.py` (e.g. `parity_tail_model.py`) that imports `domain`, `coefficients` and
`seeds` and copies `phase_words` with a parameter `j_tail` (default 24): rows j ≥ j_tail use the parity update on
the full divisor (version A, no truncation), then on per-row truncated odd divisors E_j with nested intervals
(version B, as the incumbent's `Aa` rows). The tail digits enter the phase word with weights `(2 + 1/E0)·2^-j`
scaled, E0 a fitted constant. Refit the 34 angles by the exact LP (the phase is linear in the angles: HiGHS via
`scipy.optimize.linprog`, constraint generation over all 3,133,452 inputs, exact integer recheck) and report
failures, min margin and worst error/tolerance against the unmodified model under the same LP, for
j_tail ∈ {24, 26, 28, 30}. *Go* if some j_tail ≥ 24 keeps the margin ≥ 50% of the reference (rotation cost
≤ +4.6 Toffoli-eq); record the curve (rows freed vs margin) either way. Complexity O(27·N) per evaluation (about 1 s
vectorised), memory ≈ 0.5 GB.

**E2 — scaled kickmix divider, width-capped (S2.20, 6–10 h, seconds to verify).** A small generator
(`generators/invsqrt/team/proto_parity_divider.py`) with the kickmix builder API of `nr_reference_circuit.py`:
Q of 6–7 bits (E = 2Q + 1), a 6-bit fed dividend in 6 non-restoring rows, then 6 parity rows; 12 quotient digits
rotated with fixed angles 2^-(j+1) turn (k = 0, no ring). Baseline: the same with 12 non-restoring rows and full
history. F1: HMR the 6 tail digits after the rotations, reverse them history-free with conditioned Z repayment,
mirror the feed rows. Check every (Q, dividend) pair × 256 HMR seeds with `km.Simulator`: inputs restored, scratch
zero, phase equal to `Σ θ_j q_j` to the 2^-64 unit. Then run both under width caps (full, −2, −4 wires) with one
adder policy (Gidney adder with n − 1 clean carries when wires are free, a carry-free in-place adder at about 2n
Toffolis otherwise) and report clean-up Toffolis. *Go* if, at a cap where the baseline pays a premium, F1's
clean-up recovers ≥ 30% of it.

**E3 — entry to M3.2, not part of S2.20.** The whole datapath (normalise, root, table, product, hybrid divider,
ring, rotations) on the official verifier's scaled instance (`--n-p 5`: x ≤ 2,700 on 12 qubits, or `--n-p 6`:
x ≤ 11,532 on 14 qubits, with `--precision-bits` scaled), built by our own generator with and without the parity
tail. This is the first circuit of our own invsqrt line and doubles as the GEN-A13 harness.

### Rank 2: no second family

None of F2–F7 deserves prototype time: each loses by ≥ 9% for a structural reason that a prototype would only
confirm. The second slot of S2.20 is better spent on the **scaled-proxy harness (GEN-A13)** that E3 and M3.2 need:
a parametrised version of the exact integer model (`domain()` and `phase_words` generalised to n_p, the precision
bits and the table size), plus kickmix circuits checked by the official verifier at `--n-p 5/6`. In M3.2, F1
should be stacked with F1b and the width moves INV-A21/A22/A24 on our own 44-qubit base; that combination is the
realistic −5% to −9%.

What would change this ranking: a public invsqrt below 44 Q or below ≈ 1,400 mean Toffolis (read its family before
anything else); E1 failing on the margin (F1 then shrinks to F1b or to the few rows that fit); or an INV-B01 search
that produces a datapath not covered here.

## 7. Proposed register updates (not applied in this step)

- `Analysis/IDEAS.md`: new **INV-B10** "history-free (parity) divider rows + measured-out tail digits" (gain 2.1–3.6%,
  P 0.3, 0.6 on our own divider build, effort 50 h; TODO, top of the B list); **INV-B06** superseded by INV-B10
  (restoring rows stay priced above; parity rows cost a non-restoring row); INV-B02, INV-B07, INV-B08 → NULL on paper
  (F4, F2; reasons above); `PLAN.md` §11 backlog item "direct synthesis of the phase bits" → NULL (F7).
- `experiments/NULLS.md` (our research): F2 (polynomial, ≥ 50 Q, +17%), F3 (absolute-frame streaming, +19%),
  F4 (table-free recurrences, ≥ 74 Q), F5 (phase-additive split, +9% to +13%), F6 (LNS, ≥ +120%), F7 (Boolean
  synthesis, ≥ +490%), and the lemma of 5.2 (feed rows need history).
- `experiments/LEDGER.md`: one row for S2.18 (paper, verdict: F1 to prototype, six nulls).
- `changelog/`: an S2.18 entry; `Analysis/PLAN.md`: S2.18 DONE, S2.19 = E1 + harness, S2.20 = E2.

## Appendix A. Exhaustive check of the parity row (sections 5.1 and 5.2)

Plain Python 3, no dependencies; all assertions pass in 0.7 s (run 2026-11-19 with CPython 3.11). It checks: the
parity row is a bijection of [−d−1, d] with the digit equal to the next remainder's low bit and the quotient
identity of 5.1, for every odd d in 3..199; no interval/threshold is invariant under a one-bit feed for odd d in
3..61 (5.2); 20,000 random hybrids (15 non-restoring feed rows, then 12 parity rows) reverse their tail without the
stored digits.

```python
import random

def parity_row(r, d):                      # d odd; interval [-d-1, d]
    p = 1 if r >= 0 else 0
    return (2 * r - d if p else 2 * r + d + 1), p

def check(d, J=12):
    I = range(-d - 1, d + 1)
    assert sorted(parity_row(r, d)[0] for r in I) == list(I)          # bijection of the interval
    for r in I:
        r2, p = parity_row(r, d)
        assert p == (r2 & 1)                                           # digit = low bit of next remainder
        assert (r2 - (-d if p else d + 1)) // 2 == r                   # history-free reverse
    for r0 in I:                                                       # quotient identity
        r, P = r0, 0
        for _ in range(J):
            r, p = parity_row(r, d)
            P = 2 * P + p
        assert 2**J * r0 + (d + 1) * (2**J - 1) - r == (2 * d + 1) * P

def feed_invariant(d):                     # one-bit feed b in {0,1}: search intervals and thresholds
    for L in range(-d - 3, -d + 3):
        for U in range(d - 3, d + 3):
            for t in range(L, U + 1):
                if all(L <= (2*r + b - d if r >= t else 2*r + b + d + 1) <= U
                       for r in range(L, U + 1) for b in (0, 1)):
                    return (L, U, t)
    return None

for d in range(3, 200, 2):
    check(d)
assert not any(feed_invariant(d) for d in range(3, 63, 2))

random.seed(1)
for _ in range(20000):                     # hybrid: 15 non-restoring feed rows, then 12 parity rows
    d = random.randrange(2049, 4095, 2)
    r = random.randrange(-d, d)
    for b in [random.randrange(2) for _ in range(15)]:
        r = 2 * r + b - (d if r >= 0 else -d)
        assert -d <= r < d
    rm, tail = r, []
    for _ in range(12):
        r, p = parity_row(r, d)
        tail.append(p)
    for p in reversed(tail):               # reverse without the stored digits
        q = r & 1
        assert q == p
        r = (r - (-d if q else d + 1)) // 2
    assert r == rm
```
