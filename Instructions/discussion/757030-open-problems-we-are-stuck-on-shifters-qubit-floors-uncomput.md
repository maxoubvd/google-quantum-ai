# Open problems we are stuck on: shifters, qubit floors, uncompute and loader bounds

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757030 — copied 2026-11-15 (raw Markdown, unedited).

**Rustam Bazarbayev** · 2026-11-08 16:07 UTC · 4 votes

Hi all,

This competition has turned into a nice collection of small, hard, well-defined combinatorics problems. Several of them have stopped us cold, and I suspect some of you have thought about them too. Below are the ones I find most interesting, posed as open questions. I'm not sharing circuits here, only the problems. Partial answers, pointers to papers, or "we proved this is impossible" are all welcome. Negative results save everyone time (the mu_4(6) = 14 thread, 756476, was a great example).

---

## 1. Minimum controlled-SWAP count for a normalizing barrel shifter

**Setting.** A w-bit register holds x with x != 0. We want to shift it left by r = (w - bitlen(x)) so the leading one lands in the top position. The binary digits of r are already available on control qubits (computed earlier, free to read). Each CSWAP costs 1 Toffoli. Controls may be any single available control bit (or a cheap AND of them that is already present).

**Question.** What is the minimum number of CSWAPs for w around 21 and about 22 possible shift amounts? The textbook log-depth barrel shifter is far from optimal, because many positions are known to be zero for a given r. Public circuits are now in the high 40s per shifter.

What I know:
- A SAT encoding (fixed N gates, choose controls and swap pairs, require correctness on every r) solves small cases quickly. With free controls, optimal counts I got are n=6: 7, n=7: <=9, n=8: <=11, n=9: <=13, roughly 2n - 5.
- At w = 21, N in the high 40s, a CaDiCaL run of about 11 hours gives neither SAT nor UNSAT.
- A simple counting argument gives a lower bound in the high 30s, so the gap is large.

Is there a closed form, a better encoding (symmetry breaking on swap order? layered formulation?), or a known result from sorting-network or permutation-network literature that applies when the set of permutations is a fixed list of shifts with zero padding?

## 2. Is about 62 qubits a real floor for the inverse square root task?

Every strong public invsqrt circuit I've looked at sits near 62 qubits. When I measure liveness over the input domain, the plateau is wide, not a single peak: it lasts most of the circuit, and the scratch bits look like dense data (about half ones), not mostly-zero flags.

Things that did *not* help in my tests: bipartite/multipartite tables, higher-degree polynomials with fewer segments (the partial-product count grows faster than the table shrinks), and more measurement-based uncompute (multi-bit sums need recompute, which eats the saving).

**Question.** Is there an argument (information-theoretic or from reversible pebbling) that a Newton/polynomial-style invsqrt with this output precision needs about 3*w live bits? Or has someone found a scheme that trades a few hundred extra Toffolis for, say, 8 fewer qubits? Under the Q*C cost model one qubit is worth about 25 Toffolis here, so the trade only pays if it's very lopsided.

## 3. Cheap uncompute for long chains of in-place field multiplications

Many decoder-style circuits here update a large register in place, multiplying by field scalars that are themselves produced along a recurrence. Running such a block backwards (to clean up) is often not as cheap as running it forwards, because the inverse step needs values the forward step never had to keep.

**Question.** Are there general techniques, from reversible computing or from classical Reed-Solomon / BCH hardware, for structuring a chain of in-place GF(2^m) updates so that its inverse is as cheap as the forward pass? For example, a shared evaluation sequence, a different normal form of the key equation, or a pebbling schedule suited to this kind of recurrence. Pointers to papers are very welcome.

## 4. Multiplication rank: closed, thanks to the community

For completeness: GF(2^12) as GF(4)^6 needs 3 * 14 = 42 bilinear multiplications, and thread 756476 settled mu_4(6) = 14 (asymmetric included, with a Lean proof). I had a flip-graph search running for rank 13. It can stop. Thank you to the people who did that work.

The remaining multiplier question I'd still like answered: is there any *non-bilinear* trick (for example, reusing one operand across many products when one side is a fixed scalar chain) that beats 42 Toffolis per GF(2^12) product in a reversible setting?

## 5. Chemistry: loader cost vs. a ~2*sqrt(m*L) bound, and co-designing the words

For the DF/THC task the angle loader dominates. A table of L words of m bits loaded with a standard select-swap / basis-plus-products construction costs on the order of 2*sqrt(m*L) Toffolis, and public loaders are already close to that. An exact Mobius-style (ANF) synthesis of the same table is much worse, because the words look random.

Two observations:
1. The words are not fixed. They come from training, and each is a floor/ceil choice near a continuous optimum. So the real problem is *joint*: pick roundings so that the table is cheap to load **and** the energy error stays under the gate.
2. The task metric has hard floors (E_gap and Lambda constraints). Good circuits already sit almost exactly on them, so there is very little lambda left to win. The gain has to come from cheaper loading.

**Question.** Has anyone tried loader-aware rounding, for example adding a penalty for the GF(2) rank of the products a deep XAG synthesizer would need, or a differentiable proxy for loader Toffoli count? Or is there a loader construction that beats the select-swap bound for tables with a known structure (low-rank, sign-gauge freedom)?

---

If you have a proof that any of these is optimal as it stands, I'd love to see it. Otherwise, happy to discuss approaches in the comments. Good luck everyone.
