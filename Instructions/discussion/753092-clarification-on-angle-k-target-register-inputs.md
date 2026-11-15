# Clarification on ANGLE_k target-register inputs

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753092 — copied 2026-11-15 (raw Markdown, unedited).

**Chat8903** · 2026-10-27 05:30 UTC · 4 votes

Must every ANGLE_k satisfy the XOR transformation in Eq. (72) for arbitrary target-register inputs, or only for states reachable in the specified forward/reverse traversal?


---

## Comments (4 threads)

Replies are indented with `>`, one level per reply depth.

**manwithacat** · 2026-10-28 01:23 UTC · 2 votes

Thanks to Tanuj for the ruling, and to everyone who pushed on the ANGLE_k question. Following up on my earlier comment with something checkable. As before, this was done with AI assistance (Claude wrote the checker and ran the experiments under my direction), and my linear algebra is still mostly vibes, so please treat it as a tool to poke at, not an authority.

**What I ran.** First the sampled strict checker (random values in the word registers, not just the previous batch): no violations on any of the chemistry circuits I hold. Then a proof-style check of the *contract*, for all measurement outcomes rather than the 64 the verifier samples.

**Result.** 17 distinct chemistry circuit sets are proven (11 taken from public notebook outputs, the organisers' sample submission, and 5 of my own research builds), for every address, in both traversal directions and for all measurement outcomes, under both readings of the contract: the reachable-only one the organisers ruled is official (z = batch k-1 forward, batch k in reverse), and the stricter "XOR for every z" one. That's 34 checks with no counterexamples. The checker does see real faults: a deliberately shortcut toy (a loader that assumes z) passes reachable but fails strict, as it should, and removing one conditioned CZ from a real public family gives a counterexample at 142 addresses in each direction, which replays in kickmix.

**Method.** For each address and circuit, the kickmix semantics are encoded as XOR/AND over GF(2): one variable per bit of the target register (strict only) and one per measurement outcome. The check asks whether words, workspace or phase can differ from the specification. On every real circuit this reduced to zero by exact symbolic (ANF) evaluation after fixing the address, so the "proof" is an exact algebraic evaluation. A SAT solver (CryptoMiniSat) is the fallback and was exercised on the toys and on mutated circuits.

**Limits, stated plainly.**
- Trust rests on my encoder of kickmix's semantics. It agreed with kickmix on 223,040 random cases (0 mismatches, including the measurement outcomes, by replaying kickmix's RNG), and the mutation controls behave as they should. That is strong evidence, not a proof that the encoder is right.
- It covers these 17 circuit sets only, at the addresses and word/batch structure the verifier uses.
- It models a zero initial workspace and an empty condition stack; Z_POW is unsupported.
- It says nothing about the salt. Whether final scoring re-seeds the sampling is something I do not know.
- Proving the contract does not make a circuit cheaper or better; it only says the contract holds.

**Try it.** The checker, the results table, the toy families and a pinned solver wheel are in a dataset, [vqcc-angle-k-formal-check](https://www.kaggle.com/datasets/manwithacat/vqcc-angle-k-formal-check), and [this notebook](https://www.kaggle.com/code/manwithacat/prove-angle-k-for-every-measurement-outcome) reruns it on the organisers' sample family or on your own folder, so you can check the checker (it takes a while: about half an hour on Kaggle for the sample family). A short technical note on Zenodo will follow. Prior work on SAT/BDD equivalence checking of reversible circuits (Wille et al. and others) is the foundation; I think only the application to kickmix's contract is new here, and I'd welcome pointers if it already exists.

Corrections very welcome, particularly from anyone who can break the encoder.

> **manwithacat** · 2026-10-28 01:28 UTC · 0 votes
>
> A quick note on *why* you might want this, because it's a subtle one.
>
> The chemistry verifier tries 64 random measurement outcomes per address, seeded from your circuit's hash and a salt. The scoring README says the final grader must use a private, unpredictable salt, so the outcomes drawn at final scoring may not be the ones you've seen. A circuit with a fault that only shows up on a rare combination of outcomes (say a missing fix-up after a measurement-based uncompute, or a correction conditioned on the wrong bit) can pass every local run and the public LB, then fail when the outcomes are re-drawn, and a failed task sinks the whole submission.
>
> An illustrative example (hypothetical numbers): a fault that triggers on 1 in 1,000 outcomes at a single address slips past 64 samples about 94% of the time, and a fresh draw hits it about 6% of the time. Sampling can only make that unlikely; the formal check says no outcome breaks the contract at any address.
>
> So it's worth running once before you pick your final submissions if you've edited chemistry circuits, especially measurement-based uncomputation or hand-optimised loaders, or if you're building cheaper loaders on the reachable-only ruling (you want `reachable` PROVEN; failing `strict` is fine). It isn't a faster verifier or a score check: it takes minutes to half an hour and says nothing about cost. For what it's worth, it found no such faults in any of the 17 sets above.

**Tanuj Khattar** — **Host** · 2026-10-27 19:12 UTC · 2 votes

The ANGLE_k circuits only need to satisfy the transformation on states reachable in the specified forward and reverse traversals (for valid addresses 0 <= a < N + R*B), not on arbitrary target-register inputs. This is also what the official verifier enforces.

The explanation by @manwithacat above is correct.

> **manwithacat** · 2026-10-27 19:40 UTC · 0 votes
>
> Interesting. That implies a possible gain in chemistry. Circuits that only work on the inputs the algorithm actually produces are now confirmed legal. A word-loading circuit can therefore assume z holds the previous batch and use that to save Toffolis. @egortrushin (28th) says such circuits already exist. That makes it a concrete idea for a cheaper chemistry loader, and possibly part of how the private leaders are ahead.

**manwithacat** · 2026-10-27 17:47 UTC · 0 votes

Good question. I went through the paper's block encoding (with AI assistance) to see whether the strict reading is
ever needed. A caveat up front: my linear algebra is not particularly strong, and I'm partly going on vibes here, so
someone with a bit more youth and mathematical rigour may well offer a better analysis. Corrections very welcome.

**Where ANGLE_k is used.** In the block-encoding chain (App. A.4, Eq. 100), the word register z starts at 0. In the
forward chain, ANGLE_k is followed by Givens batch k for k = 0 … K−1, then Majorana select. The reverse chain runs
Givens batch k⁻¹ followed by ANGLE_k for k = K−1 … 0, then Inner†, reflect and square, and Outer†. Sec. III.3.6 and
App. A.1 describe exactly this: ANGLE_0 loads batch 0, each later ANGLE_k updates the register in place to batch k, and
the reverse traversal clears it again. The Givens rotations are controlled *by* z and do not change it, and nothing
else in the chain touches z. So whenever ANGLE_k runs, z holds batch k−1 of the current address on the way forward,
or batch k on the way back. In both cases the address register may be in superposition.

**Why that's enough.** These circuits act on basis states as classical reversible logic plus phases that are multiples
of π. A circuit that is exactly right (words, clean workspace, unflipped phase) on every basis state |a⟩|ϑ^(k−1)(a)⟩
and |a⟩|ϑ^(k)(a)⟩ is right on every superposition of them, by linearity. Those are precisely the inputs
`recover_words_from_circuits` simulates: every address, z reset to the previous batch going forward and to the current
batch going back. The "self-inverse" property in Sec. III.3.6 is then only ever used on those states, and the cost
model (one forward plus one reverse traversal) charges exactly those uses. As far as I can see, the XOR form of Eq. (72)
for arbitrary z is a convenient specification that the algorithm never needs.

**If the organisers prefer the strict reading,** it can be checked classically and cheaply: write random values into
the word registers and check z ⊕ f_k(a) on a few thousand shots per address. Every basis input is one classical
trajectory in kickmix, so no quantum hardware is involved.

Since circuits of both kinds apparently pass today, a ruling would help everyone. @ryanbabbush

**Egor Trushin** · 2026-10-27 17:32 UTC · 0 votes

I would also appreciate clarification on this question, since solutions do exist that pass official verification and are correct throughout the prescribed forward/reverse traversal, but fail the XOR requirement for some other target inputs. Are such solutions considered valid submissions?
@ryanbabbush
