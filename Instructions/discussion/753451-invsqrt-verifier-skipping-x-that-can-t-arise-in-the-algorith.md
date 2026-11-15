# invsqrt verifier: skipping x that can't arise in the algorithm?

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753451 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-28 05:46 UTC · 3 votes

A small suggestion about verifier cost, not about the contract for scoring purposes.

In the stopping-power paper the task comes from (Rubin et al., arXiv:2308.12352, §II.4), the inverse square root is applied to the sum of squared x, y, z coordinate differences between two particles. With n_p = 10, its input is always

**x = a² + b² + c², with integers |a|, |b|, |c| ≤ 1022 and x ≥ 1**,

which is also where the verifier's max_x = 3·1022² = 3,133,452 comes from. That set is exactly defined and cheap to compute: enumerate the 1,023² sums a² + b², add each c², and mark the results. It takes under a second with numpy and gives a 391,682-byte bitmask (bit k−1 for x = k, little-endian). We get **2,395,723** reachable values, **76.5%** of 1..3,133,452. The rest are:
- **522,239** of the form 4ᵃ(8b + 7), never a sum of three squares (Legendre);
- **215,490** sums of three squares that need a coordinate above 1022. They all lie at x ≥ 2²⁰, starting with 2²⁰ = 1024² itself, since powers of 4 have only the trivial representation.

Our bitmask has SHA-256 `2d649ab768c29d3d034d565cf9ce138b88618e864084ee5ea1d820908f6bd258` if anyone wants to compare.

We checked whether a circuit could exploit these gaps. It essentially can't: they're scattered among valid neighbours and 1/√x is smooth, so tables, normalisation and arithmetic don't change. The only effect is a little slack in the rotation budget, about 0.08% of invsqrt cost for the current public circuits. So this isn't a request to relax the problem for scores.

The reason to raise it is compute. Every full invsqrt verification (locally, on Kaggle and at scoring) spends about a quarter of its 200 million shots on inputs the algorithm never produces. Restricting the verifier's x to that fixed, published set would cut that by about 24%, or about 17% with just the Legendre bit test, without changing what the benchmark measures.

We appreciate that changing the verifier mid-competition has a fairness cost, and the current all-integers contract is simple and well defined, so this is just an option for the hosts, perhaps for a future version. (Done with AI assistance.)

cc @quantumbaba @ryanbabbush, as co-authors of the stopping-power paper.
