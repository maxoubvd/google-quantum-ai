# Two exact finite-field multiplication complexities — reproducible proofs

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756476 — copied 2026-11-15 (raw Markdown, unedited).

**田村** · 2026-11-06 13:55 UTC · 11 votes

While working on this competition, we obtained computer-assisted proofs of two exact values:

- $\mu_4(5)=\mu_4^{\mathrm{sym}}(5)=11$: no rank-10 algorithm exists, symmetric or asymmetric.
- $\mu_4^{\mathrm{sym}}(6)=14$: no symmetric rank-13 algorithm exists.

The [public notebook](https://www.kaggle.com/code/xolotlmictlan/two-exact-finite-field-multiplication-complexities) includes the
proof reductions, complete exhaustive-search code, matching constructions, and verified outputs. All checks passed on Kaggle.

The degree-six result is symmetric-only; asymmetric rank 13 is not excluded. Independent verification and corrections are welcome!


---

## Comments (2 threads)

Replies are indented with `>`, one level per reply depth.

**pjmathematician** · 2026-11-06 16:01 UTC · 3 votes

Nice work, and thanks for sharing the code. I also independently arrived at both of your values, and can add the degree-six asymmetric case: μ_4(6) = 14 exactly, asymmetric algorithms included (rank 12 is excluded too), so no 13-term GF(4)-bilinear algorithm for GF(4096) exists at all.

Sketch: since no [13,6,7]_4 code exists (shortening gives a [10,3,7]_4, impossible as m_3(2,4) = 9), any rank-13 scheme has a tight slice, which after scaling is a frame: the projections of a GF(4)-basis of GF(4^6). V = W + span(frame) then has dimension 11 and holds no other rank-one, so the remaining 7 terms must lie on one projective line of M/V. Enumerating all frames up to Frobenius and Singer symmetry (4,485 configurations, 11.28M frames, 2.38e9 candidate lines) gives 0 hits, which also excludes rank 12; Griesmer covers rank ≤ 11 and CCMA gives the upper bound 14. I confirmed it with a brute force over all 63.4M bases, and a separate completion that recovers 763/763 planted rank-13 schemes yet finds none on the real tensor.

For μ_4(5) = 11 the same tight-frame method applies (no [10,5,6]_4 left code exists): all 15,597,604 tight-frame orbits give 0 hits. For μ_4^sym(6) = 14 I enumerated all C(455,4) windows of Veronese points, with 0 solutions.

> **田村** — **Topic author** · 2026-11-07 02:14 UTC · 1 votes
>
> Thanks for sharing this! The asymmetric degree-six result closes the case we had left open. Really appreciate the independent checks and the detailed sketch.

**Dyst-0** · 2026-11-07 17:16 UTC · 1 votes

Thank you both. We reached the same values independently and have **formal Lean 4 + Mathlib proofs** of both. They include the asymmetric degree-six case, for every extension L/K with |K| = 4:

- μ₄(6) = 14, as `Mu4.Final.mu4_6_general`: `bilinRank K (LinearMap.mul K L) = 14` when [L : K] = 6.
- μ₄(5) = 11, as `Mu4.N5.mu4_5_general`: the same with [L : K] = 5 and rank 11.

Here `bilinRank` is the least r with f(u, v) = Σᵢ aᵢ(u) bᵢ(v) cᵢ, for linear forms aᵢ, bᵢ and vectors cᵢ. The upper bounds are explicit CRT schemes checked by `decide`. The lower bounds go:
- substitution to a hyperplane;
- a pencil structure;
- a normal form, with kernel-checked coverage tables;
- an exhaustive search (about 2.4M candidate triangles for n = 6, none survives), whose completeness is proved in Lean.

**Trust base:** `#print axioms` gives `propext`, `Classical.choice` and `Quot.sound`, plus `Lean.ofReduceBool` and `Lean.trustCompiler`. The last two come only from `native_decide` in three finite searches that are too large for the kernel. Everything else is kernel-checked, and there is no `sorry` or `axiom` declaration.

**Archive:** sha256 `58e138356aba88aa5af1763aa1aa0d6a506cab0c9024b812d1f6498745accd51`, posted on 11-06: https://x.com/mathbbA/status/2118814683062730872. The repo has the archive, its unpacked sources and the proof as precise lemmas: https://github.com/Dyst-0/mu4-lean-proofs

Lean v4.28.0 and Mathlib v4.28.0. To build: `lake exe cache get && lake build && lake env lean Axioms.lean`.

Credit to 田村 for the first public proofs, and to pjmathematician for the asymmetric sketch. Happy to compare notes.

> **田村** — **Topic author** · 2026-11-08 11:28 UTC · 0 votes
>
> Thank you for the formalization and the credit. Great to have machine-checked proofs for both values.
