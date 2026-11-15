# Clarification request — Chemistry ANGLE XOR contract and reachable input states

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757522 — copied 2026-11-15 (raw Markdown, unedited).

**Michael O'Flaherty** · 2026-11-09 22:38 UTC · 0 votes

Hello,

We’re reviewing the quantum Chemistry ANGLE subroutine specification and would appreciate clarification on its required input domain.

Equation (72) appears to define an XOR transformation for arbitrary incoming word-register states:

|a⟩|z⟩ → |a⟩|z ⊕ Δₖ(a)⟩

Our independent analysis of an accepted ANGLE packet found that it behaves correctly on the forward/reverse states exercised by the official verifier, but does not preserve the specified transformation for certain other initial word-register values. In particular, one declared word bit is reset, producing a different output and measurement-dependent phase.

We would like to clarify whether compliance requires the unrestricted XOR transformation for every incoming word state, or whether the competition explicitly permits circuits that are correct only on the reachable states of the full Chemistry block-encoding sequence.

We don’t want to infer a weaker mathematical contract merely because the current verifier accepts the circuit.

We can provide the exact circuit identity, reproducible counterexample and independent analysis if helpful.

Thank you for clarifying the intended requirements.


---

## Comments (4 threads)

Replies are indented with `>`, one level per reply depth.

**World Champ** · 2026-11-12 03:23 UTC · 0 votes

Following up on this — the question is now deciding work across the field, not just one packet. With the Nov 20 milestone snapshot nine days out, teams are choosing which chemistry loaders to build January packages on, and readings (A) and (B) point in opposite directions.

Hosts — even a one-line answer would help: for prize eligibility, is (A) strict Eq. (72) for arbitrary word-register states, or (B) reachable-states-only as the verifier checks, the rule? Thank you.

**ABHISHEK KUMAR** · 2026-11-11 13:16 UTC · 0 votes

+1, we'd also like a ruling on this. It decides which chemistry ANGLE packets are eligible, both on the leaderboard and for the single-task prize, so an answer now would save everyone work in the wrong direction.

To make it quick to answer, here are the two readings as we understand them:

**(A) Strict Eq. (72).** Every ANGLE packet must implement |a⟩|z⟩ → |a⟩|z ⊕ Δₖ(a)⟩ for every incoming word-register state z, superpositions included. That means no input-dependent scratch cleared by measurement or reset, and no extra phase.

**(B) Reachable states only.** A packet only has to be correct on the states the full block-encoding sequence actually feeds it, which is what the verifier's forward/reverse trajectories check.

Three short questions:

1. Is (A) or (B) the rule?
2. Does the same rule apply to the leaderboard and to prize eligibility?
3. If (A): will prize review test packets on arbitrary word inputs, and will accepted (B)-only packets be rescored or ruled ineligible?

Thank you!

**Michael O'Flaherty** — **Topic author** · 2026-11-09 22:41 UTC · 0 votes

Example: In one accepted Chemistry ANGLE packet, at selector address 0, the required five-word XOR delta is [10, 7, 17, 1, 2]. Starting from [0, 0, 0, 0, 0], the circuit correctly produces [10, 7, 17, 1, 2]. However, starting from [0, 0, 0, 0, 128], the unrestricted XOR specification requires [10, 7, 17, 1, 130], while the actual circuit produces [10, 7, 17, 1, 2]. Independent analysis also identified a measurement-dependent phase difference. The circuit passes the verifier's reachable-state trajectories. Are such reachable-state-specific implementations permitted, or must every ANGLE packet satisfy Equation (72) for arbitrary initial word-register states?

**(deleted user)** · 2026-11-15 08:04 UTC · 0 votes

*(deleted)*
