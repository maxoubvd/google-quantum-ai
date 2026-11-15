# Is resetting workspace with R allowed in place of uncomputation?

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755302 — copied 2026-11-15 (raw Markdown, unedited).

**Jyoti Negi** · 2026-11-03 06:21 UTC · 5 votes

A question for the hosts about circuit validity.
kickmix has the R instruction, which resets a qubit to |0> in the Z basis without storing an outcome. The verifiers simulate basis-state inputs only, so a circuit that computes garbage into scratch qubits and then clears them with R passes the "scratch returned to |0>" check, the same as a circuit that uncomputes them properly or clears them with HMR plus phase corrections.
On real hardware these are different. If the scratch value depends on the input, an R reset measures it and collapses the superposition, so the circuit no longer implements the intended phase oracle (for invsqrt, |x>|0> -> exp(2 pi i / sqrt(x)) |x>|0> on superpositions of x).
Could you clarify:
1. Are submissions expected to be correct on superposition inputs, or only on the basis states the verifier runs?
2. Is using R to clear input-dependent scratch allowed, or will such circuits be disqualified in the final review?
Thanks.


---

## Comments (4 threads)

Replies are indented with `>`, one level per reply depth.

**Satvik Verma** · 2026-11-05 03:36 UTC · 1 votes

Following up on this question: could the hosts confirm whether resetting scratch qubits with R (instead of uncomputing them) is allowed when the scratch value depends on the input, provided the official verifier passes? A ruling either way would help everyone stay within the rules. Thank you!

**Hassan Mohamed** · 2026-11-04 14:35 UTC · 1 votes

+1. A ruling before the 20 Nov milestone would help everyone plan, because the answer changes the cost scale a lot.

One concrete data point: in the public DQI decoders, the uncompute half of the circuit is about half of all Toffolis. A circuit that clears input-dependent workspace with R instead of uncomputing it (or using HMR plus phase fix-ups) would pass the published verifier, since it only runs basis-state inputs. On a superposition input, though, it would no longer implement the phase oracle. The same applies to invsqrt's uncompute stage.

Could the hosts confirm two things?
1. Must final submissions be correct on superposition inputs, i.e. implement the stated unitary or phase oracle, not just pass on basis states?
2. Will the final review check for input-dependent R (for example with a randomised-input or symbolic check)? Or is the published verifier the whole spec?

We're not alleging anything about any submission. We just want to know which of the two design spaces we should be optimising in. Thanks!

**Syed Asad Ali** · 2026-11-06 09:59 UTC · 0 votes

Data point on the R question: the verifier already catches input-dependent R resets. kickmix's R randomizes the phase when the qubit is |1> ("if the target qubit is in the |1> state, this operation randomizes the phase"), and each x is simulated 64 times with seeded outcomes. Test: a passing invsqrt circuit plus "CX q0 q61; R q61" at the end fails with "x=5: phase error 5.000e-01 turns"; the same circuit plus "CX q0 q61; CX q0 q61" passes. So clearing input-dependent scratch with R does not pass; HMR plus phase fix-ups (or real uncomputation) is still needed.

**David Pilný** · 2026-11-05 18:42 UTC · 0 votes

Adding a concrete case to this question, since a ruling before the 20 Nov milestone affects how everyone budgets the next weeks.

In all three tasks, the verifier runs computational-basis inputs only, so it cannot distinguish (a) clearing an input-dependent scratch value by uncomputation or HMR + phase correction from (b) clearing it with R. On superposition inputs, (b) measures the scratch value and destroys the phase oracle. The size of the difference is large: in today's public DQI decoders, the reverse (uncompute) passes are roughly half of the mean Toffoli count, and in the public invsqrt circuits roughly a quarter. So the answer changes the achievable cost scale by much more than any recent public improvement.

Could the hosts confirm:

1. Must submissions implement the stated oracle on superposition inputs (so R on input-dependent scratch is disallowed), or is passing the published verifier sufficient?

2. If superposition correctness is required, will the final evaluation check it (by review or an additional test), and will circuits that rely on R resets be removed from the final leaderboard? (This connects to the open question in 755427 about the final evaluation scope.)

I am not using R resets on input-dependent values myself and am asking only so that everyone optimizes against the same rules. Thanks!
