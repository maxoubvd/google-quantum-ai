# SymQ: checking circuit rewrites, ancillas and measurement corrections — feedback welcome

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754815 — copied 2026-11-15 (raw Markdown, unedited).

**yoch** · 2026-11-01 20:32 UTC · 0 votes

Hi everyone,

When optimizing a circuit, how do you check that a cheaper replacement preserves the outputs, restores borrowed wires, and respects the assumptions where it is inserted?

I’ve been developing **[SymQ](https://github.com/yoch/symq)**, an open-source Python toolkit for checking these properties. I’d like to share it with people working on this competition and get feedback from real use.

**What it can help with today**

• Check equivalence of reversible X/CX/CCX/MCX circuits using symbolic Boolean representations or SAT.
• Find and replay counterexamples when supported checks detect a difference.
• State input assumptions and restoration requirements explicitly, including dirty ancillas whose initial values are unknown.
• Verify a local substitution in its circuit context, rather than checking two isolated blocks and assuming the surrounding conditions hold.
• Check small, restricted measurement-assisted rewrites, including phase corrections that classical output tests can miss.

The CLI supports locating an exact block, checking a replacement, saving a report and replaying the check. Locating a block is only a suggestion; it does not prove that the replacement is valid.

**A small example to try**

The [public demo](https://github.com/yoch/symq/blob/main/docs/demo.md) computes a temporary AND into a clean ancilla, then uses X-basis measurement/reset and a conditional CZ correction.

With the correction, the declared quantum channel matches the reference. Remove it, and the classical output bits can still look right while quantum coherence is lost. SymQ detects the difference and replays a witness.

The demo runs with the pure-Python core, without a native compiler or quantum service. The [contextual rewrite guide](https://github.com/yoch/symq/blob/main/docs/guide/rewrite.md) also provides workflows using source-pinned public KMB and small Qualtran examples.

**Limits worth knowing**

SymQ is experimental research software, and its API is evolving. Classical equivalence and coherent-channel equivalence are different guarantees. Assumptions and supported operations matter; a timeout or unsupported case is not a proof.

The coherent workflow is deliberately bounded. This does not replace the official competition verifier or certify an entire large competition circuit. Its gate-cost estimates are separate from Kaggle’s cost model.

**Help make it useful on real circuits**

I’d welcome:
• Small circuit blocks and proposed optimizations to test.
• Bug reports, counterexamples, timeouts and confusing diagnostics.
• Cases involving dirty ancillas, compute/uncompute or measurement corrections.
• Contributions to importers, examples, independent checks and documentation.

For a bug report, please include the commit, a minimal circuit pair, intended assumptions, command and observed result. Please share only circuits you are allowed to publish, with their provenance and license.

[Issues](https://github.com/yoch/symq/issues) and [contributions](https://github.com/yoch/symq/blob/main/CONTRIBUTING.md) are welcome. If you try it, tell me what helped, what failed, and what you would need to use it in your own optimization workflow.
