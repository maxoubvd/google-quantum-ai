# Clarification on Apache-2.0 circuit reuse and winner reproducibility

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753316 — copied 2026-11-15 (raw Markdown, unedited).

**Michael O'Flaherty** · 2026-10-27 20:40 UTC · 2 votes

Hello organizers,

I have a few questions about combining independently developed work with publicly released, Apache-2.0-licensed quantum circuits in this competition.

1. Third-party circuit reuse: If a participant incorporates an unchanged circuit from another participant's publicly released Apache-2.0 work, preserves attribution, and contributes other components and integration work, is that permitted under the original-work warranty in Section C(14)?

2.  Reproducibility: If the original author published the executable circuit but not its generator, would providing the exact publicly available binary, source references, hashes, and complete code and instructions to assemble and verify the resulting submission satisfy Sections B(5), B(7), and C(6)? Or must the participant also provide the original author's generator?

3. LICENSE/NOTICE: Since the evaluation format specifies a ZIP containing the benchmark circuit files, where should the required Apache-2.0 license and attribution notices be provided? May they accompany the submission in a separate public repository or documentation, or must they be included inside the scoring ZIP?

The objective is to reuse publicly licensed work responsibly, credit its original authors, and ensure that both submissions and any potential winning solution comply with the rules.

Thank you for clarifying.
