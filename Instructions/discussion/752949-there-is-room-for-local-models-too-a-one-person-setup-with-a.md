# There is room for local models too: a one-person setup with a limited AI budget

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752949 — copied 2026-11-15 (raw Markdown, unedited).

**Valerio Dolci** · 2026-10-26 19:51 UTC · 6 votes

A thread earlier today asked whether the leaderboard is just about who spends the most on API tokens. My experience says not necessarily, so here is how I work: one person plus Claude, with a limited AI budget.

- **Brain (one subscription).** A single Claude subscription does the "head" work: reading the public write-ups, designing exact algebraic changes, and reviewing what has been tried.
- **Brute force (local, free).** An open-weight model, Qwen3.8-Flash-Next, runs on two consumer GPUs (2× RTX 5070 Ti) through my own open-source MoE offload engine, [capy-flash](https://github.com/ValerioDolci/capy-flash). It runs inside a custom harness built for this competition: build, verify, keep or discard, write a journal entry, try the next idea, session after session through the night. Every "best" it reports is re-checked with the untouched official verifier before I trust it.

Most of my chemistry gains came from the local model, one exact, measured change at a time. Tonight, in its first two sessions on a new table, it beat the best public chemistry circuit by about 5 %.

The edge in this competition is how fast your build → verify loop runs, and that loop does not need frontier-model tokens at every step. A local model inside a good harness can run it all night for the cost of electricity.
