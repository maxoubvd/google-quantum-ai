# Some early thoughts on setting up a two stage search loop for kickmix

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753741 — copied 2026-11-15 (raw Markdown, unedited).

**Raitul Hasan Priyo** · 2026-10-29 06:28 UTC · 2 votes

Hey everyone! While exploring the kickmix baselines, I've been experimenting with a two-stage search loop: first optimizing ancilla qubit scheduling especially for dqi where the Q^2 penalty dominates, and then running a local peephole search on subcircuits to squeeze down Toffoli and T-gate counts. Using fast local kickmix fuzz testing inside a sliding-window mutation loop seems much more reliable than rewriting entire .kmx files at once. For those pushing up the leaderboard, are you seeing bigger gains from early ancilla uncomputation or algebraic gate rewrites?
