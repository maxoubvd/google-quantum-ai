# Grid search vs. exact trigonometry

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755525 — copied 2026-11-15 (raw Markdown, unedited).

**Pablo Diego Pereira** · 2026-11-03 21:07 UTC · 0 votes

Title: Grid search vs. exact trigonometry: how much FeMoCo L1 norm is still sitting on the 45-column table?
When you open the FeMoCo tensors in the DF-THC-SOS task, even after compressing the inner rank from 60 down to the 45-column compact-group table (which already drops width to 427 qubits and total Toffolis to 15,321), the natural way to trim block-encoding cost without touching the spectral errors is to sweep Givens angles across a fixed grid, but a discretized grid misses the exact spots where the remaining gains sit. Since the two-electron integral reconstruction only sees the augmented Gram product of each slice, multiplying on the right by any 45-by-45 orthogonal matrix preserves the spectral errors and Givens angles down to machine precision while changing the column L1 norms, yet that L1 objective has sharp kinks every time a matrix entry crosses zero. It is like rotating rigid crates inside a shipping container before strapping them down: the cargo delivered at the port is identical down to the last bolt, but the angle that uses the least strapping tape almost always sits right where a crate edge touches the wall, which a coarse angle grid steps right over.

Rather than brute-forcing a finer grid, we used the fact that between any two adjacent sign flips all 144 row signs stay frozen, turning the objective on each subinterval into a single smooth harmonic whose exact minimum across all 288 breakpoints vectorizes in seconds in NumPy. Combining that exact O(45) trigonometric and Cayley gauge polish (which pushed Lambda_2 on the compact table from 270.230986 down to 270.230012 with max Gram drift at 1.95e-14) with a 62-qubit inverse-square-root circuit at 97,921.54 cost, an 891-qubit cached-product DQI decoder averaging 396,874.2 Toffolis, and a bit-packed verifier across all 1,220 QROM addresses lifted our geometric mean from 22.99252 to 24.69664 on the public leaderboard (rank 94 out of 2,309 teams), breaking the tie just ahead of the 24.69663 cluster. I am not sure whether anyone will squeeze another rank reduction out of FeMoCo without blowing past the 20 mHa spectral cap before the deadline, so we shared the vectorized solver and verification plots in our [Google Quantum Compilation Starter Notebook](https://www.kaggle.com/code/pablodiegoo/google-quantum-compilation-starter-notebook) — how low have you managed to push Lambda_2 on the 45-column factorization while keeping exact Gram invariance?

![Exact piecewise-trigonometric gauge optimization and three-task compilation diagnostics](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F11792984%2Ff965985cc0f91218bb5a191aca79530f%2Fquantum_compilation_diagnostics.png?generation=1793837643870017&alt=media)


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**Pablo Diego Pereira** — **Topic author** · 2026-11-05 00:29 UTC · 0 votes

Updated... bit better numbers.
