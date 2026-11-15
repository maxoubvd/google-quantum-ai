# DQI on real hw 

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755738 — copied 2026-11-15 (raw Markdown, unedited).

**Peter Babulik** · 2026-11-04 13:40 UTC · 1 votes

All experiments were executed on IBM Quantum's 156-qubit Heron processor.

FanoDQI is an open-source quantum compilation and experimental verification framework that downscales and accelerates Decoded Quantum Interferometry (DQI)—the core algorithm of the Google Quantum AI Verifiable Circuit Compilation Challenge—onto physical superconducting hardware.
By mapping the finite field structure 
 directly onto the Fano projective plane, this framework achieves:

Sub-second Hardware Transpilation: Exploiting the shared bipartite graph symmetry between the 8D hypercube (Okteract) and IBM's Heavy-Hex lattice eliminates heuristic routing latency.
Super-Resolution Metrological Filtration: Syndrome-based post-selection and error-filtering during the DQI recombination stage deliver a +22.98 Signal-to-Noise Ratio (SNR) boost, surpassing the theoretical unpurified peak limit (20.63 vs 18.75).
Full-Chip Scalability across 154 Physical Qubits: Concurrent execution of 22 parallel 7-qubit Fano-DQI blocks simultaneously on IBM Quantum's 156-qubit Heron processor (ibm_fez), demonstrating macroscopic quantum coherence at 98.7 processor capacity.

[](https://github.com/peterbabulik/FanoDQI)


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**Peter Babulik** — **Topic author** · 2026-11-04 13:41 UTC · 0 votes

https://github.com/peterbabulik/FanoDQI
