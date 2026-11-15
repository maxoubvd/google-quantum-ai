# ML method: Quantum Circuit Optimization with AlphaTensor

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754466 — copied 2026-11-15 (raw Markdown, unedited).

**hengck23** · 2026-10-31 07:18 UTC · 3 votes

paper: Quantum Circuit Optimization with AlphaTensor  
https://arxiv.org/abs/2402.14396  

The reference code is available in the repo. It is probably useful for the LUT table circuit optimization in dfthc, which my estimated gain is 2 to 4% (e.g. 10x to 10.4). Once I go into it, the effort spent will be quite large. Any one has suggestion on its performance? is it worth the effort?


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**Zero** · 2026-10-31 15:06 UTC · 1 votes

Google DeepMind hasn't released the trained model weights. Wouldn't training AlphaTensor-Quantum from scratch be prohibitively expensive?

> **hengck23** — **Topic author** · 2026-10-31 15:27 UTC · 0 votes
>
> Yes. That is why i am asking if anyone has experience on it etc. we have 3 months to go and this is one of the marathon tricks if it works ( other will not be able to duplicate ) but i need to prove it will wirk by doing some random walk on the circuit etc
