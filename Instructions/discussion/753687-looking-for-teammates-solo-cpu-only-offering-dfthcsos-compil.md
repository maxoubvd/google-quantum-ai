# Looking for teammates (solo, CPU-only): offering dfthcsos compile-and-verify tooling, seeking optimization/chemistry skills

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753687 — copied 2026-11-15 (raw Markdown, unedited).

**Wan Xiaolong** · 2026-10-29 01:17 UTC · 3 votes

Hi all! I'm a solo participant working with AI assistants on a 6-core CPU laptop (no GPU), currently at 22.63 on the public LB. I have **no unique score gain yet**, so I'm not offering a better circuit. What I can bring is engineering support, mostly on the dfthcsos task:

- an angle-table → ANGLE_k planner (reusing SomeJoint's public emitter). On SomeJoint's v1 table its official cost is within +0.16% of theirs (4 salts PASS); on v2 it is still +1.0% (4,262 vs 4,218 Toffoli per pass);
- a fast dfthcsos model with analytic gradients that matches the official evaluator to ~1e-13, an exact incremental error_1 evaluator, and a rank proxy for lookup cost (Spearman 0.83 against 28 real compiles);
- a verify + multi-salt pipeline, daily public-notebook harvesting, and logged negative results: 7-bit two-body QAT by local search plateaus at error_2 ≈ 0.11; STE, soft-to-hard and progressive freezing all failed on my side.

I'm looking for people with continuous-optimization / QAT or chemistry-factorization experience (a GPU is a plus) who could use compile-and-verify support, or anyone strong on dqi / invsqrt who wants a second pair of hands on the chemistry side. Happy to swap code and notes first and to coordinate submissions (I have used only 1 submission so far, so a merge is well within the limit).

Reply here or send a team invitation to **wanxiaolong** if interested.


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**Wan Xiaolong** — **Topic author** · 2026-10-29 06:07 UTC · 0 votes

Update (Oct 29) — a verified number, not a reproduction kit yet.

Starting from Dmitrii Gluzdov's public v30 chemistry file (angle words = SomeJoint's xor2 table), I refit only the two-body weight arrays W / H2_I under the official constraints; U, angle_words, tree_edges, brot, H1_I and Hsym_lower stay byte-identical, so heliosli's 15 public angle circuits are reused as-is. Official verifier, default seed + salts 11/22/33 (all identical): lambda_eff 48.63243 -> 48.48773, T 16,296, Q 430, cost 340,781,040 -> 339,767,112 (-0.298 %, about +0.10 % on the geometric mean). Input sha256 26974888..., output df888cf2....

Reproducibility: the same procedure on three other public tables with different angle words (SomeJoint v1 R8/B143/C60, heliosli R8/B143/C62, Michael Rustad's R10/C46 four-word table) gives -0.15 to -0.24 % each, verified with their own public circuits (default + 3 salts PASS). Negative result: with words fixed, no rank-59/58 table reached 20 mHa (error_2 floors 0.0201 / 0.0206).

Method write-up, code and the refitted file are not public yet. Still looking for teammates with continuous-optimization / QAT experience on the chemistry factorization; I bring the compile-and-verify side. Thanks to Dmitrii Gluzdov, SomeJoint, heliosli and the authors they credit.
