# dfthcsos: why λ_eff can't go much below ~46.8 

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753633 — copied 2026-11-15 (raw Markdown, unedited).

**Lynn Langit** · 2026-10-28 20:13 UTC · 5 votes

λ_eff hasn't moved much in public work (~48.7) while everything else has. 

Alex Chilton's nulls thread (https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753185) notes that every measured fit sits at Λ ≈ 290–308, and mrmeows (https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752757) showed that refits can't improve λ_eff by more than ~0.3% on a fixed factorization. Here's a complementary bound that holds for any factorization: with the official verifier's fixed symmetry shift, λ_eff can't go below ~46.8, about 4% under today's public value. The Λ > 49 floor in the task statement never binds.

**Setup** (notation as in the task paper and dfthc_math.py). The verifier fixes S and β₁. With w_B = W + s, where s = Σ_b w_b,

- Λ = λ_h1 + λ_h2, with λ_h1 = Σ|eig(h₁′)| and λ_h2 = ¼ Σ_{ρ,c} (Σ_b|w_b| + |W + s|)²
- E_gap = const + (λ_h1 − tr h₁′) + ½ ΣW², with const = E_upper − e_core + η β₁ (fixed)

**Step 1** (the gap floor forces the offsets). E_gap > 3.914898 requires 2·Σ(negative eigenvalues of h₁′) + ½ΣW² ≥ D := 3.914898 − const. For FeMoCo76o113e, const = −559.132, so D = 563.047 Ha.

**Step 2** (the triangle inequality). Σ_b|w_b| ≥ |s| and |s| + |W + s| ≥ |W|, so λ_h2 ≥ ¼ΣW².

**Step 3** (combine). Λ ≥ Σ(positive eigenvalues) + Σ(negative eigenvalues) + ¼ΣW² ≥ Σ(positive eigenvalues) + D/2 ≥ D/2 = 281.52 Ha. So at the gap floor, λ_eff = √(E_gap(2Λ − E_gap)) ≥ 46.79.

**Where the current public build stands**((8, 143, 60) 8-bit factorization by the SomeJoint team, with Dmitrii Gluzdov's gauge; official verifier formulas). I checked these numbers against the official verifier:

|   | Value |
| --- | --- |
|  Λ  | 304.82 (λ_h1 12.81 + λ_h2 292.01) |
|  ¼ΣW²  | 270.27 |
|  Bound D/2 + Σ(positive eigenvalues)  | 283.08 |
|  λ_eff now / lower bound  | 48.70 / 46.79–46.92 |


All of the remaining gap (21.7 Ha of Λ) is the triangle-inequality slack in λ_h2. Closing it would require every w_b in a (ρ, c) column to share one sign, with W aligned against s, while still fitting the two-body tensor within 20 mHa. Even closing it completely gains at most 3.9% on λ_eff, about +1.3% on the score. In practice the offsets W are already optimal: re-optimizing them alone with the floor enforced changes λ_eff by under 10⁻⁶ %.

**Takeaway:** chemistry gains have to come from T × Q (fewer addresses, narrower words, fewer copies), not from λ_eff.

Happy to be corrected if I've missed a degree of freedom. The derivation only uses the official verifier's definitions.


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**Lynn Langit** — **Topic author** · 2026-10-30 00:52 UTC · 0 votes

Update (Oct 29): I re-checked this against the current best public chemistry build: Dmitrii Gluzdov’s v37 loaders on Atakan Aldemir’s word table, λ_eff 48.459. The bound itself doesn’t change, because it holds for any factorization. Updated numbers:

![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F311896%2Fe2776f6402f810f7532cf74ceef12058%2FScreenshot%202026-10-29%20at%207.52.36PM.png?generation=1793321574842384&alt=media)

One correction for the new build: re-optimizing the offsets W alone, with the floor enforced, does find a little. λ_eff drops by 0.014%, worth about +0.005% on the score. That is negligible, but it is not the “under 10⁻⁶ %” I measured on the earlier build.
