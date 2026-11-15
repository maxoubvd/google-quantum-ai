# Predictions for the final top leaderboard score?

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752514 — copied 2026-11-15 (raw Markdown, unedited).

**Yan Zhou** · 2026-10-25 07:08 UTC · 1 votes

Curious what everyone thinks the final winning score will be.

30.0? 40.0? 50.0+?

Any guesses?


---

## Comments (8 threads)

Replies are indented with `>`, one level per reply depth.

**Mark Cooper** · 2026-10-25 17:31 UTC · 2 votes

I’d be surprised to see anything near 70.
The score is the geometric mean of three ratios, so 70 needs a product of 343,000, which is about 70× on every task.

Each benchmark has hard floors in its cost formula.

DQI’s 70 syndrome registers alone are 840 qubits, and Q enters squared.

Chemistry’s λ_eff and the fixed state-prep cost model put a floor under that task too, which realistically keeps it in single digits.
If chemistry tops out around 10×, the other two would each need about 185×.

For DQI that means the whole 35-error decode in about 160k Toffolis at the absolute minimum qubit count.

Today’s best public decoders use about 1M, and Gu et al.’s paper needs 440k at more than twice the qubits.

For invsqrt it means a cost around 11k, against roughly 130–140k for the best circuits now.

Scores will keep rising, but unless someone finds a loophole in the verifiers, I’d expect the winner somewhere in the 25s, maybe low 30s with real breakthroughs.

**D Chakraborty** · 2026-10-25 09:52 UTC · 2 votes

6️⃣9️⃣.4️⃣2️⃣0️⃣

**yash bhaskar** · 2026-10-25 11:42 UTC · 1 votes

My Prediction:

If gold Median is x

Gold : x - 0.01 ,
Silver : x - 0.1 ,
Bronze : x - 1

> **Yan Zhou** — **Topic author** · 2026-10-25 12:17 UTC · 0 votes
>
> Wow, curious why you think the gap would be so small.

> **Marko Rupnik** · 2026-10-26 10:54 UTC · 0 votes
>
> Seems the most realistic to me, previous golf competitions had a similar shape in the end. 3 months is a long time to saturate the top.

**Navneet** · 2026-10-25 07:24 UTC · -5 votes

Thank you for the Predictions for the final top @yanzhou06

**mikelou1** · 2026-10-25 07:19 UTC · -1 votes

My predictions

1st: 59.492
Gold: 45.204
Silver: 40.592
Bronze: 37.592

> **OverfitOracle** · 2026-10-25 08:35 UTC · 2 votes
>
> Wow what a precise prediction.😂.Hope you are not a time traveler.

**hengck23** · 2026-10-25 15:10 UTC · 0 votes

i have a different question. Assume we are not using quantum computing. Do you think AI agent can make classical circuits that run really fast?

 i think AI can exploit any cirucit as long as there is a verifier?

> **Jannes Stubbemann** · 2026-10-25 15:16 UTC · 0 votes
>
> Up to some limits for sure. Not every quantum circuit is difficult to simulate and some algorithms have been shown to be simulatable efficiently with classical compute with fairly natural assumptions. AI can probably drive the bar for quantum advantage higher. At the same time, it can drive it lower by reducing the resources for circuits that will still be inherently hard to run

**NT** · 2026-10-25 11:56 UTC · 0 votes

70 for prize,
60 for gold,
50 for silver,
40 for bronze

**OverfitOracle** · 2026-10-25 08:36 UTC · 0 votes

My prediction is about 55 for 1st and for gold 50 and and then 40 for silver and 30+ for bronze
