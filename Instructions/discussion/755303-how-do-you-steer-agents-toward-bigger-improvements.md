# How do you steer agents toward bigger improvements?

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755303 — copied 2026-11-15 (raw Markdown, unedited).

**uctpexx123** · 2026-11-03 06:41 UTC · 1 votes

I’ve been experimenting with an RSI loop: propose → implement → verify → revise.

Two problems keep appearing:
- Agents converge on small optimizations around the current baseline.
- One failed formulation leads them to abandon an entire direction, even though other formulations remain untested.
Adding more agents often produces more variations of the same ideas.

For those finding substantial structural improvements, how do your agents discover promising directions and organize exploration? What helps agent distinguish “this implementation failed” from “this idea is unlikely to work”?

Practical examples of your discovery process would be greatly appreciated!


---

## Comments (3 threads)

Replies are indented with `>`, one level per reply depth.

**Agniv Banerjee** · 2026-11-03 18:29 UTC · 3 votes

I don't think is RSI loop  propose → implement → verify → revise.

You need some kind of DAG structure, some kind of map, with diffrent directions so the probability of stucking is lesser.

Something like this https://x.com/KingBootoshi/status/2105011580870242328

**ranmiao** · 2026-11-04 10:26 UTC · 0 votes

I am attempting to combine coding agents with population evolution.

**xz** · 2026-11-03 18:55 UTC · 0 votes

"Claude, make no mistake."
