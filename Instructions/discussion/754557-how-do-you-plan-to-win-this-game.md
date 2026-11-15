# How do you plan to win this game? 

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754557 — copied 2026-11-15 (raw Markdown, unedited).

**yoyobar** · 2026-10-31 14:18 UTC · 4 votes

I think using AI smartly and efficiently is the key to this challenge.

![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F8293010%2F8304785d7de67ae832997212a0e89489%2F692be7ae-6e81-47a5-9a04-1dd4b6563b7b.png?generation=1793456124108850&alt=media)


or 

![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F8293010%2F4697d125957b37cef3931aabfabcea08%2Fmoney-is-all-you-need-deep-learning-ai.jpg?generation=1793456285716830&alt=media)


---

## Comments (6 threads)

Replies are indented with `>`, one level per reply depth.

**pjmathematician** · 2026-10-31 17:16 UTC · 3 votes

bulk experiment with ai -> tune it yourself -> repeat

> **Jack** · 2026-11-01 15:49 UTC · 1 votes
>
> Wym tune yourself? I'll talk things through with my agents to help guide direction but I'm not tuning anything myself

**Tom** · 2026-11-01 11:15 UTC · 2 votes

One-shot winning

> **Sangram Patil** · 2026-11-01 14:00 UTC · 0 votes
>
> interesting

**Zero** · 2026-10-31 15:30 UTC · 1 votes

I'm planning to try a math-based reinforcement learning approach. LLM-driven research has reached a point of diminishing returns for me, so I'm looking for a fundamentally new way to explore the search space.

AlphaEvolve and CodeEvolve-style methods are also difficult to apply efficiently here, since relying on Google's full official validator for every iteration is too slow and cumbersome.

> **manwithacat** · 2026-10-31 16:09 UTC · 1 votes
>
> there's a faster rust-based version of the evaluator. Cut memory usage down to about 10% on some problems . Check the RezOps briefing for link.

> > **Zero** · 2026-11-01 04:03 UTC · 0 votes
> >
> > Thanks! I've checked out your accelerated validator, and it looks great! This could improve things quite a lot.

> **Mohit** · 2026-10-31 17:20 UTC · 0 votes
>
> Yo i had the same in mind

> **hengck23** · 2026-10-31 17:58 UTC · 0 votes
>
> Problem with training method is ground truth or exploration. I wonder how google did it

> > **Agniv Banerjee** · 2026-10-31 18:18 UTC · 0 votes
> >
> > A lot of compute. I tried something like this https://www.anthropic.com/research/designing-efficient-quantum-circuits but it has problem with a cost and scope even with rust speedup. Maybe if try on opus 5.5 it will be better, i mostly tried sol or terra.

> **(deleted user)** · 2026-11-01 04:12 UTC · 0 votes
>
> *(deleted)*

**D Chakraborty** · 2026-11-01 04:12 UTC · -2 votes

I had major breakthrough after switching to RL with the LLM's confidence as the reward function. No more validator bottleneck, unlimited exploration and every iteration is more confident than last time. Rank 1 is now just a matter of compute.

(i hope someone's agent picks this up and breaks the "plateau")

**WOOSUNG YOON** · 2026-11-01 00:03 UTC · 0 votes

Well... as long as I post my score on the leaderboard and never look at the final standings,

shouldn't I technically be winning by probability?

**hengck23** · 2026-10-31 17:59 UTC · 0 votes

It seems that there are no better mathematical formulation? So the next 2.7 months will be wiring optimization?

> **Agniv Banerjee** · 2026-10-31 18:22 UTC · 1 votes
>
> I think with releases like fable 5.5 and astra 6.2 we will definitely will have some better math formulation. According to rumors release is possibly   in next week

> **sebewe** · 2026-11-01 12:14 UTC · 1 votes
>
> Openai recently broke the exponential time hypothesis, foundational math as we know it is subject to change. Throw enough tokens at it, with the correct approach, and everything we know to be certain is subject to change.
