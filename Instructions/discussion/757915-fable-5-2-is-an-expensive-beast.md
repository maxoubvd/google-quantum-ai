# Fable 5.2 is an expensive beast

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757915 — copied 2026-11-15 (raw Markdown, unedited).

**pjmathematician** · 2026-11-11 06:14 UTC · 9 votes

I tried the new fable on my current setup and its impressive and surprising how quickly it is able to adapt and perform better than the models i previously tested. 
But it uses a loooot of tokens, and shoots up the cost super quickly. Currently I am making fable just plan out stuff and letting kimi do most of the coding. But still need to figure out an even more efficient way.


---

## Comments (7 threads)

Replies are indented with `>`, one level per reply depth.

**hengck23** · 2026-11-11 06:46 UTC · 4 votes

Kaggle competition is not the the same as before. You need to make BIG investment to win.

It is now a fight on how to reduce cost.

> **Dungdv** · 2026-11-11 07:54 UTC · 4 votes
>
> Its always like that bro.
>
> Before, you still need GPUs to train, more GPUs mean you can do your experiement faster.
>
> it's about optimization

> > **Andre** · 2026-11-11 08:58 UTC · 1 votes
> >
> > You guys using GPUs? 😂 (poor me)

> **GrantPrater** · 2026-11-11 21:28 UTC · 2 votes
>
> You're talking to the current 2nd place. I think they know exactly what kind of investment is needed. It's also less an investment and more about constraining the search space as much as possible and meta optimizing everything. It doesn't matter if you use 100x the compute if you're representing the problem in slightly the wrong way that leads to a combinatorial explosion

**Asterios Panoras** · 2026-11-12 01:01 UTC · 1 votes

I tried Fable earlier and found it genuinely useful, especially for high-level stuff, but I ended up dropping it because the token cost was just too steep. For day-to-day work, Sonnet (medium effort), Opus, and Kimi already handle things really well. That said, bringing Fable back strictly for analyzing past runs, spotting better circuit moves, or planning major architectural shifts sounds like a really smart way to get its value without burning through credits. Might give that workflow a shot.

**yuanzhe zhou** · 2026-11-11 16:48 UTC · 1 votes

Kimi3 is good. I tried to use gpt 6 and a lot of gpu, but it does not seem to be helpful ... gpt 6 is stupid or I asked stupid questions to it  👀. It kept saying things like "sorry I did not make any improvement and I will work hard".

> **Rustam Bazarbayev** · 2026-11-11 16:53 UTC · 0 votes
>
> What are you plan Plus or Pro?

> > **yuanzhe zhou** · 2026-11-11 16:59 UTC · 0 votes
> >
> > pro with codex

> > **Rustam Bazarbayev** · 2026-11-11 17:27 UTC · 0 votes
> >
> > Pro with Codex is not enough for me

> **pjmathematician** — **Topic author** · 2026-11-11 16:53 UTC · 2 votes
>
> Gpu hasnt been very helpful for me as well for this competition

> > **yuanzhe zhou** · 2026-11-11 16:59 UTC · 0 votes
> >
> > It is able to transform everything to GPU, searching very fast, but no improvement

> **hengck23** · 2026-11-11 17:26 UTC · 3 votes
>
> I haven't used kimi for this competition yet. But in other competition kimi3 is cheaper and is about 95% as accurate as Claude and ChatGPT.
>
> Here, my experience with the two:
> - ChatGPT: mr yesman. Says everything can be done. But in the end, it cannot be done, but he is not honest. E.g. I ask him to read a public notebook and reimplement the solution (eg lattice+LLL quantization). He ends up copying the byte-by-byte circuit solution.
> - Calude: mr angryman. likes to argue with you things cannot be done or his solution is better. E.g. some task above, he kept insisting current progressive quantization pq is better and refuses to implement lattice+LLL. Do a lot of tests and measurements to convince me not to do it. Even if you force him to do it, he will give "sloppy or half-baked work". But onec you prove that it can gives better results (e.g. experiment results or forum post) he suddenly give much better quality work and try more.
>
> ----
>
> But the problem may not be the agent. I feel that i have not used them correctly.  Because in NeuroGolf, when I know the code and the problem, they can perform much better when i give them better prompt and template.

**GrantPrater** · 2026-11-11 21:31 UTC · 0 votes

Yeah fable is not intended for anything other than planning and maybe orchestration and review, unless you're just trying to burn money

**xz** · 2026-11-11 17:22 UTC · 0 votes

Do you use max and ultracode?

**Adarsh** · 2026-11-11 14:41 UTC · 0 votes

![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F14670471%2F47059033c17fbc47e26532272170c48c%2Ffable%205.2.png?generation=1794408045401917&alt=media)

My Fable has been cooking for 40 mins xD

**baithunter635** · 2026-11-11 14:19 UTC · 0 votes

Ah the Plutocrats strikes again.
