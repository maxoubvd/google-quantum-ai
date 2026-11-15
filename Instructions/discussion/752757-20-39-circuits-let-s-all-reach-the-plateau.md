# 20.39 circuits, let's all reach the plateau!

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752757 — copied 2026-11-15 (raw Markdown, unedited).

**mrmeows** · 2026-10-26 07:28 UTC · 21 votes

I've published my current best submission, together with the methods and the dead ends that led to it.

Notebook (reproduces the 14.49721 submission): https://www.kaggle.com/code/mrmeows/qcc-14-50-everything-i-learned-open-circuits

Dataset (circuits, generators, audits; Apache-2.0): https://www.kaggle.com/datasets/mrmeows/qcc-1375-open-circuits

  ***Why.*** There are three months left. I think this competition is most interesting if it's decided by new circuit architectures rather than by assembling public components. Releasing a strong, fully verified baseline early removes the easy gains for everyone, including me. I'd like the remaining progress to come from genuinely new ideas.

What's included.

invsqrt, 68 qubits, 17.27×. Storage and uncomputation changes on the public 78-qubit lineage. The polynomial and rotation angles are unchanged, and every identity is checked on all 3,133,452 inputs.

dfthcsos, 6.435×. 8-bit angle words take the state-preparation workspace from 607 to about 300 qubits (total Q 742, 448), with accuracy recovered through discrete refits and a larger dictionary.

dqi. Dmitrii Gluzdov's round12 decoder, unchanged and credited.

The methods, including a bound showing that refits can't improve λ_eff by more than 0.3% on a fixed factorization, decoder width lower bounds for BM and extended-Euclid, and two Bézoutian identities.

Negative results with measurements: 7- and 9-bit words, paired rounding, Newton and cubic invsqrt variants, and more, so nobody has to repeat them.

Practical notes on the verifiers: Windows runs them serially, kickmix allocates hidden ancillas, and short dqi runs inflate R.

All component costs come from the unmodified official verifiers at full shot counts, and the hashes are pinned. Credits to busyaprime, Phuc Thai (DePie), Dmitrii Gluzdov, marsyas, Agniv Banerjee (SomeJoint), ViberTeam, manwithacat, Ahmed Berat Özer, marculera and Michael Rustad, whose public work this builds on.

If you find something better, I hope you'll publish it too. Questions and corrections are welcome.

EDIT: 20.39 circuits now up! https://www.kaggle.com/code/mrmeows/20-39-circuits-let-s-all-reach-the-plateau/notebook

EDIT: 15.75 circuits now up

EDIT: 15.05 circuits now up!


---

## Comments (5 threads)

Replies are indented with `>`, one level per reply depth.

**mrmeows** — **Topic author** · 2026-10-26 08:56 UTC · 3 votes

Just shared my 15.05 circuits! https://www.kaggle.com/code/mrmeows/15-05-circuits-let-s-all-reach-the-plateau

**WOOSUNG YOON** · 2026-10-26 08:01 UTC · 3 votes

Thanks for sharing!

I’ve been diving into some books on quantum computing recently, and I get the feeling that the current theories in this space aren't fully set in stone yet. It seems like we just don't have enough experimental data or solid evidence to back up a single, unified theory.

Right now, most textbooks still rely on classical, Turing-machine-based models to explain everything. But every time I read them, I find myself questioning if building quantum gates this way is even the right approach. It honestly feels like the current theories are a bit outdated, and it makes me wonder if we need to rethink the whole foundation from scratch.

Instead of trying to control every single internal state using classical logic, what if we treated it like a machine learning system? A sort of black-box model where we don't necessarily manage the exact internal mechanics, but focus on training the system to deliver the right outputs based on the inputs. What are your thoughts on this? 😆

> **mrmeows** — **Topic author** · 2026-10-26 08:16 UTC · 3 votes
>
> Thanks, glad it helped! I'd push back a bit on the premise though. The physics itself is really well tested. What's still hard is the engineering of building big reliable machines. The gate model also isn't a classical idea at heart. It's just one of several equivalent ways to describe quantum computation, and people use it because it's the easiest to reason about.
>
> Your black box idea is actually pretty close to things people already do. Variational algorithms like VQE and QAOA train circuit parameters with a classical optimizer, and on real hardware the control pulses are often learned rather than designed by hand.
>
> The problem in this competition is that the verifier checks every single legal input and wants the exact phase, so something trained to be right most of the time would just fail. Where learning and search do help is in finding the circuits. A lot of my gains came from automated search over small subcircuits, and then I proved the result exact before using it. Search to find, proof to check. That combo feels like the most promising direction to me 😄

> **hengck23** · 2026-10-26 11:38 UTC · 2 votes
>
> Optimization requires experience. I think some will experience in c code, asm code or bit hacking. But very few have experience in quantum circuits.  The inv sqrt is something i did in the past but the solution in the public notebook are so unfamiliar to me.

**mrmeows** — **Topic author** · 2026-10-27 02:33 UTC · 2 votes

Added 21

Taking a break,

Will continue in few days. Hope I helped someone 😮‍💨

ok goodnight

> **hengck23** · 2026-10-27 02:45 UTC · 3 votes
>
> thanks. Will take a look and make a tutorial on how to start and progress from 1 -> 5 -> 10 -> 20 based on your work.
> Once it is saturated, the real fun begins, where only ai agent can further improve

**hengck23** · 2026-10-26 11:31 UTC · 2 votes

That is king. I have a feeling that it will be top10 against the rest

The top10 will have some secret source and the rest will be mixing components from each other

> **shanzhong8** · 2026-10-26 13:03 UTC · 1 votes
>
> For competition like this, maybe lots of people hiding their score in local. Another game burns tokens and time.........

> **Kim Jinseo** · 2026-10-26 23:37 UTC · 0 votes
>
> Do you have any clue of how the private scoring would work? Wouldn't that possibly cause a shakeup?

> **(deleted user)** · 2026-10-27 10:14 UTC · 0 votes
>
> *(deleted)*

**mrmeows** — **Topic author** · 2026-10-28 00:21 UTC · 0 votes

Hello world!

Have to do some IRL stuff right now but I am still here!

When I come back I will release up to a 25x solution for general public, this will most likely be last of my public dumps unless I realize that there is much more headroom.

I am more than willing to team up with anyone who is serious and has compute (i have compute as well, not a one ended relationship)

thanks

> **D Chakraborty** · 2026-11-11 17:03 UTC · 0 votes
>
> It's been 15 days. Our lord and savior mrmeows promised 25x and vanished into the real world 😭
>
> Please return, my lord. The peasants are still stuck at the plateau, the agents are burning tokens and we're running out of public circuits to copy.
