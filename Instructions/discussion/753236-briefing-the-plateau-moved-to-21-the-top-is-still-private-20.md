# Briefing: the plateau moved to 21, the top is still private — 2026-10-27

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753236 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-27 17:07 UTC · 1 votes

Prepared with AI assistance. Numbers are labelled: *observed* (leaderboard), *author-reported* (the author's
official-verifier run), or *ours*. Window: 26 Oct 17:45 UTC to 27 Oct 17:00 UTC. Corrections welcome.

## 1. Leaderboard (27 Oct, 16:57 UTC, observed)

| | 26 Oct | 27 Oct |
| --- | --- | --- |
| Top score | 23.35 | **24.53** (Breaking1800) |
| Teams ≥ 23 | — | 7 |
| Teams ≥ 20 | 10 | **481** |
| Teams ≥ 15 | 300 | 623 |
| Teams | 1,124 | 1,275 |

The top ten are now 22.86–24.53 (Breaking1800, Schrödinger's cat lived, Thinking..., Zyy7390, Jack, Distortion!!,
Schrödinger's ZZZ, Slingshot Effect, shyjin, Riemann).

Why ≥ 20 went from 10 teams to 481 in a day: several complete public portfolios now score 21.2–21.5, and many teams
submitted them. The top moved by only about 1.2 points, and everything public sits below 21.6, so the leaders' extra
three points come from work that isn't public.

## 2. Public frontier by task (author-reported, official verifiers)

| Task | Best public | Qubits | Cost | Ratio to reference | Source |
| --- | --- | --- | --- | --- | --- |
| invsqrt | mrmeows' circuit (used unchanged by everyone) | 63 | 103,546 (mean 1,348.5 Toffolis) | 20.06 | [21.0 circuits, let's all reach the plateau](https://www.kaggle.com/code/mrmeows/21-0-circuits-let-s-all-reach-the-plateau) — Apache-2.0, generators and certificates included |
| DQI | Ahmed Berat Özer v9, two register reclamations on SomeJoint's 889-qubit fused decoder | 877 | 4.421 × 10¹¹ (mean 482,543 Toffolis) | 55.28 | [Ahmed: DF-THC Group Compression](https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression); [SomeJoint: 21.34 Local GM: Fused DQI](https://www.kaggle.com/code/somejoint/quantum-circuits-21-34-local-gm-fused-dqi) |
| Chemistry | Ahmed v9: Michael Rustad's R10/B129/C46 8-bit factors, heliosli's five-word angle circuits, finer quantisation and a λ loop | 431 | 3.552 × 10⁸ (17,330 Toffolis, λ_eff 47.55, error 19.998 mHa) | 9.45 | [Ahmed](https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression); earlier [heliosli J5](https://www.kaggle.com/code/heliosli/qcc-j5-chemistry-loaders-deeper-lookup-plans) (8.91), [Rustad](https://www.kaggle.com/code/michaelrustad/vqcc-r10-c46-four-word-chemistry-21-41-lb) (9.09) |

Author-reported local geometric mean for Ahmed's portfolio: 21.88. Public notebooks claiming leaderboard scores
range from 21.17 (Alex Chilton v44) to 21.53 (leoprovorov / saddamaditya).

## 3. What's new

- **Narrower DQI (Ahmed Berat Özer).** Two register reclamations take SomeJoint's 889-qubit decoder to 877 qubits
  for +1.7% Toffolis, a net −1.0% on Q²T, with zero value and phase errors in the author's run.
- **Fused DQI (SomeJoint / Agniv Banerjee).** A fused confluent-Hermite decoder with 17 jet nodes at 889 qubits,
  474,308 mean Toffolis.
- **Chemistry at 431 qubits.** The gains combine 8-bit factors from quantisation-aware training (Rustad) with
  five-word angle circuits (heliosli) and Ahmed's finer quantisation schedule.
- **marculera** merges DF-THC groups before dropping rows: 22% fewer rows at that stage.
- **avazaku** gives a lower bound for coefficient-only 7-bit refits: 156–199 mHa against the 20 mHa gate. 7-bit
  rounding can't work without changing the dictionary itself.
- **Alex Chilton's "What didn't work"** ([753185](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753185))
  lists about 30 measured nulls across all three tasks. Worth reading before starting any of them.
- **Valerio Dolci** ([752949](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752949))
  describes one person plus Claude plus a local open-weight model on two consumer GPUs running a build–verify–keep loop,
  and reports beating the best public chemistry circuit by about 5% (no numbers posted yet).

## 4. From the forum

- **"How relevant is bit-aware training for chemistry? Can we use 9-bit BROT?"** (hengck23, in 752895). The public
  evidence says bit-aware training is what makes low bit widths work at all: re-rounding public 10-bit words to 8 bits
  without retraining gives 162–201 mHa against the 20 mHa gate, and 7- and 6-bit rotations give 350 and 795 mHa
  (Alex's null list); the 8-bit frontier factors come from quantisation-aware training. We haven't measured 9-bit
  ourselves; it should be easier to fit than 8, at the price of the qubits and Toffolis the extra bit costs.
- **Are 63 qubits (invsqrt) and 877 (DQI) floors?** For invsqrt nobody has found a 62-qubit circuit (Alex's list; our
  notebook [Is the 63-qubit inverse square root optimal?](https://www.kaggle.com/code/manwithacat/is-the-63-qubit-inverse-square-root-optimal)
  looks at why). For DQI, Alex measures one qubit as worth about 1,100 Toffolis at this size, which is the exchange
  rate any further narrowing has to beat.
- **DQI phase errors near 8 × 10⁻⁵** ([752344](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752344)):
  two phase errors cost about +0.85% on the cost ratio (xhsmsuper); shawncsx reports a clean 0 / 102,400 with a
  margin-1 parametrisation.
- **Unanswered and worth an answer:** whether ANGLE_k must satisfy Eq. (72) for all inputs or only reachable ones
  (Chat8903); whether private scoring will reshuffle the board (752757).

## 5. Tools and notes from us

- **Verifying big circuits locally** ([753230](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753230)):
  a kickmix build that takes local DQI verification from 7.8 to 2.0 GiB peak at one worker (9.0 to 5.6 GiB at six), with
  identical results, plus a way to use all cores on Windows. Since that post, two more changes are merged (runs of CX
  gates without per-op dispatch, cheaper counting of conditioned ops): full local verification of all three tasks on
  a Mac is now about 67 s (about 115 s before the profiling-driven changes; DQI alone about 1.46× faster than the post's build), still
  field-for-field identical. Install: `pip install "git+https://github.com/manwithacat/kickmix-macos.git#subdirectory=kickmix"`.
- [Is the 63-qubit inverse square root optimal?](https://www.kaggle.com/code/manwithacat/is-the-63-qubit-inverse-square-root-optimal)

## 6. What we think it means (opinion)

The public plateau is now about 21.5, shared by hundreds of teams, and the remaining public gains are small
(a percent here and there). The gap to 24.5 is about 14% in the geometric mean, so the leaders are probably ahead on
more than one task. The tools that made the last two days possible are verifier-in-the-loop searches, and the
bottleneck for most people is how many candidate circuits they can verify. That's why we're putting effort into
making verification cheap on ordinary hardware.


---

## Comments (2 threads)

Replies are indented with `>`, one level per reply depth.

**D Chakraborty** · 2026-10-27 18:29 UTC · 0 votes

Meanwhile, at Alliance headquarters...

![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F5819143%2F25db4280ee974c2da40dcfc86d4a3397%2F550e23b7-4a82-4e96-82fb-19106859474b.png?generation=1793125665689231&alt=media)

Currently at 0 MiB used, 0% utilization , the Alliance hasn't even drawn its weapon yet.


![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F5819143%2F7f489f2294c233e0ff67a325f0b7e4f7%2FScreenshot%202026-10-11%20at%201.10.57PM.png?generation=1793125748880371&alt=media)

Stay tuned.

**Valerio Dolci** · 2026-10-27 17:29 UTC · 0 votes

Well the big jump was on the second problem a notebook that wasn't supposed to be released open, that was made public by mistake. Is a standard mistake from agents now days. Right now i think that there is some optimisation left in the 3' problem but it need lots of gpu hours and it needs some smart training, and maybe some luck. Still i feel that the real surprise will come from the first problem where we all are close to theoretical optimum, so some agent will find a complete new way forward surprising everyone, this ot all this challenge will be saturated in the first week and we will have the same scores for months
