# What are you guys scoring in each task?

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752550 — copied 2026-11-15 (raw Markdown, unedited).

**mikelou1** · 2026-10-25 10:40 UTC · 7 votes

My score in each task

| Task | Score |
| --- | --- |
|  `invsqrt` (Fusion Energy) | 20.2× |
|  `dqi` (Combinatorial Optimization) | 30.5× |
| `dfthcsos` (Quantum Chemistry & Clean Energy) | 8.25× |

Please share how much you are scoring in each task!

I have updated it to reflect my performance recently


---

## Comments (8 threads)

Replies are indented with `>`, one level per reply depth.

**Mehrab Hosain Robin** · 2026-10-26 13:49 UTC · 2 votes

A practical note for Windows users running the dqi/invsqrt verifiers locally with many workers: under memory pressure the parallel run can die silently (empty stdout/stderr) or raise MemoryError, which looks like a FAIL. Check that free commit memory is at least about 8 GB before a full run, and treat an empty-output failure as an infrastructure retry, not a circuit failure. It saved us several false alarms.

**Scott Willis** · 2026-10-25 15:49 UTC · 2 votes

| Task |Score  |
| --- | --- |
| invsqrt (Fusion Energy) | 22.68× |
| dqi (Combinatorial Optimization) | 40.06× |
| dfthcsos (Quantum Chemistry & Clean Energy) | 7.10× |

**Muhammed Ammaar Yılmazoğlu** · 2026-10-26 10:31 UTC · 1 votes

| Task | Score |
|---|---|
| invsqrt (Fusion Energy) | 19.73× |
| dqi (Combinatorial Optimization) | 36.87× |
| dfthcsos (Quantum Chemistry & Clean Energy) | 7.99× |

**Subh** · 2026-10-26 07:13 UTC · 1 votes

| Task | Verified Score |
|---|---:|
| invsqrt (Fusion Energy) | **15.9957×** |
| DQI (Combinatorial Optimization) | **36.2249×** |
| Chemistry (Quantum Chemistry and Clean Energy) | **8.1357×** |

> **DePie** · 2026-10-26 09:59 UTC · 0 votes
>
> Well. If you get public solution for invsqrt. Can you get top 5 ?

> > **Subh** · 2026-10-26 13:31 UTC · 1 votes
> >
> > I tried to improve my previous inv and chem task scores but ended up improving DQI further 😅
> >
> > Updated Scores:
> >
> > | Task | Verified Score |
> > |---|---:|
> > | Invsqrt (Fusion Energy) | **18.19723×** |
> > | DQI (Combinatorial Optimization) | **50.78074×** |
> > | Chemistry (Quantum Chemistry and Clean Energy) | **8.40553×** |

> > **Mehrab Hosain Robin** · 2026-10-26 13:55 UTC · 0 votes
> >
> > 50.78× on DQI is a huge jump, congrats! That’s about 1.6× beyond anything else posted here. Out of curiosity, did most of the gain come from qubits or from Toffolis? No details needed if you’re keeping it close.

> > **Subh** · 2026-10-26 16:47 UTC · 0 votes
> >
> > Thanks! 😄  A bit of both, but reducing qubits made the bigger difference initially. Pushing it more I was actually able to trade extra Toffolis for fewer qubits - that gave a small bump over 50

> > **Yair Be'ery** · 2026-10-27 05:09 UTC · 0 votes
> >
> > We achieved same DQI score :)
> >
> > Moving forward..
> >
> > Note: Currently working only on the DQI task

> > **Subh** · 2026-10-27 07:17 UTC · 0 votes
> >
> > That's amazing! I was able to reach 55.14x in DQI

> > **Mehrab Hosain Robin** · 2026-10-27 07:19 UTC · 0 votes
> >
> > Very good job Man

> > **Yair Be'ery** · 2026-10-27 19:16 UTC · 0 votes
> >
> > We achieved 55.78× in DQI :)

**Meet Mendapara** · 2026-10-26 17:25 UTC · 0 votes

| Task | Score |
| --- | --- |
| **invsqrt** (Fusion Energy)   |  19.13x |
| **dqi** (Combinatorial Optimization)  |  49.0x |
| **dfthcsos** (Quantum Chemistry & Clean Energy)  |  6.78x |

**DePie** · 2026-10-26 10:03 UTC · 0 votes

| Task | Score |
| --- | --- |
| invsqrt (Fusion Energy) | 16.10× |
| dqi (Combinatorial Optimization) | 27.19× |
| dfthcsos (Quantum Chemistry & Clean Energy) | 7.924× |

Still updating

**venks** · 2026-10-25 16:33 UTC · 0 votes

interesting, i'm at 14.77× invsqrt, 25.36× DQI, 4.09× chemistry

**shanzhong8** · 2026-10-25 11:43 UTC · 0 votes

got 5+ on task dfthcsos
