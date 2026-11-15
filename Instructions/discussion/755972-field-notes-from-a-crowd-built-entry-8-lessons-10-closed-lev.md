# Field notes from a crowd-built entry: 8 lessons + 10 closed levers (visual summary)

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755972 — copied 2026-11-15 (raw Markdown, unedited).

**Lynn Langit** · 2026-11-05 01:59 UTC · 0 votes

We reached the top 5% largely by verifying and blending public work, so here’s something back: a one-page visual summary (partial result shown below) of what two weeks of re-verifying, profiling and testing taught us.

[Circuit Cost Field Notes](https://claude.ai/artifact/EJstp5r8E4AcsZTFHZHNS7)

![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F311896%2F59de634923987483853306b7f705add5%2FScreenshot%202026-11-04%20at%208.01.47PM.png?generation=1793844131933833&alt=media)

What’s in the full summary:

- What each task actually charges (Q × C, R × Q² × T, λ_eff × T × Q) and where that points the digging.
- 8 principles with evidence, e.g. one failed DQI shot in 102,400 costs as much as ~1,840 Toffolis, and the 892-qubit decoder’s peak is ~200 brief moments, not one bottleneck.
- A negative-results ledger: 10 levers we measured and closed (sub-42 GF(2¹²) multipliers, pebbling the DQI peak, dropping rare branches, deleting never-firing invsqrt gates, chemistry weight refits, loader product search, 7-bit words, and more), with the tool used and the result.
- Three questions to ask before reaching for a new tool: does it cut what’s charged, does the circuit have that kind of slack, and is it safe under an unknown seed?

The most expensive lesson: our invsqrt gate deletions passed every input on the default seed, then failed the grader with a 0.5-turn phase error. Those gates were phase fixes for rare measurement outcomes. Check measurement-dependent edits on several salts.

Thank you to everyone whose Apache-2.0 work this rests on, especially Ehan Mallick, Kirill Tymchishin, SomeJoint, busyaprime, yoyobar, Syed Asad Ali, Ahmed Berat Özer, Lera, Ayodeji, heliosli, Dmitrii Gluzdov, Atakan Aldemir, Isabe11aM00re, Alex Chilton, and manwithacat for the daily briefings, blend notebook and fast verifier.

Made with Claude. Corrections welcome.
