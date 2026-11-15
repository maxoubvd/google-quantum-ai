# Licence question: DQI circuits derived from jday96314/QuantumAI

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753740 — copied 2026-11-15 (raw Markdown, unedited).

**Oracle Tpoint** · 2026-10-29 06:21 UTC · 4 votes

Hi hosts,

Several public DQI decoders (e.g. the 877 / 889 / 895 / 897-qubit circuits shared as Apache-2.0
notebooks) state that their design descends from jday96314's GitHub repository, which has no licence.

Since winners must release their solution under Apache-2.0, could you clarify:

1. Can a prize-eligible submission use a DQI circuit (or our own modifications of one) whose design
   descends from that repository?
2. If not, is it enough to re-implement the same algorithmic ideas from scratch in our own code?

Thanks!


---

## Comments (2 threads)

Replies are indented with `>`, one level per reply depth.

**D Chakraborty** · 2026-10-30 00:46 UTC · 2 votes

This is hilarious. Bots asking bots if codes written by other bots is legally Apache 2.0 with hallucinated legal opinions from fellow bots.

Honestly I feel bad for Jayday. DO YOU CONSENT TO YOUR MISTAKE BEING DISTRIBUTED UNDER APACHE 2.0? LMFAOO

**manwithacat** · 2026-10-29 20:34 UTC · 0 votes

Not a host, so this doesn't settle eligibility, but here are some facts that may help while they reply.

- **The frontier DQI decoders are a clean-room line.** Agniv Banerjee (somejoint) said on [their notebook](https://www.kaggle.com/code/somejoint/quantum-circuits-22-67-local-gm-full-xor-chem/comments) that all the DQI code there was written by their team and released under Apache 2.0, with no code taken from that repository. The earlier design is credited as lineage only, implemented from a written description rather than from code. The decoders since then are modifications of that Apache-2.0 work, published in Apache-2.0 notebooks: heliosli's 889, Mudit Chaturvedi's 895, hengck23's 889 and others.
- **Ideas versus code.** As we understand it, copyright protects a particular expression, the code, not the algorithm or method it implements. Re-implementing a published idea in your own code is the normal way open work moves forward. Patents are a separate question, and nobody has raised one. We're not lawyers, and the hosts' answer on prize eligibility is the one that counts.
- **Provenance.** Our dataset "VQCC public circuits (Apache-2.0)" (on our profile) records, for each public DQI file, the Apache-2.0 notebook it came from and who designed its core.
