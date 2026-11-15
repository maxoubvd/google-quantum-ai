#   Permission request: reuse of the QuantumAI DQI decoder (jday96314)

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753660 — copied 2026-11-15 (raw Markdown, unedited).

**Mukkai Krishnamoorthy** · 2026-10-28 22:08 UTC · 5 votes

Hi @jday96314, 
Many public DQI circuits in 
this competition (the Schur / continued-fraction lineage by Alex Chilton, Ahmed Berat
Özer, SomeJoint, heliosli and others) credit your QuantumAI decoder as their origin. heliosli's release notes that your repository stated no license and asks people to contact you before reusing the code. The repository now seems to be unavailable.

Our team would like to submit a compiled DQI circuit from that lineage (heliosli's public 889-qubit build), with full credit to you. Because of the winner source-delivery rules, we'd rather ask than assume. Would you be willing to:
1. grant permission to use circuits derived from your decoder in competition submissions, with attribution, or
2. state a license (e.g. Apache-2.0, like most public work here)?

If you'd prefer that derived work not be used, we'll respect that and build our own decoder instead. Thank you for publishing it in the first place.


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**James Day** · 2026-10-29 22:57 UTC · 1 votes

1. I don't mind you using it. I've already accepted that people will use it regardless of what I say.
2. If the prize contenders depend on it and the hosts are concerned or I rank in the gold medal range then my full repository will be made public again with an apache 2.0 license included in the repository. The repository will contain an archive of the version people are using + any optimizations I've made since then.

> **D Chakraborty** · 2026-10-30 01:24 UTC · 1 votes
>
> 🤖 DO YOU CONSENT TO YOUR MISTAKE BEING DISTRIBUTED UNDER APACHE 2.0??? 🤖
>
> Man, I actually feel bad for you at this point LMAOO

> **hengck23** · 2026-10-30 01:52 UTC · 3 votes
>
> in summary:
> - jday96314 has the right and option to disclose or not his source code
> - the host  has the right and option to decide whether any solutions without jday96314's source code are valid solutions eglieble for prize money or not
> - the host can ask jday96314 to disclose his source code after the prize award, but it is up to jday96314 to decide
> - the host can yet set an additional prize if he seems fit  (e.g. milestone opensource)
> - kaggler community can suggest to the host for the additional prize mentioned above
>
> ---
>
> There is also a technical consideration: KMB files already contain the circuit wiring, so it may be possible to reconstruct a functionally equivalent implementation without access to the original source code. This makes the distinction between sharing source code and sharing a compiled circuit somewhat less clear.
>
> ---
>
> ### Most importantly, if you appreciate jday96314's work, go beyond simply giving credit ... help him earn his gold! For example, share your experimental results and improvements based on his work.
