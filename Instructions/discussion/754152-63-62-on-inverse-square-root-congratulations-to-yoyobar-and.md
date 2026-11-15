# 63 → 62 on inverse square root: congratulations to yoyobar and Ahmed Berat Özer (and our notes on why it was hard)

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754152 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-30 15:28 UTC · 1 votes

*Corrected 30 Oct: the first public 62-qubit circuit is **yoyobar**'s; my first version of this post wrongly credited heliosli, who republished it with credit to yoyobar. Apologies to yoyobar.*

Congratulations to **yoyobar**, whose `a568bc793481` was the first public 62-qubit inverse square root ([LB 22.99 | invsqrt 101 - 62 qubits](https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits)). Congratulations also to **Ahmed Berat Özer**, whose `6dbb98b34c3b` (11:12 UTC today) builds on it and saves one more Toffoli at 62 qubits. Thanks to **heliosli** for the careful re-verified republication. Both circuits are in the public circuits dataset with credits, and the Friends notebook now picks Özer's.

We had spent a while convinced that 63 was a wall. Here is what we found, and how yoyobar got past it.

**Why 63 looked tight.** In the 63-qubit design (mrmeows's, re-synthesised by busyaprime into `9a8627c81283`), one wire, q23, holds the flag [x < 2²⁰] alone for most of the circuit. It is not a function of the other 62 live wires: there are 612,723 conflicting states over the 3,133,452 legal inputs. So it can't be parked, recomputed or measured away. For that op stream, the width census said 63 → 62 was blocked. That analysis is in our note *Is the 63-qubit inverse square root optimal?* ([doi:10.5281/zenodo.23781241](https://doi.org/10.5281/zenodo.23781241)).

**How yoyobar did it.** a568 doesn't touch q23. It keeps the same angles and the same first 239 operations, and re-synthesises **every** 63-wide stage one wire narrower:
- The coefficient lookup builds 27 of its 32 AND products directly in the output words, so only 5 temporaries remain. The rest are erased with a block of measurements and Clifford fix-ups.
- The later stages replace fresh helper qubits with classical parity logic and linear forms XORed onto wires that are already live.

The cost is about +4.5 Toffoli-equivalents for one qubit. A qubit is worth about 26 at this size, so the trade is clearly worth it. Özer's edit replaces one controlled swap whose AND term equals an affine form of live wires: a free wire, two CX, a measurement and a CZ fix-up, for −1 Toffoli.

**What we checked afterwards (exact, over all inputs).** We swept all 1,501 Toffoli-class gates in 6dbb for further affine or single-AND replacements of this kind. There are none: every apparent hit sat inside a measurement-conditioned fix-up block and was refuted on other branches. So 6dbb looks locally tight to that kind of edit. Our note's "optimal for its architecture" claim was too strong on cost. A v0.3 is coming that narrows it, credits yoyobar and Özer, and records the re-synthesis route that worked.

**Quantum Voodoo.** Elsewhere I joked that I was delaying publishing because Claude was trying to explain it to me and it looked like voodoo. That notebook is now public: [VQCC Quantum Voodoo](https://www.kaggle.com/code/manwithacat/vqcc-quantum-voodoo). It is my attempt to understand, with animations, what is going on:
- circuit surgery on the DQI decoder (889, then 888 qubits by measurement-based uncomputation);
- why DQI scores depend on a random draw;
- a fix for a chance-zero failure mechanism, so the decoder fails on fewer draws instead of relying on luck;
- a constrained refit of the chemistry weights;
- a section on the 63 → 62 story above.

Its version history is the record of when each piece existed. Everything is Apache 2.0, and the circuits are in the [VQCC Quantum Voodoo dataset](https://www.kaggle.com/datasets/manwithacat/vqcc-private-dqi-candidates). It scored 23.13505 on the public board.

**What's next.** We're doing further work on our visualisation and circuit-surgery tool set, so that anyone (or their agent) can make animations like these from a .kmb file in a couple of commands. It will go into the open kickmix-macos fork, and I'll post here when it's ready.

*AI assistance: the analysis and tooling were built with Claude (Anthropic) under my direction; the conclusions and any errors are mine.*


---

## Comments (2 threads)

Replies are indented with `>`, one level per reply depth.

**manwithacat** — **Topic author** · 2026-10-30 18:04 UTC · 0 votes

Update: the renderer mentioned above is out. [`tools/circuit-viewer/render.py`](https://github.com/manwithacat/kickmix-macos/blob/main/tools/circuit-viewer/render.py) in the open kickmix-macos fork turns any `.kmb` into explanatory PNGs and animated GIFs, one command per picture:

- `width`: live qubits stage by stage, two circuits side by side (e.g. 63 vs 62, 889 vs 888)
- `diff`: what changed between two circuits
- `trace`: one shot's register values through the circuit
- `cost`: Q, Toffolis and cost per stage from the official formulas
- `ledger`: a double-entry bookkeeping view (values opened and settled, measurement-based uncomputation as accruals, width as capital, recomputed values as paying twice)

It needs only numpy, matplotlib and Pillow (kickmix only for DQI traces) and prints a short JSON summary, so it works well driven by an agent. Quickstart, examples (including yoyobar's 62-qubit circuit against the 63-qubit design) and a Kaggle cell are in [RENDER.md](https://github.com/manwithacat/kickmix-macos/blob/main/tools/circuit-viewer/RENDER.md). Apache-2.0; issues and suggestions welcome.

**yoyobar** · 2026-10-30 15:39 UTC · 0 votes

Hi~ I’m a little sorry, but did you and your AI read the reference in their notebook? Where they did cited my work… QAQ


You hurt my heart… QQQQQQQ QQQQQ
(https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits)

> **yoyobar** · 2026-10-30 15:46 UTC · 0 votes
>
> Is your Claude really from Anthropic? Doesn’t it actually come from a proxy or agency that routes your requests to Doubao?
>
> ![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F8293010%2Fbc8a59667a2ecadca66b1be4ad6e4217%2Fimages.jfif?generation=1793375153243197&alt=media)

> > **manwithacat** — **Topic author** · 2026-10-30 15:51 UTC · 0 votes
> >
> > Yes. It's Claude Code with Anthropic's own models, no proxy.

> **manwithacat** — **Topic author** · 2026-10-30 15:51 UTC · 0 votes
>
> You're absolutely right, and I'm sorry. The first public 62-qubit circuit, `a568bc793481`, is yours: [LB 22.99 | invsqrt 101 - 62 qubits](https://www.kaggle.com/code/yoyobar/lb-22-99-invsqrt-101-62-qubits). heliosli's notebook credits you for it, and Ahmed Berat Özer's −1 Toffoli version builds on it. Our scan picked up heliosli's copy and credited it as the origin. That was our error: a stale timestamp in our pipeline made the replication study look earlier than your release. I'm correcting the credit in the public circuits dataset, and I've edited the post above. Congratulations on the 62, and thank you for the correction.
