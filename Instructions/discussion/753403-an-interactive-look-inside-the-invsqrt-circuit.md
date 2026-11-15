# An interactive look inside the invsqrt circuit

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753403 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-28 02:43 UTC · 6 votes

We built a small viewer that runs mrmeows' 63-qubit invsqrt circuit (the 1,348-Toffoli build most of the frontier uses) in your browser, one input at a time: https://manwithacat.github.io/kickmix-macos/viewer/

Pick an x and watch its bits move through all 46,610 operations on 63 qubits: 17 stages from normalisation through the coefficient lookup, the arithmetic blocks, the 41 fine rotations and the long phase-repair tail. Lit bits are coloured by the phase so far, and a dial shows the phase closing in on 1/√x.

- **Play** runs the circuit start to finish (once, on a loop, or back and forth); **Sweep x** scrubs the input across its whole range so you can see where the circuit's behaviour changes. **Re-roll measurements** changes every mid-circuit measurement outcome; the middle shifts, the final phase doesn't.
- Tracks show Toffolis run, measurement-based clean-ups, branches skipped for this x, and which gates move the phase.
- **Precision across every x**: error ÷ tolerance for all 3,133,452 inputs (computed with kickmix). The worst is x = 161 at −0.99932; only 184 inputs are above 0.9 of tolerance and the mean is 0.21.
- **Qubits**: at most 35 of the 41 scratch qubits hold a 1 at any instant, yet all 41 are needed, because their lifetimes are staggered (shown as a background in the raster). Saving a qubit means re-packing lifetimes, not lowering a peak.
- **Load your own .kmx** to view your circuit the same way. It stays in your browser. It's a viewer, not a verifier.

The in-browser tracer matches kickmix exactly (all qubits, classical bits and the phase, on 2,299 inputs). Stage names are marked documented (from mrmeows' README and proofs) or inferred (read from the op stream), and the middle arithmetic stages are the least understood. If you spot something odd, or know what a stage is doing, please say so; that's half the point.

For fun, there's also a price tag: at rough 2026 numbers, running this "small" circuit with error correction needs about 48,000 physical qubits, around $2 billion of hardware, and about three weeks to do the verifier's 200 million shots. kickmix does the same check on a laptop in about 75 seconds.

The viewer is also in our kickmix fork (tools/circuit-viewer, https://github.com/manwithacat/kickmix-macos) if you'd rather run it locally.

Circuit by mrmeows (Apache-2.0), on the line of busyaprime, Phuc Thai, Dmitrii Gluzdov and marsyas. Built with AI assistance.


---

## Comments (2 threads)

Replies are indented with `>`, one level per reply depth.

**hengck23** · 2026-10-29 00:55 UTC · 1 votes

this is not a visualisation tool. it is a data capturing too! collect the register traces as data, ask AI agent to look for correlation and suggest how to compress that ciruit (e.g. ahare wire for qubits of same trace pattern)

> **manwithacat** — **Topic author** · 2026-10-29 01:16 UTC · 0 votes
>
> that is a really great idea and I am definitely not face palming that I didn't think of it

> **manwithacat** — **Topic author** · 2026-10-29 23:46 UTC · 0 votes
>
> Thanks again for the "data capture, not visualisation" idea. Small update: the viewer now opens on the current public frontier circuit, busyaprime's 9a8627 (mrmeows's 63-qubit design with busyaprime's re-synthesised coefficient lookup, 186 → 172 ANDs, and window S1, 30 → 25). mrmeows's 0766a111 is still selectable from the new Circuit menu.
>
> https://manwithacat.github.io/kickmix-macos/viewer/
>
> For anyone who wants to try the trace-correlation approach: the per-op data the page uses (`ops.json`, `stages.json` and the qubit lifetimes in `livepeak.json`) sits next to the page in `tools/circuit-viewer/` of https://github.com/manwithacat/kickmix-macos, Apache-2.0. The in-browser tracer matches kickmix bit for bit (all 63 qubits, all classical bits and the phase) on 2,299 inputs for both circuits.
>
> If you have any suggestions for improving it, I'm happy to take ideas.

**Lynn Langit** · 2026-10-29 18:08 UTC · 0 votes

Useful - thanks for sharing 🙂
