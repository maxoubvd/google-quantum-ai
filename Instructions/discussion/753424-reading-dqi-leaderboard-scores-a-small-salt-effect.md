# Reading DQI leaderboard scores: a small salt effect

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753424 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-28 03:43 UTC · 4 votes

The DQI decoders in the current public lineage (ours included) fail on a few inputs in every 100,000, by design. The rules allow this: the repetition factor R prices it, and it never fails a submission. But each DQI score is a random draw. The test inputs are seeded from the circuit's hash plus a salt, any change to the circuit re-rolls the public draw, and final scoring uses a private salt.

We measured the public DQI circuits under 9 salts each ([notebook](https://www.kaggle.com/code/manwithacat/dqi-salt-check-public-frontier)). The 877- and 889-qubit decoders fail on about 1 in 100,000 shots and draw zero failures only about a third of the time. The 883- and 897-qubit decoders fail on about 1 in 300,000, and a 1,893-qubit decoder never failed. The expected cost of a "lucky" zero-failure public draw is 0.05–0.16% of the overall score. Our own 22.00949 was such a draw; its expected score is about 21.998.

That's small, but if you're choosing final submissions between near-ties on public score, compare expected DQI cost rather than the observed one. The notebook runs weekly with the unchanged official verifier. Corrections welcome. (Done with AI assistance.)


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**manwithacat** — **Topic author** · 2026-10-28 05:26 UTC · 0 votes

A follow-up on *what* these rare DQI failures are, since we now have enough of them to look at.

We logged 47 failing shots across four public decoders (the 877-, 883-, 889- and 897-qubit circuits; 22 of them from 10.24M shots of the 897Q). All 47 behave the same way:

- They're deterministic in the input: replaying a failing input with fresh simulator seeds fails every time, and passing inputs never do. So they're properties of particular error patterns, not measurement luck.
- They're all phase errors. The syndrome and ancillas come back clean.
- Probing the circuit with single-bit masks shows which error word it "believed" it decoded. In every case that word has the right syndrome but shares none of the true error positions; it fills positions 0..69 exactly, occasionally plus a few extras.

Positions 0..69 are exactly the span of the classical-mask projection these decoders share (the 70-position Vandermonde term that subtracts the masks' known linear phase). Its contribution equals the phase of the unique word on positions 0..69 with that syndrome. So the simplest reading is that, on these rare inputs, the decoder finds no error positions (presumably an internal degenerate step, like a zero pivot), and what's left in the phase is just the projection's term.

That's a black-box reading; we haven't looked inside any decoder. If you're building or modifying one with the projection, flagging the "found nothing" case explicitly would at least make these failures visible and countable. The trigger rates differ between designs (roughly 1 in 100,000 to 1 in 500,000 shots), so it's worth measuring yours under several salts. Logs and probe scripts are on request. (Done with AI assistance.)
