# The mandatory "let's rewrite it in Rust" post: kickmix-engine, DQI verification in a fraction of the time and memory

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753642 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-28 21:08 UTC · 3 votes

Yes, this is that post. Some of you will remember our kickmix fork ([753230](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753230)),
which made local verification leaner by patching kickmix itself. This week we went one step further and
re-implemented kickmix's simulator in Rust: **[kickmix-engine](https://github.com/manwithacat/kickmix-engine)**
(Apache-2.0, v0.4.0). It runs the competition's official verifier logic for invsqrt and DQI, with the engine as the
simulator, and returns the same summary JSON. (i.e. I typed. "rewrite this in rust don't make any mistakes" into Claude).

**Same answers first.** It's an independent re-implementation, so we checked it against kickmix shot by shot before
measuring any speed: every qubit, every classical bit, the exact phase and the op counts, on both frontier DQI
decoders and on invsqrt. That includes the shots that fail, so failures match too. Full verifications match the
official verifier field for field (including `reaction_depth`), and so do the FAIL cases we built: dirty ancillas, x
not restored, dirty scratch. Add `--cross-check` and it runs the official verifier alongside and diffs every field.
The official verifier stays the reference.

**On a Kaggle CPU notebook** (4 vCPUs, 4 workers each), from our public quickstart notebook, everything installed
offline from wheels. Every build returned identical results on every circuit:

| circuit | official kickmix | our kickmix fork | kickmix-engine |
| --- | --- | --- | --- |
| sample DQI | 220 s, 9.0 GB | 172 s, 3.4 GB | **84 s, 0.6 GB** |
| frontier DQI (heliosli's 889-qubit `4ee7a17b4cbc`) | 283 s, 4.6 GB (separate run) | 83 s, 1.9 GB | **29 s, 0.36 GB** |
| sample invsqrt | | 803 s | **317 s** |

On the frontier DQI that's about 10x official kickmix and 2.8x the fork, at a fraction of the memory. (Memory here is
the process tree's summed PSS; the official frontier-DQI row comes from a separate private benchmark on the same
machine type, so treat that ratio as approximate.)

**Why not just release the GIL?** That was our first guess too, and it was mostly wrong. What we measured:
- Threads help less than you'd hope. With the verifier's batches on threads and the GIL released around the
  simulator, DQI got 2.7x and invsqrt only about 1.2x. That's why the official verifier uses separate processes, which sidestep the GIL
  completely. So for speed, the GIL isn't the bottleneck you hit.
- Its real cost is memory. Processes mean every worker holds its own copy of the circuit and simulator state. That's
  the 9 GB on the sample DQI: the price of working around the GIL, not of the simulation itself.
- Profiling invsqrt (on a Mac, one worker) showed where the time actually goes: 46% simulating; 38% inside the
  simulator's `read_shot_phase`, which builds a Python `Fraction` object for every shot's phase; 6% converting those
  Fractions back to integers; about 10% NumPy checks and sampling.
- So compiling the Python side can't do much. We tried Numba (nogil) on everything it could reach: 1.1x, because
  the real per-shot Python is only 16-18% of the time. The expensive part is the *boundary*: one Python object per
  shot crossing from C++ into Python, and Python orchestrating every batch.

**Where the speed actually comes from**, roughly in order of size:
1. **Removing the per-shot boundary.** The engine never makes a Python object per shot. The sampler, the simulation,
   the phase readout (as plain 64-bit integers) and the checks all happen in one native call per batch. We first
   proved this inside our kickmix fork: a batched integer phase read-out plus compiled checks made invsqrt 1.7x faster
   than the fork without them (Mac). The engine takes it all the way, including a draw-for-draw Rust restatement of the DQI sampler
   (proven identical on 2,400 batches).
2. **One process, real threads, one circuit.** With no Python in the loop, the GIL is simply released for the whole
   batch. Threads then scale (7.2x on 10 Apple cores) and share a single compiled copy of the circuit (257 MiB for the
   frontier DQI). That's where 9 GB becomes 0.6 GB, and why the engine doesn't need processes at all.
3. **CPU kernels.** On x86 the wheel carries baseline, AVX2 and AVX-512 builds of the hot loop and picks the best one
   at start-up. On Kaggle's Xeon, AVX2 was 1.42x faster than baseline at 4 threads. A compiled circuit is also cached,
   so repeat runs skip the ~1 s compile.
4. **A faster inner loop, but less than you'd think.** The bit-sliced core (256 shots per word, as in kickmix) with a
   fixed-width compiled op stream and no per-op dispatch on CX runs is 0.80x kickmix's time per batch on one core.
   Useful, but it's the smallest of the four.

So what Rust bought us is mostly a place to put the whole per-shot loop where Python can't see it. Any language that
let us do that would have worked; Rust made it safe to do with threads, and easy to ship as one wheel per platform.

invsqrt gains less (about 2.5x the fork) because its per-shot checks still run in NumPy on the main thread: point 1
is only half done there. Moving that into Rust is next on the roadmap.

**Try it:**
- Copy & Edit **[kickmix quickstart](https://www.kaggle.com/code/manwithacat/kickmix-quickstart)**. It runs offline,
  shows all three builds side by side, and has a cell to verify your own circuit.
- Or in any notebook: attach [kickmix-engine-wheels](https://www.kaggle.com/datasets/manwithacat/kickmix-engine-wheels)
  plus the competition data, then
  `pip install --no-index --no-deps /kaggle/input/kickmix-engine-wheels/kickmix_engine-*-manylinux_2_28_x86_64.whl`
  and `kickmix-engine verify dqi your.kmb --verifier <competition>/verifiers/dqi/verify_submission_dqi.py --cross-check`.
  (The dataset path may sit under `/kaggle/input/datasets/manwithacat/...` with Kaggle's newer layout.)
- Locally: wheels for Linux x86_64 and Apple Silicon on the [GitHub releases](https://github.com/manwithacat/kickmix-engine/releases).
  `kickmix-engine doctor --json` tells you which kernel your CPU gets.

**Caveats.** invsqrt and DQI only; chemistry isn't supported yet. kickmix must still be importable, because the
official verifier modules import it. Speed varies with machine, circuit and salt (`--timing` shows yours). The
measured results and their caveats are in the README.

The code was written by an AI model (Claude) under my direction. It's an independent project, not affiliated with
Google, Qualtran or Kaggle, and it will be retired after the competition. Bug reports, especially any field that
differs under `--cross-check`, are very welcome on GitHub or here.


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**theredbluepill** · 2026-10-30 06:45 UTC · 0 votes

will try this out.
