# Run-time optimization rule

- **Date:** 2026-11-17 20:53 UTC
- **Plan step:** none (working method), during S1.5
- **Files:** `CLAUDE.md` (new rule "Optimize run time"), `Analysis/PLAN.md` (section 5, "Speed is a
  deliverable")

## What was done

- `CLAUDE.md`: our code must minimize run time: best algorithmic complexity first (complexity stated
  in the docstring of hot paths), NumPy or bit-parallel vectorization, caching, the 12 cores via
  `multiprocessing` with `fork` in WSL (sized to the 8 GB cap), profiling before tuning, measured run
  time and peak memory recorded in each changelog entry, and the fastest sufficient check during
  iteration (full official verifiers only at gates).
- `Analysis/PLAN.md`, section 5: run time is part of every tool's "done when".

## Why

The user asked for code optimized for compile and run time (computational complexity). The S1.5
verification showed why: a full invsqrt run takes about 1 h, while the reduced run takes about 1 min
and the generator's own exact proof already covers correctness.

## Results and analysis

Measured so far in S1.5 (WSL, 12 cores): dfthcsos verifier 26 s (953 MB peak); dqi verifier at
`--max-workers 2` 1 min 50 s per salt (4.7 GB peak), so 4 workers would need the WSL cap raised
above 8 GB.

## Next

Finish S1.5 (reduced invsqrt run, 7 extra dqi salts).
