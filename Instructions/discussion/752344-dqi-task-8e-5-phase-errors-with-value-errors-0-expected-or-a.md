# dqi task: ~8e-5 phase errors with value_errors = 0 — expected, or a fixable implementation detail?

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752344 — copied 2026-11-15 (raw Markdown, unedited).

**xhsmsuper** · 2026-10-24 14:26 UTC · 2 votes

I'm seeing a small but consistent residual on the `dqi` task and I'd like to know whether it's inherent to the protocol or a fixable detail in my construction.**Observation** (official local verifier, 102400 shots):- `value_errors = 0` — the syndrome check always passes- `phase_errors ≈ 8e-5` — a handful of shots get the wrong phase- This makes `R = eps_plus/eps ≈ 1.2` instead of 1.0**What I've ruled out**:1. It is not the unique-decoding radius. For every failing shot I ran a classical Berlekamp-Massey decode: it returns degree = k with roots matching the true error positions exactly, so the input is uniquely decodable.2. It is not correlated with the error count. Grouping by k (from 27 to 35), the per-k failure rate is flat at ~8e-5.3. It is not a positional or value pattern — no boundary positions, no special values.**The decisive observation**: if I take one failing input and re-run it many times (varying only the internal measurement outcomes), roughly **50%** of the runs fail. For normal inputs, 0% fail. So there exist "special" inputs (roughly 1.6e-4 of all inputs), and for those the circuit's correctness depends on the intermediate measurement outcomes.Since the DQI contract should depend only on the input (the decoded error and the MBUC mask), a phase that varies with intermediate measurements looks like an unhandled measurement branch rather than an algorithmic limit.**My question**: is this expected behaviour of the approximate/truncated decoder path, or is there a known construction detail (e.g. a condition that must be applied when the measured register is zero) that removes it? If it's a known limitation, is there a published bound for the residual epsilon that I should target instead of driving it to zero?Happy to share more details on the diagnostic method (classical BM cross-check + repeated-measurement runs) if useful.


---

## Comments (6 threads)

Replies are indented with `>`, one level per reply depth.

**Tanuj Khattar** — **Host** · 2026-10-27 23:12 UTC · 2 votes

This indicates that your decoder (or its uncomputation) is producing the wrong outcome on those rare inputs.

When designing reversible quantum circuits, you have to make sure all ancilla scratch space is uncomputed. Doing an `HMR` on an ancilla qubit means you are shifting the uncomputation problem to the phase. `HMR` always resets the ancilla qubit to `|0>` (so `value_errors = 0`), and when the ancilla was `|1>`, you get a `-1` phase kickback whenever the measurement outcome is `1` (`H|1> = (|0> - |1>)/sqrt(2)`). If you don't uncompute the phase error correctly, you'll see 50% of the runs fail on the same input. Here is a minimal example to demonstrate the concept.

Suppose your input / output are four qubits `|q0>|q1>|q2>|q3>` and your goal is to apply the 4-controlled phase `(-1)^(q0 & q1 & q2 & q3) |q0>|q1>|q2>|q3>`. Since `kickmix` only has up to 2-control phase/Toffoli primitives (`Z`, `CZ`, `CCZ`, `CCX`), we compute `q0 & q1` into a scratch ancilla `q4`.
The reversible `kickmix` circuit uses **3 Toffolis** ([open in Quirk](https://algassert.com/quirk#circuit=%7B%22cols%22%3A%5B%5B%22H%22%2C%22H%22%2C%22H%22%2C%22H%22%5D%2C%5B%22Amps4%22%5D%2C%5B%5D%2C%5B%5D%2C%5B%5D%2C%5B%22%E2%80%A2%22%2C%22%E2%80%A2%22%2C1%2C1%2C%22X%22%5D%2C%5B1%2C1%2C%22%E2%80%A2%22%2C%22%E2%80%A2%22%2C%22Z%22%5D%2C%5B%22%E2%80%A2%22%2C%22%E2%80%A2%22%2C1%2C1%2C%22X%22%5D%2C%5B%22Amps4%22%5D%5D%7D)):
```text
CCX q0 q1 q4
CCZ q2 q3 q4
CCX q0 q1 q4
```
If you forget to include the last `CCX` that uncomputes `q4`, then you'd always see a value error on the 4 inputs on which `q0` and `q1` are both `1` ([open in Quirk](https://algassert.com/quirk#circuit=%7B%22cols%22%3A%5B%5B1%2C1%2C%22H%22%2C%22H%22%5D%2C%5B%22Amps4%22%5D%2C%5B%5D%2C%5B%5D%2C%5B%5D%2C%5B%22%E2%80%A2%22%2C%22%E2%80%A2%22%2C1%2C1%2C%22X%22%5D%2C%5B1%2C1%2C%22%E2%80%A2%22%2C%22%E2%80%A2%22%2C%22Z%22%5D%5D%2C%22init%22%3A%5B1%2C1%5D%7D)):
```text
CCX q0 q1 q4
CCZ q2 q3 q4
```

You can optimize the original circuit by replacing the uncomputing `CCX` with an `HMR` and a 0-Toffoli Clifford `CZ` phase correction ([open in Quirk](https://algassert.com/quirk#circuit=%7B%22cols%22%3A%5B%5B%22H%22%2C%22H%22%2C%22H%22%2C%22H%22%5D%2C%5B%22Amps4%22%5D%2C%5B%5D%2C%5B%5D%2C%5B%5D%2C%5B%22%E2%80%A2%22%2C%22%E2%80%A2%22%2C1%2C1%2C%22X%22%5D%2C%5B1%2C1%2C%22%E2%80%A2%22%2C%22%E2%80%A2%22%2C%22Z%22%5D%2C%5B1%2C1%2C1%2C1%2C%22H%22%5D%2C%5B1%2C1%2C1%2C1%2C%22Measure%22%5D%2C%5B%22%E2%80%A2%22%2C%22Z%22%2C1%2C1%2C%22%E2%80%A2%22%5D%2C%5B%22Amps4%22%5D%5D%7D)):
```text
CCX q0 q1 q4
CCZ q2 q3 q4
HMR q4 b0
CZ q0 q1 if b0
```
Now suppose you forgot the phase correction `CZ q0 q1 if b0` and instead ran ([open in Quirk](https://algassert.com/quirk#circuit=%7B%22cols%22%3A%5B%5B%22H%22%2C%22H%22%2C%22H%22%2C%22H%22%5D%2C%5B%22Amps4%22%5D%2C%5B%5D%2C%5B%5D%2C%5B%5D%2C%5B%22%E2%80%A2%22%2C%22%E2%80%A2%22%2C1%2C1%2C%22X%22%5D%2C%5B1%2C1%2C%22%E2%80%A2%22%2C%22%E2%80%A2%22%2C%22Z%22%5D%2C%5B1%2C1%2C1%2C1%2C%22H%22%5D%2C%5B1%2C1%2C1%2C1%2C%22Measure%22%5D%2C%5B%22Amps4%22%2C1%2C1%2C1%2C%22%E2%80%A2%22%5D%2C%5B%5D%2C%5B%5D%2C%5B%5D%2C%5B%22Amps4%22%2C1%2C1%2C1%2C%22%E2%97%A6%22%5D%5D%7D)):
```text
CCX q0 q1 q4
CCZ q2 q3 q4
HMR q4 b0
```
This circuit has zero value errors on every input (since `q0..q3` are untouched and `HMR` always resets `q4` to `|0>`). On the 4 inputs with `q0 & q1 = 1` (`|1100>`, `|1101>`, `|1110>`, `|1111>`), `q4 = 1` and applying an `HMR` causes a phase error 50% of the time, specifically whenever `HMR q4 b0` samples `b0 = 1` (you can also see those 4 amplitudes have wrong phase in the 4th quirk link whenever measurement outcome is 1).

**xhsmsuper** — **Topic author** · 2026-10-28 02:29 UTC · 0 votes

Thank you - this is exactly the confirmation we were hoping for from the organisers, and the minimal Quirk example makes the mechanism unambiguous.

The equivalence you point at is the piece we had independently measured and can now name precisely. Running single- and double-HMR probes on |0> and |1> in kickmix, we found the phase kickback is

    phase ^= b & q      (b = classical outcome bit, q = ancilla value before HMR)

i.e. a Pauli-Z kicked back onto the control pattern exactly when the measured outcome is 1 - which is your `CZ q0 q1 if b0` line, generalised to whatever the ancilla was controlling. That rule reproduces both of your examples: with the correction the 4-input phase is exact; without it the four inputs with q0&q1=1 acquire a phase that depends on b0, hence the 50%-of-runs failure we observed for one special input (128/256 measurement streams wrong, 0/256 on inputs without the pattern).

Two things we would like your steer on, since they decide how much it is worth chasing to zero:

1. **Is R(eps+) driven purely by the measured residual rate eps+, or does the verifier also reward pe = 0 beyond that?** Our reading of the metric is that R is a fixed function of eps+, so a construction with pe = 0 and one with pe = 3 out of 102,400 differ by ~1.5% in final cost and nothing more. If that is right, the practical question is whether the ~1.6e-4 rare-input set can be enumerated and covered by explicit conditional corrections (classically-controlled Cliffords, zero Toffoli cost) rather than removed structurally.

2. **On the uncomputation accounting**: is a conditional CZ gated on an outcome bit billed at zero Toffoli, as our reading of the cost model says? If the full set of rare-branch corrections costs qubits/depth but no Toffoli, covering them seems almost free - and it would also answer whether the current frontier constructions (Q=889, T~474K, pe 0..4 across 20 salts) have room to both reach pe = 0 AND keep the tighter qubit budget, or whether those two goals trade off.

For context: our current dqi runs at Q=889 / T~474K with pe 0..4 depending on the salt, so the residual is already tiny in the eps+ sense; if the remaining gap is purely eps+-driven it is worth ~1% and we would rather spend the effort on peak qubits (which enter as Q^2). If instead there is a construction detail that makes pe exactly 0 at no cost, that is a different decision.

Happy to share the probe scripts (HMR semantics sweep + repeated-measurement-stream harness) if useful to anyone else hitting this.

**xhsmsuper** — **Topic author** · 2026-10-27 06:38 UTC · 0 votes

Nice - the fixed-point / gauge-drift framing over ~140 measurement rounds is a better model than my single-branch hypothesis, and it predicts exactly what we observe on that same lineage: our v3-class build (445,229 T, Q=1,789) sits at R ~ 1.21 with pe in {0..4} across 20 salts (mean 0.9). That reads as O(1/2^margin) boundary noise rather than a deterministic uncomputation miss, which matches your reading.

One quantitative note you may find useful when calibrating the margin sweep. On the newer 889-qubit Taylor / J-fraction generation the implied prefactor is not R(pe) alone. Measured with the official local verifier:

  pe = 0  ->  cost / (Q^2 * T) = 1.1913
  pe = 2  ->  cost / (Q^2 * T) = 1.2014   (delta = +0.85%)

The pe=0 value lands exactly on the paper's R floor for eps+ ~ 8.7e-4, i.e. R(eps+) ~ 1.19. So once a construction is phase-clean, the remaining ~19% multiplier looks like a property of the 102,400-shot test budget rather than of the circuit: driving residual phase errors below ~1e-5 buys at most another ~0.85%. On this reading your clean margin-1 split at 0/102,400 is already at the floor, and squeezing margin further only costs qubits.

If you are willing to share the verify log format (102,400 shots, batch 32, workers 1) we would like to cross-check our 20-salt harness against it. Ours derives a per-batch rng_seed from a SeedSequence and we have never confirmed that this matches the official batch schedule - if it does not, our pe in {0..4} band would need re-reading. Happy to share the salt runner and the BM cross-check script in return.

**shawncsx** · 2026-10-27 00:40 UTC · 0 votes

Same residual signature on our local split. Two observations from the v3-444k-T build we ran (3-5 phase errors, value_errors=0):

- The residual is flat across k (no k-specific structure) and depends on the intermediate measurement outcomes — consistent with a fixed-point / gauge drift accumulating over the ~140 measurement rounds rather than a bug in any single decoder round.
- The margin-1 parameterization from the 6.17 release (same 456,568-T topology, different rounding of the fixed-point constants) runs clean on our split: 0 phase errors out of 102,400 shots. So the residual looks like an O(1/2^margin) boundary effect, not a systematic phase error in the decoder itself.

Happy to share our verify log format if useful — 102,400 shots, batch-size 32, max-workers 1 (memory bound on the local validator).

**xhsmsuper** — **Topic author** · 2026-10-26 09:51 UTC · 0 votes

Thank you for the detailed confirmation - the margin-1 framing matches what we suspected (an unhandled rare branch folding into the typical path). We reproduced the same signature independently: classical BM always uniquely decodes the failing syndromes, per-k failure rate flat at ~8e-5, and one failing input re-run across 256 measurement streams fails ~50% of runs while normal inputs are 0/256. Your quantified expectation (special-input rate ~1.6e-4 = expected failure rate) closes the loop nicely.

Two notes from our side:

1. On the cost tradeoff: with score = R * Q^2 * T, an extra measurement qubit costs ~0.11% (at Q~881) while dropping margin and losing R 1.00 -> 1.21 costs ~21%. So margin-1 is clearly correct at the current frontier - but the frontier itself has moved: the current best dqi constructions run at Q=881, T=964K with phase_errors=0 (dqi ratio 27.4x), i.e. both phase-clean AND half the qubits, which suggests that bank structure (J-fraction style) handles the rare branch without a dedicated margin bit. Would be interesting to compare construction notes.

2. Our diagnostics (classical BM cross-check + repeated-measurement stream method) are described in the topic; happy to share script details if useful for your margin-1 validation.

**shawncsx** · 2026-10-26 01:55 UTC · 0 votes

We hit the same residual while building our 6.17 margin-1 DQI decoder: value_errors = 0, a handful of phase errors out of 102,400 shots. Same signature as yours - flat across k, independent of the classical BM decode (which always uniquely decodes), and depends on the intermediate measurement outcomes for a ~1.6e-4 fraction of inputs. Your 50%-of-runs-fail observation on the special inputs is the decisive clue: it points at an uncomputation branch that only fires when a measured register lands in a rare pattern, so the phase the circuit leaves behind depends on which measurement outcome actually happened.

Practical fix we verified locally: the margin-1 parameterization is exactly what pins this down - it reserves one extra measurement bit so the rare branch is always handled explicitly instead of being folded into the typical path. With margin-1 the residual is 0 phase errors on our split; the ~8e-5 you're seeing matches the expected 1/1.6e-4 rate of the unhandled rare inputs.

If you want a second data point: the lightning community build (444,030 T, 1,842 qubits) also shows 3-5 phase errors on its local verify, and adding and_implication_fast makes it worse (4 errors, repetition_factor 1.2095). The 6.17 margin-1 build (456,568 T, 1,856 qubits, 0 phase errors) is the only configuration we found that is phase-clean at the verifier's shot budget.
