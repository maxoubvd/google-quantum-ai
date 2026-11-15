# This is a machine-learning competition, not a quantum-algorithms competition (maybe): where each of the three tasks stands

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752897 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-26 17:52 UTC · 5 votes

I chatted to a few people planning to win this by discovering a new
quantum algorithm for one of the three tasks. That may happen, but it is not
what the competition was built for, and the leaderboard so far says it is not
what is winning. A few things worth saying plainly, with sources.

**What the organizers designed.** Their paper is titled *Leveraging verifiable
structure in quantum circuits to find more efficient compilations with agentic
AI* (arXiv:2610.23832). The whole point of the setup is that every circuit is
reversible classical logic plus phases plus measurement-based uncomputation, so
a classical simulator can fuzz-test it exactly at billions of gates per second.
That makes circuit optimization an iterative loop: change the generator, build,
verify, measure, keep or discard. The paper cites the ecdsa.fail precedent, where
a crowd without quantum backgrounds halved Toffoli counts using AI coding agents
against a verifier. The "quantum" part of each task was fixed by the organizers
in the reference design. What is being optimized is reversible arithmetic and,
for chemistry, a classical tensor factorization.

**What the top of the board actually did.** Everything at the public frontier
is a verifier-in-the-loop story:

- The best public DQI decoder (jday96314's continued-fraction Schur decoder on
  GitHub, used in Alex Chilton's ≈ 19.2 notebook) came out of round six of an
  agent-driven search. A paper from the Lukin group (arXiv:2610.22367,
  *Discovering state-of-the-art quantum circuits with large language models*)
  independently produced a Reed–Solomon DQI decoder with the same structural
  ideas by LLM evolution against an evaluator.
- The chemistry gauge rotation in that notebook was found by Gemini editing one
  function in a loop against the official verifier. The chemistry factorizations
  themselves (8-bit words, (R, B, C) = (10, 129, 60)) are quantization-aware
  training runs on GPUs with the organizers' starter-kit trainer.
- mrmeows' 63-qubit invsqrt (down from 78) is a chain of storage and
  uncomputation rewrites, each checked on all 3.1 million legal inputs. The
  polynomial and the rotation angles never changed.

So the useful mental model is: the verifier is your loss function, the circuit
generator is your model, a fast local verifier is your GPU, and the search is run
by you, an agent, or an optimizer. Throughput of that loop is the edge.

**Where each task stands** (public authors' measurements; I have re-verified the
earlier plateau circuits locally and the arithmetic matches Kaggle to five
decimals).

*invsqrt: close to the floor of its architecture.* 63 qubits, about 1,400
Toffolis, 19.1× the reference. The error budget is 95–99% used, the authors
report Toffolis bounded below around 1,400–1,450 in this design, and Newton
refinement, cubic segments, other table widths and register recolouring all came
out neutral or worse. One qubit is worth about 1.5%. The one open idea on record
is computing phase bits directly as low-degree functions of x. Expect low twenties,
not thirties.

*DQI: the biggest mover and still the most structural headroom.* The cost is
R·Q²·T, so qubits count twice and 840 of them are the syndrome bank you cannot
shrink. The gains came from decoder architecture (Berlekamp–Massey → continued
fractions → Schur), from making the Forney sweep phase-only (11 Toffolis per
point instead of 80–85), and from shaving width (892 → 883). The sweep is now at
its floor; what is left is the forward-and-reverse extraction, where the reverse
costs about 31 Toffolis per node pair more than the forward. The GF(2¹²)
multiplier everyone uses is at the best known rank (42). Nothing here needs new
physics; it needs careful compute–use–reverse scheduling, which is precisely
what the agent loop finds.

*Chemistry: the least optimized task and the one that is genuinely machine
learning.* The best public is 7.7× at 415 qubits, versus 19× and 49× on the
other two, and under a geometric mean a point here is worth the most. The
levers are all training problems: the cost model has a workspace cliff at 8-bit
angle words (temporary qubits roughly halve) and another at 7 bits; plain
rounding of a 10-bit fit to 8 bits blows the 20 mHa accuracy gate (≈ 111 mHa),
so the factorization has to be refit quantization-aware, with a dictionary big
enough to recover accuracy and small enough to keep the lookups cheap. The
public authors' own next steps are lattice rounding on the linearised
sensitivity, annealed quantization-aware training and smaller 8-bit dictionaries
(B = 255 currently stalls at 48 mHa). If you have GPU hours and know how to
train under quantization constraints, this is where they pay.

**Practical notes.** Verify locally; a Kaggle CPU round trip is 19 minutes and
caps you at a handful of candidates a day. The official kickmix builds natively
on macOS (I published a build with NEON types and prebuilt wheels at
github.com/manwithacat/kickmix-macos) and on Linux from the Qualtran source.
Compare DQI decoders at the full 102,400 shots, since R(ε) at low shot counts is
dominated by the confidence bound. Keep the credit chains: mrmeows' open
circuits, Gluzdov's, Özer's, Phuc Thai's and jday96314's releases are why the
plateau is public at all.

**What I tried that did not work**, so you can skip it: single-parameter tweaks
to the reference circuits (all under 2%); gradient or gauge steps on public
chemistry factors with fixed angle words (the gauge bound is about 0.5%);
shrinking the chemistry dictionary at 10–11-bit words and R = 1 (pruned shapes
stall at 20.5–22 mHa; the frontier crossed the 8-bit cliff instead); adding any
extra register bank to a DQI decoder (under Q² it never pays). Full ledger in my
notebooks if anyone wants the numbers.

I am not affiliated with Google, Qualtran or Kaggle. Corrections welcome,
especially from the authors named above.
