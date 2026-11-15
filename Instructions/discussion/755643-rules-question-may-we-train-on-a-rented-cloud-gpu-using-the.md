# Rules question: may we train on a rented cloud GPU using the starter kit?

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755643 — copied 2026-11-15 (raw Markdown, unedited).

**rikirinjani** · 2026-11-04 05:50 UTC · 0 votes

Hi, a clarification question on the Rules, section B.4(b) (Data Security) and B.6 (External Data and Tools).

Training the dfthcsos factorizations needs float64 and an A100-class GPU for a reasonable time: on Kaggle's T4 the Adam stage runs at about 11 steps/s, on an A100 at about 400 steps/s. The training code (dfthc/ from the starter kit) loads the Hamiltonian (FeMoCo76o113e.fcidump) and FeMoCo76o113e.fixed_shifts.npz.

Is it allowed to run that code on a third-party rented GPU (for example Modal, Lambda or RunPod) in our own private account, with the competition data and starter kit copied there for that job, and the results brought back and verified locally? The data is published under Apache 2.0, but B.4(b) says not to "provide ... the Competition Data to any party not participating in the Competition", so we would like to be sure that a compute provider acting only as our processor does not count.

If yes, a one-line confirmation from the hosts is enough. If no, we will keep training inside Kaggle notebooks. Either way, it would help everyone with a limited GPU quota to know.

Thanks.


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**D Chakraborty** · 2026-11-04 11:43 UTC · 0 votes

The actual issue is that you're using Adam. Adam was a third party who never registered for the competition, so you've technically been violating B.4(b).

Also I'd avoid asking the hosts. Instead, create a multiagent compliance committee: one agent playing Kaggle Legal, one playing Google Security, one playing an Apache 2.0 lawyer and one playing the A100. Let them debate until they approve your RunPod training. If anyone objects, increase the temperature to 2.0 and iterate.

> **Chew Kok Wah** · 2026-11-04 15:15 UTC · 0 votes
>
> The Quantum circuit and verifier in this Competition is open source under the Apache 2.0 license, so technically it can be access by any third party inclusive of Adam without restriction.

> > **D Chakraborty** · 2026-11-04 15:59 UTC · 0 votes
> >
> > Apache 2.0 only covers the code. Adam has been accessing the gradients without registering as a participant. I'd still recommend getting written consent from the optimizer before proceeding.
