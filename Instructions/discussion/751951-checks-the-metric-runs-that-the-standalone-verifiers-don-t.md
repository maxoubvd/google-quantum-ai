# Checks the metric runs that the standalone verifiers don't

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/751951 — copied 2026-11-15 (raw Markdown, unedited).

**Elias T.-H. Shieh** · 2026-10-23 06:50 UTC · 4 votes

Going through the metric notebook, I noticed it checks a few things before it runs any verifier, and the standalone verifiers never look at them. If one fails, the whole submission is rejected, so a local `"status": "PASS"` doesn't rule them out.

The one I'd watch most is the id cap: each circuit's largest qubit or bit id, `max(num_qubits, num_bits) - 1`, has to be below 699,050 (`MAX_QUBIT_OR_BIT_ID = 2**30 // 1536`). kickmix will parse much larger ids, and none of the verifiers check this.

The others I found:

- exactly one `invsqrt` file, one `dqi` file and one `dfthcsos/` folder, searched recursively through the zip, so a leftover backup copy anywhere counts;
- `angle_k_0` … `angle_k_<n-1>` with no gaps;
- every register needs a name;
- `submission.npz` arrays have to be floats (`tree_edges` an integer type) and little-endian. The dfthcsos verifier accepts an integer `H1_I = 0` or a big-endian array; the metric doesn't.

Also worth knowing: each benchmark is verified once per leaderboard split, each with its own secret salt, so a local run sees different tests.

A question for the hosts: the comment next to `MAX_QUBIT_OR_BIT_ID` says it can be raised if a benchmark needs more. Should we treat 699,049 as the fixed maximum for the rest of the competition?

I wrote these up as a [notebook](https://www.kaggle.com/code/b05204047elias/vqcc-submission-pitfalls-what-the-metric-checks) that runs the checks on a submission without simulating, in case it saves anyone a submission. If I've misread the metric, or you've run into another metric-only rejection, please say so.


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**Виктор Донькевич** · 2026-11-01 21:40 UTC · 0 votes

здравраствуйте, (верификация - затем подтверждение, закончится тогда, когда некому и нечего будет добавить,  П.С - удовлетворяет, работает, улучшить - добавить: улучшают)
