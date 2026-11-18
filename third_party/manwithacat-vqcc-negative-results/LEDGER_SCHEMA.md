# VQCC Negative Results ledger: schema

`ledger.csv` (UTF-8, RFC 4180 quoting, one row per idea) and `ledger.parquet` (same rows and columns) record ideas
tried on the three tasks of the Google Quantum AI verifiable circuit compilation competition, and how they failed.
Rows come from one team's preregistered work cards (`source = ours`), from the same team's second, ChatGPT-based
reference-research pathway (`source = R`), and from other entrants' public posts and notebooks
(`source = <Kaggle name>`, always credited and linked). It is a research log, not a submission artefact:
no circuits, no scores claimed.

Rows are append-only. IDs are stable (`NR-0001`, `NR-0002`, ...); a row is never renumbered or deleted. A correction
is a new row whose `supersedes` names the old ID. When a later public circuit or a later result of ours shows that a
row's conclusion was wrong, or broader than its evidence, the row gets a post-hoc annotation in the four `later_*`
columns; its original fields are never edited (see "Later findings" below).

## Columns

| Column | Type | Allowed values / format | Meaning |
| --- | --- | --- | --- |
| `id` | string | `NR-` + 4 digits | Stable identifier |
| `task` | string | `invsqrt`, `dqi`, `chemistry`, `all` | Competition task (`chemistry` = `dfthcsos`); `all` for tooling or cross-task ideas |
| `mechanism` | string | short phrase | What was tried, e.g. "drop the BM gap flag", "7-bit angle words with refit" |
| `family_target` | string | free text | Circuit family, public baseline or component targeted (e.g. "public 63-q lookup family", "hengck23-style BM decoder", "R=1 DF-THC factorization") |
| `date_tested` | date | `YYYY-MM-DD` | Date the result was recorded (card log date, or post date) |
| `source` | string | `ours`, `R` or a Kaggle display name | Who did it (`R` = the team's ChatGPT research pathway) |
| `link` | string | URL or empty | Public notebook / discussion / dataset / DOI where the evidence can be read; empty if none is public |
| `evidence_tier` | string | see below | Strongest evidence behind the outcome |
| `outcome` | string | `kill`, `null`, `partial`, `open idea`, `observation` | See below |
| `key_numbers` | string | free text | The numbers that decided the outcome, with units |
| `why_it_failed` | string | free text | The structural or measured reason, in one or two sentences |
| `reopen_if` | string | free text | The condition under which the idea is worth revisiting (from the card's next decision or the author's own words) |
| `related_cards` | string | `;`-separated | Work-card file names (ours) or related ledger IDs |
| `dataset_files` | string | `;`-separated | Raw evidence files shipped in this dataset (paths relative to the dataset root), if any |
| `keywords` | string | `;`-separated, lower case | Search terms (e.g. `qubits;uncompute;bm;gf-multiply`) |
| `licence` | string | e.g. `Apache-2.0` | Licence of the underlying evidence; `author's post (Kaggle)` when only a forum post |
| `added_utc` | string | ISO 8601 UTC | When the row was added to the ledger |
| `supersedes` | string | ledger ID or empty | Row this one corrects |
| `later_status` | string | empty, `contradicted`, `qualified`, `superseded` | Post-hoc annotation; empty = the row stands. See "Later findings" |
| `later_by` | string | free text | What showed it: public circuit SHA-256 prefix (12 hex) + author + notebook URL, or our card / ledger ID |
| `later_date` | date | `YYYY-MM-DD` or empty | When the later evidence appeared (not before `date_tested`) |
| `later_note` | string | free text | One or two sentences: what the row got wrong, and the lesson (e.g. which model assumption failed) |

### `outcome`

- `kill`: a preregistered kill/stop condition was met, or a structural argument rules the idea out.
- `null`: tried and measured; no gain (or a gain within noise), no structural proof that none exists.
- `partial`: some gain or a positive sub-result, but the line was closed short of its goal.
- `open idea`: proposed or untested (typically from the forum); listed so it is not lost or duplicated.
- `observation`: a measured fact about the circuits or the task that constrains future work but was not an idea
  tried (e.g. inputs that every decoder of a line fails). For these rows `why_it_failed` states why the fact
  matters, and `reopen_if` the finding that would lift the constraint. Supporting data sits under `observations/`.

### `evidence_tier`

- `official verifier` — the unmodified official verifier (local accepted build or hosted Kaggle) ran on a circuit.
- `measured by us` — our local analysis using the official equations or exact algebra, no verifier run.
- `local prescreen` — the project's own name (from its evidence vocabulary) for a local check with the official
  equations, a native simulator or exact models, no official verifier run. It is the same tier as `measured by us`;
  rows from 2026-11-09 onward use this name.
- `static model` — arithmetic projection or resource model only.
- `paper only` — a written argument or design study, nothing run.
- `author-reported` — another entrant's own claim, not re-verified by us.
- `observed public score` — read from the Kaggle leaderboard.

### Later findings (`later_*`, added 2026-11-17)

A negative result can be overturned. These columns record that without rewriting the original finding:

- `contradicted`: a later public circuit or a later result of ours shows the conclusion was wrong for the scope the row
  itself states (e.g. a static model priced a family at 3x and a public circuit of that family later beat the pick).
- `qualified`: the conclusion is still true for its exact scope (one file, one schedule, one grammar), but its wording
  or how it was used was broader than the evidence, and that misled later work.
- `superseded`: replaced by a newer row; `later_by` names that row's ID.

Blank means no later evidence against the row is known; it is not a re-verification. The annotation is set with
`negative_results_add.py later --id NR-.... --status ... --by ... --date ... --note ...`, which refuses to change any
other field. Rows with a `later_status` are still counted under their original outcome.

## Normalisation at staging

The master ledger is append-only, so rows written with non-schema spellings are not rewritten in place. Instead the
staging tool (`negative_results_add.py stage`) applies this explicit map to the published copy, and `add` refuses
the raw spellings for new rows:

| Column | Raw value in the master | Published value |
| --- | --- | --- |
| `task` | `dfthcsos` | `chemistry` |
| `task` | `invsqrt;chemistry` | `all` |
| `outcome` | `negative` | `kill` |
| `evidence_tier` | `static model / projection` | `static model` |
| `evidence_tier` | `static model / projection (census: local prescreen)` | `static model` (the lower tier) |
| `source` | `manwithacat (RezOps)`, `ours (census) on R's child` | `ours` |
| `licence` | `ours`, or empty on an `ours`/`R` row | `n/a (own work)` |

It also clears a `link` that is a repository path rather than a public URL (the path moves to `related_cards`),
joins `keywords` with `;` in lower case, and drops `dataset_files` entries that are not shipped in the dataset. A
stage may hold back named rows (`--exclude-ids`); an ID missing from a published `ledger.csv` is such a row, not a
deleted one. `negative_results_add.py lint` lists what staging will change.

## Querying

```python
import pandas as pd
L = pd.read_csv("/kaggle/input/vqcc-chemistry-negative-results/ledger.csv")
# Has X been tried on DQI?
L[(L.task == "dqi") & L[["mechanism", "keywords", "why_it_failed"]].apply(
    lambda c: c.str.contains("gap flag", case=False, na=False)).any(axis=1)]
```

New rows are added with `scripts/negative_results_add.py` (validates columns, assigns the next ID, refuses
duplicates on `(task, source, normalised mechanism)` and on an identical `link` + `mechanism`).
