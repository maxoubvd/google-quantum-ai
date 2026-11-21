# VQCC ThemBones milestone generators

This dataset holds the generator code behind the three circuits of team ThemBones' milestone set **s108** (public
leaderboard 31.67924, scored on the upload of 2026-11-20 04:24 UTC) in Google Quantum AI's Verifiable Quantum Circuit
Compilation competition. The circuits themselves are in the companion dataset `vqcc-thembones-milestone-set`, and the
notebook checks them by sha256. Licence: Apache-2.0 (`LICENSE`, `NOTICE`, `third_party__*`).

It is source code and small data files only. The competition data (kickmix, the verifiers and the task inputs) is
**not** included: every script expects kickmix at `<HOME>/qcc/data/kickmix`, taken from the competition's data page.

**Paths.** The scripts were written for the team's machines. Absolute home paths have been replaced by the literal
`<HOME>` (for example `<HOME>/qcc/ws/...`, `<HOME>/Quantum/tools/...`). Before you run anything, substitute your own
home directory:
`grep -rl '<HOME>' . | xargs sed -i "s#<HOME>#$HOME#g"`. The scripts also expect the team's folder layout
(`~/qcc/ws/<workstream>/...`, `~/qcc/.venv/bin/python`, the repo at `~/Quantum`). Only those path strings were changed.
No other line of code was edited.

**Where the folders go** (the drivers name these locations):

| this dataset | expected at |
|---|---|
| `*/tools/` (all three merged) | `~/Quantum/tools/` (`Q` in the scripts = `<HOME>/Quantum`) |
| `invsqrt/src_s18_srcM/` | `~/qcc/ws/invsqrt_fc/s18/srcM/` |
| `invsqrt/cfg/cfg_M.json`, `cfg_b312g.json` | `~/qcc/ws/invsqrt_fc/s18/cfg/` |
| `invsqrt/data/r59b0/` | `~/qcc/ws/invsqrt_fc/best/r59b0/` |
| `invsqrt/prog/q18_p7.json` | `~/qcc/ws/invsqrt_fc/s19/q18/q18_p7.json` |
| `dqi/kit/kv_unfused5/`, `dqi/kit/kv_unfused4_z2/` | `~/qcc/ws/dqi_zn/src/` (the line's kit copy); the `dqi_dq3_*.py` tools name the same tree at `~/qcc/ws/dqi_dq3/gen/zc2p/home/qcc/ws/dqi_zn/src/`; drivers that pass `R12_SOURCE=$D/pinvar/$N/home/qcc/ws/dqi_zn/src/kv_unfused5` (a run's pinned copy) take it too |
| `dqi/kit/z9C_oracle.json` | `~/qcc/ws/dqi_int/zptrace/z9C_oracle.json` (the chain's `B5_ZP_ORACLE`) |
| kickmix (competition data) | `~/qcc/data/kickmix/` |
| a Python venv | `~/qcc/.venv/bin/python` |

**What "the record" means below.** `[M]` = measured in the team's logs (the log is named), `[R]` = read in the named
file. Files under `*/record/` are the team's own build records (SOURCE.txt and logs), copied with the same path
redaction.

**Over-inclusion.** Each `tools/` folder was collected automatically: starting from the entry scripts, the collector
followed every import and every file name a script mentions, comments included. A folder therefore holds more files
than any single run needs. Scripts that only work on the team's machines (box sync and remote-run wrappers, card
start / stop, watchers, launch queues) were left out, so some references inside the drivers point to files that are
not here. Those references are listed under "Missing" for each task.

| task | circuit (sha256 of the set's file) | entry | byte-for-byte rebuild on record |
|---|---|---|---|
| invsqrt | `invsqrt.kmb` c639a045edc8... (b196lend) | `invsqrt/tools/f9_i1_stack.sh` (or `s5_inv_h1end_build.sh`) | **YES** (two independent runs) |
| DQI | file A 2c043a9afc55...; ticket r46zd20Aa 859ff092df1a... | `dqi/tools/sh/dqi_int_zc5b.sh` on the kit `dqi/kit/kv_unfused5`, then `dqi_int_p901pub3.sh`; tickets with `midprize_nb_reroll.py` | tickets from file A **YES**; the kit from `kv_unfused4_z2` **YES** (353/353 files); file A from its chain: **not re-run** |
| chemistry | 15 loaders + `submission.npz` e05fb18f808e... | `chem/tools/chem_c3_entry_raw.py` + `lead_gauge.py` | loaders from their `.kmx` **YES** (15/15); the `.kmx` and the table **NO** (pinned as files) |

## Version 2: what is now included, and what is still missing

- **Added in version 2 (2026-11-20): `dqi/kit/`** (708 files, 3.5 MB): the DQI kit the line runs on (`kv_unfused5/`),
  the tree it is made from (`kv_unfused4_z2/`), the zero-trace oracle `z9C_oracle.json`, and `KIT_IDENTITY.txt`
  (what is proven about them). Version 1 said these were only on the DQI build machine. That was wrong: they were on
  the packing machine one folder level deeper than its search went. See the DQI section for what is proven.
- **Still missing:** (1) the competition data: kickmix, the verifiers and the task inputs (take them from the
  competition's data page); (2) the DQI parent chain zd9's file A (693 MB), used only by the twelve-salt pair check,
  not by the build; (3) the team's machine-only wrappers (box sync, remote run, run pinning, card start / stop,
  watchers, queues), so some driver references point to absent files (listed per task); (4) the history behind the
  carried trees: invsqrt's `srcM` and the DQI kit are carried in their final state, and the patches that led there
  from the public sources are not replayed; (5) the chemistry's search seeds and the GPU fit race (the `.kmx` and the
  table are carried as files instead).

## Start here: one paragraph per task

**invsqrt.** Put the folders where the table above says. Then run the two commands of the invsqrt section below
(`tools/invsqrt_fc_s19stage.py`, then `bash tools/f9_i1_stack.sh q18 dq18c prog/dq18c.staged.json 3` with `VERIFY`
unset). It needs Python 3.13 with numpy, scipy and numba, kickmix, one CPU core and about 20 minutes. The output's
sha256 must be c639a045..., the set's `invsqrt.kmb` (b196lend). Rebuilt byte for byte: **yes, twice**, by the team's
builder seat on 2026-11-20 (runs dq18c and dd21a; `invsqrt/record/stack_dq18c.log`).

**DQI.** Put `dqi/kit/kv_unfused5/` and `dqi/kit/kv_unfused4_z2/` in `~/qcc/ws/dqi_zn/src/` and
`dqi/kit/z9C_oracle.json` in `~/qcc/ws/dqi_int/zptrace/`. Then run the chain command of the DQI section below
(`tools/sh/dqi_int_zc5b.sh` with its eleven module files, then `tools/sh/dqi_int_p901pub3.sh` for the re-rolls).
Make the ticket with `python -I tools/midprize_nb_reroll.py append <fileA.kmb> r46zd20Aa <out.kmb>`. It needs Python 3
with numpy and numba, kickmix, 2 cores for the chain and about 40 minutes. The pair check needs 2 more cores and the
parent's file A, which is not included. Rebuilt byte for byte: **ticket <-> file A** with the re-roll tool (the
team's record for set s088; the companion notebook checks s108's every time it runs), and **the kit `kv_unfused5`
from `kv_unfused4_z2`** with `tools/dqi_zn_mkunf5.py`, all 353 files (this dataset's packer, 2026-11-20 15:37 UTC).
**Not rebuilt:** file A from its chain (run once, no second run on record).

**Chemistry.** Merge `chem/tools/` into the tools folder. Then run
`RAW=table/all_c.npz python -I tools/chem_c3_entry_raw.py <BASE> <OUT> WORDS=table/asm_words.npz 0=loaders_kmx/angle_k_0.kmx ... 14=loaders_kmx/angle_k_14.kmx`
and then the gauge `tools/lead_gauge.py`, as in the chemistry section. It needs Python 3 with numpy and scipy, and
kickmix. The fit tools also need JAX and optax on a GPU. Rebuilt byte for byte: **the 15 loaders from their `.kmx`,
15/15** (this dataset's packer, 2026-11-20 14:21 UTC). **Not rebuilt:** the `.kmx` and the table (randomized searches
and a GPU race, carried as files).

---

## invsqrt (`invsqrt/`)

**What made the set's file** [R `invsqrt/record/b196lend_SOURCE.txt`]:
```
python -I tools/s5_inv_h1end_genpatch.py <copy of ~/qcc/ws/invsqrt_fc/s18/srcM>
bash tools/s5_inv_h1end_build.sh s5lend1 11    # FC's s19hold pipeline unchanged, S5LEND=1, cfg_b312g.json, R59's angles
```
**The rebuild on record.** `tools/f9_i1_stack.sh` is the same pipeline in one command. It copies the source
`src_s18_srcM` (FC's s18/srcM), applies the generator patch `s5_inv_h1end_genpatch.py` (q11 lent in H1's word un-add),
writes the staged q18 program, and then runs emit, width step, all-input check, P1 pass, post-chain, R91 b312 carry and
the exact all-input proof against b196lend:
```
python -I tools/invsqrt_fc_s19stage.py <PROG> prog/q18_p7.json prog/dq18c.staged.json   # staged = today's exit
bash tools/f9_i1_stack.sh q18 dq18c prog/dq18c.staged.json 3                          # VERIFY unset: local only
```
- Run dq18c (Nov 20, 04:53:43-05:11:53 UTC, one core) and run dd21a (the digit-21 path, 1,845 s) each produced
  sha256 c639a045... = b196lend, **byte for byte** [M `record/stack_dq18c.log`; the F9-I1 seat's log "ITEM 1"]. S5's
  own build made the file in the first place [R SOURCE.txt].
- `src_dq18c/` is the source tree as it stood after that run (srcM + the patch + the staged program), kept so it can
  be read and diffed. `cfg/cfg_dq18c.json` is its config: FC's `cfg_M.json` + the q18 stub key + cfg_b312g's
  "r91_splitdebt". `data/r59b0/` holds R59's router and angles. `prog/` holds the staged program and q18_p7.
- Needs: Python 3.13 (team venv 3.13.18), numpy 2.5.4, scipy 1.18.1; numba (imported by the generator source; the
  driver sets NUMBA_CACHE_DIR); kickmix; one CPU core. python-sat and ortools are imported only by 8 SAT / search tools
  in `tools/`, none of them in the source trees. Wall time 1,088 s on record
  (401-1,845 s depending on load), plus about 190 s a salt for the official verifier.
- Missing: the `VERIFY=1` branch of `f9_i1_stack.sh` calls the team's box lane (`box_sync.sh`, `box_verify.sh`, left
  out), so run the official verifier from the competition data yourself. The history of patches that produced srcM
  from the public generator is not replayed here: srcM is carried in its final state.

## DQI (`dqi/`)

**What made the set's file** [R `dqi/record/best_q31zc14_SOURCE.txt`, `..._repro901_zd20.txt`]: best_q31zc13's line
with this folder's `dqi_int_mon1.py` (sha256 e7bcbe33..., the line's own copy; the `tools/` copy here is that version)
and `DQI_INT_MON1_MERGE=1`, run as chain zd20 on a 2-core slot:
```
B5_INV17_MAP=dqi_int_inv17map3.py B5_ZC_PARENT=zd9 B5_ZC_PTAG=zd9A bash ~/Quantum/tools/sh/dqi_int_zc5b.sh <fresh> <c1,c2> <c1-c2> <pair cores> \
  dqi_dq3_kpn.py,dqi_dq3_inv24.py,dqi_dq3_om36.py,dqi_dq5_c12.py,dqi_int_div6k.py,dqi_dq5_c12n3.py,dqi_int_fdots.py,dqi_int_prek.py,r81_dqi_fold_k22_K154.txt,r81_dqi_fold_k23_K224.txt,dqi_int_mon1.py \
  dqi_dq3_kpn,dqi_dq3_inv24,dqi_dq3_om36,dqi_dq5_c12,dqi_int_div6k,dqi_dq5_c12n3,dqi_int_fdots,dqi_int_prek,dqi_int_mon1 \
  "DQI_INT_FDOTS_OUT=<json> DQI_INT_PREK_OUT=<json> DQI_INT_MHR_EARLYF=1 DQI_INT_MHR_EARLY_LEVELS=11,13 DQI_INT_MON1_OUT=<json> DQI_INT_MON1_MERGE=1"
bash ~/Quantum/tools/sh/dqi_int_p901pub3.sh <fresh> <cores> 2      # re-rolls
```
The fully expanded environment of the run, including all 29 extension modules of `DQI_INT_EXT`, is the first line of
`record/best_q31zc14_repro901_zd20.txt`. Every module it names is in `tools/`.
- Tickets: `python -I tools/midprize_nb_reroll.py append <fileA.kmb> r46zd20Aa <out.kmb>` gives 859ff092... and
  `... strip <ticket.kmb> <fileA.kmb>` gives file A 2c043a9a... (standard library only; the companion notebook does
  this and checks every sha256 against the verified ticket). **Byte for byte**: the same tool rebuilt the tickets of
  the earlier set s088 byte for byte, in both directions [M team record REPRODUCE_s088]. For s108 the notebook makes
  this check every time it runs.
- File A from its chain: run once (zd20, 01:01:30 UTC start, file A at 01:27:01, done at 01:40:40 [R repro901]). The
  re-rolls b, c and d were later remade byte for byte [M SOURCE.txt "RE-ROLLS"]. **No second run of the chain is on
  record**, so the rebuild of file A from its chain is not shown.
- Needs: Python 3 with numpy and numba (22 files import numba); kickmix; 2 cores for the chain plus 2 for the pair
  check; about 40 min wall time on record.
- **The kit and the oracle (INCLUDED since version 2, `dqi/kit/`; details and the commands to recheck in
  `dqi/kit/KIT_IDENTITY.txt`).** `kv_unfused5/` (353 files: 322 `.py`, 31 `.json`; 1.76 MB) is the team's working
  copy of the public DQI bundle credited in NOTICE, as the line runs it. `kv_unfused4_z2/` is the tree it is made
  from: `python3 tools/dqi_zn_mkunf5.py <copy of kv_unfused4_z2>` changes one file (`research/production/selected.py`:
  a14's reverse fusion is not installed, a12's reverse_root_lifetime is), and the result equals `kv_unfused5`, **all
  353 files byte for byte** [M 2026-11-20 15:37 UTC, the packer]. `z9C_oracle.json` (sha256 d08ab947...) is the
  zero-trace oracle the chain reads. What is proven: the packaged trees are the snapshot of generator line zc2p
  (pinned 2026-11-18 03:42 UTC), copied unchanged (no path needed redaction), and byte-identical file by file to the
  other copies on the packing machine (snapshots q15b and r37, the variant `kv_unfused5_t54cp`, a run's pinned copy).
  The team recorded hashes of the kit on 2026-11-14 ("py tree 404a19ae" in the kit builder's log, "py tree
  f744011163c23292" in the pin record of run q15b, "tree sha 18c06fb8" for kv_unfused4), and the oracle's hash
  (d08ab947, best_q27c9f's record). All four hashes come out the same on these files. **Not proven:** that the copies the zc14
  chain read on its build machine on 2026-11-19 were these bytes. No hash of them from that day is on record, and
  that machine is off and was not re-checked.
- **Still missing:** (1) kickmix and the verifiers (competition data). (2) The parent chain zd9's file A (693 MB),
  used only by the twelve-salt pair check. (3) The driver `dqi_int_pinvar.sh` / `pinvar.sh`, left out as
  machine-only: on the team's machines it wrote a run's pinned copy and, with `KIT_MK` set, made `kv_unfused5` inside
  it from `kv_unfused4_z2` with the same maker. Without it, point the drivers at your copy of `dqi/kit/kv_unfused5`
  (for example `R12_SOURCE=<path>` in the drivers that take it). (4) The steps that made `kv_unfused4_z2` from the
  public bundle are not replayed: the tree is carried as it was.

## Chemistry (`chem/`)

**What made the set's files** [R `chem/record/c4_t11r10ag2caa_g_SOURCE.txt`, `record/asm_c4_t11r10ag2caa_asm.log`]:
1. The table `t11_r10ag2ca`: a GPU fit race (JAX / optax; the race tickets are in `tools/sh/chem_t11_races.txt`),
   then the table step (rounding with refit, `chem_rk7_table7c.py` / `chem_wr7_table7c.py`; `table/fin.log`,
   `table/t7.log`). Result: `table/all_c.npz` (lambda_eff 46.933606, `table/pcheck.json`, layout `table/rel.txt`).
2. The 15 loaders were compiled by the loader port (`tools/sh/chem_c5_port.sh` -> `chem_c4_port4.sh` LAYOUT 10, plus
   C3's k0 pipeline for k0) into the `.kmx` in `loaders_kmx/` (`record/asm_c4_t11r10ag2caa_parts.txt` names the
   build folder of each loader).
3. Assembly: `RAW=table/all_c.npz python -I tools/chem_c3_entry_raw.py <BASE> <OUT> WORDS=table/asm_words.npz 0=loaders_kmx/angle_k_0.kmx ... 14=loaders_kmx/angle_k_14.kmx`
   (BASE lends only its folder layout: any folder with `dfthcsos/angle_k_*.kmb` and `submission.npz`).
4. The channel gauge `tools/lead_gauge.py`, which rotates only W and H2_I of `submission.npz` and leaves the loaders as
   they are.
- **Byte for byte:** all 15 `.kmx` written through kickmix's kmb writer give the set's 15 loaders, **15/15** [M
  2026-11-20 14:21 UTC, the packer's check: `angle_k_0..14.kmb` sha256 equal to the set's]. The gauged entry's
  loaders are byte-identical to the assembled entry's [M cmp].
- **Not rebuilt:** the `.kmx` come from randomized loader searches whose seeds and timing are not all on file. The
  table comes from a GPU fit and a rounding race whose seeds are timing-dependent and cannot be replayed bit for bit.
  The tools are here. The exact `.kmx` and the exact table are carried as files (`loaders_kmx/`, `table/all_c.npz`,
  `table/asm_words.npz`). The final `submission.npz` is in the companion set dataset.
- Needs: Python 3 with numpy and scipy; JAX and optax on a GPU for the fit (41 files import jax); numba for some loader
  tools; kickmix. On record the port ran on 20 cores for about 3 hours and the table step took about 4 hours.
- Missing: the box-run wrappers of the port and the table landings (`chem_t11_par*.sh`, `chem_t10_s6par.sh`,
  `chem_c3_chain.sh`, `*box*.sh`, remote / put / queue scripts) were left out as machine-only. The multi-GB fit dumps
  and the per-loader search folders are not included.

## Credits

See `NOTICE`. Upstream licences and credits are in the `third_party__*` files. Built by our team with the help of AI
coding agents (the competition's Description invites them).
