#!/bin/bash
# chem_own_build.sh: ONE COMMAND from a table of angle words to a complete, verified chemistry entry.
#
#   Inside WSL:      bash tools/sh/chem_own_build.sh <table.npz> <label> [options of tools/chem_own_entry.py]
#   From PowerShell: wsl -d Ubuntu -e bash -lc "export ARGS='<table.npz> <label> [options]'; tr -d '\r' <
#                    <HOME>/Quantum/tools/sh/chem_own_build.sh | bash"
#                    (add LOG=<name> to run it in the background with output in ~/qcc/ws/chem_own/logs/<name>.log)
#
#   <table.npz>  angle_words, tree_edges, W, H2_I, U  (the words and tree are the table; W gives the shape (R,B,C)
#                and, unless --weights scratch, the starting weights of the refit)
#   <label>      the entry is published to ~/qcc/ws/chem_own/best_<label>/ with MANIFEST.json (sha256 of every file)
#                and both official verifier summaries (salt "" and salt alt-salt-1)
#   options      --weights refit|scratch|keep   --effort N   --seed S   --rounds N   --workers 3   --quick-refit
#                --relax (last-block don't-cares)   --spec plan_spec.json (seed, rounds, relaxed-from round per
#                circuit, written by tools/chem_own_pick.py; tools/sh/chem_own_publish.sh = pick + this script)
#
# Pipeline: chem_own_refit.py (weights) -> chem_own_compile.py search (plans) -> chem_own_compile.py build
# (15 circuits, every address checked) -> official verifier x 2 -> publish.  At most 3 processes, niced;
# CPUS=5,7 pins everything to those cores (Level B: two pinned cores).
set -e
T=<HOME>/Quantum/tools
W=~/qcc/ws/chem_own
if [ $# -gt 0 ]; then ARGS="$*"; fi
if [ -z "$ARGS" ]; then sed -n '2,17p' "$0" 2>/dev/null; echo "usage: chem_own_build.sh <table.npz> <label> [options]"; exit 2; fi
mkdir -p $W/logs $W/build $W/tmp
cd $W
if [ ! -x venv/bin/python ]; then
  ~/.local/bin/uv venv venv --python 3.13 && ~/.local/bin/uv pip install --python venv/bin/python numpy scipy numba
fi
export NUMBA_CACHE_DIR=$W/tmp/numba_cache
# one BLAS thread: the refit is then bit-for-bit repeatable (thread count changes summation order)
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=${THREADS:-1} OPENBLAS_NUM_THREADS=${THREADS:-1} MKL_NUM_THREADS=${THREADS:-1}
if [ -n "$LOG" ]; then
  setsid nohup ${CPUS:+taskset -c $CPUS} nice -n 5 venv/bin/python $T/chem_own_entry.py $ARGS > logs/$LOG.log 2>&1 < /dev/null &
  echo "started: chem_own_entry.py $ARGS -> $W/logs/$LOG.log (pid $!)"
  sleep 3
else
  ${CPUS:+taskset -c $CPUS} nice -n 5 venv/bin/python $T/chem_own_entry.py $ARGS
fi
