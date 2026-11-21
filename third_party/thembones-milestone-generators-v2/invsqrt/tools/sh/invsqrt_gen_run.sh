#!/bin/bash
# invsqrt_gen helper: run one of our python tools inside WSL with the 3.13 venv, niced, at most 3 threads.
# usage (from PowerShell):
#   wsl -d Ubuntu -e bash -lc "tr -d '\r' < /mnt/c/.../tools/sh/invsqrt_gen_run.sh | bash -s -- <tool.py> [args]"
export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3
T=<HOME>/Quantum/tools
cd ~/qcc/ws/invsqrt_gen || exit 1
tool="$1"; shift
exec nice -n 5 ~/qcc/.venv/bin/python "$T/$tool" "$@"
