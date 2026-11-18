#!/usr/bin/env bash
# Rebuild all three circuits from source into one folder, ready for tools/package.py.
#
#   bash build/build_all.sh OUT_DIR [TASK...]      # TASK in: invsqrt dqi dfthcsos (default: all three)
#
# Run inside WSL. Current base: Syed Asad Ali's v12 builders (Apache-2.0, see NOTICE and
# third_party/PROVENANCE.md); each build checks the SHA-256 of its output. Needs CPython 3.13 with
# numpy 2.2.4 (venv ~/venvs/gqai313, override with PY). Work trees go to ~/work/build (Linux
# filesystem, much faster than /mnt/c), override with WORK_ROOT.
# Measured (12 cores, WSL): dfthcsos 5 s, invsqrt 69 s, dqi 103 s; peak memory 1.0 GB (dqi).
set -euo pipefail
[ $# -ge 1 ] || { sed -n '2,10p' "$0"; exit 2; }
OUT=$1; shift
TASKS=${*:-dfthcsos invsqrt dqi}
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SRC=$REPO/third_party/syed3000-vqcc-own-circuit-portfolio-v12
ENG=$REPO/google-verifiable-quantum-circuit-compilation/kickmix
VER=$REPO/google-verifiable-quantum-circuit-compilation/verifiers
export PY=${PY:-$HOME/venvs/gqai313/bin/python}
WORK_ROOT=${WORK_ROOT:-$HOME/work/build}
mkdir -p "$OUT"
for t in $TASKS; do
  case $t in invsqrt|dqi|dfthcsos) ;; *) echo "unknown task $t" >&2; exit 2;; esac
  export WORK=$WORK_ROOT/$t
  rm -rf "$WORK"; mkdir -p "$WORK"
  echo "== $t"
  T0=$SECONDS
  bash "$SRC/$t/${t}_source/build.sh" "$ENG" "$VER" "$OUT" > "$WORK_ROOT/$t.log" 2>&1 \
    || { tail -20 "$WORK_ROOT/$t.log" >&2; echo "build of $t failed, log $WORK_ROOT/$t.log" >&2; exit 1; }
  grep "^BUILT" "$WORK_ROOT/$t.log"
  echo "   $((SECONDS - T0)) s"
done
