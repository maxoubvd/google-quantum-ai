#!/usr/bin/env bash
# New file added by the team (Maxime Bendavid), 2026-11-21, plan step S2.26.
# Rebuilds the board dqi file 4b572141 (Tony Li / tonylica, on Team ThemBones' designs, upstream Syed Asad Ali) with
# tonylica's own builder (dataset tonylica/quantum-ai-32-85306-task-generators, Apache-2.0), which checks the digest
# itself. The ThemBones dqi chain cannot run as published (three inputs are missing: idle_l9r4dp_902.txt,
# spec_k1.txt, jt_tensors_r54c.json). Run inside WSL:  bash generators/thembones/build_dqi_tonylica.sh SRC_DIR OUT_DIR
# SRC_DIR (the unpacked dataset) is copied to ext4 first ($WORK, default ~/tl); the builder runs from that copy.
# Needs CPython 3.13 + numpy 2.5.4 (~/venvs/tb313) and the competition kickmix.so (sha256 c66b3216..., the version
# tonylica tested). About 1.5 GiB of memory and 2 GB of disk.
set -eu
R=$(cd "$(dirname "$0")/../.." && pwd)
SRC=$(readlink -f "$1"); OUT=$(readlink -m "$2"); WORK=${WORK:-$HOME/tl}
rm -rf "$WORK"; mkdir -p "$WORK" "$OUT"; cp -r "$SRC/." "$WORK/"
export PYTHONPATH=$R/google-verifiable-quantum-circuit-compilation/kickmix
T0=$(date +%s)
(cd "$WORK/dqi" && ~/venvs/tb313/bin/python build.py "$OUT/dqi.kmb")
echo "wall $(( $(date +%s) - T0 )) s; $(sha256sum "$OUT/dqi.kmb")"
