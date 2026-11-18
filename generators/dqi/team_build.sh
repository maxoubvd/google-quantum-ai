#!/usr/bin/env bash
# New file added by the team (Maxime Bendavid), 2026-11-18, derived from build.sh (Copyright 2026 Syed Asad Ali,
# Apache-2.0). Builds the dqi decoder through the team layer `src/team/build_team.py` with a chosen extractor
# cutoff (idea DQI-A03).
#
#   bash team_build.sh J_MIN ENGINE_DIR VERIFIERS_DIR OUT_DIR [EXPECTED_SHA256]
#
# Same host-file handling and pins as build.sh (kickmix.so and the official dqi verifier are copied in and their
# SHA-256 checked). With J_MIN = 15 the output must be byte-identical to build.sh's (7cbe8ecf...): that run is the
# neutrality check of the team layer. Environment: PY (CPython 3.13, not -O, with numpy; default python3),
# WORK (scratch folder on a fast disk; default a new temp folder).
set -euo pipefail
die() { echo "team_build.sh: ERROR: $*" >&2; exit 1; }
[ "$#" -ge 4 ] || die "usage: bash team_build.sh J_MIN ENGINE_DIR VERIFIERS_DIR OUT_DIR [EXPECTED_SHA256]"
J_MIN=$1; ENGINE_DIR=$2; VERIFIERS_DIR=$3; OUT_DIR=$4; EXPECTED=${5:-}
PY=${PY:-python3}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ENGINE_SHA=c66b3216c4029e0cdfd09054306aa0feadf53ab2f50c1164e37f752874e5005e
VERIFIER_SHA=f14274bc0943fb9769ddb01131be6063191a15009fe1caf1d4bdd9fd25194e65
unset PYTHONOPTIMIZE || true
export PYTHONDONTWRITEBYTECODE=1
"$PY" -c 'import sys; sys.exit(0 if (sys.implementation.name == "cpython" and sys.version_info[:2] == (3, 13) and __debug__) else 1)' \
  || die "PY=$PY must be CPython 3.13 without optimization"
sha256() { "$PY" -c 'import hashlib,sys
h=hashlib.sha256()
with open(sys.argv[1],"rb") as f:
    for b in iter(lambda: f.read(1<<20), b""): h.update(b)
print(h.hexdigest())' "$1"; }
WORK=${WORK:-$(mktemp -d "${TMPDIR:-/tmp}/dqi_team.XXXXXX")}
mkdir -p "$WORK"
RUN=$(mktemp -d "$WORK/dqi_team_j${J_MIN}.XXXXXX")
ROOT=$RUN/root
mkdir -p "$ROOT"
cp -R "$HERE/src/." "$ROOT/"
chmod -R u+w "$ROOT"
OFFICIAL=$ROOT/pro1_controller_bundle/1_task/official
mkdir -p "$OFFICIAL" "$ROOT/pro1_controller_bundle/1_task/verifiers/dqi"
cp "$ENGINE_DIR/kickmix.so" "$OFFICIAL/kickmix.so"
cp "$VERIFIERS_DIR/dqi/verify_submission_dqi.py" "$ROOT/pro1_controller_bundle/1_task/verifiers/dqi/verify_submission_dqi.py"
[ "$(sha256 "$OFFICIAL/kickmix.so")" = "$ENGINE_SHA" ] || die "kickmix.so is not the official engine"
[ "$(sha256 "$ROOT/pro1_controller_bundle/1_task/verifiers/dqi/verify_submission_dqi.py")" = "$VERIFIER_SHA" ] || die "not the official verifier"
OUT_KMB=$RUN/out/dqi.kmb
START=$SECONDS
(cd "$ROOT" && "$PY" team/build_team.py --j-min "$J_MIN" --output "$OUT_KMB" --log-directory "$RUN/logs") > "$RUN/build.stdout" 2>&1 \
  || { tail -n 30 "$RUN/build.stdout" >&2; die "build failed after $((SECONDS - START)) s; log $RUN/build.stdout"; }
GOT=$(sha256 "$OUT_KMB")
mkdir -p "$OUT_DIR"
mv -f "$OUT_KMB" "$OUT_DIR/dqi.kmb"
tail -n 1 "$RUN/build.stdout"
echo "BUILT dqi.kmb $GOT (j_min $J_MIN, $((SECONDS - START)) s)"
rm -rf "$ROOT"
if [ -n "$EXPECTED" ] && [ "$GOT" != "$EXPECTED" ]; then die "sha256 $GOT differs from expected $EXPECTED"; fi
