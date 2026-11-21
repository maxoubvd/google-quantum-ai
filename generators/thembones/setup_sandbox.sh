#!/usr/bin/env bash
# New file added by the team (Maxime Bendavid), 2026-11-21, plan step S2.26.
# Lays out team ThemBones' milestone generators (Apache-2.0, third_party/thembones-milestone-generators-v2/, pristine)
# in a sandbox "home" in WSL, the way their README's "Where the folders go" table asks, and creates the Python 3.13
# venv with the recorded versions. Run inside WSL:  bash generators/thembones/setup_sandbox.sh
# Writes only under $TBH (default ~/tb) and ~/venvs/tb313. Time: about 1 min (plus the venv install the first time).
set -eu
R=$(cd "$(dirname "$0")/../.." && pwd)
DS=$R/third_party/thembones-milestone-generators-v2
TBH=${TBH:-$HOME/tb}
rm -rf $TBH
mkdir -p $TBH/Quantum/tools $TBH/qcc/ws/invsqrt_fc/s18/cfg $TBH/qcc/ws/invsqrt_fc/s18/files/M $TBH/qcc/ws/invsqrt_fc/best \
  $TBH/qcc/ws/invsqrt_fc/s19/q18 $TBH/qcc/ws/dqi_zn/src $TBH/qcc/ws/dqi_int/zptrace $TBH/qcc/data \
  $TBH/qcc/ws/s5/inv_h1end/best/b196lend
cp -r $DS $TBH/ds
# the only edit their README asks for: the redacted home path
grep -rl '<HOME>' $TBH/ds | xargs sed -i "s#<HOME>#$TBH#g"
for t in invsqrt dqi chem; do cp -r $TBH/ds/$t/tools/. $TBH/Quantum/tools/; done
cp -r $TBH/ds/invsqrt/src_s18_srcM $TBH/qcc/ws/invsqrt_fc/s18/srcM
cp $TBH/ds/invsqrt/cfg/cfg_M.json $TBH/ds/invsqrt/cfg/cfg_b312g.json $TBH/qcc/ws/invsqrt_fc/s18/cfg/
cp -r $TBH/ds/invsqrt/data/r59b0 $TBH/qcc/ws/invsqrt_fc/best/r59b0
cp $TBH/ds/invsqrt/prog/q18_p7.json $TBH/qcc/ws/invsqrt_fc/s19/q18/q18_p7.json
cp -r $TBH/ds/dqi/kit/kv_unfused5 $TBH/ds/dqi/kit/kv_unfused4_z2 $TBH/qcc/ws/dqi_zn/src/
cp $TBH/ds/dqi/kit/z9C_oracle.json $TBH/qcc/ws/dqi_int/zptrace/
cp -r $R/google-verifiable-quantum-circuit-compilation/kickmix $TBH/qcc/data/kickmix
# Python 3.13 venv: numpy 2.5.4 and scipy 1.18.1 as recorded; numba 0.68.1 (llvmlite 0.50.0) is what pip resolved on
# 2026-11-21 (their record does not name a numba version). The venv itself must sit at ~/qcc/.venv (a symlinked
# interpreter alone loses the venv's site-packages).
if [ ! -x ~/venvs/tb313/bin/python ]; then
  ~/venvs/gqai313/bin/python -m venv ~/venvs/tb313
  ~/venvs/tb313/bin/pip install -q numpy==2.5.4 scipy==1.18.1 numba==0.68.1
fi
ln -sfn ~/venvs/tb313 $TBH/qcc/.venv
# the base the invsqrt build is proved against: the set's own invsqrt.kmb (c639a045 = their b196lend), as text
B=$TBH/qcc/ws/s5/inv_h1end/best/b196lend
cp $DS/set_invsqrt/invsqrt.kmb $B/invsqrt.kmb
$TBH/qcc/.venv/bin/python -c "
import sys; sys.path.insert(0, '$TBH/qcc/data/kickmix'); import kickmix as km
km.Circuit.from_file('$B/invsqrt.kmb', format='kmb').to_file('$B/invsqrt.kmx', format='kmx')"
echo "sandbox ready at $TBH ($(du -sh $TBH | cut -f1))"
