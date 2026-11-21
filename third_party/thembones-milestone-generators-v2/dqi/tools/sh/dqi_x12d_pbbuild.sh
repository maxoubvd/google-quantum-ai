#!/usr/bin/env bash
# Workstream dqi_x12d: build builder J's 24-bit decoder (ip4c907's recipe: J's entry tools/dqi_junc_r8_entry.py, the build
# line's patches and J's patches, read-only) on a given source tree, into ~/qcc/ws/dqi_x12d/builds/<name> (nothing is
# written in J's tree), then print the build line and the stage peaks.
#   bash tools/sh/dqi_x12d_pbbuild.sh <name> <source_dir> [core=9] [extra DQI_JUNC_R8P patches, comma list]
set -u
W=~/qcc/ws/dqi_x12d; T=<HOME>/Quantum/tools; PY=~/qcc/.venv/bin/python
NAME=${1:?name}; SRC=${2:?source dir}; CORE=${3:-9}; EXTRA=${4:-}
export DQI_BUILD_SOURCE=$SRC DQI_BUILD_PATCHES=prefix5,t1,r4jt,r4k,r5w,r6n,peak DQI_BUILD_FIRST=3
export DQI_BUILD_R4JT=$HOME/qcc/ws/dqi_build/r8/jt_tensors_lens.json
export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem${EXTRA:+,$EXTRA}
export DQI_JUNC_FWD_LEVELS= DQI_JUNC_REV_KINDS= DQI_JUNC_R8_CEIL=907
mkdir -p $W/builds/$NAME $W/logs
cd $W/builds/$NAME
nice -n 10 taskset -c "$CORE" $PY -B $T/dqi_convert_build.py $T/dqi_junc_r8_entry.py $W/builds/$NAME/dqi.kmb > $W/logs/build_$NAME.log 2>&1
tail -n 1 $W/logs/build_$NAME.log | cut -c1-200
[ -f $W/builds/$NAME/dqi.kmb ] || { grep -a -B3 -A14 'Traceback' $W/logs/build_$NAME.log | tail -n 40 | cut -c1-260; exit 1; }
echo "$NAME source=$SRC recipe=ip4c907 (prefix5,t1,r4jt,r4k,r5w,r6n,peak + ipjunc,nofullcache,cachelevels,latem${EXTRA:+,$EXTRA}; ceiling 907)" > $W/builds/$NAME/SOURCE.txt
nice -n 10 taskset -c "$CORE" $PY -B $T/dqi_build_stagepeak.py $W/builds/$NAME/dqi.kmb | tee $W/logs/peak_$NAME.txt
