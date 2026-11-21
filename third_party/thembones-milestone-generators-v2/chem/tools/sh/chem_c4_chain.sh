#!/bin/bash
# chem_c4_chain.sh: builder C4 (holder 4). ON A BOX: rebuild C4's wired loaders in order on a table, each with its lever,
# passing the dead-toggled table of every winning build on to the next loader (a batch read by two loaders).
#   CL="c1 c2 ..." bash chem_c4_chain.sh TABLE WORKDIR "k:mode k:mode ..."
# modes: plain24 (24 outputs, C4_DROP), row24 (24 outputs + 4 dead-paid left-kernel vectors in group 3),
#        dead23 (23 outputs, 4 dead-paid relations + one cubic block in group 4, annealed group 3)
# Result: WORKDIR/L<k>/best/{angle_k_<k>.kmx,table_dead.npz}; the last table in WORKDIR/table_final.npz; WORKDIR/chain.log
T=$1; WD=$2; LIST=$3
cd ~/Quantum/tools || exit 1
mkdir -p "$WD"
export CHEM_HB=77 PY=$HOME/qcc/.venv_work/bin/python C4_DROP=1 MAXP=${MAXP:-48} SAMPLES=${SAMPLES:-150}
cur=$T
for km in $LIST; do
  k=${km%%:*}; mode=${km##*:}
  rm -rf "$WD/L$k"
  case $mode in
    plain24) env C4_T1=149 bash sh/chem_c4_loader2.sh $k 0 0 "$cur" "$WD" ;;
    row24)   env C4_T1=148 C4_DEAD_ROW=4:2 C4_DEAD_OUT=EMIT bash sh/chem_c4_loader2.sh $k 0 0 "$cur" "$WD" ;;
    row24b)  env C4_T1=148 C4_DEAD_ROW=3:2:72 C4_DEAD_ROW_TRIALS=100 C4_DEAD_OUT=EMIT \
               bash sh/chem_c4_loader2.sh $k 0 0 "$cur" "$WD" ;;   # holder 4, 18:38Z: 3 vectors, rank 72 (42 equations)
    dead23)  env C4_T1=144 C4_T2=143 C4_DEAD=4 C4_CUBIC=3 C4_CUBIC_TRIES=1000 C4_DEAD_OUT=EMIT \
               bash sh/chem_c4_loader2.sh $k 0 0 "$cur" "$WD" ;;
  esac
  if [ ! -f "$WD/L$k/best/angle_k_$k.kmx" ]; then         # the lever missed: rebuild at the plain target (old circuit is invalid)
    mv "$WD/L$k" "$WD/L${k}_miss"
    case $mode in
      row24|row24b)  env C4_T1=149 bash sh/chem_c4_loader2.sh $k 0 0 "$cur" "$WD" ;;
      dead23) env C4_T1=146 C4_T2=145 C4_DEAD=4 C4_DEAD_OUT=EMIT bash sh/chem_c4_loader2.sh $k 0 0 "$cur" "$WD" ;;
      plain24) env C4_T1=151 bash sh/chem_c4_loader2.sh $k 0 0 "$cur" "$WD" ;;
    esac
  fi
  echo "$(date -u +%H:%M:%S) k$k $mode: $(grep -h -E 'SUCCESS|end without|kept' $WD/L$k/run.log | tail -n 1 | cut -c1-120)" >> "$WD/chain.log"
  if [ -f "$WD/L$k/best/table_dead.npz" ]; then cur="$WD/L$k/best/table_dead.npz"; fi
done
cp "$cur" "$WD/table_final.npz"
echo "$(date -u +%H:%M:%S) done; final table from $cur" >> "$WD/chain.log"
