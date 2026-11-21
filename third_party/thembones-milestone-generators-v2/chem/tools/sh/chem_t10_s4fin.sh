#!/bin/bash
# T10 seat 4 (LOCAL WSL, detached): THE END OF ONE RELATION TICKET (steps 1-3 of seat 3's chem_t10_pboxfin.sh), then a
# READY marker for the table-step queue (tools/sh/chem_t10_s4queue.sh, highest count first, two box slots).
#  1 wait for runs/<RUN>/DONE on its card; read the END (error_2^2) and the line 'relations broken on C2's STORED words'
#  2 END over MAXEND (401.2): stop here (written to the marker as OVER); else the weights hold on the card (220 steps)
#  3 pull held.npz (else rounded.npz), rel_used_g.txt, log.txt to ~/qcc/ws/chem_t10/pfin/in/; READY_<RUN> with count, END
#   RUN=t10_q2aa setsid nohup bash chem_t10_s4fin.sh > ~/qcc/ws/chem_t10/s4/fin_<RUN>.log 2>&1 &
: ${RUN:?}
T=<HOME>/Quantum/tools/sh/chem_t10_s4tickets.txt
read -r K G NAME LAY SEED CNT KILL FAC TC TS TCC <<< "$(grep -v '^#' $T | awk -v n="$RUN" '$3 == n')"
[ -z "$NAME" ] && { echo "no ticket $RUN"; exit 1; }
if [ "$K" = 19 ]; then B=$((32 + 8 * G)); else B=$((8 * G)); fi
MAXEND=${MAXEND:-$(python3 -c "print(round(399.34 + (3511 - int('${CNT#E}')) / 6.3, 3))")};   # the lead's gate (17:14Z)
HC=$B,$((B + 64)); [ -n "$TC" ] && HC=$TC
P=~/qcc/ws/chem_t10/pfin; S=~/qcc/ws/chem_t10/s4
M='s/[0-9]{1,3}(\.[0-9]{1,3}){3}/<addr>/g; s/port [0-9]+/port <n>/g'
mkdir -p $P/in $S/q; cd ~/qcc/ws/chem_wr7 || exit 1
echo "$(date -u +%H:%M:%SZ) waiting for the END of $RUN (card $K GPU $G, layout $LAY, count $CNT)"
until K=$K timeout 120 ./card.sh "test -e ~/qcc/ws/chem_wr7/runs/$RUN/DONE && test -s ~/qcc/ws/chem_wr7/runs/$RUN/rounded.npz && echo yes" < /dev/null 2>/dev/null | grep -q yes; do
  K=$K timeout 120 ./card.sh "test -e ~/qcc/ws/chem_wr7/runs/$RUN/DONE && echo stopped" < /dev/null 2>/dev/null | grep -q stopped && { echo "$(date -u +%H:%M:%SZ) $RUN DONE without rounded.npz (stopped)"; echo "$RUN $CNT 999 STOPPED" > $S/q/OVER_$RUN; exit 3; }
  grep -q STOP $S/kill_$RUN.log 2>/dev/null && { echo "$(date -u +%H:%M:%SZ) $RUN stopped by its kill line: $(grep STOP $S/kill_$RUN.log)"; echo "$RUN $CNT 999 KILLED" > $S/q/OVER_$RUN; exit 3; }
  sleep 60
done
sleep 5
K=$K timeout 120 ./card.sh "grep -E 'all frozen|relp-e \(t10\): relations broken|fp: rows placed' ~/qcc/ws/chem_wr7/runs/$RUN/log.txt | tail -n 3" < /dev/null 2>&1 | sed -E "$M" | cut -c1-250
EV=$(K=$K timeout 120 ./card.sh "grep 'all frozen' ~/qcc/ws/chem_wr7/runs/$RUN/log.txt | tail -n 1" < /dev/null 2>/dev/null | grep -o -E "error_2 [0-9.]+" | awk '{print $2 * $2}')
echo "$(date -u +%H:%M:%SZ) END of $RUN: ${EV:-unknown} mHa^2 (MAXEND $MAXEND)"
case $CNT in E*) echo "measurement ticket (count $CNT not counted by C3's builder): no table step"; echo "$RUN $CNT ${EV:-999} MEASURE" > $S/q/OVER_$RUN; exit 0 ;; esac
if ! python3 -c "import sys; sys.exit(0 if '${EV}' and float('${EV:-999}') <= float('$MAXEND') else 1)"; then
  echo "over MAXEND: no table step"; echo "$RUN $CNT ${EV:-999} OVER" > $S/q/OVER_$RUN; exit 3; fi
SRCN=rounded
K=$K timeout 3000 ./card.sh "cd ~/qcc/ws/chem_wr7 && mkdir -p pol/$RUN && test -e tools/chem_rk7_whold.py && (unset CUDA_VISIBLE_DEVICES JAX_PLATFORMS; export XLA_PYTHON_CLIENT_PREALLOCATE=false OMP_NUM_THREADS=2 QCC_GPU=$G; taskset -c $HC ~/qcc/ws/chem_own/venv_gpu/bin/python tools/chem_rk7_whold.py runs/$RUN/rounded.npz pol/$RUN/held.npz --steps ${HSTEPS:-220} --chunk 20 --seconds ${HSEC:-1800} > logs/whold_$RUN.log 2>&1); tail -n 2 logs/whold_$RUN.log | cut -c1-200; test -s pol/$RUN/held.npz && echo HELD" < /dev/null 2>&1 | sed -E "$M" | tee $P/in/${RUN}_hold.txt
grep -q HELD $P/in/${RUN}_hold.txt && SRCN=held
# seat 4: the relation-keeping polish after the hold (tools/chem_t10_ppol.py, the card's GPU, PMIN minutes); its output is
# taken only if it ends with 'relations broken on the stored words 0' (the tool asserts it) and a lower error
if [ $SRCN = held ] && [ "${POLISH:-1}" = 1 ]; then
  K=$K timeout 5400 ./card.sh "cd ~/qcc/ws/chem_wr7 && test -e tools/chem_t10_ppol.py && (unset CUDA_VISIBLE_DEVICES JAX_PLATFORMS; export XLA_PYTHON_CLIENT_PREALLOCATE=false OMP_NUM_THREADS=2 QCC_GPU=$G; taskset -c $HC ~/qcc/ws/chem_own/venv_gpu/bin/python tools/chem_t10_ppol.py pol/$RUN/held.npz runs/$RUN/rel_used_g.txt pol/$RUN/ppol.npz --minutes ${PMIN:-40} --four --triples --dp 3 --dp-minutes 15 > logs/ppol_$RUN.log 2>&1); tail -n 1 logs/ppol_$RUN.log | cut -c1-260; grep -q 'relations broken on the stored words 0' logs/ppol_$RUN.log && test -s pol/$RUN/ppol.npz && echo POLISHED" < /dev/null 2>&1 | sed -E "$M" | tee $P/in/${RUN}_ppol.txt
  grep -q POLISHED $P/in/${RUN}_ppol.txt && SRCN=ppol
fi
if [ $SRCN = ppol ]; then K=$K ./card.sh --pull pol/$RUN/ppol.npz $P/in/${RUN}_src.npz 2>&1 | sed -E "$M"
elif [ $SRCN = held ]; then K=$K ./card.sh --pull pol/$RUN/held.npz $P/in/${RUN}_src.npz 2>&1 | sed -E "$M"
else K=$K ./card.sh --pull runs/$RUN/rounded.npz $P/in/${RUN}_src.npz 2>&1 | sed -E "$M"; fi
K=$K ./card.sh --pull runs/$RUN/rel_used_g.txt $P/in/${RUN}_rel_used_g.txt 2>&1 | sed -E "$M"
K=$K ./card.sh --pull runs/$RUN/log.txt $P/in/${RUN}_log.txt 2>&1 | sed -E "$M"
cp ~/qcc/ws/chem_t10/k2/rel/t10s4/$LAY.txt $P/in/${RUN}_layout.txt
HV=$(grep -o -E "error_2 [0-9.]+ mHa|-> [0-9.]+ mHa\^2|[0-9]{3}\.[0-9]+ mHa\^2" $P/in/${RUN}_hold.txt | tail -n 1)
test -s $P/in/${RUN}_src.npz && test -s $P/in/${RUN}_rel_used_g.txt || { echo "PULL FAILED"; exit 4; }
echo "$RUN $CNT $EV $SRCN" > $S/q/READY_$RUN
echo "$(date -u +%H:%M:%SZ) $RUN READY for the table step: count $CNT, END $EV, source $SRCN ($(md5sum $P/in/${RUN}_src.npz | cut -c1-8)); hold: $HV"
