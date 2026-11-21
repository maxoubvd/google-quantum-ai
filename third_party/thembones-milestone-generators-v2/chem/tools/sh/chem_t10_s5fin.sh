#!/bin/bash
# T10 seat 5, v2 = + W1 after the polish (LOCAL WSL, detached): THE END OF ONE CARD-30 TICKET OF tools/sh/chem_t10_s5tickets.txt (= seat 4's
# chem_t10_s4fin.sh + the wire-capable END repair):
#  1 wait for runs/<RUN>/DONE + rounded.npz; read the END and 'relations broken on C2's STORED words'
#  2 END over MAXEND (R85's gate, tools/r85_chem_budget_t1_gate.py <count>): stop (OVER marker); else the weights hold on
#    the run's own GPU (chem_rk7_whold.py, 220 steps)
#  3 broken > 0: the END repair (tools/chem_t10_wrep.py) on the held table (the card's CPU, cores HC)
#  4 the relation-keeping polish (tools/chem_t10_ppol.py --four --triples --dp 3) on the run's GPU
#  5 pull the polished source, rel_used_g.txt, log.txt to ~/qcc/ws/chem_t10/pfin/in/; READY_<RUN> for seat 4's queue
#   RUN=t10_y46a setsid nohup bash chem_t10_s5fin.sh > ~/qcc/ws/chem_t10/s5/fin_<RUN>.log 2>&1 &
: ${RUN:?}
T=${TICKETS:-<HOME>/Quantum/tools/sh/chem_t10_s5tickets.txt}   # race: TICKETS=<race file>
Q=<HOME>/Quantum
read -r K G NAME LAY PAR PMD5 CNT KL TC TS TCC <<< "$(grep -v '^#' $T | awk -v n="$RUN" '$3 == n')"
[ -z "$NAME" ] && { echo "no ticket $RUN"; exit 1; }
MAXEND=${MAXEND:-$(python3 $Q/tools/r85_chem_budget_t1_gate.py $CNT)}
HC=${HC:-$TC}
P=~/qcc/ws/chem_t10/pfin; S=~/qcc/ws/chem_t10/s4
M='s/[0-9]{1,3}(\.[0-9]{1,3}){3}/<addr>/g; s/port [0-9]+/port <n>/g'
mkdir -p $P/in $S/q; cd ~/qcc/ws/chem_wr7 || exit 1
echo "$(date -u +%H:%M:%SZ) waiting for the END of $RUN (card $K GPU $G, layout $LAY, parent $PAR, count $CNT, MAXEND $MAXEND)"
until K=$K timeout 120 ${CARDSH:-./card.sh} "test -e ~/qcc/ws/chem_wr7/runs/$RUN/DONE && test -s ~/qcc/ws/chem_wr7/runs/$RUN/rounded.npz && echo yes" < /dev/null 2>/dev/null | grep -q yes; do
  K=$K timeout 120 ${CARDSH:-./card.sh} "test -e ~/qcc/ws/chem_wr7/runs/$RUN/DONE && sleep 8 && test ! -s ~/qcc/ws/chem_wr7/runs/$RUN/rounded.npz && echo stopped" < /dev/null 2>/dev/null | grep -q stopped && { echo "$(date -u +%H:%M:%SZ) $RUN DONE without rounded.npz (stopped)"; echo "$RUN $CNT 999 STOPPED" > $S/q/OVER_$RUN; exit 3; }
  grep -q STOP ~/qcc/ws/chem_t10/s5/kill_$RUN.log 2>/dev/null && { echo "$(date -u +%H:%M:%SZ) $RUN stopped by its kill line"; echo "$RUN $CNT 999 KILLED" > $S/q/OVER_$RUN; exit 3; }
  sleep 60
done
sleep 5
K=$K timeout 120 ${CARDSH:-./card.sh} "grep -E 'all frozen|relp-e \(t10\): relations broken|relp-thaw' ~/qcc/ws/chem_wr7/runs/$RUN/log.txt | tail -n 3" < /dev/null 2>&1 | sed -E "$M" | cut -c1-250
EV=$(K=$K timeout 120 ${CARDSH:-./card.sh} "grep 'all frozen' ~/qcc/ws/chem_wr7/runs/$RUN/log.txt | tail -n 1" < /dev/null 2>/dev/null | grep -o -E "error_2 [0-9.]+" | awk '{print $2 * $2}')
NB=$(K=$K timeout 120 ${CARDSH:-./card.sh} "grep 'relations broken on C2' ~/qcc/ws/chem_wr7/runs/$RUN/log.txt | tail -n 1" < /dev/null 2>/dev/null | grep -o -E "STORED words [0-9]+" | grep -o -E "[0-9]+$")
echo "$(date -u +%H:%M:%SZ) END of $RUN: ${EV:-unknown} mHa^2 (MAXEND $MAXEND); relations broken on the stored words ${NB:-unknown}"
if ! python3 -c "import sys; sys.exit(0 if '${EV}' and float('${EV:-999}') <= float('$MAXEND') else 1)"; then
  echo "over MAXEND: no table step"; echo "$RUN $CNT ${EV:-999} OVER" > $S/q/OVER_$RUN; exit 3; fi
GPUENV="unset CUDA_VISIBLE_DEVICES JAX_PLATFORMS; export XLA_PYTHON_CLIENT_PREALLOCATE=false OMP_NUM_THREADS=2 QCC_GPU=$G"
PY='~/qcc/ws/chem_own/venv_gpu/bin/python'      # literal: the remote shell expands it (a local ~ would name WSL's home)
K=$K timeout 3000 ${CARDSH:-./card.sh} "cd ~/qcc/ws/chem_wr7 && mkdir -p pol/$RUN && ($GPUENV; taskset -c $HC $PY tools/chem_rk7_whold.py runs/$RUN/rounded.npz pol/$RUN/held.npz --steps ${HSTEPS:-220} --chunk 20 --seconds ${HSEC:-1800} > logs/whold_$RUN.log 2>&1); tail -n 2 logs/whold_$RUN.log | cut -c1-200; test -s pol/$RUN/held.npz && echo HELD" < /dev/null 2>&1 | sed -E "$M" | tee $P/in/${RUN}_hold.txt
grep -q HELD $P/in/${RUN}_hold.txt || { echo "HOLD FAILED"; exit 4; }
POLIN=pol/$RUN/held.npz
if [ "${NB:-1}" != 0 ]; then
  K=$K timeout 3000 ${CARDSH:-./card.sh} "cd ~/qcc/ws/chem_wr7 && (export JAX_PLATFORMS=cpu OMP_NUM_THREADS=8; taskset -c ${RC:-48-55} nice -n 5 $PY tools/chem_t10_wrep.py pol/$RUN/held.npz runs/$RUN/rel_used_g.txt pol/$RUN/wrep.npz > logs/wrep_$RUN.log 2>&1); grep -v Warning logs/wrep_$RUN.log | tail -n 2 | cut -c1-250" < /dev/null 2>&1 | sed -E "$M" | tee $P/in/${RUN}_wrep.txt
  grep -q "relations broken on the stored words 0" $P/in/${RUN}_wrep.txt || { echo "REPAIR FAILED"; exit 4; }
  POLIN=pol/$RUN/wrep.npz
fi
K=$K timeout 5400 ${CARDSH:-./card.sh} "cd ~/qcc/ws/chem_wr7 && ($GPUENV; taskset -c $HC $PY tools/chem_t10_ppol.py $POLIN runs/$RUN/rel_used_g.txt pol/$RUN/ppol.npz --minutes ${PMIN:-40} --four --triples --dp 3 --dp-minutes ${DPMIN:-15} > logs/ppol_$RUN.log 2>&1); tail -n 1 logs/ppol_$RUN.log | cut -c1-260; grep -q 'relations broken on the stored words 0' logs/ppol_$RUN.log && test -s pol/$RUN/ppol.npz && echo POLISHED" < /dev/null 2>&1 | sed -E "$M" | tee $P/in/${RUN}_ppol.txt
SRCN=held; SRCF=$POLIN
grep -q POLISHED $P/in/${RUN}_ppol.txt && { SRCN=ppol; SRCF=pol/$RUN/ppol.npz; }
# seat 5 v2: the wire-capable row polish (tools/chem_t10_wpol.py, W1) after ppol, on the run's GPU; taken only if it
# ends with 0 relations broken and writes its table
if [ $SRCN = ppol ] && [ "${WPOL:-1}" = 1 ]; then
  K=$K timeout 3600 ${CARDSH:-./card.sh} "cd ~/qcc/ws/chem_wr7 && ($GPUENV; taskset -c $HC $PY tools/chem_t10_wpol.py pol/$RUN/ppol.npz runs/$RUN/rel_used_g.txt pol/$RUN/wpol.npz --passes ${WPASS:-4} --minutes ${WMIN:-30} > logs/wpol_$RUN.log 2>&1); grep END logs/wpol_$RUN.log | cut -c1-300; grep -q 'relations broken on the stored words 0' logs/wpol_$RUN.log && test -s pol/$RUN/wpol.npz && echo WPOLISHED" < /dev/null 2>&1 | sed -E "$M" | tee $P/in/${RUN}_wpol.txt
  grep -q WPOLISHED $P/in/${RUN}_wpol.txt && { SRCN=wpol; SRCF=pol/$RUN/wpol.npz; }
fi
[ $SRCN = held ] && [ "$POLIN" != pol/$RUN/held.npz ] && SRCN=wrep
K=$K ${CARDSH:-./card.sh} --pull $SRCF $P/in/${RUN}_src.npz 2>&1 | sed -E "$M"
K=$K ${CARDSH:-./card.sh} --pull runs/$RUN/rel_used_g.txt $P/in/${RUN}_rel_used_g.txt 2>&1 | sed -E "$M"
K=$K ${CARDSH:-./card.sh} --pull runs/$RUN/log.txt $P/in/${RUN}_log.txt 2>&1 | sed -E "$M"
cp ~/qcc/ws/chem_t10/k2/rel/t10s4/$LAY.txt $P/in/${RUN}_layout.txt
test -s $P/in/${RUN}_src.npz && test -s $P/in/${RUN}_rel_used_g.txt || { echo "PULL FAILED"; exit 4; }
EPF=$P/in/${RUN}_ppol.txt; [ $SRCN = wpol ] && EPF=$P/in/${RUN}_wpol.txt
EP=$(grep -o -E "\-> [0-9.]+ mHa\^2 \(" $EPF | tail -1 | grep -o -E "[0-9]+\.[0-9]+")
[ "${READY:-1}" = 1 ] && echo "$RUN $CNT ${EP:-$EV} $SRCN" > $S/q/READY_$RUN
echo "$(date -u +%H:%M:%SZ) $RUN READY for the table step: count $CNT, END $EV, source $SRCN ($(md5sum $P/in/${RUN}_src.npz | cut -c1-8)) at ${EP:-?} mHa^2"
