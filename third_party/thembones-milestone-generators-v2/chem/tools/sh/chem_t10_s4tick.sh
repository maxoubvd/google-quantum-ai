#!/bin/bash
# T10 seat 4 (LOCAL WSL): start ONE relation ticket of the table tools/sh/chem_t10_s4tickets.txt on its card and GPU.
# = seat 3's chem_t10_ptick.sh (frelgc + R67's cond + R53's carrier rule, from G9s2's 4% head, step-capped lever) with
#   the jittered tool tools/chem_t10_frelgcj.py and the layout rel/t10s4/<layout>.txt; cores by the CARD18.md map
#   (card 19: base 32 + 8 g; cards 18 / 28 / 29: base 8 g): CPUS b,b+64 / SDYN b+1 / CONDC b+2..b+7,b+65.
#   usage: bash chem_t10_s4tick.sh <name>      (e.g. t10_q2ab)
T=<HOME>/Quantum/tools/sh/chem_t10_s4tickets.txt
read -r K G NAME LAY SEED CNT KILL FAC TC TS TCC <<< "$(grep -v '^#' $T | awk -v n="$1" '$3 == n')"
[ -z "$NAME" ] && { echo "no ticket $1"; exit 1; }
if [ "$K" = 19 ]; then B=$((32 + 8 * G)); else B=$((8 * G)); fi
CPUS=$B,$((B + 64)); SDYN=$((B + 1)); CONDC=$((B + 2))-$((B + 7)),$((B + 65))
[ -n "$TCC" ] && { CPUS=$TC; SDYN=$TS; CONDC=$TCC; }      # the table's own cores (card 30)
JIT=""; [ "$SEED" != "-" ] && JIT="--order-jitter ${EPS:-0.1} --order-seed $SEED"
cd ~/qcc/ws/chem_wr7 || exit 1
K=$K ./card.sh "cd ~/qcc/ws/chem_wr7 && test ! -e runs/$NAME/log.txt && test -s tools/${TOOLF:-chem_t10_frelgcj.py} && test -s rel/t10s4/$LAY.txt && mkdir -p runs/$NAME logs && cp ~/qcc/setup/t10_g9s2_state.npz runs/$NAME/state.npz && cp ~/qcc/setup/t10_g9s2_sdyn.npy runs/$NAME/sdyn.npy && \
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 XLA_PYTHON_CLIENT_PREALLOCATE=false QCC_GPU=$G T10_ASSIGN_EVERY=${AEVERY:-3} && unset CUDA_VISIBLE_DEVICES JAX_PLATFORMS && \
( setsid nohup bash -c 'taskset -c $CPUS ~/qcc/ws/chem_own/venv_gpu/bin/python tools/${TOOLF:-chem_t10_frelgcj.py} pts/c46c.npz runs/$NAME --bits 7 --codes 0:1:8-14 --carriers none --no-relations --no-pairs --schedule fine4 --iters 10000 --late-iters 4000 --late-from 0.25 --final 20000 --kappa 2e-3 --rank sdyn --sdyn-cpu $SDYN --edges tree/order_g9k2m1s0_tree.npz --fp --relp @rel/t10s4/$LAY.txt --relp-g --relp-g-rounds 1 --smart-protect --relp-e exact --lm-until 0.1 --lm-lam-start 1e-3 --lm-stop-lam 0.1 --lm-steps 15 --lm-seconds 300 --rank-cond 0.01 --cond-cpus $CONDC --stop-at 2.0 --keep-states 0.01,0.04,0.1,0.25,0.4,0.6,0.75,0.8 --max-minutes 900 $JIT ${XTRA:-} >> logs/$NAME.log 2>&1; touch runs/$NAME/DONE' > /dev/null 2>&1 < /dev/null & ) ; sleep 2; echo started $NAME GPU $G" 2>&1
echo "$(date -u +%H:%M:%SZ) t10 s4 ticket $NAME: card $K GPU $G layout $LAY seed $SEED cpus $CPUS sdyn $SDYN cond $CONDC" | tee -a ~/qcc/ws/chem_t10/s4/starts.log
