#!/bin/bash
# T10 seat 5 (LOCAL WSL): start ONE rounding ticket of tools/sh/chem_t10_s5tickets.txt on card 30 FROM ITS PARENT, from
# stage 0 (no head: the parent is new), = seat 4's chem_t10_s4tick.sh (frelgcjl = frelgcj + the late LM window; R67's cond
# order; R53's carrier rule; step-capped lever 15 / 300) with the late window of the recipe (--lm2-from 0.75 --lm2-until
# 1.0 --lm2-steps 30 --lm2-seconds 600) and the parent as the init. The parent must already sit in card 30's
# ~/qcc/ws/chem_wr7/pts/<PARENT>.npz with the md5 of the ticket line (checked here).
#   usage: bash chem_t10_s5tick.sh <name>
#   ticket line: K G NAME LAYOUT PARENT PARENT_MD5_12 COUNT KILL_LINES CPUS SDYN CONDC
T=${TICKETS:-<HOME>/Quantum/tools/sh/chem_t10_s5tickets.txt}   # race: TICKETS=<race file>
read -r K G NAME LAY PAR PMD5 CNT KL TC TS TCC <<< "$(grep -v '^#' $T | awk -v n="$1" '$3 == n')"
[ -z "$NAME" ] && { echo "no ticket $1"; exit 1; }
cd ~/qcc/ws/chem_wr7 || exit 1
M='s/[0-9]{1,3}(\.[0-9]{1,3}){3}/<addr>/g; s/port [0-9]+/port <n>/g'
K=$K timeout 120 ${CARDSH:-./card.sh} "cd ~/qcc/ws/chem_wr7 && test ! -e runs/$NAME/log.txt && test -s tools/${TOOLF:-chem_t10_frelgcjl.py} && test -s rel/t10s4/$LAY.txt && test \"\$(md5sum < pts/$PAR.npz | cut -c1-12)\" = $PMD5 && mkdir -p runs/$NAME logs && \
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 XLA_PYTHON_CLIENT_PREALLOCATE=false QCC_GPU=$G T10_ASSIGN_EVERY=${AEVERY:-3} && unset CUDA_VISIBLE_DEVICES JAX_PLATFORMS && \
( setsid nohup bash -c 'taskset -c $TC ~/qcc/ws/chem_own/venv_gpu/bin/python tools/${TOOLF:-chem_t10_frelgcjl.py} pts/$PAR.npz runs/$NAME --bits 7 --codes 0:1:8-14 --carriers none --no-relations --no-pairs --schedule fine4 --iters 10000 --late-iters 4000 --late-from 0.25 --final 20000 --kappa 2e-3 --rank sdyn --sdyn-cpu $TS --edges tree/${EDGES:-order_g9k2m1s0_tree.npz} --fp --relp @rel/t10s4/$LAY.txt --relp-g --relp-g-rounds 1 --smart-protect --relp-e exact --lm-until 0.1 --lm-lam-start 1e-3 --lm-stop-lam 0.1 --lm-steps 15 --lm-seconds 300 --rank-cond 0.01 --cond-cpus $TCC --stop-at 2.0 --keep-states 0.01,0.04,0.1,0.25,0.4,0.6,0.75,0.8,0.9,0.95 --max-minutes 900 --lm2-from 0.75 --lm2-until 1.0 --lm2-steps 30 --lm2-seconds 600 ${XTRA:-} >> logs/$NAME.log 2>&1; touch runs/$NAME/DONE' > /dev/null 2>&1 < /dev/null & ) ; sleep 2; echo started $NAME GPU $G from $PAR" < /dev/null 2>&1 | sed -E "$M"
echo "$(date -u +%H:%M:%SZ) t10 s5 ticket $NAME: card $K GPU $G layout $LAY parent $PAR ($PMD5) cpus $TC sdyn $TS cond $TCC" | tee -a ~/qcc/ws/chem_t10/s4/starts.log
