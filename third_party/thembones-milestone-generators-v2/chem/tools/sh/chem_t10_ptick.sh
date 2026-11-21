#!/bin/bash
# T10 seat 3 (LOCAL WSL): start a P ticket on card K, GPU G: R65's frelg + R67's cond order + R53's carrier rule
# (tools/chem_t10_frelgc.py --relp-e exact), from a staged head, step-capped lever. ONE run a GPU.
#   K=18 G=3 NAME=t10_pA1 REL=rel/G10P53m.txt HEADDIR=t10/g10g04 EDGES=tree/order_g10k3s0_tree.npz INIT=pts/c46c.npz \
#     CPUS=24,88 SDYN=25 CONDC=26-31,89 [EMODE=exact] bash chem_t10_ptick.sh
#   HEADDIR: a card folder with state.npz + sdyn.npy; or HEADSETUP=g9s2 (~/qcc/setup/t10_g9s2_state.npz / _sdyn.npy).
cd ~/qcc/ws/chem_wr7 || exit 1
: ${K:?} ${G:?} ${NAME:?} ${REL:?} ${EDGES:?} ${INIT:?} ${CPUS:?} ${SDYN:?} ${CONDC:?}
if [ -n "$HEADSETUP" ]; then HC="cp ~/qcc/setup/t10_${HEADSETUP}_state.npz runs/$NAME/state.npz && cp ~/qcc/setup/t10_${HEADSETUP}_sdyn.npy runs/$NAME/sdyn.npy"
else HC="cp ${HEADDIR:?}/state.npz ${HEADDIR}/sdyn.npy runs/$NAME/"; fi
K=$K ./card.sh "cd ~/qcc/ws/chem_wr7 && test ! -e runs/$NAME/log.txt && mkdir -p runs/$NAME logs && $HC && \
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 XLA_PYTHON_CLIENT_PREALLOCATE=false QCC_GPU=$G T10_ASSIGN_EVERY=${AEVERY:-1} && unset CUDA_VISIBLE_DEVICES JAX_PLATFORMS && \
setsid nohup bash -c 'taskset -c $CPUS ~/qcc/ws/chem_own/venv_gpu/bin/python tools/chem_t10_frelgc.py $INIT runs/$NAME --bits 7 --codes 0:1:8-14 --carriers none --no-relations --no-pairs --schedule fine4 --iters 10000 --late-iters 4000 --late-from 0.25 --final 20000 --kappa 2e-3 --rank sdyn --sdyn-cpu $SDYN --edges $EDGES --fp --relp @$REL --relp-g --relp-g-rounds 1 --smart-protect --relp-e ${EMODE:-exact} --lm-until 0.1 --lm-lam-start 1e-3 --lm-stop-lam 0.1 --lm-steps 15 --lm-seconds 300 --rank-cond 0.01 --cond-cpus $CONDC --stop-at 2.0 --keep-states 0.01,0.04,0.1,0.25,0.4 --max-minutes 900 ${XTRA:-} >> logs/$NAME.log 2>&1; touch runs/$NAME/DONE' > /dev/null 2>&1 < /dev/null & sleep 2; echo started $NAME GPU $G" 2>&1
echo "$(date -u +%H:%M:%SZ) t10 P ticket $NAME: card $K GPU $G rel $REL head ${HEADDIR:-setup:$HEADSETUP} edges $EDGES cpus $CPUS sdyn $SDYN cond $CONDC"
