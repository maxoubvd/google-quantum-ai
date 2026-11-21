#!/usr/bin/env bash
# dqi_int_mhr_go.sh -- B5 (21st holder, Nov 17; runs ON box 6): best_q29kcgABC's whole recipe (its SOURCE.txt reproduce
# command, unchanged) + tools/dqi_int_mhr.py installed LAST (the kit's restore with hi made by the conv), into FRESH
# folders through dqi_int_t54.sh. The module's call record goes to ~/qcc/ws/dqi_int/mhr/<fresh>.json.
#   bash dqi_int_mhr_go.sh <fresh> <c1,c2> <c1-c2> [DQI_INT_MHR_LEVELS=a,b] [VAR=V ...]
set -u
N=${1:?fresh}; CC=${2:?c1,c2}; PC=${3:?c1-c2}; shift 3
D=$HOME/qcc/ws/dqi_int; mkdir -p $D/mhr
EXT=dqi_int_tail901,dqi_zn_pdl2,dqi_zn_kpark,dqi_zn_alias,dqi_int_aliaslog,dqi_z2_pdr,r46_dqi_syn_tail,r46_dqi_ext_fgret,r48_dqi_ret2_lateflag,dqi_int_root4,dqi_int_merge49b,dqi_int_zpk2,dqi_int_kcg,dqi_int_mhr
B5_DWX=dqi_int_dwxz.sh B5_FILES_ADD=dqi_int_kcg.py,dqi_int_mhr.py${B5_MHR_FILES:+,$B5_MHR_FILES} \
  bash $HOME/Quantum/tools/sh/dqi_int_t54.sh $N $CC $PC 100 $D/src/jt_tensors_r54c.json \
  DQI_INT_EXT=$EXT${B5_MHR_EXT:+,$B5_MHR_EXT} DQI_INT_KCG_A=1 DQI_INT_KCG_C=1 DQI_INT_MHR_OUT=$D/mhr/$N.json "$@"
