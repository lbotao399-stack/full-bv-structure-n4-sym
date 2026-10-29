#!/bin/zsh
cd /Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/k_transfer_numeric/kn
export KN_EPS=1e-1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4
export KN_T=1e1,1e2,1e3,1e4,1e5,1e6
for s in AW_AW AV_AV c_As l_l; do
  python3 kn_run_K.py defect $s main > log_defect_$s.txt 2>&1
done
for s in AW_AsW AV_AsV c_cs l_ls; do
  python3 kn_run_K.py flow $s main > log_flow_$s.txt 2>&1
done
echo ALLDONE > log_ab_done.txt
