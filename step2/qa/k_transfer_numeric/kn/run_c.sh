#!/bin/zsh
cd /Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/k_transfer_numeric/kn
while [ ! -f log_ab_done.txt ]; do sleep 10; done
# locality (scaling) test of the defect along K_b: p1 -> 2 p1
for s in AW_AW c_As; do
  KN_P1=2,4,0,2 KN_EPS=1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5 KN_T=1e2,1e4 python3 kn_run_K.py defect $s p2x > log_defect_${s}_p2x.txt 2>&1
done
# quadrature convergence: per_decade 14 vs 10
KN_PD=14 KN_EPS=1e-2,1e-3,1e-4 KN_T=1e2,1e4,1e6 python3 kn_run_K.py defect AW_AW pd14 > log_defect_AW_AW_pd14.txt 2>&1
KN_PD=14 KN_EPS=1e-2,1e-3,1e-4 KN_T=1e2,1e4,1e6 python3 kn_run_K.py flow AW_AsW pd14 > log_flow_AW_AsW_pd14.txt 2>&1
export KN_EPS=1e-1,1e-2,1e-3
export KN_T=1e2,1e3,1e4
export KN_VALS=1e-2:1e2,1e-3:1e2
export KN_PD=8
export KN_FULL4=1
for s in AW_AsW c_cs l_ls AV_AsV; do
  KN_M=00,01,02,03,10,11,12,13 python3 kn_eta.py $s mA > log_eta_${s}_mA.txt 2>&1
  KN_M=20,21,22,23,30,31,32,33 python3 kn_eta.py $s mB > log_eta_${s}_mB.txt 2>&1
done
echo ALLDONE > log_c_done.txt
