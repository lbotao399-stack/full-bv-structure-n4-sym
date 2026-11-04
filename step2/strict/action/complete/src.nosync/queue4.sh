#!/bin/zsh
cd /Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/strict/action/complete/src.nosync
python3 -u xieta.py B + > ../out.nosync/xieta_Bp.log 2>&1
python3 -u kb1.py B > ../out.nosync/kb1_B.log 2>&1
python3 -u xieta.py B - > ../out.nosync/xieta_Bm.log 2>&1
python3 -u xieta.py A + > ../out.nosync/xieta_Ap.log 2>&1
python3 -u kb1.py A > ../out.nosync/kb1_A.log 2>&1
python3 -u xieta.py A - > ../out.nosync/xieta_Am.log 2>&1
echo QUEUE4 DONE >> ../out.nosync/queue4.done
