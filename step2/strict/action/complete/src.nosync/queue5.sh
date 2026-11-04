#!/bin/zsh
cd /Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/strict/action/complete/src.nosync
until [ -f ../out.nosync/queue4.done ]; do sleep 20; done
python3 -u xieta.py B + > ../out.nosync/xieta_Bp.log 2>&1
python3 -u kb1.py B > ../out.nosync/kb1_B.log 2>&1
python3 -u xieta.py A + > ../out.nosync/xieta_Ap.log 2>&1
python3 -u xieta.py B - > ../out.nosync/xieta_Bm2.log 2>&1
echo QUEUE5 DONE >> ../out.nosync/queue5.done
