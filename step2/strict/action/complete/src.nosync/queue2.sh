#!/bin/zsh
cd /Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/strict/action/complete/src.nosync
python3 -u eta1.py B - > ../out.nosync/eta1_Bm.log 2>&1
python3 -u eta1.py A + > ../out.nosync/eta1_Ap.log 2>&1
python3 -u eta1.py A - > ../out.nosync/eta1_Am.log 2>&1
python3 -u xieta.py B + > ../out.nosync/xieta_Bp.log 2>&1
python3 -u xieta.py B - > ../out.nosync/xieta_Bm.log 2>&1
echo QUEUE2 DONE >> ../out.nosync/queue2.done
