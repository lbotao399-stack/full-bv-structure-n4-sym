#!/bin/zsh
cd /Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/strict/action/complete/src.nosync
python3 -u dxi.py B - > ../out.nosync/dxi_Bm.log 2>&1
python3 -u eta1.py B - > ../out.nosync/eta1_Bm.log 2>&1
python3 -u eta1.py A + > ../out.nosync/eta1_Ap.log 2>&1
python3 -u eta1.py A - > ../out.nosync/eta1_Am.log 2>&1
python3 -u dxi.py A + > ../out.nosync/dxi_Ap.log 2>&1
python3 -u dxi.py A - > ../out.nosync/dxi_Am.log 2>&1
echo QUEUE1 DONE >> ../out.nosync/queue1.done
