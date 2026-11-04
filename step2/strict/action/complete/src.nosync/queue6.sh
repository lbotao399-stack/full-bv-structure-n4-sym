#!/bin/zsh
cd /Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/strict/action/complete/src.nosync
until [ -f ../out.nosync/queue5.done ]; do sleep 20; done
python3 -u kb1.py A > ../out.nosync/kb1_A.log 2>&1
echo QUEUE6 DONE >> ../out.nosync/queue6.done
