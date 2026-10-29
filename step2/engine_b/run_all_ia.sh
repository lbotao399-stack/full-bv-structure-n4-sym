#!/bin/zsh
cd /Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/engine_b
python3 run_ia.py cApAp 6 > ia_cApAp.log 2>&1
python3 run_ia.py clamlam 8 > ia_clamlam.log 2>&1
python3 run_ia.py cAVAV 11 > ia_cAVAV.log 2>&1
echo DONE > ia_all.done
