#!/bin/zsh
cd /Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/composite_qsk/qsk/engine_b/src.nosync
for s in "qb_fsp.py" "qb_bil.py" "qb_tq3.py" "qb_etak.py eta" "qb_etak.py K" "qb_numeps.py"; do
  n=$(echo $s | tr ' ' '_')
  t0=$(date +%s)
  /usr/bin/time -l python3 ${=s} > ../out.nosync/log_$n.txt 2> ../out.nosync/err_$n.txt
  rc=$?
  t1=$(date +%s)
  mem=$(grep "maximum resident" ../out.nosync/err_$n.txt | awk '{print int($1/1048576)}')
  echo "$(date '+%Y-%m-%d %H:%M:%S') analysis $s: exit $rc, $((t1-t0)) s, peak ${mem} MB" >> ../PROGRESS.txt
done
