#!/bin/zsh
cd /Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/composite_qsk/qsk/engine_b/src.nosync
while pgrep -f "python3 qb_" > /dev/null; do sleep 5; done
for s in "qb_etak.py eta:2" "qb_etak.py eta:3" "qb_etak.py eta:5" "qb_fsp.py" "qb_gmisc.py" "qb_numeps.py"; do
  n=$(echo $s | tr ' :' '__')
  t0=$(date +%s)
  /usr/bin/time -l python3 ${=s} > ../out.nosync/log_$n.txt 2> ../out.nosync/err_$n.txt
  rc=$?
  t1=$(date +%s)
  mem=$(grep "maximum resident" ../out.nosync/err_$n.txt | awk '{print int($1/1048576)}')
  echo "$(date '+%Y-%m-%d %H:%M:%S') analysis(v2) $s: exit $rc, $((t1-t0)) s, peak ${mem} MB" >> ../PROGRESS.txt
done
python3 qb_tq2.py 1500 > ../out.nosync/log_tq2d.txt 2>&1
echo "$(date '+%Y-%m-%d %H:%M:%S') TQ2 driver finished (three-output leg sets, sign v2)" >> ../PROGRESS.txt
