#!/bin/zsh
cd /Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/composite_qsk/composite/qa/work.nosync/scripts
for spec in "AV|AV:AV" "AV|AV:AV AV" "AV|AV:AW AW" "AV|AW:AV AW" "AW|AW:AW AW" "AV|AV:AsV c" "AV|AW:AsW c"; do
  pair=${spec%%:*}; ls=${spec#*:}
  echo "=== $pair / $ls  $(date +%H:%M:%S)"
  python3 -B b_at_a_points.py "$pair" "$ls" 99
done
echo "ALL DONE $(date +%H:%M:%S)"
