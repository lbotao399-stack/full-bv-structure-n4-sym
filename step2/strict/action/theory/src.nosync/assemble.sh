#!/bin/sh
# Assemble action_theory.tex from the parts and compile into _build (pdflatex, three passes).
cd "$(dirname "$0")/.." || exit 1
P=src.nosync/parts
cat $P/part1.tex $P/part2.tex $P/part2b.tex $P/part3.tex $P/part4.tex $P/part5.tex \
 | sed -E 's/\\eqref\{eq:(c[0-9]+-[A-Za-z0-9]+)\}/(\1)/g' > action_theory.tex
mkdir -p _build
for i in 1 2 3; do pdflatex -interaction=nonstopmode -halt-on-error -output-directory=_build action_theory.tex > _build/run$i.log 2>&1 || exit 1; done
grep -c 'undefined' _build/action_theory.log
