#!/bin/sh
# Replace chapters 9-11 and append chapter 12 in a reader at port $1, then re-place all supplements.
set -e
cd "$(dirname "$0")/.."
python3 reader_build/merge_into_live.py $1 reader_build/ch9_blocks.json --replace-prefix c9 | tail -2
python3 reader_build/merge_into_live.py $1 reader_build/ch10_blocks.json --replace-prefix c10 | tail -2
python3 reader_build/merge_into_live.py $1 reader_build/ch11_blocks.json --replace-prefix c11 | tail -2
python3 reader_build/merge_into_live.py $1 reader_build/ch12_blocks.json | tail -2
python3 reader_build/insert_supplements.py $1 | tail -3
