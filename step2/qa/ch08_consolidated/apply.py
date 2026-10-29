#!/usr/bin/env python3
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from edits import E
root = Path(__file__).resolve().parents[3]
src = root / 'step2/qa/ch08_consolidated/ch08_before_consolidation.tex'
dst = root / 'manuscript/ch08_definitions.tex'
t = src.read_text()
log = []
for e in E:
    n = t.count(e['old'])
    if n != 1:
        print('FAIL', e['id'], e['label'], 'count', n); sys.exit(1)
    t = t.replace(e['old'], e['new'])
    log.append(e['id'])
if '--dry' in sys.argv:
    print('dry ok', len(log)); sys.exit(0)
dst.write_text(t)
print('applied', len(log), 'edits;', len(t), 'chars')
