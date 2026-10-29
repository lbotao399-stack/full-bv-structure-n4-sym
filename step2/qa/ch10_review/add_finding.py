"""Append one finding (JSON on stdin) to findings.json after checking that `replace` occurs exactly once
in the file named by `file` (default: the chapter)."""
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
FILES = {'chapter': ROOT / 'manuscript' / 'ch10_two_letter_kernels.tex',
         'supplement': ROOT / 'manuscript' / 'supplements' / 'ch10.tex'}
f = json.loads(sys.stdin.read())
target = FILES[f.get('file', 'chapter')]
text = target.read_text(encoding='utf-8')
n = text.count(f['replace'])
if n != 1:
    sys.exit(f"replace occurs {n} times in {target.name}: {f['replace'][:80]!r}")
path = HERE / 'findings.json'
data = json.loads(path.read_text())
data = [d for d in data if d['location'] != f['location'] or d['replace'] != f['replace']]
data.append(f)
path.write_text(json.dumps(data, indent=1, ensure_ascii=False))
print('ok', len(data))
