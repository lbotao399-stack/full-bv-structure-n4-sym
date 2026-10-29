"""Append findings (JSON list or object on stdin) to step2/qa/ch11_review/findings.json after checking that each
`replace` occurs exactly once in the file named by `file` ('chapter' or 'supplement'); replaces an entry with the
same location and replace string."""
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
REV = HERE.parent
ROOT = HERE.parents[3]
FILES = {'chapter': ROOT / 'manuscript' / 'ch11_composite_and_odd_kernels.tex',
         'supplement': ROOT / 'manuscript' / 'supplements' / 'ch11.tex'}
inp = json.loads(sys.stdin.read())
if isinstance(inp, dict):
    inp = [inp]
path = REV / 'findings.json'
data = json.loads(path.read_text()) if path.exists() else []
for f in inp:
    assert f['severity'] in ('fatal', 'major', 'minor', 'markup'), f['severity']
    for k in ('location', 'severity', 'file', 'claim', 'evidence', 'replace', 'by'):
        assert k in f and f[k], k
    target = FILES[f['file']]
    text = target.read_text(encoding='utf-8')
    n = text.count(f['replace'])
    if n != 1:
        sys.exit(f"replace occurs {n} times in {target.name}: {f['replace'][:100]!r}")
    assert f['replace'] != f['by']
    f['file_path'] = str(target.relative_to(ROOT))
    data = [d for d in data if not (d['location'] == f['location'] and d['replace'] == f['replace'])]
    data.append(f)
path.write_text(json.dumps(data, indent=1, ensure_ascii=False))
print('ok', len(data))
