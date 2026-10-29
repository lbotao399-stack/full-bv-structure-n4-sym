"""Run reader_build/convert.py's converter on all chapters with output redirected to this folder."""
import json, re, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / 'reader_build'))
import ch00_guard
mode = ch00_guard.install()
import convert
convert.OUT = HERE
files = sorted(convert.MS.glob('ch[0-9][0-9]_*.tex'))
c = convert.Converter(); c.figdir = HERE / 'figures.nosync'
blocks = c.run(files)
doc = {'id': 'check', 'title': 'check', 'blocks': blocks}
(HERE / 'document_check.json').write_text(json.dumps(doc, ensure_ascii=False, indent=1))
unresolved = sorted(set(re.findall(r'\?\?[\w:.-]+', json.dumps(doc))))
ch11 = [b for b in blocks if b['type'] == 'heading' and b.get('level') == 1]
print(json.dumps({'ch00': mode, 'blocks': len(blocks), 'labels': len(c.labels), 'unresolved': unresolved[:40],
                  'chapters': [b['source'][:60] for b in ch11]}, indent=1))
