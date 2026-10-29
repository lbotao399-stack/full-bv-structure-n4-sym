"""Run reader_build/convert.py's Converter on all chapters exactly as its __main__ does, with the output
(document JSON, figures) redirected to step2/qa/ch11_review so reader_build/document.json stays untouched."""
import json, re, sys, subprocess
from pathlib import Path
HERE = Path(__file__).resolve().parent
REV = HERE.parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / 'reader_build'))
import convert
src = (ROOT / 'reader_build' / 'convert.py').read_text()
main = src[src.index("if __name__"):]
print(main[:1200], file=sys.stderr)
convert.OUT = REV
files = sorted(convert.MS.glob('ch[0-9][0-9]_*.tex'))
c = convert.Converter(); c.figdir = REV / 'figures.nosync'
blocks = c.run(files)
doc = {'id': 'check', 'title': 'check', 'blocks': blocks}
(REV / 'document_check.json').write_text(json.dumps(doc, ensure_ascii=False, indent=1))
unresolved = sorted(set(re.findall(r'\?\?[\w:.-]+', json.dumps(doc))))
print(json.dumps({'files': [f.name for f in files], 'blocks': len(blocks), 'labels': len(c.labels), 'unresolved': unresolved[:40]}, indent=1))
