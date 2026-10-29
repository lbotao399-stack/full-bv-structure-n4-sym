"""Apply every finding to scratch copies (step2/qa/ch11_review/patched.nosync/), check that the replacements are
disjoint (each replace still occurs once when applied in sequence), then convert the patched chapter with the
manuscript chapters 0-10 and check refs; the manuscript files are not touched."""
import json, re, sys, shutil
from pathlib import Path
HERE = Path(__file__).resolve().parent
REV = HERE.parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / 'reader_build'))
import convert
F = json.loads((REV / 'findings.json').read_text())
src = {'chapter': (ROOT / 'manuscript' / 'ch11_composite_and_odd_kernels.tex').read_text(),
       'supplement': (ROOT / 'manuscript' / 'supplements' / 'ch11.tex').read_text()}
txt = dict(src)
for f in F:
    n = txt[f['file']].count(f['replace'])
    assert n == 1, (f['location'], n)
    txt[f['file']] = txt[f['file']].replace(f['replace'], f['by'])
out = REV / 'patched.nosync'; out.mkdir(exist_ok=True)
(out / 'ch11_composite_and_odd_kernels.tex').write_text(txt['chapter'])
(out / 'ch11_supplement.tex').write_text(txt['supplement'])
files = sorted(convert.MS.glob('ch[0-9][0-9]_*.tex'))
files = [f for f in files if not f.name.startswith('ch11')] + [out / 'ch11_composite_and_odd_kernels.tex']
convert.OUT = REV
c = convert.Converter(); c.figdir = REV / 'figures.nosync'
blocks = c.run(files)
unresolved = sorted(set(re.findall(r'\?\?[\w:.-]+', json.dumps(blocks))))
labels = set(c.labels)
sup_refs = set(re.findall(r'\\(?:eq)?ref\{([^}]*)\}', txt['supplement']))
print(json.dumps({'findings': len(F), 'blocks': len(blocks), 'unresolved': unresolved,
                  'supplement_missing_refs': sorted(r for r in sup_refs if r not in labels)}, indent=1))
