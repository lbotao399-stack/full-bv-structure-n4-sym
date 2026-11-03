"""Run reader_build/convert.py's Converter on all chapters, chapter 11 taken from a given file (default the draft),
outputs in this folder only; then render every math string of chapter 11 and of a given supplement file with KaTeX."""
import json, re, subprocess, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'reader_build'))
import convert
ch11 = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / 'ch11_new.tex'
supp = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / 'manuscript' / 'supplements' / 'ch11.tex'
files = [f for f in sorted(convert.MS.glob('ch[0-9][0-9]_*.tex')) if not f.name.startswith('ch11_')]
files.append(ch11 if re.match(r'ch11_', ch11.name) else ch11)
c = convert.Converter(); c.figdir = HERE / 'figures.nosync'
blocks = c.run(files)
doc = {'id': 'check', 'title': 'check', 'blocks': blocks}
(HERE / 'document_check.json').write_text(json.dumps(doc, ensure_ascii=False, indent=1))
unresolved = sorted(set(re.findall(r'\?\?[\w:.-]+', json.dumps(doc))))
heads = [b['source'][:70] for b in blocks if b['type'] == 'heading' and b.get('level') in (1, 2) and b['source'].startswith('11')]
print(json.dumps({'blocks': len(blocks), 'labels': len(c.labels), 'unresolved': unresolved[:40], 'ch11_headings': heads}, indent=1))
# KaTeX
bl = blocks
start = [i for i, b in enumerate(bl) if b['type'] == 'heading' and b.get('level') == 1 and b['source'].startswith('11.')][0]
items = []
for b in bl[start:]:
    for key in ('source', 'caption'):
        s = b.get(key, '')
        if not s:
            continue
        if b['type'] == 'equation' and key == 'source':
            items.append(('display', s, 'ch11'))
        else:
            for m in re.finditer(r'(?<!\\)\$([^$]*)\$', s):
                items.append(('inline', m.group(1), 'ch11'))
def tex_math(text, tag):
    text = '\n'.join(l for l in text.splitlines() if not l.lstrip().startswith('%'))
    out = []
    for m in re.finditer(r'\\\[(.*?)\\\]', text, flags=re.S):
        out.append(('display', m.group(1), tag))
    for m in re.finditer(r'\\begin\{equation\}(.*?)\\end\{equation\}', text, flags=re.S):
        out.append(('display', re.sub(r'\\label\{[^}]*\}', '', m.group(1)), tag))
    t2 = re.sub(r'\\\[(.*?)\\\]', ' ', text, flags=re.S)
    t2 = re.sub(r'\\begin\{equation\}(.*?)\\end\{equation\}', ' ', t2, flags=re.S)
    for m in re.finditer(r'(?<!\\)\$([^$]*)\$', t2):
        out.append(('inline', m.group(1), tag))
    return out
items += tex_math(supp.read_text(), 'supplement')
out = HERE / 'math_items.json'
out.write_text(json.dumps([[a, b] for a, b, c in items]))
(HERE / 'math_tags.json').write_text(json.dumps([c for a, b, c in items]))
NODE = '/Users/libotao/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node'
katex = ROOT / 'reader' / 'notion-reader' / 'app' / 'vendor' / 'katex' / 'katex.js'
js = r"""
const katex = require(process.argv[1]);
const items = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));
const tags = JSON.parse(require('fs').readFileSync(process.argv[3], 'utf8'));
let bad = 0, warn = 0;
items.forEach(([kind, s], i) => {
  try { katex.renderToString(s, {displayMode: kind === 'display', throwOnError: true,
        strict: (code, msg) => { warn++; console.log('WARN', tags[i], i, code, msg.slice(0,120), '::', s.slice(0,100)); return 'ignore'; }}); }
  catch (e) { bad++; console.log('ERR', tags[i], i, kind, String(e.message).slice(0, 200), '::', s.slice(0, 160)); }
});
console.log(JSON.stringify({items: items.length, errors: bad, warnings: warn}));
"""
r = subprocess.run([NODE, '-e', js, str(katex), str(out), str(HERE / 'math_tags.json')], capture_output=True, text=True, timeout=300)
print(r.stdout[-6000:]); print(r.stderr[-2000:], file=sys.stderr)
