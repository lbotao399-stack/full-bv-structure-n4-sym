"""Render every math string of chapters 9-11 (reader blocks of reader_build/document.json) and of
manuscript/supplements/ch09-ch11.tex with the vendored KaTeX (throwOnError, strict warn)."""
import json, re, subprocess, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NODE = '/Users/libotao/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node'
doc = json.loads((ROOT / 'reader_build' / 'document.json').read_text())
bl = doc['blocks']
heads = [i for i, b in enumerate(bl) if b['type'] == 'heading' and b.get('level') == 1]
start = [i for i in heads if bl[i]['source'].startswith('9.')][0]
after = [i for i in heads if bl[i]['source'].startswith('12.')]
end = after[0] if after else len(bl)
items = []
for b in bl[start:end]:
    for key in ('source', 'caption'):
        s = b.get(key, '')
        if not s:
            continue
        if b['type'] == 'equation' and key == 'source':
            items.append(('display', s, 'blk'))
        else:
            for m in re.finditer(r'(?<!\\)\$([^$]*)\$', s):
                items.append(('inline', m.group(1), 'blk'))
    for row in b.get('rows', []) or []:
        for cell in row:
            for m in re.finditer(r'(?<!\\)\$([^$]*)\$', str(cell)):
                items.append(('inline', m.group(1), 'table'))


def tex_math(text, tag):
    text = '\n'.join(l for l in text.splitlines() if not l.lstrip().startswith('%'))
    out = []
    for m in re.finditer(r'\\\[(.*?)\\\]', text, flags=re.S):
        out.append(('display', m.group(1), tag))
    t2 = re.sub(r'\\\[(.*?)\\\]', ' ', text, flags=re.S)
    for m in re.finditer(r'(?<!\\)\$([^$]*)\$', t2):
        out.append(('inline', m.group(1), tag))
    return out


for c in ('ch09', 'ch10', 'ch11'):
    items += tex_math((ROOT / 'manuscript' / 'supplements' / f'{c}.tex').read_text(), 'supp-' + c)
out = HERE / 'katex_items.json'
out.write_text(json.dumps([[a, b] for a, b, c in items]))
(HERE / 'katex_tags.json').write_text(json.dumps([c for a, b, c in items]))
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
r = subprocess.run([NODE, '-e', js, str(katex), str(out), str(HERE / 'katex_tags.json')], capture_output=True, text=True, timeout=300)
print(r.stdout[-6000:]); print(r.stderr[-2000:], file=sys.stderr)
