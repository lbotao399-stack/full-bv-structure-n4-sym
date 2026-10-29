"""Render every math string of chapter 8 (reader blocks) with the vendored KaTeX (displayMode for equations, throwOnError true, strict 'ignore')."""
import json, re, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
NODE = '/Users/libotao/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node'
doc = json.loads((ROOT / 'reader_build' / 'document.json').read_text())
bl = doc['blocks']
h1 = [i for i, b in enumerate(bl) if b['type'] == 'heading' and b.get('level') == 1]
start = [i for i in h1 if bl[i]['source'].startswith('8.')][0]
end = min([i for i in h1 if i > start] + [len(bl)])
items = []
for b in bl[start:end]:
    s = b.get('source', '')
    if b['type'] == 'equation':
        items.append(('display', s))
    else:
        for m in re.finditer(r'(?<!\\)\$([^$]*)\$', s):
            items.append(('inline', m.group(1)))
out = Path(__file__).with_name('math_items_ch08.json')
out.write_text(json.dumps(items))
katex = ROOT / 'reader' / 'notion-reader' / 'app' / 'vendor' / 'katex' / 'katex.min.js'
js = r"""
const katex = require(process.argv[1]);
const items = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));
let bad = 0;
items.forEach(([kind, s], i) => {
  try { katex.renderToString(s, {displayMode: kind === 'display', throwOnError: true, strict: 'ignore'}); }
  catch (e) { bad++; console.log('ERR', i, kind, String(e.message).slice(0, 200), '::', s.slice(0, 160)); }
});
console.log(JSON.stringify({items: items.length, errors: bad}));
"""
r = subprocess.run([NODE, '-e', js, str(katex), str(out)], capture_output=True, text=True, timeout=120)
print('chapter 8 blocks', start, end, bl[start]['source'][:60])
print(r.stdout[-8000:]); print(r.stderr[-2000:], file=sys.stderr)
