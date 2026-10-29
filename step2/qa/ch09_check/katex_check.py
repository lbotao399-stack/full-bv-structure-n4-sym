"""Render every math string of chapter 9 (reader blocks) with the vendored KaTeX; report parse errors."""
import json, re, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
doc = json.loads((ROOT / 'reader_build' / 'document.json').read_text())
bl = doc['blocks']
start = [i for i, b in enumerate(bl) if b['type'] == 'heading' and b.get('level') == 1 and b['source'].startswith('9.')][0]
items = []
for b in bl[start:]:
    s = b.get('source', '')
    if b['type'] == 'equation':
        items.append(('display', s))
    else:
        for m in re.finditer(r'(?<!\\)\$([^$]*)\$', s):
            items.append(('inline', m.group(1)))
out = Path(__file__).with_name('math_items.json')
out.write_text(json.dumps(items))
katex = ROOT / 'reader' / 'notion-reader' / 'app' / 'vendor' / 'katex' / 'katex.js'
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
r = subprocess.run(['node', '-e', js, str(katex), str(out)], capture_output=True, text=True, timeout=120)
print(r.stdout[-6000:]); print(r.stderr[-2000:], file=sys.stderr)
