"""Render every math string of the Chapter 6-8 physics-notation supplements with the vendored KaTeX."""
import json, re, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[4]
items = []
for ch in ('ch06', 'ch07', 'ch08'):
    s = (ROOT / 'manuscript' / 'supplements' / f'{ch}.tex').read_text()
    s = re.sub(r'^%.*$', '', s, flags=re.M)
    for m in re.finditer(r'\\\[(.*?)\\\]', s, flags=re.S):
        items.append((ch, 'display', m.group(1)))
    rest = re.sub(r'\\\[(.*?)\\\]', ' ', s, flags=re.S)
    for m in re.finditer(r'(?<!\\)\$([^$]*)\$', rest):
        items.append((ch, 'inline', m.group(1)))
    # bracket titles of remarks contain math too (already caught by the $ scan)
out = Path(__file__).with_name('math_items_supp.json')
out.write_text(json.dumps([[k, s] for _, k, s in items]))
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
