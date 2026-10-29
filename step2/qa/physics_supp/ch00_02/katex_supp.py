"""Render every math string of the Chapter 0-2 supplements with the vendored KaTeX of the reader."""
import json, re, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[4]
items = []
for ch in ('ch00', 'ch01', 'ch02'):
    s = (ROOT / 'manuscript' / 'supplements' / f'{ch}.tex').read_text()
    for m in re.finditer(r'\\\[(.*?)\\\]', s, flags=re.S):
        items.append((ch, 'display', m.group(1).strip()))
    t = re.sub(r'\\\[(.*?)\\\]', ' ', s, flags=re.S)
    for m in re.finditer(r'(?<!\\)\$([^$]*)\$', t):
        items.append((ch, 'inline', m.group(1)))
out = Path(__file__).with_name('katex_items.json')
out.write_text(json.dumps(items))
katex = ROOT / 'reader' / 'notion-reader' / 'app' / 'vendor' / 'katex' / 'katex.js'
js = r"""
const katex = require(process.argv[1]);
const items = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));
let bad = 0;
items.forEach(([ch, kind, s], i) => {
  try { katex.renderToString(s, {displayMode: kind === 'display', throwOnError: true, strict: 'ignore'}); }
  catch (e) { bad++; console.log('ERR', ch, i, kind, String(e.message).slice(0, 200), '::', s.slice(0, 160)); }
});
console.log(JSON.stringify({items: items.length, errors: bad}));
"""
r = subprocess.run(['node', '-e', js, str(katex), str(out)], capture_output=True, text=True, timeout=120)
print(r.stdout[-6000:]); print(r.stderr[-2000:], file=sys.stderr)
