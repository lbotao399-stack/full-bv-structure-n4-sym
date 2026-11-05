#!/usr/bin/env python3
"""Re-point the two Chapter-9 links inside the Chapter-7 remark block (added by fix_c7_internal.py) to the current
targets of thm:c9-P1 and prop:c9-scheme. Usage: fix_c7_links.py PORT. Revision-guarded; only that block changes."""
import json, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import merge_into_live as mil, convert, chapter_ids
port = int(sys.argv[1])
full = convert.Converter(); full.run(sorted(convert.MS.glob('ch[0-9][0-9]_*.tex')))
remap = chapter_ids.make_remap(full.blocks)
tgt = {l: (remap(b), n) for l, (b, n) in full.labels.items()}
st, state = mil.req(port, 'GET'); rev, doc = state['revision'], state['document']
A = 'sentence-8f1edb9425823d2cf3d04581'
blocks = [dict(b) for b in doc['blocks']]
k = next(i for i, b in enumerate(blocks) if b['id'] == A)
assert not blocks[k].get('formatting')
src = blocks[k]['source']
src = re.sub(r'Theorem <a href="#c9-\d+">[\d.]+</a>', f'Theorem <a href="#{tgt["thm:c9-P1"][0]}">{tgt["thm:c9-P1"][1]}</a>', src)
src = re.sub(r'Proposition <a href="#c9-\d+">[\d.]+</a>', f'Proposition <a href="#{tgt["prop:c9-scheme"][0]}">{tgt["prop:c9-scheme"][1]}</a>', src)
blocks[k]['source'] = src
assert [b for i, b in enumerate(blocks) if i != k] == [b for i, b in enumerate(doc['blocks']) if i != k]
d = dict(doc); d['blocks'] = blocks
st, res = mil.req(port, 'PUT', {'document': d, 'revision': rev}); assert st == 200, (st, res)
print({'revision_after': res['revision'], 'links': re.findall(r'href="#(c9-\d+)">([^<]+)<', src)})
