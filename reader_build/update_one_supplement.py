#!/usr/bin/env python3
"""Replace the blocks of ONE supplement entry (ID prefix, e.g. ph02-02-) in a reader, keeping every other block
byte-identical; learner highlights on a replaced block are carried to the new block of the same ID as a whole-block
highlight of the same style. Usage: update_one_supplement.py PORT PREFIX [--dry-run]"""
import json, sys, re
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import merge_into_live as mil

def plain(src):
    """Reader formatting text: math -> U+FFFC, markdown links -> their label, tags stripped."""
    maths = []
    def m(mt):
        maths.append(mt.group(1)); return '￼'
    t = re.sub(r'(?<!\\)\$([^$]+?)(?<!\\)\$', m, src)
    t = re.sub(r'\[([^\]]*)\]\(#[^)]*\)', r'\1', t)
    t = re.sub(r'<[^>]+>', '', t)
    return t, maths

def main():
    port, prefix = int(sys.argv[1]), sys.argv[2]; dry = '--dry-run' in sys.argv
    entries = json.loads((HERE / 'supplement_blocks.json').read_text())
    new = [b for e in entries for b in e['blocks'] if b['id'].startswith(prefix)
           or str(b.get('sentenceSplit', {}).get('originId', '')).startswith(prefix)]
    st, state = mil.req(port, 'GET'); assert st == 200
    rev, doc = state['revision'], state['document']; old = doc['blocks']
    (HERE.parent / f'reader/qa/backup_port{port}_rev{rev:06d}_before_{prefix.rstrip("-")}.json').write_text(json.dumps(state, ensure_ascii=False))
    mine = lambda b: b['id'].startswith(prefix) or str(b.get('sentenceSplit', {}).get('originId', '')).startswith(prefix)
    idx = [k for k, b in enumerate(old) if mine(b)]
    assert idx == list(range(idx[0], idx[-1] + 1)), 'entry blocks not contiguous'
    prev = {b['id']: b for b in old if mine(b)}
    for b in new:
        p = prev.get(b['id'])
        if p and p.get('formatting'):
            if p.get('source') == b.get('source'):
                b['formatting'] = p['formatting']
            else:
                style = p['formatting']['ranges'][0]['style']
                text, maths = plain(b['source'])
                text = text.rstrip('\n'); offs = [k for k, ch in enumerate(text) if ch == '￼']
                b['formatting'] = {'version': 1, 'text': text + '\n', 'ranges': [{'start': 0, 'end': len(text), 'quote': text,
                    'before': '', 'after': '\n', 'math': [{'offset': o, 'tex': t} for o, t in zip(offs, maths)], 'style': style}]}
    lost = [i for i, p in prev.items() if p.get('formatting') and not any(b['id'] == i and b.get('formatting') for b in new)]
    assert not lost, f'formatting would be lost on {lost}'
    blocks = old[:idx[0]] + new + old[idx[-1] + 1:]
    assert [b for b in blocks if not mine(b)] == [b for b in old if not mine(b)]
    print({'revision': rev, 'replaced': len(idx), 'new': len(new), 'carried_formatting': [b['id'] for b in new if b.get('formatting')]})
    if dry: return
    d = dict(doc); d['blocks'] = blocks
    st, res = mil.req(port, 'PUT', {'document': d, 'revision': rev}); assert st == 200, (st, res)
    print({'revision_after': res['revision']})

if __name__ == '__main__':
    main()
