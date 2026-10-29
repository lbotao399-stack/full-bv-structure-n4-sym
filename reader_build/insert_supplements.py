#!/usr/bin/env python3
"""Insert the prepared "Physics notation" supplements into a reader document with a revision-guarded PUT.

Usage: insert_supplements.py PORT [--dry-run]
Anchors: end_of_section <n> -> before the next heading of the same or higher level after heading "<n> ...";
after_label <label> -> right after the labelled block and its subtree (and after a following proof callout).
Supplement blocks already present (IDs ph*) are replaced; the script refuses if any of them carries learner
formatting. Every other block must come out byte-identical and in the same order.
"""
from __future__ import annotations
import json, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import merge_into_live as mil  # noqa: E402  (req helper)
import convert  # noqa: E402
import chapter_ids  # noqa: E402


def is_supp(b):
    return any(str(v).startswith('ph') for v in (b['id'], b.get('parentId', ''), b.get('sentenceSplit', {}).get('originId', '')))


def main():
    port = int(sys.argv[1]); dry = '--dry-run' in sys.argv
    entries = json.loads((HERE / 'supplement_blocks.json').read_text())
    st, state = mil.req(port, 'GET'); assert st == 200, state
    rev, doc = state['revision'], state['document']
    backup = ROOT / f'reader/qa/backup_port{port}_rev{rev:06d}_before_supplements.json'
    backup.write_text(json.dumps(state, ensure_ascii=False))
    old = doc['blocks']
    stale = [b for b in old if is_supp(b)]
    if any(b.get('formatting') for b in stale):
        raise SystemExit('learner formatting on existing supplement blocks; resolve by hand')
    base = [b for b in old if not is_supp(b)]
    full = convert.Converter(); full.run(sorted(convert.MS.glob('ch[0-9][0-9]_*.tex')))
    remap_id = chapter_ids.make_remap(full.blocks)
    labels = {lab: remap_id(bid) for lab, (bid, _) in full.labels.items()}
    idx = {b['id']: i for i, b in enumerate(base)}
    parent = {b['id']: b.get('parentId') for b in base}

    def descends(bid, root):
        seen = 0
        while bid and seen < 50:
            if bid == root:
                return True
            bid = parent.get(bid); seen += 1
        return False

    def position(anchor):
        kind, val = anchor['kind'], anchor['value']
        if kind == 'end_of_section':
            level = val.count('.') + 1
            h = next((i for i, b in enumerate(base) if b['type'] == 'heading' and b.get('level') == level
                      and re.match(re.escape(val) + r'\s', b['source'])), None)
            if h is None:
                raise SystemExit(f'heading {val} not found')
            j = h + 1
            while j < len(base) and not (base[j]['type'] == 'heading' and base[j].get('level', 9) <= level and not base[j].get('parentId')):
                j += 1
            return j
        if kind == 'after_label':
            bid = labels.get(val)
            if bid is None or bid not in idx:
                raise SystemExit(f'label {val} not found in the document')
            j = idx[bid] + 1
            while j < len(base) and descends(base[j].get('parentId'), bid):
                j += 1
            if j < len(base) and base[j]['type'] == 'callout' and base[j].get('calloutKind') == 'proof':
                pid = base[j]['id']; j += 1
                while j < len(base) and descends(base[j].get('parentId'), pid):
                    j += 1
            return j
        raise SystemExit(f'unknown anchor {anchor}')

    inserts = {}
    for e in entries:
        inserts.setdefault(position(e['anchor']), []).extend(e['blocks'])
    blocks = []
    for i, b in enumerate(base):
        blocks.extend(inserts.get(i, []))
        blocks.append(b)
    blocks.extend(inserts.get(len(base), []))
    ids = [b['id'] for b in blocks]
    assert len(ids) == len(set(ids)), 'duplicate IDs'
    kept = [b for b in blocks if not is_supp(b)]
    assert kept == base, 'a pre-existing block changed or moved'
    summary = {'revision_before': rev, 'blocks_before': len(old), 'blocks_after': len(blocks), 'replaced_old_supplement_blocks': len(stale),
               'supplements': len(entries), 'learner_formatting_blocks': sum(1 for b in old if b.get('formatting')), 'backup': str(backup)}
    if dry:
        print(json.dumps({**summary, 'dry_run': True}, indent=1)); return
    doc = dict(doc); doc['blocks'] = blocks
    st, res = mil.req(port, 'PUT', {'document': doc, 'revision': rev}); assert st == 200, (st, res)
    st, back = mil.req(port, 'GET')
    bb = back['document']['blocks']
    print(json.dumps({**summary, 'revision_after': res['revision'], 'readback_ids_ok': [b['id'] for b in bb] == ids,
                      'non_supplement_blocks_identical': [b for b in bb if not is_supp(b)] == base,
                      'learner_formatting_blocks_after': sum(1 for b in bb if b.get('formatting'))}, indent=1))


if __name__ == '__main__':
    main()
