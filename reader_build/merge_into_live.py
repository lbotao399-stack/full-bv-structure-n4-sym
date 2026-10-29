#!/usr/bin/env python3
"""Append prepared chapter blocks to the live reader document with a revision-guarded PUT.

Usage: merge_into_live.py PORT BLOCKS_JSON [--replace-prefix PREFIX] [--dry-run]
- Fetches the latest saved document and revision, backs it up under reader/qa/,
- appends the prepared blocks at the end (or, with --replace-prefix, replaces the blocks whose IDs
  start with PREFIX- in place, keeping any learner formatting on blocks whose source is unchanged),
- checks that every pre-existing block is untouched, PUTs with the fetched revision, reads back.
"""
from __future__ import annotations
import json, sys, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def req(port, method, body=None):
    data = None if body is None else json.dumps(body, ensure_ascii=False).encode()
    r = urllib.request.Request(f'http://127.0.0.1:{port}/api/document', data=data, method=method,
                               headers={'Content-Type': 'application/json'})
    try:
        with urllib.request.urlopen(r, timeout=60) as f:
            return f.status, json.loads(f.read())
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read() or b'{}')


def main():
    port = int(sys.argv[1]); new = json.loads(Path(sys.argv[2]).read_text())['blocks']
    prefix = None; dry = '--dry-run' in sys.argv
    if '--replace-prefix' in sys.argv:
        prefix = sys.argv[sys.argv.index('--replace-prefix') + 1]
    st, state = req(port, 'GET')
    assert st == 200, state
    rev, doc = state['revision'], state['document']
    backup = ROOT / f'reader/qa/backup_port{port}_rev{rev:06d}_before_{Path(sys.argv[2]).stem}.json'
    backup.write_text(json.dumps(state, ensure_ascii=False))
    old = doc['blocks']
    def mine(b):
        return any(str(v).startswith(prefix + '-') for v in (b['id'], b.get('parentId', ''),
                                                             b.get('sentenceSplit', {}).get('originId', '')))
    if prefix:
        keep = [b for b in old if not mine(b)]
        prev = {b['id']: b for b in old if mine(b)}
        for b in new:   # carry learner formatting over where the block text is unchanged
            p = prev.get(b['id'])
            if p and p.get('source') == b.get('source') and p.get('formatting'):
                b['formatting'] = p['formatting']
        lost = [i for i, p in prev.items() if p.get('formatting') and not any(b['id'] == i and b.get('formatting') for b in new)]
        if lost:
            raise SystemExit(f'learner formatting would be lost on {lost}; resolve by hand')
        first = next((k for k, b in enumerate(old) if mine(b)), len(old))
        head = [b for b in old[:first]]
        tail = [b for b in old[first:] if not mine(b)]
        blocks = head + new + tail
        untouched = keep
    else:
        ids = {b['id'] for b in old}
        clash = [b['id'] for b in new if b['id'] in ids]
        assert not clash, f'ID clash {clash[:5]}'
        blocks = old + new
        untouched = old
    after = {b['id']: b for b in blocks}
    assert all(after[b['id']] == b for b in untouched), 'a pre-existing block would change'
    assert len(after) == len(blocks), 'duplicate IDs'
    assert 1 <= len(blocks) <= 10000
    summary = {'revision_before': rev, 'blocks_before': len(old), 'blocks_after': len(blocks), 'backup': str(backup),
               'learner_formatting_blocks': sum(1 for b in old if b.get('formatting'))}
    if dry:
        print(json.dumps({**summary, 'dry_run': True}, indent=1)); return
    doc = dict(doc); doc['blocks'] = blocks
    st, res = req(port, 'PUT', {'document': doc, 'revision': rev})
    assert st == 200, (st, res)
    st, back = req(port, 'GET')
    ok = back['revision'] == res['revision'] and [b['id'] for b in back['document']['blocks']] == [b['id'] for b in blocks]
    formatting_kept = all({b['id']: b for b in back['document']['blocks']}[b['id']].get('formatting') == b.get('formatting')
                          for b in untouched if b.get('formatting'))
    print(json.dumps({**summary, 'revision_after': res['revision'], 'readback_ok': ok, 'learner_formatting_kept': formatting_kept}, indent=1))


if __name__ == '__main__':
    main()
