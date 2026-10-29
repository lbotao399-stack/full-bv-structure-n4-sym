#!/usr/bin/env python3
"""Prepare the blocks of newly added manuscript chapters for appending to the live reader.

Usage: prep_new_chapter.py CHAP LIVE_IDS_JSON OUT_JSON. Input: reader_build/document.json (full conversion of
ch00..chNN). Output: the blocks of chapter CHAP with stable IDs (b-NNNNN -> cCHAP-NNNNN, chapter_ids.py),
internal links to chapters >= 8 remapped the same way, sentence-split with the reader skill's
idempotent splitter, and every formula checked with the reader's KaTeX options.
"""
from __future__ import annotations
import json, os, re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import chapter_ids  # noqa: E402
SKILL = Path.home() / '.claude/skills/academic-audio-lecture/scripts'
NODE = '/Users/libotao/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node'
# the Documents/ copy of node_modules is partly evicted by iCloud; the npx cache copy comes first
NODE_PATH = '/Users/libotao/.npm/_npx/4b4c857f6efdfb61/node_modules:/Users/libotao/Documents/Codex/2026-10-22/new-chat/work/ch33_lecture/notion-reader/node_modules'
KATEX = ROOT / 'reader/notion-reader/app/vendor/katex/katex.min.js'


def main():
    chap = int(sys.argv[1])               # chapter number of the full conversion, >= 8
    live_ids_file = Path(sys.argv[2])     # JSON list of block IDs of the live document
    out = Path(sys.argv[3])
    doc = json.loads((ROOT / 'reader_build/document.json').read_text())
    ch = chapter_ids.chapter_of(doc['blocks'])
    remap_id = chapter_ids.make_remap(doc['blocks'])
    new = [b for b in doc['blocks'] if ch[b['id']] == chap]
    assert new and chap >= chapter_ids.FIRST_PREFIXED
    live_ids = set(json.loads(live_ids_file.read_text()))
    remap = {b['id']: remap_id(b['id']) for b in new}

    for b in new:
        b['id'] = remap[b['id']]
        if b.get('parentId'):
            b['parentId'] = remap_id(b['parentId'])
        for k in ('source', 'caption'):
            if k in b:
                b[k] = chapter_ids.fix_links(b[k], remap_id)
    targets = set(re.findall(r'\(#([\w-]+)\)', json.dumps([b.get('source', '') for b in new])))
    missing = sorted(t for t in targets if t not in live_ids and t not in set(remap.values()))
    clash = sorted(set(remap.values()) & live_ids)

    # sentence split with the skill's splitter (idempotent; never writes to a server)
    env = dict(os.environ, PATH=str(Path(NODE).parent) + ':' + os.environ['PATH'], NODE_PATH=NODE_PATH)
    os.environ.update(env)
    sys.path.insert(0, str(SKILL))
    import sentence_blocks as sb
    split, report = sb.transform_document({'id': 'tmp', 'title': 'tmp', 'blocks': new})
    again, report2 = sb.transform_document(split)
    assert not report2['changed'], 'splitter not idempotent'
    blocks = split['blocks']

    # KaTeX check with the reader's options
    maths = []
    for b in blocks:
        s = b.get('source', '')
        if b['type'] == 'equation':
            maths.append((b['id'], s, True))
        else:
            for m in re.finditer(r'(?<!\\)\$([^$]+?)(?<!\\)\$|\\\((.+?)\\\)', s, flags=re.S):
                maths.append((b['id'], m.group(1) or m.group(2), False))
    js = ("const k=require(process.argv[1]);const items=JSON.parse(require('fs').readFileSync(0,'utf8'));"
          "const bad=[];for(const [id,tex,d] of items){try{k.renderToString(tex,{displayMode:d,throwOnError:true,strict:'ignore',trust:false,output:'htmlAndMathml'})}catch(e){bad.push([id,tex.slice(0,120),String(e.message).slice(0,200)])}}"
          "process.stdout.write(JSON.stringify(bad))")
    r = subprocess.run([NODE, '-e', js, str(KATEX)], input=json.dumps(maths), capture_output=True, text=True, env=env)
    bad = json.loads(r.stdout or '[]')
    out.write_text(json.dumps({'blocks': blocks}, ensure_ascii=False, indent=1))
    print(json.dumps({'new_blocks_before_split': len(new), 'after_split': len(blocks),
                      'split_changed': len(report['changed']), 'split_skipped': report['skipped'],
                      'formulas_checked': len(maths), 'katex_errors': bad[:20], 'n_katex_errors': len(bad),
                      'missing_link_targets': missing, 'id_clashes': clash, 'node_stderr': r.stderr[-300:]}, indent=1))


if __name__ == '__main__':
    main()
