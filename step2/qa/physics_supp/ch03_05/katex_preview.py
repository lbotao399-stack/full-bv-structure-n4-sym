#!/usr/bin/env python3
"""Convert manuscript/supplements/*.tex ("Physics notation" remarks with % ANCHOR lines) into reader blocks.

Each supplement file is a sequence of blocks:
    % ANCHOR: end_of_section <number>      or      % ANCHOR: after_label <label>
    \\begin{remark}[Physics notation: <topic>] ... \\end{remark}
Output: reader_build/supplement_blocks.json = [{"anchor": {...}, "file": ..., "blocks": [...]}, ...] with stable IDs
ph<file>-<k>-<n>, refs resolved against the full manuscript (chapter-8 IDs remapped to c8-), callout headers
"Physics notation (<topic>).", sentence-split with the reader skill's splitter, every formula KaTeX-checked.
"""
from __future__ import annotations
import json, os, re, subprocess, sys
from pathlib import Path

HERE = Path('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/reader_build')
OUTDIR = Path('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/physics_supp/ch03_05')
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import convert  # noqa: E402

SKILL = Path.home() / '.claude/skills/academic-audio-lecture/scripts'
NODE = '/Users/libotao/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node'
NODE_PATH = '/Users/libotao/Documents/Codex/2026-10-22/new-chat/work/ch33_lecture/notion-reader/node_modules'
KATEX = ROOT / 'reader/notion-reader/app/vendor/katex/katex.min.js'
CH8_FIRST = 1399  # first block number of chapter 8 in the full conversion; live IDs are c8-NNNNN


def remap_id(bid):
    m = re.fullmatch(r'b-(\d{5})', bid)
    return f'c8-{m.group(1)}' if m and int(m.group(1)) >= CH8_FIRST else bid


def parse_file(path):
    text = path.read_text(encoding='utf-8')
    parts = re.split(r'^%\s*ANCHOR:\s*(\S+)\s+(\S+)\s*$', text, flags=re.M)
    out = []
    for k in range(1, len(parts), 3):
        kind, value, tex = parts[k], parts[k + 1], parts[k + 2]
        tex = '\n'.join(l for l in tex.splitlines() if not l.lstrip().startswith('%')).strip()
        if tex:
            out.append({'kind': kind, 'value': value, 'tex': tex})
    return out


def main():
    files = sorted(convert.MS.glob('ch[0-9][0-9]_*.tex'))
    full = convert.Converter()
    full.run(files)
    labels = {lab: (remap_id(bid), num) for lab, (bid, num) in full.labels.items()}
    sup_dir = convert.MS / 'supplements'
    order = ['ch03.tex', 'ch04.tex', 'ch05.tex']
    entries, problems = [], []
    for name in order:
        p = sup_dir / name
        if not p.exists():
            continue
        tag = 'key' if name == 'key.tex' else name[2:4]
        for k, sup in enumerate(parse_file(p)):
            c = convert.Converter()
            c.labels = dict(labels)
            c.chap = 0
            counter = {'n': 0}

            def bid(prefix='b', tag=tag, k=k, counter=counter):
                counter['n'] += 1
                return f'ph{tag}-{k:02d}-{counter["n"]:03d}'
            c.bid = bid
            if '\\item[' in sup['tex'] or '\\label' in sup['tex'] or '\\begin{equation' in sup['tex']:
                problems.append({'file': name, 'k': k, 'problem': 'forbidden markup (\\item[, \\label or equation env)'})
            c.parse_body(sup['tex'])
            for b in c.blocks:
                for key in ('source', 'caption'):
                    if key in b:
                        b[key] = c.resolve_refs(b[key])
                if b['type'] == 'callout' and b.get('calloutKind') == 'remark':
                    b['source'] = re.sub(r'^<strong>Remark [\d.]+</strong> \(Physics notation: (.*)\)\.$',
                                         r'<strong>Physics notation</strong> (\1).', b['source'])
            unresolved = re.findall(r'\?\?[\w:.-]+', json.dumps(c.blocks))
            if unresolved:
                problems.append({'file': name, 'k': k, 'problem': f'unresolved refs {unresolved}'})
            entries.append({'file': name, 'anchor': {'kind': sup['kind'], 'value': sup['value']}, 'blocks': c.blocks})
    # sentence split and KaTeX check
    env = dict(os.environ, PATH=str(Path(NODE).parent) + ':' + os.environ['PATH'], NODE_PATH=NODE_PATH)
    os.environ.update(env)
    sys.path.insert(0, str(SKILL))
    import sentence_blocks as sb
    maths = []
    for e in entries:
        split, _ = sb.transform_document({'id': 'tmp', 'title': 'tmp', 'blocks': e['blocks']})
        e['blocks'] = split['blocks']
        for b in e['blocks']:
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
    ids = [b['id'] for e in entries for b in e['blocks']]
    if len(ids) != len(set(ids)):
        problems.append({'problem': 'duplicate supplement block IDs'})
    (OUTDIR / 'katex_preview_blocks.json').write_text(json.dumps(entries, ensure_ascii=False, indent=1))
    print(json.dumps({'supplements': len(entries), 'blocks': len(ids), 'formulas': len(maths), 'katex_errors': bad[:15],
                      'n_katex_errors': len(bad), 'problems': problems,
                      'anchors': [f"{e['file']}:{e['anchor']['kind']} {e['anchor']['value']}" for e in entries]}, indent=1))


if __name__ == '__main__':
    main()
