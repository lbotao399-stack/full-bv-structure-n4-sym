#!/usr/bin/env python3
"""Block-level rewrite driver for the manuscript chapters.

A chapter file is split into top-level units: theorem-like environments, proofs, tables (tabular + following
\\tabcaption line), top-level displays, headings and paragraphs.  Units are addressed by their first \\label;
a proof by 'proof:<label of the preceding unit>'; a table by 'table:<label of the preceding unit>' (k-th table
after that unit: 'table2:...').  Operations (JSON list) are applied to the ORIGINAL file (orig_*.tex), so a run
is idempotent; the result is written to manuscript/ and every old->new passage is logged.

ops entries: {"op": "replace"|"insert_after"|"delete", "key": ..., "file": "<snippet .tex>"}
 replace: the unit(s) 'key' (or the range key..key_end if "key_end" given) replaced by the snippet text
 insert_after: snippet inserted after unit 'key' (after its proof if "after_proof": true)
"""
import json, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ENVS = ('definition', 'proposition', 'theorem', 'lemma', 'corollary', 'example', 'remark', 'axiom', 'proof',
        'tabular', 'equation')


def split_units(text):
    lines = text.split('\n')
    units, i, n = [], 0, len(lines)
    while i < n:
        ln = lines[i]
        m = re.match(r'\\begin\{(' + '|'.join(ENVS) + r')\}', ln)
        if m:
            env = m.group(1)
            depth, j = 0, i
            while j < n:
                depth += len(re.findall(r'\\begin\{' + env + r'\}', lines[j]))
                depth -= len(re.findall(r'\\end\{' + env + r'\}', lines[j]))
                if depth == 0:
                    break
                j += 1
            k = j + 1
            if env == 'tabular' and k < n and lines[k].startswith('\\tabcaption'):
                k += 1
            units.append({'kind': env, 'text': '\n'.join(lines[i:k])})
            i = k
        elif ln.strip() == '':
            units.append({'kind': 'blank', 'text': ln})
            i += 1
        elif ln.startswith('\\chapter') or ln.startswith('\\section') or ln.startswith('\\subsection'):
            units.append({'kind': 'heading', 'text': ln})
            i += 1
        elif ln.startswith('\\['):
            j = i
            while j < n and '\\]' not in lines[j]:
                j += 1
            units.append({'kind': 'display', 'text': '\n'.join(lines[i:j + 1])})
            i = j + 1
        else:
            j = i
            while j < n and lines[j].strip() != '' and not re.match(r'\\begin\{(' + '|'.join(ENVS) + r')\}', lines[j]) \
                    and not lines[j].startswith('\\section') and not lines[j].startswith('\\['):
                j += 1
            units.append({'kind': 'para', 'text': '\n'.join(lines[i:j])})
            i = j
    # keys
    last, tabcount, parcount = None, {}, {}
    for u in units:
        lab = re.search(r'\\label\{([^}]*)\}', u['text'])
        if u['kind'] == 'proof':
            u['key'] = f'proof:{last}'
        elif u['kind'] == 'tabular' and not (lab and lab.group(1).startswith('tab:')):
            tabcount[last] = tabcount.get(last, 0) + 1
            c = tabcount[last]
            u['key'] = f'table:{last}' if c == 1 else f'table{c}:{last}'
        elif u['kind'] in ('para', 'heading'):
            parcount[last] = parcount.get(last, 0) + 1
            c = parcount[last]
            u['key'] = f'{u["kind"]}:{last}' if c == 1 else f'{u["kind"]}{c}:{last}'
        elif lab:
            u['key'] = lab.group(1)
            last = lab.group(1)
        else:
            u['key'] = None
        if u['kind'] in ('tabular',) and lab and lab.group(1).startswith('tab:'):
            u['key'] = lab.group(1)
    return units


def find(units, key):
    for i, u in enumerate(units):
        if u.get('key') == key:
            return i
    raise KeyError(key)


def apply(chapter, orig, ops, log):
    text = (HERE / orig).read_text()
    units = split_units(text)
    assert '\n'.join(u['text'] for u in units) == text, 'split not lossless'
    for op in ops:
        snip = (HERE / op['file']).read_text().rstrip('\n') if 'file' in op else ''
        if op['op'] == 'replace':
            i = find(units, op['key'])
            j = find(units, op['key_end']) if 'key_end' in op else i
            old = '\n'.join(u['text'] for u in units[i:j + 1])
            units[i:j + 1] = [{'kind': 'new', 'text': snip, 'key': op.get('newkey', op['key'])}]
            log.append({'chapter': chapter, 'label': op['key'] + (('..' + op['key_end']) if 'key_end' in op else ''),
                        'op': 'replace', 'old': old, 'new': snip})
        elif op['op'] == 'insert_after':
            i = find(units, op['key'])
            if op.get('after_proof'):
                i = find(units, 'proof:' + op['key'])
            units[i + 1:i + 1] = [{'kind': 'blank', 'text': ''}, {'kind': 'new', 'text': snip, 'key': op.get('newkey')}]
            log.append({'chapter': chapter, 'label': op.get('newkey', 'after ' + op['key']), 'op': 'insert_after ' + op['key'],
                        'old': '', 'new': snip})
        elif op['op'] == 'delete':
            i = find(units, op['key'])
            log.append({'chapter': chapter, 'label': op['key'], 'op': 'delete', 'old': units[i]['text'], 'new': ''})
            del units[i]
    return '\n'.join(u['text'] for u in units)


def labels(text):
    return set(re.findall(r'\\label\{([^}]*)\}', text))


def main():
    plan = json.loads((HERE / 'plan.json').read_text())
    log = []
    for ch in plan:
        new = apply(ch['chapter'], ch['orig'], ch['ops'], log)
        lost = labels((HERE / ch['orig']).read_text()) - labels(new)
        if lost:
            print('LOST LABELS', ch['chapter'], sorted(lost))
            sys.exit(1)
        (ROOT / ch['target']).write_text(new)
        print(ch['chapter'], 'ops', len(ch['ops']), 'labels', len(labels(new)))
    out = ROOT / 'step2' / 'strict' / 'rewrite_log.json'
    out.write_text(json.dumps(log, ensure_ascii=False, indent=1))
    print('log entries', len(log))


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'keys':
        for u in split_units((HERE / sys.argv[2]).read_text()):
            if u['key']:
                print(u['kind'], u['key'])
    else:
        main()
