#!/usr/bin/env python3
"""Apply the string replacements of supp_plan.json to the original supplement files (orig_supp_chNN.tex), write
manuscript/supplements/chNN.tex and append one log entry per replaced block (anchor, old block, new block) to
step2/strict/rewrite_log.json.  Run after blocks.py, which rewrites the log."""
import json, re
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def blocks_of(text):
    parts = re.split(r'(?m)^(?=% ANCHOR:)', text)
    return [p for p in parts if p.strip()]


def main():
    plan = json.loads((HERE / 'supp_plan.json').read_text())
    logp = ROOT / 'step2' / 'strict' / 'rewrite_log.json'
    log = json.loads(logp.read_text())
    log = [e for e in log if not e.get('chapter', '').startswith('supplements/')]
    for f, reps in plan.items():
        text = (HERE / f'orig_supp_{f}.tex').read_text()
        old_blocks = blocks_of(text)
        new = text
        for old, rep in reps:
            assert new.count(old) == 1, (f, old[:80])
            new = new.replace(old, rep)
        new_blocks = blocks_of(new)
        assert len(old_blocks) == len(new_blocks)
        for a, b in zip(old_blocks, new_blocks):
            if a != b:
                anchor = a.splitlines()[0].replace('% ANCHOR: ', '')
                log.append({'chapter': f'supplements/{f}', 'label': anchor, 'op': 'replace', 'old': a.rstrip('\n'),
                            'new': b.rstrip('\n')})
        (ROOT / 'manuscript' / 'supplements' / f'{f}.tex').write_text(new)
        print(f, 'replacements', len(reps))
    logp.write_text(json.dumps(log, ensure_ascii=False, indent=1))
    print('log entries', len(log))


if __name__ == '__main__':
    main()
