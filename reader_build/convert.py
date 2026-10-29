#!/usr/bin/env python3
"""Convert manuscript/chNN_*.tex (restricted markup, see manuscript/FORMAT.md)
into (1) reader document JSON blocks and (2) a compilable LaTeX book."""
from __future__ import annotations
import json, re, shutil, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / 'manuscript'
OUT = ROOT / 'reader_build'
THM = ['definition', 'proposition', 'theorem', 'lemma', 'corollary', 'example', 'remark', 'axiom', 'proof']
THM_NAME = {'definition': 'Definition', 'proposition': 'Proposition', 'theorem': 'Theorem', 'lemma': 'Lemma',
            'corollary': 'Corollary', 'example': 'Example', 'remark': 'Remark', 'axiom': 'External input',
            'proof': 'Proof'}


def braced(s, i):
    """Return (content, index after closing brace) for s[i] == '{'."""
    assert s[i] == '{', s[i:i + 30]
    depth, j = 0, i
    while j < len(s):
        if s[j] == '\\':
            j += 2
            continue
        if s[j] == '{':
            depth += 1
        elif s[j] == '}':
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
        j += 1
    raise ValueError('unbalanced braces near: ' + s[i:i + 80])


class Converter:
    def __init__(self):
        self.blocks, self.labels, self.n = [], {}, 0
        self.chap = 0
        self.counters = {}
        self.figdir = OUT / 'figures'
        self.pending_refs = []

    def bid(self, prefix='b'):
        self.n += 1
        return f'{prefix}-{self.n:05d}'

    # ---------- inline text ----------
    def inline(self, t):
        t = ' '.join(t.split())
        maths = []
        def stash(m):
            maths.append(m.group(0))
            return f'\x00{len(maths)-1}\x00'
        t = re.sub(r'(?<!\\)\$[^$]*\$', stash, t)
        t = self._inline_text(t)
        return re.sub(r'\x00(\d+)\x00', lambda m: maths[int(m.group(1))], t)

    def _inline_text(self, t):
        t = re.sub(r'\\emph\{([^{}]*)\}', r'<em>\1</em>', t)
        t = re.sub(r'\\textbf\{([^{}]*)\}', r'<strong>\1</strong>', t)
        t = re.sub(r'\\textit\{([^{}]*)\}', r'<em>\1</em>', t)
        t = t.replace('~', '\u00a0').replace('---', '\u2014').replace('--', '\u2013')
        t = t.replace("``", '\u201c').replace("''", '\u201d')
        return t

    def resolve_refs(self, t):
        def eqref(m):
            lab = m.group(1)
            if lab not in self.labels:
                return '(??' + lab + ')'
            bid, num = self.labels[lab]
            return f'[({num})](#{bid})'

        def ref(m):
            lab = m.group(1)
            if lab not in self.labels:
                return '??' + lab
            bid, num = self.labels[lab]
            return f'[{num}](#{bid})'
        t = re.sub(r'\\eqref\{([^}]*)\}', eqref, t)
        t = re.sub(r'\\ref\{([^}]*)\}', ref, t)
        return t

    def add(self, typ, source, parent=None, **kw):
        b = {'id': kw.pop('id', None) or self.bid(), 'type': typ, 'source': source}
        if parent:
            b['parentId'] = parent
        b.update(kw)
        self.blocks.append(b)
        return b

    def next_num(self, kind):
        key = 'eq' if kind == 'eq' else ('fig' if kind == 'fig' else ('tab' if kind == 'tab' else 'thm'))
        self.counters[key] = self.counters.get(key, 0) + 1
        return f'{self.chap}.{self.counters[key]}'

    # ---------- block parser ----------
    def parse_body(self, body, parent=None):
        """Split body into paragraphs / displays / envs / figs / tables."""
        i = 0
        para = []

        def flush():
            text = ' '.join(para).strip()
            para.clear()
            if text:
                self.add('paragraph', self.inline(text), parent)
        lines = body
        pos = 0
        token = re.compile(r'\\begin\{(equation\*?|aligned|tabular|itemize|' + '|'.join(THM) + r')\}|\\\[|\\fig\{|\\section\{|\\subsection\{|\\chapter\{|\n\s*\n')
        while pos < len(lines):
            m = token.search(lines, pos)
            if not m:
                para.append(lines[pos:])
                break
            para.append(lines[pos:m.start()])
            tok = m.group(0)
            if tok.strip() == '':
                flush()
                pos = m.end()
                continue
            flush()
            if tok.startswith('\\chapter'):
                title, e = braced(lines, m.end() - 1)
                self.chap += 1
                self.counters = {}
                self.add('heading', f'{self.chap}. ' + self.inline(title), None, level=1)
                pos = e
            elif tok.startswith('\\section'):
                title, e = braced(lines, m.end() - 1)
                self.counters['sec'] = self.counters.get('sec', 0) + 1
                self.counters['sub'] = 0
                self.add('heading', f"{self.chap}.{self.counters['sec']} " + self.inline(title), None, level=2)
                pos = e
            elif tok.startswith('\\subsection'):
                title, e = braced(lines, m.end() - 1)
                self.counters['sub'] = self.counters.get('sub', 0) + 1
                self.add('heading', f"{self.chap}.{self.counters.get('sec', 0)}.{self.counters['sub']} " + self.inline(title), None, level=3)
                pos = e
            elif tok == '\\[':
                e = lines.index('\\]', m.end())
                self.add('equation', lines[m.end():e].strip(), parent)
                pos = e + 2
            elif tok.startswith('\\fig'):
                path, e = braced(lines, m.end() - 1)
                cap, e2 = braced(lines, e)
                num = self.next_num('fig')
                src = ROOT / path.strip()
                dest = self.figdir / src.name
                self.figdir.mkdir(parents=True, exist_ok=True)
                if src.exists():
                    shutil.copy2(src, dest)
                else:
                    print('MISSING FIGURE', src, file=sys.stderr)
                cap_html = f'<strong>Figure {num}.</strong> ' + self.inline(cap)
                b = self.add('image', 'figures/' + src.name, parent, image='figures/' + src.name, caption=cap_html)
                labm = re.search(r'\\label\{([^}]*)\}', cap)
                if labm:
                    self.labels[labm.group(1)] = (b['id'], num)
                    b['caption'] = re.sub(r'\\label\{[^}]*\}', '', b['caption'])
                pos = e2
            else:
                env = m.group(1)
                end_tok = '\\end{' + env + '}'
                # find matching end (no nesting of same env assumed)
                e = lines.index(end_tok, m.end())
                inner = lines[m.end():e]
                pos = e + len(end_tok)
                if env.startswith('equation'):
                    tex = inner.strip()
                    labm = re.search(r'\\label\{([^}]*)\}', tex)
                    tex = re.sub(r'\\label\{[^}]*\}', '', tex).strip()
                    num = self.next_num('eq')
                    b = self.add('equation', tex + f' \\tag{{{num}}}', parent)
                    if labm:
                        self.labels[labm.group(1)] = (b['id'], num)
                elif env == 'tabular':
                    spec, k = braced(inner, inner.index('{')) if inner.lstrip().startswith('{') else ('', 0)
                    rows_src = inner[k:]
                    rows = [r.strip() for r in re.split(r'\\\\', rows_src) if r.strip()]
                    html, first = ['<table>'], True
                    for r in rows:
                        r = r.replace('\\hline', '').strip()
                        if not r:
                            continue
                        cells = [self.inline(c.strip()) for c in re.split(r'(?<!\\)&', r)]
                        tag = 'th' if first else 'td'
                        html.append('<tr>' + ''.join(f'<{tag}>{c}</{tag}>' for c in cells) + '</tr>')
                        first = False
                    html.append('</table>')
                    capm = re.match(r'\s*\\tabcaption\{', lines[pos:])
                    cap = None
                    if capm:
                        cap, pos = braced(lines, pos + capm.end() - 1)
                    tb = self.add('table', '\n'.join(html), parent)
                    if cap:
                        num = self.next_num('tab')
                        labm = re.search(r'\\label\{([^}]*)\}', cap)
                        if labm:
                            self.labels[labm.group(1)] = (tb['id'], num)
                        cap = re.sub(r'\\label\{[^}]*\}', '', cap)
                        self.add('paragraph', f'<strong>Table {num}.</strong> ' + self.inline(cap), parent)
                elif env == 'itemize':
                    for item in [x for x in re.split(r'\\item\b', inner) if x.strip()]:
                        start = len(self.blocks)
                        self.parse_body(item.strip(), parent)
                        for b in self.blocks[start:]:
                            if b['type'] == 'paragraph' and b.get('parentId') == parent:
                                b['source'] = '\u2022 ' + b['source']
                                break
                        else:
                            self.add('paragraph', '\u2022', parent)
                else:  # theorem-like
                    title = ''
                    s = inner.lstrip()
                    if s.startswith('['):
                        depth, j = 0, 0
                        for j, ch in enumerate(s):
                            if ch == '[':
                                depth += 1
                            elif ch == ']':
                                depth -= 1
                                if depth == 0:
                                    break
                        title, s = s[1:j], s[j + 1:]
                    labm = re.search(r'\\label\{([^}]*)\}', s[:200])
                    if env == 'proof':
                        head = '<strong>Proof.</strong>' + (f' ({self.inline(title)})' if title else '')
                        num = None
                    else:
                        num = self.next_num('thm')
                        head = f'<strong>{THM_NAME[env]} {num}</strong>' + (f' ({self.inline(title)}).' if title else '.')
                    cb = self.add('callout', head, parent, calloutKind=env)
                    if labm and num:
                        self.labels[labm.group(1)] = (cb['id'], num)
                        s = s.replace(labm.group(0), '', 1)
                    start = len(self.blocks)
                    self.parse_body(s, parent=cb['id'])
                    if env == 'proof':
                        # append end-of-proof mark to last child paragraph or add one
                        kids = [b for b in self.blocks[start:] if b.get('parentId') == cb['id']]
                        if kids and kids[-1]['type'] == 'paragraph':
                            kids[-1]['source'] += ' \u220e'
                        else:
                            self.add('paragraph', '\u220e', cb['id'])
        flush()

    def run(self, files):
        for f in files:
            text = f.read_text(encoding='utf-8')
            text = '\n'.join(l for l in text.splitlines() if not l.lstrip().startswith('%'))
            mnum = re.match(r'ch(\d\d)_', f.name)
            if mnum:
                # \chapter increments the counter: chNN becomes chapter NN, as in the PDF
                self.chap = int(mnum.group(1)) - 1
            self.parse_body(text)
        for b in self.blocks:
            for key in ('source', 'caption'):
                if key in b and b['type'] != 'equation':
                    b[key] = self.resolve_refs(b[key])
            if b['type'] == 'equation':
                b['source'] = self.resolve_refs(b['source'])
        return self.blocks


def main():
    files = sorted(MS.glob('ch[0-9][0-9]_*.tex'))
    if len(sys.argv) > 1:
        files = [MS / f for f in sys.argv[1:]]
    c = Converter()
    blocks = c.run(files)
    doc = {'id': 'n4-costello-bv-step1', 'title': 'N=4 super Yang–Mills in Costello’s formalism — Step 1', 'blocks': blocks}
    (OUT / 'document.json').write_text(json.dumps(doc, ensure_ascii=False, indent=1))
    unresolved = sorted(set(re.findall(r'\?\?[\w:.-]+', json.dumps(doc))))
    print(json.dumps({'blocks': len(blocks), 'labels': len(c.labels), 'unresolved': unresolved[:40]}, indent=1))


if __name__ == '__main__':
    main()
