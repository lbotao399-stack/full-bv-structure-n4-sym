import re, os, glob, sys
MS = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/manuscript'
labels = set()
for f in glob.glob(os.path.join(MS, 'ch0*_*.tex')):
    labels |= set(re.findall(r'\\label\{([^}]*)\}', open(f).read()))
sections = {'0': ['0.1', '0.2', '0.3', '0.4'], '1': [f'1.{k}' for k in range(1, 8)], '2': [f'2.{k}' for k in range(1, 6)]}
allowed = {'begin', 'end', 'item', 'ref', 'eqref'}
ok = True
for ch in ('ch00', 'ch01', 'ch02'):
    s = open(os.path.join(MS, 'supplements', ch + '.tex')).read()
    lines = s.split('\n')
    blocks = re.split(r'^% ANCHOR: ', s, flags=re.M)[1:]
    print(ch, len(blocks), 'blocks')
    for b in blocks:
        head, body = b.split('\n', 1)
        kind, arg = head.split(' ', 1)
        if kind == 'end_of_section':
            good = arg in sections[ch[-1]]
        elif kind == 'after_label':
            good = arg in labels
        else:
            good = False
        envs = re.findall(r'\\begin\{(\w+)\}', body)
        top = body.strip()
        if not top.startswith('\\begin{remark}[Physics notation:') or not top.endswith('\\end{remark}'):
            good = False; print('  bad remark wrapper', head)
        if top.count('\\begin{remark}') != 1:
            good = False; print('  more than one remark', head)
        for e in envs:
            if e not in ('remark', 'itemize', 'aligned', 'pmatrix'):
                good = False; print('  bad env', e, head)
        for r in re.findall(r'\\(?:eq)?ref\{([^}]*)\}', body):
            if r not in labels:
                good = False; print('  missing label', r, head)
        for bad in ('$$', '\\label', '\\tag', '\\item[', '\\begin{equation}', '\\text', '\\newcommand', '\\cite', '\\footnote'):
            if bad in body:
                good = False; print('  forbidden', bad, head)
        # brace balance and $ parity outside \[ \]
        depth = 0
        for chh in body.replace('\\{', '').replace('\\}', ''):
            depth += (chh == '{') - (chh == '}')
            if depth < 0: break
        if depth != 0:
            good = False; print('  unbalanced braces', head)
        txt = re.sub(r'\\\[.*?\\\]', '', body, flags=re.S)
        if txt.replace('\\$', '').count('$') % 2:
            good = False; print('  odd $', head)
        if body.count('\\[') != body.count('\\]'):
            good = False; print('  unbalanced \\[ \\]', head)
        nlines = len([l for l in body.strip().split('\n') if l.strip()])
        print(f"  {'OK ' if good else 'BAD'} {kind} {arg}  ({nlines} source lines)  {re.search(r'Physics notation: ([^]]*)', body).group(1)}")
        ok &= good
print('ALL OK' if ok else 'PROBLEMS')
