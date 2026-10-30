"""Extract every math string of chapter 20 and its supplement (display equations and inline $...$),
as in step2/qa/ch10_review/katex_check.py, into items.json for the KaTeX run (katex_check.js, run with
node or with the macOS JavaScriptCore shell jsc)."""
import json, re, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..', '..'))
files = {'ch20': os.path.join(ROOT, 'manuscript', 'ch20_scalar_two_point.tex'),
         'supp20': os.path.join(ROOT, 'manuscript', 'supplements', 'ch20.tex')}
def tex_math(text, tag):
    text = '\n'.join(l for l in text.splitlines() if not l.lstrip().startswith('%'))
    out = []
    for m in re.finditer(r'\\\[(.*?)\\\]', text, flags=re.S):
        out.append(('display', m.group(1), tag))
    for m in re.finditer(r'\\begin\{equation\}(.*?)\\end\{equation\}', text, flags=re.S):
        out.append(('display', re.sub(r'\\label\{[^}]*\}', '', m.group(1)), tag))
    t2 = re.sub(r'\\\[(.*?)\\\]', ' ', text, flags=re.S)
    t2 = re.sub(r'\\begin\{equation\}(.*?)\\end\{equation\}', ' ', t2, flags=re.S)
    for m in re.finditer(r'(?<!\\)\$([^$]*)\$', t2):
        out.append(('inline', m.group(1), tag))
    return out
items = []
for tag, p in files.items():
    items += tex_math(open(p).read(), tag)
json.dump(items, open(os.path.join(HERE, 'items.json'), 'w'))
print(len(items), 'math items')
