import json, re, sys, os, subprocess
from pathlib import Path
ROOT = Path('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27')
sys.path.insert(0, str(ROOT/'reader_build'))
import convert
draft = Path(sys.argv[1])
files = sorted(p for p in convert.MS.glob('ch[0-9][0-9]_*.tex') if not p.name.startswith('ch12'))
files.append(draft)
c = convert.Converter()
c.figdir = ROOT/'step2/strict/ch12.nosync/check/fig'
blocks = c.run(files)
unres = sorted(set(re.findall(r'\?\?[\w:.-]+', json.dumps(blocks))))
# chapter 12 blocks
start = None
for i,b in enumerate(blocks):
    if b['type']=='heading' and b.get('level')==1 and b['source'].startswith('12.'):
        start=i
ch = blocks[start:]
labs = [l for l in c.labels if l.startswith(('def:c12','thm:c12','lem:c12','prop:c12','cor:c12','rem:c12','eq:c12','tab:c12'))]
maths=[]
for b in ch:
    s=b.get('source','')
    if b['type']=='equation': maths.append((b['id'],s,True))
    else:
        for m in re.finditer(r'(?<!\\)\$([^$]+?)(?<!\\)\$', s, flags=re.S): maths.append((b['id'],m.group(1),False))
    if b.get('caption'):
        for m in re.finditer(r'(?<!\\)\$([^$]+?)(?<!\\)\$', b['caption'], flags=re.S): maths.append((b['id'],m.group(1),False))
NODE='/Users/libotao/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node'
KATEX=ROOT/'reader/notion-reader/app/vendor/katex/katex.min.js'
js=("const k=require(process.argv[1]);const items=JSON.parse(require('fs').readFileSync(0,'utf8'));"
    "const bad=[];for(const [id,tex,d] of items){try{k.renderToString(tex,{displayMode:d,throwOnError:true,strict:'ignore',trust:false,output:'htmlAndMathml'})}catch(e){bad.push([id,tex.slice(0,100),String(e.message).slice(0,160)])}}"
    "process.stdout.write(JSON.stringify(bad))")
r=subprocess.run([NODE,'-e',js,str(KATEX)],input=json.dumps(maths),capture_output=True,text=True,timeout=300)
bad=json.loads(r.stdout or '[]')
types={}
for b in ch: types[b['type']]=types.get(b['type'],0)+1
print(json.dumps({'ch12_blocks':len(ch),'types':types,'ch12_labels':len(labs),'unresolved':unres,'formulas':len(maths),'katex_errors':bad[:20],'n_katex':len(bad),'stderr':r.stderr[-300:]},indent=1))
