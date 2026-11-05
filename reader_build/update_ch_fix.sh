#!/bin/bash
# Replace chapters 8, 11, 12 in the reader at port $1 and update the changed supplement entries one by one.
set -e -o pipefail
cd "$(dirname "$0")/.."
for c in 8 9 10 11 12; do python3 reader_build/merge_into_live.py $1 reader_build/ch${c}_blocks.json --replace-prefix c$c | tail -2; done
python3 - "$1" <<'PY'
import json,sys,subprocess,urllib.request
port=sys.argv[1]
s=json.loads(urllib.request.urlopen(f'http://127.0.0.1:{port}/api/document').read()); b=s['document']['blocks']
def org(x): return str(x.get('sentenceSplit',{}).get('originId') or x['id'])
new=json.load(open('reader_build/supplement_blocks.json'))
for e in new:
    pre=e['blocks'][0]['id'].rsplit('-',1)[0]+'-'
    live=[x for x in b if org(x).startswith(pre) or x['id'].startswith(pre)]
    lv=[{k:v for k,v in x.items() if k!='formatting'} for x in live]
    if len(live)!=len(e['blocks']) or any(l!=n for l,n in zip(lv,e['blocks'])):
        if not live: print('NEW entry (needs insert_supplements):',pre); continue
        print(subprocess.run(['python3','reader_build/update_one_supplement.py',port,pre],capture_output=True,text=True).stdout.strip())
PY
