#!/usr/bin/env python3
from pathlib import Path
import argparse,json,os,subprocess
root=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--port',type=int,required=True);a=p.parse_args()
runtime=json.loads((root/'qa/runtime.json').read_text());env=os.environ.copy()
paths=runtime.get('node_modules',[]);paths=[paths] if isinstance(paths,str) else paths
env['NODE_PATH']=os.pathsep.join(paths)
raise SystemExit(subprocess.call([runtime['python'],str(root/'notion-reader/server.py'),'--port',str(a.port)],cwd=root,env=env))
