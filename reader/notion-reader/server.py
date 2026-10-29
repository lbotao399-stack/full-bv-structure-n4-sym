#!/usr/bin/env python3
"""Local editable Notion manuscript, measured narration, and cross-window playback."""
from pathlib import Path
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlsplit, unquote
import argparse, copy, hashlib, itertools, json, math, os, queue, re, shutil, subprocess, threading, time
from authoring.selection import selection
from authoring.source_language import source_language
from authoring.dictionary import lookup as dictionary_lookup
from image_store import MAX_JSON_BYTES, capabilities as image_capabilities, store_image

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data';DATA.mkdir(exist_ok=True)
HISTORY=DATA/'history';HISTORY.mkdir(exist_ok=True)
JOBS=ROOT/'audio/jobs';JOBS.mkdir(parents=True,exist_ok=True)
RUNTIME=json.loads((ROOT.parent/'qa/runtime.json').read_text())
NODE_ENV=os.environ.copy()
node_paths=RUNTIME.get('node_modules',[])
if isinstance(node_paths,str):node_paths=[node_paths]
if node_paths:NODE_ENV['NODE_PATH']=os.pathsep.join(node_paths)
SYNTHESIS_IDENTITY={name:hashlib.sha256(Path(path).read_bytes()).hexdigest() for name,path in {'model':RUNTIME['model'],'voices':RUNTIME['voices'],'worker':ROOT/'authoring/tts_worker.py'}.items()}
SYNTHESIS_IDENTITY['speechText']=hashlib.sha256((ROOT/'app/speech-text.js').read_bytes()).hexdigest()
LOCK=threading.RLock();PLAN_LOCK=threading.Lock();JOB_LOCK=threading.RLock()
STATE_PATH=DATA/'document.json'
state=json.loads(STATE_PATH.read_text()) if STATE_PATH.exists() else {'document':json.loads((ROOT/'document.seed.json').read_text()),'revision':0}
player={'owner':None,'blockId':None,'playing':False,'time':0,'updatedAt':0}
jobs={}; work=queue.PriorityQueue(); order=itertools.count(); planner=None; synthesizer=None
LOG=open(ROOT/'qa/backend-workers.log','a',buffering=1)

def atomic(path,value):
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2));tmp.replace(path)

def validate_document(doc):
    if not isinstance(doc,dict) or not isinstance(doc.get('title'),str) or not isinstance(doc.get('blocks'),list):raise ValueError('Expected a document title and blocks.')
    if not 1<=len(doc['blocks'])<=10000:raise ValueError('The document must contain 1 to 10000 blocks.')
    seen=set()
    for b in doc['blocks']:
        validate_block(b)
        if b['id'] in seen:raise ValueError('Block IDs must be unique.')
        seen.add(b['id'])

def validate_block(b):
    if not isinstance(b,dict) or not isinstance(b.get('id'),str) or not re.fullmatch(r'[\w-]{1,100}',b['id']):raise ValueError('Invalid block ID.')
    if b.get('type') not in {'paragraph','heading','equation','image','callout','table'}:raise ValueError('Unknown block type.')
    if not isinstance(b.get('source'),str) or len(b['source'])>100000:raise ValueError('Invalid block source.')
    if b.get('speechOverride') is not None and (not isinstance(b['speechOverride'],str) or len(b['speechOverride'])>100000):raise ValueError('Invalid spoken reading.')
    if b['type']=='image' and not (ROOT/b['source']).resolve().is_relative_to(ROOT/'figures'):raise ValueError('Images must come from the local figure library.')

def make_plan(block):
    global planner
    with PLAN_LOCK:
        if planner is None or planner.poll() is not None:
            planner=subprocess.Popen([RUNTIME['node'],str(ROOT/'authoring/speech_plan.cjs')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=LOG,text=True,bufsize=1,env=NODE_ENV)
        planner.stdin.write(json.dumps({'block':block,'language':source_language(block,state['document'])})+'\n');planner.stdin.flush()
        line=planner.stdout.readline()
        if not line:raise RuntimeError('The formula reader stopped. Retry this block.')
        result=json.loads(line)
        if not result.get('ok'):raise ValueError(result.get('error','Unable to read this formula.'))
        return result

def prepare(block, priority=0):
    validate_block(block)
    plan=make_plan(block)
    key=hashlib.sha256(json.dumps({'segments':plan['segments'],'language':plan['language'],'version':'source-language-kokoro-v2','synthesis':SYNTHESIS_IDENTITY},sort_keys=True).encode()).hexdigest()[:32]
    with JOB_LOCK:
        if key not in jobs:
            saved=JOBS/(key+'.json')
            if saved.exists() and (ROOT/'audio'/(key+'.wav')).exists():jobs[key]=json.loads(saved.read_text())
            else:
                jobs[key]={'status':'queued','speech':plan['speech'],'formulaReadings':plan['formulaReadings'],'jobId':key,'queuedAt':time.time()}
                work.put((priority,next(order),key,plan))
        elif jobs[key]['status']=='queued' and priority==0:
            # A foreground request may jump ahead of background prewarming.
            work.put((0,next(order),key,plan))
        elif jobs[key]['status']=='error':
            jobs[key].update(status='queued');work.put((priority,next(order),key,plan))
    return key

def worker():
    global synthesizer
    while True:
        _,_,key,plan=work.get()
        with JOB_LOCK:
            if jobs[key]['status']!='queued':work.task_done();continue
            jobs[key]['status']='running'
        try:
            if synthesizer is None or synthesizer.poll() is not None:
                synthesizer=subprocess.Popen([RUNTIME['python'],str(ROOT/'authoring/tts_worker.py')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=LOG,text=True,bufsize=1)
            synthesizer.stdin.write(json.dumps({'key':key,'segments':plan['segments'],'language':plan['language']})+'\n');synthesizer.stdin.flush()
            line=synthesizer.stdout.readline()
            if not line:raise RuntimeError('Narration worker stopped. Retry this block.')
            result=json.loads(line)
            if not result.pop('ok',False):raise RuntimeError(result.get('error','Narration failed.'))
            with JOB_LOCK:
                jobs[key].update(result,status='ready');atomic(JOBS/(key+'.json'),jobs[key])
        except Exception as e:
            with JOB_LOCK:jobs[key].update(status='error',error=str(e))
        finally:work.task_done()

def player_state():
    result=copy.copy(player)
    if time.time()*1000-result['updatedAt']>10000:result['playing']=False
    return result

class Handler(SimpleHTTPRequestHandler):
    def __init__(self,*args,**kwargs):super().__init__(*args,directory=str(ROOT),**kwargs)
    def log_message(self,fmt,*args):
        if '/api/jobs/' not in str(args) and '/api/player-state' not in str(args):super().log_message(fmt,*args)
    def allowed(self):
        expected={f'127.0.0.1:{self.server.server_port}',f'localhost:{self.server.server_port}'}
        return self.headers.get('Host') in expected
    def json_response(self,value,status=200):
        body=json.dumps(value,ensure_ascii=False).encode();self.send_response(status);self.send_header('Content-Type','application/json; charset=utf-8');self.send_header('Content-Length',str(len(body)));self.send_header('Cache-Control','no-store');self.end_headers()
        if self.command!='HEAD':self.wfile.write(body)
    def body(self,max_size=8*1024*1024):
        if not self.allowed():raise PermissionError('Invalid local host.')
        origin=self.headers.get('Origin')
        if origin and origin not in {f'http://127.0.0.1:{self.server.server_port}',f'http://localhost:{self.server.server_port}'}:raise PermissionError('Use the local reader to edit this document.')
        if self.headers.get_content_type()!='application/json':raise ValueError('Expected JSON.')
        size=int(self.headers.get('Content-Length','0'))
        if not 0<size<=max_size:raise ValueError('Request is too large or empty.')
        return json.loads(self.rfile.read(size))
    def do_GET(self):
        if not self.allowed():return self.json_response({'error':'Invalid local host.'},403)
        route=urlsplit(self.path).path
        if route=='/api/health':return self.json_response({'ok':True,'service':'notion-reader','blocks':len(state['document']['blocks']),'jobs':{s:sum(j['status']==s for j in list(jobs.values())) for s in ['queued','running','ready','error']}})
        if route=='/api/images':return self.json_response(image_capabilities())
        if route=='/api/document':
            with LOCK:return self.json_response(state)
        if route=='/api/player-state':
            with LOCK:return self.json_response(player_state())
        if route.startswith('/api/jobs/'):
            with JOB_LOCK:job=copy.copy(jobs.get(route.rsplit('/',1)[-1],{}))
            return self.json_response(job or {'error':'Unknown narration job.'},200 if job else 404)
        if route.startswith('/api/'):return self.json_response({'error':'Not found.'},404)
        return super().do_GET()
    def do_HEAD(self):
        if urlsplit(self.path).path.startswith('/api/'):return self.do_GET()
        return super().do_HEAD()
    def do_PUT(self):
        try:
            payload=self.body()
            if urlsplit(self.path).path!='/api/document':return self.json_response({'error':'Not found.'},404)
            validate_document(payload.get('document'))
            with LOCK:
                if payload.get('revision')!=state['revision']:return self.json_response({'error':'This document changed in another window. Reload before saving.','revision':state['revision']},409)
                atomic(HISTORY/f"revision-{state['revision']:06d}.json",state)
                state['document']=payload['document'];state['revision']+=1;atomic(STATE_PATH,state)
                for old in sorted(HISTORY.glob('revision-*.json'))[:-50]:old.unlink()
                return self.json_response({'revision':state['revision']})
        except PermissionError as e:return self.json_response({'error':str(e)},403)
        except (ValueError,TypeError) as e:return self.json_response({'error':str(e)},400)
        except Exception as e:return self.json_response({'error':str(e)},500)
    def do_POST(self):
        try:
            route=urlsplit(self.path).path
            payload=self.body(MAX_JSON_BYTES if route=='/api/images' else 8*1024*1024)
            if route=='/api/images':return self.json_response(store_image(ROOT,payload))
            if route=='/api/dictionary':return self.json_response(dictionary_lookup(payload.get('term'),RUNTIME.get('dictionary_path')))
            if route=='/api/selection':
                with LOCK:return self.json_response(selection(payload,state['document']))
            if route=='/api/prepare':
                key=prepare(payload.get('block'),100 if payload.get('background') else 0)
                with JOB_LOCK:return self.json_response({'jobId':key,**jobs[key]})
            if route=='/api/plan':return self.json_response(make_plan(payload['block']))
            if route=='/api/player-state':
                with LOCK:
                    owner=payload.get('owner')
                    if not isinstance(owner,str) or len(owner)>200:raise ValueError('Invalid player owner.')
                    if payload.get('claim') or player['owner']==owner:
                        t=float(payload.get('time',0))
                        if not math.isfinite(t) or t<0:raise ValueError('Invalid playback time.')
                        player.update(owner=owner,blockId=payload.get('blockId'),playing=bool(payload.get('playing')),time=t,updatedAt=time.time()*1000)
                    return self.json_response(player_state())
            return self.json_response({'error':'Not found.'},404)
        except PermissionError as e:return self.json_response({'error':str(e)},403)
        except (ValueError,TypeError,KeyError) as e:return self.json_response({'error':str(e)},400)
        except Exception as e:return self.json_response({'error':str(e)},500)
    def send_head(self):
        self.remaining=None
        if not self.allowed():self.send_error(403);return None
        route=unquote(urlsplit(self.path).path)
        if route=='/':self.send_response(302);self.send_header('Location','/app/index.html');self.send_header('Content-Length','0');self.end_headers();return None
        path=Path(self.translate_path(self.path)).resolve()
        if not any(path.is_relative_to(ROOT/d) for d in ['app','figures','audio']) or not path.is_file() or '/jobs/' in route or '/segments/' in route:
            self.send_error(404);return None
        requested=self.headers.get('Range')
        if not requested:return super().send_head()
        size=path.stat().st_size;match=re.fullmatch(r'bytes=(\d*)-(\d*)',requested)
        if match and any(match.groups()):
            left,right=match.groups()
            if left:start,end=int(left),min(int(right),size-1) if right else size-1
            else:start,end=max(0,size-int(right)),size-1
        else:start,end=size,size-1
        if not 0<=start<=end<size:self.send_response(416);self.send_header('Content-Range',f'bytes */{size}');self.send_header('Content-Length','0');self.end_headers();return None
        stream=path.open('rb');stream.seek(start);self.remaining=end-start+1
        self.send_response(206);self.send_header('Content-Type',self.guess_type(str(path)));self.send_header('Content-Range',f'bytes {start}-{end}/{size}');self.send_header('Content-Length',str(self.remaining));self.end_headers();return stream
    def end_headers(self):
        self.send_header('Accept-Ranges','bytes');self.send_header('X-Content-Type-Options','nosniff');self.send_header('Cache-Control','no-cache');super().end_headers()
    def copyfile(self,source,outputfile):
        if self.remaining is None:return super().copyfile(source,outputfile)
        while self.remaining:
            data=source.read(min(65536,self.remaining))
            if not data:break
            outputfile.write(data);self.remaining-=len(data)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--port',type=int,default=8770);args=parser.parse_args()
    threading.Thread(target=worker,daemon=True).start()
    server=ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    print(f'Notion Reader: http://127.0.0.1:{args.port}/app/index.html',flush=True)
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:
        server.server_close()
        for process in [planner,synthesizer]:
            if process and process.poll() is None:process.terminate()
