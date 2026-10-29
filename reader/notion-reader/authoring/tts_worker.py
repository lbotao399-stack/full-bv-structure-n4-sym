"""Persistent local Kokoro worker. Stdout is JSONL; captions use exact PCM boundaries."""
from pathlib import Path
import sys, json, hashlib, wave, os, contextlib, traceback, re
import numpy as np
import onnxruntime as ort
import kokoro_onnx

ROOT = Path(__file__).resolve().parents[1]
RUNTIME = json.loads((ROOT.parent / 'qa/runtime.json').read_text())
CACHE = ROOT / 'audio/segments'
CACHE.mkdir(parents=True, exist_ok=True)
CONFIG = dict(voice='af_heart', lang='en-us', speed=.95, sentence_pause=.34, clause_pause=.15)
ZH_CONFIG = dict(voice='zf_xiaobei', lang='cmn', speed=.95, sentence_pause=.34, clause_pause=.15)
IDENTITY = {**CONFIG, 'model':hashlib.sha256(Path(RUNTIME['model']).read_bytes()).hexdigest(), 'voices':hashlib.sha256(Path(RUNTIME['voices']).read_bytes()).hexdigest(), 'pcm':'mono-int16-gain-cap-0.92-v1'}
options=ort.SessionOptions(); options.intra_op_num_threads=2; options.inter_op_num_threads=1
kokoro_onnx.create_session=lambda p:ort.InferenceSession(p,sess_options=options,providers=['CPUExecutionProvider'])
model=None
g2p=None


def mandarin_phonemes(text):
    global g2p
    if g2p is None:
        deps=Path(RUNTIME.get('tts_deps',''))
        if not RUNTIME.get('tts_deps') or not deps.is_dir():raise RuntimeError('Set tts_deps to the installed Mandarin phonemizer directory in qa/runtime.json.')
        sys.path.insert(0,str(deps))
        from misaki.zh import ZHG2P
        g2p=ZHG2P()
    punctuation=str.maketrans({'，':',','。':'.','；':';','：':':','？':'?','！':'!','、':',','“':'"','”':'"','‘':"'",'’':"'",'（':'(', '）':')','—':',','–':',','−':'减','…':',','①':'一','②':'二','③':'三','④':'四','⑤':'五','⑥':'六','⑦':'七','⑧':'八','⑨':'九','§':'节'})
    text=text.translate(punctuation)
    parts=[]
    for part in re.split(r'([A-Za-z]+)',text):
        if not part:continue
        if re.fullmatch(r'[A-Za-z]+',part):
            spoken=' '.join(part) if part.isupper() and len(part)>1 else part
            parts.append(model.tokenizer.phonemize(spoken,lang='en-us'))
        else:parts.append(g2p(part))
    phonemes=' '.join(parts)
    unknown=set(phonemes)-set(model.tokenizer.vocab)
    if unknown:raise ValueError(f'Unsupported Mandarin speech symbols: {sorted(unknown)}')
    return phonemes


def write_wav(path, data, rate):
    temporary=path.with_suffix('.part.wav')
    with wave.open(str(temporary),'wb') as w:
        w.setnchannels(1); w.setsampwidth(2); w.setframerate(rate); w.writeframes(data.astype('<i2').tobytes())
    temporary.replace(path)

def synthesize(req):
    global model
    language='zh-CN' if req.get('language')=='zh-CN' else 'en-US'
    config=ZH_CONFIG if language=='zh-CN' else CONFIG
    identity={**IDENTITY,**config,'language':language,'phonemizer':'misaki-0.7.4-zh-english-letters-v1' if language=='zh-CN' else 'espeak-en-us'}
    parts=[]; cues=[]; offset=0; rate=24000
    for text in req['segments']:
        key=hashlib.sha256(json.dumps({'text':text,**identity},sort_keys=True).encode()).hexdigest()
        wav=CACHE/(key+'.wav')
        if wav.exists():
            with wave.open(str(wav),'rb') as w:
                rate=w.getframerate(); data=np.frombuffer(w.readframes(w.getnframes()),dtype='<i2').copy()
        else:
            with contextlib.redirect_stdout(sys.stderr):
                if model is None:model=kokoro_onnx.Kokoro(RUNTIME['model'],RUNTIME['voices'])
                if language=='zh-CN':samples,rate=model.create(mandarin_phonemes(text),**config,is_phonemes=True)
                else:samples,rate=model.create(text,**config)
            samples=np.asarray(samples)
            if samples.ndim!=1 or not len(samples) or not np.isfinite(samples).all():raise ValueError('Invalid synthesized samples')
            peak=float(np.abs(samples).max())
            if peak==0:raise ValueError('Silent synthesis')
            data=np.rint(samples*min(1.,.92/peak)*32767).astype('<i2')
            write_wav(wav,data,rate)
        start=offset/rate; parts.append(data); offset+=len(data)
        cues.append({'start':start,'end':offset/rate,'text':text})
        silence=np.zeros(round(rate*.34),dtype='<i2');parts.append(silence);offset+=len(silence)
    data=np.concatenate(parts)
    key=req['key']; dest=ROOT/'audio'/(key+'.wav');write_wav(dest,data,rate)
    return {'audio':'/audio/'+dest.name, 'segments':cues,'duration':len(data)/rate,'sampleRate':rate,'frames':len(data),'pcmPeak':int(np.abs(data.astype(np.int32)).max()),'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'voice':config,'language':language,'timing':'Exact synthesized clause PCM boundaries; 0.34 seconds of silence between segments','actualListening':'unverified'}

for line in sys.stdin:
    try: print(json.dumps({'ok':True,**synthesize(json.loads(line))}),flush=True)
    except Exception as e:
        traceback.print_exc(file=sys.stderr);print(json.dumps({'ok':False,'error':str(e)}),flush=True)
