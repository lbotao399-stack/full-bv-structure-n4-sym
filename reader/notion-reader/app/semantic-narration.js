/* Source-locked explanation tracks; no source edits or synthesis on typing. */
(function(root,factory){
 'use strict';const helpers=factory();
 if(typeof module==='object'&&module.exports)module.exports=helpers;
 if(root)root.ReaderSemanticNarration=helpers.create();
})(typeof window!=='undefined'?window:typeof globalThis!=='undefined'?globalThis:null,function(){
 'use strict';
 const clone=value=>JSON.parse(JSON.stringify(value));
 const validMode=value=>value==='exact'?'exact':'guided';
 const sourceKey=b=>JSON.stringify([b?.id,b?.type,b?.source||'',b?.speechOverride||'',...(b?.type==='image'?[b.caption||'',b.alt||'',b.image||null]:[])]);
 const imageSource=b=>typeof b?.image==='string'?b.image:b?.image?.src||b?.source||'';
 const imageAlt=b=>b?.image?.alt||b?.alt||'';
 const hashShape=value=>/^[a-f0-9]{64}$/.test(value||'');
 const validRect=rect=>Array.isArray(rect)&&rect.length===4&&rect.every(Number.isFinite)&&rect[0]>=0&&rect[1]>=0&&rect[2]>0&&rect[3]>0&&rect[0]+rect[2]<=1&&rect[1]+rect[3]<=1;
 function imageGeometry(box,naturalWidth,naturalHeight,fit='contain',position='50% 50%'){
  if(!(naturalWidth>0&&naturalHeight>0&&box.width>0&&box.height>0))return null;
  let width=box.width,height=box.height;
  const contain=Math.min(box.width/naturalWidth,box.height/naturalHeight),cover=Math.max(box.width/naturalWidth,box.height/naturalHeight);
  if(fit!=='fill'){const scale=fit==='cover'?cover:fit==='none'?1:fit==='scale-down'?Math.min(1,contain):contain;width=naturalWidth*scale;height=naturalHeight*scale;}
  let values=String(position||'50% 50%').trim().split(/\s+/);if(values.length===1)values.push('50%');if(values.length!==2)return null;
  if(/^(top|bottom)$/.test(values[0]))values=[values[1],values[0]];
  function offset(value,free){if(value==='left'||value==='top')return 0;if(value==='right'||value==='bottom')return free;if(value==='center')return free/2;if(/^-?[\d.]+%$/.test(value))return free*parseFloat(value)/100;if(/^-?[\d.]+(?:px)?$/.test(value))return parseFloat(value);return free/2;}
  return {left:box.left+offset(values[0],box.width-width),top:box.top+offset(values[1],box.height-height),width,height,clip:box};
 }
 function targetGeometry(geometry,rect){
  if(!geometry||!validRect(rect))return null;const raw={left:geometry.left+rect[0]*geometry.width,top:geometry.top+rect[1]*geometry.height,width:rect[2]*geometry.width,height:rect[3]*geometry.height},clip=geometry.clip;
  const left=Math.max(raw.left,clip.left),top=Math.max(raw.top,clip.top),right=Math.min(raw.left+raw.width,clip.left+clip.width),bottom=Math.min(raw.top+raw.height,clip.top+clip.height);
  return right>left&&bottom>top?{left,top,width:right-left,height:bottom-top}:null;
 }
 function normalizeSpoken(value,equationLabels=[]){
  const speechText=typeof module==='object'&&module.exports?require('./speech-text.js'):window.ReaderSpeechText;
  const references=typeof module==='object'&&module.exports?require('./speech-references.js'):window.ReaderSpeechReferences;
  return references.normalize(speechText.normalize(value),{equationLabels})
   .replace(/[ \t]+/g,' ').replace(/\s*\n\s*/g,'\n').replace(/;\s*\./g,'.').replace(/\s+([,.!?;:])/g,'$1').replace(/\.{2,}/g,'.').replace(/\s+/g,' ').trim();
 }
 function mapMeasured(metadata,record,equationLabels=[]){
  if(!Array.isArray(metadata?.segments)||!Array.isArray(record?.segments))return [];
  const measured=record.segments;let cursor=0;const result=[];
  for(const segment of metadata.segments){
   const expected=normalizeSpoken(segment.text,equationLabels);if(!expected)return [];
   const from=cursor;let accumulated='';
   while(cursor<measured.length){
    const next=normalizeSpoken(measured[cursor]?.text,equationLabels);if(!next)return [];
    accumulated+=(accumulated?' ':'')+next;cursor++;
    if(accumulated===expected)break;
    if(!expected.startsWith(accumulated))return [];
   }
   if(accumulated!==expected)return [];
   result.push({from,to:cursor,targets:segment.targets||[]});
  }
  return cursor===measured.length?result:[];
 }
 function create(environment={}){
  let api={},control=null,loaded=null,loadedAt=0,loadedOK=false,lastCue=null,highlighted=[],focusNode=null,imageState=null;
  const modes=new Map(),hashes=new Map(),mapped=new WeakMap(),imageRequests=new Map();
  const dom=()=>environment.document||document;
  const storage=()=>environment.storage||localStorage;
  const network=(...args)=>(environment.fetch||fetch)(...args);
  const now=()=>environment.now?environment.now():Date.now();
  const abort=signal=>{if(signal?.aborted)throw new DOMException('Narration preparation cancelled.','AbortError');};
  function documentID(){return api.getDocument?.()?.id||null;}
  function mode(){
   const id=documentID();if(!id)return 'guided';
   if(!modes.has(id)){let value;try{value=storage().getItem(`notion-reader:${id}:narration-mode`);}catch(_e){}modes.set(id,validMode(value));}
   return modes.get(id);
  }
  function cacheTag(block){return block?._readerMode==='exact'||block?._readerMode==='guided'?block._readerMode:mode();}
  function refreshControl(){if(control)control.value=mode();}
  function setMode(value){
   const next=validMode(value),id=documentID();if(!id)return;
   modes.set(id,next);refreshControl();
   try{storage().setItem(`notion-reader:${id}:narration-mode`,next);}catch(_e){api.toast?.('Narration style works in this tab, but browser storage is unavailable.');}
   api.onModeChange?.(next);
  }
  async function sha256(text){
   text=String(text);if(hashes.has(text))return hashes.get(text);
   const promise=(async()=>{const cryptoAPI=environment.crypto||crypto;const digest=await cryptoAPI.subtle.digest('SHA-256',new TextEncoder().encode(text));return [...new Uint8Array(digest)].map(n=>n.toString(16).padStart(2,'0')).join('');})();
   hashes.set(text,promise);if(hashes.size>128)hashes.delete(hashes.keys().next().value);
   try{return await promise;}catch(error){hashes.delete(text);throw error;}
  }
  function validPlan(plan){
   const type=plan?.type;
   return plan?.reviewed===true&&typeof plan.blockId==='string'&&['table','equation','image','paragraph','heading','callout'].includes(type)&&hashShape(plan.sourceSha256)&&hashShape(plan.originalSpeechOverrideSha256)&&['en-US','zh-CN'].includes(plan.language)&&Array.isArray(plan.segments)&&plan.segments.length>0&&plan.segments.every(segment=>typeof segment?.text==='string'&&segment.text.trim()&&Array.isArray(segment.targets)&&segment.targets.every(target=>target&&typeof target==='object'&&(type==='table'?Number.isInteger(target.row)&&target.row>=0&&Number.isInteger(target.col)&&target.col>=0:type==='image'?validRect(target.rect)&&typeof target.label==='string'&&target.label.trim():type==='equation'&&typeof target.tex==='string'&&target.tex.trim())))&&Array.isArray(plan.context||[])&&(plan.context||[]).every(ref=>typeof ref?.blockId==='string'&&hashShape(ref.sourceSha256))&&(type!=='image'||hashShape(plan.imageSha256)&&hashShape(plan.captionSha256)&&hashShape(plan.altSha256));
  }
  async function verifyImage(b,plan,signal){
   if(await sha256(b.caption||'')!==plan.captionSha256||await sha256(imageAlt(b))!==plan.altSha256)return null;
   abort(signal);const src=imageSource(b);if(!src)return null;const url=api.asset?.(src)||new URL(src,dom().baseURI).href,key=url+'|'+plan.imageSha256;
   let pending=imageRequests.get(key);
   if(!pending){pending=(async()=>{const response=await network(url,{cache:'no-cache'});if(!response.ok)throw Error('Unavailable source figure.');const data=await response.arrayBuffer(),cryptoAPI=environment.crypto||crypto;const digest=await cryptoAPI.subtle.digest('SHA-256',data);const hash=[...new Uint8Array(digest)].map(n=>n.toString(16).padStart(2,'0')).join('');return hash===plan.imageSha256?{url,sha256:hash}:null;})();imageRequests.set(key,pending);pending.then(()=>{if(imageRequests.get(key)===pending)imageRequests.delete(key);},()=>{if(imageRequests.get(key)===pending)imageRequests.delete(key);});}
   return await waitFor(pending,signal);
  }
  async function manifest(){
   // A failed optional sidecar is retried later, without turning background playback into a request loop.
   if(loaded&&(loadedOK||now()-loadedAt<60000))return loaded;
   loadedAt=now();loaded=(async()=>{try{
    const url=api.manifestURL||api.asset?.('app/narration-plans.json')||new URL('narration-plans.json',dom().baseURI).href;
    const response=await network(url,{cache:'no-store'});if(!response.ok)throw Error('Unavailable narration plans.');
    const payload=await response.json();if(payload?.version!==1||!Array.isArray(payload.plans))throw Error('Unsupported narration plans.');
    loadedOK=true;return payload.plans.filter(validPlan);
   }catch(_error){loadedOK=false;return [];}})();return loaded;
  }
  function waitFor(value,signal){
   if(!signal)return value;abort(signal);
   return new Promise((resolve,reject)=>{const cancelled=()=>{signal.removeEventListener('abort',cancelled);reject(new DOMException('Narration preparation cancelled.','AbortError'));};signal.addEventListener('abort',cancelled,{once:true});Promise.resolve(value).then(result=>{signal.removeEventListener('abort',cancelled);resolve(result);},error=>{signal.removeEventListener('abort',cancelled);reject(error);});});
  }
  async function resolve(original,{signal}={}){
   abort(signal);const selected=cacheTag(original),b=clone(original);b._readerMode=selected;delete b._readerNarration;
   if(!['table','equation','image','paragraph','heading','callout'].includes(b.type))return b;
   if(selected==='exact'){delete b.speechOverride;return b;}
   try{
    const plans=await waitFor(manifest(),signal);abort(signal);
    const candidates=plans.filter(p=>p.blockId===b.id&&p.type===b.type);if(!candidates.length)return b;
    const sourceHash=await sha256(b.source||''),overrideHash=await sha256(b.speechOverride||'');abort(signal);
    for(const plan of candidates){
     if(plan.sourceSha256!==sourceHash||plan.originalSpeechOverrideSha256!==overrideHash)continue;
     const imageVerified=b.type==='image'?await verifyImage(b,plan,signal):null;if(b.type==='image'&&!imageVerified)continue;abort(signal);
     const currentDoc=api.getDocument?.(),current=currentDoc?.blocks?.find(item=>item.id===b.id);
     if(!current||sourceKey(current)!==sourceKey(b))continue;
     const context=[];let matches=true;
     for(const ref of plan.context||[]){
      const item=currentDoc.blocks.find(candidate=>candidate.id===ref.blockId);
      if(!item||await sha256(item.source||'')!==ref.sourceSha256){matches=false;break;}
      context.push({blockId:item.id,source:item.source||''});abort(signal);
     }
     if(!matches)continue;
     const latest=api.getDocument?.();
     if(sourceKey(latest?.blocks?.find(item=>item.id===b.id))!==sourceKey(b)||context.some(ref=>latest?.blocks?.find(item=>item.id===ref.blockId)?.source!==ref.source))continue;
     b.speechOverride=plan.segments.map(item=>item.text.trim()).join('\n');
     b._readerNarration={version:1,mode:'guided',language:plan.language,blockId:b.id,sourceKey:sourceKey(original),context,segments:clone(plan.segments),sourceSha256:plan.sourceSha256,...(imageVerified?{imageVerified}:{})};
     return b;
    }
   }catch(error){if(error?.name==='AbortError')throw error;/* Optional explanation failure keeps existing source reading usable. */}
   return b;
  }
  function clear(){
   if(imageState){imageState.observer?.disconnect();imageState.image.removeEventListener?.('load',imageState.update);(environment.window||globalThis).removeEventListener?.('resize',imageState.update);imageState.overlay.remove();imageState=null;}
   for(const node of highlighted)node.classList.remove('semantic-reading-cue','semantic-equation-cue');highlighted=[];
   focusNode?.remove();focusNode=null;lastCue=null;
  }
  function contextCurrent(metadata){
   const blocks=api.getDocument?.()?.blocks||[];
   return (metadata.context||[]).every(ref=>blocks.find(item=>item.id===ref.blockId)?.source===ref.source);
  }
  function showImageCue(block,metadata,cue,wrapper){
   const image=wrapper.querySelector('.block-content img'),targets=cue.targets.filter(target=>validRect(target?.rect));
   if(!image||!targets.length||!metadata.imageVerified)return;
   const shown=image.currentSrc||image.src;
   if(shown&&new URL(shown,dom().baseURI).href!==metadata.imageVerified.url)return;
   const overlay=dom().createElement('div');overlay.className='semantic-image-overlay';overlay.setAttribute('aria-hidden','true');wrapper.append(overlay);
   focusNode=dom().createElement('aside');focusNode.className='narration-focus narration-image-focus';focusNode.setAttribute('aria-label',metadata.language==='zh-CN'?'正在讲解的图片区域':'Current figure panels');wrapper.append(focusNode);
   function update(){
    if(imageState?.overlay!==overlay)return;
    const width=image.naturalWidth,height=image.naturalHeight;if(!(width>0&&height>0)){overlay.hidden=true;focusNode.hidden=true;return;}
    const imageBox=image.getBoundingClientRect(),wrapBox=wrapper.getBoundingClientRect(),style=(environment.getComputedStyle||(typeof getComputedStyle==='function'?getComputedStyle:()=>({})))(image),px=key=>parseFloat(style[key])||0;
    const borderLeft=px('borderLeftWidth'),borderRight=px('borderRightWidth'),borderTop=px('borderTopWidth'),borderBottom=px('borderBottomWidth'),paddingLeft=px('paddingLeft'),paddingRight=px('paddingRight'),paddingTop=px('paddingTop'),paddingBottom=px('paddingBottom');
    const box={left:imageBox.left-wrapBox.left-(wrapper.clientLeft||0)+borderLeft+paddingLeft,top:imageBox.top-wrapBox.top-(wrapper.clientTop||0)+borderTop+paddingTop,width:imageBox.width-borderLeft-borderRight-paddingLeft-paddingRight,height:imageBox.height-borderTop-borderBottom-paddingTop-paddingBottom};
    const geometry=imageGeometry(box,width,height,style.objectFit||'contain',style.objectPosition||'50% 50%');
    overlay.replaceChildren();overlay.hidden=!geometry;
    for(const target of targets){
     const bounds=targetGeometry(geometry,target.rect);if(!bounds)continue;
     const rectangle=dom().createElement('div');rectangle.className='semantic-image-cue';rectangle.title=target.label||'';
     for(const key of ['left','top','width','height'])rectangle.style[key]=bounds[key]+'px';overlay.append(rectangle);
    }
    focusNode.replaceChildren();focusNode.hidden=false;
    const label=dom().createElement('span');label.className='narration-focus-label';label.textContent=metadata.language==='zh-CN'?'当前区域':'Current panel';focusNode.append(label);
    const panels=dom().createElement('div');panels.className='narration-image-panels';focusNode.append(panels);
    const ns='http://www.w3.org/2000/svg';
    for(const target of targets){
     const [x,y,w,h]=target.rect,cropWidth=w*width,cropHeight=h*height,displayWidth=Math.min(400,320*cropWidth/cropHeight),displayHeight=displayWidth*cropHeight/cropWidth;
     const figure=dom().createElement('figure');figure.className='narration-image-panel';figure.style.maxWidth=displayWidth+'px';
     const svg=dom().createElementNS(ns,'svg');svg.setAttribute('viewBox',`${x*width} ${y*height} ${cropWidth} ${cropHeight}`);svg.setAttribute('width',String(displayWidth));svg.setAttribute('height',String(displayHeight));svg.setAttribute('role','img');svg.setAttribute('aria-label',target.label||'Current panel');svg.setAttribute('preserveAspectRatio','xMidYMid meet');
     const pixels=dom().createElementNS(ns,'image');pixels.setAttribute('href',metadata.imageVerified.url);pixels.setAttribute('width',String(width));pixels.setAttribute('height',String(height));pixels.setAttribute('preserveAspectRatio','none');svg.append(pixels);figure.append(svg);
     const caption=dom().createElement('figcaption');caption.textContent=target.label||'';figure.append(caption);panels.append(figure);
    }
   }
   imageState={image,overlay,update,observer:null};
   const Observer=environment.ResizeObserver||(typeof ResizeObserver==='function'?ResizeObserver:null);
   if(Observer){imageState.observer=new Observer(update);imageState.observer.observe(image);imageState.observer.observe(wrapper);}
   else (environment.window||globalThis).addEventListener?.('resize',update);
   image.addEventListener?.('load',update);update();
  }
  function isCurrent(block,record){
   const metadata=record?.readerNarration;
   return !metadata||metadata.blockId===block?.id&&metadata.sourceKey===sourceKey(block)&&contextCurrent(metadata);
  }
  function sync(block,record,index,wrapper){
   const metadata=record?.readerNarration;
   if(!wrapper||!metadata||metadata.blockId!==block?.id||metadata.sourceKey!==sourceKey(block)||!contextCurrent(metadata)){clear();return;}
   let cues=mapped.get(record);if(!cues){cues=mapMeasured(metadata,record,(typeof module==='object'&&module.exports?require('./speech-references.js'):window.ReaderSpeechReferences).labels(api.getDocument?.()));mapped.set(record,cues);}
   const cue=cues.find(item=>index>=item.from&&index<item.to);
   if(!cue){clear();return;}
   if(lastCue?.record===record&&lastCue?.cue===cue&&lastCue?.wrapper===wrapper)return;
   clear();lastCue={record,cue,wrapper};
   if(block.type==='table'){
    const table=wrapper.querySelector('table');if(!table)return;
    for(const target of cue.targets){
     if(!target||!Number.isInteger(target.row)||target.row<0||!Number.isInteger(target.col)||target.col<0)continue;
     const cell=table.rows?.[target.row]?.cells?.[target.col];if(cell&&!highlighted.includes(cell)){cell.classList.add('semantic-reading-cue');highlighted.push(cell);}
    }
    return;
   }
   if(block.type==='image'){showImageCue(block,metadata,cue,wrapper);return;}
   const source=block.source||'',terms=[...new Set(cue.targets.filter(t=>typeof t?.tex==='string'&&t.tex.trim()&&source.includes(t.tex)).map(t=>t.tex))];
   if(!terms.length)return;
   const equation=wrapper.querySelector('.math-display')||wrapper.querySelector('.katex-display')||wrapper.querySelector('.block-content');
   if(equation){equation.classList.add('semantic-equation-cue');highlighted.push(equation);}
   if(typeof api.renderMath!=='function')return;
   // This supplement is outside editable source: it identifies the exact term without claiming token coordinates.
   focusNode=dom().createElement('aside');focusNode.className='narration-focus';focusNode.setAttribute('aria-label',metadata.language==='zh-CN'?'正在讲解的公式部分':'Current formula terms');
   const label=dom().createElement('span');label.className='narration-focus-label';label.textContent=metadata.language==='zh-CN'?'当前项':'Current term';focusNode.append(label);
   for(const tex of terms){const term=dom().createElement('div');term.className='narration-focus-term';term.innerHTML=api.renderMath(tex,true);focusNode.append(term);}
   wrapper.append(focusNode);
  }
  function init(options){
   api=options||{};if(control){refreshControl();return;}
   const speed=dom().querySelector('#speed');if(!speed)return;
   control=dom().createElement('select');control.id='narration-mode';control.className='narration-mode';control.setAttribute('aria-label','Narration style');control.title='Guided explanations where reviewed; Exact reads displayed tables, formulas, or figure captions.';
   for(const [value,text]of [['guided','Guided'],['exact','Exact']]){const option=dom().createElement('option');option.value=value;option.textContent=text;control.append(option);}
   control.onchange=()=>setMode(control.value);speed.parentNode.insertBefore(control,speed);refreshControl();
  }
  return {init,resolve,sync,clear,mode,cacheTag,refreshControl,setMode,isCurrent};
 }
 return {create,normalizeSpoken,mapMeasured,imageGeometry,targetGeometry};
});
