/* Local source manuscript editor and measured-audio reader. */
(function(){'use strict';
const $=s=>document.querySelector(s),$$=(s,p=document)=>Array.from(p.querySelectorAll(s));
const ROOT=new URL('../',document.baseURI),floating=new URLSearchParams(location.search).get('floating')==='1';
const owner='reader-'+(crypto.randomUUID?.()||Math.random().toString(36).slice(2));
const audio=$('#audio'),cache=new Map(),historyStack=[],redoStack=[];
const md=window.markdownit?.({html:true,breaks:false,linkify:false,typographer:false});
let doc=null,revision=0,selectedId=null,activeId=null,editorId=null,busy=false,playIntent=false,continuous=true;
let activeRecording=null,figureCueNode=null;
let serial=0,controller=null,segments=[],lastCaption=-1,raf=0,saveTimer=0,saving=false,pendingSave=false,dirty=false;
let follow=false,lastUndo=0,lastHeartbeat=0,remoteOwner=null,playerAPI=true,loadingDocument=false,toastTimer=0,blockNodes=new Map();
let controlsOnly=floating&&new URLSearchParams(location.search).get('controls')==='1',temporaryBlock=null,selectionQueue=null,selectionRequest=0,documentReady=null;
const channel=typeof BroadcastChannel!=='undefined'?new BroadcastChannel('notion-manuscript-reader'):null;
if(floating)document.body.classList.add('floating');
document.body.classList.toggle('controls-only',controlsOnly);
const clone=x=>JSON.parse(JSON.stringify(x));
const asset=s=>new URL(s||'',ROOT).href;
const block=id=>temporaryBlock?.id===id?temporaryBlock:doc?.blocks.find(b=>b.id===id);
const index=id=>doc?.blocks.findIndex(b=>b.id===id)??-1;
const fingerprint=b=>JSON.stringify([b?.id,b?.type,b?.level,b?.source,b?.image,b?.alt,b?.caption,b?.speechOverride]);
const narrationFingerprint=b=>fingerprint(b)+'|'+(window.ReaderSemanticNarration?.cacheTag(b)||'guided')+'|'+(window.ReaderVoiceTransport?.cacheTag()||'legacy');
const narrationCurrent=(b,record)=>window.ReaderSemanticNarration?.isCurrent(b,record)!==false;
const renderedBlocks=new Map();
let draftTimer=0,draftDeadline=0,draftPending=false,outlineTimer=0,outlineDirty=true;
let pages=[],pageIndex=0,pageAnchor=null,pageDocumentId=null,pageChanging=false,pageOptionsKey='',printing=false,printPosition=null;
const renderKey=b=>JSON.stringify([b.type,b.level,b.blockStyle,b.source,b.image,b.alt,b.caption,b.formatting,b.sourceBlockUrl,Boolean(b.figureCues),b.syncGroup]);
const storeKey=()=>`notion-reader:${doc?.id||'document'}`;
// Small preferences save independently of manuscript edits and narration.
const playbackSpeeds=[0.5,0.6,0.7,0.75,0.8,0.85,0.9,1,1.1,1.15,1.2,1.25,1.3,1.4,1.5,1.6,1.75,2,2.25,2.5,2.75,3];
function readSettings(){
 let value;try{value=JSON.parse(localStorage.getItem(storeKey()+':settings')||'null');}catch(_e){}
 return {speed:playbackSpeeds.includes(value?.speed)?value.speed:1,follow:value?.follow===true,outline:value?.outline===true};
}
function setPlaybackSpeed(value){
 const rate=playbackSpeeds.includes(Number(value))?Number(value):1;
 $('#speed').value=String(rate);enforcePlaybackSpeed();
}
function enforcePlaybackSpeed(){
 const selected=Number($('#speed').value),rate=playbackSpeeds.includes(selected)?selected:1;
 if(audio.defaultPlaybackRate!==rate)audio.defaultPlaybackRate=rate;
 if(audio.playbackRate!==rate)audio.playbackRate=rate;
}
function setFollow(enabled,{reveal=false}={}){
 follow=Boolean(enabled)&&!controlsOnly;$('#follow').setAttribute('aria-pressed',String(follow));$('#follow').textContent=follow?'Follow on':'Follow off';
 if(follow&&reveal)revealBlock(activeId||selectedId);
}
function setOutline(open){
 open=Boolean(open)&&!floating;$('#outline').hidden=!open;document.body.classList.toggle('outline-open',open);$('#outline-toggle').setAttribute('aria-expanded',String(open));
 if(open&&doc&&outlineDirty)drawOutline();
}
function restoreSettings(){const settings=readSettings();setPlaybackSpeed(settings.speed);setFollow(settings.follow);setOutline(settings.outline);}
function saveSettings(patch){
 if(!doc)return;
 try{localStorage.setItem(storeKey()+':settings',JSON.stringify({...readSettings(),...patch}));}
 catch(_e){toast('Settings work in this tab, but browser storage is unavailable.');}
}
// Narration ignores empty layout blocks; editing and their positions are preserved.
function hasReadableText(value){
 let text=String(value||'').replace(/<!--[\s\S]*?-->/g,'').replace(/<(script|style)\b[^>]*>[\s\S]*?<\/\1\s*>/gi,'').replace(/<[^>]*>/g,'');
 text=text.replace(/&#(?:x([\da-f]+)|(\d+));/gi,(_m,hex,decimal)=>{const n=parseInt(hex||decimal,hex?16:10);return n>=0&&n<=0x10ffff?String.fromCodePoint(n):'';});
 text=text.replace(/&(?:nbsp|ensp|emsp|thinsp|ZeroWidthSpace|zwnj|zwj);/gi,' ');
 text=text.replace(/[\s\u200b-\u200f\u202a-\u202e\u2060-\u206f\ufeff]/g,'').replace(/\\[()[\]]/g,'').replace(/[$*_`~#]/g,'');
 return Boolean(text)&&!/^[-|:]+$/.test(text);
}
function isNarratable(b){
 if(!b)return false;
 if(b.type==='image')return Boolean(typeof b.image==='string'?b.image.trim():b.image?.src?.trim()||b.source?.trim())||hasReadableText(b.caption||b.alt||b.image?.alt);
 return hasReadableText(b.source);
}
function readableNeighbor(id,delta=1,includeCurrent=false){
 const ids=selectionQueue||doc?.blocks.map(b=>b.id)||[],start=ids.indexOf(id);
 if(start<0)return null;
 for(let i=start+(includeCurrent?0:delta);i>=0&&i<ids.length;i+=delta)if(isNarratable(block(ids[i])))return ids[i];
 return null;
}
// End empty-block narration helpers.
function time(n){n=Number.isFinite(n)?Math.max(0,Math.floor(n)):0;return `${Math.floor(n/60)}:${String(n%60).padStart(2,'0')}`;}
function status(s){$('#play-status').textContent=s;}
function toast(s){$('#toast').textContent=s;$('#toast').classList.add('visible');clearTimeout(toastTimer);toastTimer=setTimeout(()=>$('#toast').classList.remove('visible'),4000);}
function saved(s,error=false){$('#save-status').textContent=s;$('#save-status').classList.toggle('error',error);}
function node(tag,cls,text){const n=document.createElement(tag);if(cls)n.className=cls;if(text!==undefined)n.textContent=text;return n;}
function escape(s){return String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));}
function clean(s){return window.DOMPurify?DOMPurify.sanitize(s,{ADD_ATTR:['target','rel'],ADD_TAGS:['annotation','semantics']}):escape(s);}
function renderMath(s,display){try{return katex.renderToString(String(s),{displayMode:display,throwOnError:true,strict:'ignore',trust:false,output:'htmlAndMathml'});}catch(e){return '<span class="math-error" title="'+escape(e.message)+'">'+escape(s)+'</span>';}}
function markup(raw){
 let source=String(raw||''),tokens=[];
 const protect=(tex,display)=>{const key=`MATHPLACEHOLDER${tokens.length}END`;tokens.push([key,`<${display?'div':'span'} class="${display?'math-display':'math-inline'}">${renderMath(tex,display)}</${display?'div':'span'}>`]);return display?'\n\n'+key+'\n\n':key;};
 source=source.replace(/\$\$([\s\S]*?)\$\$/g,(_m,t)=>protect(t,true));
 source=source.replace(/\\\[([\s\S]*?)\\\]/g,(_m,t)=>protect(t,true));
 source=source.replace(/\$`([\s\S]*?)`\$/g,(_m,t)=>protect(t,false));
 source=source.replace(/\\\(([\s\S]*?)\\\)/g,(_m,t)=>protect(t,false));
 source=source.replace(/(?<!\\)\$([^$\n]+?)\$/g,(_m,t)=>protect(t,false));
 source=source.replace(/<span color="yellow_bg">/g,'<span class="source-highlight">');
 let result=md?md.render(source):'<p>'+escape(source).replace(/\n/g,'<br>')+'</p>';
 for(const [key,value]of tokens)result=result.replaceAll('<p>'+key+'</p>',value).replaceAll(key,value);
 const box=document.createElement('template');box.innerHTML=clean(result);
 $$('a',box.content).forEach(a=>{if((a.getAttribute('href')||'').startsWith('#')){a.removeAttribute('target');a.removeAttribute('rel');}else{a.target='_blank';a.rel='noopener';}});
 $$('img',box.content).forEach(im=>{const s=im.getAttribute('src');if(s)im.src=asset(s);});
 return box.innerHTML;
}
function headingMarkup(source){
 const holder=document.createElement('div');holder.innerHTML=markup(source);
 // Preserve every paragraph and formula when a multi-line block becomes a heading.
 for(const el of [...holder.querySelectorAll('p,div,h1,h2,h3,h4,h5,h6,blockquote,pre,ul,ol,li')].reverse()){
  const span=document.createElement('span');for(const a of el.attributes)span.setAttribute(a.name,a.value);
  span.append(...el.childNodes);const end=document.createElement('br');el.replaceWith(span,end);
 }
 while(holder.lastChild?.nodeName==='BR')holder.lastChild.remove();return holder.innerHTML;
}
function contentFor(b){
 if(b.type==='equation')return '<div class="math-display">'+renderMath(b.source||'',true)+'</div>';
 if(b.type==='image'){
  const src=typeof b.image==='string'?b.image:(b.image?.src||b.source),alt=b.image?.alt||b.alt||'';
  const image=`<img src="${escape(asset(src))}" alt="${escape(alt)}">`;
  const visual=b.figureCues?`<div class="figure-stage">${image}<div class="figure-cue" hidden aria-hidden="true"></div></div>`:image;
  return `<figure>${visual}<figcaption>${markup(b.caption||'')}</figcaption></figure>`;
 }
 if(b.type==='heading'){const level=Math.min(3,Math.max(1,Number(b.level)||2));return `<h${level}>${headingMarkup(b.source||'')}</h${level}>`;}
 return markup(b.source||'')||'<p class="empty-block">Empty passage · double-click to write</p>';
}
function blockLabel(b){if(!b)return '';if(b.type==='image')return String(b.caption||'Figure').replace(/<[^>]*>/g,'').slice(0,100);return String(b.source||b.type).replace(/<[^>]*>/g,'').replace(/\s+/g,' ').slice(0,110);}
function renderBlock(b){
 const wrap=node('section',`source-block type-${b.type}`);wrap.classList.toggle('block-quote',b.type==='paragraph'&&b.blockStyle==='quote');wrap.dataset.blockId=b.id;if(b.sentenceSplit?.originId)wrap.dataset.paragraph=b.sentenceSplit.originId;wrap.tabIndex=0;wrap.setAttribute('aria-label',`${b.type}: ${blockLabel(b)}`);
 const cursor=node('span','reading-cursor');cursor.setAttribute('aria-hidden','true');wrap.append(cursor);
 const content=node('div','block-content');content.innerHTML=contentFor(b);window.ReaderFormatting?.decorate(b,content);wrap.append(content);if(!floating){window.ReaderDirectEdit?.bind(b,content);window.ReaderBlocks?.bind(b,wrap);window.ReaderBlockSelection?.bind(b,wrap);window.ReaderImages?.bind(b,wrap);}
 wrap.onclick=e=>{if(e.target.closest('a,button,input,textarea,select,summary'))return;selectBlock(b.id,{reveal:false});};
 wrap.ondblclick=e=>{if(e.target.closest('a,button,[contenteditable]')||floating)return;if(!window.getSelection()?.toString())openEditor(b.id);};
 wrap.onkeydown=e=>{if(e.target!==wrap)return;if(e.key==='Enter'){e.preventDefault();floating?playFrom(b.id,true):openEditor(b.id);}};
 return wrap;
}
// Pagination changes the rendered view, never the full manuscript model.
function rebuildPages(preferred=pageAnchor){
 pages=window.ReaderPagination.partition(doc.blocks,{chapterLevel:doc.pagination?.chapterLevel,documentTitle:doc.title,frontMatterLabel:doc.pagination?.frontMatterLabel});
 let next=window.ReaderPagination.pageFor(pages,doc.blocks,preferred);
 if(next<0)next=Math.min(pageIndex,Math.max(0,pages.length-1));
 pageIndex=next;pageAnchor=pages[pageIndex]?.anchor||null;pageDocumentId=doc.id;
}
function currentPageContains(id){const p=pages[pageIndex],at=index(id);return Boolean(p&&at>=p.start&&at<p.end);}
function persistPage(){
 if(!doc||printing)return;
 try{localStorage.setItem(storeKey()+':page',JSON.stringify({anchor:pageAnchor}));}catch(_e){}
}
function plainPageLabel(source){const holder=node('span');holder.innerHTML=headingMarkup(source);for(const mathml of holder.querySelectorAll('.katex-mathml'))mathml.remove();return holder.textContent.replace(/\s+/g,' ').trim();}
function refreshPageControls(){
 const signature=JSON.stringify(pages.map(p=>[p.anchor,p.label]));
 if(signature!==pageOptionsKey){
  pageOptionsKey=signature;
  for(const select of $$('[data-page-select]')){
   const options=document.createDocumentFragment();
   pages.forEach((page,i)=>{const option=node('option');option.value=String(i);const full=plainPageLabel(page.label),label=full.length>90?full.slice(0,87)+'…':full;option.textContent=label;options.append(option);});select.replaceChildren(options);
  }
 }
 for(const select of $$('[data-page-select]'))select.value=String(pageIndex);
 for(const button of $$('[data-page-prev]'))button.disabled=pageIndex<=0;
 for(const button of $$('[data-page-next]'))button.disabled=pageIndex>=pages.length-1;
 for(const label of $$('[data-page-status]'))label.textContent=`${pageIndex+1} / ${pages.length}`;
 for(const nav of $$('.reader-pagination'))nav.hidden=printing||controlsOnly||pages.length<=1;
 const heading=$('.paper-heading');if(heading)heading.hidden=!printing&&pageIndex!==0;
}
function visibleBlocks(){if(printing)return doc.blocks;const p=pages[pageIndex];return [...(temporaryBlock?[temporaryBlock]:[]),...(p?doc.blocks.slice(p.start,p.end):[])];}
function renderPage(){
 window.ReaderInsertion?.clear();window.ReaderSemanticNarration?.clear();clearFigureCue();
 const holder=$('#blocks'),blocks=visibleBlocks(),ids=new Set(blocks.map(b=>b.id));
 for(const [id,n]of blockNodes){if(!ids.has(id)){n.remove();blockNodes.delete(id);renderedBlocks.delete(id);}}
 let next=holder.firstChild;
 for(const b of blocks){
  let n=blockNodes.get(b.id);const key=renderKey(b);
  if(!n||renderedBlocks.get(b.id)!==key){const fresh=renderBlock(b);if(n){if(next===n)next=fresh;n.replaceWith(fresh);}n=fresh;blockNodes.set(b.id,n);renderedBlocks.set(b.id,key);}
  if(n!==next)holder.insertBefore(n,next);next=n.nextSibling;window.ReaderBlocks?.updateLabel(b,n,index(b.id));
 }
 window.ReaderBlocks?.refreshStructure();window.ReaderInsertion?.refresh();window.ReaderBlockSelection?.refresh();refreshPageControls();updateSelection();
 if(activeRecording)syncAudio();persistPage();
}
function preparePageChange(){
 if(window.ReaderInsertion?.preparePageChange?.()===false)return false;
 if(window.ReaderDirectEdit?.preparePageChange){if(window.ReaderDirectEdit.preparePageChange()===false)return false;}
 else window.ReaderDirectEdit?.flush();
 flushDraft();if(pendingSave)saveDocument();return true;
}
function showPage(next,{scroll=false,updateURL=false}={}){
 if(!doc||printing)return false;
 next=Math.max(0,Math.min(pages.length-1,Number(next)||0));
 if(next===pageIndex)return true;
 if(pageChanging)return false;
 pageChanging=true;
 try{
  if(!preparePageChange()){refreshPageControls();return false;}
  pageIndex=next;pageAnchor=pages[pageIndex]?.anchor||null;
  if(updateURL&&!playIntent&&!busy&&audio.paused)selectedId=pageAnchor;
  renderPage();
  if(updateURL&&pageAnchor)try{history.replaceState(null,'','#'+encodeURIComponent(pageAnchor));}catch(_e){}
  if(scroll){const first=blockNodes.get(pageAnchor),top=first?.getBoundingClientRect().top;const inset=($('.topbar')?.getBoundingClientRect().height||0)+($('.page-navigation-top')?.getBoundingClientRect().height||0)+14;if(Number.isFinite(top))window.scrollTo({top:Math.max(0,scrollY+top-inset),behavior:'instant'});}
  return true;
 }finally{pageChanging=false;}
}
function ensureBlockPage(id){
 if(!doc||temporaryBlock?.id===id||printing)return true;
 // Structural edits may have inserted a new block since the previous partition.
 if(pageDocumentId!==doc.id||!pages.some(p=>p.anchor===pageAnchor)||!currentPageContains(id)&&window.ReaderPagination.pageFor(pages,doc.blocks,id)<0)rebuildPages();
 const target=window.ReaderPagination.pageFor(pages,doc.blocks,id);return target>=0&&showPage(target,{scroll:true});
}
function drawDocument(){
 $('#document-name').textContent=doc.title;$('#document-title').textContent=doc.title;document.title=(floating?'Read · ':'')+doc.title;
 $('#source-link').href=doc.sourceUrl||'#';$('#source-link').hidden=!doc.sourceUrl;
 rebuildPages(selectedId||pageAnchor);renderPage();drawOutline();
 $('#undo').disabled=!historyStack.length;$('#redo').disabled=!redoStack.length;
}
// Structural edits retain every unaffected visible block and the complete model.
function reconcileDocument(){rebuildPages(selectedId||pageAnchor);renderPage();}
window.ReaderPageNavigation={
 ensureBlock:ensureBlockPage,
 visibleIds:()=>visibleBlocks().map(b=>b.id),
 pageInfo:()=>({current:pageIndex,total:pages.length,anchor:pageAnchor,start:pages[pageIndex]?.start,end:pages[pageIndex]?.end}),
 goToPage:page=>showPage(page,{scroll:true,updateURL:true})
};
function updateBlock(id){const b=block(id),old=blockNodes.get(id);if(!b||!old)return;const n=renderBlock(b);old.replaceWith(n);blockNodes.set(id,n);renderedBlocks.set(id,renderKey(b));window.ReaderBlocks?.refreshStructure();window.ReaderInsertion?.refresh();window.ReaderBlockSelection?.refresh();updateSelection();}
function updateSelection(){
 for(const [id,n]of blockNodes){n.classList.toggle('selected',id===selectedId);n.classList.toggle('reading',id===activeId);if(id===activeId)n.setAttribute('aria-current','true');else n.removeAttribute('aria-current');}
 const i=index(activeId||selectedId),queueIndex=selectionQueue?.indexOf(activeId||selectedId)??-1;$('#position').textContent=queueIndex>=0?`${queueIndex+1} / ${selectionQueue.length} selected`:i>=0?`${i+1} / ${doc.blocks.length}`:'Selection';$('#previous').disabled=(controlsOnly&&!selectionQueue)||!readableNeighbor(selectedId,-1);$('#next').disabled=(controlsOnly&&!selectionQueue)||!readableNeighbor(selectedId,1);
 $$('#outline-items button').forEach(n=>n.classList.toggle('active',n.dataset.id===selectedId));
 $('#full-editor-link').href=new URL('index.html#'+encodeURIComponent(selectedId||''),document.baseURI).href;
 if(floating&&selectedId){const b=block(selectedId);$('#source-link').href=b?.sourceBlockUrl||doc.sourceUrl||'#';}
}
function selectBlock(id,{reveal=false,broadcast=true,show=true}={}){
 if(!block(id))return;if(show&&!ensureBlockPage(id))return;selectedId=id;updateSelection();try{if(show||currentPageContains(id))history.replaceState(null,'','#'+encodeURIComponent(id));if(index(id)>=0)localStorage.setItem(storeKey()+':position',JSON.stringify({id,time:Date.now(),owner}));}catch(_e){}
 if(broadcast)channel?.postMessage({type:'select',owner,document:doc.id,id});
 if(reveal)revealBlock(id);if(show&&editorId&&editorId!==id&&!floating)openEditor(id,{focus:false});
}
function revealBlock(id){if(!controlsOnly&&ensureBlockPage(id))blockNodes.get(id)?.scrollIntoView({behavior:'instant',block:floating?'start':'center'});}
function scheduleOutline(){
 outlineDirty=true;clearTimeout(outlineTimer);
 if(!$('#outline').hidden)outlineTimer=setTimeout(drawOutline,250);
}
function outlineLabel(b){
 const label=node('span','outline-label'),source=b.type==='image'?(b.caption||b.image?.alt||b.alt||'Figure'):b.source||'';
 // Render complete TeX first. Visual overflow never truncates formula source.
 label.innerHTML=b.type==='equation'?'<span class="math-inline">'+renderMath(source,false)+'</span>':headingMarkup(source);
 for(const el of label.querySelectorAll('img,figure,video,audio,button,input,select,textarea'))el.remove();
 for(const link of label.querySelectorAll('a')){const span=node('span');span.append(...link.childNodes);link.replaceWith(span);}
 return label;
}
function drawOutline(){
 clearTimeout(outlineTimer);outlineDirty=true;if($('#outline').hidden)return;outlineDirty=false;
 const holder=$('#outline-items');holder.replaceChildren();const term=$('#outline-search').value.toLowerCase();
 for(const b of doc.blocks){if(term?!`${b.source||''} ${b.caption||''}`.toLowerCase().includes(term):b.type!=='heading')continue;const n=node('button','outline-entry');n.append(outlineLabel(b));n.dataset.id=b.id;n.dataset.level=String(b.type==='heading'?Math.min(3,Math.max(1,Number(b.level)||2)):1);n.classList.toggle('active',b.id===selectedId);n.onclick=()=>selectBlock(b.id,{reveal:true});holder.append(n);}
 if(!holder.children.length)holder.append(node('p','help','No matching passages.'));
}
function undoPoint(){if(Date.now()-lastUndo>650||!historyStack.length){historyStack.push(clone(doc));if(historyStack.length>50)historyStack.shift();}lastUndo=Date.now();redoStack.length=0;}
// Snapshot after a short idle interval, with a deadline during continuous typing.
function flushDraft(){
 clearTimeout(draftTimer);clearTimeout(draftDeadline);draftTimer=0;draftDeadline=0;
 if(!draftPending||!doc)return;draftPending=false;
 try{localStorage.setItem(storeKey()+':draft',JSON.stringify({document:doc,revision,dirty,updatedAt:Date.now()}));}catch(_e){}
}
function persistDraft(){
 draftPending=true;clearTimeout(draftTimer);draftTimer=setTimeout(flushDraft,250);
 if(!draftDeadline)draftDeadline=setTimeout(flushDraft,2000);
}
function changed(id,{structural=false,presentation=false,render=true}={}){
 const synced=window.ReaderBlocks?.syncFrom(id)||[];
 if(window.ReaderBlocks?.normalizeStructure()?.length)structural=true;
 dirty=true;pendingSave=true;if(!presentation)for(const changedId of [id,...synced]){cache.delete(changedId);nearbyNarration.invalidate(changedId);}
 if(activeId&&!block(activeId))stopPlayback({message:'The reading passage was removed. Select another passage to continue.'});
 else if(!presentation&&(activeId===id||synced.includes(activeId))&&(playIntent||busy||audio.getAttribute('src')))status('Edit saved for the next reading. The current recording continues.');
 persistDraft();saved('Saving local draft…');clearTimeout(saveTimer);saveTimer=setTimeout(saveDocument,550);
 if(render){structural?reconcileDocument():updateBlock(id);}else if(blockNodes.has(id)){renderedBlocks.set(id,renderKey(block(id)));}
 if(!structural)for(const peerId of synced)updateBlock(peerId);
 if(block(id)?.type==='heading'&&!structural){rebuildPages(pageAnchor);refreshPageControls();}
 $('#undo').disabled=!historyStack.length;$('#redo').disabled=!redoStack.length;if(!presentation||structural)scheduleOutline();
}
function formattingChanged(ids){
 lastUndo=0;
 for(const id of window.ReaderBlocks?.syncMany(ids)||ids)updateBlock(id);
 dirty=true;pendingSave=true;persistDraft();saved('Saving local draft…');clearTimeout(saveTimer);saveTimer=setTimeout(saveDocument,550);
 $('#undo').disabled=!historyStack.length;$('#redo').disabled=!redoStack.length;
}
function clearStaleSpeech(b){
 if(b.speechOverride){delete b.speechOverride;if(editorId===b.id)$('#speech-override').value='';}
}
function sourceTool(kind){
 const b=block(editorId),field=$('#block-source');if(!b||b.type==='equation'||b.type==='image')return;
 const start=field.selectionStart,end=field.selectionEnd,value=field.value,selected=value.slice(start,end);
 const delimiters={bold:['**','**'],italic:['*','*'],underline:['<u>','</u>'],math:['$','$'],equation:['\n\n$$\n','\n$$\n\n'],link:['[','](https://)']};
 const [left,right]=delimiters[kind];let next,selectionStart,selectionEnd;
 if(start>=left.length&&value.slice(start-left.length,start)===left&&value.slice(end,end+right.length)===right){next=value.slice(0,start-left.length)+selected+value.slice(end+right.length);selectionStart=start-left.length;selectionEnd=end-left.length;}
 else{next=value.slice(0,start)+left+selected+right+value.slice(end);selectionStart=start+left.length;selectionEnd=end+left.length;}
 editField('source',next,{keepSpeech:['bold','italic','underline'].includes(kind)});field.value=next;field.focus();field.setSelectionRange(selectionStart,selectionEnd);
}

async function saveDocument(){
 if(saving||!pendingSave||!doc)return;saving=true;pendingSave=false;const sent=JSON.stringify(doc),sentRevision=revision;
 try{const r=await fetch('/api/document',{method:'PUT',headers:{'Content-Type':'application/json'},body:JSON.stringify({document:JSON.parse(sent),revision:sentRevision})});const data=await r.json();
  if(!r.ok)throw Error(r.status===409?'Another view saved a newer draft. Your edits are kept in this browser; keep this view open to recover them.':data.error||'Could not save the draft.');
  revision=data.revision??revision+1;dirty=JSON.stringify(doc)!==sent;pendingSave=dirty;saved(dirty?'Saving local draft…':'Saved locally');persistDraft();flushDraft();channel?.postMessage({type:'saved',owner,document:doc.id,revision});
 }catch(e){saved('Draft kept in this browser',true);status(e.message);pendingSave=false;persistDraft();flushDraft();}finally{saving=false;if(pendingSave)setTimeout(saveDocument,0);}
}
function openEditor(id,{focus=true}={}){
 const b=block(id);if(!b||!ensureBlockPage(id))return;if(index(id)<0){toast('Copied selections are temporary and do not edit the saved manuscript.');return;}if(floating){window.open(new URL('index.html#'+encodeURIComponent(id),document.baseURI).href,'_blank');return;}
 editorId=id;selectedId=id;$('#editor').hidden=false;document.body.classList.add('editing');$('#editor-title').textContent=`Passage ${index(id)+1}`;$('#block-type').value=b.type||'paragraph';$('#heading-level').value=String(b.level||2);$('#block-source').value=b.source||'';$('#image-caption').value=b.caption||'';$('#speech-override').value=b.speechOverride||'';$('#generated-speech').hidden=true;
 $('#block-source-link').href=b.sourceBlockUrl||doc.sourceUrl||'#';$('#block-source-link').hidden=!(b.sourceBlockUrl||doc.sourceUrl);editorFields(b);updateSelection();if(focus)$('#block-source').focus();
}
function editorFields(b){$('#source-tools').hidden=['equation','image'].includes(b.type);$('#heading-level-row').hidden=b.type!=='heading';$('#image-caption-row').hidden=b.type!=='image';$('#source-label').textContent=b.type==='equation'?'LaTeX equation':b.type==='image'?'Local figure path':'Markdown and formulas';$('#editing-help').textContent=b.type==='equation'?'Type LaTeX without delimiters. The rendered equation updates on the page.':b.type==='image'?'Use an existing path in figures/. The caption supports formulas.':'Use $…$ for inline math and $$…$$ for a displayed equation. The page updates as you type.';}
function closeEditor(){$('#editor').hidden=true;document.body.classList.remove('editing');editorId=null;saveDocument();window.ReaderFormatting?.clearSelection();}
function editField(field,value,{keepSpeech=false}={}){const b=block(editorId);if(!b||b[field]===value)return;undoPoint();b[field]=value;if(field==='source'&&b.type==='image'&&b.image){delete b.image;}if(!keepSpeech&&(field==='source'||field==='caption'))clearStaleSpeech(b);changed(b.id,{structural:field==='level'});}
function addBlock(type){undoPoint();const b={id:'local-'+(crypto.randomUUID?.()||Date.now()),type,source:type==='equation'?'x^a = x^i (e_i)^a':'',level:2};const i=index(selectedId),parent=window.ReaderBlocks?.insertionParent(selectedId,doc.blocks[i+1]?.id);if(parent)b.parentId=parent;doc.blocks.splice(i+1,0,b);selectedId=b.id;changed(b.id,{structural:true});revealBlock(b.id);if(window.ReaderDirectEdit&&!floating)window.ReaderDirectEdit.focusBlock(b.id);else openEditor(b.id);}
function moveBlock(delta){const id=editorId;if(!block(id))return;window.ReaderBlocks?.moveBy(id,delta);openEditor(id,{focus:false});}
function undo(redo=false){const source=redo?redoStack:historyStack,target=redo?historyStack:redoStack;if(!source.length)return;target.push(clone(doc));doc=source.pop();lastUndo=0;if(!block(selectedId))selectedId=doc.blocks[0]?.id;changed(selectedId,{structural:true});if(editorId&&block(editorId))openEditor(editorId,{focus:false});else closeEditor();}
function download(name,content,type){const url=URL.createObjectURL(new Blob([content],{type})),a=node('a');a.href=url;a.download=name;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),5000);}
function markdownExport(){return '# '+doc.title+'\n\n'+doc.blocks.map(b=>b.type==='equation'?'$$\n'+b.source+'\n$$':b.type==='heading'?'#'.repeat(Math.min(6,Math.max(1,b.level||2)))+' '+b.source:b.type==='image'?'!['+(b.image?.alt||b.alt||'')+']('+(typeof b.image==='string'?b.image:b.image?.src||b.source)+')\n\n'+(b.caption||''):b.source||'').join('\n\n');}
function setPlayUI(){const playing=!audio.paused&&!audio.ended;document.body.classList.toggle('playing',playing);$('#play').classList.toggle('busy',busy);$('#play').textContent=busy?'…':playing?'Ⅱ':'▶';$('#play').setAttribute('aria-label',busy?'Cancel narration preparation':playing?'Pause narration':'Play narration');}
function stopPlayback({message='',remote=false,keepNearby=false}={}){
 window.ReaderSemanticNarration?.clear();
 if(!keepNearby)nearbyNarration.stop();
 clearFigureCue();const wasOwner=playIntent||busy||(!audio.paused&&!audio.ended);serial++;playIntent=false;busy=false;controller?.abort();controller=null;cancelAnimationFrame(raf);audio.pause();audio.removeAttribute('src');audio.load();activeRecording=null;segments=[];lastCaption=-1;activeId=null;$('#seek').disabled=true;$('#seek').value=0;$('#time').textContent='0:00 / 0:00';$('#caption').textContent='Choose a passage and press Play.';delete $('#caption').dataset.segment;setPlayUI();if(doc)updateSelection();if(message)status(message);if(wasOwner&&!remote)postPlayer(false,false);
}
async function postPlayer(playing,claim=false){
 if(!playerAPI||!doc)return;
 try{const r=await fetch('/api/player-state',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({owner,blockId:activeId||selectedId,playing,time:audio.currentTime||0,claim})});if(r.status===404){playerAPI=false;return;}if(r.ok){const state=await r.json();if(state.owner&&state.owner!==owner&&playIntent)stopPlayback({message:'Another window controls the reading. Press Play to continue here.',remote:true});}}catch(_e){}
}
async function pollPlayer(){
 if(!doc||!playerAPI||(document.hidden&&!playIntent&&!busy))return;
 if((playIntent||busy)&&Date.now()-lastHeartbeat>2000){lastHeartbeat=Date.now();await postPlayer(true,false);}
 try{const r=await fetch('/api/player-state',{cache:'no-store'});if(r.status===404){playerAPI=false;return;}if(!r.ok)return;const s=await r.json();remoteOwner=s.owner||null;
  if(s.owner&&s.owner!==owner&&s.playing){const show=floating||follow||currentPageContains(activeId);if(playIntent||busy||!audio.paused)stopPlayback({message:'Reading is active in another window.',remote:true});if(!editorId&&!window.ReaderDirectEdit?.isEditing()&&!controlsOnly&&block(s.blockId)){activeId=s.blockId;selectBlock(s.blockId,{broadcast:false,show,reveal:follow});status('Reading in another window · press Play to continue here.');}}
  else if(s.owner!==owner&&!s.playing&&!playIntent&&!busy&&audio.paused&&activeId){activeId=null;updateSelection();status('Choose a passage or press Play to read in this window.');}
 }catch(_e){}
}
// Prepare only a small playback neighborhood; typing never starts synthesis.
const warmedAudio=new Map();
function nearbyCandidates(id){
 const result=[];let next=id;
 for(let i=0;i<20;i++){next=readableNeighbor(next,1);if(!next)break;result.push(block(next));}
 const previous=readableNeighbor(id,-1);if(previous)result.push(block(previous));
 return result;
}
function trimWarmedAudio(){
 let bytes=[...warmedAudio.values()].reduce((sum,item)=>sum+item.bytes,0);
 for(const [key,item]of warmedAudio){
  if(warmedAudio.size<=24&&bytes<=96*1024*1024)break;
  if(audio.getAttribute('src')===item.url)continue;
  URL.revokeObjectURL(item.url);warmedAudio.delete(key);bytes-=item.bytes;
 }
}
async function warmRecording(snapshot,record,signal){
 if(fingerprint(block(snapshot.id))!==fingerprint(snapshot)||signal.aborted)return;
 const key=asset(record.audio);if(warmedAudio.has(key))return;
 let blob=await window.ReaderSavedAudio?.get(key);
 if(!blob){
  const response=await fetch(key,{signal});if(!response.ok)return;
  const length=Number(response.headers.get('Content-Length'));
  if(length>32*1024*1024){await response.body?.cancel();return;}
  blob=await response.blob();
  if(!signal.aborted)void window.ReaderSavedAudio?.put(key,blob);
 }
 if(signal.aborted||blob.size>32*1024*1024||fingerprint(block(snapshot.id))!==fingerprint(snapshot))return;
 if(!warmedAudio.has(key))warmedAudio.set(key,{url:URL.createObjectURL(blob),bytes:blob.size});
 trimWarmedAudio();
}
function recordingSource(record){
 const key=asset(record.audio),item=warmedAudio.get(key);
 if(item){warmedAudio.delete(key);warmedAudio.set(key,item);return item.url;}
 return key;
}
const sharedPreparation=window.ReaderSharedPreparation.create({
 keyOf:narrationFingerprint,request:requestRecording,
 onPromote:async(b,signal)=>{const spoken=window.ReaderNarration?await window.ReaderNarration.resolve(b,{signal}):b;await window.ReaderVoiceTransport.start(spoken,{signal,background:false});}
});
function prepare(b,signal,options={}){
 if(signal?.aborted)return Promise.reject(new DOMException('Aborted','AbortError'));
 b={...clone(b),_readerMode:window.ReaderSemanticNarration?.cacheTag(b)||'guided'};
 const key=narrationFingerprint(b),entry=cache.get(b.id);
 if(!options.refresh&&entry?.key===key&&narrationCurrent(b,entry.record))return (options.background?Promise.resolve():warmRecording(b,entry.record,signal||new AbortController().signal).catch(error=>{if(error.name==='AbortError')throw error;})).then(()=>entry.record);
 return sharedPreparation.prepare(b,signal,options);
}
const nearbyNarration=window.ReaderNearbyNarration.create({
 getCandidates:nearbyCandidates,prepare,fingerprint,
 isCurrent:id=>playIntent&&activeId===id&&!audio.paused,
 onReady:warmRecording
});
window.addEventListener('reader-voice-change',()=>{nearbyNarration.stop();cache.clear();status('Voice saved. It applies to the next passage you start.');});
async function requestRecording(b,signal,{background=false,refresh=false}={}){
 const key=narrationFingerprint(b);if(!refresh&&cache.has(b.id)&&cache.get(b.id).key===key&&narrationCurrent(b,cache.get(b.id).record))return cache.get(b.id).record;
 const spoken=window.ReaderNarration?await window.ReaderNarration.resolve(b,{signal,background,refresh}):b;
 const voiceRequest=await window.ReaderVoiceTransport.start(spoken,{signal,background});let record=voiceRequest.record;
 if(record.jobId&&!record.audio&&!record.result?.audio){const jobId=record.jobId;const started=Date.now();for(;;){if(signal.aborted)throw new DOMException('Aborted','AbortError');const r=await fetch(voiceRequest.jobURL||('/api/jobs/'+encodeURIComponent(jobId)),{signal,cache:'no-store'});record=await r.json();if(!r.ok||record.status==='error')throw Error(record.error||'Narration preparation failed.');if(record.status==='ready'||record.audio||record.result?.audio){record=record.result||record;break;}if(Date.now()-started>600000)throw Error('Narration is taking longer than expected. Try this passage again.');await new Promise((resolve,reject)=>{const t=setTimeout(resolve,550);signal.addEventListener('abort',()=>{clearTimeout(t);reject(new DOMException('Aborted','AbortError'));},{once:true});});}}
 record=record.result||record;
 if(!record.audio)throw Error('The reading is ready but has no audio file.');
 if(spoken._readerNarration)record={...record,readerNarration:spoken._readerNarration};
 if(!narrationCurrent(block(b.id),record))throw Error('The explanation context changed. Press Play to read the latest source.');
 if(block(b.id)&&fingerprint(block(b.id))===fingerprint(b))cache.set(b.id,{key:narrationFingerprint(b),record});
 if(!background)try{await warmRecording(b,record,signal);}catch(error){if(error.name==='AbortError')throw error;}
 return record;
}
async function playFrom(id,keepGoing=true,claim=true){
 if(selectionQueue&&!selectionQueue.includes(id))selectionQueue=null;
 if(!block(id))return;
 id=readableNeighbor(id,1,true);
 if(!id){stopPlayback({message:selectionQueue?'No readable passages remain in this selection.':'No readable passages remain. Select an earlier passage to replay.'});return;}
 const show=claim||(!editorId&&!window.ReaderDirectEdit?.isEditing()&&(follow||currentPageContains(activeId)));const b=block(id);stopPlayback({remote:true,keepNearby:true});activeId=id;selectBlock(id,{reveal:follow,broadcast:claim,show});continuous=keepGoing;playIntent=true;busy=true;const token=serial;controller=new AbortController();setPlayUI();status('Preparing source-language reading…');$('#caption').textContent='Preparing this passage and its formulas…';if(claim)channel?.postMessage({type:'claim',owner,document:doc.id,id});await postPlayer(true,claim);
 try{if(token!==serial||!playIntent||!controller)return;const snapshot={...clone(b),_readerMode:window.ReaderSemanticNarration?.mode()||'guided'},pending=prepare(snapshot,controller.signal);nearbyNarration.stop();const record=await pending;if(token!==serial||!playIntent)return;if(fingerprint(block(id))!==fingerprint(snapshot)){busy=false;playIntent=false;setPlayUI();status('This passage changed while preparing. Press Play to read the latest text.');postPlayer(false,false);return;}activeRecording={key:fingerprint(snapshot),narrationKey:narrationFingerprint(snapshot),record,figureCues:validatedFigureCues(snapshot,record)};segments=record.segments||[];lastCaption=-1;audio.src=recordingSource(record);setPlaybackSpeed($('#speed').value);audio.load();enforcePlaybackSpeed();trimWarmedAudio();$('#seek').disabled=false;busy=false;setPlayUI();
  if(editorId===id&&record.speech){$('#generated-speech').textContent='Generated reading: '+(Array.isArray(record.speech)?record.speech.join(' '):record.speech);$('#generated-speech').hidden=false;}
  try{enforcePlaybackSpeed();await audio.play();if(token!==serial){audio.pause();return;}status(selectionQueue?`Reading selection ${selectionQueue.indexOf(id)+1} / ${selectionQueue.length} · stops after this selection.`:'Reading · '+(floating?'Source-language narration':'Click a passage to choose where to continue.'));}catch(_e){playIntent=false;await postPlayer(false,false);status('Narration ready. Press Play to begin.');setPlayUI();}syncAudio();
 }catch(e){if(token!==serial)return;busy=false;playIntent=false;setPlayUI();status(e.name==='AbortError'?'Reading preparation cancelled.':e.message);$('#caption').textContent='This passage remains available to read and edit.';postPlayer(false,false);}
}
async function togglePlay(){
 if(!doc?.blocks.length){status('Add or select a passage before starting narration.');return;}
 if(controlsOnly&&!selectionQueue&&!audio.getAttribute('src')&&!busy){status('Copy a passage in Codex Notion, then choose Read copied selection.');return;}
 if(busy){stopPlayback({message:'Preparation paused. Press Play to try again.'});return;}
 if(!audio.paused&&!audio.ended){playIntent=false;audio.pause();status('Paused · text and formulas remain selectable.');postPlayer(false,false);return;}
 if(audio.ended&&selectionQueue?.length>1){playFrom(selectionQueue[0],true);return;}
 if(audio.getAttribute('src')&&activeId===selectedId&&activeRecording&&activeRecording.key===fingerprint(block(activeId))&&activeRecording.narrationKey===narrationFingerprint(block(activeId))&&narrationCurrent(block(activeId),activeRecording.record)&&(!cache.has(activeId)||cache.get(activeId).record===activeRecording.record)){playIntent=true;const token=serial;channel?.postMessage({type:'claim',owner,document:doc.id,id:activeId});await postPlayer(true,true);if(token!==serial||!playIntent)return;if(audio.ended)audio.currentTime=0;enforcePlaybackSpeed();audio.play().catch(()=>{playIntent=false;postPlayer(false,false);status('Press Play to begin narration.');});}
 else playFrom(selectedId||doc.blocks[0]?.id,true);
}
// Cue ranges address measured synthesis segments, not estimated word durations.
function validatedFigureCues(b,record){
 const data=b.figureCues,actual=record.segments||[];
 if(!data||!Array.isArray(data.segments)||data.segments.length!==actual.length||!data.segments.every((text,i)=>text===actual[i].text)||!Array.isArray(data.cues))return [];
 return data.cues.filter(c=>Number.isInteger(c.from)&&Number.isInteger(c.to)&&c.from>=0&&c.to>c.from&&c.to<=actual.length&&Array.isArray(c.rect)&&c.rect.length===4&&c.rect.every(Number.isFinite)&&c.rect[0]>=0&&c.rect[1]>=0&&c.rect[2]>0&&c.rect[3]>0&&c.rect[0]+c.rect[2]<=1&&c.rect[1]+c.rect[3]<=1);
}
function clearFigureCue(){if(figureCueNode){figureCueNode.hidden=true;figureCueNode=null;}}
function syncFigureCue(i){
 const b=block(activeId),overlay=blockNodes.get(activeId)?.querySelector('.figure-cue');
 const cue=activeRecording?.key===fingerprint(b)?activeRecording?.figureCues?.find(c=>i>=c.from&&i<c.to):null;
 if(!overlay||!cue){clearFigureCue();return;}
 if(figureCueNode!==overlay){clearFigureCue();figureCueNode=overlay;}
 const [x,y,w,h]=cue.rect;overlay.style.left=`${100*x}%`;overlay.style.top=`${100*y}%`;overlay.style.width=`${100*w}%`;overlay.style.height=`${100*h}%`;overlay.dataset.feature=cue.target;overlay.title=cue.label||'';overlay.hidden=false;
}
function syncAudio(){
 const t=audio.currentTime||0,d=Number.isFinite(audio.duration)?audio.duration:0;$('#time').textContent=`${time(t)} / ${time(d)}`;$('#seek').value=d?Math.min(1000,1000*t/d):0;
 let i=-1;for(let k=0;k<segments.length;k++){if(t+.001>=Number(segments[k].start))i=k;else break;}
 if(i!==lastCaption){lastCaption=i;$('#caption').textContent=i>=0?segments[i].text||'':'';$('#caption').dataset.segment=String(i);}
 syncFigureCue(i);window.ReaderSemanticNarration?.sync(block(activeId),activeRecording?.record,i,blockNodes.get(activeId));setPlayUI();
}
function tick(){syncAudio();if(!audio.paused&&!audio.ended)raf=requestAnimationFrame(tick);}
function nextPassage(delta,autoplay){if(controlsOnly&&!selectionQueue)return;const id=readableNeighbor(selectedId||activeId,delta);if(!id)return;if(autoplay)playFrom(id,true);else{stopPlayback({message:'Passage selected. Press Play to listen.'});selectBlock(id,{reveal:!floating});}}
async function loadDocument({focusRefresh=false}={}){
 if(loadingDocument||focusRefresh&&(dirty||saving||editorId||window.ReaderDirectEdit?.isEditing()||playIntent||busy))return;loadingDocument=true;
 try{const r=await fetch('/api/document',{cache:'no-store'});if(!r.ok)throw Error('Could not open the local manuscript.');const data=await r.json();if(focusRefresh&&(dirty||saving||editorId||window.ReaderDirectEdit?.isEditing()||playIntent||busy))return;if(focusRefresh&&Number(data.revision)===revision)return;const next=data.document||data;if(!Array.isArray(next.blocks))throw Error('The document has no editable blocks.');if(focusRefresh){historyStack.length=0;redoStack.length=0;lastUndo=0;window.ReaderFormatting?.clearSelection();}doc=next;revision=data.revision??0;
  let recoveredConflict=false;
  if(!focusRefresh){let draft;try{draft=JSON.parse(localStorage.getItem(storeKey()+':draft')||'null');}catch(_e){}if(draft?.dirty&&Array.isArray(draft.document?.blocks)){const migrated=await window.ReaderSentenceMigration?.recoverDraft(doc,revision,draft);if(migrated){doc=migrated;dirty=true;pendingSave=true;toast('Recovered your edits and the sentence layout update.');setTimeout(saveDocument,800);}else{doc=draft.document;dirty=true;if(draft.revision===revision){pendingSave=true;toast('Recovered your unsaved local draft.');setTimeout(saveDocument,800);}else{revision=draft.revision;pendingSave=false;recoveredConflict=true;toast('Recovered conflicting edits. They remain in this view and in the browser recovery draft.');}}}}
  if(!focusRefresh)restoreSettings();window.ReaderSemanticNarration?.refreshControl();
  let hash;try{hash=decodeURIComponent(location.hash.slice(1));}catch(_e){}let remembered;try{remembered=JSON.parse(localStorage.getItem(storeKey()+':position')||'null')?.id;}catch(_e){}
  let rememberedPage;try{rememberedPage=JSON.parse(localStorage.getItem(storeKey()+':page')||'null')?.anchor;}catch(_e){}
  selectedId=block(hash)?hash:block(selectedId)?selectedId:block(rememberedPage)?rememberedPage:block(remembered)?remembered:doc.blocks[0]?.id;drawDocument();$('#loading').hidden=true;$('#manuscript').hidden=false;saved(recoveredConflict?'Conflicting draft kept locally':dirty?'Unsaved local draft':'Saved locally',recoveredConflict);if(recoveredConflict)status('A newer saved copy exists. Your conflicting edits are kept in this browser; keep this view open to recover them.');selectBlock(selectedId,{broadcast:false,reveal:!focusRefresh&&Boolean(block(hash))});
 }catch(e){$('#loading').hidden=false;$('#loading').textContent=e.message;status('Open this reader through its local server.');}finally{loadingDocument=false;}
}
window.readerSetControlsMode=function(enabled){
 controlsOnly=Boolean(enabled);document.body.classList.toggle('controls-only',controlsOnly);
 if(controlsOnly){follow=false;$('#follow').setAttribute('aria-pressed','false');$('#follow').textContent='Follow off';}
 if(doc){updateSelection();refreshPageControls();}return controlsOnly;
};
window.readerReadCopiedSelection=async function(selection={}){
 const text=typeof selection.text==='string'?selection.text:'',html=typeof selection.html==='string'?selection.html:'';
 if(!text.trim()&&!html.trim()){status('Nothing copied. Select a passage in Codex Notion, copy it, then choose Read copied selection.');return {ok:false,error:'No copied selection.'};}
 const request=++selectionRequest;stopPlayback({message:'Opening copied selection…'});
 try{
  await documentReady;if(!doc)throw Error('The local manuscript has not loaded.');
  await loadDocument({focusRefresh:true});if(request!==selectionRequest)return {ok:false,cancelled:true};
  const response=await fetch('/api/selection',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({text,html})});const result=await response.json();
  if(request!==selectionRequest)return {ok:false,cancelled:true};if(!response.ok)throw Error(result.error||'This copied selection could not be read.');
  if(result.kind==='source'){
   const ids=Array.isArray(result.blockIds)&&result.blockIds.length?[...new Set(result.blockIds)]:[result.startBlockId];
   if(ids[0]!==result.startBlockId||ids.some(id=>!doc.blocks.some(b=>b.id===id)))throw Error('The matched source is unavailable. Copy the passage again; no nearby passage was substituted.');
   temporaryBlock=null;selectionQueue=ids.filter(id=>isNarratable(block(id)));
  }else if(result.kind==='selection'&&result.block&&['paragraph','equation'].includes(result.block.type)&&typeof result.block.source==='string'&&result.block.source.trim()){
   temporaryBlock=clone(result.block);selectionQueue=[temporaryBlock.id];
  }else throw Error('The copied selection did not provide a readable passage.');
  if(!selectionQueue.some(id=>isNarratable(block(id)))){stopPlayback({message:'The selection contains only empty blocks.'});return {ok:true,kind:result.kind,blockIds:[],skipped:true};}
  drawDocument();selectBlock(selectionQueue[0],{broadcast:false});await playFrom(selectionQueue[0],true,true);
  return {ok:true,kind:result.kind,blockIds:[...selectionQueue]};
 }catch(e){if(request===selectionRequest){status(e.message);$('#caption').textContent='Copy the intended passage again. The saved manuscript has not changed.';}return {ok:false,error:e.message};}
};
$('#manuscript').hidden=true;
for(const button of $$('[data-page-prev]'))button.onclick=()=>showPage(pageIndex-1,{scroll:true,updateURL:true});
for(const button of $$('[data-page-next]'))button.onclick=()=>showPage(pageIndex+1,{scroll:true,updateURL:true});
for(const select of $$('[data-page-select]'))select.onchange=()=>showPage(Number(select.value),{scroll:true,updateURL:true});
// Printing explicitly expands the full document, then restores its small view.
window.addEventListener('beforeprint',()=>{if(!doc||printing)return;if(!preparePageChange())return;printPosition={anchor:pageAnchor,x:scrollX,y:scrollY};printing=true;renderPage();});
window.addEventListener('afterprint',()=>{if(!printing)return;printing=false;rebuildPages(printPosition?.anchor);renderPage();if(printPosition)scrollTo(printPosition.x,printPosition.y);printPosition=null;});
function closePastePanel(){$('#paste-selection-panel').hidden=true;$('#paste-selection-toggle').setAttribute('aria-expanded','false');}
$('#paste-selection-toggle').onclick=()=>{const open=$('#paste-selection-panel').hidden;$('#paste-selection-panel').hidden=!open;$('#paste-selection-toggle').setAttribute('aria-expanded',String(open));if(open)$('#paste-selection-text').focus();};
$('#paste-selection-close').onclick=closePastePanel;
$('#read-pasted-selection').onclick=async()=>{const text=$('#paste-selection-text').value;if(!text.trim()){status('Paste the passage you want to hear.');$('#paste-selection-text').focus();return;}closePastePanel();$('#read-pasted-selection').disabled=true;try{await window.readerReadCopiedSelection({text,html:''});}finally{$('#read-pasted-selection').disabled=false;}};
$('#outline-toggle').onclick=()=>{setOutline($('#outline').hidden);saveSettings({outline:!$('#outline').hidden});};$('#outline-search').oninput=drawOutline;
$$('[data-source-tool]').forEach(button=>{button.onmousedown=e=>e.preventDefault();button.onclick=()=>sourceTool(button.dataset.sourceTool);});
$('#block-source').addEventListener('keydown',e=>{if((e.metaKey||e.ctrlKey)&&['b','i','u'].includes(e.key.toLowerCase())&&!['equation','image'].includes(block(editorId)?.type)){e.preventDefault();sourceTool({b:'bold',i:'italic',u:'underline'}[e.key.toLowerCase()]);}});
$('#editor-close').onclick=closeEditor;$('#block-source').oninput=e=>editField('source',e.target.value);$('#image-caption').oninput=e=>editField('caption',e.target.value);$('#speech-override').oninput=e=>editField('speechOverride',e.target.value);$('#heading-level').onchange=e=>editField('level',Number(e.target.value));
$('#block-type').onchange=e=>{const b=block(editorId);if(!b)return;undoPoint();const textTypes=['paragraph','heading','callout'],presentation=textTypes.includes(b.type)&&textTypes.includes(e.target.value);b.type=e.target.value;delete b.blockStyle;if(!presentation)clearStaleSpeech(b);changed(b.id,{structural:true,presentation});editorFields(b);};$('#read-edited').onclick=()=>playFrom(editorId,true);$('#edit-selected').onclick=()=>openEditor(selectedId);$('#add-text').onclick=()=>addBlock('paragraph');$('#add-equation').onclick=()=>addBlock('equation');$('#move-up').onclick=()=>moveBlock(-1);$('#move-down').onclick=()=>moveBlock(1);
$('#duplicate-block').onclick=()=>{const id=window.ReaderBlocks?.duplicateBlock(editorId);if(id)openEditor(id);};
$('#delete-block').onclick=()=>window.ReaderBlocks?.deleteBlock(editorId);
$('#undo').onclick=()=>undo();$('#redo').onclick=()=>undo(true);$('#export-menu-toggle').onclick=()=>{const open=$('#export-menu').hidden;$('#export-menu').hidden=!open;$('#export-menu-toggle').setAttribute('aria-expanded',String(open));};$('#export-json').onclick=()=>{download('notion-manuscript.json',JSON.stringify(doc,null,2),'application/json');$('#export-menu').hidden=true;};$('#export-markdown').onclick=()=>{download('notion-manuscript.md',markdownExport(),'text/markdown');$('#export-menu').hidden=true;};$('#print-document').onclick=()=>{if(!preparePageChange())return;$('#export-menu').hidden=true;window.print();};
$('#reload-saved').onclick=()=>{if(!doc||saving){toast('Wait for the current save to finish.');return;}if(dirty){download('notion-manuscript-unsaved-backup.json',JSON.stringify(doc,null,2),'application/json');try{localStorage.setItem(storeKey()+':conflict-backup',JSON.stringify({document:doc,revision,updatedAt:Date.now()}));}catch(_e){}}clearTimeout(saveTimer);clearTimeout(draftTimer);clearTimeout(draftDeadline);draftTimer=0;draftDeadline=0;draftPending=false;pendingSave=false;dirty=false;try{localStorage.removeItem(storeKey()+':draft');}catch(_e){}historyStack.length=0;redoStack.length=0;closeEditor();stopPlayback();$('#export-menu').hidden=true;loadDocument();};
$('#play').onclick=togglePlay;$('#previous').onclick=()=>nextPassage(-1,playIntent||busy||!audio.paused);$('#next').onclick=()=>nextPassage(1,playIntent||busy||!audio.paused);$('#position').onclick=()=>revealBlock(activeId||selectedId);$('#follow').onclick=()=>{setFollow(!follow,{reveal:true});if(!controlsOnly)saveSettings({follow});};
$('#speed').onchange=()=>{setPlaybackSpeed($('#speed').value);saveSettings({speed:Number($('#speed').value)});};$('#seek').oninput=e=>{if(Number.isFinite(audio.duration)){audio.currentTime=audio.duration*Number(e.target.value)/1000;syncAudio();}};
audio.onplay=()=>{enforcePlaybackSpeed();playIntent=true;cancelAnimationFrame(raf);tick();status(selectionQueue?`Reading selection ${selectionQueue.indexOf(activeId)+1} / ${selectionQueue.length} · stops after this selection.`:'Reading · Original-language narration.');};audio.onpause=()=>{cancelAnimationFrame(raf);if(!playIntent)nearbyNarration.stop();syncAudio();};audio.onloadedmetadata=()=>{enforcePlaybackSpeed();syncAudio();};audio.onseeking=audio.onseeked=audio.ontimeupdate=syncAudio;
audio.addEventListener('loadstart',enforcePlaybackSpeed);audio.addEventListener('ratechange',enforcePlaybackSpeed);
audio.addEventListener('playing',()=>{enforcePlaybackSpeed();if(playIntent)nearbyNarration.update(activeId);});
audio.onended=()=>{cancelAnimationFrame(raf);syncAudio();const next=readableNeighbor(activeId,1);if(playIntent&&continuous&&next)playFrom(next,true,false);else{playIntent=false;nearbyNarration.stop();postPlayer(false,false);status(selectionQueue?'Selection complete. Copy another passage or replay this selection.':'Reading complete. Replay or select another passage.');}};
audio.onerror=()=>{if(audio.getAttribute('src')){busy=false;playIntent=false;nearbyNarration.stop();setPlayUI();status('The audio could not load. Try this passage again.');postPlayer(false,false);}};
document.addEventListener('keydown',e=>{if(!doc||e.defaultPrevented||e.repeat)return;if(e.key==='Escape'&&!$('#paste-selection-panel').hidden){e.preventDefault();closePastePanel();return;}const typing=e.target.closest('textarea,input,select,[contenteditable]');if((e.metaKey||e.ctrlKey)&&!e.shiftKey&&e.key.toLowerCase()==='e'&&!floating){e.preventDefault();openEditor(selectedId);return;}if(typing)return;if((e.metaKey||e.ctrlKey)&&e.key.toLowerCase()==='z'){e.preventDefault();undo(e.shiftKey);return;}if(e.metaKey||e.ctrlKey||e.altKey)return;if(window.getSelection()?.toString())return;if(e.key==='Escape'){closeEditor();return;}if((e.code==='Space'||e.key===' ')&&!e.target.closest('button,a,summary')){e.preventDefault();togglePlay();}else if(e.key==='ArrowRight'){e.preventDefault();nextPassage(1,playIntent);}else if(e.key==='ArrowLeft'){e.preventDefault();nextPassage(-1,playIntent);}});
window.addEventListener('hashchange',()=>{let id;try{id=decodeURIComponent(location.hash.slice(1));}catch(_e){}if(block(id))selectBlock(id,{reveal:!floating,broadcast:false});});
window.addEventListener('focus',()=>loadDocument({focusRefresh:true}));window.addEventListener('pagehide',flushDraft);document.addEventListener('visibilitychange',()=>{if(document.hidden)flushDraft();});window.addEventListener('beforeunload',e=>{flushDraft();if(dirty){e.preventDefault();e.returnValue='';}if(playIntent)postPlayer(false,false);});
if(channel)channel.onmessage=e=>{const m=e.data;if(m.owner===owner||!doc||m.document!==doc.id)return;if(m.type==='claim'){const show=floating||follow||currentPageContains(activeId);stopPlayback({message:'Reading is active in another window.',remote:true});if(!editorId&&!window.ReaderDirectEdit?.isEditing()&&!controlsOnly&&block(m.id)){activeId=m.id;selectBlock(m.id,{broadcast:false,show,reveal:follow});}}else if(m.type==='select'&&!editorId&&!window.ReaderDirectEdit?.isEditing()&&!playIntent&&!controlsOnly&&block(m.id))selectBlock(m.id,{broadcast:false,show:floating||follow,reveal:follow});else if(m.type==='saved')loadDocument({focusRefresh:true});};
new ResizeObserver(()=>document.documentElement.style.setProperty('--player',Math.ceil($('.player').getBoundingClientRect().height)+'px')).observe($('.player'));
const pageNavigation=$('.page-navigation-top');if(pageNavigation)new ResizeObserver(()=>document.documentElement.style.setProperty('--page-navigation-height',Math.ceil(pageNavigation.getBoundingClientRect().height)+'px')).observe(pageNavigation);
window.ReaderFormatting?.init({getDocument:()=>doc,onBeforeChange:()=>{lastUndo=0;undoPoint();},onChange:formattingChanged,onEdit:id=>openEditor(id)});

function commitDirectBlock(id,patch,{render=false}={}){
 const b=block(id);if(!b)return;
 const speechChanged=['source','caption','type'].some(k=>k in patch&&patch[k]!==b[k]);
 const layoutChanged=['type','level','parentId'].some(k=>k in patch&&patch[k]!==b[k]);
 for(const [key,value] of Object.entries(patch)){if(value===null)delete b[key];else b[key]=value;}
 if(speechChanged)clearStaleSpeech(b);
 changed(id,{render,structural:layoutChanged,presentation:!speechChanged});
}
function structuralDirectEdit(id,mutate){
 mutate();selectedId=id;changed(id,{structural:true});
}
window.ReaderSemanticNarration?.init({getDocument:()=>doc,asset,toast,renderMath,onModeChange:()=>{nearbyNarration.stop();if(playIntent&&!audio.paused)nearbyNarration.update(activeId);toast('Narration style applies to the next reading.');}});
window.ReaderNarration?.init({
 getDocument:()=>doc&&({...doc,blocks:doc.blocks.filter(isNarratable)}),getSelectedId:()=>selectedId,getBlock:block,fingerprint,prepare,
 select:id=>selectBlock(id,{reveal:true,broadcast:false})
});
window.ReaderImages?.init({
 getDocument:()=>doc,getBlock:block,getSelectedId:()=>selectedId,toast,
 select:id=>selectBlock(id,{broadcast:false}),
 beforeChange:()=>{lastUndo=0;undoPoint();},
 insertAfter:(anchor,blocks,{focus=true}={})=>{const i=index(anchor);if(i<0)return;const parent=window.ReaderBlocks?.insertionParent(anchor,doc.blocks[i+1]?.id);if(parent)for(const item of blocks)item.parentId=parent;doc.blocks.splice(i+1,0,...blocks);const id=blocks.at(-1).id;if(focus)selectedId=id;changed(id,{structural:true,presentation:true});if(focus)blockNodes.get(id)?.querySelector('[data-reader-image]')?.focus({preventScroll:true});}
});
window.ReaderBlocks?.init({
 getDocument:()=>doc,getNode:id=>blockNodes.get(id),toast,edit:openEditor,
 beforeChange:()=>{lastUndo=0;undoPoint();},changed,afterChange:()=>{lastUndo=0;},
 select:id=>{if(editorId)closeEditor();selectBlock(id,{broadcast:false});}
});
window.ReaderInsertion?.init({
 getDocument:()=>doc,getNode:id=>blockNodes.get(id),
 beforeChange:()=>{lastUndo=0;undoPoint();},afterChange:()=>{lastUndo=0;},
 insert:(beforeId,afterId,next)=>{let at=afterId?index(afterId):doc.blocks.length;if(at<0)return;doc.blocks.splice(at,0,next);selectedId=next.id;changed(next.id,{structural:true,presentation:true});selectBlock(next.id,{broadcast:false});}
});
window.ReaderBlockSelection?.init({
 getDocument:()=>doc,getNode:id=>blockNodes.get(id),getSelectedId:()=>selectedId,
 select:id=>selectBlock(id,{broadcast:false}),deleteBlocks:ids=>window.ReaderBlocks?.deleteBlocks(ids),
 undo:redo=>undo(redo),toast
});
window.ReaderDirectEdit?.init({
 getDocument:()=>doc,getBlock:block,renderMath,toast,
 select:id=>{if(editorId)closeEditor();selectBlock(id,{broadcast:false});},
 onBeforeChange:(_id,options={})=>{if(options.group!=='typing')lastUndo=0;undoPoint();},
 commitBlock:commitDirectBlock,
 renderBlock:updateBlock,
 insertAfter:(id,next)=>{const i=index(id);if(i<0)return;delete next.sourceBlockUrl;const parent=window.ReaderBlocks?.insertionParent(id,doc.blocks[i+1]?.id);if(parent)next.parentId=parent;structuralDirectEdit(next.id,()=>doc.blocks.splice(i+1,0,next));},
 split:(id,patch,next)=>{if(index(id)<0)return;structuralDirectEdit(next.id,()=>{const result=window.ReaderBlocks.split(id,patch,next);for(const changedId of result.changed)cache.delete(changedId);});},
 merge:(previousId,currentId,patch)=>{if(index(currentId)<0||!block(previousId))return;structuralDirectEdit(previousId,()=>{const result=window.ReaderBlocks.merge(previousId,currentId,patch);for(const changedId of result.changed)cache.delete(changedId);if(result.detached)toast('Merged locally. Other synced copies remain independent of this merge.');});},
 undo:redo=>{const selection=window.ReaderDirectEdit.snapshotSelection();undo(redo);window.ReaderDirectEdit.restoreSelection(selection);}
});
window.ReaderWorkspaceEditor={
 ready:()=>Boolean(doc),documentId:()=>doc?.id,selected:()=>selectedId,block,blocks:()=>doc?.blocks||[],
 busy:()=>Boolean(dirty||saving||editorId||window.ReaderDirectEdit?.isEditing()),
 flush:()=>window.ReaderDirectEdit?.flush(),notice:toast,
 patch:(id,patch)=>{if(!block(id))throw Error('Block was removed.');lastUndo=0;undoPoint();commitDirectBlock(id,patch,{render:true});},
 insert:(anchor,next)=>{window.ReaderDirectEdit?.flush();const i=index(anchor);if(i<0)throw Error('Select an existing block first.');lastUndo=0;undoPoint();const parent=window.ReaderBlocks?.insertionParent(anchor,doc.blocks[i+1]?.id);if(parent)next.parentId=parent;doc.blocks.splice(i+1,0,next);selectedId=next.id;changed(next.id,{structural:true});return next.id;}
};
setInterval(pollPlayer,1000);documentReady=loadDocument();
})();
