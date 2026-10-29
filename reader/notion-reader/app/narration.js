/* On-demand source-language readings and explicit background preparation. No typing hooks. */
(function(){'use strict';
 let api={},panel=null,button=null,progress=null,run=null;
 const clone=b=>JSON.parse(JSON.stringify(b));
 const readable=b=>b.type==='image'||Boolean(String(b.source||'').replace(/<[^>]*>/g,'').trim());
 const needsPreparation=b=>!b.speechOverride?.trim() || /[\u3400-\u9fff]/.test(b.type==='image'?(b.caption||b.alt||''):b.source||'')&&!/[\u3400-\u9fff]/.test(b.speechOverride);
 // The local planner resolves source language and ignores translation overrides.
 // Preparation never contacts a translation service or changes the manuscript.
 async function resolve(b,{signal}={}){
  if(signal?.aborted)throw new DOMException('Cancelled','AbortError');
  return window.ReaderSemanticNarration?await window.ReaderSemanticNarration.resolve(b,{signal}):b;
 }
 function toggle(){panel.hidden=!panel.hidden;button.setAttribute('aria-expanded',String(!panel.hidden));if(!panel.hidden)counts();}
 function counts(){if(!panel||run)return;const blocks=api.getDocument?.()?.blocks||[],edited=blocks.filter(b=>readable(b)&&needsPreparation(b));panel.querySelector('[data-prepare="edited"]').textContent=`Prepare edited blocks (${edited.length})`;}
 function update(text){progress.textContent=text;progress.title=text;}
 async function prepare(mode){
  if(run)return;window.ReaderDirectEdit?.flush();
  const doc=api.getDocument?.();if(!doc)return;
  const blocks=doc.blocks.filter(b=>readable(b)&&(mode==='all'||mode==='selected'&&b.id===api.getSelectedId?.()||mode==='edited'&&needsPreparation(b))).map(clone);
  if(!blocks.length){update('No blocks need preparation.');return;}
  const controller=new AbortController();run={controller,done:0,failed:[],stale:0,total:blocks.length};const current=run;
  panel.querySelectorAll('[data-prepare]').forEach(b=>b.disabled=true);panel.querySelector('[data-cancel]').hidden=false;
  update(`Preparing 0 / ${blocks.length}…`);let cursor=0;
  async function worker(){
   while(cursor<blocks.length&&!controller.signal.aborted){const b=blocks[cursor++];
    try{await api.prepare?.(b,controller.signal,{background:true,refresh:mode==='selected'});if(api.fingerprint?.(api.getBlock?.(b.id))!==api.fingerprint?.(b))current.stale++;}
    catch(e){if(e.name==='AbortError')break;current.failed.push({id:b.id,error:e.message});}
    current.done++;update(`Prepared ${current.done-current.failed.length-current.stale} / ${current.total}${current.failed.length?` · ${current.failed.length} failed`:''}${current.stale?` · ${current.stale} changed while preparing`:''}`);
   }
  }
  // A small bounded queue prepares local recordings without typing hooks.
  await Promise.all(Array.from({length:Math.min(3,blocks.length)},worker));
  const label=controller.signal.aborted?'Stopped':'Finished';
  update(`${label}: ${current.done-current.failed.length-current.stale} ready${current.failed.length?`, ${current.failed.length} failed`:''}${current.stale?`, ${current.stale} changed; prepare again`:''}.`);
  const errors=panel.querySelector('.narration-errors');errors.replaceChildren();
  for(const failure of current.failed){const row=document.createElement('button');row.className='narration-failure';row.textContent=`Show failed passage: ${failure.error}`;row.onclick=()=>api.select?.(failure.id);errors.append(row);}
  run=null;panel.querySelectorAll('[data-prepare]').forEach(b=>b.disabled=false);panel.querySelector('[data-cancel]').hidden=true;counts();
 }
 function init(options){api=options;const toolbar=document.querySelector('.toolbar-actions');if(!toolbar||document.body.classList.contains('floating'))return;
  button=document.createElement('button');button.id='narration-toggle';button.className='quiet';button.textContent='Narration';button.setAttribute('aria-expanded','false');button.setAttribute('aria-controls','narration-preparation');button.onclick=toggle;toolbar.append(button);
  panel=document.createElement('section');panel.id='narration-preparation';panel.className='narration-preparation';panel.hidden=true;panel.setAttribute('aria-label','Prepare narration');
  const title=document.createElement('strong');title.textContent='Prepare narration';const help=document.createElement('p');help.textContent='Passages are read in their original language: English or Mandarin. The next three readable blocks and previous block preload during playback. You can also prepare recordings here in advance.';panel.append(title,help);
  for(const [mode,text] of [['selected','Refresh selected block'],['edited','Prepare edited blocks'],['all','Prepare all blocks']]){const b=document.createElement('button');b.type='button';b.dataset.prepare=mode;b.textContent=text;b.onclick=()=>void prepare(mode);panel.append(b);}
  const cancel=document.createElement('button');cancel.type='button';cancel.dataset.cancel='';cancel.textContent='Stop preparation';cancel.hidden=true;cancel.onclick=()=>run?.controller.abort();panel.append(cancel);
  progress=document.createElement('p');progress.className='narration-progress';progress.setAttribute('role','status');progress.textContent='Existing recordings are reused when their text is unchanged.';const errors=document.createElement('div');errors.className='narration-errors';panel.append(progress,errors);toolbar.closest('.document-toolbar').after(panel);
 }
 window.ReaderNarration={resolve,init,testing:{prepare,readable,needsPreparation}};
})();
