/* Saved Library files, mixed-language dictation, shared blocks and editable ink. */
(() => {
 'use strict';
 const $=s=>document.querySelector(s), el=(tag,text)=>{const n=document.createElement(tag);if(text)n.textContent=text;return n;};
 const uid=()=>crypto.randomUUID();
 const remote=!!window.ReaderWorkspace;
 const base=window.ReaderLibraryEndpoint||(remote?'/workspace-api/':'http://127.0.0.1:8840/api/');
 const home=remote?location.origin+'/':new URL('../',base).href;
 const api=()=>window.ReaderWorkspaceEditor;
 let anchor=null,recognition=null,session=0,strokes=[],stroke=null,busy=false,conflicts=new Map();
 const fields=['type','source','level','blockStyle','formatting'];
 const sharedContent=b=>Object.fromEntries(fields.filter(k=>b[k]!==undefined).map(k=>[k,b[k]]));
 const signature=b=>JSON.stringify(sharedContent(b));
 async function request(path,body){
  const r=await fetch(base+path,body?{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)}:{cache:'no-store'});
  const value=await r.json();if(!r.ok){const e=Error(value.error||'Library request failed.');e.status=r.status;throw e;}return value;
 }
 function library(){
  if(parent!==window){let origin;try{origin=new URL(document.referrer).origin;}catch(_){}
   if(origin){parent.postMessage({type:'reader:library'},origin);return;}}
  // Keep this editor and its audio alive in the original tab.
  window.open(home,'_blank','noopener');
 }
 const dialog=el('dialog');dialog.id='workspace-insert';dialog.setAttribute('aria-label','Insert into this file');
 dialog.innerHTML=`<header><strong>Insert</strong><button type="button" data-close aria-label="Close insert panel">Close</button></header>
 <nav aria-label="Insertion tools"><button data-tab="dictation">Dictation</button><button data-tab="ink">Handwriting</button><button data-tab="shared">Shared blocks</button><button data-tab="file">New file</button></nav>
 <p id="workspace-message" role="status"></p>
 <section data-panel="dictation"><label>Language <select id="dictation-language"><option value="zh-CN">Mandarin</option><option value="en-US">English</option></select></label>
 <p class="workspace-help">Uses your browser's speech recognition service. Start requests microphone access. You can also use keyboard dictation here.</p>
 <button id="dictation-start">Start microphone</button><button id="dictation-stop" disabled>Stop</button>
 <label><input type="checkbox" id="dictation-live"> Insert final phrases live</label>
 <label for="dictation-text">Transcript / formula source</label><textarea id="dictation-text" rows="5" placeholder="Speak or type. Review before inserting."></textarea><p id="dictation-interim" aria-live="polite"></p>
 <label>Insert as <select id="dictation-kind"><option value="paragraph">Text</option><option value="equation">Equation (LaTeX)</option></select></label>
 <button id="dictation-insert">Insert transcript</button></section>
 <section data-panel="ink" hidden><p class="workspace-help">Write or draw with Apple Pencil. Strokes remain editable; this does not convert handwriting to LaTeX.</p>
 <div class="ink-actions"><label>Pen <select id="ink-width"><option value="2">Fine</option><option value="4" selected>Medium</option><option value="7">Bold</option></select></label><label><input type="checkbox" id="ink-finger"> Draw with finger</label><button id="ink-undo">Undo stroke</button><button id="ink-clear">Clear</button></div>
 <canvas id="workspace-ink" width="1200" height="600" aria-label="Handwriting canvas"></canvas><button id="ink-save">Insert drawing</button></section>
 <section data-panel="shared" hidden><p class="workspace-help">Shared text and formulas update across files while they are open. Closed files catch up when reopened. Conflicting edits are kept for review.</p>
 <button id="shared-publish">Share selected block</button><button id="shared-unlink">Unlink selected block</button><button id="shared-refresh">Refresh list</button><div id="shared-conflict" hidden><p>Both this file and the shared block changed.</p><button id="shared-use-remote">Use shared version</button><button id="shared-use-local">Publish my version</button></div><div id="shared-list"></div></section>
 <section data-panel="file" hidden><label for="workspace-file-title">Title</label><input id="workspace-file-title" maxlength="160" placeholder="Untitled"><label>Language <select id="workspace-file-language"><option value="en-US">English</option><option value="zh-CN">Mandarin</option></select></label><button id="workspace-create">Create file in Library</button><a id="workspace-created" hidden target="_blank" rel="noopener">Open new file</a></section>`;
 document.body.append(dialog);
 const message=t=>$('#workspace-message').textContent=t;
 const run=fn=>async()=>{try{await fn();}catch(e){message(e.message);}};
 function stop(){session++;if(recognition){recognition.onend=null;recognition.abort();recognition=null;}$('#dictation-start').disabled=false;$('#dictation-stop').disabled=true;$('#dictation-interim').textContent='';}
 function tab(name){dialog.querySelectorAll('[data-panel]').forEach(n=>n.hidden=n.dataset.panel!==name);dialog.querySelectorAll('[data-tab]').forEach(n=>n.setAttribute('aria-pressed',String(n.dataset.tab===name)));message('');if(name==='shared')void refreshShared().catch(e=>message(e.message));}
 function open(name='dictation'){
  api()?.flush();anchor=api()?.selected();if(!api()?.block(anchor)){message('Select a block first.');return;}
  if(!dialog.open)dialog.showModal();tab(name);
  if(name==='ink'){
   const b=api().block(anchor);strokes=b.ink?.version===1?JSON.parse(JSON.stringify(b.ink.strokes)):[];
   dialog.dataset.inkEdit=b.ink?.version===1?anchor:'';$('#ink-save').textContent=dialog.dataset.inkEdit?'Save drawing':'Insert drawing';draw();
  }
 }
 dialog.querySelector('[data-close]').onclick=()=>dialog.close();dialog.addEventListener('close',stop);dialog.addEventListener('cancel',stop);
 dialog.querySelectorAll('[data-tab]').forEach(n=>n.onclick=()=>{stop();if(n.dataset.tab==='ink')open('ink');else tab(n.dataset.tab);});
 const bar=$('header.topbar')||document.body;
 const back=el('button','Library');back.id='reader-library-home';back.className='quiet';back.onclick=library;bar.prepend(back);
 const trigger=el('button','Insert');trigger.id='reader-insert-tools';trigger.className='quiet';trigger.onclick=()=>open();bar.append(trigger);
 function insert(source,type='paragraph',extra={}){if(!source.trim())throw Error('There is no content to insert.');const id=api().insert(anchor,{id:'local-'+uid(),type,source,...extra});anchor=id;return id;}
 function draft(){try{localStorage.setItem('workspace:dictation:'+api().documentId(),$('#dictation-text').value);}catch(_){}}
 $('#dictation-text').addEventListener('input',draft);
 trigger.addEventListener('click',()=>{try{$('#dictation-text').value=localStorage.getItem('workspace:dictation:'+api().documentId())||$('#dictation-text').value;}catch(_){}});
 $('#dictation-insert').onclick=run(()=>{insert($('#dictation-text').value,$('#dictation-kind').value);$('#dictation-text').value='';draft();message('Inserted. Autosaving in this file.');});
 $('#dictation-stop').onclick=stop;
 $('#dictation-start').onclick=run(()=>{
  const Recognizer=window.SpeechRecognition||window.webkitSpeechRecognition;
  if(!Recognizer)throw Error('Live speech recognition is unavailable in this browser. Use Safari dictation on the keyboard, or type in the transcript.');
  stop();const token=++session;recognition=new Recognizer();recognition.lang=$('#dictation-language').value;recognition.continuous=true;recognition.interimResults=true;
  recognition.onresult=e=>{if(token!==session)return;let interim='';
   for(let i=e.resultIndex;i<e.results.length;i++){const result=e.results[i],text=result[0].transcript;if(result.isFinal){
    if($('#dictation-live').checked){try{insert(text,'paragraph');message('Final phrase inserted. Listening…');}catch(err){$('#dictation-text').value+=text;message(err.message);}}
    else $('#dictation-text').value+=($('#dictation-text').value?' ':'')+text;
   }else interim+=text;}
   $('#dictation-interim').textContent=interim;draft();};
  recognition.onerror=e=>{message('Microphone: '+e.error+'. Your transcript is retained.');stop();};
  recognition.onend=()=>{if(token!==session)return;stop();message('Listening stopped. Review the transcript or start again.');};
  recognition.start();$('#dictation-start').disabled=true;$('#dictation-stop').disabled=false;message('Listening…');
 });
 const canvas=$('#workspace-ink'),ctx=canvas.getContext('2d');
 function draw(){ctx.clearRect(0,0,1200,600);ctx.fillStyle='#ffffff';ctx.fillRect(0,0,1200,600);ctx.strokeStyle='#20252b';ctx.fillStyle='#20252b';ctx.lineCap='round';ctx.lineJoin='round';
  for(const s of strokes){for(let i=0;i<s.points.length;i++){const p=s.points[i],prev=s.points[Math.max(0,i-1)];ctx.lineWidth=s.width*(0.5+p[2]);ctx.beginPath();ctx.moveTo(prev[0],prev[1]);ctx.lineTo(p[0]+.01,p[1]);ctx.stroke();}}}
 function point(e){const r=canvas.getBoundingClientRect();return [(e.clientX-r.left)*1200/r.width,(e.clientY-r.top)*600/r.height,e.pressure||0.5];}
 canvas.addEventListener('pointerdown',e=>{if(e.pointerType==='touch'&&!$('#ink-finger').checked)return;e.preventDefault();canvas.setPointerCapture(e.pointerId);stroke={pointer:e.pointerId,width:Number($('#ink-width').value),points:[point(e)]};strokes.push(stroke);draw();});
 canvas.addEventListener('pointermove',e=>{if(!stroke||e.pointerId!==stroke.pointer)return;e.preventDefault();for(const p of e.getCoalescedEvents?.()||[e])stroke.points.push(point(p));draw();});
 const endInk=e=>{if(stroke?.pointer===e.pointerId)stroke=null;};canvas.addEventListener('pointerup',endInk);canvas.addEventListener('pointercancel',endInk);
 $('#ink-undo').onclick=()=>{strokes.pop();draw();};$('#ink-clear').onclick=()=>{strokes=[];draw();};
 $('#ink-save').onclick=run(async()=>{
  if(!strokes.length)throw Error('Draw something first.');const button=$('#ink-save');button.disabled=true;
  try{const savedStrokes=JSON.parse(JSON.stringify(strokes));const r=await fetch('/api/images',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({dataUrl:canvas.toDataURL('image/png')})});const data=await r.json();if(!r.ok)throw Error(data.error||'Drawing could not be saved.');
   const patch={type:'image',source:data.source,alt:'Handwritten notes',ink:{version:1,width:1200,height:600,strokes:savedStrokes}};
   if(dialog.dataset.inkEdit){api().patch(dialog.dataset.inkEdit,patch);message('Drawing updated.');}else{const id=insert(data.source,'image',patch);dialog.dataset.inkEdit=id;button.textContent='Save drawing';message('Drawing inserted. Select it and open Handwriting to edit strokes.');}
  }finally{button.disabled=false;}
 });
 let createId=null;
 $('#workspace-file-title').oninput=()=>createId=null;
 $('#workspace-create').onclick=run(async()=>{const button=$('#workspace-create');button.disabled=true;try{createId||=uid();const result=await request('create-file',{title:$('#workspace-file-title').value,language:$('#workspace-file-language').value,requestId:createId});const link=$('#workspace-created');link.href=home+'#reader='+result.id;link.hidden=false;message('Created '+result.title+'.');}finally{button.disabled=false;}});
 async function refreshShared(){const data=await request('shared');const list=$('#shared-list');list.replaceChildren();
  for(const row of data.blocks){const button=el('button',row.title+' · Insert synced block');button.onclick=run(async()=>{const value=await request('shared/get',{id:row.id});insert(value.content.source,value.content.type,{...value.content,librarySync:{id:value.id,revision:value.revision,base:JSON.stringify(value.content)}});message('Synced block inserted.');});list.append(button);}
  if(!data.blocks.length)list.textContent='Share a selected block to start.';
  $('#shared-conflict').hidden=!conflicts.has(anchor);
 }
 function applyRemote(id,value){const patch=Object.fromEntries(fields.map(k=>[k,value.content[k]??null]));patch.librarySync={id:value.id,revision:value.revision,base:JSON.stringify(value.content)};api().patch(id,patch);conflicts.delete(id);}
 $('#shared-publish').onclick=run(async()=>{api().flush();const b=api().block(anchor);if(!b)throw Error('Select a saved block first.');if(b.librarySync){message('This block is already shared.');return;}
  const start=signature(b);const value=await request('shared/create',{id:uid(),title:b.source.replace(/<[^>]*>/g,'').slice(0,80),content:sharedContent(b)});
  if(!api().block(anchor))throw Error('Selected block was removed. The shared copy is in the list.');
  api().patch(anchor,{librarySync:{id:value.id,revision:value.revision,base:JSON.stringify(value.content)}});message(start===signature(api().block(anchor))?'Block shared. Insert it from another file’s Shared blocks list.':'Block shared; newer local edits will sync next.');await refreshShared();});
 $('#shared-unlink').onclick=run(()=>{api().patch(anchor,{librarySync:null});conflicts.delete(anchor);message('This copy is independent.');$('#shared-conflict').hidden=true;});
 $('#shared-refresh').onclick=run(refreshShared);
 $('#shared-use-remote').onclick=run(async()=>{const b=api().block(anchor);const value=await request('shared/get',{id:b.librarySync.id});applyRemote(anchor,value);await refreshShared();message('Shared version applied. Undo is available.');});
 $('#shared-use-local').onclick=run(async()=>{const b=api().block(anchor),value=await request('shared/get',{id:b.librarySync.id}),snapshot=sharedContent(b);
  const result=await request('shared/update',{id:value.id,revision:value.revision,content:snapshot});api().patch(anchor,{librarySync:{id:result.id,revision:result.revision,base:JSON.stringify(result.content)}});conflicts.delete(anchor);await refreshShared();message('Local version published.');});
 async function sync(){if(busy||!api()?.ready()||api().busy()||document.hidden)return;busy=true;
  try{for(const initial of api().blocks().filter(b=>b.librarySync)){
   if(api().busy())break;const id=initial.id,binding={...initial.librarySync};if(conflicts.has(id))continue;const start=signature(initial);
   try{const value=await request('shared/get',{id:binding.id});const current=api().block(id);if(!current||signature(current)!==start||api().busy()||current.librarySync?.id!==binding.id)continue;
    const localChanged=start!==binding.base;
    if(value.revision!==binding.revision){if(localChanged&&JSON.stringify(value.content)!==start){conflicts.set(id,value);api().notice('Shared block conflict. Local edits retained; open Insert → Shared blocks.');continue;}applyRemote(id,value);}
    else if(localChanged){const result=await request('shared/update',{id:binding.id,revision:binding.revision,content:sharedContent(current)});const latest=api().block(id);if(latest?.librarySync?.id===binding.id)api().patch(id,{librarySync:{id:result.id,revision:result.revision,base:JSON.stringify(result.content)}});}
   }catch(e){if(e.status===409)conflicts.set(id,true);api().notice(e.message);break;}
  }}finally{busy=false;}}
 setInterval(()=>void sync(),2500);
 window.ReaderWorkspaceTools={open,library,sync,testing:{draw,getStrokes:()=>strokes,conflicts}};
})();
