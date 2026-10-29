/* Block handles, reversible reordering, and document-local synced content. */
(function(){'use strict';
 const SHARED=['type','level','blockStyle','source','image','alt','caption','formatting','speechOverride','figureCues','tableCues'];
 const clone=value=>value===undefined?undefined:JSON.parse(JSON.stringify(value));
 const equal=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
 const uid=()=> 'local-'+(crypto.randomUUID?.()||Date.now()+'-'+Math.random().toString(36).slice(2));
 const groupId=()=> 'sync-'+(crypto.randomUUID?.()||Date.now()+'-'+Math.random().toString(36).slice(2));
 const find=(doc,id)=>doc.blocks.find(b=>b.id===id);
 const peers=(doc,b)=>b?.syncGroup?doc.blocks.filter(x=>x.syncGroup===b.syncGroup):b?[b]:[];
 function assign(target,patch){for(const [key,value]of Object.entries(patch)){if(value==null)delete target[key];else target[key]=clone(value);}}
 function copyShared(source,target){let changed=false;for(const key of SHARED){if(equal(source[key],target[key]))continue;changed=true;if(source[key]===undefined)delete target[key];else target[key]=clone(source[key]);}return changed;}
 function syncFrom(doc,id){const source=find(doc,id);if(!source?.syncGroup)return [];const updated=[];for(const target of peers(doc,source))if(target!==source&&copyShared(source,target))updated.push(target.id);return updated;}
 function syncMany(doc,ids){const updated=new Set(ids),seen=new Set();for(const id of [...ids].reverse()){const b=find(doc,id);if(!b?.syncGroup||seen.has(b.syncGroup))continue;seen.add(b.syncGroup);for(const peer of syncFrom(doc,id))updated.add(peer);}return [...updated];}
 const isContainer=b=>Boolean(b&&(b.type==='callout'||b.type==='paragraph'&&b.blockStyle==='quote'));
 // A container owns a contiguous run in reading order; structural metadata is
 // deliberately excluded from SHARED, so linked content never reparents peers.
 function childList(doc,id){const at=doc.blocks.findIndex(b=>b.id===id),result=[];if(at<0)return result;for(let i=at+1;i<doc.blocks.length&&doc.blocks[i].parentId===id;i++)result.push(doc.blocks[i]);return result;}
 function emptyHeader(doc,b){return Boolean(isContainer(b)&&childList(doc,b.id).length&&!/<(?:img|video|audio|svg|math)\b/i.test(b.source||'')&&!String(b.source||'').replace(/<[^>]*>|&nbsp;|&#(?:160|x[aA]0);|[\s\u200b]/g,''));}
 function canMerge(previous,current){return Boolean(previous&&current&&!isContainer(current)&&((!isContainer(previous)&&(previous.parentId||null)===(current.parentId||null))||current.parentId===previous.id));}
 function normalizeStructure(doc){
  if(!doc?.blocks)return [];
  const changed=[];let parent=null;
  // Slash commands can turn a child into a new quote/callout without going
  // through setType. Lift that new container after its old siblings first.
  for(const b of [...doc.blocks])if(b.parentId&&isContainer(b)){const oldParent=find(doc,b.parentId);if(isContainer(oldParent)&&!oldParent.parentId&&childList(doc,oldParent.id).includes(b)){doc.blocks.splice(doc.blocks.indexOf(b),1);const at=doc.blocks.indexOf(oldParent)+childList(doc,oldParent.id).length+1;delete b.parentId;doc.blocks.splice(at,0,b);changed.push(b.id);}}
  for(const b of doc.blocks){
   if(b.parentId&&(!parent||b.parentId!==parent.id||isContainer(b))){delete b.parentId;changed.push(b.id);}
   if(!b.parentId)parent=isContainer(b)?b:null;
  }
  return changed;
 }
 function parentOf(doc,b){if(!b?.parentId)return null;const p=find(doc,b.parentId);return isContainer(p)&&!p.parentId&&childList(doc,p.id).includes(b)?p:null;}
 function insertionParent(doc,beforeId,afterId){
  const before=find(doc,beforeId),after=find(doc,afterId);
  if(before&&after&&doc.blocks.indexOf(after)!==doc.blocks.indexOf(before)+1)return null;
  const parent=parentOf(doc,before)||(isContainer(before)&&!before.parentId?before:null);
  return parent&&after?.parentId===parent.id?parent.id:null;
 }
 function packet(doc,id){const b=find(doc,id);return b?[b,...(!b.parentId&&isContainer(b)?childList(doc,id):[])]:[];}
 function copyPacket(doc,id,newId,newGroup,synced){
  const source=find(doc,id);if(!source)return null;const originals=packet(doc,id),copies=[];
  for(let i=0;i<originals.length;i++){
   const item=originals[i];if(synced)item.syncGroup=item.syncGroup||(i?groupId():newGroup);
   const next=clone(item);next.id=i?uid():newId;delete next.sourceBlockUrl;if(!synced){delete next.syncGroup;delete next.librarySync;}if(i)next.parentId=newId;copies.push(next);
  }
  doc.blocks.splice(doc.blocks.indexOf(originals.at(-1))+1,0,...copies);return copies[0];
 }
 function createCopy(doc,id,newId=uid(),newGroup=groupId()){return copyPacket(doc,id,newId,newGroup,true);}
 function duplicate(doc,id,newId=uid()){return copyPacket(doc,id,newId,null,false);}
 function unlink(doc,id){const b=find(doc,id);if(!b?.syncGroup&&!b?.librarySync)return false;delete b.syncGroup;delete b.librarySync;return true;}
 function move(doc,id,targetId,after=false){
  const source=find(doc,id),target=find(doc,targetId);if(!source||!target||id===targetId)return false;
  const moving=packet(doc,id);if(moving.includes(target))return false;
  const oldOrder=doc.blocks.map(b=>b.id).join('|'),oldParent=source.parentId;
  let anchor=target,parent=parentOf(doc,target);
  // Containers are never nested. Moving a whole group across a child crosses
  // that child's complete container instead of breaking the contiguous run.
  if(isContainer(source)&&parent){anchor=parent;parent=null;}
  for(const item of moving)doc.blocks.splice(doc.blocks.indexOf(item),1);
  let at=doc.blocks.indexOf(anchor);
  if(after)at+=packet(doc,anchor.id).length;
  if(parent&&!isContainer(source))source.parentId=parent.id;else delete source.parentId;
  doc.blocks.splice(at,0,...moving);normalizeStructure(doc);
  return oldOrder!==doc.blocks.map(b=>b.id).join('|')||oldParent!==source.parentId;
 }
 function moveInto(doc,id,parentId){
  const source=find(doc,id),parent=find(doc,parentId);
  if(!source||!parent||source===parent||isContainer(source)||!isContainer(parent)||parent.parentId)return false;
  const children=childList(doc,parent.id);if(source.parentId===parentId&&children.at(-1)===source)return false;
  doc.blocks.splice(doc.blocks.indexOf(source),1);const at=doc.blocks.indexOf(parent)+childList(doc,parentId).length+1;source.parentId=parentId;doc.blocks.splice(at,0,source);normalizeStructure(doc);return true;
 }
 function indentCandidate(doc,id){const b=find(doc,id),previous=doc.blocks[doc.blocks.indexOf(b)-1];if(!b||b.parentId||isContainer(b)||!previous)return null;return parentOf(doc,previous)||(isContainer(previous)&&!previous.parentId?previous:null);}
 function indent(doc,id){const parent=indentCandidate(doc,id);return parent?moveInto(doc,id,parent.id):false;}
 function outdent(doc,id){const b=find(doc,id),parent=parentOf(doc,b);if(!parent)return false;doc.blocks.splice(doc.blocks.indexOf(b),1);const at=doc.blocks.indexOf(parent)+childList(doc,parent.id).length+1;delete b.parentId;doc.blocks.splice(at,0,b);normalizeStructure(doc);return true;}
 function removeBlock(doc,id){const b=find(doc,id);if(!b)return false;for(const child of childList(doc,id))delete child.parentId;doc.blocks.splice(doc.blocks.indexOf(b),1);normalizeStructure(doc);return true;}
 function removeMany(doc,ids){
  const requested=new Set(ids),removed=doc.blocks.filter(b=>requested.has(b.id)).map(b=>b.id);
  if(!removed.length)return null;
  const first=doc.blocks.findIndex(b=>requested.has(b.id));
  doc.blocks=doc.blocks.filter(b=>!requested.has(b.id));
  let createdId=null;
  if(!doc.blocks.length){const next={id:uid(),type:'paragraph',source:''};doc.blocks.push(next);createdId=next.id;}
  normalizeStructure(doc);
  return {removedIds:removed,selectedId:doc.blocks[Math.min(first,doc.blocks.length-1)].id,createdId};
 }
 const TEXT_TYPES=new Set(['paragraph','heading','callout']);
 const typeValue=b=>b?.type==='heading'?'heading-'+Math.min(3,Math.max(1,Number(b.level)||2)):b?.type==='paragraph'&&b.blockStyle==='quote'?'quote':b?.type;
 function setType(doc,id,value){
  const b=find(doc,id);if(!b||!TEXT_TYPES.has(b.type)||!['paragraph','quote','callout','heading-1','heading-2','heading-3'].includes(value)||typeValue(b)===value)return false;
  delete b.blockStyle;delete b.level;
  if(value==='quote'){b.type='paragraph';b.blockStyle='quote';}else if(value==='paragraph'||value==='callout')b.type=value;else{b.type='heading';b.level=Number(value.slice(-1));}
  if(isContainer(b)&&b.parentId)outdent(doc,id);normalizeStructure(doc);return true;
 }
 function split(doc,id,patch,next){
  const source=find(doc,id);if(!source)return {changed:[]};
  const instances=peers(doc,source),continuationGroup=source.syncGroup?groupId():null;
  assign(source,patch);delete source.speechOverride;
  const changed=[source.id,...syncFrom(doc,id)];
  for(const current of instances){const continuation=clone(next);if(current.id!==id)continuation.id=uid();delete continuation.sourceBlockUrl;delete continuation.syncGroup;delete continuation.librarySync;if(continuationGroup)continuation.syncGroup=continuationGroup;if(isContainer(current))continuation.parentId=current.id;else if(current.parentId)continuation.parentId=current.parentId;doc.blocks.splice(doc.blocks.indexOf(current)+1,0,continuation);changed.push(continuation.id);}
  changed.push(...normalizeStructure(doc));return {changed};
 }
 function merge(doc,previousId,currentId,patch){
  const previous=find(doc,previousId),current=find(doc,currentId);if(!canMerge(previous,current))return {changed:[],detached:false};
  let pairs=[[previous,current]],detached=false;
  if(previous.syncGroup||current.syncGroup){
   const left=peers(doc,previous),right=peers(doc,current);
   const matching=previous.syncGroup&&current.syncGroup&&previous.syncGroup!==current.syncGroup&&left.length===right.length&&left.every(b=>{const next=doc.blocks[doc.blocks.indexOf(b)+1];return next?.syncGroup===current.syncGroup&&canMerge(b,next);});
   if(matching)pairs=left.map(b=>[b,doc.blocks[doc.blocks.indexOf(b)+1]]);
   else {detached=Boolean(previous.syncGroup||current.syncGroup);delete previous.syncGroup;delete current.syncGroup;}
  }
  for(const [left,right]of pairs){if(left.librarySync||right.librarySync){delete left.librarySync;delete right.librarySync;detached=true;}}
  assign(previous,patch);delete previous.speechOverride;
  for(const [left,right]of pairs){for(const child of childList(doc,right.id)){if(isContainer(left)&&!left.parentId)child.parentId=left.id;else delete child.parentId;}}
  const removed=new Set(pairs.map(pair=>pair[1].id));doc.blocks=doc.blocks.filter(b=>!removed.has(b.id));
  const changed=[previous.id,...syncFrom(doc,previous.id),...removed,...normalizeStructure(doc)];
  return {changed,detached};
 }
 let api={},menu=null,menuId=null,dragId=null,dropTarget=null,pointerDrag=null,suppressClickUntil=0,initialized=false,scrollFrame=0;
 const documentModel=()=>api.getDocument?.();
 const blockNode=id=>api.getNode?.(id);
 const button=(text,action,disabled=false)=>{const b=document.createElement('button');b.type='button';b.textContent=text;b.disabled=disabled;b.setAttribute('role','menuitem');b.onclick=()=>{closeMenu();action();};return b;};
 function closeMenu({focus=false}={}){const id=menuId;menu?.remove();menu=null;menuId=null;document.querySelectorAll('.block-handle[aria-expanded="true"]').forEach(n=>n.setAttribute('aria-expanded','false'));if(focus)blockNode(id)?.querySelector('.block-handle')?.focus({preventScroll:true});}
 function before(){window.ReaderDirectEdit?.flush();api.beforeChange?.();}
 function finish(id,{focus=false,message=''}={}){normalizeStructure(documentModel());api.changed?.(id,{structural:true,presentation:true});api.afterChange?.();if(focus)blockNode(id)?.querySelector('.block-handle')?.focus({preventScroll:true});if(message)api.toast?.(message);}
 function reorder(id,targetId,after=false,inside=false){const doc=documentModel();if(!doc||!find(doc,id)||!find(doc,targetId)||id===targetId)return false;before();if(inside?moveInto(doc,id,targetId):move(doc,id,targetId,after)){api.select?.(id);finish(id,{focus:true});return true;}return false;}
 function changeNesting(id,remove=false){const doc=documentModel(),b=find(doc,id);if(!b||!(remove?parentOf(doc,b):indentCandidate(doc,id)))return false;before();const changed=remove?outdent(doc,id):indent(doc,id);if(changed){api.select?.(id);finish(id,{focus:true,message:remove?'Moved out of container.':'Moved into the preceding quote or callout.'});}return changed;}
 function moveBy(id,delta){const doc=documentModel();if(!doc||!delta)return false;const index=doc.blocks.findIndex(b=>b.id===id);if(index<0)return false;const target=doc.blocks[index+(delta>0?packet(doc,id).length:-1)];return target?reorder(id,target.id,delta>0):false;}
 function duplicateBlock(id){const doc=documentModel();if(!doc||!find(doc,id))return null;before();const next=duplicate(doc,id);api.select?.(next.id);finish(next.id,{focus:true});return next.id;}
 function deleteBlock(id){const doc=documentModel(),b=doc&&find(doc,id);if(!b||doc.blocks.length<2)return false;before();const at=doc.blocks.indexOf(b);removeBlock(doc,id);const selected=doc.blocks[Math.min(at,doc.blocks.length-1)]?.id;api.select?.(selected);finish(id,{message:isContainer(b)?'Container deleted; its contents were kept. Undo restores the container.':'Block deleted. Undo restores it.'});return selected;}
 function deleteBlocks(ids){
  const doc=documentModel();if(!doc||!ids.some(id=>find(doc,id)))return false;
  before();const result=removeMany(doc,ids);if(!result)return false;
  api.select?.(result.selectedId);finish(result.removedIds[0],{message:result.removedIds.length+' block'+(result.removedIds.length===1?'':'s')+' deleted. Undo restores them.'});
  return result;
 }
 function menuFor(id,handle){
  if(menuId===id){closeMenu({focus:true});return;}closeMenu();const doc=documentModel(),b=find(doc,id);if(!b)return;
  api.select?.(id);menuId=id;menu=document.createElement('div');menu.className='block-action-menu';menu.setAttribute('role','menu');menu.setAttribute('aria-label','Block actions');
  const header=document.createElement('div');header.className='block-menu-heading';header.textContent=b.syncGroup?`${peers(doc,b).length} synced copies · this document`:'Block';menu.append(header);window.ReaderBlockSelection?.addMenuItems(menu,id,button);
  if(TEXT_TYPES.has(b.type)){
   const row=document.createElement('label');row.className='block-type-picker';const label=document.createElement('span');label.textContent='Turn into';
   const picker=document.createElement('select');picker.setAttribute('aria-label','Block type');
   for(const [value,text]of [['paragraph','Text'],['heading-1','Heading 1'],['heading-2','Heading 2'],['heading-3','Heading 3'],['quote','Quote'],['callout','Callout']]){const option=document.createElement('option');option.value=value;option.textContent=text;picker.append(option);}
   picker.value=typeValue(b);picker.onchange=()=>{const value=picker.value;closeMenu();before();if(setType(doc,id,value))finish(id,{focus:true});};row.append(label,picker);menu.append(row);
  }
  if(parentOf(doc,b))menu.append(button('Move out of container',()=>changeNesting(id,true)));
  else if(indentCandidate(doc,id))menu.append(button('Move into preceding quote or callout',()=>changeNesting(id)));
  menu.append(button('Move up',()=>moveBy(id,-1),doc.blocks[0]===b),button('Move down',()=>moveBy(id,1),doc.blocks.at(-1)===packet(doc,id).at(-1)));
  menu.append(button('Create synced copy below',()=>{before();const copy=createCopy(doc,id);api.select?.(copy.id);finish(copy.id,{focus:true,message:'Synced copies share text, formulas and colors.'});}));
  menu.append(button('Share across Library files…',()=>{api.select?.(id);window.ReaderWorkspaceTools?.open('shared');}));
  menu.append(button('Dictate below…',()=>{api.select?.(id);window.ReaderWorkspaceTools?.open('dictation');}));
  menu.append(button(b.ink?'Edit handwriting…':'Draw below…',()=>{api.select?.(id);window.ReaderWorkspaceTools?.open('ink');}));
  if(b.syncGroup)menu.append(button('Unlink this copy',()=>{before();unlink(doc,id);finish(id,{focus:true,message:'This block is now independent.'});}));
  menu.append(button('Duplicate independently',()=>duplicateBlock(id)));
  window.ReaderImages?.addMenuItems(menu,id,button);
  if(emptyHeader(doc,b))menu.append(button('Add text at start',()=>window.ReaderInsertion?.insert(b.id,childList(doc,b.id)[0].id)));
  menu.append(button('Source and narration…',()=>api.edit?.(id)));
  if(b.sourceBlockUrl){const link=document.createElement('a');link.textContent='Open original source ↗';link.href=b.sourceBlockUrl;link.target='_blank';link.rel='noopener';link.setAttribute('role','menuitem');link.onclick=()=>closeMenu();menu.append(link);}
  const remove=button(isContainer(b)?'Remove container, keep contents':'Delete this block',()=>deleteBlock(id),doc.blocks.length<2);remove.className='danger';menu.append(remove);
  document.body.append(menu);handle.setAttribute('aria-expanded','true');const rect=handle.getBoundingClientRect();menu.style.left=Math.max(8,Math.min(rect.left,innerWidth-menu.offsetWidth-8))+'px';menu.style.top=Math.max(8,Math.min(rect.bottom+5,innerHeight-menu.offsetHeight-8))+'px';menu.querySelector('button:not(:disabled)')?.focus({preventScroll:true});
 }
 function clearDrop(){dropTarget?.classList.remove('drop-before','drop-after','drop-inside');dropTarget=null;document.querySelectorAll('.block-dragging').forEach(n=>n.classList.remove('block-dragging'));}
 function startPointer(e,id,source){
  if(e.button!==0||e.isPrimary===false||pointerDrag)return;
  pointerDrag={id,x:e.clientX,y:e.clientY,pointerId:e.pointerId,handle:source,active:false,lastX:e.clientX,lastY:e.clientY};
  source.setPointerCapture?.(e.pointerId);
  // A browser-native image drag inserts HTML into a contenteditable destination.
  // Keep the gesture entirely within the existing reversible block move path.
  if(source.matches('img,.block-drag-equation'))e.preventDefault();
 }
 function dragTarget(x,y){
  const top=document.querySelector('.topbar')?.getBoundingClientRect().bottom||0,bottom=document.querySelector('.player')?.getBoundingClientRect().top||innerHeight;
  if(y<top||y>bottom)return null;
  const target=document.elementFromPoint(x,y)?.closest('.source-block');
  if(target)return target;
  const holder=document.querySelector('#blocks'),bounds=holder?.getBoundingClientRect();
  if(!bounds||x<bounds.left-55||x>bounds.right+30||y<bounds.top-20||y>bounds.bottom+20)return null;
  let nearest=null,distance=Infinity;
  for(const block of holder.querySelectorAll(':scope > .source-block')){
   const rect=block.getBoundingClientRect(),gap=y<rect.top?rect.top-y:y>rect.bottom?y-rect.bottom:0;
   if(gap<distance){nearest=block;distance=gap;}
  }
  return nearest;
 }
 function markDrop(x,y){
  const target=dragTarget(x,y);
  if(dropTarget&&dropTarget!==target)dropTarget.classList.remove('drop-before','drop-after','drop-inside');dropTarget=null;
  if(target&&target.dataset.blockId!==(pointerDrag?.id||dragId)){
   dropTarget=target;const r=target.getBoundingClientRect(),b=find(documentModel(),target.dataset.blockId),source=find(documentModel(),pointerDrag?.id||dragId);
   // The gutter belongs to the same drop lane as the block: a straight
   // vertical handle drag must reach the callout, not a no-op reorder.
   const inside=isContainer(b)&&!b.parentId&&!isContainer(source)&&y>r.top+Math.min(10,r.height*.15)&&y<r.bottom-Math.min(10,r.height*.15)&&x>=r.left-55&&x<=r.right+12;
   target.classList.toggle('drop-inside',inside);target.classList.toggle('drop-before',!inside&&y<=r.top+r.height/2);target.classList.toggle('drop-after',!inside&&y>r.top+r.height/2);
  }
 }
 function scrollDrag(){
  scrollFrame=0;const state=pointerDrag;if(!state?.active)return;
  const top=document.querySelector('.topbar')?.getBoundingClientRect().bottom||0,bottom=document.querySelector('.player')?.getBoundingClientRect().top||innerHeight;
  const speed=state.lastY<top+45?-14:state.lastY>bottom-45?14:0;
  if(speed){window.scrollBy(0,speed);markDrop(state.lastX,state.lastY);}
  scrollFrame=requestAnimationFrame(scrollDrag);
 }
 function updateLabel(b,wrap,position){const handle=wrap.querySelector('.block-handle');if(handle){handle.setAttribute('aria-label',`Move or manage passage ${position+1}`);handle.title='Drag to move · click for actions · Tab to move into a quote/callout · Shift+Tab to move out · Option ↑/↓ to reorder';}const label=wrap.querySelector('.synced-block-badge');if(label)label.title=`${peers(documentModel(),b).length} linked copies in this document; edits update all copies.`;}
 function refreshStructure(){
  const doc=documentModel();if(!doc)return;
  for(const b of doc.blocks){
   const wrap=blockNode(b.id);if(!wrap)continue;const parent=parentOf(doc,b),children=!b.parentId&&isContainer(b)?childList(doc,b.id):[];
   for(const name of ['block-container','has-block-children','block-child','container-quote','container-callout','container-last-child'])wrap.classList.remove(name);
   delete wrap.dataset.parentId;
   wrap.classList.toggle('container-empty-header',emptyHeader(doc,b));
   wrap.classList.toggle('container-first-child',Boolean(parent&&emptyHeader(doc,parent)&&childList(doc,parent.id)[0]===b));
   if(isContainer(b)&&!b.parentId){wrap.classList.add('block-container',b.type==='callout'?'container-callout':'container-quote');wrap.classList.toggle('has-block-children',children.length>0);}
   if(parent){wrap.classList.add('block-child',parent.type==='callout'?'container-callout':'container-quote');wrap.classList.toggle('container-last-child',childList(doc,parent.id).at(-1)===b);wrap.dataset.parentId=parent.id;}
  }
 }
 function bind(b,wrap){
  if(!documentModel()?.blocks.includes(b))return;
  const handle=document.createElement('button');handle.type='button';handle.className='block-handle';handle.draggable=false;handle.textContent='⠿';handle.setAttribute('aria-haspopup','menu');handle.setAttribute('aria-expanded','false');handle.setAttribute('data-direct-ignore','');handle.onclick=e=>{e.stopPropagation();if(performance.now()<suppressClickUntil)return;menuFor(b.id,handle);};
  handle.onpointerdown=e=>startPointer(e,b.id,handle);
  handle.onkeydown=e=>{if(e.altKey&&['ArrowUp','ArrowDown'].includes(e.key)){e.preventDefault();e.stopPropagation();moveBy(b.id,e.key==='ArrowUp'?-1:1);}else if(e.key==='Escape'){closeMenu({focus:true});}};
  handle.ondragstart=e=>{window.ReaderDirectEdit?.flush();closeMenu();dragId=b.id;e.dataTransfer.effectAllowed='move';e.dataTransfer.setData('text/plain',b.id);wrap.classList.add('block-dragging');};
  handle.ondragend=()=>{dragId=null;clearDrop();};wrap.append(handle);
  if(b.syncGroup||b.librarySync){wrap.classList.add('synced-block');const badge=document.createElement('button');badge.type='button';badge.className='synced-block-badge';badge.textContent='↔';badge.setAttribute('aria-label','Synced block actions');badge.setAttribute('data-direct-ignore','');badge.onclick=e=>{e.stopPropagation();menuFor(b.id,handle);};wrap.append(badge);}
  wrap.ondragover=e=>{if(!dragId||dragId===b.id)return;e.preventDefault();e.dataTransfer.dropEffect='move';markDrop(e.clientX,e.clientY);if(e.clientY<90)window.scrollBy(0,-18);else if(e.clientY>innerHeight-100)window.scrollBy(0,18);};
  wrap.ondrop=e=>{if(!dragId)return;e.preventDefault();e.stopPropagation();const from=dragId,after=wrap.classList.contains('drop-after'),inside=wrap.classList.contains('drop-inside');dragId=null;clearDrop();if(from!==b.id)reorder(from,b.id,after,inside);};
  for(const image of wrap.querySelectorAll('.block-content img')){
   image.draggable=false;image.contentEditable='false';image.classList.add('block-drag-image');image.tabIndex=0;
   image.title=b.type==='image'?'Drag to move this image block':'Drag to move this passage with its image';
   image.addEventListener('pointerdown',e=>startPointer(e,b.id,image));
   image.addEventListener('dragstart',e=>{e.preventDefault();e.stopImmediatePropagation();});
   image.addEventListener('keydown',e=>{if(e.altKey&&['ArrowUp','ArrowDown'].includes(e.key)){e.preventDefault();e.stopPropagation();moveBy(b.id,e.key==='ArrowUp'?-1:1);}});
  }
  if(b.type==='equation')for(const equation of wrap.querySelectorAll('.block-content > .math-display')){
   equation.draggable=false;equation.classList.add('block-drag-equation');
   equation.title='Drag to move this equation · click to edit';
   equation.addEventListener('pointerdown',e=>startPointer(e,b.id,equation));
   equation.addEventListener('dragstart',e=>{e.preventDefault();e.stopImmediatePropagation();});
  }
  // A completed/cancelled drag must not also open the equation editor.
  wrap.addEventListener('click',e=>{if(e.target.closest?.('.block-drag-image,.block-drag-equation')&&performance.now()<suppressClickUntil){e.preventDefault();e.stopImmediatePropagation();}},true);
  updateLabel(b,wrap,documentModel().blocks.indexOf(b));
 }
 function pointerMove(e){
  const state=pointerDrag;if(!state||e.pointerId!==state.pointerId)return;
  if(!state.active&&Math.hypot(e.clientX-state.x,e.clientY-state.y)<5)return;
  if(!state.active){state.active=true;window.ReaderDirectEdit?.flush();closeMenu();dragId=state.id;blockNode(state.id)?.classList.add('block-dragging');document.body.classList.add('block-pointer-dragging');scrollFrame=requestAnimationFrame(scrollDrag);}
  e.preventDefault();state.lastX=e.clientX;state.lastY=e.clientY;markDrop(e.clientX,e.clientY);
 }
 function pointerEnd(e,cancel=false){
  const state=pointerDrag;if(!state||e.pointerId!==state.pointerId)return;
  const target=dropTarget?.dataset.blockId,after=dropTarget?.classList.contains('drop-after'),inside=dropTarget?.classList.contains('drop-inside');pointerDrag=null;dragId=null;clearDrop();cancelAnimationFrame(scrollFrame);scrollFrame=0;document.body.classList.remove('block-pointer-dragging');
  if(state.handle.hasPointerCapture?.(state.pointerId))state.handle.releasePointerCapture(state.pointerId);
  if(state.active){suppressClickUntil=performance.now()+300;e.preventDefault();if(!cancel&&target)reorder(state.id,target,after,inside);}
 }
 function cancelPointer(){if(pointerDrag)pointerEnd({pointerId:pointerDrag.pointerId,preventDefault(){}},true);}
 function init(options){api=options;if(initialized)return;initialized=true;document.addEventListener('pointermove',pointerMove,{passive:false});document.addEventListener('pointerup',e=>pointerEnd(e));document.addEventListener('pointercancel',e=>pointerEnd(e,true));
  document.addEventListener('lostpointercapture',e=>{if(pointerDrag?.handle===e.target)pointerEnd(e,true);});
  document.addEventListener('beforeinput',e=>{if((dragId||pointerDrag?.active)&&['insertFromDrop','deleteByDrag'].includes(e.inputType)){e.preventDefault();e.stopImmediatePropagation();}},true);
  document.addEventListener('drop',e=>{if(pointerDrag){e.preventDefault();e.stopImmediatePropagation();cancelPointer();}},true);
  window.addEventListener('blur',cancelPointer);document.addEventListener('visibilitychange',()=>{if(document.hidden)cancelPointer();});
  document.addEventListener('pointerdown',e=>{if(menu&&!menu.contains(e.target)&&!e.target.closest('.block-handle,.synced-block-badge'))closeMenu();});document.addEventListener('keydown',e=>{if(e.key==='Escape'&&pointerDrag){e.preventDefault();e.stopImmediatePropagation();cancelPointer();return;}if(!menu){if(e.key==='Tab'&&!e.metaKey&&!e.ctrlKey&&!e.altKey&&!e.isComposing&&e.target.matches?.('.block-handle,.direct-math-atom,.block-drag-image,.source-block')&&(!e.target.matches('.direct-math-atom')||e.target.closest('.source-block')?.classList.contains('type-equation'))){const id=e.target.closest('.source-block')?.dataset.blockId;if(changeNesting(id,e.shiftKey)){e.preventDefault();e.stopImmediatePropagation();}}return;}if(e.target.closest?.('select')&&e.key!=='Escape')return;if(e.key==='Escape'){e.preventDefault();e.stopImmediatePropagation();closeMenu({focus:true});}else if(['ArrowDown','ArrowUp','Home','End'].includes(e.key)){e.preventDefault();const items=[...menu.querySelectorAll('[role="menuitem"]')].filter(n=>!n.disabled),at=items.indexOf(document.activeElement);items[e.key==='Home'?0:e.key==='End'?items.length-1:(at+(e.key==='ArrowDown'?1:-1)+items.length)%items.length]?.focus();}},true);window.addEventListener('resize',()=>closeMenu());window.addEventListener('scroll',()=>closeMenu(),{passive:true});}
 window.ReaderBlocks={deleteBlocks,init,bind,updateLabel,refreshStructure,moveBy,duplicateBlock,deleteBlock,normalizeStructure:()=>normalizeStructure(documentModel()),insertionParent:(before,after)=>insertionParent(documentModel(),before,after),indent:id=>changeNesting(id),outdent:id=>changeNesting(id,true),syncFrom:id=>syncFrom(documentModel(),id),syncMany:ids=>syncMany(documentModel(),ids),split:(id,patch,next)=>split(documentModel(),id,patch,next),merge:(left,right,patch)=>merge(documentModel(),left,right,patch),testing:{removeMany,isContainer,childList,normalizeStructure,parentOf,insertionParent,packet,moveInto,indentCandidate,indent,outdent,removeBlock,SHARED,copyShared,syncFrom,syncMany,createCopy,duplicate,unlink,move,setType,typeValue,split,merge}};
 Object.assign(window.ReaderBlocks,{isContainer,emptyHeader:id=>emptyHeader(documentModel(),find(documentModel(),id)),children:id=>childList(documentModel(),id),canMerge:(left,right)=>canMerge(find(documentModel(),left),find(documentModel(),right))});
})();
