/* Explicit block selection is separate from text selection and narration. */
(function(){'use strict';
 let api={},initialized=false,toolbar=null,countLabel=null,removeButton=null,selected=new Set(),anchorId=null,focusId=null,documentId=null,recentDeletion=false;
 const model=()=>api.getDocument?.();
 const orderedIds=()=>model()?.blocks?.map(b=>b.id)||[];
 const node=id=>api.getNode?.(id);
 const picker=id=>node(id)?.querySelector(':scope > .block-select-toggle');
 const isTyping=target=>Boolean(target?.closest?.('input,textarea,select,[contenteditable="true"],[contenteditable="plaintext-only"],.direct-math-editor,[role="dialog"]'));
 const isControl=target=>Boolean(target?.closest?.('.block-select-toggle,.block-selection-toolbar,.block-action-menu,.block-handle,.synced-block-badge'));
 const rangeIds=(ids,a,b)=>{const left=ids.indexOf(a),right=ids.indexOf(b);return left<0||right<0?[]:ids.slice(Math.min(left,right),Math.max(left,right)+1);};
 const validSelected=()=>orderedIds().filter(id=>selected.has(id));
 function focusPicker(id,{reveal=false}={}){if(window.ReaderPageNavigation?.ensureBlock?.(id)===false)return false;const button=picker(id);if(!button)return false;button.focus({preventScroll:true});if(reveal)button.scrollIntoView({block:'nearest'});return true;}
 function refresh(){
  const doc=model(),ids=doc?.blocks?.map(b=>b.id)||[],valid=new Set(ids);
  if(documentId!=null&&doc?.id!==documentId){selected.clear();anchorId=null;focusId=null;recentDeletion=false;}
  documentId=doc?.id??null;
  selected=new Set([...selected].filter(id=>valid.has(id)));
  if(anchorId&&!valid.has(anchorId))anchorId=validSelected()[0]||null;
  if(focusId&&!valid.has(focusId))focusId=validSelected().at(-1)||null;
  for(let i=0;i<ids.length;i++){
   const id=ids[i],wrap=node(id);if(!wrap)continue;const picked=selected.has(id),toggle=wrap.querySelector(':scope > .block-select-toggle');
   wrap.classList.toggle('block-picked',picked);wrap.toggleAttribute('data-block-selected',picked);
   if(toggle){toggle.setAttribute('aria-checked',String(picked));toggle.setAttribute('aria-label',`${picked?'Deselect':'Select'} passage ${i+1}`);toggle.title='Select block · Shift-click for a range';}
  }
  if(toolbar){const count=selected.size;toolbar.hidden=count===0;countLabel.textContent=`${count} ${count===1?'block':'blocks'} selected`;removeButton.disabled=count===0;}
  document.body.classList.toggle('block-selection-active',selected.size>0);
 }
 function clear({focus=false}={}){const last=focusId;selected.clear();anchorId=null;focusId=null;refresh();if(focus)focusPicker(last);}
 function select(id,{shift=false,additive=false,toggle=false,focus=true}={}){
  const ids=orderedIds();if(!ids.includes(id))return false;
  window.ReaderDirectEdit?.flush();window.ReaderFormatting?.clearSelection?.();window.getSelection()?.removeAllRanges();recentDeletion=false;
  if(shift&&anchorId&&ids.includes(anchorId)){
   const range=rangeIds(ids,anchorId,id);selected=additive?new Set([...selected,...range]):new Set(range);
  }else{
   if(toggle||additive){if(selected.has(id))selected.delete(id);else selected.add(id);}
   else selected=new Set([id]);
   anchorId=id;
  }
  focusId=id;api.select?.(id);refresh();if(focus)focusPicker(id);return true;
 }
 function selectAll(){const ids=orderedIds();if(!ids.length)return;selected=new Set(ids);anchorId=ids[0];focusId=focusId&&selected.has(focusId)?focusId:ids.at(-1);refresh();}
 function removeSelected(){
  const ids=validSelected();if(!ids.length)return false;
  window.ReaderDirectEdit?.flush();const result=api.deleteBlocks?.(ids);
  if(result===false||result===undefined){refresh();return false;}
  const survivors=new Set(orderedIds());
  if(ids.every(id=>survivors.has(id))){refresh();return false;}
  selected.clear();anchorId=null;focusId=null;recentDeletion=true;refresh();
  const next=typeof result==='string'?result:result?.selectedId||api.getSelectedId?.();if(next)focusPicker(next);
  return true;
 }
 function button(text,action,className=''){const el=document.createElement('button');el.type='button';el.textContent=text;el.className=className;el.onclick=action;return el;}
 function ensureToolbar(){
  if(toolbar)return;toolbar=document.createElement('div');toolbar.className='block-selection-toolbar';toolbar.hidden=true;toolbar.setAttribute('role','toolbar');toolbar.setAttribute('aria-label','Selected block actions');toolbar.dataset.directIgnore='';
  countLabel=document.createElement('span');countLabel.className='block-selection-count';countLabel.setAttribute('role','status');countLabel.setAttribute('aria-live','polite');countLabel.setAttribute('aria-atomic','true');
  removeButton=button('Delete',removeSelected,'block-selection-delete');removeButton.setAttribute('aria-label','Delete selected blocks');removeButton.title='Delete selected blocks · Backspace or Delete · Undo restores them';
  const all=button('Select all',selectAll),cancel=button('Clear',()=>clear({focus:true}));all.setAttribute('aria-label','Select all blocks');cancel.setAttribute('aria-label','Clear block selection');
  toolbar.append(countLabel,all,removeButton,cancel);document.body.append(toolbar);
 }
 function bind(block,wrap){
  if(document.body.classList.contains('floating')||!model()?.blocks?.some(b=>b.id===block.id)||wrap.querySelector(':scope > .block-select-toggle'))return;
  wrap.classList.add('block-pickable');const toggle=document.createElement('button');toggle.type='button';toggle.className='block-select-toggle';toggle.setAttribute('role','checkbox');toggle.setAttribute('aria-checked',String(selected.has(block.id)));toggle.setAttribute('aria-label','Select block');toggle.dataset.directIgnore='';toggle.title='Select block · Shift-click for a range';
  toggle.innerHTML='<svg viewBox="0 0 18 18" aria-hidden="true" focusable="false"><rect x="2.5" y="2.5" width="13" height="13" rx="2"></rect><path d="m5.5 9 2.3 2.5 4.8-5"></path></svg>';
  toggle.addEventListener('pointerdown',e=>{e.stopPropagation();});
  toggle.addEventListener('click',e=>{e.preventDefault();e.stopPropagation();select(block.id,{shift:e.shiftKey,additive:e.metaKey||e.ctrlKey,toggle:!e.shiftKey});});
  wrap.append(toggle);
 }
 function addMenuItems(menu,id,makeButton){menu.append(makeButton(selected.has(id)?'Deselect block':'Select block',()=>select(id,{toggle:selected.has(id)})));if(selected.size>0&&anchorId&&anchorId!==id)menu.append(makeButton('Select through this block',()=>select(id,{shift:true})));}
 function handleKeydown(e){
  if(e.defaultPrevented||e.isComposing||document.body.classList.contains('floating')||isTyping(e.target))return;
  const key=e.key.toLowerCase(),mod=e.metaKey||e.ctrlKey;
  if(e.target.closest?.('[role="menu"]'))return;
  if(mod&&key==='z'&&recentDeletion&&api.undo){e.preventDefault();e.stopImmediatePropagation();api.undo(e.shiftKey);recentDeletion=false;refresh();return;}
  if(!selected.size)return;
  if(e.key==='Escape'){e.preventDefault();e.stopImmediatePropagation();clear({focus:true});return;}
  if((e.key==='Backspace'||e.key==='Delete')&&!mod&&!e.altKey){e.preventDefault();e.stopImmediatePropagation();removeSelected();return;}
  if(mod&&key==='a'&&!e.altKey){e.preventDefault();e.stopImmediatePropagation();selectAll();return;}
  if(['ArrowUp','ArrowDown'].includes(e.key)&&!mod&&!e.altKey&&e.target.closest?.('.block-select-toggle,.block-selection-toolbar')){
   const ids=orderedIds(),at=ids.indexOf(focusId||api.getSelectedId?.()),next=ids[at+(e.key==='ArrowDown'?1:-1)];
   e.preventDefault();e.stopImmediatePropagation();if(next){select(next,{shift:e.shiftKey});focusPicker(next,{reveal:true});}return;
  }
 }
 function init(options){
  api=options||{};if(initialized)return;initialized=true;ensureToolbar();
  document.addEventListener('keydown',handleKeydown,true);
  document.addEventListener('pointerdown',e=>{
   if(!selected.size||isControl(e.target))return;
   if(e.target.closest?.('[data-direct-edit-root],.direct-math-atom,.block-insertion-gap,.block-content img'))clear();
  },true);
  // Tabbed entry into text editing must also stop an old block selection from
  // remaining armed, while leaving the caret and native text deletion alone.
  document.addEventListener('focusin',e=>{if(selected.size&&isTyping(e.target)&&!isControl(e.target))clear();});
 }
 window.ReaderBlockSelection={init,bind,refresh,clear,select,selectAll,removeSelected,addMenuItems,selectedIds:validSelected,testing:{rangeIds,isTyping}};
})();
