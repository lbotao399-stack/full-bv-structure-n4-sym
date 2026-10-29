/* Empty inter-block caret stops become real blocks only when the user writes. */
(function(){'use strict';
 let api={},restoreGap=null,materializing=false,composingGap=null;
 const uid=()=> 'local-'+(crypto.randomUUID?.()||Date.now()+'-'+Math.random().toString(36).slice(2));
 const esc=text=>String(text).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
 const gaps=()=>[...document.querySelectorAll('.block-insertion-gap')];
 const selectionEnd=el=>{el.focus({preventScroll:true});const r=document.createRange();r.selectNodeContents(el);r.collapse(false);const s=getSelection();s.removeAllRanges();s.addRange(r);el.scrollIntoView({block:'nearest'});};
 function clear(){
  const focused=document.activeElement?.closest?.('.block-insertion-gap');
  restoreGap=focused&&!materializing?{before:focused.dataset.beforeId||null,after:focused.dataset.afterId||null}:null;
  for(const el of gaps())el.remove();
 }
 function insert(beforeId,afterId,text=''){
  const blocks=api.getDocument?.()?.blocks||[];
  if(afterId&&!blocks.some(b=>b.id===afterId)||beforeId&&!blocks.some(b=>b.id===beforeId))return null;
  const next={id:uid(),type:'paragraph',source:text?'<p>'+esc(text).replace(/\r\n?|\n/g,'<br>')+'</p>':''};
  const parent=window.ReaderBlocks?.insertionParent?.(beforeId,afterId);if(parent)next.parentId=parent;
  materializing=true;try{window.ReaderDirectEdit?.flush();api.beforeChange?.();api.insert?.(beforeId,afterId,next);window.ReaderDirectEdit?.caret?.(next.id,text.length);window.ReaderDirectEdit?.refreshCommands?.(next.id);api.afterChange?.();}finally{materializing=false;restoreGap=null;}
  return next.id;
 }
 function materialize(el,text=el.textContent||''){if(!el.isConnected)return null;return insert(el.dataset.beforeId||null,el.dataset.afterId||null,text);}
 function makeGap(beforeId,afterId,position){
  const el=document.createElement('div');el.className='block-insertion-gap';el.dataset.beforeId=beforeId||'';el.dataset.afterId=afterId||'';el.dataset.directIgnore='';el.contentEditable='true';el.tabIndex=0;el.spellcheck=false;el.setAttribute('role','textbox');el.setAttribute('aria-label',afterId?'Insert block before passage '+(position+1):'Insert block at end');el.setAttribute('aria-multiline','false');el.setAttribute('aria-description','Click, type, or press Enter to insert a block.');
  // The empty seam remains clickable; reveal its control only at the left edge.
  const gutterHover=e=>{const x=e.clientX-el.getBoundingClientRect().left;el.classList.toggle('gutter-active',x>=-24&&x<=24);};
  el.addEventListener('pointermove',gutterHover);
  el.addEventListener('pointerenter',gutterHover);
  el.addEventListener('pointerleave',()=>el.classList.remove('gutter-active'));
  let composing=false;
  el.addEventListener('click',e=>{e.preventDefault();e.stopPropagation();materialize(el);});
  el.addEventListener('keydown',e=>{
   if(e.isComposing||composing||e.metaKey||e.ctrlKey||e.altKey)return;
   if(e.key==='ArrowUp'||e.key==='ArrowDown'){
    e.preventDefault();e.stopImmediatePropagation();const dir=e.key==='ArrowUp'?-1:1,id=dir<0?beforeId:afterId;
    if(id)window.ReaderDirectEdit?.navigateToBlock?.(id,dir);return;
   }
   if(e.key==='Enter'){e.preventDefault();e.stopImmediatePropagation();materialize(el);}
   else if(e.key==='Backspace'||e.key==='Delete'){if(!el.textContent){e.preventDefault();e.stopImmediatePropagation();const id=e.key==='Backspace'?beforeId:afterId;if(id)window.ReaderDirectEdit?.navigateToBlock?.(id,e.key==='Backspace'?-1:1);}}
   else if(e.key==='Escape'){e.preventDefault();el.blur();}
  });
  el.addEventListener('compositionstart',()=>{composing=true;composingGap=el;});
  el.addEventListener('compositionend',()=>{composing=false;if(composingGap===el)composingGap=null;if(el.textContent)materialize(el);});
  el.addEventListener('input',e=>{if(!composing&&!e.isComposing&&el.textContent)materialize(el);});
  el.addEventListener('paste',e=>{const text=e.clipboardData?.getData('text/plain');if(text){e.preventDefault();e.stopPropagation();materialize(el,text);}});
  return el;
 }
 function refresh(){
  if(document.body.classList.contains('floating')||!api.getDocument?.())return;
  const holder=document.querySelector('#blocks'),blocks=api.getDocument().blocks;if(!holder)return;
  // Block wrappers stay direct children so caret, drag and narration lookup remain stable.
  if(!gaps().length){
   const visible=[...holder.querySelectorAll(':scope > .source-block')].map(n=>n.dataset.blockId),indices=visible.map(id=>blocks.findIndex(b=>b.id===id)).filter(i=>i>=0);
   for(const i of indices){const before=blocks[i-1]?.id||null,after=blocks[i].id,anchor=api.getNode?.(after);if(anchor?.parentElement===holder)holder.insertBefore(makeGap(before,after,i),anchor);}
   // Page seams still address the adjacent blocks of the complete document.
   if(indices.length){const i=indices.at(-1)+1;holder.append(makeGap(blocks[i-1]?.id||null,blocks[i]?.id||null,i));}
  }
  for(const gap of gaps()){
   const parent=window.ReaderBlocks?.insertionParent?.(gap.dataset.beforeId||null,gap.dataset.afterId||null);
   gap.classList.toggle('inside-container',Boolean(parent));
  }
  if(restoreGap){const x=restoreGap;restoreGap=null;focusBetween(x.before,x.after);}
 }
 function focusBetween(beforeId,afterId){const gap=gaps().find(el=>(el.dataset.beforeId||null)===(beforeId||null)&&(el.dataset.afterId||null)===(afterId||null));if(!gap)return false;window.ReaderDirectEdit?.flush();selectionEnd(gap);return true;}
 window.ReaderInsertion={init:options=>{api=options||{};},clear,refresh,focusBetween,insert,preparePageChange:()=>!composingGap?.isConnected};
})();
