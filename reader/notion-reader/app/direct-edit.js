/* Direct manuscript editing. Rendered formula nodes remain indivisible TeX atoms. */
(function () {
  'use strict';
  let api = {}, initialized = false, active = null, mathEditor = null, slash = null, closingMath = false;
  const sessions = new WeakMap();
  let verticalGoalX = null;
  const EDITABLE_TYPES = new Set(['paragraph','heading','callout','table']);
  const ATOM = '\ufffc';
  const $ = (s, root = document) => root.querySelector(s);
  const idOf = root => root?.closest('.source-block')?.dataset.blockId;
  const rootOf = node => (node?.nodeType === 1 ? node : node?.parentElement)?.closest('[data-direct-edit-root]');
  const blockOf = root => api.getBlock?.(idOf(root));
  const rootById = id => [...document.querySelectorAll('[data-direct-edit-root]')].find(root => idOf(root) === id);
  const clone = x => JSON.parse(JSON.stringify(x));
  const uid = () => 'local-' + (crypto.randomUUID?.() || Date.now() + '-' + Math.random().toString(36).slice(2));
  const escapeText = text => String(text).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/[\\$*_\[\]`|]/g,c=>'&#'+c.charCodeAt(0)+';').replace(/\u200b/g,'');
  const escapeAttr = text => escapeText(text).replace(/"/g,'&quot;');
  const same = (a,b) => JSON.stringify(a) === JSON.stringify(b);

  function texOf(node) { return node.dataset.directTex ?? node.querySelector('annotation[encoding="application/x-tex"]')?.textContent ?? node.querySelector('.math-error')?.textContent ?? '';  }
  function safeStyle(node) {
    if (node.matches('[data-reader-format],.reader-format-mark')) return '';
    const properties=['color','background-color','font-weight','font-style','text-decoration-line','text-align'];
    return properties.map(key => {const value=node.style?.getPropertyValue(key)||'';return value && !/url|expression|[<>]/i.test(value) ? key+':'+value : '';}).filter(Boolean).join(';');
  }
  function serializeNode(node) {
    if (node.nodeType === Node.TEXT_NODE) return escapeText(node.data);
    if (node.nodeType !== Node.ELEMENT_NODE && node.nodeType !== Node.DOCUMENT_FRAGMENT_NODE) return '';
    if (node.nodeType === Node.DOCUMENT_FRAGMENT_NODE) return [...node.childNodes].map(serializeNode).join('');
    if (node.matches('[data-direct-ignore],.katex-mathml,script,style,.empty-block')) return '';
    if (node.matches('.math-inline,.math-display')) {
      const tex=texOf(node);return node.matches('.math-display') ? '\n\n$$\n'+tex+'\n$$\n\n' : '\\('+tex+'\\)';
    }
    let tag=node.tagName.toLowerCase(),body=[...node.childNodes].map(serializeNode).join('');
    if (node.hasAttribute('data-reader-format')) return body;
    if (tag==='br') return '<br>';
    if (tag==='img') {
      const src=node.getAttribute('src')||'';return /^(?:https?:|figures\/|\.\.?\/)/i.test(src) ? '<img src="'+escapeAttr(src)+'" alt="'+escapeAttr(node.getAttribute('alt')||'')+'">' : '';
    }
    if (!['p','div','span','strong','b','em','i','u','s','del','code','pre','a','mark','sub','sup','blockquote','ul','ol','li','table','thead','tbody','tfoot','tr','td','th','h1','h2','h3','h4','h5','h6'].includes(tag)) return body;
    if (tag==='div') tag='p';
    let attrs='';
    if (tag==='a') {
      const href=node.getAttribute('href')||'';
      if (/^(?:https?:|mailto:|#|\/|\.\.?\/)/i.test(href)) attrs+=' href="'+escapeAttr(href)+'"';
      if (node.getAttribute('title')) attrs+=' title="'+escapeAttr(node.getAttribute('title'))+'"';
    }
    if ((tag==='td'||tag==='th')) for(const key of ['colspan','rowspan']) {const n=Number(node.getAttribute(key));if(Number.isInteger(n)&&n>1&&n<100)attrs+=' '+key+'="'+n+'"';}
    if(tag==='ol'&&node.hasAttribute('start')){const start=Number(node.getAttribute('start'));if(Number.isInteger(start))attrs+=' start="'+start+'"';}
    let style=safeStyle(node);if(node.classList.contains('source-highlight'))style+=(style?';':'')+'background-color:#fff0a7';if(style)attrs+=' style="'+escapeAttr(style)+'"';
    return '<'+tag+attrs+'>'+body+'</'+tag+'>';
  }
  function serializeEditable(root, block = {}) {
    const copy=root.cloneNode(true);
    for(const n of copy.querySelectorAll('[data-direct-ignore]'))n.remove();
    if(block.type==='heading') {
      const head=copy.querySelector('h1,h2,h3,h4,h5,h6');
      if(head) return '<p>'+[...head.childNodes].map(serializeNode).join('')+'</p>';
    }
    let result=[...copy.childNodes].map(serializeNode).join('');
    if(!result.replace(/<[^>]*>|&nbsp;|\s/g,''))return '';
    if(!/^\s*<(?:p|h[1-6]|table|ul|ol|blockquote|pre)\b/.test(result))result='<p>'+result+'</p>';
    return result;
  }
  function linearProjection(root) {
    const units=[];let text='';
    function walk(n) {
      if(n.nodeType===Node.TEXT_NODE){if(n.data){units.push({node:n,start:text.length,end:text.length+n.data.length,type:'text'});text+=n.data;}return;}
      if(n.nodeType!==Node.ELEMENT_NODE&&n.nodeType!==Node.DOCUMENT_FRAGMENT_NODE)return;
      if(n.nodeType===Node.ELEMENT_NODE&&n.matches('[data-direct-ignore],.katex-mathml,script,style'))return;
      if(n.nodeType===Node.ELEMENT_NODE&&n.matches('.math-inline,.math-display')){units.push({node:n,start:text.length,end:text.length+1,type:'math',tex:texOf(n)});text+=ATOM;return;}
      for(const child of n.childNodes)walk(child);
    }walk(root);return {text,units};
  }
  function extractFormatting(root) {
    const tools=window.ReaderFormatting?.testing;
    if(!tools?.projection||!tools?.anchor)return null;
    const p=tools.projection(root),runs=[];
    for(const unit of p.units){
      const chain=[];let n=unit.type==='math'?unit.node:unit.node.parentElement;
      while(n&&n!==root){if(n.matches?.('.reader-format-mark'))chain.unshift(n);n=n.parentElement;}
      const style={};
      for(const el of chain){
        if(el.dataset.readerColor)style.color=el.dataset.readerColor;
        if(el.dataset.readerHighlight)style.highlight=el.dataset.readerHighlight;
        if(Number(el.style.fontWeight)>=600||el.style.fontWeight==='bold')style.bold=true;
        if(el.style.fontStyle==='italic')style.italic=true;
        if(el.style.textDecorationLine.includes('underline'))style.underline=true;
      }
      const cleaned=tools.cleanStyle(style);if(!Object.keys(cleaned).length)continue;
      const prev=runs[runs.length-1];if(prev&&prev.end===unit.start&&same(prev.style,cleaned))prev.end=unit.end;else runs.push({start:unit.start,end:unit.end,style:cleaned});
    }
    return runs.length?{version:1,text:p.text,ranges:runs.map(run=>tools.anchor(p,run.start,run.end,run.style))}:null;
  }
  function pointOffset(root,node,offset) {
    try {const range=document.createRange();range.selectNodeContents(root);range.setEnd(node,offset);return linearProjection(range.cloneContents()).text.length;}catch(_){return 0;}
  }
  function domPoint(root,offset) {
    const p=linearProjection(root);offset=Math.max(0,Math.min(offset,p.text.length));
    for(const u of p.units){
      if(offset<=u.end){if(u.type==='text')return [u.node,Math.max(0,offset-u.start)];const parent=u.node.parentNode,index=[...parent.childNodes].indexOf(u.node);return [parent,index+(offset>u.start?1:0)];}
    }
    return [root,root.childNodes.length];
  }
  function snapshotSelection() {
    const s=window.getSelection();if(!s?.rangeCount)return null;
    const a=rootOf(s.anchorNode),f=rootOf(s.focusNode);if(!a||!f)return null;
    return {anchor:{id:idOf(a),offset:pointOffset(a,s.anchorNode,s.anchorOffset)},focus:{id:idOf(f),offset:pointOffset(f,s.focusNode,s.focusOffset)}};
  }
  function restoreSelection(snapshot,{focus=true}={}) {
    if(!snapshot||window.ReaderPageNavigation?.ensureBlock?.(snapshot.anchor.id)===false)return false;const a=rootById(snapshot.anchor.id),f=rootById(snapshot.focus.id);if(!a||!f)return false;
    const ap=domPoint(a,snapshot.anchor.offset),fp=domPoint(f,snapshot.focus.offset);if(focus)a.focus({preventScroll:true});
    const s=window.getSelection();s.removeAllRanges();if(s.setBaseAndExtent)s.setBaseAndExtent(...ap,...fp);else{const r=document.createRange();r.setStart(...ap);r.setEnd(...fp);s.addRange(r);}return true;
  }
  function preserveSelection(fn){const snapshot=snapshotSelection();const result=fn();restoreSelection(snapshot);return result;}
  function caret(id,offset=0){restoreSelection({anchor:{id,offset},focus:{id,offset}});}

  function selectionRect(range) {
    const rects=[...range.getClientRects()].filter(r=>r.height>0);
    if(rects.length)return rects[0];
    if(!range.collapsed)return null;
    const container=range.startContainer,at=range.startOffset;
    if(container.nodeType!==Node.ELEMENT_NODE)return null;
    const after=at>=container.childNodes.length,node=container.childNodes[after?at-1:at];
    if(!node)return null;
    const probe=document.createRange();
    if(node.nodeType===Node.TEXT_NODE){if(!node.length)return null;probe.setStart(node,after?node.length-1:0);probe.setEnd(node,after?node.length:1);}
    else probe.selectNode(node);
    const rect=[...probe.getClientRects()].find(r=>r.height>0);
    return rect?{left:after?rect.right:rect.left,right:after?rect.right:rect.left,top:rect.top,bottom:rect.bottom,height:rect.height,width:0}:null;
  }
  function selectionInRoot(root) {
    const selection=window.getSelection();
    return selection?.rangeCount && selection.isCollapsed && root.contains(selection.anchorNode) ? selection : null;
  }
  function placePoint(root,node,offset) {
    // Hit testing on a rendered formula may return one of KaTeX's internal nodes.
    const element=node.nodeType===1?node:node.parentElement;
    const atom=element?.closest('.direct-math-atom');
    if(atom&&root.contains(atom)){
      const r=atom.getBoundingClientRect(),parent=atom.parentNode,index=[...parent.childNodes].indexOf(atom);
      offset=index+(verticalGoalX!=null&&verticalGoalX>r.left+r.width/2?1:0);node=parent;
    }
    const selection=window.getSelection(),range=document.createRange();
    range.setStart(node,offset);range.collapse(true);selection.removeAllRanges();selection.addRange(range);
  }
  function revealVerticalPoint(y,height=20) {
    const nav=document.querySelector('.page-navigation-top');
    const top=Math.max(document.querySelector('.topbar')?.getBoundingClientRect().bottom||0,nav&&!nav.hidden?nav.getBoundingClientRect().bottom:0)+12;
    const bottom=(document.querySelector('.player')?.getBoundingClientRect().top||innerHeight)-12;
    if(y<top)window.scrollBy(0,y-top);
    else if(y+height>bottom)window.scrollBy(0,y+height-bottom);
  }
  function edgeGeometry(root,direction) {
    const range=document.createRange();range.selectNodeContents(root);
    const rects=[...range.getClientRects()].filter(r=>r.height>0);
    const rect=direction>0?rects[0]:rects.at(-1),box=root.getBoundingClientRect();
    const lineHeight=parseFloat(getComputedStyle(root).lineHeight)||parseFloat(getComputedStyle(root).fontSize)*1.5||20;
    if(!rect)return {x:box.left,y:direction>0?box.top:Math.max(box.top,box.bottom-lineHeight),height:lineHeight};
    return {x:rect.left,y:direction>0?rect.top:rect.bottom-Math.min(rect.height,lineHeight),height:Math.min(rect.height,lineHeight)};
  }
  function navigateToBlock(id,direction=1,preferredX=null) {
    if(window.ReaderBlocks?.emptyHeader?.(id)){
      const blocks=api.getDocument?.().blocks||[],at=blocks.findIndex(b=>b.id===id);
      let next=at+direction;
      while(blocks[next]&&window.ReaderBlocks.emptyHeader(blocks[next].id))next+=direction;
      const target=blocks[next]||blocks[at+1];
      if(target)return navigateToBlock(target.id,direction,preferredX);
    }
    if(window.ReaderPageNavigation?.ensureBlock?.(id)===false)return false;
    const b=api.getBlock?.(id),wrap=[...document.querySelectorAll('.source-block')].find(n=>n.dataset.blockId===id);
    if(!b||!wrap)return false;
    flush();hideSlash();if(Number.isFinite(preferredX))verticalGoalX=preferredX;
    const root=rootById(id);
    if(!root){
      const target=wrap.querySelector('.direct-math-atom,.block-drag-image')||wrap;
      target.focus({preventScroll:true});window.getSelection()?.removeAllRanges();api.select?.(id);
      const rect=target.getBoundingClientRect();revealVerticalPoint(rect.top,Math.min(rect.height,60));return true;
    }
    root.focus({preventScroll:true});
    // Start from a valid editable point; native hit testing chooses the nearest
    // glyph on the first/last visual line, including proportional fonts and CJK.
    caret(id,direction>0?0:linearProjection(root).text.length);
    let edge=edgeGeometry(root,direction);revealVerticalPoint(edge.y,edge.height);edge=edgeGeometry(root,direction);
    const box=root.getBoundingClientRect(),x=Math.max(box.left+1,Math.min(box.right-1,verticalGoalX??edge.x));
    const y=edge.y+Math.max(1,edge.height/2);
    let node=null,offset=0;
    if(document.caretPositionFromPoint){const point=document.caretPositionFromPoint(x,y);node=point?.offsetNode;offset=point?.offset??0;}
    else if(document.caretRangeFromPoint){const point=document.caretRangeFromPoint(x,y);node=point?.startContainer;offset=point?.startOffset??0;}
    if(node&&root.contains(node))placePoint(root,node,offset);
    api.select?.(id);return true;
  }
  function adjacentVertical(id,direction) {
    const blocks=api.getDocument?.().blocks||[],index=blocks.findIndex(b=>b.id===id);
    const next=blocks[index+direction];
    return next?navigateToBlock(next.id,direction,verticalGoalX):false;
  }
  function verticalKey(event,root) {
    const direction=event.key==='ArrowDown'?1:-1;
    if(event.metaKey||event.ctrlKey||event.altKey||event.shiftKey||event.isComposing)return false;
    const wrap=event.target.closest?.('.source-block'),id=wrap?.dataset.blockId;
    if(!id)return false;
    const b=api.getBlock?.(id);
    if(!root){
      if(!['equation','image'].includes(b?.type))return false;
      const rect=event.target.getBoundingClientRect();if(verticalGoalX==null)verticalGoalX=rect.left;
      const blocks=api.getDocument?.().blocks||[],index=blocks.findIndex(block=>block.id===id);
      const previous=blocks[index-1]?.id||null,next=blocks[index+1]?.id||null;
      const gap=window.ReaderInsertion?.focusBetween?.(direction>0?id:previous,direction>0?next:id);
      if((gap!==false&&gap!==undefined)||adjacentVertical(id,direction)){
        event.preventDefault();event.stopImmediatePropagation();return true;
      }return false;
    }
    if(sessions.get(root)?.composing)return false;
    const focusedAtom=event.target.matches?.('.direct-math-atom')?event.target:null;
    if(focusedAtom){
      const r=document.createRange();direction>0?r.setStartAfter(focusedAtom):r.setStartBefore(focusedAtom);r.collapse(true);
      const s=window.getSelection();s.removeAllRanges();s.addRange(r);root.focus({preventScroll:true});
    }
    const selection=selectionInRoot(root);if(!selection)return false;
    const original=selection.getRangeAt(0).cloneRange(),originalRect=selectionRect(original);
    if(verticalGoalX==null)verticalGoalX=originalRect?.left??root.getBoundingClientRect().left;
    let nativeMoved=false;
    if(typeof selection.modify==='function'){
      selection.modify('move',direction>0?'forward':'backward','line');
      const now=selection.rangeCount?selection.getRangeAt(0):null,newRect=now?selectionRect(now):null;
      // Native movement inside a different visual line is already exactly what
      // the editor needs. Do not replace it with logical text-offset arithmetic.
      nativeMoved=Boolean(now&&root.contains(now.startContainer)&&newRect&&originalRect&&
        (direction>0?newRect.top>=originalRect.bottom-2:newRect.bottom<=originalRect.top+2));
      if(nativeMoved){event.preventDefault();event.stopImmediatePropagation();revealVerticalPoint(newRect.top,newRect.height);return true;}
      selection.removeAllRanges();selection.addRange(original);
    }else{
      const edge=edgeGeometry(root,direction<0?1:-1);
      if(originalRect&&(direction>0?originalRect.bottom<edge.y+edge.height-3:originalRect.top>edge.y+3))return false;
    }
    if(adjacentVertical(id,direction)){
      event.preventDefault();event.stopImmediatePropagation();return true;
    }
    // At the document's ends retain the browser's usual first/last-line action.
    return false;
  }

  function before(root,group='typing'){api.onBeforeChange?.(idOf(root),{group});}
  function commit(root) {
    const session=sessions.get(root),b=blockOf(root);if(!session||!b||session.composing||!session.dirty)return false;
    session.dirty=false;const source=serializeEditable(root,b),formatting=extractFormatting(root);
    const patch=session.field==='caption'?{caption:source,formatting}:{source,formatting};
    api.commitBlock?.(b.id,patch,{render:false});session.html=root.innerHTML;
    if(window.ReaderBlocks?.isContainer?.(b))window.ReaderBlocks.refreshStructure();
    return true;
  }
  function flush(){if(active)commit(active);}
  function preparePageChange(){if(sessions.get(active)?.composing)return false;if(mathEditor&&!commitMath())return false;flush();hideSlash();return true;}
  function changed(root){const session=sessions.get(root);if(!session)return;session.dirty=true;if(!session.composing){commit(root);updateSlash(root);}}
  function isEditing(){return Boolean(mathEditor||document.activeElement?.closest?.('[data-direct-edit-root]'));}
  function singleSelection(root){const s=window.getSelection();if(!s?.rangeCount)return null;const r=s.getRangeAt(0);return root.contains(r.startContainer)&&root.contains(r.endContainer)?r:null;}
  function mutableSelection(root) {
    const r=singleSelection(root);if(r)return r;api.toast?.('Edit one passage at a time. A selection across passages can still be highlighted.');return null;
  }
  function insertNodeAtRange(root,range,node) {
    range.deleteContents();range.insertNode(node);const tail=document.createTextNode('');node.after(tail);
    const next=document.createRange();next.setStart(tail,0);next.collapse(true);const selection=window.getSelection();selection.removeAllRanges();selection.addRange(next);root.focus({preventScroll:true});
  }
  function insertBreak(root) {
    const r=mutableSelection(root);if(!r)return;before(root,'structural');
    const br=document.createElement('br');insertNodeAtRange(root,r,br);if(!br.nextSibling?.nextSibling)br.after(document.createElement('br'));changed(root);
  }
  function splitBlock(root) {
    const b=blockOf(root),r=mutableSelection(root);if(!b||!r)return;
    if(b.type==='table'||root.closest('figcaption')||r.startContainer.parentElement?.closest('td,th,pre')){insertBreak(root);return;}
    before(root,'structural');r.deleteContents();
    const left=document.createRange();left.selectNodeContents(root);left.setEnd(r.startContainer,r.startOffset);
    const right=document.createRange();right.selectNodeContents(root);right.setStart(r.startContainer,r.startOffset);
    const lf=document.createElement('div'),rf=document.createElement('div');lf.append(left.cloneContents());rf.append(right.cloneContents());
    const leftPatch={source:serializeEditable(lf,b),formatting:extractFormatting(lf)};
    const next={id:uid(),type:'paragraph',source:serializeEditable(rf,b),formatting:extractFormatting(rf)};
    if(!next.formatting)delete next.formatting;
    hideSlash();sessions.get(root).dirty=false;api.split?.(b.id,leftPatch,next);caret(next.id,0);
  }
  function deleteEmptyContainerLine(root,direction) {
    const b=blockOf(root),selection=singleSelection(root);
    if(!b||!selection?.collapsed||serializeEditable(root,b))return false;
    const children=window.ReaderBlocks?.children?.(b.id)||[];
    if(window.ReaderBlocks?.isContainer?.(b)){
      if(!children.length)return false;
      flush();navigateToBlock(children[0].id,1);return true;
    }
    if(!b.parentId||b.type!=='paragraph')return false;
    const siblings=window.ReaderBlocks.children(b.parentId),at=siblings.findIndex(x=>x.id===b.id);
    const target=siblings[at+direction]||siblings[at-direction]||api.getBlock(b.parentId);
    if(!target||!window.ReaderBlocks.deleteBlock(b.id))return false;
    navigateToBlock(target.id,direction);return true;
  }
  function mergePrevious(root) {
    const b=blockOf(root),s=window.getSelection(),r=singleSelection(root);if(!b||!r||!s.isCollapsed||pointOffset(root,r.startContainer,r.startOffset)!==0)return false;
    const blocks=api.getDocument?.().blocks||[],i=blocks.findIndex(x=>x.id===b.id),previous=blocks[i-1];
    // Text deletion cannot consume a container or cross its outer boundary.
    if(window.ReaderBlocks?.isContainer?.(b)||previous&&!window.ReaderBlocks?.canMerge?.(previous.id,b.id))return true;
    if(!previous||!['paragraph','heading','callout'].includes(previous.type)||!['paragraph','heading','callout'].includes(b.type))return false;
    before(root,'structural');
    if(window.ReaderPageNavigation?.ensureBlock?.(previous.id)===false)return false;
    const earlier=rootById(previous.id);if(!earlier)return false;
    const combined=document.createElement('div');
    function contents(source){if(!serializeEditable(source))return [];const copy=source.cloneNode(true);const wrapper=copy.children.length===1&&copy.firstElementChild.matches('p,h1,h2,h3,h4,h5,h6')?copy.firstElementChild:copy;return [...wrapper.childNodes].map(x=>x.cloneNode(true));}
    const p=document.createElement('p');p.append(...contents(earlier));const offset=linearProjection(p).text.length;p.append(...contents(root));combined.append(p);
    api.merge?.(previous.id,b.id,{source:serializeEditable(combined,previous),formatting:extractFormatting(combined)});caret(previous.id,offset);return true;
  }
  function markAtoms(root) {
    for(const math of root.querySelectorAll('.math-inline,.math-display')){
      math.contentEditable='false';math.tabIndex=0;math.setAttribute('role','button');math.setAttribute('aria-label','Edit formula: '+texOf(math));math.classList.add('direct-math-atom');
    }
  }
  function mathMarkup(tex,display){return api.renderMath?.(tex,display)||window.katex.renderToString(tex,{displayMode:display,throwOnError:false});}
  function validMath(tex,display=false){try{window.katex.renderToString(tex,{throwOnError:true,displayMode:display,strict:'ignore',trust:false});return null;}catch(e){return e.message;}}
  function mathAtom(tex,display=false){const span=document.createElement(display?'div':'span');span.className=display?'math-display':'math-inline';span.dataset.directTex=tex;span.innerHTML=mathMarkup(tex,display);markAtoms({querySelectorAll:()=>[span]});return span;}
  function closeMath({cancel=false}={}) {
    if(!mathEditor||closingMath)return;closingMath=true;const current=mathEditor;mathEditor=null;current.panel.remove();
    if(cancel&&current.cancelPatch){api.commitBlock?.(current.id,current.cancelPatch,{render:true});const restored=rootById(current.id);if(restored)caret(current.id,linearProjection(restored).text.length);}
    else if(cancel&&current.snapshot)restoreSelection(current.snapshot);
    else if(cancel&&current.displayBlock)current.node?.focus({preventScroll:true});
    closingMath=false;
  }
  function commitMath() {
    if(!mathEditor)return;const current=mathEditor;
    if(current.input.value===current.tex&&!current.cancelPatch&&(current.node||current.displayBlock)){const snapshot=current.snapshot;closeMath();if(snapshot)restoreSelection(snapshot);else current.node?.focus({preventScroll:true});return true;}
    const tex=current.input.value.trim(),error=validMath(tex,current.displayBlock);
    if(error){current.error.textContent=error;current.input.focus();return false;}
    if(!tex){current.error.textContent='Enter a LaTeX expression.';current.input.focus();return false;}
    if(current.displayBlock){
      api.onBeforeChange?.(current.id,{group:'math'});api.commitBlock?.(current.id,{source:tex,type:'equation'},{render:true});closeMath();
      const wrap=[...document.querySelectorAll('.source-block')].find(n=>n.dataset.blockId===current.id);wrap?.querySelector('.direct-math-atom')?.focus({preventScroll:true});return true;
    }
    const root=rootById(current.id);if(!root||!current.range?.startContainer?.isConnected){closeMath();api.toast?.('The passage changed. Select the formula again.');return false;}
    before(root,'math');const atom=mathAtom(tex,false);
    if(current.node?.isConnected){
      // Annotation styles are presentation data, retained on the replacement atom.
      for(const key of ['readerColor','readerHighlight'])if(current.node.dataset[key])atom.dataset[key]=current.node.dataset[key];
      if(current.node.classList.contains('reader-format-mark')){atom.classList.add('reader-format-mark');atom.style.cssText=current.node.style.cssText;}
      const range=document.createRange();range.selectNode(current.node);insertNodeAtRange(root,range,atom);
    }else insertNodeAtRange(root,current.range,atom);
    closeMath();changed(root);return true;
  }
  function openMath({root,node=null,range=null,tex='',displayBlock=false,id=null,cancelPatch=null}) {
    if(mathEditor&&!commitMath())return;flush();hideSlash();closeMath();id=id||idOf(root);if(!id)return;
    const snapshot=snapshotSelection();const panel=document.createElement('div');panel.className='direct-math-editor';panel.setAttribute('role','dialog');panel.setAttribute('aria-label',displayBlock?'Edit display equation':'Edit inline equation');panel.dataset.directIgnore='';
    const label=document.createElement('label');label.textContent=displayBlock?'Display equation':'Inline equation';
    const input=document.createElement('textarea');input.setAttribute('aria-label',displayBlock?'Display equation LaTeX':'Inline equation LaTeX');input.spellcheck=false;input.value=tex;input.rows=displayBlock?3:2;
    const hint=document.createElement('p');hint.className='direct-edit-hint';hint.textContent='Enter to render · Shift+Enter for a new line · Esc to cancel';
    const preview=document.createElement('div');preview.className='direct-math-preview';const error=document.createElement('div');error.className='direct-math-error';error.setAttribute('role','status');
    const controls=document.createElement('div');controls.className='direct-math-actions';const done=document.createElement('button');done.type='button';done.textContent='Done';done.onclick=commitMath;const cancel=document.createElement('button');cancel.type='button';cancel.textContent='Cancel';cancel.onclick=()=>closeMath({cancel:true});controls.append(cancel,done);
    panel.append(label,input,hint,preview,error,controls);document.body.append(panel);
    const targetRect=node?.getBoundingClientRect()||range?.getBoundingClientRect()||root?.getBoundingClientRect();const width=Math.min(520,window.innerWidth-24);panel.style.width=width+'px';const rect=targetRect||{left:window.innerWidth/2,top:80,bottom:100};panel.style.left=Math.max(12,Math.min(window.innerWidth-width-12,rect.left))+'px';panel.style.top=Math.max(70,Math.min(window.innerHeight-320,rect.bottom+8))+'px';
    mathEditor={id,root,node,range:range?.cloneRange(),tex,displayBlock,panel,input,error,snapshot,cancelPatch};
    input.oninput=()=>{preview.innerHTML=mathMarkup(input.value,displayBlock);error.textContent='';};input.oninput();
    panel.addEventListener('keydown',event=>{event.stopPropagation();if(event.key==='Escape'){event.preventDefault();closeMath({cancel:true});}else if(event.key==='Enter'&&!event.shiftKey&&!event.isComposing){event.preventDefault();commitMath();}});
    input.focus();input.setSelectionRange(input.value.length,input.value.length);
  }
  function inlineMath() {
    const selection=window.getSelection(),root=rootOf(selection?.anchorNode);if(!root)return false;const range=mutableSelection(root);if(!range)return true;
    const fragment=range.cloneContents(),selectedAtoms=fragment.querySelectorAll('.math-inline,.math-display');
    if(selectedAtoms.length===1&&linearProjection(fragment).text.replace(ATOM,'').trim()===''){const original=[...root.querySelectorAll('.math-inline,.math-display')].find(n=>range.intersectsNode(n));if(original){const atomRange=document.createRange();atomRange.selectNode(original);openMath({root,node:original,range:atomRange,tex:texOf(original)});return true;}}
    const selected=range.toString().trim();if(!range.collapsed&&selected){const tex=selected.replace(/^\$([\s\S]*)\$$/,'$1').replace(/^\\\(([\s\S]*)\\\)$/,'$1');
      if(!validMath(tex)){before(root,'math');insertNodeAtRange(root,range,mathAtom(tex));changed(root);return true;}
      openMath({root,range,tex});return true;
    }
    openMath({root,range,tex:''});return true;
  }
  function hideSlash(){if(slash){slash.panel.remove();slash=null;}}
  const slashCommands=[
    {key:'math',label:'Math equation',hint:'/math',aliases:['math','equation']},
    {key:'quote',label:'Quote',hint:'/quote',aliases:['quote','q','ut']},
    {key:'callout',label:'Callout',hint:'/callout',aliases:['callout','out']},
    {key:'heading-1',label:'Heading 1',hint:'/h1',aliases:['h1','heading1']},
    {key:'heading-2',label:'Heading 2',hint:'/h2',aliases:['h2','heading2']},
    {key:'heading-3',label:'Heading 3',hint:'/h3',aliases:['h3','heading3']},
    {key:'paragraph',label:'Text',hint:'/text',aliases:['text','paragraph']}
  ];
  function makeDisplay(root) {
    const b=blockOf(root);if(!b)return;before(root,'structural');hideSlash();sessions.get(root).dirty=false;
    const cancelPatch={type:b.type,level:b.level??null,blockStyle:b.blockStyle??null,source:b.source,formatting:b.formatting?clone(b.formatting):null};
    api.commitBlock?.(b.id,{type:'equation',level:null,blockStyle:null,source:'',formatting:null},{render:true});openMath({root:null,id:b.id,tex:'',displayBlock:true,cancelPatch});
  }
  function applySlash(root,command){
    if(!command)return;if(command.key==='math'){makeDisplay(root);return;}
    const b=blockOf(root);if(!b)return;before(root,'structural');hideSlash();sessions.get(root).dirty=false;
    const patch={type:'paragraph',level:null,blockStyle:null,source:'',formatting:null};
    if(command.key==='callout')patch.type='callout';
    else if(command.key==='quote')patch.blockStyle='quote';
    else if(command.key.startsWith('heading-')){patch.type='heading';patch.level=Number(command.key.slice(-1));}
    api.commitBlock?.(b.id,patch,{render:true});caret(b.id,0);
  }
  function selectSlash(index){
    if(!slash)return;slash.index=(index+slash.commands.length)%slash.commands.length;
    [...slash.panel.querySelectorAll('[role=option]')].forEach((item,i)=>{item.setAttribute('aria-selected',String(i===slash.index));if(i===slash.index)item.scrollIntoView({block:'nearest'});});
  }
  function updateSlash(root) {
    const b=blockOf(root);if(!b||!['paragraph','heading','callout'].includes(b.type)){hideSlash();return;}
    const text=root.textContent.replace(/\u200b/g,'').trim(),match=text.match(/^\/([a-z0-9-]*)$/i);if(!match){hideSlash();return;}
    const query=match[1].toLowerCase(),commands=slashCommands.filter(c=>c.aliases.some(a=>a.startsWith(query)));
    if(!commands.length){hideSlash();return;}if(slash?.root===root&&slash.query===query)return;
    const previous=slash?.commands[slash.index]?.key;hideSlash();const panel=document.createElement('div');panel.className='direct-slash-menu';panel.dataset.directIgnore='';panel.setAttribute('role','listbox');panel.setAttribute('aria-label','Block commands');
    commands.forEach((command,index)=>{const item=document.createElement('button');item.type='button';item.setAttribute('role','option');item.setAttribute('aria-selected','false');const label=document.createElement('span');label.textContent=command.label;const hint=document.createElement('small');hint.textContent=command.hint;item.append(label,hint);item.onpointerdown=e=>e.preventDefault();item.onpointermove=()=>selectSlash(index);item.onclick=()=>applySlash(root,command);panel.append(item);});document.body.append(panel);
    const selection=window.getSelection(),r=(selection?.rangeCount?selection.getRangeAt(0).getBoundingClientRect():null)||root.getBoundingClientRect();panel.style.left=Math.max(12,Math.min(window.innerWidth-250,r.left))+'px';panel.style.top=Math.max(70,Math.min(window.innerHeight-panel.offsetHeight-14,r.bottom+9))+'px';slash={panel,root,query,commands,index:Math.max(0,commands.findIndex(c=>c.key===previous))};selectSlash(slash.index);
  }
  function performUndo(redo){
    flush();const snapshot=snapshotSelection(),blocks=api.getDocument?.().blocks||[],index=blocks.findIndex(b=>b.id===snapshot?.anchor.id),previous=blocks[index-1]?.id,next=blocks[index+1]?.id;
    api.undo?.(redo);if(restoreSelection(snapshot))return;
    const id=api.getBlock?.(previous)?previous:api.getBlock?.(next)?next:null;
    if(id&&window.ReaderPageNavigation?.ensureBlock?.(id)!==false){const root=rootById(id);if(root)caret(id,linearProjection(root).text.length);}
  }
  function handleKeydown(event) {
    if(event.defaultPrevented||event.isComposing)return false;
    if(!['ArrowUp','ArrowDown'].includes(event.key))verticalGoalX=null;
    const key=event.key.toLowerCase(),mod=event.metaKey||event.ctrlKey;
    if(event.target.closest?.('input,textarea,select'))return false;
    if(event.target.matches?.('.direct-math-atom')&&['Enter',' '].includes(event.key)){event.preventDefault();event.stopImmediatePropagation();event.target.click();return true;}
    if(mod&&event.shiftKey&&key==='e'){if(inlineMath()){event.preventDefault();event.stopImmediatePropagation();return true;}return false;}
    const root=rootOf(event.target);
    if(mod&&key==='z'&&(root||event.target.closest?.('.source-block'))){event.preventDefault();event.stopImmediatePropagation();performUndo(event.shiftKey);return true;}
    if(!root){if(['ArrowUp','ArrowDown'].includes(event.key))return verticalKey(event,null);return false;}
    if(slash?.root===root&&!mod&&!event.altKey){
      if(['ArrowDown','ArrowUp'].includes(event.key)){event.preventDefault();event.stopPropagation();selectSlash(slash.index+(event.key==='ArrowDown'?1:-1));return true;}
      if(event.key==='Escape'){event.preventDefault();event.stopPropagation();hideSlash();return true;}
      if((event.key==='Enter'||event.key==='Tab')&&!event.shiftKey){event.preventDefault();event.stopPropagation();applySlash(root,slash.commands[slash.index]);return true;}
      if(event.key===' '){const exact=slash.commands.find(c=>c.aliases.includes(slash.query));if(exact){event.preventDefault();event.stopPropagation();applySlash(root,exact);return true;}}
    }
    if(['ArrowUp','ArrowDown'].includes(event.key)&&verticalKey(event,root))return true;
    if(event.key==='Escape'){hideSlash();flush();root.blur();return false;}
    if(event.key==='Enter'&&!mod&&!event.altKey){event.preventDefault();event.stopPropagation();if(event.shiftKey)insertBreak(root);else splitBlock(root);return true;}
    if(['Backspace','Delete'].includes(event.key)&&!mod&&!event.altKey){
      if(deleteEmptyContainerLine(root,event.key==='Backspace'?-1:1)||event.key==='Backspace'&&mergePrevious(root)){event.preventDefault();event.stopPropagation();return true;}
    }
    return false;
  }
  function bind(block,content) {
    if(document.body.classList.contains('floating'))return;
    markAtoms(content);
    content.addEventListener('click',event=>{
      const atom=event.target.closest('.math-inline,.math-display');
      if(atom&&content.contains(atom)){event.preventDefault();event.stopPropagation();api.select?.(block.id);const range=document.createRange();range.selectNode(atom);openMath({root:rootOf(atom),node:atom,range,tex:texOf(atom),id:block.id,displayBlock:block.type==='equation'});return;}
      if(event.target.closest('a')&&!event.metaKey&&!event.ctrlKey)event.preventDefault();
    });
    content.addEventListener('keydown',event=>{if(event.target.matches('.direct-math-atom')&&['Enter',' '].includes(event.key)){event.preventDefault();event.stopPropagation();event.target.click();}});
    let root=content,field='source';
    if(block.type==='image'){root=content.querySelector('figcaption');field='caption';if(!root)return;}
    else if(!EDITABLE_TYPES.has(block.type))return;
    const placeholder=root.querySelector('.empty-block');if(placeholder)placeholder.replaceWith(document.createElement('br'));
    root.contentEditable='true';root.dataset.directEditRoot='';root.setAttribute('role','textbox');root.setAttribute('aria-multiline','true');root.setAttribute('aria-label',field==='caption'?'Edit image caption':'Edit manuscript passage');root.spellcheck=true;
    sessions.set(root,{id:block.id,field,html:root.innerHTML,dirty:false,composing:false});
    root.addEventListener('focusin',()=>{if(active&&active!==root)commit(active);active=root;api.select?.(block.id);});
    root.addEventListener('beforeinput',event=>{
      if(event.inputType==='historyUndo'||event.inputType==='historyRedo'){event.preventDefault();performUndo(event.inputType==='historyRedo');return;}
      if(!event.isComposing&&['deleteContentBackward','deleteContentForward'].includes(event.inputType)){
        if(deleteEmptyContainerLine(root,event.inputType==='deleteContentBackward'?-1:1)||event.inputType==='deleteContentBackward'&&mergePrevious(root)){event.preventDefault();return;}
      }
      const selection=window.getSelection();if(selection?.rangeCount&&rootOf(selection.anchorNode)!==rootOf(selection.focusNode)){event.preventDefault();api.toast?.('Edit one passage at a time. Selections across passages can be highlighted.');return;}
      if(!sessions.get(root).composing)before(root,'typing');
    });
    root.addEventListener('input',()=>changed(root));
    root.addEventListener('compositionstart',()=>{before(root,'typing');sessions.get(root).composing=true;hideSlash();});
    root.addEventListener('compositionend',()=>{sessions.get(root).composing=false;changed(root);});
    root.addEventListener('blur',()=>{commit(root);setTimeout(()=>{if(slash?.root===root&&!slash.panel.contains(document.activeElement))hideSlash();},0);});
    root.addEventListener('paste',event=>{
      const text=event.clipboardData?.getData('text/plain');if(text===undefined)return;event.preventDefault();const range=mutableSelection(root);if(!range)return;before(root,'typing');range.deleteContents();const fragment=document.createDocumentFragment();const lines=text.replace(/\r\n?/g,'\n').split('\n');lines.forEach((line,i)=>{if(i)fragment.append(document.createElement('br'));fragment.append(document.createTextNode(line));});const last=fragment.lastChild;range.insertNode(fragment);if(last){range.setStartAfter(last);range.collapse(true);const selection=window.getSelection();selection.removeAllRanges();selection.addRange(range);}changed(root);
    });
  }
  function refreshCommands(id){const root=rootById(id);if(root)updateSlash(root);}
  function focusBlock(id,offset=0){if(window.ReaderPageNavigation?.ensureBlock?.(id)===false)return false;const b=api.getBlock?.(id);if(b?.type==='equation'){const atom=[...document.querySelectorAll('.source-block')].find(n=>n.dataset.blockId===id)?.querySelector('.direct-math-atom');openMath({root:null,node:atom,tex:b.source||'',id,displayBlock:true});return;}caret(id,offset);}
  function init(options) {
    api=options||{};if(initialized)return;initialized=true;
    document.addEventListener('keydown',handleKeydown,true);
    document.addEventListener('pointerdown',event=>{
      verticalGoalX=null;
      if(mathEditor&&!mathEditor.panel.contains(event.target)&&!event.target.closest('.direct-math-atom')){if(!commitMath()){event.preventDefault();return;}}
      if(slash&&!slash.panel.contains(event.target)&&!slash.root.contains(event.target))hideSlash();
    },true);
    window.addEventListener('beforeunload',flush);
  }
  window.ReaderDirectEdit={init,bind,isEditing,flush,preparePageChange,focusBlock,navigateToBlock,caret,refreshCommands,handleKeydown,snapshotSelection,restoreSelection,preserveSelection,testing:{serializeEditable,serializeNode,extractFormatting,linearProjection,texOf,edgeGeometry}};
})();
