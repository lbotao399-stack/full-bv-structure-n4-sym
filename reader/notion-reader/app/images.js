/* Clipboard pictures are ordinary persistent image blocks; audio is independent. */
(function(){'use strict';
 const MAX_BYTES=16*1024*1024,MAX_PIXELS=32000000;
 let api={},initialized=false,endpoint=null,queue=Promise.resolve(),pendingCopy=null,preparedCopy=null;
 const pendingAnchors=new Map();
 const uid=()=> 'local-'+crypto.randomUUID();
 const source=b=>typeof b.image==='string'?b.image:(b.image?.src||b.source);
 const absolute=s=>new URL(s,new URL('../',document.baseURI)).href;
 const selected=()=>api.getSelectedId?.();
 const get=id=>api.getBlock?.(id);
 const escape=s=>String(s||'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
 const fail=e=>api.toast?.(e?.message||'The image could not be inserted.');
 function clipboardHTML(b,pixels){return `<figure data-reader-image-copy="1"><img src="${escape(pixels||absolute(source(b)))}" alt="${escape(b.image?.alt||b.alt||'')}"><figcaption>${escape(b.caption||'')}</figcaption></figure>`;}
 function clipboardImage(html){
  if(!html)return null;const parsed=new DOMParser().parseFromString(html,'text/html');
  const figure=parsed.querySelector('figure[data-reader-image-copy="1"]'),image=figure?.querySelector('img')||parsed.querySelector('img');
  if(!image)return null;
  if(!figure){const copy=parsed.body.cloneNode(true);copy.querySelectorAll('img,meta,style,script').forEach(n=>n.remove());if(copy.textContent.trim())return null;}
  const src=image.getAttribute('src')||'';
  if(!/^(?:https?:\/\/|data:image\/(?:png|jpeg|webp|gif|bmp);base64,)/i.test(src))return null;
  return {src,alt:(image.getAttribute('alt')||'').slice(0,1000),caption:figure?.querySelector('figcaption')?.textContent.slice(0,100000)||''};
 }
 function readDataURL(blob){return new Promise((resolve,reject)=>{const reader=new FileReader();reader.onload=()=>resolve(reader.result);reader.onerror=()=>reject(Error('The clipboard image could not be read.'));reader.readAsDataURL(blob);});}
 async function png(blob){
  if(!blob?.size)throw Error('The clipboard contains an empty image.');
  if(blob.size>MAX_BYTES)throw Error('Choose an image smaller than 16 MB.');
  if(!/^image\/(?:png|jpeg|jpg|webp|gif|bmp|avif)$/i.test(blob.type))throw Error('Paste a PNG, JPEG, WebP, GIF, BMP or AVIF image.');
  const url=URL.createObjectURL(blob),image=new Image();
  try{image.src=url;await image.decode();const width=image.naturalWidth,height=image.naturalHeight;
   if(!width||!height||width*height>MAX_PIXELS)throw Error('Choose an image with at most 32 million pixels.');
   const canvas=document.createElement('canvas');canvas.width=width;canvas.height=height;
   const ctx=canvas.getContext('2d');if(!ctx)throw Error('Image decoding is unavailable.');ctx.drawImage(image,0,0);
   const result=await new Promise(resolve=>canvas.toBlob(resolve,'image/png'));
   if(!result||result.size>MAX_BYTES)throw Error('The decoded image exceeds 16 MB; use a smaller image.');return result;
  }catch(e){if(e.name==='EncodingError')throw Error('This image format could not be decoded. Try copying it as PNG.');throw e;}
  finally{URL.revokeObjectURL(url);}
 }
 async function imageBlob(src){const r=await fetch(src);if(!r.ok)throw Error('The copied image could not be loaded.');return r.blob();}
 function prepareCopy(b){
  const url=absolute(source(b));if(preparedCopy?.url===url)return preparedCopy;
  const record={url,blob:imageBlob(url).then(png),data:null};
  record.ready=record.blob.then(readDataURL).then(value=>{record.data=value;return value;});record.ready.catch(()=>{});preparedCopy=record;return record;
 }
 async function uploadEndpoint(){
  if(endpoint)return endpoint;
  for(const url of ['/api/images']){
   try{const r=await fetch(url,{cache:'no-store'});if(r.ok&&r.headers.get('content-type')?.includes('application/json')){const d=await r.json();if(d.supported===true){endpoint=url;return url;}}}catch(_e){}
  }
  throw Error('Image storage is unavailable. Start the reader image service, then paste again.');
 }
 async function upload(blob){const url=await uploadEndpoint();let response;
  try{response=await fetch(url,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({dataUrl:await readDataURL(blob)})});}catch(_e){endpoint=null;throw Error('Image storage disconnected. Paste again after restarting the reader.');}
  const result=await response.json();if(!response.ok)throw Error(result.error||'The image could not be saved.');
  if(!/^figures\/pasted\/[a-f0-9]+\.png$/.test(result.source||''))throw Error('Image storage returned an invalid file path.');return result;
 }
 function add(items,anchor=selected()){
  if(!get(anchor))return Promise.resolve();window.ReaderDirectEdit?.flush();
  const predecessor=pendingAnchors.get(anchor);
  const run=async()=>{
   const preceding=predecessor?await predecessor.catch(()=>null):null;
   const at=preceding?.at(-1)?.id||anchor;
   if(!get(at))throw Error('The insertion block was removed. Select a passage and paste again.');
   api.toast?.('Saving image…');const blocks=[];
   for(const item of items){const blob=await png(item.blob||await imageBlob(item.src));const saved=await upload(blob);
    blocks.push({id:uid(),type:'image',source:saved.source,alt:item.alt||'',caption:item.caption||''});
   }
   if(!get(at))throw Error('The insertion block was removed. Select a passage and paste again.');
   window.ReaderDirectEdit?.flush();api.beforeChange?.();api.insertAfter?.(at,blocks,{focus:selected()===anchor||selected()===at});
   api.toast?.(blocks.length===1?'Image inserted. Click its caption to edit.':`${blocks.length} images inserted.`);return blocks;
  };
  const promise=queue.then(run);pendingAnchors.set(anchor,promise);
  promise.finally(()=>{if(pendingAnchors.get(anchor)===promise)pendingAnchors.delete(anchor);}).catch(()=>{});
  queue=promise.catch(fail);return promise;
 }
 function isPasteTarget(target){return !document.body.classList.contains('floating')&&!target.closest?.('input,textarea,select,.direct-math-editor,.editor')&&(target.closest?.('#manuscript')||target===document.body);}
 function paste(event){
  if(!isPasteTarget(event.target)||!get(selected()))return;
  const data=event.clipboardData;if(!data)return;
  const files=[...data.items||[]].filter(item=>item.kind==='file'&&item.type.startsWith('image/')).map(item=>item.getAsFile()).filter(Boolean);
  const copied=clipboardImage(data.getData('text/html'));
  if(!files.length&&!copied)return;
  event.preventDefault();event.stopImmediatePropagation();
  const anchor=event.target.closest?.('.source-block')?.dataset.blockId||selected();
  const items=files.length?files.map((blob,i)=>({blob,alt:i===0?copied?.alt||'':'',caption:i===0?copied?.caption||'':''})):[copied];
  void add(items,anchor).catch(()=>{});
 }
 function copyEvent(event){const b=pendingCopy?.block||get(document.activeElement?.dataset.readerImage);if(!b||b.type!=='image'||!event.clipboardData)return;
  let pixels=pendingCopy?.data||(preparedCopy?.url===absolute(source(b))?preparedCopy.data:null);
  if(!pixels){
   // A native copy event is synchronous, including immediately after image focus.
   const image=[...document.querySelectorAll('img[data-reader-image]')].find(n=>n.dataset.readerImage===b.id);
   try{if(image?.complete&&image.naturalWidth&&image.naturalWidth*image.naturalHeight<=MAX_PIXELS){
    const canvas=document.createElement('canvas');canvas.width=image.naturalWidth;canvas.height=image.naturalHeight;canvas.getContext('2d').drawImage(image,0,0);
    const data=canvas.toDataURL('image/png');if(data.length<=4*Math.ceil(MAX_BYTES/3)+22)pixels=data;
   }}catch(_e){}
  }
  if(!pixels){event.preventDefault();api.toast?.('The image is preparing. Press Command+C again, or use Copy image in its menu.');return;}
  event.preventDefault();event.clipboardData.setData('text/html',clipboardHTML(b,pixels));event.clipboardData.setData('text/plain',absolute(source(b)));
 }
 async function copy(id){
  const b=get(id);if(b?.type!=='image')return;
  const record=prepareCopy(b),text=absolute(source(b));
  try{
   if(navigator.clipboard?.write&&window.ClipboardItem){
    // Promise-valued PNG keeps clipboard permission tied to this user gesture.
    await navigator.clipboard.write([new ClipboardItem({'image/png':record.blob,'text/html':record.ready.then(data=>new Blob([clipboardHTML(b,data)],{type:'text/html'})),'text/plain':new Blob([text],{type:'text/plain'})})]);
   }else{pendingCopy={block:b,data:await record.ready};let ok;try{ok=document.execCommand('copy');}finally{pendingCopy=null;}if(!ok)throw Error('Click the image and press Command+C to copy it.');}
   api.toast?.('Image copied. Paste it into any reader.');
  }catch(e){fail(e.name==='NotAllowedError'?Error('Use the native shortcut instead: click the image and press Command+C.'):e);}
 }
 function choose(anchor=selected()){
  if(!get(anchor))return;const input=document.createElement('input');input.type='file';input.accept='image/png,image/jpeg,image/webp,image/gif,image/bmp,image/avif';input.multiple=true;input.hidden=true;document.body.append(input);
  input.onchange=()=>{const items=[...input.files].map(blob=>({blob}));input.remove();if(items.length)void add(items,anchor).catch(()=>{});};input.oncancel=()=>input.remove();input.click();
 }
 function bind(b,wrap){
  if(b.type!=='image')return;const image=wrap.querySelector('.block-content img');if(!image)return;
  image.tabIndex=0;image.dataset.readerImage=b.id;image.setAttribute('aria-label',(b.alt||'Image')+' · Drag to move this block · Command+C to copy');
  image.addEventListener('focus',()=>prepareCopy(b));
  image.addEventListener('click',e=>{e.stopPropagation();getSelection()?.removeAllRanges();api.select?.(b.id);image.focus({preventScroll:true});prepareCopy(b);});
 }
 function addMenuItems(menu,id,button){if(get(id)?.type==='image')menu.append(button('Copy image',()=>void copy(id)));menu.append(button('Insert image below…',()=>choose(id)));}
 function init(options){api=options;if(initialized)return;initialized=true;document.addEventListener('paste',paste,true);document.addEventListener('copy',copyEvent,true);
  const toolbar=document.querySelector('.toolbar-actions');if(toolbar){const b=document.createElement('button');b.id='add-image';b.className='quiet';b.textContent='+ Image';b.onclick=()=>choose();toolbar.append(b);}
 }
 window.ReaderImages={init,bind,addMenuItems,copy,choose,testing:{clipboardHTML,clipboardImage,png,upload,add}};
})();
