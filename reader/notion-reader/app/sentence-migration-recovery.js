/* Rebase a pre-migration local draft only over the exact sentence-only revision. */
(function(root){'use strict';
 const clone=value=>JSON.parse(JSON.stringify(value));
 const canonical=value=>JSON.stringify(sort(value));
 function sort(value){if(Array.isArray(value))return value.map(sort);if(value&&typeof value==='object')return Object.fromEntries(Object.keys(value).sort().map(k=>[k,sort(value[k])]));return value;}
 function recover(saved,revision,draft,receipt){
  if(!draft?.dirty||!Array.isArray(draft.document?.blocks)||!receipt||receipt.version!==1)return null;
  if(draft.document.id!==saved.id||receipt.before?.document?.id!==saved.id||receipt.after?.document?.id!==saved.id)return null;
  if(draft.revision!==receipt.before.revision||revision!==receipt.after.revision||revision!==draft.revision+1)return null;
  if(canonical(saved)!==canonical(receipt.after.document))return null;
  const old=new Map(receipt.before.document.blocks.map(b=>[b.id,b]));
  const current=new Map(saved.blocks.map(b=>[b.id,b]));
  const result=clone(draft.document);const ids=new Set();result.blocks=[];
  for(const b of draft.document.blocks){
   const mapped=receipt.fragments?.[b.id];
   const unchanged=old.has(b.id)&&canonical(b)===canonical(old.get(b.id));
   const pieces=unchanged&&Array.isArray(mapped)&&mapped.length?mapped.map(id=>current.get(id)):[b];
   if(pieces.some(p=>!p))return null;
   for(const p of pieces){if(ids.has(p.id))return null;ids.add(p.id);result.blocks.push(clone(p));}
  }
  if(!result.blocks.length)return null;
  return result;
 }
 async function recoverDraft(saved,revision,draft){
  if(!draft?.dirty||draft.revision===revision)return null;
  try{const response=await fetch(new URL('sentence-migration.json',document.baseURI),{cache:'no-store'});if(!response.ok)return null;return recover(saved,revision,draft,await response.json());}catch(_e){return null;}
 }
 const api={recoverDraft,recover};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.ReaderSentenceMigration=api;
})(typeof window!=='undefined'?window:globalThis);
