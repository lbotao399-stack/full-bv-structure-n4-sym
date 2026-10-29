/* Share pending narration work without sharing a caller's cancellation. */
(function(root,factory){
 const api=factory();if(typeof module==='object'&&module.exports)module.exports=api;
 if(root)root.ReaderSharedPreparation=api;
})(typeof window!=='undefined'?window:typeof globalThis!=='undefined'?globalThis:null,function(){'use strict';
 function cancelled(){if(typeof DOMException==='function')return new DOMException('Narration preparation cancelled.','AbortError');const error=new Error('Narration preparation cancelled.');error.name='AbortError';return error;}
 function create({request,keyOf,onPromote}){
  const pending=new Map();
  function forget(entry){if(pending.get(entry.key)===entry)pending.delete(entry.key);}
  function finish(consumer,ok,value){
   if(consumer.done)return;consumer.done=true;
   consumer.signal?.removeEventListener('abort',consumer.abort);consumer.entry.consumers.delete(consumer);
   if(ok)consumer.resolve(value);else consumer.reject(value);
  }
  function settle(entry,ok,value){
   entry.settled=true;forget(entry);entry.promotion?.abort();
   for(const consumer of [...entry.consumers])finish(consumer,ok,value);
  }
  function start(block,key,options){
   const entry={block,key,options,controller:new AbortController(),consumers:new Set(),settled:false,promoted:false,promotion:null};
   // Attach rejection handling immediately, including synchronous request failures.
   Promise.resolve().then(()=>{
    if(entry.controller.signal.aborted)throw cancelled();
    return request(entry.block,entry.controller.signal,entry.options);
   }).then(value=>settle(entry,true,value),error=>settle(entry,false,error));
   return entry;
  }
  function promote(entry){
   if(!entry.options.background||entry.promoted||entry.settled||entry.controller.signal.aborted)return;
   entry.promoted=true;if(typeof onPromote!=='function')return;
   // Promotion has its own cancellation; its failure never invalidates the recording.
   const controller=new AbortController();entry.promotion=controller;
   Promise.resolve().then(()=>{
    if(controller.signal.aborted||entry.settled)return;
    return onPromote(entry.block,controller.signal);
   }).catch(()=>{});
  }
  function prepare(block,signal,{background=false,refresh=false}={}){
   if(signal?.aborted)return Promise.reject(cancelled());
   let key;try{key=keyOf(block);}catch(error){return Promise.reject(error);}
   let entry=refresh?null:pending.get(key);
   if(!entry||entry.settled||entry.controller.signal.aborted){
    entry=start(block,key,{background,refresh});if(!refresh)pending.set(key,entry);
   }else if(!background)promote(entry);
   return new Promise((resolve,reject)=>{
    const consumer={entry,signal,resolve,reject,done:false,abort:null};
    consumer.abort=()=>{
     finish(consumer,false,cancelled());
     if(!entry.settled&&!entry.consumers.size){forget(entry);entry.controller.abort();entry.promotion?.abort();}
    };
    entry.consumers.add(consumer);signal?.addEventListener('abort',consumer.abort,{once:true});
    if(signal?.aborted)consumer.abort();
   });
  }
  return {prepare};
 }
 return {create};
});
