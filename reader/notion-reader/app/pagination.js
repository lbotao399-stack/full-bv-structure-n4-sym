/* Whole source chapters. The complete document remains the save/export model. */
(function(root,factory){
 'use strict';const api=factory();if(typeof module==='object'&&module.exports)module.exports=api;if(root)root.ReaderPagination=api;
})(typeof window!=='undefined'?window:null,function(){
 'use strict';
 function title(block){return String(block?.source||'').replace(/<[^>]*>/g,'').replace(/\s+/g,' ').trim();}
 function partition(blocks,{chapterLevel=1,documentTitle='',frontMatterLabel='Front matter'}={}){
  if(!Array.isArray(blocks)||!blocks.length)return [];
  chapterLevel=Number(chapterLevel);if(!Number.isInteger(chapterLevel)||chapterLevel<1||chapterLevel>6)chapterLevel=1;
  const firstHeading=blocks.findIndex(b=>b.type==='heading'&&!b.parentId);
  const normalizedTitle=title({source:documentTitle});
  const documentHeading=firstHeading>=0&&(blocks[firstHeading].sourceLabel==='title'||normalizedTitle&&title(blocks[firstHeading])===normalizedTitle)?firstHeading:-1;
  // Subheadings stay inside their chapter. Container headings cannot split it.
  const starts=[];
  for(let i=0;i<blocks.length;i++){
   const b=blocks[i];
   if(b.type==='heading'&&!b.parentId&&Number(b.level||2)===chapterLevel&&i!==documentHeading)starts.push(i);
  }
  if(!starts.length)return [{start:0,end:blocks.length,anchor:blocks[0].id,label:normalizedTitle||(documentHeading>=0?title(blocks[documentHeading]):'')||'Manuscript'}];
  const frontMatter=starts[0]>0;if(frontMatter)starts.unshift(0);
  // Consecutive chapter headings delimit separate units, even an empty appendix.
  return starts.map((start,i)=>({start,end:starts[i+1]??blocks.length,anchor:blocks[start].id,label:frontMatter&&i===0?frontMatterLabel:title(blocks[start])||'Untitled chapter'}));
 }
 function pageFor(pages,blocks,id){const at=blocks.findIndex(b=>b.id===id);if(at<0)return -1;return pages.findIndex(p=>at>=p.start&&at<p.end);}
 return {partition,pageFor};
});
