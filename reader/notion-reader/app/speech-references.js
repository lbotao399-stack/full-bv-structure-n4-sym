/* Source-aware spoken references, shared by planning and measured cue matching. */
(function(root,factory){const api=factory();if(typeof module==='object'&&module.exports)module.exports=api;if(root)root.ReaderSpeechReferences=api;})(typeof window!=='undefined'?window:null,function(){
  'use strict';
  const number='\\d+(?:[.-]\\d+)+[a-z]?';
  const item='[（(]?'+number+'[)）]?';
  const names={eq:'Equation',eqs:'Equations',equation:'Equation',equations:'Equations',fig:'Figure',figs:'Figures',figure:'Figure',figures:'Figures',sec:'Section',secs:'Sections',sect:'Section',sects:'Sections',section:'Section',sections:'Sections','§':'Section','§§':'Sections',table:'Table',tables:'Tables',theorem:'Theorem',lemma:'Lemma',proposition:'Proposition',definition:'Definition',corollary:'Corollary'};
  const explicit=new RegExp('(?:\\b(eq(?:uation)?s?|fig(?:ure)?s?|sec(?:tion|t)?s?|tables?|theorem|lemma|proposition|definition|corollary)\\.?|([§]{1,2}))\\s*('+item+'(?:\\s*(?:,|and|or|to|[–—])\\s*'+item+')*)','gi');
  function labels(document){
    const found=new Set();
    for(const block of document?.blocks||[])for(const m of String(block.source||'').matchAll(/\\tag\*?\s*\{\s*([^{}]+?)\s*\}/g))found.add(m[1]);
    return [...found];
  }
  function spoken(number,language){return number.replace(/\./g,language==='zh-CN'?' 点 ':' point ').replace(/-/g,language==='zh-CN'?'，':', ');}
  function normalize(value,{equationLabels=[],language='en-US'}={}){
    const saved=[],known=new Set(equationLabels),hold=value=>{const id=saved.length;saved.push(value);return ' READERREFPROTECTED'+id+'END ';};
    let text=String(value||'').replace(/https?:\/\/[^\s<>]+|`[^`]*`|\$\$[\s\S]*?\$\$|\$[^$\n]*\$|\\\[[\s\S]*?\\\]|\\\([\s\S]*?\\\)/g,hold);
    text=text.replace(explicit,(_whole,label,symbol,values)=>{
      const name=names[(label||symbol).toLowerCase()];
      return hold(name+' '+values.replace(/[()（）]/g,'').replace(new RegExp(number,'gi'),value=>spoken(value,language)).replace(/[–—]/g,' to '));
    });
    text=text.replace(new RegExp('[（(]('+number+')[)）]','gi'),(whole,value,offset)=>{
      if(!known.has(value))return whole;
      // A known equation identifier still need not be a citation in a numeric example.
      const before=text.slice(Math.max(0,offset-70),offset),after=text.slice(offset+whole.length);
      if(/\b(?:value|number|decimal|measurement|temperature|mass|voltage|length|width|height|pH|version)\s*$/i.test(before)||/^\s*(?:kg|cm|mm|km|Hz|MHz|GHz|°C|volts?|seconds?)\b/.test(after))return whole;
      return (language==='zh-CN'?'式 ':'Equation ')+spoken(value,language);
    });
    return text.replace(/\s*READERREFPROTECTED(\d+)END\s*/g,(_m,id)=>' '+saved[Number(id)]+' ').replace(/[ \t]+/g,' ').replace(/\s+([,.;:!?])/g,'$1').trim();
  }
  return {version:'source-references-v1',labels,normalize};
});
