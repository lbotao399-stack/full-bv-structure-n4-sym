/* Resolve reviewed prose annotations before mathematical speech, without editing TeX. */
const {DOMParser,XMLSerializer}=require('@xmldom/xmldom');
const speechText=require('../app/speech-text.js');

function groupAt(tex,start){
  let depth=0;
  for(let i=start;i<tex.length;i++){
    if(tex[i]==='\\'){i++;continue;}
    if(tex[i]==='{')depth++;
    else if(tex[i]==='}'&&--depth===0)return {body:tex.slice(start+1,i),end:i+1};
  }
  throw Error('Unclosed formula annotation.');
}
function tagOf(tex){
  const match=/\\tag\*?\s*\{/.exec(tex);if(!match)return null;
  const group=groupAt(tex,match.index+match[0].length-1);
  return {label:group.body.trim(),raw:tex.slice(match.index,group.end),expression:tex.slice(0,match.index)+tex.slice(group.end)};
}
function romanLabel(body){
  const match=/^\s*\(\s*([A-Za-z\s~]+)\s*\)\s*$/.exec(body);if(!match)return null;
  const words=match[1].split(/~+/).flatMap(part=>{
    const tokens=part.trim().split(/\s+/).filter(Boolean);
    return tokens.length>1&&tokens.every(token=>token.length===1)?[tokens.join('')]:tokens;
  });
  // This is a review candidate, not permission to turn arbitrary roman products into words.
  return words.length>=2&&words.every(word=>word.length>=3)?words.join(' ').toLowerCase():null;
}
function textAnnotations(tex,labels={}){
  let result='',at=0;const command=/\\(text|textrm|textnormal|textsf|texttt|textbf|textit|mathrm)\s*\{/g;let match;
  while((match=command.exec(tex))){
    const group=groupAt(tex,match.index+match[0].length-1);
    result+=tex.slice(at,match.index);
    if(match[1]!=='mathrm'){
      // Text commands already declare prose; normalize prose hyphens before SRE sees them.
      result+=tex.slice(match.index,match.index+match[0].length)+speechText.normalize(group.body)+'}';
    }else{
      const label=romanLabel(group.body);
      if(label){
        const reviewed=labels[label];
        if(typeof reviewed!=='string'||reviewed.toLowerCase()!==label||!/^[A-Za-z ]+$/.test(reviewed))throw Error('Review the prose label "'+label+'" in formula-labels.json or supply a source-bound formula reading before narration.');
        result+='\\text{('+speechText.normalize(reviewed)+')}';
      }else result+=tex.slice(match.index,group.end);
    }
    at=group.end;command.lastIndex=group.end;
  }
  return result+tex.slice(at);
}
function normalizeMathText(math){
  const doc=new DOMParser().parseFromString(math,'application/xml');
  for(const node of Array.from(doc.getElementsByTagName('mtext')))node.textContent=speechText.normalize(node.textContent);
  return new XMLSerializer().serializeToString(doc.documentElement);
}
module.exports={tagOf,textAnnotations,normalizeMathText};
