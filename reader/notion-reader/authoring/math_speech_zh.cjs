/* Structural Mandarin reading of KaTeX MathML; no translation of source prose.
 * row/mew are speech-only aliases for the English names rho/mu in Xiaoxiao clauses. */
const {scriptContext, referenceNumber}=require('./math_reading_style.cjs');
const words={
 'α':'阿尔法','β':'贝塔','γ':'伽马','δ':'德尔塔','ε':'艾普西隆','ϵ':'艾普西隆','ζ':'泽塔','η':'伊塔','θ':'西塔','ϑ':'西塔','ι':'约塔','κ':'卡帕','λ':'lambda','μ':'mew','ν':'纽','ξ':'克西','π':'派','ρ':'row','ϱ':'row','σ':'西格玛','τ':'陶','υ':'宇普西隆','φ':'斐','ϕ':'斐','χ':'凯','ψ':'普西','ω':'欧米伽',
 'Γ':'伽马','Δ':'德尔塔','Θ':'西塔','Λ':'兰姆达','Ξ':'克西','Π':'派','Σ':'西格玛','Φ':'斐','Ψ':'普西','Ω':'欧米伽',
 '+':'加','−':'减','-':'减','±':'正负','∓':'负正','=':'等于','≠':'不等于','<':'小于','>':'大于','≤':'小于等于','≥':'大于等于','≡':'恒等于','≅':'同构于','≈':'约等于','∼':'相似于',
 '×':'乘','⋅':'乘','·':'乘','⊗':'张量积','⊕':'直和','∈':'属于','∉':'不属于','⊂':'包含于','⊆':'包含于或等于','⊃':'包含','∪':'并','∩':'交','∅':'空集','∞':'无穷',
 '∑':'求和','∏':'连乘','∫':'积分','∮':'闭路积分','∂':'偏导','∇':'纳布拉','→':'趋于','↦':'映到','⇒':'推出','⇔':'当且仅当','↔':'双向对应','∀':'任意','∃':'存在','¬':'非','∧':'且','∨':'或',
 '(':'左括号',')':'右括号','[':'左方括号',']':'右方括号','{':'左花括号','}':'右花括号','|':'竖线','∣':'竖线','‖':'双竖线','⟨':'左尖括号','⟩':'右尖括号',
 ',':'，',';':'；',':':'：','.':'。','…':'省略号','⋯':'省略号','⋮':'竖省略号','⋱':'斜省略号','!':'阶乘','′':'撇','″':'双撇','*':'星号','∗':'星号','†':'厄米共轭','□':'证毕','■':'证毕',
 'sin':'正弦','cos':'余弦','tan':'正切','cot':'余切','sinh':'双曲正弦','cosh':'双曲余弦','exp':'指数函数','log':'对数','ln':'自然对数','lim':'极限','det':'行列式','tr':'迹','Tr':'迹','dim':'维数','ker':'核','Hom':'态射空间','Re':'实部','Im':'虚部',
 'irreps':'不等价不可约表示','classes':'共轭类','counting':'计数测度','states':'态','cosets':'陪集','weights':'权重'
};
function token(text){
 const s=String(text||'').trim();if(/^[\u2061-\u2064\u200b]*$/.test(s))return '';if(words[s])return words[s];
 return [...s.normalize('NFKD')].map(c=>words[c]||c).join(' ');
}
function children(n){return [...Array(n.childNodes.length)].map((_,i)=>n.childNodes[i]).filter(x=>x.nodeType===1);}
function speech(n,context){
 if(!n)return '';
 const c=children(n),name=n.localName||n.nodeName;
 const read=x=>speech(x,context);
 const join=()=>{
  const bar=x=>(x.localName||x.nodeName)==='mo'&&/^[|∣]$/.test(x.textContent.trim());
  const paired=c.filter(bar).length>0&&c.filter(bar).length%2===0;
  const parts=[];
  for(let i=0;i<c.length;i++){
   if(paired&&bar(c[i])){
    let end=i+1;while(end<c.length&&!bar(c[end]))end++;
    const inner=c.slice(i+1,end).map(read).filter(Boolean).join(' ');
    parts.push(c.slice(i+1,end).length===1?`${inner}的绝对值`:`绝对值，${inner}，绝对值结束`);i=end;
   }else parts.push(read(c[i]));
  }
  return parts.filter(Boolean).join(' ');
 };
 switch(name){
  case 'annotation':case 'annotation-xml':case 'mspace':case 'mphantom':return '';
  case 'semantics':return read(c[0]);
  case 'mi':case 'mn':case 'mo':return token(n.textContent);
  case 'mtext':return n.textContent;
  case 'mfrac':return `分式，分子 ${read(c[0])}，分母 ${read(c[1])}，分式结束`;
  case 'msqrt':return `平方根，${join()}，根号结束`;
  case 'mroot':return `${read(c[1])} 次方根，${read(c[0])}，根号结束`;
  case 'msub':return context.simple(c[1])&&!context.ambiguous.has(c[0]?.textContent)?`${read(c[0])} ${read(c[1])}`:`${read(c[0])} 下标 ${read(c[1])}`;
  case 'msup':return c[1]?.textContent==='2'?`${read(c[0])} 的平方`:c[1]?.textContent==='3'?`${read(c[0])} 的立方`:`${read(c[0])} 上标 ${read(c[1])}`;
  case 'msubsup':return `${read(c[0])} 下标 ${read(c[1])}，上标 ${read(c[2])}`;
  case 'munder':return `${read(c[0])}，下方 ${read(c[1])}，`;
  case 'mover':{
   const accent=c[1]?.textContent;const label={'^':'帽','ˆ':'帽','~':'波浪号','˜':'波浪号','¯':'横线','→':'矢量箭头','˙':'点','¨':'双点'}[accent];
   return label?`${read(c[0])} 上加${label}`:`${read(c[0])}，上方 ${read(c[1])}，`;
  }
  case 'munderover':return `${read(c[0])}，下限 ${read(c[1])}，上限 ${read(c[2])}，`;
  case 'mtable':{
   // KaTeX wraps tagged equations in a full-width one-row layout table.
   if(n.getAttribute('width')==='100%'&&c.length===1){
    const cells=children(c[0]).filter(x=>x.textContent.trim());
    const last=cells.at(-1),tag=last?.textContent.trim();
    if(cells.length>1&&/^\([\d.-]+\)$/.test(tag))return `式 ${referenceNumber(tag.slice(1,-1),'zh-CN')}，${cells.slice(0,-1).map(read).filter(Boolean).join(' ')}`;
    return cells.map(read).filter(Boolean).join(' ');
   }
   const aligned=/right left/.test(n.getAttribute('columnalign')||'');
   return c.map((r,i)=>(aligned?'':`第 ${i+1} 行，`)+read(r)).join('；');
  }
  case 'mtr':case 'mlabeledtr':return c.map(read).filter(Boolean).join(' ');
  case 'menclose':return `括起 ${join()}`;
  default:return c.length?join():token(n.textContent);
 }
}
module.exports=mathml=>{const context=scriptContext(mathml);return speech(context.root,context).replace(/\s+/g,' ').replace(/\s*([，；。])\s*/g,'$1').replace(/，；/g,'；').replace(/；。/g,'。').trim();};
