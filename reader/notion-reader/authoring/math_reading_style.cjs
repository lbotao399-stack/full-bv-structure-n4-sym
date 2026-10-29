/* Speech-only conventions. Displayed TeX is never changed. */
const {DOMParser} = require('@xmldom/xmldom');
const digits = {
  'en-US': ['zero','one','two','three','four','five','six','seven','eight','nine'],
  'zh-CN': ['零','一','二','三','四','五','六','七','八','九']
};
function referenceNumber(value, language) {
  const names = digits[language] || digits['en-US'];
  return String(value).split(/([-.])/).map(part => part === '-' ? (language === 'zh-CN' ? '，' : ', ')
    : part === '.' ? (language === 'zh-CN' ? ' 点 ' : ' point ')
    : [...part].map(c => names[Number(c)] || c).join(language === 'zh-CN' ? '' : ' ')).join('');
}
function references(text, language) {
  return text.replace(/\b(equations?|eq\.?|theorem|proposition|lemma|corollary|definition|figure|fig\.?|table|section|sec\.?)\s*\(?(\d+(?:[-.]\d+)+)\)?/gi,
    (_, label, number) => ({eq:'Equation',fig:'Figure',sec:'Section'}[label.toLowerCase().replace('.','')]||label) + ' ' + referenceNumber(number, 'en-US'))
    .replace(/(公式|方程|定理|命题|引理|推论|定义|式|图|表|节)\s*[（(]?(\d+(?:[-.]\d+)+)[）)]?/g,
    (_, label, number) => label + ' ' + referenceNumber(number, 'zh-CN'));
}
function elements(node) {
  return Array.from({length: node?.childNodes?.length || 0}, (_, i) => node.childNodes[i]).filter(n => n.nodeType === 1);
}
function scriptContext(math) {
  const root = typeof math === 'string' ? new DOMParser().parseFromString(math, 'application/xml').documentElement : math;
  const positions = new Map(), subscripts = [];
  function visit(node) {
    const kind = node.localName || node.nodeName, c = elements(node);
    if (['msub','msup','msubsup'].includes(kind)) {
      const base = c[0]?.textContent || '';
      if (!positions.has(base)) positions.set(base, new Set());
      positions.get(base).add(kind);
      if (kind === 'msub') subscripts.push({base, index: c[1]});
    }
    c.forEach(visit);
  }
  visit(root);
  const ambiguous = new Set([...positions].filter(([, kinds]) => kinds.size > 1 || kinds.has('msubsup')).map(([base]) => base));
  const simple = node => /^[A-Za-zα-ωΑ-Ω0-9]$/.test((node?.textContent || '').trim());
  return {root, ambiguous, simple,
    compactAllSubscripts: subscripts.length > 0 && subscripts.every(x => simple(x.index) && !ambiguous.has(x.base))};
}
module.exports = {referenceNumber, references, scriptContext};
