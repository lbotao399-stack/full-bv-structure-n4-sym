/* Remove spoken position cues from labels, preserving the original math tree. */
const {DOMParser, XMLSerializer} = require('@xmldom/xmldom');
const namespace = 'http://www.w3.org/1998/Math/MathML';
const children = node => Array.from(node.childNodes || []).filter(n => n.nodeType === 1);
function label(node) {
  if (!node) return false;
  const kind = node.localName, parts = children(node);
  if (kind === 'mi') return /^[\p{L}]+$/u.test(node.textContent.trim());
  if (kind === 'mn') return /^\d+$/.test(node.textContent.trim());
  if (kind === 'mo') return node.textContent.trim() === ',';
  if (kind === 'mspace') return true;
  if (kind === 'mover') return label(parts[0]) && /^[ˉ¯‾]$/.test(parts[1]?.textContent.trim());
  return ['mrow','mstyle'].includes(kind) && parts.length > 0 && parts.every(label);
}
module.exports = function naturalIndices(math, speak) {
  const document = new DOMParser().parseFromString(math, 'application/xml');
  const serializer = new XMLSerializer(), positions = new Map(), scripts = [];
  for (const kind of ['msub','msup','msubsup']) for (const node of Array.from(document.getElementsByTagName(kind))) {
    const parts = children(node), base = parts[0]?.textContent || '';
    if (!positions.has(base)) positions.set(base, new Set());
    positions.get(base).add(kind);
    if (kind !== 'msup') scripts.push({node, parts, base, kind});
  }
  const replacements = [];
  for (const {node, parts, base, kind} of scripts) {
    const original = speak('<math xmlns="'+namespace+'">'+serializer.serializeToString(node)+'</math>');
    if (!/\b(?:subscript|sub)\b/i.test(original)) continue;
    const compact = kind === 'msub' && label(parts[1]) && positions.get(base).size === 1;
    const marker = compact ? '' : label(parts[1]) ? 'lower ' : 'index ';
    let replacement = original.replace(/\b(?:subscript|sub)\s*/gi, marker);
    if (!compact && !label(parts[1]) && kind === 'msub') replacement += ', end index';
    replacements.push({original, replacement});
  }
  // Match generated subexpression speech, never delete words globally from
  // prose. Rendering the unchanged full tree retains brackets and operator roles.
  let text = speak(math);
  for (const {original, replacement} of replacements.sort((a,b) => b.original.length-a.original.length)) text = text.split(original).join(replacement);
  return text;
};
