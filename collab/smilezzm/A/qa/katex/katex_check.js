// Render every item of items.json with the vendored KaTeX: throwOnError true, strict 'ignore',
// displayMode for display items.  Run: jsc katex_check.js  (macOS JavaScriptCore) or node katex_check.js
var isNode = (typeof require !== 'undefined');
var root = isNode ? require('path').resolve(__dirname, '../../../../../') : '/Users/okongoyango/Projects/full-bv-structure-n4-sym/';
var katexPath = root + '/reader/notion-reader/app/vendor/katex/katex.min.js';
var itemsPath = (isNode ? __dirname : '/Users/okongoyango/Projects/full-bv-structure-n4-sym/collab/smilezzm/A/qa/katex') + '/items.json';
var katex, items;
if (isNode) { katex = require(katexPath); items = JSON.parse(require('fs').readFileSync(itemsPath, 'utf8')); }
else { load(katexPath); katex = this.katex || globalThis.katex; items = JSON.parse(readFile(itemsPath)); }
var fails = [], n = 0;
for (var i = 0; i < items.length; i++) {
  var it = items[i]; n++;
  try { katex.renderToString(it[1], {displayMode: it[0] === 'display', throwOnError: true, strict: 'ignore'}); }
  catch (e) { fails.push({tag: it[2], kind: it[0], tex: it[1].slice(0, 200), error: String(e.message || e)}); }
}
print(JSON.stringify({rendered: n, failures: fails}, null, 1));
