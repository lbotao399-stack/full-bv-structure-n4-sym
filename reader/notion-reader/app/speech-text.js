/* Speech-only prose normalization. Never write this output into the manuscript. */
(function (root, factory) {
  'use strict';
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  if (root) root.ReaderSpeechText = api;
})(typeof window !== 'undefined' ? window : globalThis, function () {
  'use strict';
  function prose(text) {
    return text.replace(/\u00ad/g, '')
      .replace(/\b[A-Za-z]{2,}(?:[-\u2010\u2011][A-Za-z]{2,})+\b/g, word => word.replace(/[-\u2010\u2011]/g, ' '))
      .replace(/\bn[-\u2010\u2011]point\b/g, 'n point');
  }
  function normalize(value) {
    // Protect TeX, code, and URLs even when called on an authored speech override.
    const protectedParts = /(\$\$[\s\S]*?\$\$|\$[^$\n]*?\$|\\\[[\s\S]*?\\\]|\\\([\s\S]*?\\\)|`[^`]*`|https?:\/\/\S+)/g;
    return String(value || '').split(protectedParts).map((part, i) => i % 2 ? part : prose(part)).join('');
  }
  return {normalize, version: 'prose-hyphens-v1'};
});
