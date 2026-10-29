/* Speech-only prose normalization. Never write this output into the manuscript. */
(function (root, factory) {
  'use strict';
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  if (root) root.ReaderSpeechText = api;
})(typeof window !== 'undefined' ? window : globalThis, function () {
  'use strict';
  const digits = ['zero','one','two','three','four','five','six','seven','eight','nine'];
  function frequencies(value, language) {
    const chinese = language === 'zh-CN';
    const units = {Hz:['hertz','赫兹'],kHz:['kilohertz','千赫兹'],MHz:['megahertz','兆赫兹'],GHz:['gigahertz','吉赫兹'],THz:['terahertz','太赫兹']};
    const words = Object.fromEntries(Object.values(units).map(pair=>[pair[0],pair]));
    let text = value.replace(/(?<![A-Za-z])(?:THz|GHz|MHz|kHz|Hz|[Tt]erahertz|[Gg]igahertz|[Mm]egahertz|[Kk]ilohertz|[Hh]ertz)(?![A-Za-z])/g, token => (units[token] || words[token.toLowerCase()])[chinese?1:0]);
    // The English guide can treat a raw decimal point as a sentence stop.
    // Normalize quantities before synthesis, retaining all fractional digits.
    return text.replace(/(?<![\d.])([+−-]?)(\d+)\.(\d+)\s*(?=(?:terahertz|gigahertz|megahertz|kilohertz|hertz|太赫兹|吉赫兹|兆赫兹|千赫兹|赫兹)(?![A-Za-z]))/g, (_,sign,integer,fraction,offset) => {
      const mathIsland = /\b(?:equals|times|over|plus|minus)\b[^。！？;；]*$/i.test(text.slice(0,offset));
      const zh = chinese && !mathIsland;
      const digitWords = number => [...number].map(d=>zh?'零一二三四五六七八九'[Number(d)]:digits[Number(d)]).join(zh?'':' ');
      const whole = integer.length===1 || integer[0]==='0' ? digitWords(integer) : integer;
      const prefix = sign==='+' ? (zh?'正':'plus ') : sign ? (zh?'负':'minus ') : '';
      return prefix+whole+(zh?'点':' point ')+digitWords(fraction)+' ';
    });
  }
  function percentages(text, language) {
    const chinese = language ? language === 'zh-CN' : /[\u3400-\u9fff]/.test(text);
    return text.replace(/([+−-]?)(\d+(?:\.\d+)?|\.\d+)\s*(?:%|percent\b)/gi, (_, sign, number) => {
      // Work on characters, never floating point: retain every decimal digit.
      const [integer, fraction] = number.split('.');
      const digitWords = value => [...value].map(d => chinese ? '零一二三四五六七八九'[Number(d)] : digits[Number(d)]).join(chinese ? '' : ' ');
      // Leave a multi-digit integer as one token for ordinary cardinal reading;
      // explicitly spell the fractional digits so segmentation cannot regroup them.
      let spoken = integer && integer.length > 1 && integer[0] !== '0' ? integer : digitWords(integer || '0');
      if (fraction !== undefined) spoken += (chinese ? '点' : ' point ') + digitWords(fraction);
      const prefix = sign === '+' ? (chinese ? '正' : 'plus ') : sign ? (chinese ? '负' : 'minus ') : '';
      return chinese ? prefix + '百分之' + spoken : prefix + spoken + ' percent';
    });
  }
  function prose(text) {
    return text.replace(/\u00ad/g, '')
      .replace(/\b[A-Za-z]{2,}(?:[-\u2010\u2011][A-Za-z]{2,})+\b/g, word => word.replace(/[-\u2010\u2011]/g, ' '))
      .replace(/\bn[-\u2010\u2011]point\b/g, 'n point');
  }
  function normalize(value) {
    // Protect TeX, code, and URLs even when called on an authored speech override.
    const protectedParts = /(\$\$[\s\S]*?\$\$|\$[^$\n]*?\$|\\\[[\s\S]*?\\\]|\\\([\s\S]*?\\\)|`[^`]*`|https?:\/\/\S+)/g;
    const text = String(value || '');
    const language = /[\u3400-\u9fff]/.test(text) ? 'zh-CN' : 'en-US';
    return text.split(protectedParts).map((part, i) => i % 2 ? part : frequencies(percentages(prose(part), language), language)).join('');
  }
  return {normalize, percentages, frequencies, version: 'speech-frequency-v3'};
});
