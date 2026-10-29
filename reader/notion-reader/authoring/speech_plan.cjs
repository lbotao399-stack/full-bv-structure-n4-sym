// Source-language prose, natural English mathematical readings, typed technical tokens.
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');
const readline = require('readline');
const katex = require('katex');
const sre = require('speech-rule-engine');
const technicalSpeech = require('./technical_speech.cjs');
const {tagOf,textAnnotations,normalizeMathText}=require('./formula_annotations.cjs');
const {referenceNumber, references} = require('./math_reading_style.cjs');
const naturalIndices = require('./natural_indices.cjs');
const MarkdownIt = require('markdown-it');
const md = new MarkdownIt({html: true});
const speechText = require('../app/speech-text.js');
const speechReferences = require('../app/speech-references.js');
const root = path.resolve(__dirname, '..');
let referenceStamp='', referenceLabels=[];
function equationReferences() {
  const file=[path.join(root,'data/document.json'),path.join(root,'document.seed.json')].find(file=>fs.existsSync(file));
  if(!file)return [];
  const stamp=file+':'+fs.statSync(file).mtimeMs;
  if(stamp!==referenceStamp){const value=JSON.parse(fs.readFileSync(file,'utf8'));referenceLabels=speechReferences.labels(value.document||value);referenceStamp=stamp;}
  return referenceLabels;
}
const readingsPath = path.join(__dirname, 'formula-readings.json');
let readings = {}, readingTime = -1;
function reviewed() {
  const mtime = fs.statSync(readingsPath).mtimeMs;
  if (mtime !== readingTime) { readings = JSON.parse(fs.readFileSync(readingsPath, 'utf8')); readingTime = mtime; }
  return readings;
}
function formula(tex, audit, language) {
  tex = tex.trim();
  const id = 'f' + crypto.createHash('sha256').update(tex).digest('hex').slice(0,12);
  // Validation also runs for reviewed text so broken edited formulas are never silently spoken.
  const tag = tagOf(tex);
  // A tag is metadata, whereas a minus inside the expression stays subtraction.
  const expression = tag ? tag.expression : tex;
  const prefix = tag ? (language === 'zh-CN' ? '式 ' : 'Equation ') + referenceNumber(tag.label, language) + (language === 'zh-CN' ? '，' : '. ') : '';
  katex.renderToString(tex,{displayMode:true,throwOnError:true,strict:'ignore',trust:false});
  let text = reviewed()[id];
  if (text && typeof text === 'object') text = text['en-US'] || text.en || '';
  if (text && !/[\u3400-\u9fff]/.test(text)) { audit.reviewed++; return prefix+technicalSpeech.mathSpoken(text); }
  const labels=JSON.parse(fs.readFileSync(path.join(__dirname,'formula-labels.json'),'utf8'));
  const markup = katex.renderToString(technicalSpeech.mathTex(textAnnotations(expression,labels)), {output: 'mathml', displayMode: true, throwOnError: true, strict:'ignore', trust:false});
  const rawMath = markup.match(/<math[\s\S]*?<\/math>/)?.[0];
  const math=rawMath&&normalizeMathText(rawMath);
  if (!math) throw new Error('This formula could not be converted to MathML.');
  text = naturalIndices(math, value => sre.toSpeech(value));
  text = text.replace(/\blamda\b/gi,'lambda').replace(/(?<=[\u3400-\u9fff])\s+(?=[\u3400-\u9fff])/g,'');
  text = text.replace(/\bSuperscript 2(?= Baseline|$)/g,'squared').replace(/\bSuperscript 3(?= Baseline|$)/g,'cubed');
  if (!text || /\\[A-Za-z]/.test(text)) throw new Error('Please enter a spoken reading for this formula.');
  audit.automatic++;
  text = technicalSpeech.mathSpoken(text);
  return prefix + technicalSpeech.models(text).replace(/\bStartAbsoluteValue\b/g, 'the absolute value of').replace(/\bEndAbsoluteValue\b/g, ', end absolute value,').replace(/\bnormal (?=[A-Za-z]\b)/g, '').replace(/\bupper (?=[A-Za-z]\b)/g, 'capital ').replace(/\bSuperscript\b/g, 'to the power').replace(/\bBaseline\b/g, ',').replace(/\bStartFraction\b/g, 'the fraction').replace(/\bOver\b/g, 'over').replace(/\bEndFraction\b/g, ', end fraction,').replace(/\bStartRoot\b/g,'the square root of').replace(/\bEndRoot\b/g, ', end root,');
}
function decode(s) { return s.replace(/&#(x[0-9a-f]+|\d+);/gi, (_, n) => String.fromCodePoint(n[0].toLowerCase()==='x'?parseInt(n.slice(1),16):Number(n))).replace(/&lt;/g,'<').replace(/&gt;/g,'>').replace(/&amp;/g,'&').replace(/&quot;/g,'"').replace(/&#39;|&apos;/g,"'").replace(/&nbsp;/g,' '); }
function sourceLanguage(block, fallback) {
  const source=block.type==='image'?(block.caption||block.alt||''):block.source||'';
  if(block.type==='equation')return /[\u3400-\u9fff]/.test(source)||block._readerNarration?.mode==='guided'&&/[\u3400-\u9fff]/.test(block.speechOverride||'')?'zh-CN':'en-US';
  const prose=block.type==='equation'?'':source.replace(/\$\$[\s\S]*?\$\$|\$`[\s\S]*?`\$|\$[^$]*?\$|\\\[[\s\S]*?\\\]|\\\([\s\S]*?\\\)/g,' ').replace(/<[^>]*>/g,' ');
  if (/[\u3400-\u9fff]/.test(prose)) return 'zh-CN';
  if (/[A-Za-z]{2,}/.test(prose)) return 'en-US';
  return fallback==='zh-CN'?'zh-CN':'en-US';
}
function plan(block, requestedLanguage) {
  const language=sourceLanguage(block, requestedLanguage);
  const chineseOverride=/[\u3400-\u9fff]/.test(block.speechOverride||'');
  const compatibleOverride=block.speechOverride?.trim() && (block.type==='equation'
    ? (!chineseOverride || block._readerNarration?.mode==='guided')
    : (language==='zh-CN' ? chineseOverride : !chineseOverride));
  const audit = {reviewed:0, automatic:0, override:false};
  let text;
  if (compatibleOverride) { text = speechText.normalize(block.speechOverride.trim()); audit.override = true; }
  else if (block.type === 'equation') text = formula(block.source, audit, language);
  else {
    const source = block.type === 'image' ? (block.caption || block.alt || (language==='zh-CN'?'图片。':'Figure.')) : block.source;
    // Protect math from Markdown/HTML: relation signs inside TeX are not HTML tags.
    const math = [];
    const protectedSource = source.replace(/\$`([\s\S]*?)`\$|\$\$([\s\S]*?)\$\$|\\\[([\s\S]*?)\\\]|\\\(([\s\S]*?)\\\)|\$([^$]+)\$/g, (_, notion, display, bracket, paren, inline) => {const n=math.length; math.push(formula(notion || display || bracket || paren || inline, audit, language)); return ` MATHSPEECHTOKEN${n}END `;});
    let html = md.render(protectedSource);
    // Tables are spoken row by row; keep visible header order.
    html = html.replace(/<tr\b[^>]*>/gi,language==='zh-CN'?' 一行： ':' Row: ').replace(/<\/(?:td|th)>/gi,'; ').replace(/<\/tr>/gi,'. ').replace(/<br\s*\/?\s*>/gi,'\n').replace(/<\/(?:p|div|h[1-6]|li|blockquote)>/gi,'\n');
    text = technicalSpeech.prose(speechText.normalize(decode(html.replace(/<[^>]*>/g,''))),language);
    text = text.replace(/MATHSPEECHTOKEN(\d+)END/g, (_, i) => math[Number(i)]);
  }
  text = speechReferences.normalize(text,{equationLabels:equationReferences(),language});
  text = speechText.normalize(text);
  text = text.replace(/,\s*,/g, ',').replace(/[ \t]+/g,' ').replace(/\s*\n\s*/g,'\n').replace(/;\s*\./g,'.').replace(/\s+([,.!?;:])/g,'$1').replace(/\.{2,}/g,'.').trim();
  if (!text) throw new Error('This block has no text to read.');
  // Timings will come from independently synthesized clauses, never word-count interpolation.
  const segments=[];
  for (const paragraph of text.split(/\n+/)) {
    for (const sentence of paragraph.split(language==='zh-CN'?/(?<=[。！？])/u:/(?<=[.!?])\s+(?=[A-Z0-9])/)) {
      if (!sentence.trim()) continue;
      if (sentence.length > (language==='zh-CN'?110:320)) {
        // Decimal points inside English formula islands are not clause boundaries.
        const clauses=sentence.split(language==='zh-CN'?/(?<=[，；：,!?;])|(?<=[.:])(?=\s|$)/u:/(?<=[;:])\s+/); let part='';
        for(const clause of clauses) { if(part && part.length+clause.length>(language==='zh-CN'?110:320)){segments.push(part.trim());part='';} part+=(part?' ':'')+clause; }
        if(part.trim())segments.push(part.trim());
      } else segments.push(sentence.trim());
    }
  }
  return {speech:text, segments, language, formulaReadings:audit};
}
(async () => {
  await sre.setupEngine({locale:'en', domain:'mathspeak', style:'default'}); await sre.engineReady();
  const lines=readline.createInterface({input:process.stdin, terminal:false});
  for await (const line of lines) {
    try { const request=JSON.parse(line); console.log(JSON.stringify({ok:true, ...plan(request.block,request.language)})); }
    catch(e) { console.log(JSON.stringify({ok:false,error:String(e.message || e)})); }
  }
})();
