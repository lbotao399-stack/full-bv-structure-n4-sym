/* Classify technical tokens before TTS. Never alter the displayed source. */
'use strict';
const units = {
  V:['volts','伏特'], mV:['millivolts','毫伏'], kV:['kilovolts','千伏'],
  A:['amperes','安培'], mA:['milliamperes','毫安'], 'μA':['microamperes','微安'],
  Hz:['hertz','赫兹'], kHz:['kilohertz','千赫兹'], MHz:['megahertz','兆赫兹'], GHz:['gigahertz','吉赫兹'],
  m:['meters','米'], cm:['centimeters','厘米'], mm:['millimeters','毫米'], 'μm':['micrometers','微米'],
  s:['seconds','秒'], ms:['milliseconds','毫秒'], 'm/s':['meters per second','米每秒'],
  'm/s^2':['meters per second squared','米每二次方秒'],
  Ω:['ohms','欧姆'], kΩ:['kiloohms','千欧姆'], W:['watts','瓦特'], mW:['milliwatts','毫瓦'],
  F:['farads','法拉'], 'μF':['microfarads','微法'], nF:['nanofarads','纳法'], pF:['picofarads','皮法']
};
const escape=s=>s.replace(/[.*+?^${}()|[\]\\]/g,'\\$&');
const unitPattern=Object.keys(units).sort((a,b)=>b.length-a.length).map(escape).join('|');
const digits=['zero','one','two','three','four','five','six','seven','eight','nine'];
function modelName(value){
  return [...value].map(c=>/\d/.test(c)?digits[Number(c)]:c==='/'?'slash':/[-‐‑]/.test(c)?'dash':c).join(' ');
}
function models(value){
  // Multiple capital prefix letters and digits distinguish instrument identifiers
  // from x-1, A/B and genuine fractions. Retain every digit, including leading zero.
  return value.replace(/\b[A-Z]{2,6}[-‐‑]\d[A-Za-z0-9]*(?:[/-][A-Za-z0-9]+)*\b/g,modelName);
}
function prose(value,language){
  const zh=language==='zh-CN',word=u=>units[u][zh?1:0];
  let text=models(value);
  text=text.replace(new RegExp('(\\d+(?:\\.\\d+)?)(?:\\s*([—–])\\s*(\\d+(?:\\.\\d+)?))?\\s+('+unitPattern+')(?![A-Za-z0-9])','g'),
    (_,a,dash,b,u)=>a+(dash?(zh?' 到 ':' to ')+b:'')+' '+word(u));
  // A table's stated physical quantity supplies the unit role; a bare V does not.
  text=text.replace(/(电压|[Vv]oltage)\s*\/\s*(mV|kV|V)\b/g,(_,label,u)=>label+(zh?'，单位 ':' in ')+word(u));
  text=text.replace(/(刻度|波长|长度|[Ww]avelength|[Ll]ength)\s*\/\s*(mm|cm|m)\b/g,(_,label,u)=>label+(zh?'，单位 ':' in ')+word(u));
  text=text.replace(/(频率|[Ff]requency)\s*\/\s*(GHz|MHz|kHz|Hz)\b/g,(_,label,u)=>label+(zh?'，单位 ':' in ')+word(u));
  text=text.replace(/\b(GHz|MHz|kHz|Hz)\b/g,u=>word(u));
  return text.replace(/(\d)\s*[×✖]️?\s*(?=[+−-]?\d)/g,'$1 times ');
}
function mathTex(tex){
  // Upright multi-letter SI units are explicit; isolated V/A/m remain variables
  // unless a numerical quantity immediately precedes their upright annotation.
  return tex.replace(/\\(?:mathrm|text)\s*\{([^{}]+)\}|\\mathrm\s+([VAmWsF])\b/g,(full,body,single,offset)=>{
    const token=(body||single).trim();
    if(!units[token])return full;
    const numeric=/\d(?:\s|\\[,;!]|\\quad|\\qquad)*$/.test(tex.slice(0,offset));
    if(token.length===1&&!numeric)return full;
    return '\\text{'+units[token][0]+'}';
  });
}
function mathSpoken(text){
  return require('../app/speech-text.js').percentages(text,'en-US')
    .replace(/(?<!capital |upper )\b([EI])\b/g,'capital $1');
}
module.exports={prose,models,mathTex,mathSpoken,version:'technical-roles-v2'};
