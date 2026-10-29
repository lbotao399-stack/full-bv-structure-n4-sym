/* Explicit language-specific voice routes; no silent provider substitution. */
(function (root, factory) {
  'use strict';
  if (typeof module === 'object' && module.exports) module.exports = factory;
  if (root) root.ReaderVoiceTransport = factory();
})(typeof window !== 'undefined' ? window : null, function (environment) {
  'use strict';
  environment = environment || {};
  function speechText() {
    const helper = environment.speechText || (typeof module === 'object' && module.exports ? require('./speech-text.js') : window.ReaderSpeechText);
    if (!helper || typeof helper.normalize !== 'function' || typeof helper.version !== 'string') {
      throw Error('The speech helper did not load. Reload this reader; if the error persists, repair its speech-text.js script installation.');
    }
    return helper;
  }
  let configuration;
  const profileIds = {};
  const profiles = new Map();
  let voiceChoice = 'auto';
  try {
    const storage = environment.storage || localStorage;
    const saved = storage.getItem('reader-voice-choice-v2');
    const legacy = JSON.parse(storage.getItem('reader-voice-choices-v1') || '{}');
    if (saved) voiceChoice = saved;
    else if (legacy['zh-CN'] && legacy['zh-CN'] === legacy['en-US'] && legacy['zh-CN'] !== 'default') voiceChoice = legacy['zh-CN'];
  } catch (_) {}
  const selected = language => voiceChoice === 'auto' ? (language === 'zh-CN' ? 'saki' : 'satsuki') : voiceChoice;
  const request = (...args) => (environment.fetch || fetch)(...args);
  function localService(value) {
    const url = new URL(value);
    if (url.protocol !== 'http:' || !['127.0.0.1', 'localhost', '[::1]'].includes(url.hostname) || url.username || url.password || url.search || url.hash || !/^(?:\/|\/voices\/(?:sakiko|satsuki|anon|uika|mai|renako)\/(?:en-US|zh-CN))?$/.test(url.pathname)) throw Error('The voice service must use a local loopback address.');
    return url.origin + (url.pathname === '/' ? '' : url.pathname);
  }
  async function config() {
    if (!configuration) configuration = (async () => {
      const response = await request(environment.configURL || new URL('speech-service.json', document.baseURI).href, {cache: 'no-store'});
      if (response.status === 404) return null;
      if (!response.ok) throw Error('Unable to load the narration voice settings.');
      const value = await response.json();
      if (value.version !== 1) throw Error('Unsupported narration voice settings.');
      const result = {...value};
      if (value.enabled !== false && value.mandarinBackend !== 'kokoro') result.url = localService(value.url);
      if (value.englishService?.enabled) result.englishService = {...value.englishService, url: localService(value.englishService.url)};
      return result;
    })();
    try { return await configuration; }
    catch (error) { configuration = null; throw error; }
  }
  function checkCancellation(signal) {
    if (signal?.aborted) throw new DOMException('Narration preparation cancelled.', 'AbortError');
  }
  function waitFor(value, signal) {
    if (!signal) return value;
    checkCancellation(signal);
    return new Promise((resolve, reject) => {
      const cancel = () => { signal.removeEventListener('abort', cancel); reject(new DOMException('Narration preparation cancelled.', 'AbortError')); };
      signal.addEventListener('abort', cancel, {once: true});
      Promise.resolve(value).then(result => { signal.removeEventListener('abort', cancel); resolve(result); }, error => { signal.removeEventListener('abort', cancel); reject(error); });
    });
  }
  async function post(url, body, signal) {
    checkCancellation(signal);
    const response = await request(url, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(body), signal});
    const record = await response.json();
    checkCancellation(signal);
    if (!response.ok || record.ok === false) throw Error(record.error || 'Narration preparation failed.');
    return record;
  }
  async function start(block, {signal, background = false} = {}) {
    checkCancellation(signal);
    // The existing source-language/math planner remains authoritative, including
    // newly typed Mandarin in an initially English manuscript.
    const plan = await post('/api/plan', {block}, signal);
    const spokenSegments = plan.segments?.map(text => speechText().normalize(text));
    if (!['en-US', 'zh-CN'].includes(plan.language)) throw Error('The source-language planner did not identify a supported reading language.');
    const configuration = await waitFor(config(), signal);
    checkCancellation(signal);
    const mandarin = plan.language === 'zh-CN';
    const choice = selected(plan.language);
    let settings = mandarin ? configuration : configuration?.englishService;
    if (configuration?.choices?.[choice]?.available === false) throw Error('The selected voice model is not installed.');
    if (configuration?.choices && !configuration.choices[choice]?.[plan.language]) throw Error('The selected voice is unavailable for this language.');
    if (choice && configuration?.choices?.[choice]?.[plan.language]) settings = configuration.choices[choice][plan.language];
    const useService = settings?.backend !== 'kokoro' && (mandarin ? configuration?.mandarinBackend !== 'kokoro' || Boolean(choice) : settings?.enabled);
    if (useService) {
        if (!settings || settings.enabled === false) throw Error('Mandarin narration is not configured. Configure the Mandarin voice service, or explicitly choose the legacy Kokoro Mandarin backend.');
        if (!Array.isArray(plan.segments) || !plan.segments.length) throw Error('This passage has no speech segments.');
        let record;
        let profileId;
        try {
          settings = {...settings, url: localService(settings.url)};
          let profile = profiles.get(settings.url);
          if (!profile || profile.expires < Date.now()) {
            const profileResponse = await request(settings.url + '/api/profile', {cache: 'no-store', signal});
            profile = await profileResponse.json();
            if (!profileResponse.ok || typeof profile.profileId !== 'string' || !profile.profileId.trim()) throw Error(profile.error || 'Unable to identify the configured voice profile.');
            profile = {...profile, expires: Date.now() + 60000};
            profiles.set(settings.url, profile);
          }
          checkCancellation(signal);
          if (profile.language && profile.language !== plan.language) throw Error('The voice service language does not match this passage.');
          const selected = settings.voice === 'saki' ? 'sakiko' : settings.voice;
          if (selected && selected !== 'default' && profile.voice !== selected) throw Error('The voice service does not match the selected character.');
          profileId = profile.profileId;
          profileIds[plan.language] = profileId;
          record = await post(settings.url + '/api/prepare', {segments: spokenSegments, language: plan.language, background}, signal);
        } catch (error) {
          if (error.name === 'AbortError') throw error;
          throw Error((mandarin ? 'Mandarin' : 'English') + ' voice unavailable. Keep the selected voice service running. ' + error.message);
        }
        return {record, profileId, jobURL: record.jobId ? settings.url + '/api/jobs/' + encodeURIComponent(record.jobId) : null};
    }
    const normalizedBlock = spokenSegments?.some((text, i) => text !== plan.segments[i])
      ? {...block, speechOverride: spokenSegments.join('\n')} : block;
    const record = await post('/api/prepare', {block: normalizedBlock, background}, signal);
    return {record, jobURL: record.jobId ? '/api/jobs/' + encodeURIComponent(record.jobId) : null};
  }
  function selectVoice(voice) {
    voiceChoice = voice;
    try { (environment.storage || localStorage).setItem('reader-voice-choice-v2', voice); } catch (_) {}
    if (typeof window !== 'undefined') window.dispatchEvent(new Event('reader-voice-change'));
  }
  return {start, config, select: (_language, voice) => selectVoice(voice), selectVoice, selected, choice: () => voiceChoice, cacheTag: () => 'single-voice-v9:' + speechText().version + ':' + voiceChoice + ':' + JSON.stringify(profileIds)};
});
