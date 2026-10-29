(() => {
  'use strict';
  async function init() {
    const transport = window.ReaderVoiceTransport;
    const config = await transport?.config().catch(() => null);
    if (!config?.choices || document.querySelector('.reader-voice-picker')) return;
    const row = document.querySelector('.player .player-row');
    if (!row) return;
    const panel = document.createElement('label'); panel.className = 'reader-voice-picker full-only';
    const select = document.createElement('select');
    select.id = 'reader-voice'; select.setAttribute('aria-label', 'Voice');
    select.append(new Option('Default', 'auto'));
    for (const [key, value] of Object.entries(config.choices)) {
      if (!value['zh-CN'] || !value['en-US']) continue;
      const option = new Option(key === 'default' ? 'Legacy voice' : value.label.replace(' · Sakiko Togawa', ''), key);
      option.disabled = value.available === false;
      select.append(option);
    }
    select.value = transport.choice();
    if (!select.value) {
      const unavailable = new Option('Saved voice unavailable', transport.choice());
      unavailable.disabled = true; select.append(unavailable); select.value = transport.choice();
    }
    const sync = () => {
      select.value = transport.choice();
      select.title = select.value === 'auto' ? 'English: Satsuki · Mandarin: Saki' : select.selectedOptions[0]?.textContent || 'Voice';
    };
    select.onchange = () => { transport.selectVoice(select.value); sync(); };
    window.addEventListener('reader-voice-change', sync);
    panel.append('Voice', select);
    row.insertBefore(panel, row.querySelector('#speed')?.closest('label') || row.querySelector('#speed'));
    sync();
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init, {once:true}); else void init();
})();
