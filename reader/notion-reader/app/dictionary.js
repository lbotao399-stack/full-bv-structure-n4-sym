/* Automatic local definitions for selected prose; never edits or starts audio. */
(function () {
  'use strict';
  let panel, content, title, pending, previousFocus, anchorRect;
  let selectionKey = '', dismissed = false, debounce;
  const cache = new Map();
  function close() {
    clearTimeout(debounce); debounce = null; dismissed = true;
    pending?.abort(); pending = null;
    if (!panel) return;
    const restore = panel.contains(document.activeElement);
    panel.hidden = true;
    if (restore && previousFocus?.isConnected) previousFocus.focus({preventScroll: true});
  }
  function update(selection, rect, key = selection) {
    const term = String(selection || '').trim();
    const plain = term.replace(/[’]/g, "'").replace(/[‐‑]/g, '-').replace(/^[.,;:!?“”"()\[\]{}]+|[.,;:!?“”"()\[\]{}]+$/g, '');
    if (!plain || term.length > 80 || term.split(/\s+/).length > 8 || !/^[A-Za-z][A-Za-z '\-]*$/.test(plain)) {
      close(); selectionKey = ''; return;
    }
    if (key === selectionKey) {if (!dismissed && panel && !panel.hidden) position(rect); return;}
    close(); selectionKey = key; dismissed = false; anchorRect = rect;
    // Selection changes during a drag or Shift+Arrow settle before a local request.
    debounce = setTimeout(() => {debounce = null; if (!dismissed) void lookup(term, anchorRect);}, 180);
  }
  function ensurePanel() {
    if (panel) return;
    panel = document.createElement('section'); panel.id = 'reader-dictionary'; panel.hidden = true;
    panel.setAttribute('role', 'dialog'); panel.setAttribute('aria-label', 'Chinese dictionary');
    const header = document.createElement('header');
    title = document.createElement('strong');
    const button = document.createElement('button'); button.type = 'button'; button.textContent = '×';
    button.setAttribute('aria-label', 'Close dictionary'); button.addEventListener('click', close);
    header.append(title, button);
    content = document.createElement('div'); content.className = 'dictionary-content';
    content.setAttribute('role', 'status'); content.setAttribute('aria-live', 'polite');
    panel.append(header, content); document.body.append(panel);
    // Preserve the source selection and prevent reader playback shortcuts here.
    panel.addEventListener('pointerdown', event => {event.stopPropagation(); if (event.target.closest('button')) event.preventDefault();});
    panel.addEventListener('click', event => event.stopPropagation());
    panel.addEventListener('keydown', event => {event.stopPropagation(); if(event.key === 'Escape'){event.preventDefault(); close();}});
    document.addEventListener('keydown', event => {
      if (!panel.hidden && event.key === 'Escape') {event.preventDefault(); event.stopImmediatePropagation(); close();}
    }, true);
    document.addEventListener('pointerdown', event => {
      if (!panel.hidden && !panel.contains(event.target) && !event.target.closest('#selection-formatting')) close();
    });
    window.addEventListener('resize', () => {if (!panel.hidden) position();});
  }
  function position(rect) {
    if (rect) anchorRect = rect;
    rect = anchorRect;
    const width = Math.min(370, Math.max(180, window.innerWidth - 24));
    panel.style.width = width + 'px';
    const bottom = Math.min(window.innerHeight - 12, (document.querySelector('.player')?.getBoundingClientRect().top || window.innerHeight) - 8);
    panel.style.maxHeight = Math.max(120, bottom - 12) + 'px';
    panel.style.left = Math.max(12, Math.min(rect?.left || 12, window.innerWidth - width - 12)) + 'px';
    let top = rect ? rect.bottom + 8 : 80;
    if (top + panel.offsetHeight > bottom) top = (rect?.top || bottom) - panel.offsetHeight - 8;
    panel.style.top = Math.max(12, Math.min(top, bottom - panel.offsetHeight)) + 'px';
  }
  function text(tag, value, className) {
    const el = document.createElement(tag); el.textContent = value; if(className)el.className = className; return el;
  }
  async function lookup(selection, rect) {
    ensurePanel(); pending?.abort();
    const controller = new AbortController(); pending = controller;
    previousFocus = document.activeElement;
    const term = String(selection || '').trim();
    title.textContent = term || 'Chinese dictionary'; content.textContent = 'Looking up locally…';
    panel.hidden = false; position(rect);
    try {
      if (!term || term.includes('\ufffc') || term.length > 80 || term.split(/\s+/).length > 8) throw Error('Select an English word or short phrase without formulas.');
      let result = cache.get(term);
      if (!result) {
        const response = await fetch('/api/dictionary', {method:'POST', headers:{'Content-Type':'application/json'}, body:JSON.stringify({term}), signal:controller.signal});
        result = await response.json();
        if (!response.ok) throw Error(result.error || 'Local dictionary unavailable.');
        if (controller.signal.aborted) return;
        cache.set(term, result); if (cache.size > 64) cache.delete(cache.keys().next().value);
      }
      if (controller.signal.aborted) return;
      content.replaceChildren();
      if (!result.entries?.length) content.append(text('p', result.message || 'No local entry for this selection.'));
      for (const entry of result.entries || []) {
        const item = document.createElement('div'); item.className = 'dictionary-entry';
        item.append(text('b', entry.word));
        if (entry.phonetic) item.append(text('span', ' /' + entry.phonetic + '/', 'dictionary-phonetic'));
        if (entry.match === 'inflection') item.append(text('small', 'Word form of ' + entry.word));
        const translation = text('p', entry.translation || 'No Chinese definition in this entry.', 'dictionary-translation'); translation.lang = 'zh-CN'; item.append(translation);
        if (entry.definition) {const details = document.createElement('details');details.append(text('summary','English definition'),text('p',entry.definition));item.append(details);}
        item.append(text('small', entry.source)); content.append(item);
      }
      position(rect);
    } catch (error) {
      if (!controller.signal.aborted) {content.textContent = error.message; position(rect);}
    }
  }
  window.ReaderDictionary = {lookup, close, update};
})();
