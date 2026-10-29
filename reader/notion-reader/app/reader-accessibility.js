/* Reach the existing narration controls without changing audio or manuscript. */
(() => {
  'use strict';
  function init() {
    if (window.ReaderAccessibleControls) return;
    const player = document.querySelector('.player, .audio-dock');
    if (!player) return;
    const embedded = window.parent !== window;
    let parentOrigin = '';
    try { parentOrigin = new URL(document.referrer).origin; } catch (_) {}
    let expanded = false, previousFocus = null;
    player.classList.add('reader-accessible-player');
    if (!player.id) player.id = 'reader-narration-controls';
    const heading = document.createElement('div');
    heading.className = 'reader-player-heading';
    const title = document.createElement('strong');
    title.textContent = 'Narration';
    const close = document.createElement('button');
    close.type = 'button'; close.textContent = 'Close';
    close.setAttribute('aria-label', 'Close narration controls');
    heading.append(title, close); player.prepend(heading);
    const backdrop = document.createElement('div');
    backdrop.className = 'reader-player-backdrop'; backdrop.hidden = true;
    document.body.append(backdrop);
    let trigger = null;
    if (!embedded && !document.body.classList.contains('controls-only')) {
      trigger = document.createElement('button');
      trigger.type = 'button'; trigger.id = 'reader-show-narration';
      trigger.className = 'quiet'; trigger.textContent = 'Narration';
      trigger.setAttribute('aria-controls', player.id);
      trigger.setAttribute('aria-expanded', 'false');
      (document.querySelector('header.topbar, .topbar') || document.body).append(trigger);
      trigger.addEventListener('click', () => setExpanded(!expanded));
    }
    const originalRole = player.getAttribute('role');
    function announce() {
      trigger?.setAttribute('aria-expanded', String(expanded));
      if (embedded && parentOrigin) parent.postMessage({type: 'reader:narration-state', open: expanded}, parentOrigin);
    }
    function setExpanded(open) {
      if (expanded === Boolean(open)) return;
      expanded = Boolean(open);
      if (expanded) previousFocus = document.activeElement;
      document.body.classList.toggle('reader-player-expanded', expanded);
      backdrop.hidden = !expanded;
      if (expanded) {
        player.setAttribute('role', 'dialog'); player.setAttribute('aria-modal', 'true');
        close.focus({preventScroll: true});
      } else {
        player.removeAttribute('aria-modal');
        if (originalRole) player.setAttribute('role', originalRole); else player.removeAttribute('role');
        if (!embedded && previousFocus?.isConnected) previousFocus.focus({preventScroll: true});
      }
      measure(); announce();
    }
    close.addEventListener('click', () => setExpanded(false));
    backdrop.addEventListener('click', () => setExpanded(false));
    window.addEventListener('message', event => {
      if (!embedded || event.source !== parent || event.origin !== parentOrigin || event.data?.type !== 'reader:narration-controls') return;
      setExpanded(event.data.open !== false);
    });
    document.addEventListener('keydown', event => {
      if (!expanded) return;
      if (event.key === 'Escape') { event.preventDefault(); setExpanded(false); }
      if (event.key === 'Tab') {
        const controls = [...player.querySelectorAll('button,select,input,[tabindex="0"]')].filter(node => !node.disabled && node.getClientRects().length);
        const first = controls[0], last = controls.at(-1);
        if (event.shiftKey && document.activeElement === first) { event.preventDefault(); last?.focus(); }
        else if (!event.shiftKey && document.activeElement === last) { event.preventDefault(); first?.focus(); }
      }
    });
    function measure() {
      const viewport = window.visualViewport;
      const height = viewport?.height || innerHeight;
      const top = viewport?.offsetTop || 0;
      const inset = Math.max(0, innerHeight - height - top);
      document.documentElement.style.setProperty('--reader-visible-height', height + 'px');
      document.documentElement.style.setProperty('--reader-visible-top', top + 'px');
      document.documentElement.style.setProperty('--reader-bottom-inset', inset + 'px');
      if (!expanded) {
        const size = Math.ceil(player.getBoundingClientRect().height + inset);
        document.documentElement.style.setProperty('--reader-player-space', size + 'px');
      }
    }
    const observer = new ResizeObserver(measure); observer.observe(player);
    window.visualViewport?.addEventListener('resize', measure);
    window.visualViewport?.addEventListener('scroll', measure);
    window.addEventListener('resize', measure);
    window.ReaderAccessibleControls = {open: () => setExpanded(true), close: () => setExpanded(false)};
    measure();
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
