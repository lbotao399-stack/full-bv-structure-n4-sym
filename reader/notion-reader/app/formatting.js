/* Presentation-only annotations. Source text, TeX and speech stay unchanged. */
(function () {
  'use strict';
  const TEXT = {
    default: null, gray: '#69716d', red: '#b13737', amber: '#9a5c0c',
    green: '#21714d', blue: '#245faa', purple: '#7949a4'
  };
  const HIGHLIGHT = {
    none: null, yellow: '#fff0a7', green: '#d7edce', blue: '#d6e9ff',
    pink: '#f9dce5', purple: '#e8dff7'
  };
  const ATOM = '\ufffc';
  let api = {}, toolbar, targets = [], pinned = false, captureQueued = false, initialized = false, pointerSelecting = false;
  const $ = (selector, root = document) => root.querySelector(selector);
  const blockElement = id => [...document.querySelectorAll('.source-block')].find(el => el.dataset.blockId === id);
  const getBlock = id => api.getDocument?.()?.blocks?.find(block => block.id === id);
  const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);

  function cleanStyle(style) {
    const out = {};
    if (style && TEXT[style.color]) out.color = style.color;
    if (style && (HIGHLIGHT[style.highlight] || style.highlight === 'none')) out.highlight = style.highlight;
    for (const key of ['bold', 'italic', 'underline']) if (style?.[key] === true) out[key] = true;
    return out;
  }
  function projection(root) {
    const units = [], math = [];
    let text = '';
    function walk(node) {
      if (node.nodeType === Node.TEXT_NODE) {
        if (!node.data) return;
        units.push({ node, start: text.length, end: text.length + node.data.length, type: 'text' });
        text += node.data; return;
      }
      if (node.nodeType !== Node.ELEMENT_NODE) return;
      if (node.matches('.katex-mathml,script,style,[data-reader-format-ignore]')) return;
      if (node.matches('.math-inline,.math-display')) {
        const tex = node.querySelector('annotation[encoding="application/x-tex"]')?.textContent || node.textContent;
        units.push({ node, start: text.length, end: text.length + 1, type: 'math', tex });
        math.push({ offset: text.length, tex }); text += ATOM; return;
      }
      for (const child of node.childNodes) walk(child);
    }
    walk(root); return { text, units, math };
  }
  function mathWithin(p, start, end) {
    return p.math.filter(item => item.offset >= start && item.offset < end).map(item => ({ offset: item.offset - start, tex: item.tex }));
  }
  function anchor(p, start, end, style) {
    return {
      start, end, quote: p.text.slice(start, end),
      before: p.text.slice(Math.max(0, start - 24), start), after: p.text.slice(end, end + 24),
      math: mathWithin(p, start, end), style: cleanStyle(style)
    };
  }
  function resolveAnchor(item, p, originalText) {
    if (!item || typeof item.quote !== 'string' || !item.quote.length || !Number.isInteger(item.start) || !Number.isInteger(item.end)) return null;
    const length = item.quote.length;
    const matches = start => start >= 0 && p.text.slice(start, start + length) === item.quote &&
      same(mathWithin(p, start, start + length), item.math || []);
    if (p.text === originalText && matches(item.start)) return { start: item.start, end: item.start + length, style: cleanStyle(item.style) };
    const candidates = [];
    for (let at = p.text.indexOf(item.quote), count = 0; at >= 0 && count < 500; at = p.text.indexOf(item.quote, at + 1), count++) {
      if (!matches(at)) continue;
      const before = typeof item.before === 'string' ? item.before : '', after = typeof item.after === 'string' ? item.after : '';
      let score = 0;
      if (before && p.text.slice(Math.max(0, at - before.length), at) === before) score += before.length;
      if (after && p.text.slice(at + length, at + length + after.length) === after) score += after.length;
      if (!before && !after && p.text === item.quote) score = 1;
      if (score) candidates.push({ start: at, end: at + length, score, style: cleanStyle(item.style) });
    }
    candidates.sort((a, b) => b.score - a.score);
    if (!candidates.length || candidates[1]?.score === candidates[0].score) return null;
    return candidates[0];
  }
  function resolveRanges(formatting, p) {
    if (!formatting || formatting.version !== 1 || !Array.isArray(formatting.ranges)) return [];
    return formatting.ranges.slice(0, 2000).map(item => resolveAnchor(item, p, formatting.text)).filter(item => item && Object.keys(item.style).length);
  }
  function atStyle(ranges, offset) {
    const style = {};
    for (const range of ranges) if (range.start <= offset && offset < range.end) Object.assign(style, range.style);
    return cleanStyle(style);
  }
  function sourceHighlight(element) {
    if (element.matches('.reader-format-mark')) return null;
    if (element.dataset.readerSourceHighlight) return element.dataset.readerSourceHighlight;
    const background = element.style.backgroundColor;
    if (background) {
      const swatch = document.createElement('span');
      for (const [name, value] of [...Object.entries(HIGHLIGHT), ...['#fff1b7', '#ffff00', '#fbf3db', '#fdecc8'].map(value => ['yellow', value])]) {
        if (!value) continue;
        swatch.style.backgroundColor = value;
        if (background === swatch.style.backgroundColor) return name;
      }
      return null;
    }
    const notion = (element.getAttribute('color') || '').replace(/_bg$/, '');
    if (HIGHLIGHT[notion] && element.getAttribute('color').endsWith('_bg')) return notion;
    return element.matches('mark,.source-highlight') ? 'yellow' : null;
  }
  function sourceHighlights(p, flatten = false) {
    const ranges = [], elements = new Map();
    for (const unit of p.units) {
      let element = unit.type === 'math' ? unit.node : unit.node.parentElement, highlight;
      while (element && !element.matches('.block-content')) {
        const value = sourceHighlight(element);
        if (value) { highlight ||= value; elements.set(element, value); }
        element = element.parentElement;
      }
      if (!highlight) continue;
      const previous = ranges.at(-1);
      if (previous?.end === unit.start && previous.style.highlight === highlight) previous.end = unit.end;
      else ranges.push({start: unit.start, end: unit.end, style: {highlight}});
    }
    // Paint source and learner highlights through the same range layer. Keeping
    // the original carrier in the DOM lets selection toggles see the source
    // baseline; source strings themselves are never rewritten by formatting.
    if (flatten) for (const [element, highlight] of elements) {
      element.dataset.readerSourceHighlight = highlight;
      element.classList.remove('source-highlight');
      if ((element.getAttribute('color') || '').endsWith('_bg')) element.removeAttribute('color');
      element.style.removeProperty('background-color');
      if (element.tagName === 'MARK') element.style.backgroundColor = 'transparent';
    }
    return ranges;
  }
  function updateRanges(ranges, start, end, patch, length) {
    const cuts = [...new Set([0, length, start, end, ...ranges.flatMap(item => [item.start, item.end])])].filter(x => x >= 0 && x <= length).sort((a, b) => a - b);
    const out = [];
    for (let i = 0; i < cuts.length - 1; i++) {
      const left = cuts[i], right = cuts[i + 1]; if (right <= left) continue;
      let style = atStyle(ranges, left);
      if (left >= start && right <= end) {
        if (patch === null) style = {highlight: 'none'};
        else for (const [key, value] of Object.entries(patch)) {
          if (key === 'highlight' && value === null) style.highlight = 'none';
          else if (value === null || value === false) delete style[key];
          else style[key] = value;
        }
      }
      style = cleanStyle(style); if (!Object.keys(style).length) continue;
      const last = out[out.length - 1];
      if (last && last.end === left && same(last.style, style)) last.end = right;
      else out.push({ start: left, end: right, style });
    }
    return out;
  }
  function applyStyle(element, style) {
    element.classList.add('reader-format-mark');
    if (style.color) { element.style.color = TEXT[style.color]; element.dataset.readerColor = style.color; }
    if (style.highlight) { element.style.backgroundColor = HIGHLIGHT[style.highlight] || 'transparent'; element.dataset.readerHighlight = style.highlight; }
    if (style.bold) element.style.fontWeight = '700';
    if (style.italic) element.style.fontStyle = 'italic';
    if (style.underline) { element.style.textDecorationLine = 'underline'; element.style.textUnderlineOffset = '.15em'; }
  }
  function decorate(block, content) {
    const p = projection(content), ranges = [...sourceHighlights(p, true), ...resolveRanges(block.formatting, p)];
    if (!ranges.length) return;
    for (const unit of p.units) {
      if (unit.type === 'math') {
        const style = atStyle(ranges, unit.start);
        if (Object.keys(style).length) applyStyle(unit.node, style);
        continue;
      }
      const cuts = [...new Set([unit.start, unit.end, ...ranges.flatMap(range => [Math.max(unit.start, range.start), Math.min(unit.end, range.end)])])]
        .filter(x => x >= unit.start && x <= unit.end).sort((a, b) => a - b);
      if (!cuts.some(x => Object.keys(atStyle(ranges, x)).length)) continue;
      const fragment = document.createDocumentFragment();
      for (let i = 0; i < cuts.length - 1; i++) {
        const left = cuts[i], right = cuts[i + 1]; if (right <= left) continue;
        const text = unit.node.data.slice(left - unit.start, right - unit.start), style = atStyle(ranges, left);
        if (!Object.keys(style).length) fragment.append(document.createTextNode(text));
        else { const span = document.createElement('span'); span.dataset.readerFormat = ''; applyStyle(span, style); span.textContent = text; fragment.append(span); }
      }
      unit.node.replaceWith(fragment);
    }
  }
  function contentRange(unit) {
    const r = document.createRange();
    if (unit.type === 'math') r.selectNode(unit.node); else r.selectNodeContents(unit.node);
    return r;
  }
  function selectionPart(range, p) {
    let start = Infinity, end = -Infinity;
    for (const unit of p.units) {
      try {
        const r = contentRange(unit);
        // Strict overlap: touching a boundary does not select the adjacent unit.
        if (range.compareBoundaryPoints(Range.END_TO_START, r) >= 0 || range.compareBoundaryPoints(Range.START_TO_END, r) <= 0) continue;
        let left = 0, right = unit.end - unit.start;
        if (unit.type === 'text') {
          if (range.startContainer === unit.node) left = range.startOffset;
          if (range.endContainer === unit.node) right = range.endOffset;
        }
        if (right > left) { start = Math.min(start, unit.start + left); end = Math.max(end, unit.start + right); }
      } catch (_) { /* Detached selections are ignored. */ }
    }
    return end > start ? { start, end } : null;
  }
  function currentTargets() {
    const selection = window.getSelection();
    if (!selection || selection.rangeCount !== 1 || selection.isCollapsed) return [];
    const range = selection.getRangeAt(0);
    const found = [];
    for (const root of document.querySelectorAll('.source-block .block-content')) {
      const id = root.closest('.source-block')?.dataset.blockId; if (!id || !getBlock(id)) continue;
      try { if (!range.intersectsNode(root)) continue; } catch (_) { continue; }
      const p = projection(root), part = selectionPart(range, p);
      if (part) found.push({ id, text: p.text, ...anchor(p, part.start, part.end, {}) });
    }
    return found;
  }
  function rectForTargets() {
    const selection = window.getSelection();
    if (selection && !selection.isCollapsed && selection.rangeCount) {
      const rect = selection.getRangeAt(0).getBoundingClientRect(); if (rect.width || rect.height) return rect;
    }
    const root = targets[0] && blockElement(targets[0].id);
    return root?.getBoundingClientRect();
  }
  function placeToolbar(rect) {
    if (!toolbar || !rect) return;
    const pad = 12, topbar = document.querySelector('.topbar')?.getBoundingClientRect().bottom || 0;
    toolbar.hidden = false;
    const box = toolbar.getBoundingClientRect();
    const left = Math.max(pad, Math.min(window.innerWidth - box.width - pad, rect.left + rect.width / 2 - box.width / 2));
    let top = rect.top - box.height - 9;
    if (top < topbar + 8) top = Math.max(topbar + 8, rect.bottom + 9);
    const playerTop = document.querySelector('.player')?.getBoundingClientRect().top || window.innerHeight;
    top = Math.min(top, playerTop - box.height - 10); top = Math.max(pad, top);
    toolbar.style.left = left + 'px'; toolbar.style.top = top + 'px';
  }
  function capture() {
    captureQueued = false;
    if (!initialized || pointerSelecting || document.body.classList.contains('floating')) return;
    if (toolbar?.contains(document.activeElement) || document.activeElement?.matches('input,textarea,select')) return;
    const selectionNode = window.getSelection()?.anchorNode;
    if ((selectionNode?.nodeType === 1 ? selectionNode : selectionNode?.parentElement)?.closest('#reader-dictionary')) return;
    const found = currentTargets();
    if (!found.length) { if (!pinned) { hideToolbar(); window.ReaderDictionary?.update(''); } return; }
    targets = found; pinned = false; refreshToolbar(); placeToolbar(rectForTargets());
    window.ReaderDictionary?.update(targets.map(item => item.quote || '').join(' ').trim(), toolbar.getBoundingClientRect(),
      JSON.stringify(targets.map(item => [item.id, item.start, item.end, item.quote])));
  }
  function queueCapture() { if (!captureQueued) { captureQueued = true; requestAnimationFrame(capture); } }
  function hideToolbar() { if (toolbar) toolbar.hidden = true; targets = []; pinned = false; window.ReaderDictionary?.close(); }
  function clearSelection() { hideToolbar(); }
  function resolvedTargets() {
    return targets.map(target => {
      const block = getBlock(target.id), root = blockElement(target.id)?.querySelector('.block-content'); if (!block || !root) return null;
      const p = projection(root), resolved = resolveAnchor(target, p, target.text || p.text);
      return resolved ? { block, p, start: resolved.start, end: resolved.end } : null;
    }).filter(Boolean);
  }
  function rangeHasStyle(ranges, start, end, property, value = true) {
    if (end <= start) return false;
    const cuts = [...new Set([start, end, ...ranges.flatMap(x => [x.start, x.end])])]
      .filter(x => x >= start && x <= end).sort((a, b) => a - b);
    return cuts.slice(0, -1).every(x => atStyle(ranges, x)[property] === value);
  }
  function allStyled(property, value = true) {
    const resolved = resolvedTargets(); if (!resolved.length) return false;
    return resolved.every(({ block, p, start, end }) => {
      const ranges = [...sourceHighlights(p), ...resolveRanges(block.formatting, p)];
      const spans = property === 'bold' || property === 'italic' ? p.units.filter(unit => unit.type === 'text').map(unit => ({ start: Math.max(start, unit.start), end: Math.min(end, unit.end) })).filter(span => span.end > span.start) : [{ start, end }];
      return spans.every(span => rangeHasStyle(ranges, span.start, span.end, property, value));
    });
  }
  function hasSelectedProse() { return targets.some(target => target.quote.replaceAll(ATOM, '').trim().length > 0); }
  function refreshToolbar() {
    if (!toolbar) return;
    const count = targets.reduce((n, t) => n + (t.quote.match(/\S+/g)?.length || 0), 0);
    $('#format-selection-label', toolbar).textContent = targets.length > 1 ? `${targets.length} passages` : `${count} selected ${count === 1 ? 'word' : 'words'}`;
    for (const prop of ['bold', 'italic', 'underline']) {
      const button = $(`[data-format-toggle="${prop}"]`, toolbar);
      button.disabled = prop !== 'underline' && !hasSelectedProse();
      button.setAttribute('aria-pressed', String(!button.disabled && allStyled(prop)));
      button.title = button.disabled ? 'Select prose to change its font. Formulas retain their mathematical font.' : prop[0].toUpperCase() + prop.slice(1);
    }
  }
  function selectionSnapshot() {
    const selection = window.getSelection(), range = selection?.rangeCount ? selection.getRangeAt(0) : null;
    return {
      direct: window.ReaderDirectEdit?.snapshotSelection(),
      targets: targets.map(target => ({ ...target })),
      backward: Boolean(range && (selection.anchorNode !== range.startContainer || selection.anchorOffset !== range.startOffset))
    };
  }
  function pointAt(root, offset) {
    const p = projection(root); offset = Math.max(0, Math.min(offset, p.text.length));
    for (const unit of p.units) {
      if (offset > unit.end) continue;
      if (unit.type === 'text') return [unit.node, Math.max(0, offset - unit.start)];
      const parent = unit.node.parentNode, index = [...parent.childNodes].indexOf(unit.node);
      return [parent, index + (offset > unit.start ? 1 : 0)];
    }
    return [root, root.childNodes.length];
  }
  function restoreFormattingSelection(snapshot) {
    if (snapshot.direct && window.ReaderDirectEdit?.restoreSelection(snapshot.direct)) {
      const bounds = list => list.map(({ id, start, end }) => ({ id, start, end }));
      if (same(bounds(currentTargets()), bounds(snapshot.targets))) return true;
    }
    // Display equations are not editable roots; retain their selection too.
    const resolved = snapshot.targets.map(target => {
      const root = blockElement(target.id)?.querySelector('.block-content'); if (!root) return null;
      const p = projection(root), range = resolveAnchor(target, p, target.text || p.text);
      return range ? { root, ...range } : null;
    });
    if (!resolved.length || resolved.some(item => !item)) return false;
    const first = resolved[0], last = resolved[resolved.length - 1];
    const start = pointAt(first.root, first.start), end = pointAt(last.root, last.end);
    const selection = window.getSelection(); if (!selection) return false;
    if (first.root.isContentEditable) first.root.focus({ preventScroll: true });
    selection.removeAllRanges();
    if (selection.setBaseAndExtent) selection.setBaseAndExtent(...(snapshot.backward ? end : start), ...(snapshot.backward ? start : end));
    else { const range = document.createRange(); range.setStart(...start); range.setEnd(...end); selection.addRange(range); }
    return true;
  }
  function apply(patch) {
    const resolved = resolvedTargets(); if (!resolved.length) { hideToolbar(); return; }
    const changes = [];
    for (const { block, p, start, end } of resolved) {
      let runs = resolveRanges(block.formatting, p);
      if (patch && ('bold' in patch || 'italic' in patch)) {
        for (const unit of p.units) {
          if (unit.type !== 'text') continue;
          const left = Math.max(start, unit.start), right = Math.min(end, unit.end);
          if (right > left) runs = updateRanges(runs, left, right, patch, p.text.length);
        }
      } else runs = updateRanges(runs, start, end, patch, p.text.length);
      const formatting = runs.length ? { version: 1, text: p.text, ranges: runs.map(run => anchor(p, run.start, run.end, run.style)) } : undefined;
      if (!same(block.formatting, formatting)) changes.push({ block, formatting });
    }
    if (!changes.length) return;
    const rect = toolbar?.getBoundingClientRect(), caret = selectionSnapshot();
    api.onBeforeChange?.(changes.map(item => item.block.id));
    for (const { block, formatting } of changes) { if (formatting) block.formatting = formatting; else delete block.formatting; }
    pinned = true;
    window.getSelection()?.removeAllRanges();
    api.onChange?.(changes.map(item => item.block.id));
    restoreFormattingSelection(caret);
    refreshToolbar();
    if (rect) { toolbar.hidden = false; toolbar.style.left = rect.left + 'px'; toolbar.style.top = rect.top + 'px'; }
  }
  function makeButton(label, text, fn, className = '') {
    const button = document.createElement('button'); button.type = 'button'; button.className = 'format-button ' + className;
    button.setAttribute('aria-label', label); button.title = label; button.textContent = text;
    button.addEventListener('pointerdown', event => event.preventDefault());
    button.addEventListener('click', event => { event.stopPropagation(); fn(); }); return button;
  }
  function makePalette(label, values, property) {
    const row = document.createElement('div'); row.className = 'format-palette'; row.setAttribute('role', 'group'); row.setAttribute('aria-label', label);
    const name = document.createElement('span'); name.className = 'format-palette-label'; name.textContent = label; row.append(name);
    for (const [key, color] of Object.entries(values)) {
      const actionLabel = property === 'color' ? (color ? `Text color: ${key}` : 'Default text color') : (color ? `Highlight: ${key}` : 'Remove highlight');
      const button = makeButton(actionLabel, property === 'color' ? 'A' : '', () => apply({ [property]: color ? key : null }), 'format-swatch');
      button.dataset.formatProperty = property; button.dataset.formatValue = key;
      if (property === 'color') button.style.color = color || '#24332f';
      else if (color) { button.style.setProperty('--swatch', color); button.classList.add('highlight-swatch'); }
      else { button.textContent = '∅'; button.classList.add('clear-swatch'); }
      row.append(button);
    }
    return row;
  }
  function createToolbar() {
    const el = document.createElement('div'); el.id = 'selection-formatting'; el.className = 'selection-formatting'; el.hidden = true;
    el.setAttribute('role', 'toolbar'); el.setAttribute('aria-label', 'Selection formatting');
    const row = document.createElement('div'); row.className = 'format-actions';
    const label = document.createElement('span'); label.id = 'format-selection-label'; label.className = 'format-selection-label'; row.append(label);
    for (const [prop, letter] of [['bold', 'B'], ['italic', 'I'], ['underline', 'U']]) {
      const button = makeButton(prop[0].toUpperCase() + prop.slice(1), letter, () => apply({ [prop]: !allStyled(prop) }), 'format-' + prop);
      button.dataset.formatToggle = prop; button.setAttribute('aria-pressed', 'false'); row.append(button);
    }
    row.append(makeButton('Clear selection formatting', 'Clear', () => apply(null), 'format-clear'));
    row.append(makeButton('Edit selected passage', 'Edit', () => { const id = targets[0]?.id; hideToolbar(); if (id) api.onEdit?.(id); }, 'format-edit'));
    row.append(makeButton('Close formatting toolbar', '×', hideToolbar, 'format-close'));
    el.append(row, makePalette('Text', TEXT, 'color'), makePalette('Highlight', HIGHLIGHT, 'highlight'));
    el.addEventListener('pointerdown', event => event.stopPropagation());
    el.addEventListener('keydown', event => {
      // Keep arrows and Space inside this toolbar; the reader uses them for playback.
      if (event.key === 'Escape') { event.preventDefault(); hideToolbar(); }
      event.stopPropagation();
    });
    document.body.append(el); return el;
  }
  function init(options) {
    api = options || {}; if (initialized) return;
    initialized = true; toolbar = createToolbar();
    document.addEventListener('pointerup', event => { pointerSelecting = false; if (!toolbar.contains(event.target)) queueCapture(); });
    document.addEventListener('keyup', event => { if (!toolbar.contains(event.target)) queueCapture(); });
    document.addEventListener('selectionchange', queueCapture);
    document.addEventListener('pointerdown', event => {
      if (toolbar.contains(event.target)) return;
      pointerSelecting = Boolean(event.target.closest('.block-content')); hideToolbar();
    });
    document.addEventListener('keydown', event => {
      const formField=event.target.closest('input,textarea,select');
      if((event.metaKey||event.ctrlKey)&&event.shiftKey&&event.key.toLowerCase()==='h'&&!formField){
        targets=currentTargets();if(targets.length){event.preventDefault();event.stopImmediatePropagation();if(!event.repeat)apply({highlight:allStyled('highlight','yellow')?null:'yellow'});}return;
      }
      if((event.metaKey||event.ctrlKey)&&['b','i','u'].includes(event.key.toLowerCase())&&!formField)targets=currentTargets();
      if (event.key === 'Escape' && !toolbar.hidden && !event.target.closest('input,textarea,select')) { hideToolbar(); return; }
      if ((event.metaKey || event.ctrlKey) && ['b', 'i', 'u'].includes(event.key.toLowerCase()) && targets.length && !event.target.closest('input,textarea,select')) {
        const prop = { b: 'bold', i: 'italic', u: 'underline' }[event.key.toLowerCase()];
        event.preventDefault(); event.stopImmediatePropagation(); if (prop !== 'underline' && !hasSelectedProse()) return; apply({ [prop]: !allStyled(prop) });
      }
      if (event.altKey && event.shiftKey && event.key.toLowerCase() === 'a' && targets.length) { event.preventDefault(); toolbar.querySelector('button')?.focus(); }
    }, true);
    window.addEventListener('resize', () => { if (!toolbar.hidden) placeToolbar(rectForTargets()); });
    window.addEventListener('scroll', () => { if (toolbar.hidden) return; if (pinned) hideToolbar(); else placeToolbar(rectForTargets()); }, { passive: true });
  }
  window.ReaderFormatting = { init, decorate, clearSelection,
    // Pure helpers permit interval and anchor regressions without a live reader.
    testing: { cleanStyle, resolveAnchor, resolveRanges, updateRanges, anchor, atStyle, rangeHasStyle, projection, selectionPart }
  };
})();
