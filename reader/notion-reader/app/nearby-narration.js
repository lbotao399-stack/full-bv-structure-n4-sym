/* Bounded narration look-ahead. It never changes the document or starts audio. */
(function (root) {
  'use strict';
  function create({ getCandidates, prepare, fingerprint = b => JSON.stringify(b), isCurrent = () => true, onReady = () => {} }) {
    let anchor = null, stopped = true, wanted = [], running = null;
    const attempted = new Map();
    const clone = value => JSON.parse(JSON.stringify(value));
    const current = () => !stopped && isCurrent(anchor);
    const stillWanted = task => wanted.some(item => item.id === task.id && item.key === task.key);

    async function pump() {
      if (running || !current()) return;
      const next = wanted.find(item => attempted.get(item.id) !== item.key);
      if (!next) return;
      const task = { ...next, controller: new AbortController() };
      running = task;
      attempted.set(task.id, task.key);
      try {
        const result = await prepare(task.block, task.controller.signal, { background: true });
        if (!task.controller.signal.aborted && current() && stillWanted(task)) {
          await onReady(task.block, result, task.controller.signal);
        }
      } catch (_error) {
        // Failed/aborted preloads leave foreground reading available and do not retry in a loop.
      } finally {
        if (running === task) running = null;
        // A completed or cancelled task opens one slot; never create an unbounded queue.
        void pump();
      }
    }

    function update(activeId) {
      anchor = activeId;
      stopped = !activeId;
      const ids = new Set();
      wanted = [];
      if (!stopped) {
        for (const source of getCandidates(activeId) || []) {
          if (!source?.id || source.id === activeId || ids.has(source.id)) continue;
          ids.add(source.id);
          const block = clone(source);
          wanted.push({ id: block.id, key: fingerprint(block), block });
          if (wanted.length === 21) break;
        }
      }
      for (const [id, key] of attempted) {
        if (!wanted.some(item => item.id === id && item.key === key)) attempted.delete(id);
      }
      if (running && (!current() || !stillWanted(running))) running.controller.abort();
      void pump();
    }

    function stop() {
      stopped = true;
      anchor = null;
      wanted = [];
      attempted.clear();
      running?.controller.abort();
    }

    function invalidate(id) {
      // Remove the obsolete snapshot. Only a later playback update can re-add its new wording.
      wanted = wanted.filter(item => item.id !== id);
      attempted.delete(id);
      if (running?.id === id) running.controller.abort();
    }

    return { update, stop, invalidate };
  }
  root.ReaderNearbyNarration = { create };
})(typeof window === 'undefined' ? globalThis : window);
