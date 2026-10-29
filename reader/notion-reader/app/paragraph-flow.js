/* Paragraph flow.
   Import splits prose into one editable passage per sentence; each passage
   keeps the source paragraph it came from (block.sentenceSplit.originId, set on
   the element as data-paragraph by app.js). This view-only module marks a
   passage whose preceding passage belongs to the same source paragraph with
   .flow-continued, so theme.css can set the sentences as one paragraph.
   The document model, block IDs and editing behaviour are untouched. */
(function () {
  'use strict';

  const CONTINUED = 'flow-continued';

  function previousBlock(block) {
    let node = block.previousElementSibling;
    while (node && !node.classList.contains('source-block')) node = node.previousElementSibling;
    return node;
  }

  function mark(root) {
    for (const block of root.querySelectorAll('.source-block')) {
      const paragraph = block.dataset.paragraph;
      const continued = Boolean(paragraph) && previousBlock(block)?.dataset.paragraph === paragraph;
      if (block.classList.contains(CONTINUED) !== continued) block.classList.toggle(CONTINUED, continued);
    }
  }

  // Only additions or removals of whole passages can change the flow; text
  // edits inside a passage are ignored.
  const touchesBlocks = record => [...record.addedNodes, ...record.removedNodes]
    .some(node => node.nodeType === 1 && node.classList.contains('source-block'));

  function start() {
    const root = document.getElementById('blocks');
    if (!root) return;
    let queued = false;
    const schedule = () => {
      if (queued) return;
      queued = true;
      requestAnimationFrame(() => { queued = false; mark(root); });
    };
    new MutationObserver(records => { if (records.some(touchesBlocks)) schedule(); })
      .observe(root, { childList: true, subtree: true });
    schedule();
  }

  window.ReaderParagraphFlow = { refresh: () => { const root = document.getElementById('blocks'); if (root) mark(root); } };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start);
  else start();
})();
