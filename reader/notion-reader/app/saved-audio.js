/* Bounded persistent audio bytes; quota/private-mode failures remain optional. */
(() => {
  'use strict';
  const name = 'reader-audio-v1', limit = 256 * 1024 * 1024, count = 128;
  let writes = Promise.resolve();
  async function get(url) {
    try {
      const cache = await caches.open(name), response = await cache.match(url);
      if (!response || !response.ok) return null;
      return await response.blob();
    } catch (_) { return null; }
  }
  function put(url, blob) {
    if (!blob.size || blob.size > 32 * 1024 * 1024) return Promise.resolve();
    writes = writes.catch(() => {}).then(async () => {
      try {
        const cache = await caches.open(name);
        await cache.delete(url);
        await cache.put(url, new Response(blob, {headers: {'Content-Type': 'audio/wav', 'Content-Length': String(blob.size)}}));
        const keys = await cache.keys();
        let size = 0;
        for (let i = keys.length - 1; i >= 0; i--) {
          const response = await cache.match(keys[i]);
          size += Number(response.headers.get('Content-Length')) || (await response.blob()).size;
          if (size > limit || keys.length - i > count) await cache.delete(keys[i]);
        }
      } catch (_) { /* Cache eviction/quota must never prevent reading. */ }
    });
    return writes;
  }
  window.ReaderSavedAudio = {get, put};
})();
