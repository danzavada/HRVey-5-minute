/* HRVey 5-min service worker — offline app shell.
 *
 * App code (index.html, ecg-core.js, manifest) is NETWORK-FIRST so a deploy
 * reaches returning users on their next load; the cache is only the offline
 * fallback. Vendor libraries and icons are cache-first (large, rarely change) —
 * bump CACHE whenever one of them does, which also drops the old cache.
 * (v1 served ecg-core.js cache-first under a fixed name, so browsers that had
 * visited before the row-detection fix kept running the old core forever.)
 */
const CACHE = "hrvey5-v2";
const NETWORK_FIRST = ["./", "./index.html", "./ecg-core.js", "./manifest.webmanifest"];
const CACHE_FIRST = [
  "./vendor/pdf.min.js",
  "./vendor/pdf.worker.min.js",
  "./vendor/xlsx.full.min.js",
  "./icons/icon-192.png",
  "./icons/icon-512.png",
  "./icons/icon-maskable-512.png",
  "./icons/favicon-16.png",
  "./icons/favicon-32.png",
  "./icons/favicon-180.png"
];

self.addEventListener("install", e => {
  e.waitUntil(caches.open(CACHE).then(c => c.addAll(NETWORK_FIRST.concat(CACHE_FIRST))).then(() => self.skipWaiting()));
});
self.addEventListener("activate", e => {
  e.waitUntil(
    caches.keys()
      .then(keys => Promise.all(keys.filter(k => k !== CACHE).map(k => caches.delete(k))))
      .then(() => self.clients.claim())
  );
});
function cachePut(req, resp){ if (resp && resp.ok) { const c = resp.clone(); caches.open(CACHE).then(x => x.put(req, c)).catch(() => {}); } return resp; }
function isAppCode(req){
  const path = new URL(req.url).pathname;
  return req.mode === "navigate" || /\/(index\.html|ecg-core\.js|manifest\.webmanifest)$/.test(path) || path.endsWith("/");
}
self.addEventListener("fetch", e => {
  if (e.request.method !== "GET") return;
  const req = e.request;
  if (isAppCode(req)) {
    e.respondWith(
      fetch(req).then(r => cachePut(req, r))
        .catch(() => caches.match(req, { ignoreSearch: true }).then(m => m || caches.match("./index.html")))
    );
  } else {
    e.respondWith(caches.match(req).then(hit => hit || fetch(req).then(r => cachePut(req, r))));
  }
});
