// Minimal IndexedDB wrapper for the proply polar cache (no dependencies).
//
// One object store ("polars") keyed by the polar cache key
// "<foil hash>|<reynolds>|<mach>", each value {alpha, cl, cd} arrays — the
// same records the Rust PolarStore keeps.  All failures are expected to be
// handled by the caller: the cache is a pure performance artifact, so a
// missing/blocked database just means a cold design run.

const DB_NAME = 'proply-polar-cache';
const STORE = 'polars';

function openDb() {
  return new Promise((resolve, reject) => {
    const req = indexedDB.open(DB_NAME, 1);
    req.onupgradeneeded = () => {
      if (!req.result.objectStoreNames.contains(STORE)) {
        req.result.createObjectStore(STORE);
      }
    };
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error);
  });
}

/** All stored records as [key, {alpha, cl, cd}] pairs. */
export async function loadAllPolars() {
  const db = await openDb();
  return new Promise((resolve, reject) => {
    const req = db.transaction(STORE, 'readonly').objectStore(STORE).openCursor();
    const out = [];
    req.onsuccess = () => {
      const cur = req.result;
      if (cur) {
        out.push([cur.key, cur.value]);
        cur.continue();
      } else {
        resolve(out);
      }
    };
    req.onerror = () => reject(req.error);
  });
}

/** Store [key, {alpha, cl, cd}] pairs in one transaction. */
export async function putPolars(entries) {
  const db = await openDb();
  return new Promise((resolve, reject) => {
    const tx = db.transaction(STORE, 'readwrite');
    const store = tx.objectStore(STORE);
    for (const [key, polar] of entries) store.put(polar, key);
    tx.oncomplete = () => resolve();
    tx.onerror = () => reject(tx.error);
  });
}

/** Drop every cached polar. */
export async function clearPolars() {
  const db = await openDb();
  return new Promise((resolve, reject) => {
    const tx = db.transaction(STORE, 'readwrite');
    tx.objectStore(STORE).clear();
    tx.oncomplete = () => resolve();
    tx.onerror = () => reject(tx.error);
  });
}
