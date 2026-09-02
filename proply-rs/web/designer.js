// The design worker: owns the WebAssembly module, the polar-cache session
// and its IndexedDB persistence, so a design run never blocks the page
// (a cold XFOIL design takes minutes; the main thread stays interactive).

import init, { PropSession } from './pkg/proply_rs.js';
import { loadAllPolars, putPolars, clearPolars } from './idb.js';

let session = null;

// A fresh design session with the per-polar persistence hook installed:
// every polar the moment it is calculated — a good sweep or a failed one
// (cached as a degenerate marker) — is handed to the host, which writes
// it into IndexedDB immediately.  A design interrupted mid-way (the page
// closed, the worker killed) therefore keeps every completed sweep,
// exactly like the CLI's per-polar disk checkpoint.  The end-of-design
// drain in runDesign remains as the catch-up flush.
function makeSession() {
  const s = new PropSession();
  s.set_on_polar((key, alpha, cl, cd) => {
    putPolars([[key, { alpha: Array.from(alpha), cl: Array.from(cl), cd: Array.from(cd) }]])
      .catch((e) => console.warn('could not persist polar:', e));
  });
  return s;
}

async function boot() {
  await init();
  session = makeSession();
  let hydrated = 0;
  try {
    for (const [key, polar] of await loadAllPolars()) {
      session.hydrate_entry(
        key,
        Float64Array.from(polar.alpha),
        Float64Array.from(polar.cl),
        Float64Array.from(polar.cd),
      );
      hydrated++;
    }
  } catch (e) {
    console.warn('polar cache unavailable, starting cold:', e);
  }
  postMessage({ type: 'ready', hydrated });
}

async function runDesign(params) {
  const t0 = performance.now();
  try {
    const out = session.design(JSON.stringify(params));
    let newPolars = 0;
    const entries = Object.entries(JSON.parse(session.take_new_json()));
    if (entries.length > 0) {
      try {
        await putPolars(entries);
        newPolars = entries.length;
      } catch (e) {
        console.warn('could not persist new polars:', e);
      }
    }
    postMessage({
      type: 'design-complete',
      elapsed: (performance.now() - t0) / 1000,
      thrust: out.thrust,
      torque: out.torque,
      rpm: out.rpm,
      yaml: out.yaml,
      step: out.step,
      name: params.name || 'prop',
      warning: out.warning || '',
      newPolars,
      totalPolars: session.polar_count(),
    });
  } catch (e) {
    postMessage({ type: 'design-error', message: String(e) });
  }
}

self.onmessage = (ev) => {
  const msg = ev.data;
  if (session === null) {
    return; // still booting the wasm module
  }
  if (msg.type === 'design') {
    runDesign(msg.params);
  } else if (msg.type === 'clear-cache') {
    clearPolars()
      .then(() => {
        session = makeSession(); // drop the hydrated store too
        postMessage({ type: 'cache-cleared' });
      })
      .catch((e) => {
        postMessage({ type: 'design-error', message: `could not clear cache: ${e}` });
      });
  }
};

boot();
