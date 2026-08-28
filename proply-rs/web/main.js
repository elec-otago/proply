// proply-rs browser demo: design a propeller entirely in the browser.
//
// The WebAssembly module (built by `make wasm` into ../pkg/) exposes a
// long-lived PropSession: its polar cache is hydrated from IndexedDB at
// startup, designs run against it, and freshly simulated polars are
// written back after each run.  The design call is synchronous and can
// block the tab for a while (longer without "plate polars": real XFOIL
// polars are computed in-wasm).

import init, { PropSession } from '../pkg/proply_rs.js';
import { loadAllPolars, putPolars, clearPolars } from './idb.js';

const $ = (id) => document.getElementById(id);

const DEFAULT_PARAMS = {
  name: 'browser_demo',
  altitude: 0.0,
  forward_airspeed: 3.0,
  motor_Kv: 1900,
  motor_volts: 11.0,
  motor_no_load_current: 0.5,
  motor_winding_resistance: 0.405,
  blades: 3,
  thrust: '3 N',
  radius: '68 mm',
  tip_chord: '5 mm',
  scimitar_percent: 0.0,
  trailing_edge: '0.25 mm',
  center_hole: '1.5 mm',
  hub_radius: '6 mm',
  hub_depth: '6 mm',
  chord_spline_n: 3,
  bem: true,
  resolution: 10,
  n: 12,
};

let session = null;

function setStatus(text) {
  $('status').textContent = text;
}

async function boot() {
  setStatus('loading WebAssembly…');
  await init();
  session = new PropSession();

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
  setStatus(`ready — cache: ${hydrated} polars hydrated`);
  $('design').disabled = false;
}

function runDesign() {
  let params;
  try {
    params = JSON.parse($('params').value);
  } catch (e) {
    setStatus(`invalid design JSON: ${e}`);
    return;
  }
  params.plate = $('plate').checked;

  $('design').disabled = true;
  $('yaml').textContent = '';
  $('step-link').replaceChildren();
  setStatus('designing… (the tab is busy while this runs)');
  // Let the status paint before the synchronous design call blocks.
  setTimeout(async () => {
    const t0 = performance.now();
    try {
      const out = session.design(JSON.stringify(params));
      const seconds = (performance.now() - t0) / 1000;

      $('yaml').textContent = out.yaml;
      const blob = new Blob([out.step], { type: 'text/plain' });
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = `${params.name || 'prop'}.step`;
      a.textContent = `download ${params.name || 'prop'}.step (${out.step.length} bytes)`;
      $('step-link').replaceChildren(a);

      // Persist the polars this run simulated.
      let stored = 0;
      const entries = Object.entries(JSON.parse(session.take_new_json()));
      if (entries.length > 0) {
        try {
          await putPolars(entries);
          stored = entries.length;
        } catch (e) {
          console.warn('could not persist new polars:', e);
        }
      }
      setStatus(
        `designed in ${seconds.toFixed(1)} s — thrust ${out.thrust.toFixed(2)} N, ` +
        `torque ${out.torque.toFixed(3)} N·m at ${out.rpm.toFixed(0)} rpm; ` +
        `${stored} new polars cached (${session.polar_count()} total)`,
      );
    } catch (e) {
      setStatus(`design failed: ${e}`);
      console.error(e);
    } finally {
      $('design').disabled = false;
    }
  }, 50);
}

async function clearCache() {
  try {
    await clearPolars();
    session = new PropSession(); // drop the hydrated store too
    setStatus('cache cleared — next design runs cold');
  } catch (e) {
    setStatus(`could not clear cache: ${e}`);
  }
}

$('params').value = JSON.stringify(DEFAULT_PARAMS, null, 2);
$('design').addEventListener('click', runDesign);
$('clear-cache').addEventListener('click', clearCache);
boot();
