// proply-rs browser demo: design a propeller entirely in the browser.
//
// The WebAssembly designer runs in a dedicated worker (designer.js) —
// this file only sends design requests and renders the results, so the
// page stays responsive while a design runs (longer without "plate
// polars": real XFOIL polars are computed in-wasm).
//
// The 3D preview (viewer.js) loads on demand: three.js and
// occt-import-js only start downloading on the first "Show 3D preview"
// click, not on every design.

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

let worker = null;
let busyTimer = null;
let viewer = null; // viewer.js module, imported on the first preview click
let latestStep = null; // STEP text of the most recent completed design

async function viewerModule() {
  if (!viewer) {
    viewer = await import('./viewer.js');
  }
  return viewer;
}

function setStatus(text) {
  $('status-text').textContent = text;
}

/// Show/hide the spinner.  While busy, the status line ticks with the
/// elapsed time (the design itself runs silently in the worker).
function setBusy(busy, label = '') {
  $('spinner').style.visibility = busy ? 'visible' : 'hidden';
  if (busy) {
    const t0 = performance.now();
    setStatus(label);
    busyTimer = setInterval(() => {
      setStatus(`${label} — ${((performance.now() - t0) / 1000).toFixed(1)} s elapsed`);
    }, 100);
  } else if (busyTimer !== null) {
    clearInterval(busyTimer);
    busyTimer = null;
  }
}

function boot() {
  setStatus('loading WebAssembly (worker)…');
  worker = new Worker('designer.js', { type: 'module' });
  worker.onmessage = (ev) => {
    const msg = ev.data;
    if (msg.type === 'ready') {
      setStatus(`ready — cache: ${msg.hydrated} polars hydrated`);
      $('design').disabled = false;
    } else if (msg.type === 'design-complete') {
      render(msg);
    } else if (msg.type === 'design-error') {
      setBusy(false);
      setStatus(`design failed: ${msg.message}`);
      console.error(msg.message);
      const note = $('viewer-note');
      if (note) {
        note.textContent = 'Run a design to preview the propeller here.';
        note.style.display = '';
      }
      $('design').disabled = false;
    } else if (msg.type === 'cache-cleared') {
      setStatus('cache cleared — next design runs cold');
    }
  };
  worker.onerror = (e) => {
    setBusy(false);
    setStatus(`worker failed: ${e.message || e}`);
    $('design').disabled = true;
  };
}

function render(msg) {
  $('yaml').textContent = msg.yaml;
  const blob = new Blob([msg.step], { type: 'text/plain' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `${msg.name}.step`;
  a.textContent = `download ${msg.name}.step (${msg.step.length} bytes)`;
  $('step-link').replaceChildren(a);
  latestStep = msg.step;
  const note = $('viewer-note');
  if (note) {
    note.textContent = 'Click "Show 3D preview" to view the propeller in 3D.';
    note.style.display = '';
  }
  $('preview').disabled = false;
  setBusy(false);
  setStatus(
    `done in ${msg.elapsed.toFixed(1)} s — thrust ${msg.thrust.toFixed(2)} N, ` +
    `torque ${msg.torque.toFixed(3)} N·m at ${msg.rpm.toFixed(0)} rpm, ` +
    `power at the operating point ${(msg.torque * msg.rpm * 2 * Math.PI / 60).toFixed(1)} W; ` +
    `${msg.newPolars} new polars cached (${msg.totalPolars} total)`,
  );
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
  $('preview').disabled = true;
  latestStep = null;
  $('yaml').textContent = '';
  $('step-link').replaceChildren();
  if (viewer) {
    viewer.clearStep(); // the on-screen model no longer matches this run
  }
  const note = $('viewer-note');
  if (note) {
    note.textContent = 'designing…';
    note.style.display = '';
  }
  setBusy(true, 'designing… (worker busy; the page stays responsive)');
  worker.postMessage({ type: 'design', params });
}

function clearCache() {
  worker.postMessage({ type: 'clear-cache' });
}

/// Tessellate and display the most recent design's STEP in the 3D window.
/// The viewer module (and with it three.js + occt-import-js) is only
/// loaded on this first explicit request.
async function showPreview() {
  const btn = $('preview');
  btn.disabled = true;
  const note = $('viewer-note');
  note.textContent = 'tessellating the STEP model…';
  note.style.display = '';
  try {
    const v = await viewerModule();
    await v.showStep(latestStep, $('viewer'), note);
  } finally {
    btn.disabled = false;
  }
}

$('params').value = JSON.stringify(DEFAULT_PARAMS, null, 2);
$('design').addEventListener('click', runDesign);
$('clear-cache').addEventListener('click', clearCache);
$('preview').addEventListener('click', showPreview);
boot();
