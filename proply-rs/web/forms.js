// proply-rs browser demo: tabbed editors for the design JSON.
//
// The design parameters are one flat JSON document (the same schema the
// CLI and the wasm parse).  The tabs edit the descriptors of three named
// parts of that document — the propeller itself, the electric motor
// model, and an "other" motor given directly by its operating point.
// The composed JSON stays in the textarea (the source of truth for
// running a design) and is persisted in localStorage, so a design
// survives page reloads.
//
// The pure parts (TABS, composeDesign, syncForm, readForm) run anywhere
// — node included — so the round-trip logic is testable without a
// browser; buildForm touches the DOM and is only used from main.js.

export const STORAGE_KEY = 'proply-design-params';

// Field definitions, one entry per JSON descriptor.
//
// type "quantity": a dimensioned value.  The input holds the value in
// `unit` (the DEFAULT_PARAMS style); composeDesign passes the text
// through unchanged, so both suffixed strings ("68 mm") and bare
// numbers (the JSON's SI units) are accepted by the wasm.  `bareToDisplay`
// converts a bare-number JSON value into the display unit when a form is
// synced from JSON (radius: 0.068 -> "68 mm"; trailing_edge is already
// stored in mm, so it converts by 1).  `optional` fields are omitted
// from the JSON when blank, so the Rust defaults apply.
export const TABS = [
  {
    id: 'prop',
    label: 'Propeller Specifications',
    fields: [
      { key: 'name', label: 'Name', type: 'text' },
      { key: 'blades', label: 'Blades', type: 'number', min: 2 },
      { key: 'radius', label: 'Radius', type: 'quantity', unit: 'mm', bareToDisplay: 1000 },
      { key: 'thrust', label: 'Thrust', type: 'quantity', unit: 'N', bareToDisplay: 1 },
      { key: 'tip_chord', label: 'Tip chord', type: 'quantity', unit: 'mm', bareToDisplay: 1000 },
      { key: 'hub_radius', label: 'Hub radius', type: 'quantity', unit: 'mm', bareToDisplay: 1000 },
      { key: 'hub_depth', label: 'Hub depth', type: 'quantity', unit: 'mm', bareToDisplay: 1000 },
      { key: 'center_hole', label: 'Center hole', type: 'quantity', unit: 'mm', bareToDisplay: 1000, optional: true },
      { key: 'trailing_edge', label: 'Trailing edge', type: 'quantity', unit: 'mm', bareToDisplay: 1 },
      { key: 'scimitar_percent', label: 'Scimitar (%)', type: 'number' },
      { key: 'forward_airspeed', label: 'Forward airspeed (m/s)', type: 'number' },
      { key: 'altitude', label: 'Altitude (m)', type: 'number' },
    ],
  },
  {
    id: 'electric',
    label: 'Electric Motor',
    fields: [
      { key: 'motor_Kv', label: 'Kv (RPM/V)', type: 'number' },
      { key: 'motor_volts', label: 'Voltage (V)', type: 'number' },
      { key: 'motor_no_load_current', label: 'No-load current (A)', type: 'number' },
      { key: 'motor_winding_resistance', label: 'Winding resistance (Ω)', type: 'number' },
    ],
  },
  {
    id: 'other',
    label: 'Other Motor',
    fields: [
      { key: 'motor_torque', label: 'Torque (N·m)', type: 'number', optional: true },
      { key: 'motor_RPM', label: 'RPM', type: 'number', optional: true },
    ],
  },
];

function allFields() {
  return TABS.flatMap((tab) => tab.fields);
}

/// Merge the tab field values into `current` (a parsed design JSON
/// carrying anything the tabs do not own — run options, camber, ...).
/// Blank values drop their key so the Rust defaults apply; number
/// fields are written as JSON numbers.
export function composeDesign(current, values) {
  const out = { ...current };
  for (const field of allFields()) {
    const v = values[field.key];
    if (v === '' || v === null || v === undefined) {
      delete out[field.key];
      continue;
    }
    if (field.type === 'number') {
      const n = Number(v);
      if (Number.isFinite(n)) out[field.key] = n;
      else delete out[field.key];
    } else {
      out[field.key] = String(v).trim();
    }
  }
  return out;
}

function fmt(n) {
  return String(Number(n.toFixed(4)));
}

/// Populate the inputs from a parsed design JSON.  Quantity fields
/// stored as bare numbers (SI) are converted into the display unit;
/// suffixed strings are shown exactly as stored.
export function syncForm(inputs, json) {
  for (const field of allFields()) {
    const el = inputs[field.key];
    if (!el) continue;
    const v = json[field.key];
    if (v === undefined || v === null) {
      el.value = '';
    } else if (field.type === 'quantity' && typeof v === 'number') {
      el.value = `${fmt(v * field.bareToDisplay)} ${field.unit}`;
    } else {
      el.value = String(v);
    }
  }
}

/// Read every input as {key: value}; number inputs give the raw string
/// (composeDesign converts), so an empty field stays distinguishable.
export function readForm(inputs) {
  const out = {};
  for (const field of allFields()) {
    const el = inputs[field.key];
    if (el) out[field.key] = el.value;
  }
  return out;
}

/// Build the tab bar and field panels into `container` (which must be
/// empty).  Returns the input map (key -> element) for readForm/syncForm.
export function buildForm(container) {
  const bar = document.createElement('div');
  bar.className = 'tabs';
  const panels = {};
  const inputs = {};

  for (const tab of TABS) {
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'tab-btn';
    btn.textContent = tab.label;
    const panel = document.createElement('div');
    panel.className = 'tab-panel';
    const grid = document.createElement('div');
    grid.className = 'field-grid';
    for (const field of tab.fields) {
      const label = document.createElement('label');
      label.className = 'field';
      const span = document.createElement('span');
      span.textContent = field.label;
      const input = document.createElement('input');
      input.type = field.type === 'number' ? 'number' : 'text';
      if (field.type === 'number' && field.min !== undefined) {
        input.min = field.min;
      }
      inputs[field.key] = input;
      label.append(span, input);
      grid.append(label);
    }
    panel.append(grid);
    const activate = () => {
      for (const b of bar.children) b.classList.toggle('active', b === btn);
      for (const p of Object.values(panels)) p.classList.toggle('active', p === panel);
    };
    btn.addEventListener('click', activate);
    bar.append(btn);
    panels[tab.id] = panel;
    container.append(panel);
  }
  container.prepend(bar);
  bar.firstChild.classList.add('active');
  panels[TABS[0].id].classList.add('active');
  return inputs;
}

/// The last design JSON held in localStorage, or null when absent or
/// unreadable (the caller falls back to DEFAULT_PARAMS).
export function loadStored() {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (!raw) return null;
    const parsed = JSON.parse(raw);
    return parsed && typeof parsed === 'object' && !Array.isArray(parsed) ? parsed : null;
  } catch {
    return null;
  }
}

/// Persist the current design JSON.  Failures are silent: storage
/// unavailable (private mode, quota) only costs persistence, never the
/// current session.
export function saveStored(params) {
  try {
    localStorage.setItem(STORAGE_KEY, JSON.stringify(params));
  } catch {
    // ignore
  }
}
