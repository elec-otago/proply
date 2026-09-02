// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! WebAssembly bindings: proply-rs in the browser.
//!
//! The browser entry point is a long-lived [`PropSession`] holding the
//! polar cache in memory.  The host hydrates it at startup from its own
//! storage (e.g. IndexedDB) and, after each design, persists the polars
//! that were freshly simulated — there is no filesystem to save to.  The
//! design itself runs the same [`crate::pipeline`] the CLI runs.

use std::sync::{Arc, Mutex};

use js_sys::Float64Array;
use wasm_bindgen::prelude::*;
use wasm_bindgen::JsValue;

use crate::cache::{PolarStore, StoredPolar};
use crate::design_parameters::DesignParameters;
use crate::pipeline::{self, DesignOutcome};

/// One finished design: the STEP (AP242) document, the YAML summary and
/// the headline numbers.
#[wasm_bindgen(getter_with_clone)]
pub struct DesignOutput {
    pub step: String,
    pub yaml: String,
    pub thrust: f64,
    pub torque: f64,
    pub rpm: f64,
    pub power: f64,
    /// Empty when the design reached its operating point; otherwise an
    /// explicit note describing the closest achievable design.
    pub warning: String,
}

impl From<DesignOutcome> for DesignOutput {
    fn from(o: DesignOutcome) -> Self {
        DesignOutput {
            step: o.step,
            yaml: o.yaml,
            thrust: o.thrust,
            torque: o.torque,
            rpm: o.rpm,
            power: o.power,
            warning: o.warning.unwrap_or_default(),
        }
    }
}

/// A design session with a warm polar cache, kept across design calls.
#[wasm_bindgen]
pub struct PropSession {
    store: Arc<Mutex<PolarStore>>,
}

impl Default for PropSession {
    fn default() -> Self {
        Self::new()
    }
}

#[wasm_bindgen]
impl PropSession {
    #[wasm_bindgen(constructor)]
    pub fn new() -> PropSession {
        console_error_panic_hook::set_once();
        PropSession {
            store: Arc::new(Mutex::new(PolarStore::in_memory())),
        }
    }

    /// Number of cached polars (e.g. after startup hydration).
    pub fn polar_count(&self) -> usize {
        self.store.lock().unwrap().len()
    }

    /// Insert one pre-existing polar record — the bulk startup hydration
    /// path.  Records inserted here are not reported by
    /// [`PropSession::take_new_json`] (the host already has them).
    pub fn hydrate_entry(&self, key: &str, alpha: &Float64Array, cl: &Float64Array, cd: &Float64Array) {
        let polar = StoredPolar {
            alpha: alpha.to_vec(),
            cl: cl.to_vec(),
            cd: cd.to_vec(),
        };
        self.store.lock().unwrap().hydrate(key.to_string(), polar);
    }

    /// Hydrate from a full cache document ([`PropSession::cache_to_json`]
    /// format), replacing any current contents.
    pub fn hydrate_json(&self, json: &str) {
        *self.store.lock().unwrap() = PolarStore::from_json_str(json);
    }

    /// Install the host's per-polar persistence hook: called synchronously
    /// for every polar the moment it is freshly calculated — a good sweep
    /// or a degenerate failure marker — with the cache key and the
    /// (alpha, cl, cd) arrays.  The host writes each record to its
    /// IndexedDB cache immediately, so a design interrupted mid-way keeps
    /// every completed sweep, exactly like the native CLI's per-polar disk
    /// checkpoint.  The hook replaces any previously installed one.
    pub fn set_on_polar(&self, on_polar: js_sys::Function) {
        let mut store = self.store.lock().unwrap();
        store.set_on_insert(Box::new(move |key: &str, p: &StoredPolar| {
            let key = JsValue::from_str(key);
            let alpha = Float64Array::from(p.alpha.as_slice());
            let cl = Float64Array::from(p.cl.as_slice());
            let cd = Float64Array::from(p.cd.as_slice());
            let _ = on_polar.call4(&JsValue::UNDEFINED, &key, &alpha, &cl, &cd);
        }));
    }

    /// The whole cache as a JSON document (export/migration escape hatch).
    pub fn cache_to_json(&self) -> String {
        self.store
            .lock()
            .unwrap()
            .to_json_string()
            .unwrap_or_else(|| "{}".into())
    }

    /// Run one full design from JSON design parameters (the same format
    /// and validation as the CLI's `--param` file).
    pub fn design(&self, params_json: String) -> Result<DesignOutput, JsValue> {
        let param =
            DesignParameters::from_json(&params_json).map_err(|e| JsValue::from_str(&e))?;
        if !param.bem && !param.lifting_line {
            return Err(JsValue::from_str(
                "select a design loop (set `bem` or `lifting_line` in the design JSON)",
            ));
        }
        if param.cst && param.arad {
            return Err(JsValue::from_str(
                "choose one foil family (naca, cst or arad)",
            ));
        }
        let outcome =
            pipeline::run_design(&param, self.store.clone()).map_err(|e| JsValue::from_str(&e))?;
        Ok(outcome.into())
    }

    /// A JSON map `key -> {alpha, cl, cd}` of the polars simulated since
    /// the last call — what the host should persist (e.g. into IndexedDB).
    /// The session keeps every polar for future designs.
    pub fn take_new_json(&self) -> String {
        let mut store = self.store.lock().unwrap();
        let new = store.take_new_entries();
        serde_json::to_string(&new).unwrap_or_else(|_| "{}".into())
    }
}
