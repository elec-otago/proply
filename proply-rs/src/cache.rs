// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Persistent polar cache, replacing the Python sqlite `foil_simulator.db`.
//!
//! Polars are stored keyed by `"<foil hash>|<reynolds>|<mach>"` in a single
//! JSON file in the working directory.  Alpha values are stored in radians,
//! matching what the Python database kept.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// One stored polar sweep: alpha (radians), cl, cd per point.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StoredPolar {
    pub alpha: Vec<f64>,
    pub cl: Vec<f64>,
    pub cd: Vec<f64>,
}

/// An in-memory map of cached polars with lazy persistence.
#[derive(Default)]
pub struct PolarStore {
    path: String,
    foils: HashMap<String, StoredPolar>,
    dirty: bool,
    /// Keys inserted since the last [`PolarStore::take_new_entries`] drain —
    /// freshly simulated polars a host without a filesystem (the browser)
    /// persists itself.
    new_keys: Vec<String>,
    /// The per-polar persistence hook, fired synchronously by
    /// [`PolarStore::insert`] for every freshly simulated polar — good or a
    /// degenerate failure marker.  Native hosts checkpoint to disk inside
    /// `insert` and leave this unset; the browser host (wasm.rs) installs
    /// it so each new polar is handed to the IndexedDB cache the moment it
    /// is calculated, instead of only when a design finishes.  Present only
    /// where it can be set: wasm builds and the crate's own tests.
    #[cfg(any(test, target_arch = "wasm32"))]
    on_insert: Option<Box<dyn Fn(&str, &StoredPolar) + Send + Sync>>,
}

// Manual Debug (the persist hook is a closure): `Mutex<PolarStore>`'s
// `lock().unwrap()` needs it.
impl std::fmt::Debug for PolarStore {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PolarStore")
            .field("path", &self.path)
            .field("entries", &self.foils.len())
            .field("dirty", &self.dirty)
            .field("pending", &self.new_keys.len())
            .finish_non_exhaustive()
    }
}

pub fn cache_key(hash: &str, reynolds: f64, mach: f64) -> String {
    format!("{}|{}|{}", hash, reynolds, mach)
}

impl PolarStore {
    /// An empty store that never persists anywhere (in-memory sessions,
    /// e.g. the WebAssembly build).
    pub fn in_memory() -> Self {
        Self::default()
    }

    /// Load the cache from `path` (missing file = empty cache).
    pub fn load(path: &str) -> Self {
        let foils = std::fs::read_to_string(path)
            .ok()
            .and_then(|s| Self::parse(&s))
            .unwrap_or_default();
        Self {
            path: path.to_string(),
            foils,
            dirty: false,
            new_keys: Vec::new(),
            #[cfg(any(test, target_arch = "wasm32"))]
            on_insert: None,
        }
    }

    fn parse(json: &str) -> Option<HashMap<String, StoredPolar>> {
        serde_json::from_str(json).ok()
    }

    /// Hydrate from a cache document produced by [`PolarStore::to_json_string`]
    /// (the same format as the on-disk cache file).  The store keeps no path,
    /// so it never persists.  Invalid JSON yields an empty store, matching
    /// [`PolarStore::load`]'s tolerance of a corrupt cache file.
    pub fn from_json_str(json: &str) -> Self {
        Self {
            path: String::new(),
            foils: Self::parse(json).unwrap_or_default(),
            dirty: false,
            new_keys: Vec::new(),
            #[cfg(any(test, target_arch = "wasm32"))]
            on_insert: None,
        }
    }

    /// The cache contents as a JSON document (the on-disk format).  `None`
    /// only if serialization fails.
    pub fn to_json_string(&self) -> Option<String> {
        serde_json::to_string_pretty(&self.foils).ok()
    }

    /// Write the cache back to disk (only if it changed since the last save).
    pub fn save(&mut self) {
        if !self.dirty || self.path.is_empty() {
            return;
        }
        if let Some(json) = self.to_json_string() {
            let _ = std::fs::write(&self.path, json);
            self.dirty = false;
        }
    }

    pub fn get(&self, key: &str) -> Option<&StoredPolar> {
        self.foils.get(key)
    }

    /// Insert a freshly simulated polar, persisting it to disk
    /// immediately: every completed calculation is durable on its own, so
    /// an interrupted run loses at most the sweep in flight (the write of
    /// the whole JSON file takes milliseconds next to the seconds each
    /// rust-foil sweep takes).  Stores without a path — the wasm build —
    /// never persist and the save is a no-op, so the wasm host instead
    /// installs the per-polar hook ([`PolarStore::on_insert`]) and each
    /// new polar is pushed to its cache here, at calculation time.
    pub fn insert(&mut self, key: String, polar: StoredPolar) {
        #[cfg(any(test, target_arch = "wasm32"))]
        if let Some(h) = self.on_insert.as_ref() {
            h(&key, &polar);
        }
        self.foils.insert(key.clone(), polar);
        self.new_keys.push(key);
        self.dirty = true;
        self.save();
    }

    /// Install the per-polar persistence hook (the wasm host installs it;
    /// the crate's tests exercise it).  See [`PolarStore::on_insert`].
    #[cfg(any(test, target_arch = "wasm32"))]
    pub(crate) fn set_on_insert(&mut self, hook: Box<dyn Fn(&str, &StoredPolar) + Send + Sync>) {
        self.on_insert = Some(hook);
    }

    /// Insert pre-existing data (a warm-up load, not a new simulation): the
    /// entry is available for lookups but is not marked for persistence.
    pub fn hydrate(&mut self, key: String, polar: StoredPolar) {
        self.foils.insert(key, polar);
    }

    /// The polars inserted (via [`PolarStore::insert`]) since the last call —
    /// what a host without a filesystem should persist itself.  The store
    /// keeps every entry; only the new-entry tracking is drained.
    pub fn take_new_entries(&mut self) -> HashMap<String, StoredPolar> {
        let keys = std::mem::take(&mut self.new_keys);
        keys.into_iter()
            .filter_map(|k| self.foils.get(&k).cloned().map(|p| (k, p)))
            .collect()
    }

    /// Number of stored polars.
    pub fn len(&self) -> usize {
        self.foils.len()
    }

    pub fn is_empty(&self) -> bool {
        self.foils.is_empty()
    }

    /// The cache file location (for tests).
    pub fn path(&self) -> &str {
        &self.path
    }
}

/// Default cache file name, matching the role of `foil_simulator.db`.
pub fn default_cache_path() -> String {
    "foil_cache.json".to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_through_disk() {
        let dir = std::env::temp_dir().join("proply_rs_cache_test");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("cache.json");
        let _ = std::fs::remove_file(&path); // start clean
        let path_str = path.to_str().unwrap().to_string();

        {
            let mut store = PolarStore::load(&path_str);
            assert!(store.get("k").is_none());
            store.insert(
                "k".into(),
                StoredPolar {
                    alpha: vec![0.0, 0.1],
                    cl: vec![0.1, 0.2],
                    cd: vec![0.01, 0.02],
                },
            );
            store.save();
        }
        let store = PolarStore::load(&path_str);
        let polar = store.get("k").expect("cached polar");
        assert_eq!(polar.alpha, vec![0.0, 0.1]);
        assert_eq!(polar.cl, vec![0.1, 0.2]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn missing_file_is_empty() {
        let store = PolarStore::load("/nonexistent/dir/cache.json");
        assert!(store.get("anything").is_none());
    }

    #[test]
    fn insert_persists_each_polar_immediately() {
        // Every freshly simulated polar is written to disk at once: an
        // interrupted run keeps every completed calculation.
        let dir = std::env::temp_dir().join("proply_rs_cache_checkpoint");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("cache.json");
        let _ = std::fs::remove_file(&path);
        let path_str = path.to_str().unwrap().to_string();

        let mut store = PolarStore::load(&path_str);
        store.insert("k1".into(), sample(0.1));
        assert!(
            std::path::Path::new(&path_str).exists(),
            "the first insert must write to disk immediately"
        );
        assert!(!store.dirty, "insert leaves the store clean");
        store.insert("k2".into(), sample(0.2));
        assert!(!store.dirty, "the second insert persists too");

        // The on-disk copy holds every inserted polar.
        let on_disk = PolarStore::load(&path_str);
        assert_eq!(on_disk.len(), 2);
        assert!(on_disk.get("k1").is_some());
        assert!(on_disk.get("k2").is_some());
        let _ = std::fs::remove_dir_all(&dir);
    }

    fn sample(key_cl: f64) -> StoredPolar {
        StoredPolar {
            alpha: vec![0.0, 0.1],
            cl: vec![key_cl, key_cl + 0.1],
            cd: vec![0.01, 0.02],
        }
    }

    #[test]
    fn json_string_round_trips() {
        let mut store = PolarStore::in_memory();
        store.insert("k".into(), sample(0.1));
        let json = store.to_json_string().expect("serializable");
        let mut hydrated = PolarStore::from_json_str(&json);
        assert_eq!(hydrated.len(), 1);
        assert_eq!(hydrated.get("k").cloned(), Some(sample(0.1)));
        // Hydrated data is pre-existing: not pending persistence.
        assert!(hydrated.take_new_entries().is_empty());
    }

    #[test]
    fn from_bad_json_is_empty() {
        assert!(PolarStore::from_json_str("not json").is_empty());
    }

    #[test]
    fn take_new_entries_returns_only_inserted() {
        let mut store = PolarStore::in_memory();
        store.hydrate("old".into(), sample(0.1));
        store.insert("new".into(), sample(0.2));

        let drained = store.take_new_entries();
        assert_eq!(drained.len(), 1);
        assert_eq!(drained.get("new").cloned(), Some(sample(0.2)));

        // The store keeps everything; the drain is one-shot.
        assert_eq!(store.len(), 2);
        assert!(store.take_new_entries().is_empty());
    }

    #[test]
    fn in_memory_store_does_not_persist() {
        let mut store = PolarStore::in_memory();
        store.insert("k".into(), sample(0.1));
        store.save(); // no path: must not panic or mark anything
        assert_eq!(store.len(), 1);
    }

    #[test]
    fn insert_fires_the_per_polar_persist_hook() {
        // The wasm host's checkpoint: every freshly simulated polar — a
        // good sweep or a degenerate failure marker — must reach the hook
        // exactly once, at insert time; pre-existing hydrated data must
        // not (the host already has it).
        let calls = std::sync::Arc::new(std::sync::Mutex::new(Vec::<(String, usize)>::new()));
        let calls2 = calls.clone();
        let mut store = PolarStore::in_memory();
        store.set_on_insert(Box::new(move |key: &str, p: &StoredPolar| {
            calls2.lock().unwrap().push((key.to_string(), p.alpha.len()));
        }));
        store.insert("good".into(), sample(0.1));
        store.insert(
            "bad".into(),
            StoredPolar {
                alpha: Vec::new(),
                cl: Vec::new(),
                cd: Vec::new(),
            },
        );
        store.hydrate("hydrated".into(), sample(0.3));

        let got = calls.lock().unwrap();
        assert_eq!(got.len(), 2, "hook fired per insert, not per hydrate");
        assert_eq!(got[0].0, "good");
        assert_eq!(got[0].1, 2);
        assert_eq!(got[1].0, "bad");
        assert_eq!(got[1].1, 0, "failure markers reach the hook too");
    }

    #[test]
    fn failed_sweep_marker_round_trips_the_web_persistence_chain() {
        // The browser flow the wasm session follows: a failed rust-foil
        // sweep stores a *marker* (an empty polar, via the same insert as
        // a successful sweep); `PropSession::take_new_json` drains it and
        // the host stores it in IndexedDB; the next page load hydrates it
        // back.  The marker must survive that round trip as an empty
        // polar — the degenerate entry that makes `bucket_fits` take the
        // flat-plate fallback instead of ever re-running the doomed sweep.
        let marker = StoredPolar {
            alpha: Vec::new(),
            cl: Vec::new(),
            cd: Vec::new(),
        };

        // Design session: the failed sweep inserts the marker in memory.
        let mut session = PolarStore::in_memory();
        session.insert("h|bucket|mach".into(), marker.clone());
        assert_eq!(session.len(), 1, "marker cached in the session");

        // take_new_json: the marker is among the freshly simulated polars.
        let drained = session.take_new_entries();
        assert_eq!(drained.len(), 1, "marker reported for persistence");
        assert!(drained["h|bucket|mach"].alpha.is_empty());

        // The host persists the drained map; a fresh page session hydrates
        // it (hydrate_json is the same JSON document take_new_json made).
        let json = serde_json::to_string(&drained).unwrap();
        let mut reload = PolarStore::from_json_str(&json);
        let got = reload.get("h|bucket|mach").expect("marker hydrated");
        assert!(got.alpha.is_empty(), "empty marker survives the round trip");
        // Hydrated data is pre-existing: it is not reported again, so the
        // host never re-persists it on the next design.
        assert!(reload.take_new_entries().is_empty());
    }
}
