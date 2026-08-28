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
#[derive(Debug, Default)]
pub struct PolarStore {
    path: String,
    foils: HashMap<String, StoredPolar>,
    dirty: bool,
    /// Keys inserted since the last [`PolarStore::take_new_entries`] drain —
    /// freshly simulated polars a host without a filesystem (the browser)
    /// persists itself.
    new_keys: Vec<String>,
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

    /// Insert a freshly simulated polar (marked for persistence).
    pub fn insert(&mut self, key: String, polar: StoredPolar) {
        self.foils.insert(key.clone(), polar);
        self.new_keys.push(key);
        self.dirty = true;
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
}
