//! Persistent polar cache, replacing the Python sqlite `foil_simulator.db`.
//!
//! Polars are stored keyed by `"<foil hash>|<reynolds>|<mach>"` in a single
//! JSON file in the working directory.  Alpha values are stored in radians,
//! matching what the Python database kept.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// One stored polar sweep: alpha (radians), cl, cd per point.
#[derive(Debug, Clone, Serialize, Deserialize)]
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
}

pub fn cache_key(hash: &str, reynolds: f64, mach: f64) -> String {
    format!("{}|{}|{}", hash, reynolds, mach)
}

impl PolarStore {
    /// Load the cache from `path` (missing file = empty cache).
    pub fn load(path: &str) -> Self {
        let foils = std::fs::read_to_string(path)
            .ok()
            .and_then(|s| serde_json::from_str(&s).ok())
            .unwrap_or_default();
        Self {
            path: path.to_string(),
            foils,
            dirty: false,
        }
    }

    /// Write the cache back to disk (only if it changed since the last save).
    pub fn save(&mut self) {
        if !self.dirty {
            return;
        }
        if let Ok(json) = serde_json::to_string_pretty(&self.foils) {
            let _ = std::fs::write(&self.path, json);
            self.dirty = false;
        }
    }

    pub fn get(&self, key: &str) -> Option<&StoredPolar> {
        self.foils.get(key)
    }

    pub fn insert(&mut self, key: String, polar: StoredPolar) {
        self.foils.insert(key, polar);
        self.dirty = true;
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
}
