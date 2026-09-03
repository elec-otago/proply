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
    /// Inserts appended to the journal since the last whole-file rewrite
    /// ([`PolarStore::save`]): a path-backed store's per-polar durability
    /// is the journal append, so the whole-file rewrite only happens every
    /// [`FULL_SAVE_EVERY`] inserts (or at save/exit).
    pending: usize,
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

/// How many journal-appended inserts accumulate before the whole cache
/// file is rewritten ([`PolarStore::save`] compacts: full pretty JSON to a
/// temp file, atomic rename, journal dropped).  The journal makes each
/// individual polar durable with an O(one polar) append; rewriting the
/// whole (growing) file on *every* insert made the per-polar checkpoint
/// cost O(the whole cache) — the design-log bottleneck once the file
/// reaches tens of megabytes.
const FULL_SAVE_EVERY: usize = 200;

/// One journal record: a freshly inserted polar (`foil_cache.json.journal`,
/// one NDJSON line per record).
#[derive(Serialize, Deserialize)]
struct JournalEntry {
    key: String,
    alpha: Vec<f64>,
    cl: Vec<f64>,
    cd: Vec<f64>,
}

pub fn cache_key(hash: &str, reynolds: f64, mach: f64) -> String {
    format!("{}|{}|{}", hash, reynolds, mach)
}

/// The append-only journal's path: the cache path with `.journal` appended.
fn journal_path(cache_path: &str) -> String {
    format!("{cache_path}.journal")
}

impl PolarStore {
    /// An empty store that never persists anywhere (in-memory sessions,
    /// e.g. the WebAssembly build).
    pub fn in_memory() -> Self {
        Self::default()
    }

    /// Load the cache from `path` (missing file = empty cache): the whole
    /// file first, then every journal record appended since its last
    /// rewrite — journal records win over the file's copy of a key (they
    /// are newer).  A torn journal tail (a killed append) is skipped.
    pub fn load(path: &str) -> Self {
        let mut foils = std::fs::read_to_string(path)
            .ok()
            .and_then(|s| Self::parse(&s))
            .unwrap_or_default();
        if let Ok(text) = std::fs::read_to_string(journal_path(path)) {
            for line in text.lines() {
                if let Ok(rec) = serde_json::from_str::<JournalEntry>(line) {
                    foils.insert(
                        rec.key,
                        StoredPolar {
                            alpha: rec.alpha,
                            cl: rec.cl,
                            cd: rec.cd,
                        },
                    );
                }
            }
        }
        Self {
            path: path.to_string(),
            foils,
            dirty: false,
            new_keys: Vec::new(),
            pending: 0,
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
            pending: 0,
            #[cfg(any(test, target_arch = "wasm32"))]
            on_insert: None,
        }
    }

    /// The cache contents as a JSON document (the on-disk format).  `None`
    /// only if serialization fails.
    pub fn to_json_string(&self) -> Option<String> {
        serde_json::to_string_pretty(&self.foils).ok()
    }

    /// Compact the cache to the whole file (only if it changed since the
    /// last save): the full pretty JSON is written to a temp file and
    /// renamed over the cache (atomic on POSIX — a killed save leaves the
    /// previous file intact and the journal to replay), then the journal
    /// is dropped.  Called by [`PolarStore::insert`] every
    /// [`FULL_SAVE_EVERY`] inserts and at exit (the CLI's final save).
    pub fn save(&mut self) {
        if !self.dirty || self.path.is_empty() {
            return;
        }
        if let Some(json) = self.to_json_string() {
            let tmp = format!("{}.tmp", self.path);
            let ok = std::fs::write(&tmp, json).is_ok()
                && std::fs::rename(&tmp, &self.path).is_ok();
            if ok {
                // Every journaled record is now in the file: drop the
                // journal so a later load does not replay stale records
                // (harmless if this removal fails — replay of identical
                // records is idempotent).
                let _ = std::fs::remove_file(journal_path(&self.path));
                self.dirty = false;
                self.pending = 0;
            } else {
                let _ = std::fs::remove_file(tmp);
            }
        }
    }

    /// Append one freshly inserted polar to the journal: the per-polar
    /// durability step for path-backed stores (a single small line, not a
    /// rewrite of the whole cache).
    fn append_journal(&mut self, key: &str, polar: &StoredPolar) {
        let line = serde_json::to_string(&JournalEntry {
            key: key.to_string(),
            alpha: polar.alpha.clone(),
            cl: polar.cl.clone(),
            cd: polar.cd.clone(),
        });
        let Ok(line) = line else { return };
        if let Ok(mut f) = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(journal_path(&self.path))
        {
            use std::io::Write;
            let _ = writeln!(f, "{line}");
        }
    }

    pub fn get(&self, key: &str) -> Option<&StoredPolar> {
        self.foils.get(key)
    }

    /// Insert a freshly simulated polar, persisting it to disk at once:
    /// path-backed stores append the record to the journal immediately —
    /// every completed calculation is durable on its own, at O(one polar)
    /// — and rewrite the whole file only every [`FULL_SAVE_EVERY`] inserts
    /// (an interrupted run loses at most the sweep in flight; a journal
    /// line is milliseconds next to the seconds each rust-foil sweep
    /// takes).  Stores without a path — the wasm build — never journal and
    /// the save is a no-op, so the wasm host instead installs the
    /// per-polar hook ([`PolarStore::on_insert`]) and each new polar is
    /// pushed to its cache here, at calculation time.
    pub fn insert(&mut self, key: String, polar: StoredPolar) {
        #[cfg(any(test, target_arch = "wasm32"))]
        if let Some(h) = self.on_insert.as_ref() {
            h(&key, &polar);
        }
        if !self.path.is_empty() {
            self.append_journal(&key, &polar);
            self.pending += 1;
        }
        self.foils.insert(key.clone(), polar);
        self.new_keys.push(key);
        self.dirty = true;
        if self.pending >= FULL_SAVE_EVERY {
            self.save();
        }
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
    fn each_insert_is_journaled_immediately_and_reloads() {
        // Path-backed stores append every freshly simulated polar to the
        // journal at once (the per-polar durability step): a reload before
        // any whole-file rewrite must still see every insert.
        let dir = std::env::temp_dir().join("proply_rs_cache_journal");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("cache.json");
        let journal = dir.join("cache.json.journal");
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&journal);
        let path_str = path.to_str().unwrap().to_string();

        {
            let mut store = PolarStore::load(&path_str);
            store.insert("k1".into(), sample(0.1));
            assert!(
                journal.exists(),
                "the first insert must be journaled immediately"
            );
            assert!(store.dirty, "the whole file waits for the next compaction");
            store.insert("k2".into(), sample(0.2));
        }
        // The main file was never rewritten, yet a fresh load replays the
        // journal and sees both polars.
        assert_eq!(std::fs::read_to_string(&path).unwrap_or_default(), "");
        let reloaded = PolarStore::load(&path_str);
        assert_eq!(reloaded.len(), 2);
        assert_eq!(reloaded.get("k1").cloned(), Some(sample(0.1)));
        assert_eq!(reloaded.get("k2").cloned(), Some(sample(0.2)));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn journal_replay_tolerates_a_torn_tail_and_duplicates() {
        // A killed append can leave a partial final line (skipped), and a
        // compaction that crashed between the rename and the journal drop
        // leaves records already in the file (idempotent re-apply).
        let dir = std::env::temp_dir().join("proply_rs_cache_journal_torn");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("cache.json");
        let journal = dir.join("cache.json.journal");
        let path_str = path.to_str().unwrap().to_string();

        // Seed the main file with k1 and journal it again (duplicate) plus
        // a good k2 and a torn tail.
        std::fs::write(&path, r#"{"k1": {"alpha": [0.0], "cl": [0.1], "cd": [0.01]}}"#)
            .unwrap();
        std::fs::write(
            &journal,
            r#"{"key":"k1","alpha":[0.0],"cl":[0.1],"cd":[0.01]}
{"key":"k2","alpha":[0.2],"cl":[0.3],"cd":[0.02]}
{"key":"k3","alpha":[0.4],"cl":[0.5],"cd":[0.03]
"#,
        )
        .unwrap();

        let store = PolarStore::load(&path_str);
        assert_eq!(store.len(), 2, "file + journal; the torn tail line is skipped");
        assert_eq!(store.get("k1").unwrap().cl, vec![0.1], "journal re-apply is idempotent");
        assert_eq!(store.get("k2").unwrap().cl, vec![0.3]);
        assert!(store.get("k3").is_none(), "the torn line must not parse");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn compaction_rewrites_the_whole_file_and_drops_the_journal() {
        // After FULL_SAVE_EVERY inserts the store compacts: the main file
        // holds everything and the journal is gone, so later loads do not
        // replay stale records.
        let dir = std::env::temp_dir().join("proply_rs_cache_compact");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("cache.json");
        let journal = dir.join("cache.json.journal");
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&journal);
        let path_str = path.to_str().unwrap().to_string();

        let mut store = PolarStore::load(&path_str);
        for i in 0..FULL_SAVE_EVERY {
            store.insert(format!("k{i}"), sample(i as f64 * 0.01));
        }
        assert!(!store.dirty, "the threshold insert compacts");
        assert!(!journal.exists(), "the journal is dropped on compaction");

        let reloaded = PolarStore::load(&path_str);
        assert_eq!(reloaded.len(), FULL_SAVE_EVERY);
        assert_eq!(
            reloaded
                .get(&format!("k{}", FULL_SAVE_EVERY - 1))
                .cloned(),
            Some(sample((FULL_SAVE_EVERY as f64 - 1.0) * 0.01))
        );
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
