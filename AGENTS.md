# Development notes

- The two crates are developed independently.  When working on one, build
  and test only that crate (`cargo test -p proply-rs`, `cargo clippy -p
  rust-foil`, ...): do not run the other crate's test suite.  Workspace
  flags (`--workspace`, `--all`) drag in the other crate's tests.
  `proply-rs` depends on `rust-foil`, so compiling the dependency is
  unavoidable; running its tests is not. 
- The webassembly version is deployed using vercel in proply-rs/web. 
- The webassembly needs to be built before committing to github (see make wasm) in Makefile.
- The web demo shows a build label ("yyyy-mm-dd.xx", the last commit's
  date plus its per-day build number) in proply-rs/web/build.js:
  regenerate it with `make build-date` before committing web changes,
  so the deployed label matches the deployed sources.
