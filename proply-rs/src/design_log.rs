// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Optional design trace log (the CLI's `--log <file>`).
//!
//! The design loops print a detailed convergence trace — the
//! operating-point torque match, per-station BEM state, the camber scan
//! and the warnings — to the console.  When the CLI opens a trace file
//! here, every [`dprintln!`] / [`deprintln!`] line is mirrored to it as
//! well, so a run's full convergence history survives the terminal (and a
//! failed design can be inspected after the fact).
//!
//! Nothing is written, and no file is touched, when no log file is open —
//! which is the state the WebAssembly build (no filesystem) and the
//! library tests run in.  The file I/O therefore lives entirely behind
//! [`open`], which only the native CLI ever calls.

use std::fs::File;
use std::io::Write;
use std::sync::Mutex;

/// The open trace file, if any.  A `Mutex` because the lifting-line camber
/// scan runs its design passes on scoped threads; each line is one locked
/// write, so log lines never interleave.
static SINK: Mutex<Option<File>> = Mutex::new(None);

fn lock() -> std::sync::MutexGuard<'static, Option<File>> {
    SINK.lock().unwrap_or_else(|p| p.into_inner())
}

/// Open `path` for the design trace, truncating any previous contents.
/// Returns the io error when the file cannot be created — the CLI treats
/// a requested but unwritable log as fatal.
pub fn open(path: &str) -> std::io::Result<()> {
    let file = File::create(path)?;
    *lock() = Some(file);
    Ok(())
}

/// Whether a trace file is currently open.
pub fn is_open() -> bool {
    lock().is_some()
}

/// Close the trace file (the tests reset the global sink with this).
pub fn close() {
    *lock() = None;
}

/// Mirror one line to the open trace file.  A no-op when no `--log` file
/// is open (the default for the wasm build and the library tests).  Each
/// line is flushed immediately, so a crash still leaves the trace up to
/// that point.
pub fn write(line: &str) {
    if let Some(file) = lock().as_mut() {
        let _ = writeln!(file, "{line}");
        let _ = file.flush();
    }
}

/// `println!` mirrored to the trace file: prints to stdout exactly as
/// `println!` would, and appends the same line to the open `--log` file.
#[macro_export]
macro_rules! dprintln {
    ($($arg:tt)*) => {{
        let line = format!($($arg)*);
        println!("{line}");
        $crate::design_log::write(&line);
    }};
}

/// `eprintln!` mirrored to the trace file.
#[macro_export]
macro_rules! deprintln {
    ($($arg:tt)*) => {{
        let line = format!($($arg)*);
        eprintln!("{line}");
        $crate::design_log::write(&line);
    }};
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    /// The sink is process-global, so the content assertion only requires
    /// this test's own marker (other tests' design prints may land in the
    /// same file while it is open; they are harmless extra lines).
    #[test]
    fn open_mirrors_lines_and_closes() {
        let path = std::env::temp_dir().join(format!("proply-log-test-{}.txt", std::process::id()));
        close();
        open(path.to_str().unwrap()).unwrap();
        assert!(is_open());
        crate::dprintln!("proply design trace marker {}", 42);
        close();
        assert!(!is_open());
        let mut text = String::new();
        std::fs::File::open(&path)
            .unwrap()
            .read_to_string(&mut text)
            .unwrap();
        assert!(
            text.contains("proply design trace marker 42"),
            "log content: {text}"
        );
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn closed_sink_is_a_noop() {
        close();
        assert!(!is_open());
        crate::dprintln!("this goes nowhere");
        write("this goes nowhere");
    }
}
