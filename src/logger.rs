//! Global stdout log collector for generating PDF summary pages.

use std::sync::{Mutex, OnceLock};

static LOG_BUFFER: OnceLock<Mutex<Vec<String>>> = OnceLock::new();

fn get_buffer() -> &'static Mutex<Vec<String>> {
    LOG_BUFFER.get_or_init(|| Mutex::new(Vec::new()))
}

/// Records a line to the global stdout log buffer.
pub fn record_line(line: &str) {
    if let Ok(mut buf) = get_buffer().lock() {
        buf.push(line.to_string());
    }
}

/// Clears the global stdout log buffer (e.g., at the start of a run).
pub fn clear_log() {
    if let Ok(mut buf) = get_buffer().lock() {
        buf.clear();
    }
}

/// Retrieves all captured stdout log lines joined by newlines.
pub fn get_log() -> String {
    if let Ok(buf) = get_buffer().lock() {
        buf.join("\n")
    } else {
        String::new()
    }
}

/// Macro to print to stdout and record in the stdout log buffer.
#[macro_export]
macro_rules! log_println {
    () => {{
        println!();
        $crate::logger::record_line("");
    }};
    ($($arg:tt)*) => {{
        let msg = format!($($arg)*);
        println!("{}", msg);
        $crate::logger::record_line(&msg);
    }};
}
