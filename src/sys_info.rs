//! Utility module to retrieve system hardware and operating system info.

use num_format::{Locale, ToFormattedString};
use std::fs;

/// Returns the CPU model string (e.g. "AMD Ryzen 7 5700G with Radeon Graphics").
pub fn get_cpu_model() -> String {
    if let Ok(cpuinfo) = fs::read_to_string("/proc/cpuinfo") {
        for line in cpuinfo.lines() {
            if line.starts_with("model name") {
                if let Some(pos) = line.find(':') {
                    return line[pos + 1..].trim().to_string();
                }
            }
        }
    }
    format!("{} ({})", std::env::consts::ARCH, std::env::consts::OS)
}

/// Returns the total physical memory as a formatted string (e.g. "58.7 GB (60,161 MB)").
pub fn get_total_memory() -> String {
    if let Ok(meminfo) = fs::read_to_string("/proc/meminfo") {
        for line in meminfo.lines() {
            if line.starts_with("MemTotal:") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let Ok(kb) = parts[1].parse::<u64>() {
                        let mb = kb / 1024;
                        let gb = (mb as f64) / 1024.0;
                        return format!("{:.1} GB ({} MB)", gb, mb.to_formatted_string(&Locale::en));
                    }
                }
            }
        }
    }
    "Unknown".to_string()
}

/// Returns the operating system description (e.g. "Kali GNU/Linux Rolling (x86_64)").
pub fn get_operating_system() -> String {
    let mut name = String::new();
    if let Ok(os_release) = fs::read_to_string("/etc/os-release") {
        for line in os_release.lines() {
            if line.starts_with("PRETTY_NAME=") {
                let val = line["PRETTY_NAME=".len()..].trim_matches('"');
                name = val.to_string();
                break;
            }
        }
    }
    if name.is_empty() {
        name = std::env::consts::OS.to_string();
    }
    format!("{} ({})", name, std::env::consts::ARCH)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sys_info_retrieval() {
        let cpu = get_cpu_model();
        let mem = get_total_memory();
        let os = get_operating_system();

        assert!(!cpu.is_empty());
        assert!(!mem.is_empty());
        assert!(!os.is_empty());
    }
}
