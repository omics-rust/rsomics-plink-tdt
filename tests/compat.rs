//! Differential compatibility tests against PLINK 1.9 `--tdt`.
//!
//! Our binary's `.tdt` output is compared field-by-field to PLINK's: CHR/SNP/
//! A1/A2 strings exact, BP/T/U integers exact, OR/CHISQ/P numeric tokens exact
//! (PLINK prints `%.4g`, reproduced byte-for-byte). With `plink` on PATH we diff
//! live; otherwise we diff against checked-in PLINK 1.9 golden output.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-plink-tdt"))
}

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn plink_available() -> bool {
    Command::new("plink")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn fields(text: &str) -> Vec<Vec<String>> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(str::to_string).collect())
        .collect()
}

fn run_ours(prefix: &Path) -> String {
    let out = Command::new(ours())
        .arg(prefix)
        .output()
        .expect("run rsomics-plink-tdt");
    assert!(
        out.status.success(),
        "rsomics-plink-tdt failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf8")
}

fn assert_fields_equal(ours: &str, ref_text: &str) {
    let a = fields(ours);
    let b = fields(ref_text);
    assert_eq!(a.len(), b.len(), "row count differs");
    for (i, (x, y)) in a.iter().zip(&b).enumerate() {
        assert_eq!(x, y, "row {i} differs:\n ours: {x:?}\n ref:  {y:?}");
    }
}

fn matches_golden(prefix: &str) {
    let ours = run_ours(&golden_dir().join(prefix));
    let golden = std::fs::read_to_string(golden_dir().join(format!("{prefix}.tdt.golden")))
        .expect("read golden");
    // Byte-exact is the contract; assert it directly.
    assert_eq!(ours, golden, "{prefix}: output differs from PLINK golden");
}

#[test]
fn trio_matches_golden() {
    matches_golden("trio");
}

#[test]
fn withx_matches_golden() {
    matches_golden("withx");
}

#[test]
fn many_matches_golden() {
    matches_golden("many");
}

#[test]
fn header_is_plink_shape() {
    let ours = run_ours(&golden_dir().join("trio"));
    let header: Vec<&str> = ours.lines().next().unwrap().split_whitespace().collect();
    assert_eq!(
        header,
        ["CHR", "SNP", "BP", "A1", "A2", "T", "U", "OR", "CHISQ", "P"]
    );
}

/// Live differential against the upstream PLINK binary, when present.
#[test]
fn matches_live_plink() {
    if !plink_available() {
        eprintln!("plink not on PATH; skipping live differential");
        return;
    }
    for prefix in ["trio", "withx", "many"] {
        let tmp = tempfile::Builder::new()
            .prefix("plink-tdt-compat-")
            .tempdir_in(std::env::var("TMPDIR").unwrap_or_else(|_| "/tmp".into()))
            .expect("tempdir");
        let out_prefix = tmp.path().join("ref");
        let status = Command::new("plink")
            .args([
                "--bfile",
                golden_dir().join(prefix).to_str().unwrap(),
                "--tdt",
                "--allow-no-sex",
                "--out",
                out_prefix.to_str().unwrap(),
            ])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .expect("run plink");
        assert!(status.success(), "plink --tdt failed for {prefix}");
        let ref_text =
            std::fs::read_to_string(out_prefix.with_extension("tdt")).expect("read .tdt");
        assert_fields_equal(&run_ours(&golden_dir().join(prefix)), &ref_text);
    }
}
