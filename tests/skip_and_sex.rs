//! Behaviour PLINK 1.9 has that a naive port misses: it accepts the `-9`
//! (unknown-sex) code in the `.fam` sex column, and it writes *no* `.tdt` at
//! all — emitting a `Warning: Skipping --tdt …` and exiting 0 — when nothing is
//! testable. These lock both against baked fixtures (no live PLINK needed).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-plink-tdt"))
}

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

/// A minimal variant-major `.bed` (all-HomA1) sized for `n_samples` × 1 variant.
fn bed_one_variant() -> [u8; 4] {
    [0x6c, 0x1b, 0x01, 0x00]
}

fn write_fileset(dir: &Path, name: &str, fam: &str, bim: &str, bed: &[u8]) -> PathBuf {
    let prefix = dir.join(name);
    fs::write(prefix.with_extension("fam"), fam).unwrap();
    fs::write(prefix.with_extension("bim"), bim).unwrap();
    fs::write(prefix.with_extension("bed"), bed).unwrap();
    prefix
}

struct Run {
    ok: bool,
    stdout: String,
    stderr: String,
}

fn run(prefix: &Path) -> Run {
    let out = Command::new(ours()).arg(prefix).output().unwrap();
    Run {
        ok: out.status.success(),
        stdout: String::from_utf8(out.stdout).unwrap(),
        stderr: String::from_utf8(out.stderr).unwrap(),
    }
}

/// `-9` sex on the affected offspring (PLINK only rejects it on *parents*) must
/// load and yield the same report as the golden trio, since `--tdt` ignores sex.
#[test]
fn neg9_sex_on_offspring_matches_golden() {
    let dir = tempfile::tempdir().unwrap();
    let fam = "FAM1 DAD1 0 0 1 1\n\
               FAM1 MOM1 0 0 2 1\n\
               FAM1 KID1 DAD1 MOM1 -9 2\n\
               FAM2 DAD2 0 0 1 1\n\
               FAM2 MOM2 0 0 2 1\n\
               FAM2 KID2 DAD2 MOM2 -9 2\n";
    let bim = fs::read_to_string(golden_dir().join("trio.bim")).unwrap();
    let bed = fs::read(golden_dir().join("trio.bed")).unwrap();
    let prefix = write_fileset(dir.path(), "neg9", fam, &bim, &bed);

    let r = run(&prefix);
    assert!(r.ok, "stderr: {}", r.stderr);
    let golden = fs::read_to_string(golden_dir().join("trio.tdt.golden")).unwrap();
    assert_eq!(r.stdout, golden);
}

/// No offspring with both parents present → `there are no trios`, no output.
#[test]
fn skip_no_trios() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = write_fileset(
        dir.path(),
        "founders",
        "F1 A 0 0 1 2\nF1 B 0 0 2 1\n",
        "1 rs1 0 1000 A G\n",
        &bed_one_variant(),
    );
    let r = run(&prefix);
    assert!(r.ok);
    assert!(r.stdout.is_empty(), "expected no report, got: {}", r.stdout);
    assert_eq!(
        r.stderr,
        "Warning: Skipping --tdt since there are no trios.\n"
    );
}

/// A complete trio but no affected child → `no trios with an affected child`.
#[test]
fn skip_no_affected_trio() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = write_fileset(
        dir.path(),
        "unaff",
        "F1 DAD 0 0 1 1\nF1 MOM 0 0 2 1\nF1 KID DAD MOM 1 1\n",
        "1 rs1 0 1000 A G\n",
        &bed_one_variant(),
    );
    let r = run(&prefix);
    assert!(r.ok);
    assert!(r.stdout.is_empty());
    assert_eq!(
        r.stderr,
        "Warning: Skipping --tdt since there are no trios with an affected child, and no\n\
         discordant parent pairs.\n"
    );
}

/// Affected trio present but every variant is on Y → `no autosomal or Xchr data`.
#[test]
fn skip_no_autosomal_or_xchr() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = write_fileset(
        dir.path(),
        "yonly",
        "F1 DAD 0 0 1 1\nF1 MOM 0 0 2 1\nF1 KID DAD MOM 1 2\n",
        "24 rsY 0 1000 A G\n",
        &bed_one_variant(),
    );
    let r = run(&prefix);
    assert!(r.ok);
    assert!(r.stdout.is_empty());
    assert_eq!(
        r.stderr,
        "Warning: Skipping --tdt since there is no autosomal or Xchr data.\n"
    );
}

/// Skipping writes no `<out>.tdt` file (matching PLINK, which leaves none).
#[test]
fn skip_writes_no_output_file() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = write_fileset(
        dir.path(),
        "founders",
        "F1 A 0 0 1 2\nF1 B 0 0 2 1\n",
        "1 rs1 0 1000 A G\n",
        &bed_one_variant(),
    );
    let out = dir.path().join("report");
    let status = Command::new(ours())
        .arg(&prefix)
        .arg("--out")
        .arg(&out)
        .output()
        .unwrap();
    assert!(status.status.success());
    assert!(!out.with_extension("tdt").exists());
}
