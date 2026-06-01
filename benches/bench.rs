use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_pgen::Pgen;
use rsomics_plink_tdt::tdt;
use std::hint::black_box;
use std::path::PathBuf;

fn bench_tdt(c: &mut Criterion) {
    // Set RSOMICS_TDT_BFILE to a representative-large trio fileset prefix to
    // bench the hot path; otherwise the in-repo golden is used.
    let prefix = std::env::var("RSOMICS_TDT_BFILE")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/trio"));
    let pgen = Pgen::load(&prefix).expect("load fileset");

    c.bench_function("tdt", |b| {
        b.iter(|| black_box(tdt(black_box(&pgen))));
    });
}

criterion_group!(benches, bench_tdt);
criterion_main!(benches);
