# rsomics-plink-tdt

Transmission disequilibrium test (TDT) over complete trios from a PLINK1 binary
fileset — a Rust reimplementation of `plink --tdt`.

For each variant, over affected offspring whose two parents are both genotyped
and Mendel-consistent, the minor allele transmitted (`T`) and untransmitted
(`U`) by heterozygous parents is counted:

| column  | meaning                                                  |
|---------|----------------------------------------------------------|
| `CHR`   | chromosome (X/Y/XY/MT mapped to 23/24/25/26)             |
| `SNP`   | variant ID                                               |
| `BP`    | base-pair position                                       |
| `A1`    | minor allele (the counted allele)                        |
| `A2`    | major allele                                             |
| `T`     | A1 transmissions from heterozygous parents               |
| `U`     | A1 non-transmissions from heterozygous parents           |
| `OR`    | `T / U` (`NA` if `U = 0`)                                 |
| `CHISQ` | `(T − U)² / (T + U)`, 1 df (`NA` if `T + U = 0`)          |
| `P`     | upper-tail 1-df chi-square p-value                       |

`A1` is the minor allele by **founder** allele frequency, matching PLINK's
convention; this fixes the `T`/`U` orientation and the `A1`/`A2` labels. Only
affected children (`.fam` phenotype `2`) with both parents present count; the
fully-haploid chromosomes Y and MT are skipped.

## Usage

```sh
# write the .tdt table to stdout
rsomics-plink-tdt path/to/fileset

# write to <out>.tdt instead of stdout (matches plink --out)
rsomics-plink-tdt path/to/fileset --out result

# choose the worker-thread count for the per-variant pass
rsomics-plink-tdt path/to/fileset -t 8
```

`path/to/fileset` is the prefix shared by `fileset.bed`, `fileset.bim`,
`fileset.fam` (no extension), exactly as PLINK's `--bfile`.

## Compatibility

The `.tdt` output is **byte-identical** to PLINK 1.9, including the `%.4g`
numeric formatting (with the switch to scientific notation), the `NA` tokens,
the X/Y/XY/MT chromosome remapping, and the adaptive column widths.
`tests/compat.rs` diffs against PLINK live when the `plink` binary is on `PATH`,
and against checked-in PLINK 1.9 golden output otherwise. The transmission table
was validated against PLINK across all 64 (father, mother, child) genotype
combinations.

The default test is implemented. The parent-of-origin (`poo`), exact (`exact`,
`exact-midp`), and permutation (`perm`, `mperm`) variants — which write a
different report or run a separate procedure — are out of scope for this crate.

## Origin

This crate is an independent Rust reimplementation of `plink --tdt` based on:

- The published method: Spielman, McGinnis & Ewens 1993 (TDT,
  doi:10.1086/229918), Chang et al. 2015 (PLINK 1.9,
  doi:10.1186/s13742-015-0047-8), and Purcell et al. 2007 (PLINK 1,
  doi:10.1086/519795).
- The public PLINK 1.9 family-based association documentation
  (<https://www.cog-genomics.org/plink/1.9/fam_assoc>) and binary-fileset format
  spec (<https://www.cog-genomics.org/plink/1.9/formats>).
- Black-box behaviour testing against the `plink` 1.9 binary.

No source code from the GPL upstream was used as reference during
implementation. Test fixtures are independently generated.

License: MIT OR Apache-2.0.
Upstream credit: [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
(Christopher Chang et al., GPLv3).
