//! PLINK1 `--tdt`: transmission disequilibrium test over complete trios.
//!
//! For each variant, over affected offspring whose two parents are both
//! genotyped and Mendel-consistent, count the minor allele transmitted (T) and
//! untransmitted (U) by heterozygous parents. The statistic is McNemar's:
//!   CHISQ = (T-U)² / (T+U),  1 df,  OR = T/U.
//!
//! A1 is the minor allele by founder allele frequency (PLINK's convention),
//! which fixes the T/U orientation and the A1/A2 labels.

use rsomics_pgen::Pgen;
use std::io::{self, Write};

pub struct TdtRecord {
    pub chrom: String,
    pub snp: String,
    pub bp: u64,
    pub a1: String,
    pub a2: String,
    pub t: u32,
    pub u: u32,
}

/// Flat 64-entry table keyed by the 6-bit `(dad<<4)|(mom<<2)|child` genotype
/// triple. Each entry packs the transmitted (high nibble) and untransmitted
/// (low nibble) bim-A1 count; both are at most 2. Codes are PLINK 2-bit:
/// 0=HomA1, 1=Missing, 2=Het, 3=HomA2 — missing or Mendel-inconsistent trios
/// score zero.
fn build_transmit_lut() -> [u8; 64] {
    let mut lut = [0u8; 64];
    for dad in 0..4u8 {
        for mom in 0..4u8 {
            for child in 0..4u8 {
                let (t, u) = transmit(dad, mom, child);
                lut[((dad << 4) | (mom << 2) | child) as usize] = (t << 4) | u;
            }
        }
    }
    lut
}

/// A genotype code is a present (non-missing) diploid call.
fn called(code: u8) -> bool {
    matches!(code, 0 | 2 | 3)
}

fn transmit(dad: u8, mom: u8, child: u8) -> (u8, u8) {
    if !called(dad) || !called(mom) || !called(child) {
        return (0, 0);
    }
    let child_a1 = match child {
        0 => 2u8,
        2 => 1,
        _ => 0,
    };
    // One ordered (dad allele, mom allele) pair must reconstruct the child;
    // het parents then transmit a determinable allele. A1 carried as a bit.
    let dad_alleles = parent_alleles(dad);
    let mom_alleles = parent_alleles(mom);
    let mut sol = None;
    for &da in dad_alleles {
        for &ma in mom_alleles {
            if da + ma == child_a1 {
                sol = Some((da, ma));
            }
        }
    }
    let Some((da, ma)) = sol else { return (0, 0) };
    let (mut t, mut u) = (0u8, 0u8);
    if dad == 2 {
        if da == 1 { t += 1 } else { u += 1 }
    }
    if mom == 2 {
        if ma == 1 { t += 1 } else { u += 1 }
    }
    (t, u)
}

/// Distinct A1-counts a parent genotype can transmit (1 = A1, 0 = A2).
fn parent_alleles(code: u8) -> &'static [u8] {
    match code {
        0 => &[1],
        2 => &[1, 0],
        3 => &[0],
        _ => &[],
    }
}

struct Trio {
    dad: usize,
    mom: usize,
    child: usize,
}

/// Affected offspring whose parents are both present in the same family.
fn trios(pgen: &Pgen) -> Vec<Trio> {
    use std::collections::HashMap;
    let mut by_key: HashMap<(&str, &str), usize> = HashMap::new();
    for (i, s) in pgen.samples.iter().enumerate() {
        by_key.insert((s.fid.as_str(), s.iid.as_str()), i);
    }
    pgen.samples
        .iter()
        .enumerate()
        .filter(|(_, s)| s.phen == "2" && s.pid != "0" && s.mid != "0")
        .filter_map(|(child, s)| {
            let dad = *by_key.get(&(s.fid.as_str(), s.pid.as_str()))?;
            let mom = *by_key.get(&(s.fid.as_str(), s.mid.as_str()))?;
            Some(Trio { dad, mom, child })
        })
        .collect()
}

fn founder_mask(pgen: &Pgen) -> Vec<bool> {
    pgen.samples
        .iter()
        .map(|s| s.pid == "0" && s.mid == "0")
        .collect()
}

#[inline]
fn code_at(row: &[u8], s: usize) -> u8 {
    (row[s / 4] >> ((s % 4) * 2)) & 0b11
}

/// Signed A1-minus-A2 founder dosage contribution of one genotype code.
#[inline]
fn dosage_diff(code: u8) -> i32 {
    match code {
        0 => 2,
        3 => -2,
        _ => 0,
    }
}

/// PLINK's `--tdt` skips the fully-haploid chromosomes Y and MT; autosomes,
/// the unplaced chromosome 0, X (23) and XY (25) are all tested.
fn tdt_tested(chrom: &str) -> bool {
    !matches!(report_chrom(chrom).as_str(), "24" | "26")
}

/// A trio with flags marking whether its parents' founder dosage should be
/// counted here. Each founder is attributed to the first trio that names it, so
/// summing over trios counts every founder exactly once.
struct TrioWork {
    dad: usize,
    mom: usize,
    child: usize,
    count_dad: bool,
    count_mom: bool,
}

#[must_use]
pub fn tdt(pgen: &Pgen) -> Vec<TdtRecord> {
    use rayon::prelude::*;
    let lut = build_transmit_lut();
    let triples = trios(pgen);

    // The minor allele is decided over founders. A trio's parents are founders,
    // so their dosage is folded into the trio loop — but a founder may parent
    // several trios, so attribute it to its first trio to count it once.
    let founders = founder_mask(pgen);
    let mut seen = vec![false; pgen.n_samples()];
    let work: Vec<TrioWork> = triples
        .iter()
        .map(|t| {
            let count_dad = !std::mem::replace(&mut seen[t.dad], true);
            let count_mom = !std::mem::replace(&mut seen[t.mom], true);
            TrioWork {
                dad: t.dad,
                mom: t.mom,
                child: t.child,
                count_dad,
                count_mom,
            }
        })
        .collect();
    let other_founders: Vec<u32> = (0..pgen.n_samples())
        .filter(|&s| founders[s] && !seen[s])
        .map(|s| s as u32)
        .collect();

    let bpv = pgen.bytes_per_variant();
    let gt = &pgen.gt_raw;

    (0..pgen.n_variants())
        .into_par_iter()
        .filter(|&v| tdt_tested(&pgen.variants[v].chrom))
        .map(|v| {
            let row = &gt[v * bpv..v * bpv + bpv];
            let (mut t_a1, mut u_a1) = (0u32, 0u32);
            let mut diff = 0i32;
            for w in &work {
                let dad = code_at(row, w.dad);
                let mom = code_at(row, w.mom);
                let key = (dad << 4) | (mom << 2) | code_at(row, w.child);
                let packed = lut[key as usize];
                t_a1 += u32::from(packed >> 4);
                u_a1 += u32::from(packed & 0x0f);
                if w.count_dad {
                    diff += dosage_diff(dad);
                }
                if w.count_mom {
                    diff += dosage_diff(mom);
                }
            }
            for &s in &other_founders {
                diff += dosage_diff(code_at(row, s as usize));
            }
            let var = &pgen.variants[v];
            let (a1, a2, t, u) = if diff <= 0 {
                (&var.a1, &var.a2, t_a1, u_a1)
            } else {
                (&var.a2, &var.a1, u_a1, t_a1)
            };
            TdtRecord {
                chrom: var.chrom.clone(),
                snp: var.id.clone(),
                bp: var.pos,
                a1: a1.clone(),
                a2: a2.clone(),
                t,
                u,
            }
        })
        .collect()
}

/// PLINK maps the sex chromosomes and MT onto numeric codes in its reports.
fn report_chrom(chrom: &str) -> String {
    match chrom {
        "X" | "x" => "23".to_string(),
        "Y" | "y" => "24".to_string(),
        "XY" | "xy" => "25".to_string(),
        "MT" | "mt" | "M" | "m" => "26".to_string(),
        other => other.to_string(),
    }
}

struct Widths {
    chr: usize,
    snp: usize,
    a1: usize,
    a2: usize,
}

impl Widths {
    fn measure(records: &[TdtRecord]) -> Self {
        let mut chr = 0;
        let mut snp = 0;
        let mut a1 = 0;
        let mut a2 = 0;
        for r in records {
            chr = chr.max(report_chrom(&r.chrom).len());
            snp = snp.max(r.snp.len());
            a1 = a1.max(r.a1.len());
            a2 = a2.max(r.a2.len());
        }
        Self {
            chr: chr.max(2) + 2,
            snp: if snp < 5 { 5 } else { snp + 3 },
            a1: a1.max(2) + 2,
            a2: a2.max(2) + 2,
        }
    }
}

/// Write the records in PLINK's `.tdt` layout (default, non-poo).
pub fn write_tdt<W: Write>(records: &[TdtRecord], out: &mut W) -> io::Result<()> {
    let w = Widths::measure(records);
    writeln!(
        out,
        "{:>cw$}{:>sw$}{:>13}{:>a1$}{:>a2$}{:>7}{:>7}{:>13}{:>13}{:>13} ",
        "CHR",
        "SNP",
        "BP",
        "A1",
        "A2",
        "T",
        "U",
        "OR",
        "CHISQ",
        "P",
        cw = w.chr,
        sw = w.snp,
        a1 = w.a1,
        a2 = w.a2,
    )?;
    for r in records {
        let (or, chisq, p) = stats(r.t, r.u);
        writeln!(
            out,
            "{:>cw$}{:>sw$}{:>13}{:>a1$}{:>a2$}{:>7}{:>7}{:>13}{:>13}{:>13}  ",
            report_chrom(&r.chrom),
            r.snp,
            r.bp,
            r.a1,
            r.a2,
            r.t,
            r.u,
            or,
            chisq,
            p,
            cw = w.chr,
            sw = w.snp,
            a1 = w.a1,
            a2 = w.a2,
        )?;
    }
    Ok(())
}

/// OR / CHISQ / P tokens. T+U=0 → all NA; U=0 → OR NA; T=0 → OR 0.
fn stats(t: u32, u: u32) -> (String, String, String) {
    let n = t + u;
    if n == 0 {
        return ("NA".into(), "NA".into(), "NA".into());
    }
    let or = if u == 0 {
        "NA".to_string()
    } else {
        fmt_g(f64::from(t) / f64::from(u))
    };
    let diff = f64::from(t) - f64::from(u);
    let chisq = diff * diff / f64::from(n);
    let p = chisq_1df_sf(chisq);
    (or, fmt_g(chisq), fmt_g(p))
}

/// Upper-tail probability of a 1-df chi-square = the regularised upper
/// incomplete gamma `Q(1/2, x/2)`, evaluated to full precision.
fn chisq_1df_sf(x: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    gamma_q(0.5, x / 2.0)
}

/// ln Γ(z) — Lanczos approximation.
fn ln_gamma(z: f64) -> f64 {
    const C: [f64; 6] = [
        76.180_091_729_471_46,
        -86.505_320_329_416_77,
        24.014_098_240_830_91,
        -1.231_739_572_450_155,
        0.001_208_650_973_866_179,
        -0.000_005_395_239_384_953,
    ];
    let mut x = z;
    let mut tmp = z + 5.5;
    tmp -= (z + 0.5) * tmp.ln();
    let mut ser = 1.000_000_000_190_015;
    for c in C {
        x += 1.0;
        ser += c / x;
    }
    -tmp + (2.506_628_274_631_000_5 * ser / z).ln()
}

/// Regularised upper incomplete gamma `Q(a, x)` (Numerical Recipes `gser`/`gcf`).
fn gamma_q(a: f64, x: f64) -> f64 {
    if x < a + 1.0 {
        1.0 - gamma_p_series(a, x)
    } else {
        gamma_q_cf(a, x)
    }
}

fn gamma_p_series(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..400 {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * 1e-16 {
            break;
        }
    }
    sum * (-x + a * x.ln() - gln).exp()
}

fn gamma_q_cf(a: f64, x: f64) -> f64 {
    const TINY: f64 = 1e-300;
    let gln = ln_gamma(a);
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / TINY;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..400 {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < TINY {
            d = TINY;
        }
        c = b + an / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-16 {
            break;
        }
    }
    (-x + a * x.ln() - gln).exp() * h
}

/// PLINK's numeric output format: the shortest round-tripping decimal rounded
/// to 4 significant figures with round-half-to-even, then `%g`-displayed
/// (trailing zeros stripped, scientific notation outside the exponent range
/// -4..4). PLINK's dtoa differs from libc `%g` only at exact half-way ties,
/// which it breaks toward the even digit; reproducing it needs the shortest
/// decimal, not the raw binary value.
fn fmt_g(x: f64) -> String {
    const SIG: usize = 4;
    if x.is_nan() {
        return "nan".to_string();
    }
    if x == 0.0 {
        return "0".to_string();
    }
    let neg = x < 0.0;
    // Shortest round-tripping decimal digits (no sign, no point) and the
    // power-of-ten exponent of the leading digit.
    let (digits, lead_exp) = shortest_decimal(x.abs());
    let (digits, exp) = round_sig_half_even(&digits, lead_exp, SIG);

    let mut s = if !(-4..SIG as i32).contains(&exp) {
        let mant = mantissa(&digits, 1);
        format!("{mant}e{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs())
    } else if exp >= 0 {
        mantissa(&digits, (exp + 1) as usize)
    } else {
        let zeros = "0".repeat((-exp - 1) as usize);
        strip_trailing(&format!("0.{zeros}{digits}"))
    };
    if neg {
        s.insert(0, '-');
    }
    s
}

/// Digits of `x`'s shortest round-tripping decimal (Ryū via `{}`), returned as
/// a digit string with no leading zeros plus the base-ten exponent of the
/// first digit. `x` must be finite and positive.
fn shortest_decimal(x: f64) -> (String, i32) {
    let sci = format!("{:e}", x); // e.g. "9.3125e-1"
    let (mant, e) = sci.split_once('e').unwrap();
    let exp: i32 = e.parse().unwrap();
    let digits: String = mant.chars().filter(|c| c.is_ascii_digit()).collect();
    (digits, exp)
}

/// Round a digit string (first digit has power `lead_exp`) to `sig` significant
/// figures, half-to-even. Returns the rounded digit string and the power of its
/// leading digit (which may shift on carry).
fn round_sig_half_even(digits: &str, lead_exp: i32, sig: usize) -> (String, i32) {
    let bytes: Vec<u8> = digits.bytes().map(|b| b - b'0').collect();
    if bytes.len() <= sig {
        let mut d: Vec<u8> = bytes;
        while d.len() > 1 && *d.last().unwrap() == 0 {
            d.pop();
        }
        return (d.iter().map(|&b| (b + b'0') as char).collect(), lead_exp);
    }
    let mut kept: Vec<u8> = bytes[..sig].to_vec();
    let next = bytes[sig];
    let rest_nonzero = bytes[sig + 1..].iter().any(|&b| b != 0);
    let round_up = next > 5 || (next == 5 && (rest_nonzero || kept[sig - 1] % 2 == 1));
    let mut lead = lead_exp;
    if round_up {
        let mut i = sig;
        loop {
            if i == 0 {
                kept.insert(0, 1);
                lead += 1;
                kept.pop();
                break;
            }
            i -= 1;
            if kept[i] == 9 {
                kept[i] = 0;
            } else {
                kept[i] += 1;
                break;
            }
        }
    }
    while kept.len() > 1 && *kept.last().unwrap() == 0 {
        kept.pop();
    }
    (kept.iter().map(|&b| (b + b'0') as char).collect(), lead)
}

/// Format `digits` with `int_len` digits before the decimal point (zero-padding
/// on the right if short), trailing zeros stripped.
fn mantissa(digits: &str, int_len: usize) -> String {
    let padded = if digits.len() < int_len {
        format!("{digits}{}", "0".repeat(int_len - digits.len()))
    } else {
        digits.to_string()
    };
    if padded.len() <= int_len {
        padded
    } else {
        strip_trailing(&format!("{}.{}", &padded[..int_len], &padded[int_len..]))
    }
}

fn strip_trailing(s: &str) -> String {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transmit_table_matches_verified_combos() {
        // Spot-check entries validated byte-exact against plink 1.9 over all 64
        // (dad, mom, child) combos. (t_a1, u_a1) where A1 is the bim-A1 allele.
        assert_eq!(transmit(0, 2, 0), (1, 0)); // mom het, child HomA1
        assert_eq!(transmit(0, 2, 2), (0, 1)); // mom het, child Het
        assert_eq!(transmit(2, 2, 2), (1, 1)); // both het, child Het
        assert_eq!(transmit(2, 2, 0), (2, 0)); // both het, child HomA1
        assert_eq!(transmit(2, 2, 3), (0, 2)); // both het, child HomA2
        assert_eq!(transmit(0, 0, 2), (0, 0)); // Mendel-inconsistent
        assert_eq!(transmit(1, 2, 2), (0, 0)); // missing parent
        assert_eq!(transmit(2, 3, 2), (1, 0));
    }

    #[test]
    fn stats_edge_cases() {
        assert_eq!(stats(0, 0), ("NA".into(), "NA".into(), "NA".into()));
        let (or, chisq, _) = stats(3, 0);
        assert_eq!(or, "NA");
        assert_eq!(chisq, "3");
        let (or, _, _) = stats(0, 2);
        assert_eq!(or, "0");
    }

    #[test]
    fn g_formatting_matches_plink() {
        assert_eq!(fmt_g(1.273), "1.273");
        assert_eq!(fmt_g(0.8868), "0.8868");
        assert_eq!(fmt_g(0.1573), "0.1573");
        assert_eq!(fmt_g(0.0), "0");
        assert_eq!(fmt_g(1.44), "1.44");
        assert_eq!(fmt_g(3.0), "3");
        assert_eq!(fmt_g(0.0006871), "0.0006871");
    }

    #[test]
    fn g_half_ties_round_to_even() {
        // PLINK breaks exact decimal half-ties toward the even digit.
        assert_eq!(fmt_g(894.0 / 960.0), "0.9312");
        assert_eq!(fmt_g(473.0 / 400.0), "1.182");
        assert_eq!(fmt_g(696.0 / 640.0), "1.088");
        assert_eq!(fmt_g(763.0 / 800.0), "0.9538");
        assert_eq!(fmt_g(431.0 / 400.0), "1.078");
        assert_eq!(fmt_g(0.91875), "0.9188");
    }

    #[test]
    fn p_value_matches_plink() {
        // chi2.sf(x, 1) values plink prints.
        assert_eq!(fmt_g(chisq_1df_sf(2.0)), "0.1573");
        assert_eq!(fmt_g(chisq_1df_sf(5.76)), "0.0164");
        assert_eq!(fmt_g(chisq_1df_sf(0.0)), "1");
        assert_eq!(fmt_g(chisq_1df_sf(3.0)), "0.08326");
    }
}
