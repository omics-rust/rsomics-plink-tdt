use clap::Parser;
use rsomics_plink_tdt::{TdtOutput, load_fileset, tdt_report, write_tdt};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "rsomics-plink-tdt",
    about = "PLINK1 transmission disequilibrium test for trios (plink --tdt)",
    version
)]
struct Cli {
    /// Path prefix for the .bed/.bim/.fam fileset (without extension).
    bfile: PathBuf,

    /// Write the report to <OUT>.tdt instead of stdout (plink --out).
    #[arg(long)]
    out: Option<PathBuf>,

    /// Worker threads for the per-variant counting pass.
    #[arg(short = 't', long, default_value_t = num_cpus())]
    threads: usize,
}

fn num_cpus() -> usize {
    std::thread::available_parallelism().map_or(1, |n| n.get())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: Cli) -> anyhow::Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads.max(1))
        .build_global()
        .ok();

    let pgen = load_fileset(&cli.bfile)?;
    let records = match tdt_report(&pgen) {
        TdtOutput::Report(records) => records,
        // PLINK writes no report and exits 0 when nothing is testable.
        TdtOutput::Skip(reason) => {
            eprintln!("{}", reason.warning());
            return Ok(());
        }
    };

    match cli.out {
        Some(prefix) => {
            let mut w =
                BufWriter::with_capacity(1 << 20, File::create(prefix.with_extension("tdt"))?);
            write_tdt(&records, &mut w)?;
            w.flush()?;
        }
        None => {
            let stdout = io::stdout();
            let mut w = BufWriter::with_capacity(1 << 20, stdout.lock());
            write_tdt(&records, &mut w)?;
            w.flush()?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_definition_is_valid() {
        <Cli as clap::CommandFactory>::command().debug_assert();
    }
}
