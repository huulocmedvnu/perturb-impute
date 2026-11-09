//! preprocess pb <cache_dir> <out_prefix>
//!     mean CP10k per perturbation for all cells and two interleaved halves, from a cell cache
//!     written by `dosecurve cache <in.h5ad> <cache_dir>` (https://github.com/huulocmedvnu/dosecurve)
//! preprocess z --pb <out_prefix> --genes <genes.txt> --out <z_prefix>
//!     z-scores against non-targeting control cells, read by `panelimp` (<z_prefix>.z.f32, .meta.tsv, .genes.tsv)
mod cache;
mod pb;
mod zstat;
use anyhow::{bail, Result};

fn arg(a: &[String], k: &str) -> Result<String> {
    match a.iter().position(|x| x == k) { Some(i) if i + 1 < a.len() => Ok(a[i + 1].clone()), _ => bail!("missing {k}") }
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    match a.get(1).map(|s| s.as_str()) {
        Some("pb") if a.len() == 4 => pb::run(&a[2], &a[3]),
        Some("z") => zstat::run(&arg(&a, "--pb")?, &arg(&a, "--genes")?, &arg(&a, "--out")?),
        _ => bail!("usage: preprocess pb <cache_dir> <out_prefix> | preprocess z --pb <prefix> --genes <file> --out <prefix>"),
    }
}
