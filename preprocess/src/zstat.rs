//! z statistics from a pseudobulk (layout of `preprocess pb`: <pb>.f32 = 3 x P x G mean CP10k for
//! all cells, half 0 (even member positions), half 1; <pb>.groups.tsv name\tn; <pb>.genes.txt).
//! Null per-cell variance per gene from half differences, pooled over perturbations:
//! s2_j = median_p[(m0 - m1)^2 / (1/n0 + 1/n1)] / median(chi2_1).
//! z = (m_p - m_c) / sqrt(s2_j (1/n_p + 1/n_c)), for all cells and each half.
use anyhow::{bail, Result};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufWriter, Write};

pub const MED_CHI2_1: f64 = 0.454936423119572;

pub fn median(v: &mut [f64]) -> f64 {
    let n = v.len();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if n % 2 == 1 { v[n / 2] } else { 0.5 * (v[n / 2 - 1] + v[n / 2]) }
}

pub fn run(pb: &str, genes_file: &str, out: &str) -> Result<()> {
    let t = std::time::Instant::now();
    let groups: Vec<(String, usize)> = std::fs::read_to_string(format!("{pb}.groups.tsv"))?
        .lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_string(), f[1].parse().unwrap())
        })
        .collect();
    let genes: Vec<String> = std::fs::read_to_string(format!("{pb}.genes.txt"))?.lines().map(String::from).collect();
    let (p, g) = (groups.len(), genes.len());
    let raw = std::fs::read(format!("{pb}.f32"))?;
    let x: &[f32] = bytemuck::cast_slice(&raw);
    if x.len() != 3 * p * g { bail!("size mismatch") }
    let at = |k: usize, i: usize, j: usize| x[(k * p + i) * g + j] as f64;
    let gi: HashMap<&str, usize> = genes.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let want: Vec<usize> = std::fs::read_to_string(genes_file)?.lines().filter_map(|s| gi.get(s).copied()).collect();
    let Some(ci) = groups.iter().position(|g| g.0 == "control") else { bail!("no control") };
    let nc = groups[ci].1 as f64;
    let halves = |n: usize| (((n + 1) / 2) as f64, (n / 2) as f64);
    let perts: Vec<usize> = (0..p).filter(|&i| i != ci && groups[i].1 >= 2).collect();
    let s2: Vec<f64> = want
        .par_iter()
        .map(|&j| {
            let mut v: Vec<f64> = perts.iter().map(|&i| {
                let (n0, n1) = halves(groups[i].1);
                (at(1, i, j) - at(2, i, j)).powi(2) / (1.0 / n0 + 1.0 / n1)
            }).collect();
            median(&mut v) / MED_CHI2_1
        })
        .collect();
    let keep: Vec<usize> = (0..want.len()).filter(|&k| s2[k] > 0.0).collect();
    let gk = keep.len();
    let mut z = vec![0f32; 3 * perts.len() * gk];
    for k in 0..3 {
        z[k * perts.len() * gk..(k + 1) * perts.len() * gk].par_chunks_mut(gk).zip(&perts).for_each(|(row, &i)| {
            let n = groups[i].1 as f64;
            let (n0, n1) = halves(groups[i].1);
            let nn = [n, n0, n1][k];
            for (c, &kk) in keep.iter().enumerate() {
                let j = want[kk];
                row[c] = ((at(k, i, j) - at(0, ci, j)) / (s2[kk] * (1.0 / nn + 1.0 / nc)).sqrt()) as f32;
            }
        });
    }
    std::fs::write(format!("{out}.z.f32"), bytemuck::cast_slice(&z))?;
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.meta.tsv"))?);
    writeln!(w, "pert\tn\tneff_all\tneff_h0\tneff_h1")?;
    for &i in &perts {
        let n = groups[i].1 as f64;
        let (n0, n1) = halves(groups[i].1);
        let ne = |a: f64| 1.0 / (1.0 / a + 1.0 / nc);
        writeln!(w, "{}\t{}\t{}\t{}\t{}", groups[i].0, groups[i].1, ne(n), ne(n0), ne(n1))?;
    }
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.genes.tsv"))?);
    writeln!(w, "gene\ts2\tctrl_mean")?;
    for &kk in &keep {
        writeln!(w, "{}\t{}\t{}", genes[want[kk]], s2[kk], at(0, ci, want[kk]))?;
    }
    eprintln!("z: {} perturbations x {gk} genes ({} dropped, s2 = 0) in {:.1}s", perts.len(), want.len() - gk, t.elapsed().as_secs_f64());
    Ok(())
}
