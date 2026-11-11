//! perturb-impute: choose gene panels on a reference CRISPR perturbation atlas, fit a reference
//! imputation model, and impute transcriptome-wide responses of a targeted screen.
//!
//!   perturb-impute select        --ref REF [--panel-size 500] [--mode standard|low-budget]
//!                                [--umi-fraction FILE] [--gene-ids FILE] --out panel.tsv
//!   perturb-impute fit-reference --ref REF --panel panel.tsv [--groups FILE]
//!                                [--noise FILE --budget 0.02] --out model.bin
//!   perturb-impute impute        --model model.bin --input panel_matrix.tsv --out-dir DIR
//!   perturb-impute validate      --model model.bin --input-z PREFIX --out per_gene.tsv
//!   perturb-impute project       --ref REF --input matrix.tsv [--programs 160] --out projected.tsv
//!
//! REF is either a z-score prefix (PREFIX.meta.tsv, PREFIX.genes.tsv, PREFIX.z.f32 with layers
//! all cells / half 0 / half 1) or a TSV matrix (perturbations x genes, first column = name).
#[allow(dead_code)]
#[path = "../la.rs"]
mod la;

use anyhow::{anyhow, bail, Context, Result};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::time::Instant;

const ALPHAS: [f64; 5] = [1e-3, 1e-2, 1e-1, 1.0, 10.0];
const MAGIC: &[u8; 8] = b"PIMPv001";

// ------------------------------------------------------------------ arguments
struct Args(Vec<String>);
impl Args {
    fn get(&self, k: &str) -> Option<String> { self.0.iter().position(|x| x == k).and_then(|i| self.0.get(i + 1).cloned()) }
    fn req(&self, k: &str) -> Result<String> { self.get(k).ok_or_else(|| anyhow!("missing argument {k} (see --help)")) }
}

// ------------------------------------------------------------------ data
/// Profiles: perturbations x genes, one layer (all cells) or three (all, half 0, half 1).
struct Prof { perts: Vec<String>, genes: Vec<String>, layers: Vec<Vec<f32>> }
impl Prof {
    fn load(path: &str) -> Result<Prof> {
        if std::path::Path::new(&format!("{path}.z.f32")).exists() {
            let col0 = |f: String| -> Result<Vec<String>> { Ok(std::fs::read_to_string(&f).with_context(|| f.clone())?.lines().skip(1).map(|l| l.split('\t').next().unwrap().to_string()).collect()) };
            let perts = col0(format!("{path}.meta.tsv"))?;
            let genes = col0(format!("{path}.genes.tsv"))?;
            let raw = std::fs::read(format!("{path}.z.f32"))?;
            let z: Vec<f32> = bytemuck::cast_slice::<u8, f32>(&raw).to_vec();
            let pg = perts.len() * genes.len();
            if z.len() % pg != 0 { bail!("{path}.z.f32 size does not match meta and genes") }
            let layers = z.chunks(pg).map(|c| c.to_vec()).collect();
            return Ok(Prof { perts, genes, layers });
        }
        let f = std::fs::File::open(path).with_context(|| format!("open {path}"))?;
        let mut lines = BufReader::new(f).lines();
        let hdr = lines.next().ok_or_else(|| anyhow!("empty {path}"))??;
        let genes: Vec<String> = hdr.split(['\t', ',']).skip(1).map(String::from).collect();
        let (mut perts, mut v) = (vec![], vec![]);
        for l in lines {
            let l = l?;
            let mut it = l.split(['\t', ',']);
            perts.push(it.next().unwrap().to_string());
            let row: Vec<f32> = it.map(|x| x.parse::<f32>().unwrap_or(f32::NAN)).collect();
            if row.len() != genes.len() { bail!("row {} of {path} has {} values for {} genes", perts.len(), row.len(), genes.len()) }
            v.extend(row);
        }
        Ok(Prof { perts, genes, layers: vec![v] })
    }
    fn mat(&self, layer: usize, rows: &[usize], cols: &[usize]) -> Vec<f64> {
        let (g, l) = (self.genes.len(), &self.layers[layer]);
        let mut m = vec![0.0; rows.len() * cols.len()];
        m.par_chunks_mut(cols.len().max(1)).zip(rows).for_each(|(r, &p)| { for (k, &j) in cols.iter().enumerate() { let v = l[p * g + j]; r[k] = if v.is_finite() { v as f64 } else { 0.0 }; } });
        m
    }
}

fn read_map(path: &str) -> Result<HashMap<String, String>> {
    let mut m = HashMap::new();
    for l in std::fs::read_to_string(path).with_context(|| path.to_string())?.lines() {
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() >= 2 { m.insert(f[0].to_string(), f[1].to_string()); }
    }
    Ok(m)
}

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for (a, b) in x.iter().zip(y) { sxy += (a - mx) * (b - my); sxx += (a - mx).powi(2); syy += (b - my).powi(2); }
    sxy / (sxx * syy).sqrt()
}
fn phi(x: f64) -> f64 {
    let z = x.abs() / std::f64::consts::SQRT_2;
    let t = 1.0 / (1.0 + 0.5 * z);
    let r = t * (-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))))).exp();
    if x >= 0.0 { 1.0 - 0.5 * r } else { 0.5 * r }
}

/// inner holdout: 10% of groups (or of perturbations without a group), deterministic hash
fn inner_split(p: &Prof, groups: &Option<HashMap<String, String>>) -> (Vec<usize>, Vec<usize>) {
    let key = |i: usize| -> u64 {
        let k: u64 = match groups.as_ref().and_then(|g| g.get(&p.perts[i])).and_then(|s| s.parse::<i64>().ok()) { Some(v) if v >= 0 => v as u64, _ => 1_000_000 + i as u64 };
        (k.wrapping_mul(0x9E3779B97F4A7C15) >> 32) % 10
    };
    ((0..p.perts.len()).filter(|&i| key(i) != 0).collect(), (0..p.perts.len()).filter(|&i| key(i) == 0).collect())
}

// ------------------------------------------------------------------ select
/// greedy minimum mean squared error selection on Sigma = C/n - I, gains divided by weights
fn greedy(c: &[f64], n: usize, g: usize, k: usize, w: &[f64]) -> Vec<(usize, f64)> {
    let mut r: Vec<f64> = c.iter().map(|v| v / n as f64).collect();
    for j in 0..g { r[j * g + j] = (r[j * g + j] - 1.0).max(0.0); }
    let mut norm: Vec<f64> = r.par_chunks(g).map(|row| row.iter().map(|v| v * v).sum()).collect();
    let mut used = vec![false; g];
    let mut sel = vec![];
    for _ in 0..k {
        let best = (0..g).into_par_iter().filter(|&j| !used[j] && w[j].is_finite())
            .map(|j| (norm[j] / (r[j * g + j] + 1.0) / w[j], j))
            .reduce(|| (f64::NEG_INFINITY, usize::MAX), |a, b| if b.0 > a.0 || (b.0 == a.0 && b.1 < a.1) { b } else { a });
        let j = best.1;
        if j == usize::MAX { break; }
        used[j] = true;
        let col: Vec<f64> = r[j * g..(j + 1) * g].to_vec();
        let d = col[j] + 1.0;
        sel.push((j, norm[j] / d));
        r.par_chunks_mut(g).zip(norm.par_iter_mut()).enumerate().for_each(|(i, (row, nm))| {
            let f = col[i] / d;
            let mut s = 0.0;
            for (v, cc) in row.iter_mut().zip(&col) { *v -= f * cc; s += *v * *v; }
            *nm = s;
        });
    }
    sel
}

fn cmd_select(a: &Args) -> Result<()> {
    let t0 = Instant::now();
    let p = Prof::load(&a.req("--ref")?)?;
    let k: usize = a.get("--panel-size").map(|v| v.parse()).transpose()?.unwrap_or(500);
    let mode = a.get("--mode").unwrap_or("standard".into());
    let g = p.genes.len();
    let all: Vec<usize> = (0..p.perts.len()).collect();
    let c = la::gram(&p.mat(0, &all, &(0..g).collect::<Vec<_>>()), all.len(), g);
    let umi: HashMap<String, String> = match a.get("--umi-fraction") { Some(f) => read_map(&f)?, None => HashMap::new() };
    let frac: Vec<f64> = p.genes.iter().map(|s| umi.get(s).and_then(|v| v.parse().ok()).unwrap_or(f64::NAN)).collect();
    let w: Vec<f64> = match mode.as_str() {
        "standard" => vec![1.0; g],
        "low-budget" => {
            if umi.is_empty() { bail!("--mode low-budget needs --umi-fraction (gene, control-cell UMI fraction)") }
            frac.iter().map(|u| if u.is_finite() { (u + 1e-4).sqrt() } else { f64::INFINITY }).collect()
        }
        m => bail!("unknown --mode {m} (standard or low-budget)"),
    };
    let total: f64 = (0..g).map(|j| (c[j * g + j] / all.len() as f64 - 1.0).max(0.0)).sum();
    let sel = greedy(&c, all.len(), g, k, &w);
    let ids: HashMap<String, String> = match a.get("--gene-ids") { Some(f) => read_map(&f)?, None => HashMap::new() };
    let mut out = BufWriter::new(std::fs::File::create(a.req("--out")?)?);
    writeln!(out, "rank\tgene_symbol\tgene_id\tcontrol_umi_fraction\tcumulative_umi_fraction\tvariance_explained\tcumulative_variance_explained")?;
    let (mut cu, mut cv) = (0.0, 0.0);
    for (r, &(j, gain)) in sel.iter().enumerate() {
        if frac[j].is_finite() { cu += frac[j]; }
        cv += gain / total;
        let f = |x: f64| if x.is_finite() { format!("{x:.6e}") } else { "NA".into() };
        writeln!(out, "{}\t{}\t{}\t{}\t{}\t{:.6e}\t{:.6}", r + 1, p.genes[j], ids.get(&p.genes[j]).map(|s| s.as_str()).unwrap_or("NA"), f(frac[j]), if umi.is_empty() { "NA".into() } else { format!("{cu:.6}") }, gain / total, cv)?;
    }
    eprintln!("select: {} genes ({mode}) from {} perturbations x {g} genes in {:.1}s; variance explained {:.3}", sel.len(), all.len(), t0.elapsed().as_secs_f64(), cv);
    Ok(())
}

// ------------------------------------------------------------------ model
struct Model {
    n_ref: u64, alpha: f64, budget: f64,
    panel: Vec<String>, targets: Vec<String>,
    css: Vec<f64>,     // k x k reference Gram, panel block
    cus: Vec<f64>,     // nu x k reference Gram, target x panel
    extra: Vec<f64>,   // k, extra panel noise on the diagonal (Gram units)
    scale: Vec<f64>,   // k, multiplier applied to input panel z-scores
    info: Vec<f64>, slope: Vec<f64>, errvar: Vec<f64>, // nu
}
fn wstr(w: &mut impl Write, s: &str) -> Result<()> { w.write_all(&(s.len() as u32).to_le_bytes())?; w.write_all(s.as_bytes())?; Ok(()) }
fn wvec(w: &mut impl Write, v: &[f64]) -> Result<()> { w.write_all(&(v.len() as u64).to_le_bytes())?; for x in v { w.write_all(&x.to_le_bytes())?; } Ok(()) }
impl Model {
    fn save(&self, path: &str) -> Result<()> {
        let mut w = BufWriter::new(std::fs::File::create(path)?);
        w.write_all(MAGIC)?;
        w.write_all(&self.n_ref.to_le_bytes())?; w.write_all(&self.alpha.to_le_bytes())?; w.write_all(&self.budget.to_le_bytes())?;
        w.write_all(&(self.panel.len() as u64).to_le_bytes())?; for s in &self.panel { wstr(&mut w, s)?; }
        w.write_all(&(self.targets.len() as u64).to_le_bytes())?; for s in &self.targets { wstr(&mut w, s)?; }
        for v in [&self.css, &self.cus, &self.extra, &self.scale, &self.info, &self.slope, &self.errvar] { wvec(&mut w, v)?; }
        Ok(())
    }
    fn load(path: &str) -> Result<Model> {
        let mut b = vec![];
        std::fs::File::open(path).with_context(|| path.to_string())?.read_to_end(&mut b)?;
        if &b[..8] != MAGIC { bail!("{path} is not a perturb-impute model") }
        let mut o = 8;
        let u64r = |o: &mut usize| -> u64 { let v = u64::from_le_bytes(b[*o..*o + 8].try_into().unwrap()); *o += 8; v };
        let f64r = |o: &mut usize| -> f64 { let v = f64::from_le_bytes(b[*o..*o + 8].try_into().unwrap()); *o += 8; v };
        let sr = |o: &mut usize| -> String { let n = u32::from_le_bytes(b[*o..*o + 4].try_into().unwrap()) as usize; let s = String::from_utf8_lossy(&b[*o + 4..*o + 4 + n]).into_owned(); *o += 4 + n; s };
        let n_ref = u64r(&mut o); let alpha = f64r(&mut o); let budget = f64r(&mut o);
        let np = u64r(&mut o) as usize; let panel = (0..np).map(|_| sr(&mut o)).collect();
        let nt = u64r(&mut o) as usize; let targets = (0..nt).map(|_| sr(&mut o)).collect();
        let mut vr = || -> Vec<f64> { let n = u64r(&mut o) as usize; (0..n).map(|_| f64r(&mut o)).collect() };
        Ok(Model { n_ref, alpha, budget, panel, targets, css: vr(), cus: vr(), extra: vr(), scale: vr(), info: vr(), slope: vr(), errvar: vr() })
    }
    /// imputation weights W' (nu x m) for the available panel genes `av` (indices into panel)
    fn weights(&self, av: &[usize]) -> Result<Vec<f64>> {
        let (k, m, nu) = (self.panel.len(), av.len(), self.targets.len());
        let mut a = vec![0.0; m * m];
        for (i, &p) in av.iter().enumerate() { for (j, &q) in av.iter().enumerate() { a[i * m + j] = self.css[p * k + q]; } }
        let tr: f64 = (0..m).map(|i| a[i * m + i]).sum::<f64>() / m as f64;
        for (i, &p) in av.iter().enumerate() { a[i * m + i] += self.alpha * tr + self.extra[p]; }
        let mut bt = vec![0.0; nu * m];
        bt.par_chunks_mut(m).enumerate().for_each(|(u, r)| { for (j, &p) in av.iter().enumerate() { r[j] = self.cus[u * k + p]; } });
        la::posv(a, m, &mut bt, nu)?;
        Ok(bt)
    }
    /// out-of-reference score: fraction of each input panel vector outside the leading programs of the panel block
    fn oor(&self, av: &[usize], x: &[f64], n: usize) -> Result<Vec<f64>> {
        let (k, m) = (self.panel.len(), av.len());
        let mut a = vec![0.0; m * m];
        for (i, &p) in av.iter().enumerate() { for (j, &q) in av.iter().enumerate() { a[i * m + j] = self.css[p * k + q]; } }
        let kq = 160.min(m / 2).max(1);
        let (pv, _) = la::top_eig(&a, m, kq)?;
        let pr = la::gemm(false, false, n, kq, m, x, &pv);
        Ok((0..n).map(|i| { let nx: f64 = x[i * m..(i + 1) * m].iter().map(|v| v * v).sum(); if nx > 0.0 { 1.0 - pr[i * kq..(i + 1) * kq].iter().map(|v| v * v).sum::<f64>() / nx } else { f64::NAN } }).collect())
    }
}

fn cmd_fit(a: &Args) -> Result<()> {
    let t0 = Instant::now();
    let p = Prof::load(&a.req("--ref")?)?;
    let g = p.genes.len();
    let gi: HashMap<&str, usize> = p.genes.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let ptxt = std::fs::read_to_string(a.req("--panel")?)?;
    let mut plines = ptxt.lines();
    let hdr: Vec<&str> = plines.next().unwrap_or("").split('\t').collect();
    let col = hdr.iter().position(|h| *h == "gene_symbol").unwrap_or(0);
    let mut s: Vec<usize> = vec![];
    for l in std::iter::once(if hdr.iter().any(|h| *h == "gene_symbol") { "" } else { ptxt.lines().next().unwrap_or("") }).chain(plines) {
        if l.is_empty() { continue; }
        let sym = l.split('\t').nth(col).unwrap_or("");
        match gi.get(sym) { Some(&j) => { if !s.contains(&j) { s.push(j) } } None => eprintln!("fit-reference: panel gene {sym} not in the reference, skipped") }
    }
    let k = s.len();
    let inp: HashSet<usize> = s.iter().copied().collect();
    let u: Vec<usize> = (0..g).filter(|j| !inp.contains(j)).collect();
    let nu = u.len();
    let groups = match a.get("--groups") { Some(f) => Some(read_map(&f)?), None => None };
    let (itr, iho) = inner_split(&p, &groups);
    let all: Vec<usize> = (0..p.perts.len()).collect();
    let cols: Vec<usize> = (0..g).collect();
    let c = la::gram(&p.mat(0, &all, &cols), all.len(), g);
    let cin = la::gram(&p.mat(0, &itr, &cols), itr.len(), g);
    let halves = p.layers.len() >= 3;
    let (lin, ltr) = if halves { (1, 2) } else { (0, 0) };
    // optional heteroscedastic panel noise at a targeted sequencing budget
    let budget: f64 = a.get("--budget").map(|v| v.parse()).transpose()?.unwrap_or(1.0);
    let (mut extra, mut scale) = (vec![0.0; k], vec![1.0; k]);
    if let Some(nf) = a.get("--noise") {
        // columns: gene, cv2_bio, cv2_tech, cv2_tot, share
        let mut rows: HashMap<String, Vec<f64>> = HashMap::new();
        for l in std::fs::read_to_string(&nf)?.lines().skip(1) { let f: Vec<&str> = l.split('\t').collect(); rows.insert(f[0].into(), f[1..].iter().map(|x| x.parse().unwrap_or(f64::NAN)).collect()); }
        let share: f64 = s.iter().filter_map(|&j| rows.get(&p.genes[j])).map(|r| r[3]).filter(|x| x.is_finite()).sum();
        let rate = (budget / share).min(1.0);
        for (i, &j) in s.iter().enumerate() {
            if let Some(r) = rows.get(&p.genes[j]) {
                let (bio, tech, tot) = (r[0], r[1], r[2]);
                if tot.is_finite() && tot > 0.0 { let inv_f = (bio + tech / rate) / tot; extra[i] = all.len() as f64 * (inv_f - 1.0).max(0.0); scale[i] = inv_f.sqrt(); }
            }
        }
        eprintln!("fit-reference: budget {budget}, panel UMI share {share:.4}, targeted rate {rate:.4}");
    }
    let mk = |cm: &[f64], n: u64, alpha: f64| Model { n_ref: n, alpha, budget, panel: s.iter().map(|&j| p.genes[j].clone()).collect(), targets: u.iter().map(|&j| p.genes[j].clone()).collect(),
        css: s.iter().flat_map(|&i| s.iter().map(move |&j| cm[i * g + j])).collect(), cus: u.iter().flat_map(|&i| s.iter().map(move |&j| cm[i * g + j])).collect(),
        extra: extra.iter().map(|e| e * n as f64 / all.len() as f64).collect(), scale: scale.clone(), info: vec![], slope: vec![], errvar: vec![] };
    // alpha, INFO, calibration on the inner holdout
    let xh: Vec<f64> = { let x = p.mat(lin, &iho, &s); x.chunks(k).flat_map(|r| r.iter().zip(&scale).map(|(v, sc)| v * sc).collect::<Vec<_>>()).collect() };
    let th = p.mat(ltr, &iho, &u);
    let nh = iho.len();
    let av: Vec<usize> = (0..k).collect();
    let mut best = (f64::NEG_INFINITY, ALPHAS[0], vec![]);
    for &al in &ALPHAS {
        let wt = mk(&cin, itr.len() as u64, al).weights(&av)?;
        let im = la::gemm(false, true, nh, nu, k, &xh, &wt);
        let r: f64 = (0..nh).map(|i| pearson(&im[i * nu..(i + 1) * nu], &th[i * nu..(i + 1) * nu])).filter(|x| x.is_finite()).sum::<f64>() / nh as f64;
        if r > best.0 { best = (r, al, im); }
    }
    let (_, alpha, im) = best;
    let col = |m: &[f64], j: usize| -> Vec<f64> { (0..nh).map(|i| m[i * nu + j]).collect() };
    let stats: Vec<(f64, f64, f64)> = (0..nu).into_par_iter().map(|j| {
        let (x, y) = (col(&im, j), col(&th, j));
        let sxx: f64 = x.iter().map(|v| v * v).sum();
        let b = if sxx > 0.0 { x.iter().zip(&y).map(|(p, q)| p * q).sum::<f64>() / sxx } else { 0.0 };
        let mse = x.iter().zip(&y).map(|(p, q)| (q - b * p).powi(2)).sum::<f64>() / nh as f64;
        (pearson(&x, &y), b, if halves { (mse - 1.0).max(0.05) } else { mse.max(0.05) })
    }).collect();
    let mut model = mk(&c, all.len() as u64, alpha);
    model.info = stats.iter().map(|s| s.0).collect();
    model.slope = stats.iter().map(|s| s.1).collect();
    model.errvar = stats.iter().map(|s| s.2).collect();
    let out = a.req("--out")?;
    model.save(&out)?;
    let size = std::fs::metadata(&out)?.len() as f64 / 1e6;
    eprintln!("fit-reference: {} reference perturbations ({} inner training, {nh} inner holdout{}), {k} panel genes, {nu} targets, alpha {alpha}, median INFO {:.3}; {out} {size:.1} MB in {:.1}s",
        all.len(), itr.len(), if halves { ", half 0 -> half 1" } else { ", single layer: INFO is not cross-cell" }, { let mut v: Vec<f64> = model.info.iter().copied().filter(|x| x.is_finite()).collect(); v.sort_by(|a, b| a.partial_cmp(b).unwrap()); v[v.len() / 2] }, t0.elapsed().as_secs_f64());
    Ok(())
}

// ------------------------------------------------------------------ impute
/// panel matrix of the input aligned to the model; returns (perturbations, available panel indices, n x m matrix)
fn align(m: &Model, p: &Prof, layer: usize) -> (Vec<usize>, Vec<f64>) {
    let gi: HashMap<&str, usize> = p.genes.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let av: Vec<usize> = (0..m.panel.len()).filter(|&i| gi.contains_key(m.panel[i].as_str())).collect();
    let cols: Vec<usize> = av.iter().map(|&i| gi[m.panel[i].as_str()]).collect();
    let rows: Vec<usize> = (0..p.perts.len()).collect();
    let mut x = p.mat(layer, &rows, &cols);
    let mm = av.len();
    x.par_chunks_mut(mm).for_each(|r| { for (c, &i) in av.iter().enumerate() { r[c] *= m.scale[i]; } });
    (av, x)
}

fn cmd_impute(a: &Args) -> Result<()> {
    let t0 = Instant::now();
    let m = Model::load(&a.req("--model")?)?;
    let p = match a.get("--input-z") { Some(z) => Prof::load(&z)?, None => Prof::load(&a.req("--input")?)? };
    let layer: usize = a.get("--layer").map(|v| v.parse()).transpose()?.unwrap_or(0);
    let (av, x) = align(&m, &p, layer);
    if let Some(f) = a.get("--write-input") {
        let mut w = BufWriter::new(std::fs::File::create(&f)?);
        writeln!(w, "perturbation\t{}", av.iter().map(|&i| m.panel[i].as_str()).collect::<Vec<_>>().join("\t"))?;
        for i in 0..p.perts.len() { writeln!(w, "{}\t{}", p.perts[i], x[i * av.len()..(i + 1) * av.len()].iter().zip(&av).map(|(v, &g)| format!("{}", (v / m.scale[g]) as f32)).collect::<Vec<_>>().join("\t"))?; }
    }
    if av.len() < m.panel.len() / 2 { bail!("only {} of {} panel genes found in the input", av.len(), m.panel.len()) }
    let (n, mm, nu) = (p.perts.len(), av.len(), m.targets.len());
    let t1 = Instant::now();
    let wt = m.weights(&av)?;
    let mut im = la::gemm(false, true, n, nu, mm, &x, &wt);
    im.par_chunks_mut(nu).for_each(|r| for j in 0..nu { r[j] *= m.slope[j]; });
    let oor = m.oor(&av, &x, n)?;
    let tcomp = t1.elapsed().as_secs_f64();
    let dir = a.req("--out-dir")?;
    std::fs::create_dir_all(&dir)?;
    let write_mat = |name: &str, f: &(dyn Fn(usize, usize) -> f64 + Sync), prec: usize| -> Result<()> {
        let mut w = BufWriter::new(std::fs::File::create(format!("{dir}/{name}"))?);
        writeln!(w, "perturbation\t{}", m.targets.join("\t"))?;
        let rows: Vec<String> = (0..n).into_par_iter().map(|i| { let mut s = p.perts[i].clone(); for j in 0..nu { s.push('\t'); s.push_str(&format!("{:.*}", prec, f(i, j))); } s }).collect();
        for r in rows { writeln!(w, "{r}")?; }
        Ok(())
    };
    write_mat("imputed_zscores.tsv", &|i, j| im[i * nu + j], 4)?;
    write_mat("imputed_lfsr.tsv", &|i, j| phi(-im[i * nu + j].abs() / m.errvar[j].sqrt()), 4)?;
    let mut w = BufWriter::new(std::fs::File::create(format!("{dir}/perturbation_qc.tsv"))?);
    writeln!(w, "perturbation\toor_score\tpanel_norm\tpanel_genes_used")?;
    for i in 0..n { writeln!(w, "{}\t{:.6}\t{:.6}\t{}", p.perts[i], oor[i], x[i * mm..(i + 1) * mm].iter().map(|v| v * v).sum::<f64>().sqrt(), mm)?; }
    let mut w = BufWriter::new(std::fs::File::create(format!("{dir}/gene_info.tsv"))?);
    writeln!(w, "gene\tinfo\tcalibration_slope\terror_variance")?;
    for j in 0..nu { writeln!(w, "{}\t{:.6}\t{:.6}\t{:.6}", m.targets[j], m.info[j], m.slope[j], m.errvar[j])?; }
    eprintln!("impute: {n} perturbations, {mm} of {} panel genes, {nu} genes imputed; computation {tcomp:.2}s, total {:.1}s", m.panel.len(), t0.elapsed().as_secs_f64());
    if mm < m.panel.len() { eprintln!("impute: {} panel genes missing in the input, weights solved on the available genes; INFO and calibration refer to the full panel", m.panel.len() - mm); }
    Ok(())
}

/// per-gene correlation between imputed (input layer) and measured (truth layer) responses
fn cmd_validate(a: &Args) -> Result<()> {
    let m = Model::load(&a.req("--model")?)?;
    let p = Prof::load(&a.req("--input-z")?)?;
    let lin: usize = a.get("--layer").map(|v| v.parse()).transpose()?.unwrap_or(1);
    let ltr: usize = a.get("--truth-layer").map(|v| v.parse()).transpose()?.unwrap_or(2);
    let min_r: f64 = a.get("--min-r").map(|v| v.parse()).transpose()?.unwrap_or(0.2);
    let (av, x) = align(&m, &p, lin);
    let (n, mm, nu) = (p.perts.len(), av.len(), m.targets.len());
    let wt = m.weights(&av)?;
    let im = la::gemm(false, true, n, nu, mm, &x, &wt);
    let gi: HashMap<&str, usize> = p.genes.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let tj: Vec<Option<usize>> = m.targets.iter().map(|s| gi.get(s.as_str()).copied()).collect();
    let present: Vec<usize> = (0..nu).filter(|&j| tj[j].is_some()).collect();
    let rows: Vec<usize> = (0..n).collect();
    let tr = p.mat(ltr, &rows, &present.iter().map(|&j| tj[j].unwrap()).collect::<Vec<_>>());
    let np = present.len();
    // responding perturbations: imputed and measured profiles both available, measured halves agree (r >= min_r)
    let h0 = p.mat(1.min(p.layers.len() - 1), &rows, &present.iter().map(|&j| tj[j].unwrap()).collect::<Vec<_>>());
    let resp: Vec<usize> = (0..n).filter(|&i| pearson(&h0[i * np..(i + 1) * np], &tr[i * np..(i + 1) * np]) >= min_r).collect();
    let r: Vec<f64> = (0..np).into_par_iter().map(|q| { let j = present[q]; pearson(&resp.iter().map(|&i| im[i * nu + j]).collect::<Vec<_>>(), &resp.iter().map(|&i| tr[i * np + q]).collect::<Vec<_>>()) }).collect();
    let mut w = BufWriter::new(std::fs::File::create(a.req("--out")?)?);
    writeln!(w, "gene\tobserved_r")?;
    let rmap: HashMap<usize, f64> = present.iter().copied().zip(r.iter().copied()).collect();
    for j in 0..nu { writeln!(w, "{}\t{}", m.targets[j], rmap.get(&j).map(|v| format!("{v:.6}")).unwrap_or("NA".into()))?; }
    let mut v: Vec<f64> = r.iter().copied().filter(|x| x.is_finite()).collect();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    eprintln!("validate: {} perturbations ({} responding), {mm} of {} panel genes, {np} target genes measured; median per-gene r {:.3}", n, resp.len(), m.panel.len(), v[v.len() / 2]);
    Ok(())
}

/// Orthogonal projection of profiles (n x g, row-major) on the k leading eigenvectors of the
/// second-moment matrix of the reference profiles (the reference programs).
fn project(reference: &[f64], n_ref: usize, g: usize, x: &[f64], n: usize, k: usize) -> Result<Vec<f64>> {
    let c = la::gram(reference, n_ref, g);
    // exact eigenvectors for small gene sets, block subspace iteration for transcriptome-wide ones
    let (v, _) = if g <= 1500 { la::top_eig(&c, g, k)? } else { la::top_eig_sub(&c, g, k, 40)? };
    Ok(la::gemm(false, true, n, g, k, &la::gemm(false, false, n, k, g, x, &v), &v))
}

/// project: whole-transcriptome profiles denoised by projection on the reference programs.
/// Genes are matched by name; reference genes absent from the input enter as zero.
fn cmd_project(a: &Args) -> Result<()> {
    let t0 = Instant::now();
    let r = Prof::load(&a.req("--ref")?)?;
    let x = Prof::load(&a.req("--input")?)?;
    let k: usize = a.get("--programs").map(|v| v.parse()).transpose()?.unwrap_or(160);
    let (n, g) = (r.perts.len(), r.genes.len());
    if k == 0 || k > g.min(n) { bail!("--programs {k} must be between 1 and {}", g.min(n)) }
    let y = r.mat(0, &(0..n).collect::<Vec<_>>(), &(0..g).collect::<Vec<_>>());
    let gi: HashMap<&str, usize> = x.genes.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let (nx, gx) = (x.perts.len(), x.genes.len());
    let mut m = vec![0.0; nx * g];
    let mut matched = 0;
    for (j, name) in r.genes.iter().enumerate() {
        if let Some(&c) = gi.get(name.as_str()) {
            matched += 1;
            for i in 0..nx { let v = x.layers[0][i * gx + c]; if v.is_finite() { m[i * g + j] = v as f64; } }
        }
    }
    if matched < g / 2 { bail!("only {matched} of {g} reference genes found in the input") }
    let p = project(&y, n, g, &m, nx, k)?;
    let mut w = BufWriter::new(std::fs::File::create(a.req("--out")?)?);
    writeln!(w, "perturbation\t{}", r.genes.join("\t"))?;
    for i in 0..nx { writeln!(w, "{}\t{}", x.perts[i], p[i * g..(i + 1) * g].iter().map(|v| format!("{:.5}", v)).collect::<Vec<_>>().join("\t"))?; }
    eprintln!("project: {nx} profiles, {matched} of {g} reference genes matched, {k} programs, {:.1}s", t0.elapsed().as_secs_f64());
    Ok(())
}

const HELP: &str = "perturb-impute: gene panels and transcriptome-wide imputation for CRISPR perturbation screens

  select        --ref REF --out panel.tsv [--panel-size 500] [--mode standard|low-budget]
                [--umi-fraction gene_fraction.tsv] [--gene-ids symbol_id.tsv]
  fit-reference --ref REF --panel panel.tsv --out model.bin [--groups pert_group.tsv]
                [--noise gene_cv2.tsv --budget 0.02]
  impute        --model model.bin --input matrix.tsv --out-dir DIR      (or --input-z PREFIX [--layer 0])
  validate      --model model.bin --input-z PREFIX --out per_gene_r.tsv [--layer 1 --truth-layer 2]
  project       --ref REF --input matrix.tsv --out projected.tsv [--programs 160]
                (whole-transcriptome profiles projected on the leading reference programs)

REF: z-score prefix (PREFIX.meta.tsv, PREFIX.genes.tsv, PREFIX.z.f32; layers all cells, half 0, half 1)
     or a TSV/CSV matrix with perturbations in rows, genes in columns and names in the first column.";

fn main() -> Result<()> {
    let v: Vec<String> = std::env::args().collect();
    let a = Args(v.clone());
    match v.get(1).map(|s| s.as_str()) {
        Some("select") => cmd_select(&a),
        Some("fit-reference") => cmd_fit(&a),
        Some("impute") => cmd_impute(&a),
        Some("validate") => cmd_validate(&a),
        Some("project") => cmd_project(&a),
        _ => { println!("{HELP}"); Ok(()) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    /// a model saved to disk loads back bit-identically
    #[test]
    fn model_roundtrip() {
        let m = Model { n_ref: 7, alpha: 0.1, budget: 0.02, panel: vec!["A".into(), "B".into()], targets: vec!["C".into()], css: vec![2.0, 0.5, 0.5, 1.0], cus: vec![0.3, -0.2],
            extra: vec![0.0, 1.5], scale: vec![1.0, 2.0], info: vec![0.4], slope: vec![0.9], errvar: vec![0.2] };
        let p = std::env::temp_dir().join("pimp_test.bin");
        m.save(p.to_str().unwrap()).unwrap();
        let l = Model::load(p.to_str().unwrap()).unwrap();
        assert_eq!((l.n_ref, l.alpha, l.budget, &l.panel, &l.targets), (7, 0.1, 0.02, &m.panel, &m.targets));
        assert_eq!((&l.css, &l.cus, &l.extra, &l.scale, &l.info, &l.slope, &l.errvar), (&m.css, &m.cus, &m.extra, &m.scale, &m.info, &m.slope, &m.errvar));
    }
    /// projection on reference programs returns profiles inside the program space unchanged and removes the rest
    #[test]
    fn projection_keeps_program_space() {
        // reference: 6 profiles spanned by two programs over 4 genes
        let (p1, p2) = ([1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, -1.0]);
        let coef = [(3.0, 0.5), (-2.0, 1.0), (1.0, -4.0), (0.5, 2.0), (-1.5, -1.0), (2.5, 3.0)];
        let y: Vec<f64> = coef.iter().flat_map(|&(a, b)| (0..4).map(move |j| a * p1[j] + b * p2[j])).collect();
        let inside = [2.0, 2.0, -1.0, 1.0];
        let outside = [1.0, -1.0, 1.0, 1.0];
        let x: Vec<f64> = inside.iter().chain(outside.iter()).copied().collect();
        let p = project(&y, 6, 4, &x, 2, 2).unwrap();
        for j in 0..4 { assert!((p[j] - inside[j]).abs() < 1e-8, "inside {j}: {}", p[j]); assert!(p[4 + j].abs() < 1e-8, "outside {j}: {}", p[4 + j]); }
    }
    /// weights solved on a subset equal the ridge solution written out directly
    #[test]
    fn subset_weights_match_direct_solve() {
        let m = Model { n_ref: 1, alpha: 0.0, budget: 1.0, panel: vec!["A".into(), "B".into(), "C".into()], targets: vec!["T".into()],
            css: vec![4.0, 1.0, 0.5, 1.0, 3.0, 0.2, 0.5, 0.2, 2.0], cus: vec![1.0, 2.0, 3.0], extra: vec![0.0; 3], scale: vec![1.0; 3], info: vec![0.0], slope: vec![1.0], errvar: vec![1.0] };
        let w = m.weights(&[0, 2]).unwrap();
        // [[4, .5], [.5, 2]] w = [1, 3]
        let det = 4.0 * 2.0 - 0.25;
        let (w0, w1) = ((2.0 * 1.0 - 0.5 * 3.0) / det, (4.0 * 3.0 - 0.5 * 1.0) / det);
        assert!((w[0] - w0).abs() < 1e-12 && (w[1] - w1).abs() < 1e-12);
    }
}
