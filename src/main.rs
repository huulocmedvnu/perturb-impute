//! panelimp: analysis pipeline for panel imputation of perturbation profiles.
//!   panelimp cv       --z PREFIX --folds folds.tsv --hubs hub_excl.tsv --l1000 FILE --out PREFIX [--sizes 100,200,500,718] [--folds-run 0,1,..]
//!   panelimp transfer --z PREFIX --zt PREFIX --l1000 FILE --out PREFIX [--sizes ...]
//! z files: coexsc format (<p>.meta.tsv, <p>.genes.tsv, <p>.z.f32 = 3 x P x G: all, half 0, half 1).
mod cache;
mod la;
mod sim;
mod upg;
mod de;
mod folds;
mod multi;
mod vcm;
mod real;
mod causal;
mod explore;
mod explore_test;

use anyhow::{anyhow, bail, Result};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};

const NPC: usize = 10; // global programs removed for the specific metric
const ALPHAS: [f64; 5] = [1e-3, 1e-2, 1e-1, 1.0, 10.0];
const DEN_K: [usize; 6] = [20, 50, 100, 160, 300, 500];

/// per-fold reference Gram matrices and programs
struct Grams {
    c: Vec<f64>,
    c_in: Vec<f64>,
    v: Vec<f64>,  // g x NPC
    vk: Vec<f64>, // g x max(DEN_K)
    ev: Vec<f64>, // top max(DEN_K) eigenvalues, descending
}
/// Gram matrices only (no eigendecomposition), for callers that use c and c_in alone
fn grams_c(z: &Z, refs: &[usize], itr: &[usize]) -> Grams {
    let g = z.genes.len();
    let all: Vec<usize> = (0..g).collect();
    Grams { c: la::gram(&z.mat(0, refs, &all), refs.len(), g), c_in: la::gram(&z.mat(0, itr, &all), itr.len(), g), v: vec![], vk: vec![], ev: vec![] }
}
fn grams(z: &Z, refs: &[usize], itr: &[usize]) -> Result<Grams> {
    let g = z.genes.len();
    let all: Vec<usize> = (0..g).collect();
    let c = la::gram(&z.mat(0, refs, &all), refs.len(), g);
    let c_in = la::gram(&z.mat(0, itr, &all), itr.len(), g);
    let nk = *DEN_K.iter().max().unwrap();
    let (vk, ev) = la::top_eig(&c, g, nk)?;
    let v: Vec<f64> = (0..g).flat_map(|i| vk[i * nk..i * nk + NPC].to_vec()).collect();
    Ok(Grams { c, c_in, v, vk, ev })
}

fn arg(a: &[String], k: &str) -> Option<String> {
    a.iter().position(|x| x == k).and_then(|i| a.get(i + 1).cloned())
}
fn req(a: &[String], k: &str) -> Result<String> {
    arg(a, k).ok_or_else(|| anyhow!("missing {k}"))
}

struct Z {
    perts: Vec<String>,
    genes: Vec<String>,
    z: Vec<f32>,
}
impl Z {
    fn load(pre: &str) -> Result<Z> {
        let perts = std::fs::read_to_string(format!("{pre}.meta.tsv"))?.lines().skip(1).map(|l| l.split('\t').next().unwrap().to_string()).collect::<Vec<_>>();
        let genes = std::fs::read_to_string(format!("{pre}.genes.tsv"))?.lines().skip(1).map(|l| l.split('\t').next().unwrap().to_string()).collect::<Vec<_>>();
        let raw = std::fs::read(format!("{pre}.z.f32"))?;
        let z = bytemuck::cast_slice::<u8, f32>(&raw).to_vec();
        if z.len() != 3 * perts.len() * genes.len() { bail!("z size") }
        Ok(Z { perts, genes, z })
    }
    fn at(&self, layer: usize, p: usize, g: usize) -> f64 {
        self.z[(layer * self.perts.len() + p) * self.genes.len() + g] as f64
    }
    /// rows x cols matrix (row-major) for given perturbations and genes
    fn mat(&self, layer: usize, rows: &[usize], cols: &[usize]) -> Vec<f64> {
        let mut m = vec![0.0; rows.len() * cols.len()];
        m.par_chunks_mut(cols.len()).zip(rows).for_each(|(r, &p)| {
            for (k, &g) in cols.iter().enumerate() {
                r[k] = self.at(layer, p, g);
            }
        });
        m
    }
}

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for (a, b) in x.iter().zip(y) {
        sxy += (a - mx) * (b - my);
        sxx += (a - mx).powi(2);
        syy += (b - my).powi(2);
    }
    sxy / (sxx * syy).sqrt()
}

/// average ranks (1-based)
fn ranks(x: &[f64]) -> Vec<f64> {
    let mut o: Vec<usize> = (0..x.len()).collect();
    o.sort_by(|&a, &b| x[a].partial_cmp(&x[b]).unwrap());
    let mut r = vec![0.0; x.len()];
    let mut i = 0;
    while i < o.len() {
        let mut j = i;
        while j + 1 < o.len() && x[o[j + 1]] == x[o[i]] { j += 1; }
        for k in i..=j { r[o[k]] = (i + j) as f64 / 2.0 + 1.0; }
        i = j + 1;
    }
    r
}
fn spearman(x: &[f64], y: &[f64]) -> f64 { pearson(&ranks(x), &ranks(y)) }
/// standard normal CDF via erfc (Numerical Recipes erfcc, rel. error < 1.2e-7)
fn phi(x: f64) -> f64 {
    let z = x.abs() / std::f64::consts::SQRT_2;
    let t = 1.0 / (1.0 + 0.5 * z);
    let r = t * (-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))))).exp();
    if x >= 0.0 { 1.0 - 0.5 * r } else { 0.5 * r }
}
const DE_Z: f64 = 4.0;
const LFSR: f64 = 0.01;

/// AUC of score for labels (Mann-Whitney with average ranks); None if < 5 of either class.
fn auc(score: &[f64], lab: &[bool]) -> Option<f64> {
    let np = lab.iter().filter(|&&l| l).count();
    let nn = lab.len() - np;
    if np < 5 || nn < 5 { return None; }
    let mut o: Vec<usize> = (0..score.len()).collect();
    o.sort_by(|&a, &b| score[a].partial_cmp(&score[b]).unwrap());
    let mut rank = vec![0.0; score.len()];
    let mut i = 0;
    while i < o.len() {
        let mut j = i;
        while j + 1 < o.len() && score[o[j + 1]] == score[o[i]] { j += 1; }
        let r = (i + j) as f64 / 2.0 + 1.0;
        for k in i..=j { rank[o[k]] = r; }
        i = j + 1;
    }
    let rs: f64 = (0..lab.len()).filter(|&k| lab[k]).map(|k| rank[k]).sum();
    Some((rs - np as f64 * (np as f64 + 1.0) / 2.0) / (np as f64 * nn as f64))
}

/// orthonormal basis (row-major n x k) of columns of v restricted to rows `rows` (Gram-Schmidt)
fn basis_on(v: &[f64], k: usize, rows: &[usize]) -> Vec<Vec<f64>> {
    let mut b: Vec<Vec<f64>> = vec![];
    for j in 0..k {
        let mut c: Vec<f64> = rows.iter().map(|&i| v[i * k + j]).collect();
        for _ in 0..2 {
            for q in &b {
                let d: f64 = c.iter().zip(q).map(|(a, b)| a * b).sum();
                for (x, y) in c.iter_mut().zip(q) { *x -= d * y; }
            }
        }
        let nrm = c.iter().map(|x| x * x).sum::<f64>().sqrt();
        if nrm > 1e-10 { b.push(c.into_iter().map(|x| x / nrm).collect()); }
    }
    b
}
fn resid(x: &[f64], b: &[Vec<f64>]) -> Vec<f64> {
    let mut r = x.to_vec();
    for q in b {
        let d: f64 = r.iter().zip(q).map(|(a, b)| a * b).sum();
        for (v, y) in r.iter_mut().zip(q) { *v -= d * y; }
    }
    r
}

/// ridge W^T (|U| x k) from the reference Gram matrix C = Y^T Y (g x g):
/// (C_SS + alpha tr I) W = C_SU, so W^T rows = C_US
fn ridge(c: &[f64], g: usize, s: &[usize], u: &[usize], alpha: f64) -> Result<Vec<f64>> {
    let k = s.len();
    let mut a = vec![0.0; k * k];
    for (i, &si) in s.iter().enumerate() {
        for (j, &sj) in s.iter().enumerate() { a[i * k + j] = c[si * g + sj]; }
    }
    let tr: f64 = (0..k).map(|i| a[i * k + i]).sum::<f64>() / k as f64;
    for i in 0..k { a[i * k + i] += alpha * tr; }
    let mut bt = vec![0.0; u.len() * k];
    bt.par_chunks_mut(k).zip(u).for_each(|(r, &ui)| { for (j, &sj) in s.iter().enumerate() { r[j] = c[ui * g + sj]; } });
    la::posv(a, k, &mut bt, u.len())?;
    Ok(bt)
}

/// impute from layer `lin` panel of rows `rows` of zt: |rows| x |U|
fn impute(zt: &Z, rows: &[usize], s: &[usize], wt: &[f64], nu: usize, lin: usize) -> Vec<f64> {
    let x = zt.mat(lin, rows, s);
    la::gemm(false, true, rows.len(), nu, s.len(), &x, wt)
}

struct Panel {
    name: String,
    genes: Vec<usize>,
}

/// splitmix64 step
fn sm64(st: &mut u64) -> u64 {
    *st = st.wrapping_add(0x9E3779B97F4A7C15);
    let mut x = *st;
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}
const KM_DIM: usize = 160;
/// module panel: k-means++ (seed 1) + Lloyd on gene loadings V_160 * sqrt(eig); gene nearest each centroid
fn module_panel(gr: &Grams, g: usize, k: usize) -> Vec<usize> {
    let nk = *DEN_K.iter().max().unwrap();
    let d = KM_DIM;
    let x: Vec<f64> = (0..g).flat_map(|j| (0..d).map(move |c| (j, c))).map(|(j, c)| gr.vk[j * nk + c] * gr.ev[c].max(0.0).sqrt()).collect();
    let xn: Vec<f64> = (0..g).map(|j| x[j * d..(j + 1) * d].iter().map(|v| v * v).sum()).collect();
    let mut st = 1u64.wrapping_mul(0x9E3779B97F4A7C15);
    let unif = |st: &mut u64| (sm64(st) >> 11) as f64 / (1u64 << 53) as f64;
    // k-means++
    let mut cen: Vec<f64> = Vec::with_capacity(k * d);
    let first = (sm64(&mut st) % g as u64) as usize;
    cen.extend_from_slice(&x[first * d..(first + 1) * d]);
    let mut dmin: Vec<f64> = (0..g).map(|j| (0..d).map(|c| (x[j * d + c] - cen[c]).powi(2)).sum()).collect();
    for _ in 1..k {
        let tot: f64 = dmin.iter().sum();
        let mut t = unif(&mut st) * tot;
        let mut pick = g - 1;
        for (j, &v) in dmin.iter().enumerate() { if t < v { pick = j; break; } t -= v; }
        let c0 = cen.len();
        cen.extend_from_slice(&x[pick * d..(pick + 1) * d]);
        dmin.par_iter_mut().enumerate().for_each(|(j, m)| { let v: f64 = (0..d).map(|c| (x[j * d + c] - cen[c0 + c]).powi(2)).sum(); if v < *m { *m = v; } });
    }
    // Lloyd
    let mut asg = vec![usize::MAX; g];
    for _it in 0..300 {
        let cn: Vec<f64> = (0..k).map(|i| cen[i * d..(i + 1) * d].iter().map(|v| v * v).sum()).collect();
        let xc = la::gemm(false, true, g, k, d, &x, &cen);
        let new: Vec<usize> = (0..g).into_par_iter().map(|j| {
            let mut b = (f64::INFINITY, 0);
            for i in 0..k { let v = xn[j] - 2.0 * xc[j * k + i] + cn[i]; if v < b.0 { b = (v, i); } }
            b.1
        }).collect();
        if new == asg { break; }
        asg = new;
        let mut sum = vec![0.0; k * d];
        let mut cnt = vec![0usize; k];
        for j in 0..g { cnt[asg[j]] += 1; for c in 0..d { sum[asg[j] * d + c] += x[j * d + c]; } }
        for i in 0..k { if cnt[i] > 0 { for c in 0..d { cen[i * d + c] = sum[i * d + c] / cnt[i] as f64; } } }
    }
    // nearest gene to each non-empty centroid (distinct)
    let mut used = HashSet::new();
    let mut out = vec![];
    for i in 0..k {
        let mut b = (f64::INFINITY, usize::MAX);
        for j in 0..g { if asg[j] == i && !used.contains(&j) { let v: f64 = (0..d).map(|c| (x[j * d + c] - cen[i * d + c]).powi(2)).sum(); if v < b.0 { b = (v, j); } } }
        if b.1 != usize::MAX { used.insert(b.1); out.push(b.1); }
    }
    out
}

/// external panels: file fold\tsize\tgene; fold "all" applies to every fold
fn extra_panels(spec: &Option<String>, z: &Z, fold: i64) -> Result<Vec<Panel>> {
    let Some(spec) = spec else { return Ok(vec![]) };
    let gi: HashMap<&str, usize> = z.genes.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let mut out = vec![];
    for item in spec.split(',') {
        let (name, path) = item.split_once('=').ok_or_else(|| anyhow!("--extra name=path"))?;
        let mut by: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
        for l in std::fs::read_to_string(path)?.lines().skip(1) {
            let f: Vec<&str> = l.split('\t').collect();
            if f[0] == "all" || f[0].parse::<i64>().ok() == Some(fold) {
                if let Some(&j) = gi.get(f[2]) { by.entry(f[1].parse()?).or_default().push(j); }
            }
        }
        for (k, gs) in by { out.push(Panel { name: format!("{name}_{k}"), genes: gs }); }
    }
    Ok(out)
}

fn panels(z: &Z, refs: &[usize], sizes: &[usize], l1000: &HashSet<String>, nseeds: u64) -> Result<Vec<Panel>> {
    let g = z.genes.len();
    let y = z.mat(0, refs, &(0..g).collect::<Vec<_>>());
    let piv = la::qr_pivots(&y, refs.len(), g)?;
    let lm: Vec<usize> = (0..g).filter(|&j| l1000.contains(&z.genes[j])).collect();
    let mut out = vec![];
    for &k in sizes {
        out.push(Panel { name: format!("designed_{k}"), genes: piv[..k].to_vec() });
        for seed in 1..=nseeds {
            // deterministic shuffle (splitmix64)
            let mut st = seed.wrapping_mul(0x9E3779B97F4A7C15);
            let mut idx: Vec<usize> = (0..g).collect();
            for i in (1..g).rev() {
                st = st.wrapping_add(0x9E3779B97F4A7C15);
                let mut x = st;
                x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
                x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
                x ^= x >> 31;
                idx.swap(i, (x % (i as u64 + 1)) as usize);
            }
            out.push(Panel { name: format!("random{seed}_{k}"), genes: idx[..k].to_vec() });
        }
        if k == lm.len() {
            out.push(Panel { name: format!("l1000_{k}"), genes: lm.clone() });
        }
    }
    Ok(out)
}

/// top-KMAX reference neighbours per query by Pearson of panel vectors, with LS scale
const KS: [usize; 5] = [5, 10, 20, 50, 100];
const PWS: [f64; 3] = [1.0, 2.0, 4.0];
struct Nb { idx: Vec<usize>, r: Vec<f64>, sc: Vec<f64> }
fn neighbours(xq: &[f64], nq: usize, xp: &[f64], np: usize, ks: usize) -> Vec<Nb> {
    let stdz = |x: &[f64], n: usize| -> Vec<f64> {
        let mut o = x.to_vec();
        o.par_chunks_mut(ks).for_each(|r| {
            let m = r.iter().sum::<f64>() / ks as f64;
            r.iter_mut().for_each(|v| *v -= m);
            let nr = r.iter().map(|v| v * v).sum::<f64>().sqrt();
            if nr > 0.0 { r.iter_mut().for_each(|v| *v /= nr); }
        });
        let _ = n; o
    };
    let (qs, ps) = (stdz(xq, nq), stdz(xp, np));
    let c = la::gemm(false, true, nq, np, ks, &qs, &ps);
    let pn: Vec<f64> = (0..np).map(|i| xp[i * ks..(i + 1) * ks].iter().map(|v| v * v).sum()).collect();
    let kmax = KS[KS.len() - 1].min(np);
    (0..nq).into_par_iter().map(|q| {
        let row = &c[q * np..(q + 1) * np];
        let mut o: Vec<usize> = (0..np).collect();
        o.select_nth_unstable_by(kmax - 1, |&a, &b| row[b].partial_cmp(&row[a]).unwrap());
        o.truncate(kmax);
        o.sort_by(|&a, &b| row[b].partial_cmp(&row[a]).unwrap().then(a.cmp(&b)));
        let xqr = &xq[q * ks..(q + 1) * ks];
        let sc = o.iter().map(|&i| if pn[i] > 0.0 { xqr.iter().zip(&xp[i * ks..(i + 1) * ks]).map(|(a, b)| a * b).sum::<f64>() / pn[i] } else { 0.0 }).collect();
        Nb { r: o.iter().map(|&i| row[i]).collect(), idx: o, sc }
    }).collect()
}
/// base + weighted mean of scaled pool rows over the first k neighbours (weights max(r,0)^pw)
fn knn_pred(nb: &[Nb], base: Option<&[f64]>, pool: &[f64], nu: usize, k: usize, pw: f64) -> Vec<f64> {
    let mut out = vec![0.0; nb.len() * nu];
    out.par_chunks_mut(nu).zip(nb).enumerate().for_each(|(q, (o, n))| {
        if let Some(b) = base { o.copy_from_slice(&b[q * nu..(q + 1) * nu]); }
        let w: Vec<f64> = n.r.iter().take(k).map(|&r| r.max(0.0).powf(pw)).collect();
        let sw: f64 = w.iter().sum();
        if sw <= 0.0 { return; }
        for (t, &wi) in w.iter().enumerate() {
            let f = wi * n.sc[t] / sw;
            if f == 0.0 { continue; }
            let pr = &pool[n.idx[t] * nu..(n.idx[t] + 1) * nu];
            for (a, b) in o.iter_mut().zip(pr) { *a += f * b; }
        }
    });
    out
}
fn spa_de(v: &[f64], t: &[f64]) -> f64 {
    let de: Vec<usize> = (0..t.len()).filter(|&j| t[j].abs() > DE_Z).collect();
    if de.len() < 5 { return f64::NAN; }
    spearman(&de.iter().map(|&j| v[j].abs()).collect::<Vec<_>>(), &de.iter().map(|&j| t[j].abs()).collect::<Vec<_>>())
}
fn median(mut v: Vec<f64>) -> f64 {
    v.retain(|x| x.is_finite());
    if v.is_empty() { return f64::NAN; }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 { v[n / 2] } else { (v[n / 2 - 1] + v[n / 2]) / 2.0 }
}

/// evaluate one panel: choose alpha on inner holdout, INFO per gene, fit on all refs, score tests
#[allow(clippy::too_many_arguments)]
fn eval_panel(z: &Z, zt: &Z, gr: &Grams, refs: &[usize], itr: &[usize], inner_ho: &[usize], tests: &[usize], p: &Panel, fold: i64, common: Option<&[usize]>, wr: &mut impl Write, wg: &mut impl Write, wd: &mut impl Write, wk: &mut impl Write) -> Result<()> {
    let (v, vk) = (&gr.v, &gr.vk);
    let g = z.genes.len();
    let inpanel: HashSet<usize> = p.genes.iter().copied().collect();
    // evaluation genes: complement of this panel, or a common set shared by all panels
    let u: Vec<usize> = match common {
        Some(e) => { assert!(e.iter().all(|j| !inpanel.contains(j))); e.to_vec() }
        None => (0..g).filter(|j| !inpanel.contains(j)).collect(),
    };
    let nu = u.len();
    // inner holdout: alpha by mean per-perturbation r (half0 panel -> half1 truth)
    let truth_ho = z.mat(2, inner_ho, &u);
    let mut best = (f64::NEG_INFINITY, ALPHAS[0], vec![], vec![]);
    for &al in &ALPHAS {
        let wt = ridge(&gr.c_in, g, &p.genes, &u, al)?;
        let im = impute(z, inner_ho, &p.genes, &wt, nu, 1);
        let r: f64 = (0..inner_ho.len()).map(|i| pearson(&im[i * nu..(i + 1) * nu], &truth_ho[i * nu..(i + 1) * nu])).filter(|x| x.is_finite()).sum::<f64>() / inner_ho.len() as f64;
        if r > best.0 { best = (r, al, im, wt); }
    }
    let (_, alpha, im_ho, wt_in) = best;
    // INFO_g: per-gene r across inner holdout perturbations
    let col = |m: &[f64], n: usize, j: usize| -> Vec<f64> { (0..n).map(|i| m[i * nu + j]).collect() };
    let info: Vec<f64> = (0..nu).into_par_iter().map(|j| pearson(&col(&im_ho, inner_ho.len(), j), &col(&truth_ho, inner_ho.len(), j))).collect();
    // per-gene calibration slope b and imputation error variance e (inner holdout only)
    let nh = inner_ho.len();
    let (bcal, ecal): (Vec<f64>, Vec<f64>) = (0..nu).into_par_iter().map(|j| {
        let (x, y) = (col(&im_ho, nh, j), col(&truth_ho, nh, j));
        let sxx: f64 = x.iter().map(|a| a * a).sum();
        let b = if sxx > 0.0 { x.iter().zip(&y).map(|(a, c)| a * c).sum::<f64>() / sxx } else { 0.0 };
        let mse = x.iter().zip(&y).map(|(a, c)| (c - b * a).powi(2)).sum::<f64>() / nh as f64;
        (b, (mse - 1.0).max(0.05))
    }).unzip();
    // final model
    let wt = ridge(&gr.c, g, &p.genes, &u, alpha)?;
    let im = impute(zt, tests, &p.genes, &wt, nu, 1);
    let tr = zt.mat(2, tests, &u);
    let dr = zt.mat(1, tests, &u);
    let real: Vec<f64> = (0..nu).into_par_iter().map(|j| pearson(&col(&im, tests.len(), j), &col(&tr, tests.len(), j))).collect();
    for j in 0..nu {
        writeln!(wg, "{fold}\t{}\t{}\t{:.6}\t{:.6}", p.name, z.genes[u[j]], info[j], real[j])?;
    }
    // OOR: fraction of half-0 panel vector inside top-K reference panel subspace
    let k = p.genes.len();
    let kk = 160.min(k / 2);
    let mut css = vec![0.0; k * k];
    for (i, &si) in p.genes.iter().enumerate() { for (j, &sj) in p.genes.iter().enumerate() { css[i * k + j] = gr.c[si * g + sj]; } }
    let (pv, _) = la::top_eig(&css, k, kk)?;
    // denoised full-transcriptome half-0 profile: projection on top-K reference programs
    let h0full = zt.mat(1, tests, &(0..g).collect::<Vec<_>>());
    let nk = DEN_K.iter().max().copied().unwrap();
    let coef = la::gemm(false, false, tests.len(), nk, g, &h0full, vk); // nT x nk
    // denoised profiles for each K via BLAS: (nT x K) (K x g); vk row-major g x nk -> take first K cols
    let denf: Vec<Vec<f64>> = DEN_K.iter().map(|&kd| {
        let ck: Vec<f64> = (0..tests.len()).flat_map(|i| coef[i * nk..i * nk + kd].to_vec()).collect();
        let vt: Vec<f64> = (0..g).flat_map(|j| vk[j * nk..j * nk + kd].to_vec()).collect(); // g x kd
        la::gemm(false, true, tests.len(), g, kd, &ck, &vt)
    }).collect();
    let pb: Vec<Vec<f64>> = (0..kk).map(|c| (0..k).map(|i| pv[i * kk + c]).collect()).collect();
    let basis = basis_on(v, NPC, &u);
    let rows: Vec<String> = (0..tests.len())
        .into_par_iter()
        .map(|i| {
            let (yi, ti, di) = (&im[i * nu..(i + 1) * nu], &tr[i * nu..(i + 1) * nu], &dr[i * nu..(i + 1) * nu]);
            let (ry, rt, rd) = (resid(yi, &basis), resid(ti, &basis), resid(di, &basis));
            let lab: Vec<bool> = ti.iter().map(|t| t.abs() > 4.0).collect();
            let a_imp = auc(&yi.iter().map(|x| x.abs()).collect::<Vec<_>>(), &lab).unwrap_or(f64::NAN);
            let a_dir = auc(&di.iter().map(|x| x.abs()).collect::<Vec<_>>(), &lab).unwrap_or(f64::NAN);
            let x: Vec<f64> = p.genes.iter().map(|&j| zt.at(1, tests[i], j)).collect();
            let nx: f64 = x.iter().map(|v| v * v).sum();
            let inside: f64 = pb.iter().map(|q| x.iter().zip(q).map(|(a, b)| a * b).sum::<f64>().powi(2)).sum();
            let den: Vec<String> = (0..DEN_K.len()).map(|q| {
                let d: Vec<f64> = u.iter().map(|&j| denf[q][i * g + j]).collect();
                format!("{:.6}", pearson(&resid(&d, &basis), &rt))
            }).collect();
            format!("{fold}\t{}\t{}\t{alpha}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}", p.name, zt.perts[tests[i]], pearson(di, ti), pearson(yi, ti), pearson(&rd, &rt), pearson(&ry, &rt), a_dir, a_imp, 1.0 - inside / nx, den.join("\t"))
        })
        .collect();
    for r in rows { writeln!(wr, "{r}")?; }
    // DE call, direction and relative size (truth = half 1, direct = half 0 measured)
    let derows: Vec<String> = (0..tests.len()).into_par_iter().map(|i| {
        let (yi, ti, di) = (&im[i * nu..(i + 1) * nu], &tr[i * nu..(i + 1) * nu], &dr[i * nu..(i + 1) * nu]);
        let ci: Vec<f64> = (0..nu).map(|j| bcal[j] * yi[j]).collect();
        let de: Vec<usize> = (0..nu).filter(|&j| ti[j].abs() > DE_Z).collect();
        // calls: (value, called) for direct, raw imputed, calibrated imputed
        let calls = |v: &[f64], f: &dyn Fn(usize) -> bool| -> (usize, usize, usize) {
            let (mut nc, mut tp, mut disc) = (0, 0, 0);
            for j in 0..nu {
                if !f(j) { continue; }
                nc += 1;
                let same = v[j].signum() == ti[j].signum();
                if !same { disc += 1; }
                if same && ti[j].abs() > DE_Z { tp += 1; }
            }
            (nc, tp, disc)
        };
        let m = calls(di, &|j| di[j].abs() > DE_Z);
        let r = calls(yi, &|j| yi[j].abs() > DE_Z);
        let cc = |j: usize| phi(-ci[j].abs() / ecal[j].sqrt()) < LFSR;
        let c = calls(&ci, &cc);
        // matched: top-n genes by |c|, n = number of direct calls for this perturbation
        let mut o: Vec<usize> = (0..nu).collect();
        o.sort_by(|&a, &b| ci[b].abs().partial_cmp(&ci[a].abs()).unwrap());
        let top: HashSet<usize> = o[..m.0].iter().copied().collect();
        let k = calls(&ci, &|j| top.contains(&j));
        let exp_disc: f64 = (0..nu).filter(|&j| cc(j)).map(|j| phi(-ci[j].abs() / (ecal[j] + 1.0).sqrt())).sum();
        let sums = |v: &[f64]| -> (f64, f64) { (v.iter().zip(ti).map(|(a, b)| a * b).sum(), v.iter().map(|a| a * a).sum()) };
        let (sm, sr, sc) = (sums(di), sums(yi), sums(&ci));
        let sp = |v: &[f64], ab: bool| -> f64 {
            if de.len() < 5 { return f64::NAN; }
            let f = |x: f64| if ab { x.abs() } else { x };
            spearman(&de.iter().map(|&j| f(v[j])).collect::<Vec<_>>(), &de.iter().map(|&j| f(ti[j])).collect::<Vec<_>>())
        };
        format!("{fold}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}",
            p.name, zt.perts[tests[i]], de.len(), m.0, m.1, m.2, r.0, r.1, r.2, c.0, c.1, c.2, exp_disc,
            sm.0, sm.1, sr.0, sr.1, sc.0, sc.1, sp(di, false), sp(yi, false), sp(&ci, false), sp(di, true), sp(yi, true), sp(&ci, true), k.1, k.2)
    }).collect();
    for r in derows { writeln!(wd, "{r}")?; }
    // neighbour copying (knn) and ridge + neighbour residuals (hybrid); tune on inner holdout
    let ks = p.genes.len();
    let tune = |pool_rows: &[usize], wtx: &[f64], qz: &Z, qrows: &[usize], base: &[f64]| -> (Vec<Nb>, Vec<f64>, Vec<f64>) {
        let xp = z.mat(0, pool_rows, &p.genes);
        let yp = z.mat(0, pool_rows, &u);
        let mut rp = la::gemm(false, true, pool_rows.len(), nu, ks, &xp, wtx);
        rp.par_iter_mut().zip(&yp).for_each(|(a, b)| *a = b - *a);
        let xq = qz.mat(1, qrows, &p.genes);
        let _ = base;
        (neighbours(&xq, qrows.len(), &xp, pool_rows.len(), ks), yp, rp)
    };
    let (nb_in, yp_in, rp_in) = tune(itr, &wt_in, z, inner_ho, &im_ho);
    let mut bestk = [(f64::NEG_INFINITY, 0usize, 0.0f64); 2];
    for &kk in &KS {
        for &pw in &PWS {
            for mode in 0..2 {
                let pr = if mode == 0 { knn_pred(&nb_in, None, &yp_in, nu, kk, pw) } else { knn_pred(&nb_in, Some(&im_ho), &rp_in, nu, kk, pw) };
                let obj = median((0..nh).into_par_iter().map(|i| spa_de(&pr[i * nu..(i + 1) * nu], &truth_ho[i * nu..(i + 1) * nu])).collect());
                if obj > bestk[mode].0 { bestk[mode] = (obj, kk, pw); }
            }
        }
    }
    drop((nb_in, yp_in, rp_in));
    let (nb_t, yp_t, rp_t) = tune(refs, &wt, zt, tests, &im);
    let pk = knn_pred(&nb_t, None, &yp_t, nu, bestk[0].1, bestk[0].2);
    let ph = knn_pred(&nb_t, Some(&im), &rp_t, nu, bestk[1].1, bestk[1].2);
    drop((nb_t, yp_t, rp_t));
    let krows: Vec<String> = (0..tests.len()).into_par_iter().map(|i| {
        let (ti, di) = (&tr[i * nu..(i + 1) * nu], &dr[i * nu..(i + 1) * nu]);
        let nm = di.iter().filter(|x| x.abs() > DE_Z).count();
        let rt = resid(ti, &basis);
        let ntde = ti.iter().filter(|x| x.abs() > DE_Z).count();
        let mut f = format!("{fold}\t{}\t{}\t{ntde}\t{nm}", p.name, zt.perts[tests[i]]);
        for v in [&im[i * nu..(i + 1) * nu], &pk[i * nu..(i + 1) * nu], &ph[i * nu..(i + 1) * nu]] {
            let mut o: Vec<usize> = (0..nu).collect();
            o.sort_by(|&a, &b| v[b].abs().partial_cmp(&v[a].abs()).unwrap());
            let (mut tp, mut disc) = (0, 0);
            for &j in &o[..nm] { let same = v[j].signum() == ti[j].signum(); if !same { disc += 1; } if same && ti[j].abs() > DE_Z { tp += 1; } }
            let de: Vec<usize> = (0..nu).filter(|&j| ti[j].abs() > DE_Z).collect();
            let sp = if de.len() >= 5 { spearman(&de.iter().map(|&j| v[j]).collect::<Vec<_>>(), &de.iter().map(|&j| ti[j]).collect::<Vec<_>>()) } else { f64::NAN };
            f += &format!("\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{tp}\t{disc}", pearson(v, ti), pearson(&resid(v, &basis), &rt), sp, spa_de(v, ti));
        }
        f + &format!("\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}", bestk[0].1, bestk[0].2, bestk[1].1, bestk[1].2, bestk[0].0, bestk[1].0)
    }).collect();
    for r in krows { writeln!(wk, "{r}")?; }
    Ok(())
}

/// native baselines: hvgresp (variance of reference z) and module (k-means hubs)
fn native_baselines(z: &Z, refs: &[usize], gr: &Grams, sizes: &[usize]) -> Vec<Panel> {
    let g = z.genes.len();
    let var: Vec<f64> = (0..g).into_par_iter().map(|j| {
        let v: Vec<f64> = refs.iter().map(|&i| z.at(0, i, j)).collect();
        let m = v.iter().sum::<f64>() / v.len() as f64;
        v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len() - 1) as f64
    }).collect();
    let mut o: Vec<usize> = (0..g).collect();
    o.sort_by(|&a, &b| var[b].partial_cmp(&var[a]).unwrap().then(a.cmp(&b)));
    let mut out = vec![];
    for &k in sizes {
        out.push(Panel { name: format!("hvgresp_{k}"), genes: o[..k].to_vec() });
        out.push(Panel { name: format!("module_{k}"), genes: module_panel(gr, g, k) });
        out.push(Panel { name: format!("greedy_{k}"), genes: sim::greedy(gr, refs.len(), g, k, &vec![1.0; g]) });
    }
    out
}

/// genes outside every panel (same evaluation set for all panels of a fold)
fn common_set(ps: &[Panel], g: usize, on: bool) -> Option<Vec<usize>> {
    if !on { return None; }
    let used: HashSet<usize> = ps.iter().flat_map(|p| p.genes.iter().copied()).collect();
    Some((0..g).filter(|j| !used.contains(j)).collect())
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let cmd = a.get(1).cloned().unwrap_or_default();
    if cmd == "pbdump" { return sim::pbdump(&a); }
    if cmd == "upgsum" { return upg::summary(&a); }
    if cmd == "simsum" { return upg::simsum(&a); }
    if cmd == "desum" { return de::summary(&a); }
    if cmd == "cvsum" { return upg::cvsum(&a); }
    if cmd == "subset" { return multi::subset(&a); }
    if cmd == "concat" { return multi::concat(&a); }
    if cmd == "multi" { return multi::run(&a); }
    if cmd == "multisum" { return multi::summary(&a); }
    if cmd == "multitest" { return multi::multitest(&a); }
    if cmd == "vcm" { return vcm::run(&a); }
    if cmd == "vcmsum" { return vcm::summary(&a); }
    if cmd == "causal" { return causal::run(&a); }
    if cmd == "explore" { return explore::run(&a); }
    if cmd == "mtx2cache" { return real::mtx2cache(&a); }
    if cmd == "realsum" { return real::summary(&a); }
    if cmd == "cellsum" { return real::cellsum(&a); }
    let z = Z::load(&req(&a, "--z")?)?;
    let l1000: HashSet<String> = std::fs::read_to_string(req(&a, "--l1000")?)?.lines().map(String::from).collect();
    let n_lm = z.genes.iter().filter(|g| l1000.contains(*g)).count();
    let sizes: Vec<usize> = arg(&a, "--sizes").map(|s| s.split(',').map(|x| if x == "lm" { n_lm } else { x.parse().unwrap() }).collect()).unwrap_or(vec![100, 200, 500, n_lm]);
    let out = req(&a, "--out")?;
    let common = a.contains(&"--common".to_string());
    let base = a.contains(&"--baselines".to_string());
    let extra = arg(&a, "--extra");
    let nseeds: u64 = arg(&a, "--seeds").map(|v| v.parse().unwrap()).unwrap_or(3);
    let mut wr = BufWriter::new(std::fs::File::create(format!("{out}.perts.tsv"))?);
    let mut wg = BufWriter::new(std::fs::File::create(format!("{out}.genes.tsv"))?);
    let mut wk = BufWriter::new(std::fs::File::create(format!("{out}.knn.tsv"))?);
    writeln!(wk, "fold\tpanel\tpert\tn_de\tm_call\t{}\tk_knn\tp_knn\tk_hyb\tp_hyb\tobj_knn\tobj_hyb", ["ridge", "knn", "hyb"].iter().map(|m| ["r", "rperp", "sp", "spa", "ktp", "kdisc"].iter().map(|c| format!("{m}_{c}")).collect::<Vec<_>>().join("\t")).collect::<Vec<_>>().join("\t"))?;
    let mut wd = BufWriter::new(std::fs::File::create(format!("{out}.de.tsv"))?);
    writeln!(wd, "fold\tpanel\tpert\tn_de\tm_call\tm_tp\tm_disc\tr_call\tr_tp\tr_disc\tc_call\tc_tp\tc_disc\tc_exp_disc\tm_sxy\tm_sxx\tr_sxy\tr_sxx\tc_sxy\tc_sxx\tsp_m\tsp_r\tsp_c\tspa_m\tspa_r\tspa_c\tk_tp\tk_disc")?;
    writeln!(wr, "fold\tpanel\tpert\talpha\tr_ceil\tr_imp\tr_ceil_perp\tr_imp_perp\tauc_ceil\tauc_imp\toor\t{}", DEN_K.iter().map(|k| format!("r_den{k}_perp")).collect::<Vec<_>>().join("\t"))?;
    writeln!(wg, "fold\tpanel\tgene\tinfo\treal")?;
    let g = z.genes.len();
    let pi: HashMap<&str, usize> = z.perts.iter().enumerate().map(|(i, p)| (p.as_str(), i)).collect();
    // fold, group per perturbation
    let mut fold = vec![-1i64; z.perts.len()];
    let mut group = vec![-1i64; z.perts.len()];
    let mut excl: HashSet<(usize, i64)> = HashSet::new();
    if cmd == "real" || cmd == "realcells" || cmd == "cv" || cmd == "sim" || cmd == "upg" || cmd == "de" || cmd == "example" {
        for l in std::fs::read_to_string(req(&a, "--folds")?)?.lines().skip(1) {
            let f: Vec<&str> = l.split('\t').collect();
            if let Some(&i) = pi.get(f[0]) { group[i] = f[1].parse()?; fold[i] = f[2].parse()?; }
        }
        for l in std::fs::read_to_string(req(&a, "--hubs")?)?.lines().skip(1) {
            let f: Vec<&str> = l.split('\t').collect();
            if let Some(&i) = pi.get(f[0]) { excl.insert((i, f[1].parse()?)); }
        }
    }
    let inner_split = |refs: &[usize]| -> (Vec<usize>, Vec<usize>) {
        // 10% of groups (hash of group id; hubs keyed by index) as inner holdout
        let key = |i: usize| -> u64 { let k = if group[i] >= 0 { group[i] as u64 } else { 1_000_000 + i as u64 }; (k.wrapping_mul(0x9E3779B97F4A7C15) >> 32) % 10 };
        (refs.iter().copied().filter(|&i| key(i) != 0).collect(), refs.iter().copied().filter(|&i| key(i) == 0).collect())
    };
    let t0 = std::time::Instant::now();
    match cmd.as_str() {
        "cv" => {
            let runs: Vec<i64> = arg(&a, "--folds-run").map(|s| s.split(',').map(|x| x.parse().unwrap()).collect()).unwrap_or((0..5).collect());
            for f in runs {
                let refs: Vec<usize> = (0..z.perts.len()).filter(|&i| fold[i] != f && !excl.contains(&(i, f))).collect();
                let tests: Vec<usize> = (0..z.perts.len()).filter(|&i| fold[i] == f).collect();
                let (itr, iho) = inner_split(&refs);
                let gr = grams(&z, &refs, &itr)?;
                let mut ps = panels(&z, &refs, &sizes, &l1000, nseeds)?;
                if base { ps.extend(native_baselines(&z, &refs, &gr, &sizes)); ps.extend(extra_panels(&extra, &z, f)?); }
                let ev = common_set(&ps, g, common);
                for p in &ps {
                    eval_panel(&z, &z, &gr, &refs, &itr, &iho, &tests, p, f, ev.as_deref(), &mut wr, &mut wg, &mut wd, &mut wk)?;
                }
                if let Some(e) = &ev { eprintln!("fold {f}: common evaluation genes {}", e.len()); }
                eprintln!("fold {f}: refs {} tests {} panels {} at {:.1}s", refs.len(), tests.len(), ps.len(), t0.elapsed().as_secs_f64());
            }
        }
        "transfer" => {
            let zt = Z::load(&req(&a, "--zt")?)?;
            if zt.genes != z.genes { bail!("gene order differs between --z and --zt") }
            let refs: Vec<usize> = (0..z.perts.len()).collect();
            let tests: Vec<usize> = (0..zt.perts.len()).collect();
            let (itr, iho) = inner_split(&refs);
            let gr = grams(&z, &refs, &itr)?;
            let mut ps = panels(&z, &refs, &sizes, &l1000, nseeds)?;
            if base { ps.extend(native_baselines(&z, &refs, &gr, &sizes)); ps.extend(extra_panels(&extra, &z, -9)?); }
            let ev = common_set(&ps, g, common);
            for p in &ps {
                eval_panel(&z, &zt, &gr, &refs, &itr, &iho, &tests, p, -9, ev.as_deref(), &mut wr, &mut wg, &mut wd, &mut wk)?;
            }
            if let Some(e) = &ev { eprintln!("transfer: common evaluation genes {}", e.len()); }
            eprintln!("transfer: refs {} tests {} at {:.1}s", refs.len(), tests.len(), t0.elapsed().as_secs_f64());
        }
        "greedytest" => sim::greedytest(&a, &z)?,
        "de" => match arg(&a, "--zt") {
            Some(p) => { let zt = Z::load(&p)?; if zt.genes != z.genes { bail!("gene order differs") } de::run(&a, &z, &zt, &fold, &group, &excl, &out)? }
            None => de::run(&a, &z, &z, &fold, &group, &excl, &out)?,
        },
        "folds" => folds::run(&a, &z)?,
        "genevar" => multi::genevar(&a, &z)?,
        "example" => {
            // gene-level profiles of chosen test perturbations: greedy 500 panel and ridge from the fold reference
            let f: i64 = req(&a, "--fold")?.parse()?;
            let alpha: f64 = arg(&a, "--alpha").map(|v| v.parse().unwrap()).unwrap_or(1.0);
            let k: usize = arg(&a, "--panel-size").map(|v| v.parse().unwrap()).unwrap_or(500);
            let refs: Vec<usize> = (0..z.perts.len()).filter(|&i| fold[i] != f && !excl.contains(&(i, f))).collect();
            let gr = grams_c(&z, &refs, &refs[..1]);
            let s = sim::greedy(&gr, refs.len(), g, k, &vec![1.0; g]);
            let inp: HashSet<usize> = s.iter().copied().collect();
            let u: Vec<usize> = (0..g).filter(|j| !inp.contains(j)).collect();
            let w = ridge(&gr.c, g, &s, &u, alpha)?;
            let mut wx = BufWriter::new(std::fs::File::create(format!("{out}.example.tsv"))?);
            writeln!(wx, "pert\tgene\tin_panel\timputed\thalf0\thalf1\tall")?;
            for name in req(&a, "--perts")?.split(',') {
                let Some(&i) = pi.get(name) else { eprintln!("example: {name} not found"); continue };
                if fold[i] != f { eprintln!("example: {name} is in fold {}", fold[i]); continue }
                let x: Vec<f64> = s.iter().map(|&j| z.at(1, i, j)).collect();
                for (q, &j) in u.iter().enumerate() {
                    let y: f64 = (0..s.len()).map(|r| w[q * s.len() + r] * x[r]).sum();
                    writeln!(wx, "{name}\t{}\t0\t{y:.4}\t{:.4}\t{:.4}\t{:.4}", z.genes[j], z.at(1, i, j), z.at(2, i, j), z.at(0, i, j))?;
                }
                for &j in &s { writeln!(wx, "{name}\t{}\t1\t{:.4}\t{:.4}\t{:.4}\t{:.4}", z.genes[j], z.at(1, i, j), z.at(1, i, j), z.at(2, i, j), z.at(0, i, j))?; }
            }
        }
        "mp" => multi::mpcount(&a, &z)?,
        "upg" => match arg(&a, "--zt") {
            Some(p) => { let zt = Z::load(&p)?; if zt.genes != z.genes { bail!("gene order differs") } upg::run(&a, &z, &zt, &fold, &group, &excl, &out)? }
            None => upg::run(&a, &z, &z, &fold, &group, &excl, &out)?,
        },
        "realcells" => real::cells(&a, &z, &fold, &group, &excl, &l1000, &out)?,
        "real" => real::run(&a, &z, &fold, &group, &excl, &l1000, &out)?,
        "sim" => sim::run(&a, &z, &fold, &group, &excl, &l1000, &sizes, &out)?,
        _ => bail!("usage: panelimp cv|transfer|sim|pbdump ..."),
    }
    Ok(())
}
