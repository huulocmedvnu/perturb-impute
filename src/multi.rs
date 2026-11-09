//! Multi-atlas reference (Experiment 2). References: normalised second-moment matrices
//! Sigma_k = C_k / n_k of genome-wide atlases on a common gene set. Panels: greedy MMSE on one
//! atlas or on the equal-weight pan matrix. Imputation: ridge from a single atlas, the equal pan
//! matrix, or a pan matrix weighted by non-negative least squares fit of the target screen's
//! measured panel covariance. Inter-atlas disagreement: spread of the single-atlas imputations.
use crate::*;

/// z prefix restricted to a gene list (in that order); genes absent from the prefix are an error
pub fn subset(a: &[String]) -> Result<()> {
    let z = Z::load(&req(a, "--z")?)?;
    let want: Vec<String> = std::fs::read_to_string(req(a, "--genes")?)?.lines().map(String::from).collect();
    let gi: HashMap<&str, usize> = z.genes.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let cols: Vec<usize> = want.iter().map(|s| gi.get(s.as_str()).copied().ok_or_else(|| anyhow!("gene {s} missing"))).collect::<Result<_>>()?;
    let (p, g) = (z.perts.len(), z.genes.len());
    let mut out = vec![0f32; 3 * p * cols.len()];
    for l in 0..3 { out[l * p * cols.len()..(l + 1) * p * cols.len()].par_chunks_mut(cols.len()).enumerate().for_each(|(i, r)| for (k, &j) in cols.iter().enumerate() { r[k] = z.z[(l * p + i) * g + j]; }); }
    let o = req(a, "--out")?;
    std::fs::write(format!("{o}.z.f32"), bytemuck::cast_slice(&out))?;
    let src = req(a, "--z")?;
    std::fs::copy(format!("{src}.meta.tsv"), format!("{o}.meta.tsv"))?;
    let mut w = BufWriter::new(std::fs::File::create(format!("{o}.genes.tsv"))?);
    writeln!(w, "gene")?;
    for s in &want { writeln!(w, "{s}")?; }
    eprintln!("subset: {p} perturbations x {} genes", cols.len());
    Ok(())
}

/// non-negative least squares for a handful of variables by enumerating active sets
pub fn nnls_small(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let m = a.len();
    let dot = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(p, q)| p * q).sum::<f64>();
    let mut best = (f64::INFINITY, vec![0.0; m]);
    for mask in 1u32..(1 << m) {
        let idx: Vec<usize> = (0..m).filter(|&i| mask >> i & 1 == 1).collect();
        let k = idx.len();
        let mut g = vec![0.0; k * k];
        let mut r: Vec<f64> = idx.iter().map(|&i| dot(&a[i], b)).collect();
        for (p, &i) in idx.iter().enumerate() { for (q, &j) in idx.iter().enumerate() { g[p * k + q] = dot(&a[i], &a[j]); } }
        for p in 0..k { g[p * k + p] *= 1.0 + 1e-12; }
        if la::posv(g, k, &mut r, 1).is_err() || r.iter().any(|&v| v < 0.0) { continue; }
        let mut x = vec![0.0; m];
        for (p, &i) in idx.iter().enumerate() { x[i] = r[p]; }
        let res: f64 = (0..b.len()).map(|t| { let f: f64 = (0..m).map(|i| x[i] * a[i][t]).sum(); (b[t] - f).powi(2) }).sum();
        if res < best.0 { best = (res, x); }
    }
    best.1
}

struct Atlas { name: String, z: Z, sig: Vec<f64>, s2: f64 }
fn atlas(name: &str, pre: &str) -> Result<Atlas> {
    let z = Z::load(pre)?;
    let (p, g) = (z.perts.len(), z.genes.len());
    let all: Vec<usize> = (0..p).collect();
    let mut sig = la::gram(&z.mat(0, &all, &(0..g).collect::<Vec<_>>()), p, g);
    sig.par_iter_mut().for_each(|v| *v /= p as f64);
    let s2 = all.par_iter().map(|&i| (0..g).map(|j| (z.at(1, i, j) - z.at(2, i, j)).powi(2)).sum::<f64>()).sum::<f64>() / (2.0 * (p * g) as f64);
    eprintln!("atlas {name}: {p} perturbations, noise variance {s2:.3}");
    Ok(Atlas { name: name.into(), z, sig, s2 })
}
fn gr_of(c: Vec<f64>) -> Grams { Grams { c, c_in: vec![], v: vec![], vk: vec![], ev: vec![] } }

pub fn run(a: &[String]) -> Result<()> {
    let t0 = std::time::Instant::now();
    let refs: Vec<Atlas> = req(a, "--refs")?.split(',').map(|kv| { let (k, v) = kv.split_once('=').unwrap(); atlas(k, v) }).collect::<Result<_>>()?;
    let targets: Vec<(String, String)> = req(a, "--targets")?.split(',').map(|kv| { let (k, v) = kv.split_once('=').unwrap(); (k.to_string(), v.to_string()) }).collect();
    let g = refs[0].z.genes.len();
    for r in &refs { if r.z.genes != refs[0].z.genes { bail!("gene order differs for {}", r.name) } }
    let k = 500;
    let alpha = 1.0;
    let mut w = BufWriter::new(std::fs::File::create(format!("{}.multi.tsv", req(a, "--out")?))?);
    writeln!(w, "target\tpanel\tref\tpert\tr\trperp\tden\tr_ceil\trperp_ceil\tn_de\tncall\ttp\tdisc\tspa\tsxy\tsxx\toor\tdisagree")?;
    let mut wl = BufWriter::new(std::fs::File::create(format!("{}.multi.log", req(a, "--out")?))?);
    for (tname, tpre) in &targets {
        let zt = Z::load(tpre)?;
        if zt.genes != refs[0].z.genes { bail!("target {tname} gene order differs") }
        // references available for this target: all atlases except the target itself
        let rs: Vec<&Atlas> = refs.iter().filter(|r| &r.name != tname).collect();
        let nr = rs.len() as f64;
        let mut pan = vec![0.0; g * g];
        for r in &rs { pan.par_iter_mut().zip(&r.sig).for_each(|(p, s)| *p += s / nr); }
        let pan_s2 = rs.iter().map(|r| r.s2).sum::<f64>() / nr;
        // evaluation frame from the pan matrix: 10 global programs, 160 denoising programs
        // --frame target: programs from the half-0 profiles of the target screen, neutral to every reference
        let (vk, _) = if arg(a, "--frame").as_deref() == Some("target") {
            let all_t: Vec<usize> = (0..zt.perts.len()).collect();
            la::top_eig(&la::gram(&zt.mat(1, &all_t, &(0..g).collect::<Vec<_>>()), all_t.len(), g), g, 160)?
        } else { la::top_eig(&pan, g, 160)? };
        let mut panels: Vec<(String, Vec<usize>)> = rs.iter().map(|r| (format!("greedy_{}", r.name), sim::greedy(&gr_of(r.sig.iter().map(|v| v / (1.0f64).max(r.s2)).collect()), 1, g, k, &vec![1.0; g]))).collect();
        panels.push(("greedy_pan".into(), sim::greedy(&gr_of(pan.iter().map(|v| v / pan_s2).collect()), 1, g, k, &vec![1.0; g])));
        let tests: Vec<usize> = (0..zt.perts.len()).collect();
        let nt = tests.len();
        for (pname, s) in &panels {
            let inp: HashSet<usize> = s.iter().copied().collect();
            let u: Vec<usize> = (0..g).filter(|j| !inp.contains(j)).collect();
            let nu = u.len();
            let v10: Vec<f64> = (0..g).flat_map(|j| vk[j * 160..j * 160 + NPC].to_vec()).collect();
            let basis = basis_on(&v10, NPC, &u);
            let h0 = zt.mat(1, &tests, &(0..g).collect::<Vec<_>>());
            let coef = la::gemm(false, false, nt, 160, g, &h0, &vk);
            let vu: Vec<f64> = u.iter().flat_map(|&j| vk[j * 160..(j + 1) * 160].to_vec()).collect();
            let den = la::gemm(false, true, nt, nu, 160, &coef, &vu);
            let (xt, tr, dr) = (zt.mat(1, &tests, s), zt.mat(2, &tests, &u), zt.mat(1, &tests, &u));
            // unsupervised mixture weights from the measured panel covariance of the target
            let st: Vec<f64> = { let mut m = la::gemm(true, false, k, k, nt, &xt, &xt); m.iter_mut().for_each(|v| *v /= nt as f64); m };
            let mut design: Vec<Vec<f64>> = rs.iter().map(|r| s.iter().flat_map(|&i| s.iter().map(move |&j| (i, j))).map(|(i, j)| r.sig[i * g + j]).collect()).collect();
            design.push((0..k * k).map(|q| if q / k == q % k { 1.0 } else { 0.0 }).collect());
            let wts = nnls_small(&design, &st);
            let ws: f64 = wts[..rs.len()].iter().sum();
            let mut panw = vec![0.0; g * g];
            for (r, &wk) in rs.iter().zip(&wts) { if wk > 0.0 { panw.par_iter_mut().zip(&r.sig).for_each(|(p, x)| *p += wk / ws * x); } }
            writeln!(wl, "{tname} {pname}: mixture weights {}", rs.iter().zip(&wts).map(|(r, w)| format!("{} {:.3}", r.name, w / ws)).collect::<Vec<_>>().join(", "))?;
            // OOR on the pan panel block
            let mut css = vec![0.0; k * k];
            for (i, &si) in s.iter().enumerate() { for (j, &sj) in s.iter().enumerate() { css[i * k + j] = pan[si * g + sj]; } }
            let (pv, _) = la::top_eig(&css, k, 160)?;
            let proj = la::gemm(false, false, nt, 160, k, &xt, &pv);
            let oor: Vec<f64> = (0..nt).map(|i| { let nx: f64 = xt[i * k..(i + 1) * k].iter().map(|v| v * v).sum(); 1.0 - proj[i * 160..(i + 1) * 160].iter().map(|v| v * v).sum::<f64>() / nx }).collect();
            let mut preds: Vec<(String, Vec<f64>)> = vec![];
            for r in &rs { preds.push((r.name.clone(), la::gemm(false, true, nt, nu, k, &xt, &ridge(&r.sig, g, s, &u, alpha)?))); }
            preds.push(("pan_equal".into(), la::gemm(false, true, nt, nu, k, &xt, &ridge(&pan, g, s, &u, alpha)?)));
            preds.push(("pan_weighted".into(), la::gemm(false, true, nt, nu, k, &xt, &ridge(&panw, g, s, &u, alpha)?)));
            // inter-atlas disagreement: variance across single-atlas predictions / mean squared prediction
            let ns = rs.len();
            let dis: Vec<f64> = (0..nt).map(|i| {
                let (mut v, mut m2) = (0.0, 0.0);
                for j in 0..nu { let xs: Vec<f64> = (0..ns).map(|q| preds[q].1[i * nu + j]).collect(); let mu = xs.iter().sum::<f64>() / ns as f64; v += xs.iter().map(|x| (x - mu).powi(2)).sum::<f64>() / ns as f64; m2 += mu * mu; }
                v / m2.max(1e-12)
            }).collect();
            for (rname, im) in &preds {
                let rows: Vec<String> = (0..nt).into_par_iter().map(|i| {
                    let (y, t, d, dd) = (&im[i * nu..(i + 1) * nu], &tr[i * nu..(i + 1) * nu], &dr[i * nu..(i + 1) * nu], &den[i * nu..(i + 1) * nu]);
                    let rt = resid(t, &basis);
                    let nde = t.iter().filter(|x| x.abs() > DE_Z).count();
                    let nc = d.iter().filter(|x| x.abs() > DE_Z).count();
                    let mut o: Vec<usize> = (0..nu).collect();
                    o.sort_by(|&a, &b| y[b].abs().partial_cmp(&y[a].abs()).unwrap());
                    let (mut tp, mut disc) = (0, 0);
                    for &j in &o[..nc] { let same = y[j].signum() == t[j].signum(); if !same { disc += 1; } if same && t[j].abs() > DE_Z { tp += 1; } }
                    let (sxy, sxx): (f64, f64) = (y.iter().zip(t).map(|(p, q)| p * q).sum(), y.iter().map(|p| p * p).sum());
                    format!("{tname}\t{pname}\t{rname}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{nde}\t{nc}\t{tp}\t{disc}\t{:.6}\t{sxy:.6e}\t{sxx:.6e}\t{:.6}\t{:.6}",
                        zt.perts[tests[i]], pearson(y, t), pearson(&resid(y, &basis), &rt), pearson(&resid(dd, &basis), &rt), pearson(d, t), pearson(&resid(d, &basis), &rt), spa_de(y, t), oor[i], dis[i])
                }).collect();
                for row in rows { writeln!(w, "{row}")?; }
            }
            writeln!(wl, "{tname} {pname} done at {:.1}s", t0.elapsed().as_secs_f64())?;
            wl.flush()?;
        }
    }
    Ok(())
}

/// summary of {pre}.multi.tsv per target, panel and reference (responding perturbations r_ceil >= 0.2)
pub fn summary(a: &[String]) -> Result<()> {
    let pre = req(a, "--in")?;
    let txt = std::fs::read_to_string(format!("{pre}.multi.tsv"))?;
    let rows: Vec<Vec<String>> = txt.lines().skip(1).map(|l| l.split('\t').map(String::from).collect()).collect();
    let f = |r: &Vec<String>, i: usize| -> f64 { r[i].parse().unwrap_or(f64::NAN) };
    // responders: r_ceil >= 0.2 (default) or, with --spec, specific ceiling rperp_ceil >= 0.2
    let rc = if a.contains(&"--spec".to_string()) { 8 } else { 7 };
    let mut keys: Vec<(String, String, String)> = vec![];
    for r in &rows { let k = (r[0].clone(), r[1].clone(), r[2].clone()); if !keys.contains(&k) { keys.push(k); } }
    println!("target  panel             ref           n    recovery  spec_r  recall  opp    rank   slope  | AUC fail: OOR  disagree");
    for (t, p, rf) in &keys {
        let s: Vec<&Vec<String>> = rows.iter().filter(|r| &r[0] == t && &r[1] == p && &r[2] == rf).collect();
        let resp: Vec<&&Vec<String>> = s.iter().filter(|r| f(r, rc) >= 0.2).collect();
        let fails: Vec<bool> = resp.iter().map(|r| f(r, 5) / f(r, 8) < 0.5).collect();
        let sum = |i: usize| s.iter().map(|r| f(r, i)).sum::<f64>();
        println!("{t:7} {p:17} {rf:12} {:5} {:.3}     {:.3}   {:.3}   {:.3}  {:.3}  {:.3}  | {:.3} {:.3} ({})", resp.len(),
            median(resp.iter().map(|r| f(r, 5) / f(r, 6)).collect()), median(resp.iter().map(|r| f(r, 5)).collect()),
            sum(11) / sum(9), sum(12) / sum(10).max(1.0), median(s.iter().map(|r| f(r, 13)).collect()), sum(14) / sum(15),
            auc(&resp.iter().map(|r| f(r, 16)).collect::<Vec<_>>(), &fails).unwrap_or(f64::NAN), auc(&resp.iter().map(|r| f(r, 17)).collect::<Vec<_>>(), &fails).unwrap_or(f64::NAN), fails.iter().filter(|x| **x).count());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn nnls_recovers_nonnegative_weights() {
        let a = vec![vec![1.0, 0.0, 1.0, 2.0], vec![0.0, 1.0, 1.0, -1.0], vec![1.0, 1.0, 0.0, 0.0]];
        let b: Vec<f64> = (0..4).map(|t| 0.7 * a[0][t] + 0.2 * a[2][t]).collect();
        let x = nnls_small(&a, &b);
        assert!((x[0] - 0.7).abs() < 1e-9 && x[1].abs() < 1e-9 && (x[2] - 0.2).abs() < 1e-9, "{x:?}");
        let b2: Vec<f64> = (0..4).map(|t| -a[1][t]).collect();
        assert!(nnls_small(&a, &b2).iter().all(|&v| v >= 0.0));
    }
}

/// per-gene signal variance: mean z^2 over perturbations (all cells) minus the noise variance from half differences
pub fn genevar(a: &[String], z: &Z) -> Result<()> {
    let (p, g) = (z.perts.len(), z.genes.len());
    let rows: Vec<String> = (0..g).into_par_iter().map(|j| {
        let m2 = (0..p).map(|i| z.at(0, i, j).powi(2)).sum::<f64>() / p as f64;
        let nz = (0..p).map(|i| (z.at(1, i, j) - z.at(2, i, j)).powi(2)).sum::<f64>() / (2.0 * p as f64);
        format!("{}\t{m2:.6}\t{nz:.6}\t{:.6}", z.genes[j], (m2 - nz).max(0.0))
    }).collect();
    let mut w = BufWriter::new(std::fs::File::create(req(a, "--out")?)?);
    writeln!(w, "gene\tmean_z2\tnoise_var\tsignal_var")?;
    for r in rows { writeln!(w, "{r}")?; }
    Ok(())
}

/// Number of reference programs above the Marchenko-Pastur edge, with the noise variance of the
/// full-cell z taken from half differences (as in `upg::run`), over all perturbations.
pub fn mpcount(_a: &[String], z: &Z) -> Result<()> {
    let (p, g) = (z.perts.len(), z.genes.len());
    let all: Vec<usize> = (0..p).collect();
    let x = z.mat(0, &all, &(0..g).collect::<Vec<_>>());
    let s2 = all.par_iter().map(|&i| (0..g).map(|j| (z.at(1, i, j) - z.at(2, i, j)).powi(2)).sum::<f64>()).sum::<f64>() / (2.0 * (p * g) as f64);
    let c = la::gram(&x, p, g);
    drop(x);
    let sp = crate::upg::spec_from(&c, g, p, s2, 500)?;
    println!("perts {p} genes {g} s2 {s2:.4} kMP {}", sp.kmp);
    Ok(())
}

/// Pooled reference for a pan-atlas model: perturbations of several z prefixes (same gene order)
/// stacked, rows of atlas r scaled by sqrt(N / (R n_r)) so that the Gram matrix is the equal-weight
/// average of the per-atlas Gram matrices (the `pan_equal` reference of `run`). Names get an atlas prefix.
pub fn concat(a: &[String]) -> Result<()> {
    let ins: Vec<(String, Z)> = req(a, "--in")?.split(',').map(|kv| { let (k, v) = kv.split_once('=').unwrap(); Ok((k.to_string(), Z::load(v)?)) }).collect::<Result<_>>()?;
    let g = ins[0].1.genes.len();
    for (n, z) in &ins { if z.genes != ins[0].1.genes { bail!("gene order differs for {n}") } }
    let ntot: usize = ins.iter().map(|(_, z)| z.perts.len()).sum();
    let nr = ins.len() as f64;
    let o = req(a, "--out")?;
    let mut out = vec![0f32; 3 * ntot * g];
    let mut off = 0;
    let mut w = BufWriter::new(std::fs::File::create(format!("{o}.meta.tsv"))?);
    writeln!(w, "pert\tatlas\tweight")?;
    for (name, z) in &ins {
        let p = z.perts.len();
        let wt = (ntot as f64 / (nr * p as f64)).sqrt() as f32;
        for l in 0..3 { for i in 0..p { for j in 0..g { out[(l * ntot + off + i) * g + j] = z.z[(l * p + i) * g + j] * wt; } } }
        for s in &z.perts { writeln!(w, "{name}:{s}\t{name}\t{wt:.4}")?; }
        eprintln!("concat: {name} {p} perturbations, row weight {wt:.4}");
        off += p;
    }
    std::fs::write(format!("{o}.z.f32"), bytemuck::cast_slice(&out))?;
    let mut w = BufWriter::new(std::fs::File::create(format!("{o}.genes.tsv"))?);
    writeln!(w, "gene")?;
    for s in &ins[0].1.genes { writeln!(w, "{s}")?; }
    Ok(())
}

/// paired comparison of the pooled reference with each single atlas: recovery per specific
/// responder (rperp_ceil >= 0.2), Wilcoxon signed-rank P and share of perturbations improved
pub fn multitest(a: &[String]) -> Result<()> {
    let pre = req(a, "--in")?;
    let panel = arg(a, "--panel").unwrap_or("greedy_pan".into());
    let txt = std::fs::read_to_string(format!("{pre}.multi.tsv"))?;
    let mut rec: HashMap<(String, String), HashMap<String, f64>> = HashMap::new();
    let mut targets: Vec<String> = vec![];
    for l in txt.lines().skip(1) {
        let r: Vec<&str> = l.split('\t').collect();
        if r[1] != panel { continue; }
        let (rp, den, rpc): (f64, f64, f64) = (r[5].parse()?, r[6].parse()?, r[8].parse()?);
        if !(rpc >= 0.2) { continue; }
        if !targets.contains(&r[0].to_string()) { targets.push(r[0].into()); }
        rec.entry((r[0].into(), r[2].into())).or_default().insert(r[3].into(), rp / den);
    }
    println!("target  ref          vs pan_equal: median_single median_pan  share_pan_higher  P");
    for t in &targets {
        let Some(pan) = rec.get(&(t.clone(), "pan_equal".into())) else { continue };
        let mut refs: Vec<&String> = rec.keys().filter(|(tt, r)| tt == t && r != "pan_equal").map(|(_, r)| r).collect();
        refs.sort();
        for rf in refs {
            let s = &rec[&(t.clone(), rf.clone())];
            let d: Vec<(f64, f64)> = s.iter().filter_map(|(p, v)| pan.get(p).map(|w| (*v, *w))).filter(|(v, w)| v.is_finite() && w.is_finite()).collect();
            let diff: Vec<f64> = d.iter().map(|(v, w)| w - v).collect();
            let hi = diff.iter().filter(|x| **x > 0.0).count() as f64 / diff.len() as f64;
            println!("{t:7} {rf:12} {:.3} {:.3}  {hi:.2}  {:.2e}  (n {})", median(d.iter().map(|x| x.0).collect()), median(d.iter().map(|x| x.1).collect()), crate::upg::wilcoxon(&diff), diff.len());
        }
    }
    Ok(())
}
