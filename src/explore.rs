//! Additional analyses (Additional file 1: Table S1, Note S1): what ridge imputation leaves unexplained.
//!   panelimp explore --z PREFIX --folds F --hubs H --out PREFIX [--fold 0] [--tss FILE] [--umi FILE] [--stage diag|test]
//! Stage `diag` uses the training fold only (model on the inner training set, inner holdout as test):
//! reproducible share of the variance between halves, target gene and neighbours, calibration by strength,
//! and the noise shared between genes inside one half. Stage `test` is in explore_test.rs.
use crate::causal::*;
use crate::*;

pub(crate) struct Ctx { pub d: Data, pub refs: Vec<usize>, pub itr: Vec<usize>, pub iho: Vec<usize>, pub tests: Vec<usize>, pub s: Vec<usize>, pub u: Vec<usize>, pub alpha: f64 }
pub(crate) fn setup(a: &[String]) -> Result<Ctx> {
    let d = Data::open(&req(a, "--z")?)?;
    let f: i64 = arg(a, "--fold").map(|v| v.parse().unwrap()).unwrap_or(0);
    let (g, np) = (d.genes.len(), d.perts.len());
    let pi: HashMap<&str, usize> = d.perts.iter().enumerate().map(|(i, p)| (p.as_str(), i)).collect();
    let (mut fold, mut group) = (vec![-1i64; np], vec![-1i64; np]);
    for l in std::fs::read_to_string(req(a, "--folds")?)?.lines().skip(1) { let x: Vec<&str> = l.split('\t').collect(); if let Some(&i) = pi.get(x[0]) { group[i] = x[1].parse()?; fold[i] = x[2].parse()?; } }
    let mut excl: HashSet<(usize, i64)> = HashSet::new();
    for l in std::fs::read_to_string(req(a, "--hubs")?)?.lines().skip(1) { let x: Vec<&str> = l.split('\t').collect(); if let Some(&i) = pi.get(x[0]) { excl.insert((i, x[1].parse()?)); } }
    let refs: Vec<usize> = (0..np).filter(|&i| fold[i] != f && !excl.contains(&(i, f))).collect();
    let tests: Vec<usize> = (0..np).filter(|&i| fold[i] == f).collect();
    let key = |i: usize| -> u64 { let k = if group[i] >= 0 { group[i] as u64 } else { 1_000_000 + i as u64 }; (k.wrapping_mul(0x9E3779B97F4A7C15) >> 32) % 10 };
    let (itr, iho): (Vec<usize>, Vec<usize>) = (refs.iter().copied().filter(|&i| key(i) != 0).collect(), refs.iter().copied().filter(|&i| key(i) == 0).collect());
    let c = sgram(&d.rows(0, &refs), refs.len(), g);
    let s = greedy32(&c, refs.len(), g, 500);
    drop(c);
    let ins: HashSet<usize> = s.iter().copied().collect();
    let u: Vec<usize> = (0..g).filter(|j| !ins.contains(j)).collect();
    let c_in = sgram(&d.rows(0, &itr), itr.len(), g);
    let (xh, th) = (cols(&d.rows(1, &iho), g, &s), cols(&d.rows(2, &iho), g, &u));
    let mut best = (f64::NEG_INFINITY, ALPHAS[0]);
    for &al in &ALPHAS {
        let p = sgemm(false, true, iho.len(), u.len(), 500, &xh, &ridge32(&c_in, g, &s, &u, al)?);
        let r: f64 = p.par_chunks(u.len()).zip(th.par_chunks(u.len())).map(|(x, t)| pearson(&x.iter().map(|&v| v as f64).collect::<Vec<_>>(), &t.iter().map(|&v| v as f64).collect::<Vec<_>>())).filter(|v| v.is_finite()).sum::<f64>() / iho.len() as f64;
        if r > best.0 { best = (r, al); }
    }
    Ok(Ctx { d, refs, itr, iho, tests, s, u, alpha: best.1 })
}
/// TSS table with columns gene, chr, tss (header names): gene -> (chromosome, position)
pub(crate) fn read_tss(f: &str) -> Result<Vec<(String, (String, f64))>> {
    let txt = std::fs::read_to_string(f)?;
    let mut lines = txt.lines();
    let hdr: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let col = |n: &str| hdr.iter().position(|h| *h == n).ok_or_else(|| anyhow!("{f}: no column {n}"));
    let (cg, cc, ct) = (col("gene")?, col("chr")?, col("tss")?);
    Ok(lines.filter_map(|l| { let x: Vec<&str> = l.split('\t').collect(); x.get(ct).and_then(|v| v.parse::<f64>().ok()).map(|t| (x[cg].to_string(), (x[cc].to_string(), t))) }).collect())
}
/// pooled mean of products
fn mp(a: &[f32], b: &[f32]) -> f64 { a.par_iter().zip(b.par_iter()).map(|(x, y)| (*x as f64) * (*y as f64)).sum::<f64>() / a.len() as f64 }
fn sub(a: &[f32], b: &[f32]) -> Vec<f32> { a.par_iter().zip(b.par_iter()).map(|(x, y)| x - y).collect() }
/// rows with the span of the orthonormal columns of b (nu x k) removed
fn deflate(m: &[f32], nu: usize, b: &[f32]) -> Vec<f32> { let (n, k) = (m.len() / nu, b.len() / nu); sub(m, &sgemm(false, true, n, nu, k, &sgemm(false, false, n, k, nu, m, b), b)) }

pub fn run(a: &[String]) -> Result<()> {
    match arg(a, "--stage").as_deref() { Some("test") => crate::explore_test::run(a), _ => diag(a) }
}

fn diag(a: &[String]) -> Result<()> {
    let t0 = std::time::Instant::now();
    let c = setup(a)?;
    let (d, s, u, itr, iho) = (&c.d, &c.s, &c.u, &c.itr, &c.iho);
    let (g, nu, k, nh) = (d.genes.len(), u.len(), s.len(), iho.len());
    let mut log = BufWriter::new(std::fs::File::create(format!("{}.explore.diag.log", req(a, "--out")?))?);
    macro_rules! say { ($($t:tt)*) => {{ let s = format!($($t)*); println!("{s}"); writeln!(log, "{s}")?; }} }
    say!("inner training {} perturbations, inner holdout {nh}, panel {k}, genes outside the panel {nu}, ridge alpha {}", itr.len(), c.alpha);
    let c_in = sgram(&d.rows(0, itr), itr.len(), g);
    let w = ridge32(&c_in, g, s, u, c.alpha)?;
    let vbig = top_eig32(&c_in, g, NPC, 30)?;
    drop(c_in);
    let v10: Vec<f64> = (0..g).flat_map(|j| vbig[j * NPC..j * NPC + NPC].iter().map(|&v| v as f64).collect::<Vec<_>>()).collect();
    let b10: Vec<f32> = { let b = basis_on(&v10, NPC, u); (0..nu).flat_map(|q| b.iter().map(|col| col[q] as f32).collect::<Vec<_>>()).collect() };
    let (z1, z2) = (d.rows(1, iho), d.rows(2, iho));
    let (x1, t1, t2) = (cols(&z1, g, s), cols(&z1, g, u), cols(&z2, g, u));
    let pb = sgemm(false, true, nh, nu, k, &x1, &w);
    // ---- 1. how much of the variance is reproducible between the two halves of the cells
    say!("\n[1] variance of the held-out half, genes outside the panel (z units of a half)");
    for (frame, f) in [("all programs", None), ("specific (10 programs removed)", Some(&b10))] {
        let fr = |m: &[f32]| -> Vec<f32> { match f { Some(b) => deflate(m, nu, b), None => m.to_vec() } };
        let (t1f, t2f) = (fr(&t1), fr(&t2));
        let (tot, sig) = (mp(&t2f, &t2f), mp(&t1f, &t2f));
        say!("  {frame}: total variance of the held-out half {tot:.3}; reproducible between halves {sig:.3} = {:.1}% of the total", 100.0 * sig / tot);
    }
    // ---- 2. where the reproducible error sits
    say!("\n[2] target gene, neighbours, calibration and shared noise");
    let sig_tot = mp(&t1, &t2) * (nh * nu) as f64;
    // (a) the knocked-down gene itself
    let gi: HashMap<&str, usize> = d.genes.iter().enumerate().map(|(i, x)| (x.as_str(), i)).collect();
    let upos: HashMap<usize, usize> = u.iter().enumerate().map(|(q, &j)| (j, q)).collect();
    let (mut own, mut n_own, mut z_own, mut p_own) = (0.0f64, 0usize, vec![], vec![]);
    for (b, &p) in iho.iter().enumerate() {
        if let Some(&q) = gi.get(d.perts[p].as_str()).and_then(|j| upos.get(j)) { own += (t1[b * nu + q] as f64) * (t2[b * nu + q] as f64); n_own += 1; z_own.push(t2[b * nu + q] as f64); p_own.push(pb[b * nu + q] as f64); }
    }
    say!("  (a) target gene of the knockdown: outside the panel for {n_own} of {nh} perturbations; median measured z {:.1}, median ridge prediction {:.2};", median(z_own.clone()), median(p_own));
    say!("      this single gene holds {:.1}% of the reproducible variance ({} of {} gene and perturbation pairs)", 100.0 * own / sig_tot, n_own, nh * nu);
    // (b) genomic neighbours of the target
    if let Some(tf) = arg(a, "--tss") {
        let mut tss: HashMap<String, (String, f64)> = HashMap::new();
        for (name, v) in read_tss(&tf)? { tss.entry(name).or_insert(v); }
        for (lo, hi) in [(0.0, 2e3), (2e3, 2e4), (2e4, 2e5), (2e5, 1e6)] {
            let (mut n, mut zs) = (0usize, vec![]);
            for (b, &p) in iho.iter().enumerate() {
                let Some((chr, pos)) = tss.get(&d.perts[p]) else { continue };
                for (q, &j) in u.iter().enumerate() {
                    if d.genes[j] == d.perts[p] { continue; }
                    if let Some((c2, p2)) = tss.get(&d.genes[j]) { let dist = (p2 - pos).abs(); if c2 == chr && dist >= lo && dist < hi { n += 1; zs.push(t2[b * nu + q] as f64); } }
                }
            }
            say!("  (b) genes with a TSS {:.0} to {:.0} kb from the target: {n} pairs, median z {:.2}", lo / 1e3, hi / 1e3, median(zs));
        }
    }
    // (f) calibration by strength of the prediction
    for (lo, hi) in [(0.0f32, 1.0f32), (1.0, 2.0), (2.0, 4.0), (4.0, 8.0), (8.0, f32::INFINITY)] {
        let (mut xy, mut xx, mut n) = (0.0f64, 0.0f64, 0usize);
        for (p, t) in pb.iter().zip(&t2) { if p.abs() >= lo && p.abs() < hi { xy += (*p as f64) * (*t as f64); xx += (*p as f64) * (*p as f64); n += 1; } }
        say!("  (f) predictions with |z| in [{lo}, {hi}): {n} pairs, slope of the held-out half on the prediction {:.2}", xy / xx);
    }
    // (g) noise shared between genes inside one half: panel block of (C11 + C22) / 2 - C12
    {
        let (a1, a2) = (cols(&d.rows(1, itr), g, s), cols(&d.rows(2, itr), g, s));
        let n = itr.len() as f32;
        let (c11, c22, c12) = (sgram(&a1, itr.len(), k), sgram(&a2, itr.len(), k), sgemm(true, false, k, k, itr.len(), &a1, &a2));
        let nz: Vec<f64> = (0..k * k).map(|x| { let (i, j) = (x / k, x % k); ((0.5 * (c11[x] + c22[x]) - 0.5 * (c12[x] + c12[j * k + i])) / n) as f64 }).collect();
        let dg: f64 = (0..k).map(|i| nz[i * k + i]).sum::<f64>() / k as f64;
        let (fro, dsq): (f64, f64) = (nz.iter().map(|v| v * v).sum(), (0..k).map(|i| nz[i * k + i].powi(2)).sum());
        let (_, ev) = la::top_eig(&nz, k, 5)?;
        let sg: f64 = (0..k).map(|i| (0.5 * (c12[i * k + i] + c12[i * k + i]) / n) as f64).sum::<f64>() / k as f64;
        say!("  (g) noise of the 500 panel genes in a half: mean variance {dg:.2} (signal variance {sg:.2}); off-diagonal share of the noise covariance {:.1}%; leading noise eigenvalues {:.1} {:.1} {:.1} (1.0 for independent noise)", 100.0 * (fro - dsq) / fro, ev[0], ev[1], ev[2]);
    }
    say!("\nresident {:.0} MB, {:.1}s", rss_mb(), t0.elapsed().as_secs_f64());
    Ok(())
}
