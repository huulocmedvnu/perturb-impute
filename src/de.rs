//! DE calls from imputed profiles and the hybrid readout. A1 posterior-scaled calls with INFO gate, A2 Huber IRLS
//! reweighting, B communality-corrected size ranking, C unsupervised panel alignment of the
//! target screen (CORAL, reliability ratio), D hybrid measured-panel plus imputed readout.
use crate::*;

/// symmetric matrix power V diag(max(l, eps)^p) V' (n x n)
pub fn sym_pow(a: &[f64], n: usize, p: f64) -> Result<Vec<f64>> {
    let (v, l) = la::top_eig(a, n, n)?;
    let vs: Vec<f64> = (0..n).flat_map(|i| (0..n).map(move |k| (i, k))).map(|(i, k)| v[i * n + k] * l[k].max(1e-9).powf(p)).collect();
    Ok(la::gemm(false, true, n, n, n, &vs, &v))
}

/// matched calls: top `ncall` genes by score; (correct-sign DE, opposite sign)
fn calls(score: &[f64], val: &[f64], t: &[f64], ncall: usize) -> (usize, usize) {
    let mut o: Vec<usize> = (0..score.len()).collect();
    o.sort_by(|&a, &b| score[b].partial_cmp(&score[a]).unwrap().then(a.cmp(&b)));
    let (mut tp, mut disc) = (0, 0);
    for &j in &o[..ncall.min(o.len())] {
        let same = val[j].signum() == t[j].signum();
        if !same { disc += 1; }
        if same && t[j].abs() > DE_Z { tp += 1; }
    }
    (tp, disc)
}

/// recall with the correct sign among DE genes binned by |truth z|: (n_de, tp) per bin [4,6), [6,10), [10,inf)
fn calls_bins(score: &[f64], val: &[f64], t: &[f64], ncall: usize) -> [(usize, usize); 3] {
    let mut o: Vec<usize> = (0..score.len()).collect();
    o.sort_by(|&a, &b| score[b].partial_cmp(&score[a]).unwrap().then(a.cmp(&b)));
    let bin = |z: f64| if z < 6.0 { 0 } else if z < 10.0 { 1 } else { 2 };
    let mut r = [(0usize, 0usize); 3];
    for j in 0..t.len() { if t[j].abs() > DE_Z { r[bin(t[j].abs())].0 += 1; } }
    for &j in &o[..ncall.min(o.len())] { if t[j].abs() > DE_Z && val[j].signum() == t[j].signum() { r[bin(t[j].abs())].1 += 1; } }
    r
}

/// per-gene stats on a set of holdout perturbations: (INFO r, slope b, error variance e)
fn gene_stats(im: &[f64], tr: &[f64], n: usize, nu: usize) -> Vec<(f64, f64, f64)> {
    (0..nu).into_par_iter().map(|j| {
        let x: Vec<f64> = (0..n).map(|i| im[i * nu + j]).collect();
        let y: Vec<f64> = (0..n).map(|i| tr[i * nu + j]).collect();
        let sxx: f64 = x.iter().map(|v| v * v).sum();
        let b = if sxx > 0.0 { x.iter().zip(&y).map(|(p, q)| p * q).sum::<f64>() / sxx } else { 0.0 };
        let mse = x.iter().zip(&y).map(|(p, q)| (q - b * p).powi(2)).sum::<f64>() / n as f64;
        let r = pearson(&x, &y);
        (if r.is_finite() { r } else { 0.0 }, b, (mse - 1.0).max(0.05))
    }).collect()
}

/// ridge weights W' (|U| x k) from a Gram matrix (same scaling as `ridge`)
fn ridge_w(c: &[f64], g: usize, s: &[usize], u: &[usize], alpha: f64) -> Result<Vec<f64>> { ridge(c, g, s, u, alpha) }

/// Gram of rows scaled by sqrt(w)
fn wgram(z: &Z, rows: &[usize], w: &[f64], g: usize) -> Vec<f64> {
    let mut y = z.mat(0, rows, &(0..g).collect::<Vec<_>>());
    y.par_chunks_mut(g).zip(w).for_each(|(r, &wi)| { let s = wi.sqrt(); r.iter_mut().for_each(|v| *v *= s); });
    la::gram(&y, rows.len(), g)
}

#[allow(clippy::too_many_arguments)]
pub fn run(a: &[String], z: &Z, zt: &Z, fold: &[i64], group: &[i64], excl: &HashSet<(usize, i64)>, out: &str) -> Result<()> {
    let t0 = std::time::Instant::now();
    let g = z.genes.len();
    let transfer = !std::ptr::eq(z, zt);
    let f: i64 = if transfer { -9 } else { arg(a, "--fold").map(|v| v.parse().unwrap()).unwrap_or(0) };
    let k: usize = arg(a, "--panel-size").map(|v| v.parse().unwrap()).unwrap_or(500);
    let mut wl = BufWriter::new(std::fs::File::create(format!("{out}.de.log"))?);
    let hk = |i: usize| -> u64 { let kk = if group[i] >= 0 { group[i] as u64 } else { 1_000_000 + i as u64 }; kk.wrapping_mul(0x9E3779B97F4A7C15) };
    let refs: Vec<usize> = (0..z.perts.len()).filter(|&i| transfer || (fold[i] != f && !excl.contains(&(i, f)))).collect();
    let tests: Vec<usize> = if transfer { (0..zt.perts.len()).collect() } else { (0..z.perts.len()).filter(|&i| fold[i] == f).collect() };
    let itr: Vec<usize> = refs.iter().copied().filter(|&i| (hk(i) >> 32) % 10 != 0).collect();
    let iho: Vec<usize> = refs.iter().copied().filter(|&i| (hk(i) >> 32) % 10 == 0).collect();
    let ia: Vec<usize> = iho.iter().copied().filter(|&i| (hk(i) >> 20) % 2 == 0).collect();
    let ib: Vec<usize> = iho.iter().copied().filter(|&i| (hk(i) >> 20) % 2 == 1).collect();
    let gr = grams_c(z, &refs, &itr);
    let s = sim::greedy(&gr, refs.len(), g, k, &vec![1.0; g]);
    let inp: HashSet<usize> = s.iter().copied().collect();
    let u: Vec<usize> = (0..g).filter(|j| !inp.contains(j)).collect();
    let nu = u.len();
    writeln!(wl, "fold {f}: refs {} itr {} iho {} (A {}, B {}) tests {}; setup {:.1}s", refs.len(), itr.len(), iho.len(), ia.len(), ib.len(), tests.len(), t0.elapsed().as_secs_f64())?;
    let mean_r = |im: &[f64], t: &[f64], n: usize| -> f64 { (0..n).map(|i| pearson(&im[i * nu..(i + 1) * nu], &t[i * nu..(i + 1) * nu])).filter(|x| x.is_finite()).sum::<f64>() / n as f64 };
    // ridge alpha on the whole inner holdout (as before)
    let xh = z.mat(1, &iho, &s);
    let th = z.mat(2, &iho, &u);
    let mut best = (f64::NEG_INFINITY, 0.0);
    for &al in &ALPHAS {
        let wt = ridge_w(&gr.c_in, g, &s, &u, al)?;
        let m = mean_r(&la::gemm(false, true, iho.len(), nu, k, &xh, &wt), &th, iho.len());
        if m > best.0 { best = (m, al); }
    }
    let alpha = best.1;
    // inner A / B matrices
    let (xa, ta, da) = (z.mat(1, &ia, &s), z.mat(2, &ia, &u), z.mat(1, &ia, &u));
    let (xb, tb, db) = (z.mat(1, &ib, &s), z.mat(2, &ib, &u), z.mat(1, &ib, &u));
    let (na, nb) = (ia.len(), ib.len());
    let ncb: Vec<usize> = (0..nb).map(|i| db[i * nu..(i + 1) * nu].iter().filter(|v| v.abs() > DE_Z).count()).collect();
    let recall_b = |pred: &[f64], score: &dyn Fn(usize, usize, f64) -> f64| -> (f64, f64) {
        let (mut tp, mut disc, mut nde, mut nc) = (0usize, 0usize, 0usize, 0usize);
        for i in 0..nb {
            let p = &pred[i * nu..(i + 1) * nu];
            let sc: Vec<f64> = (0..nu).map(|j| score(i, j, p[j])).collect();
            let (a1, b1) = calls(&sc, p, &tb[i * nu..(i + 1) * nu], ncb[i]);
            tp += a1; disc += b1; nde += tb[i * nu..(i + 1) * nu].iter().filter(|v| v.abs() > DE_Z).count(); nc += ncb[i];
        }
        (tp as f64 / nde as f64, disc as f64 / nc.max(1) as f64)
    };
    // A2: Huber IRLS on the inner training set, delta chosen on B
    let xi = z.mat(0, &itr, &s);
    let yi = z.mat(0, &itr, &u);
    let w_in0 = ridge_w(&gr.c_in, g, &s, &u, alpha)?;
    let fit_i = la::gemm(false, true, itr.len(), nu, k, &xi, &w_in0);
    let resn: Vec<f64> = (0..itr.len()).map(|i| (0..nu).map(|j| (yi[i * nu + j] - fit_i[i * nu + j]).powi(2)).sum::<f64>().sqrt()).collect();
    let mut sorted = resn.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut best_h = (f64::NEG_INFINITY, f64::INFINITY);
    for &qq in &[f64::INFINITY, 0.9, 0.75, 0.5] {
        let wt = if qq.is_infinite() { w_in0.clone() } else {
            let d = sorted[((sorted.len() as f64 - 1.0) * qq) as usize];
            let w: Vec<f64> = resn.iter().map(|r| if *r > d { d / r } else { 1.0 }).collect();
            ridge_w(&wgram(z, &itr, &w, g), g, &s, &u, alpha)?
        };
        let pb = la::gemm(false, true, nb, nu, k, &xb, &wt);
        let (rec, _) = recall_b(&pb, &|_, _, v| v.abs());
        writeln!(wl, "  A2 huber quantile {qq}: inner-B recall {rec:.4}")?;
        if rec > best_h.0 { best_h = (rec, qq); }
    }
    let huber_q = best_h.1;
    let w_in = if huber_q.is_infinite() { w_in0.clone() } else {
        let d = sorted[((sorted.len() as f64 - 1.0) * huber_q) as usize];
        let w: Vec<f64> = resn.iter().map(|r| if *r > d { d / r } else { 1.0 }).collect();
        ridge_w(&wgram(z, &itr, &w, g), g, &s, &u, alpha)?
    };
    // gene stats from A, tuning on B
    let pa = la::gemm(false, true, na, nu, k, &xa, &w_in);
    let pb = la::gemm(false, true, nb, nu, k, &xb, &w_in);
    let st_a = gene_stats(&pa, &ta, na, nu);
    let _ = &da;
    let (base_rec0, _) = recall_b(&pb, &|_, _, v| v.abs());
    let mut best_a = (base_rec0, -1.0, -1.0);
    for &kap in &[0.0, 0.25, 0.5, 1.0] { for &gate in &[-1.0, 0.05, 0.1] {
        let (rec, _) = recall_b(&pb, &|_, j, v| if st_a[j].0 < gate { f64::NEG_INFINITY } else { (st_a[j].1 * v).abs() / st_a[j].2.powf(kap) });
        if rec > best_a.0 { best_a = (rec, kap, gate); }
    } }
    let (base_rec, _) = recall_b(&pb, &|_, _, v| v.abs());
    let spa_b = |score: &dyn Fn(usize, f64) -> f64| -> f64 { median((0..nb).map(|i| { let p = &pb[i * nu..(i + 1) * nu]; spa_de(&(0..nu).map(|j| score(j, p[j])).collect::<Vec<_>>(), &tb[i * nu..(i + 1) * nu]) }).collect()) };
    let mut best_b = (f64::NEG_INFINITY, 0.0, 0.0);
    for &tau in &[0.01, 0.05, 0.1, 0.25] { for &gam in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        let m = spa_b(&|j, v| v.abs() / (st_a[j].0.max(0.0).powi(2) + tau).powf(gam));
        if m > best_b.0 { best_b = (m, tau, gam); }
    } }
    writeln!(wl, "  alpha {alpha}; A2 huber {huber_q}; A1 kappa {} gate {} (inner-B recall {:.4} vs {:.4}); B tau {} gamma {} (inner-B size rank {:.4})", best_a.1, best_a.2, best_a.0, base_rec, best_b.1, best_b.2, best_b.0)?;
    // final model on all references (Huber weights from the reference fit if chosen)
    let w0 = ridge_w(&gr.c, g, &s, &u, alpha)?;
    let wfin = if huber_q.is_infinite() { w0.clone() } else {
        let xr = z.mat(0, &refs, &s); let yr = z.mat(0, &refs, &u);
        let fr = la::gemm(false, true, refs.len(), nu, k, &xr, &w0);
        let rn: Vec<f64> = (0..refs.len()).map(|i| (0..nu).map(|j| (yr[i * nu + j] - fr[i * nu + j]).powi(2)).sum::<f64>().sqrt()).collect();
        let mut so = rn.clone(); so.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let d = so[((so.len() as f64 - 1.0) * huber_q) as usize];
        let w: Vec<f64> = rn.iter().map(|r| if *r > d { d / r } else { 1.0 }).collect();
        ridge_w(&wgram(z, &refs, &w, g), g, &s, &u, alpha)?
    };
    let pih = la::gemm(false, true, iho.len(), nu, k, &xh, &w_in);
    let st = gene_stats(&pih, &th, iho.len(), nu); // final INFO, b, e from the whole inner holdout
    // test inputs, raw and aligned (C)
    let nt = tests.len();
    // --same-cells: input panel genes and truth from all cells of each test perturbation (layer 0)
    let same = a.contains(&"--same-cells".to_string());
    let (lin, ltr) = if same { (0, 0) } else { (1, 2) };
    let xt = zt.mat(lin, &tests, &s);
    let sr: Vec<f64> = { let mut m = vec![0.0; k * k]; for (i, &a1) in s.iter().enumerate() { for (j, &b1) in s.iter().enumerate() { m[i * k + j] = gr.c[a1 * g + b1] / refs.len() as f64; } } m };
    let stt = { let mut m = la::gemm(true, false, k, k, nt, &xt, &xt); m.iter_mut().for_each(|v| *v /= nt as f64); m };
    let mut sr_i = sr.clone(); for i in 0..k { sr_i[i * k + i] += 1.0; }
    let mut st_i = stt.clone(); for i in 0..k { st_i[i * k + i] += 1.0; }
    let coral = la::gemm(false, false, k, k, k, &sym_pow(&st_i, k, -0.5)?, &sym_pow(&sr_i, k, 0.5)?);
    let x_coral = la::gemm(false, false, nt, k, k, &xt, &coral);
    let sig_t: Vec<f64> = (0..k).map(|i| (stt[i * k + i] - 1.0).max(0.0)).collect();
    let sig_r: Vec<f64> = (0..k).map(|i| (sr[i * k + i] - 1.0).max(0.0)).collect();
    let kap_s = sig_t.iter().sum::<f64>() / sig_r.iter().sum::<f64>();
    let rel: Vec<f64> = (0..k).map(|i| { let a1 = sig_t[i] / (sig_t[i] + 1.0); let b1 = kap_s * sig_r[i] / (kap_s * sig_r[i] + 1.0); if b1 > 0.0 { (a1 / b1).min(2.0) } else { 0.0 } }).collect();
    let x_rel: Vec<f64> = xt.chunks(k).flat_map(|r| r.iter().zip(&rel).map(|(v, w)| v * w).collect::<Vec<_>>()).collect();
    writeln!(wl, "  C: signal scale target/reference {kap_s:.3}, median reliability ratio {:.3}", median(rel.clone()))?;
    let p_raw = la::gemm(false, true, nt, nu, k, &xt, &w0);
    let p_a2 = la::gemm(false, true, nt, nu, k, &xt, &wfin);
    let p_cor = la::gemm(false, true, nt, nu, k, &x_coral, &wfin);
    let p_rel = la::gemm(false, true, nt, nu, k, &x_rel, &wfin);
    // D: hybrid, panel genes measured and shrunk toward leave-one-out panel prediction
    let th_r = upg::ridge_theta(&gr.c, g, &s, alpha)?;
    let e_loo = upg::loo_resid(&xt, nt, k, &th_r);
    let wsh: Vec<f64> = (0..k).map(|i| sig_r[i] / (sig_r[i] + 1.0)).collect();
    let tr_u = zt.mat(ltr, &tests, &u);
    let dr_u = zt.mat(lin, &tests, &u);
    let tr_s = zt.mat(ltr, &tests, &s);
    let dr_s = zt.mat(lin, &tests, &s);
    let a_score = |j: usize, v: f64| -> f64 { if best_a.1 < 0.0 { v.abs() } else if st[j].0 < best_a.2 { f64::NEG_INFINITY } else { (st[j].1 * v).abs() / st[j].2.powf(best_a.1) } };
    let b_score = |j: usize, v: f64| -> f64 { v.abs() / (st[j].0.max(0.0).powi(2) + best_b.1).powf(best_b.2) };
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.dex.tsv"))?);
    writeln!(w, "fold\tpert\tmethod\tn_de\tncall\ttp\tdisc\tspa")?;
    let rows: Vec<String> = (0..nt).into_par_iter().flat_map(|i| {
        let (t, d) = (&tr_u[i * nu..(i + 1) * nu], &dr_u[i * nu..(i + 1) * nu]);
        let nde = t.iter().filter(|v| v.abs() > DE_Z).count();
        let nc = d.iter().filter(|v| v.abs() > DE_Z).count();
        let mut o = vec![];
        let mut push = |name: &str, sc: Vec<f64>, val: Vec<f64>, rk: Vec<f64>| {
            let (tp, disc) = calls(&sc, &val, t, nc);
            o.push(format!("{f}\t{}\t{name}\t{nde}\t{nc}\t{tp}\t{disc}\t{:.6}", zt.perts[tests[i]], spa_de(&rk, t)));
        };
        let row = |p: &[f64]| p[i * nu..(i + 1) * nu].to_vec();
        let (r0, r2, rc, rr) = (row(&p_raw), row(&p_a2), row(&p_cor), row(&p_rel));
        push("direct_out", d.iter().map(|v| v.abs()).collect(), d.to_vec(), d.to_vec());
        push("baseline", r0.iter().map(|v| v.abs()).collect(), r0.clone(), r0.clone());
        push("A", (0..nu).map(|j| a_score(j, r2[j])).collect(), r2.clone(), r2.clone());
        push("B", r0.iter().map(|v| v.abs()).collect(), r0.clone(), (0..nu).map(|j| b_score(j, r0[j])).collect());
        push("ABC_coral", (0..nu).map(|j| a_score(j, rc[j])).collect(), rc.clone(), (0..nu).map(|j| b_score(j, rc[j])).collect());
        push("ABC_rel", (0..nu).map(|j| a_score(j, rr[j])).collect(), rr.clone(), (0..nu).map(|j| b_score(j, rr[j])).collect());
        push("AB", (0..nu).map(|j| a_score(j, r2[j])).collect(), r2.clone(), (0..nu).map(|j| b_score(j, r2[j])).collect());
        // D: all genes = panel (measured, shrunk) + imputed (A+B, no alignment); calls by |value| / sd, size by |calibrated value|
        let (ts, ds) = (&tr_s[i * k..(i + 1) * k], &dr_s[i * k..(i + 1) * k]);
        let xs = &xt[i * k..(i + 1) * k];
        let sh: Vec<f64> = (0..k).map(|q| wsh[q] * xs[q] + (1.0 - wsh[q]) * (xs[q] - e_loo[i * k + q])).collect();
        let tall: Vec<f64> = t.iter().chain(ts.iter()).copied().collect();
        let dall: Vec<f64> = d.iter().chain(ds.iter()).copied().collect();
        let cal: Vec<f64> = (0..nu).map(|j| st[j].1 * r2[j]).collect();
        let vall: Vec<f64> = cal.iter().chain(sh.iter()).copied().collect();
        let scall: Vec<f64> = (0..nu).map(|j| a_score(j, r2[j])).chain((0..k).map(|q| sh[q].abs() / wsh[q].max(0.05).sqrt())).collect();
        let ndeall = tall.iter().filter(|v| v.abs() > DE_Z).count();
        let ncall = dall.iter().filter(|v| v.abs() > DE_Z).count();
        let (tp, disc) = calls(&scall, &vall, &tall, ncall);
        o.push(format!("{f}\t{}\thybrid_all\t{ndeall}\t{ncall}\t{tp}\t{disc}\t{:.6}", zt.perts[tests[i]], spa_de(&vall, &tall)));
        for (q, (n, tp)) in calls_bins(&scall, &vall, &tall, ncall).iter().enumerate() { o.push(format!("{f}\t{}\thybrid_all_z{}\t{n}\t0\t{tp}\t0\tNaN", zt.perts[tests[i]], ["4-6", "6-10", "10+"][q])); }
        for (q, (n, tp)) in calls_bins(&r0.iter().map(|v| v.abs()).collect::<Vec<_>>(), &r0, t, nc).iter().enumerate() { o.push(format!("{f}\t{}\tbaseline_z{}\t{n}\t0\t{tp}\t0\tNaN", zt.perts[tests[i]], ["4-6", "6-10", "10+"][q])); }
        let (tp, disc) = calls(&dall.iter().map(|v| v.abs()).collect::<Vec<_>>(), &dall, &tall, ncall);
        o.push(format!("{f}\t{}\tdirect_all\t{ndeall}\t{ncall}\t{tp}\t{disc}\t{:.6}", zt.perts[tests[i]], spa_de(&dall, &tall)));
        o
    }).collect();
    for r in rows { writeln!(w, "{r}")?; }
    writeln!(wl, "done {:.1}s", t0.elapsed().as_secs_f64())?;
    Ok(())
}

/// summary of {pre}.dex.tsv: recall with correct sign, opposite-sign rate, median size rank
pub fn summary(a: &[String]) -> Result<()> {
    let pre = req(a, "--in")?;
    let txt = std::fs::read_to_string(format!("{pre}.dex.tsv"))?;
    let mut agg: Vec<(String, f64, f64, f64, f64, Vec<f64>)> = vec![];
    for l in txt.lines().skip(1) {
        let f: Vec<&str> = l.split('\t').collect();
        let m = f[2].to_string();
        let pos = match agg.iter().position(|x| x.0 == m) { Some(p) => p, None => { agg.push((m.clone(), 0.0, 0.0, 0.0, 0.0, vec![])); agg.len() - 1 } };
        let e = &mut agg[pos];
        e.1 += f[3].parse::<f64>()?; e.2 += f[4].parse::<f64>()?; e.3 += f[5].parse::<f64>()?; e.4 += f[6].parse::<f64>()?;
        let s: f64 = f[7].parse().unwrap_or(f64::NAN);
        if s.is_finite() { e.5.push(s); }
    }
    println!("method       recall_correct_sign  opposite_sign_rate  size_rank_r  (n_rank)");
    for (m, nde, nc, tp, disc, sp) in agg { println!("{m:12} {:.4}               {:.4}              {:.4}       ({})", tp / nde, disc / nc.max(1.0), median(sp.clone()), sp.len()); }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    /// sym_pow(A, 1/2) squared gives A, and sym_pow(A, -1/2) whitens A
    #[test]
    fn sym_pow_square_root() {
        let n = 6;
        let mut st = 3u64;
        let x: Vec<f64> = (0..n * 20).map(|_| (sm64(&mut st) >> 11) as f64 / (1u64 << 53) as f64 - 0.5).collect();
        let a = la::gemm(true, false, n, n, 20, &x, &x);
        let h = sym_pow(&a, n, 0.5).unwrap();
        let back = la::gemm(false, false, n, n, n, &h, &h);
        for (p, q) in back.iter().zip(&a) { assert!((p - q).abs() < 1e-9); }
        let wh = sym_pow(&a, n, -0.5).unwrap();
        let id = la::gemm(false, false, n, n, n, &la::gemm(false, false, n, n, n, &wh, &a), &wh);
        for i in 0..n { for j in 0..n { assert!((id[i * n + j] - if i == j { 1.0 } else { 0.0 }).abs() < 1e-8); } }
    }
    /// matched calls count correct-sign DE genes and sign errors among the top-scored genes
    #[test]
    fn calls_counts() {
        let score = [0.9, 0.1, 0.8, 0.7];
        let val = [5.0, 1.0, -3.0, 2.0];
        let t = [6.0, 0.0, 4.5, -5.0];
        assert_eq!(calls(&score, &val, &t, 3), (1, 2));
    }
}
