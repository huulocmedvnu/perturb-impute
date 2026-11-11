//! Additional analyses (Additional file 1: Table S1, Note S1): alternatives to the all-cell ridge, tested per fold.
//!   panelimp explore --stage test --z PREFIX --folds F --hubs H --out PREFIX [--fold 0] [--tss FILE] [--umi FILE]
//! H1 cross-half ridge, H2 identity prior (target gene and close neighbours), H3 local experts,
//! H4 quadratic terms. Everything is fitted on the training knockdowns; choices on the inner holdout.
use crate::causal::*;
use crate::explore::*;
use crate::*;

/// cross-half ridge: predict the other half of the cells from the panel genes of one half.
/// numerator (nu x k) = mean cross-half covariance of other genes with panel genes (noise cancels),
/// denominator (k x k) = mean within-half covariance of the panel genes (signal and noise of an input half)
struct Cross { num: Vec<f32>, den: Vec<f64> }
fn cross_stats(d: &Data, tr: &[usize], s: &[usize], u: &[usize], wts: Option<&[f32]>) -> Cross {
    let (g, k, nu, n) = (d.genes.len(), s.len(), u.len(), tr.len());
    let (z1, z2) = (d.rows(1, tr), d.rows(2, tr));
    let (mut a1, mut a2) = (cols(&z1, g, s), cols(&z2, g, s));
    let (b1, b2) = (cols(&z1, g, u), cols(&z2, g, u));
    drop((z1, z2));
    if let Some(w) = wts { for (r, &wi) in a1.chunks_mut(k).zip(w) { for v in r { *v *= wi; } } for (r, &wi) in a2.chunks_mut(k).zip(w) { for v in r { *v *= wi; } } }
    let unw = |m: &[f32]| -> Vec<f32> { match wts { Some(w) => m.chunks(k).zip(w).flat_map(|(r, &wi)| r.iter().map(move |v| if wi > 0.0 { v / wi } else { 0.0 })).collect(), None => m.to_vec() } };
    // weights enter once: weighted panel block against unweighted rows
    let (a1u, a2u) = (unw(&a1), unw(&a2));
    let mut num = sgemm(true, false, nu, k, n, &b2, &a1);
    let n2 = sgemm(true, false, nu, k, n, &b1, &a2);
    num.par_iter_mut().zip(&n2).for_each(|(x, y)| *x = 0.5 * (*x + *y));
    let (d1, d2) = (sgemm(true, false, k, k, n, &a1u, &a1), sgemm(true, false, k, k, n, &a2u, &a2));
    let den: Vec<f64> = (0..k * k).map(|x| { let (i, j) = (x / k, x % k); 0.25 * (d1[x] + d1[j * k + i] + d2[x] + d2[j * k + i]) as f64 }).collect();
    Cross { num, den }
}
fn cross_w(c: &Cross, k: usize, lambda: f64) -> Result<Vec<f32>> {
    let nu = c.num.len() / k;
    let mut a = c.den.clone();
    let tr: f64 = (0..k).map(|i| a[i * k + i]).sum::<f64>() / k as f64;
    for i in 0..k { a[i * k + i] += lambda * tr; }
    let mut bt: Vec<f64> = c.num.iter().map(|&v| v as f64).collect();
    la::posv(a, k, &mut bt, nu)?;
    Ok(bt.iter().map(|&v| v as f32).collect())
}
fn mean_r(p: &[f32], t: &[f32], nu: usize) -> f64 {
    let v: Vec<f64> = p.par_chunks(nu).zip(t.par_chunks(nu)).map(|(x, y)| pearson(&x.iter().map(|&v| v as f64).collect::<Vec<_>>(), &y.iter().map(|&v| v as f64).collect::<Vec<_>>())).filter(|v| v.is_finite()).collect();
    v.iter().sum::<f64>() / v.len().max(1) as f64
}
/// k-means on unit-norm rows (n x q), deterministic start from evenly spaced rows
fn kmeans(x: &[f32], q: usize, kc: usize, iters: usize) -> Vec<f32> {
    let n = x.len() / q;
    let mut cen: Vec<f32> = (0..kc).flat_map(|c| x[(c * n / kc) * q..(c * n / kc + 1) * q].to_vec()).collect();
    for _ in 0..iters {
        let lab: Vec<usize> = x.par_chunks(q).map(|r| nearest(r, &cen, q)).collect();
        let mut sum = vec![0.0f32; kc * q]; let mut cnt = vec![0usize; kc];
        for (r, &l) in x.chunks(q).zip(&lab) { cnt[l] += 1; for (s, v) in sum[l * q..(l + 1) * q].iter_mut().zip(r) { *s += v; } }
        for c in 0..kc { if cnt[c] > 0 { let nr = sum[c * q..(c + 1) * q].iter().map(|v| v * v).sum::<f32>().sqrt().max(1e-12); for (o, s) in cen[c * q..(c + 1) * q].iter_mut().zip(&sum[c * q..(c + 1) * q]) { *o = s / nr; } } }
    }
    cen
}
fn nearest(r: &[f32], cen: &[f32], q: usize) -> usize { (0..cen.len() / q).map(|c| (cen[c * q..(c + 1) * q].iter().zip(r).map(|(a, b)| a * b).sum::<f32>(), c)).fold((f32::NEG_INFINITY, 0), |a, b| if b.0 > a.0 { b } else { a }).1 }
fn unit_rows(m: &mut [f32], q: usize) { m.par_chunks_mut(q).for_each(|r| { let nr = r.iter().map(|v| v * v).sum::<f32>().sqrt().max(1e-12); for v in r { *v /= nr; } }); }

/// simplified all-gene readout: panel genes as measured in the input half, the other genes as predicted;
/// as many calls as direct measurement of all genes makes (|z| > 4), ranked by absolute value
fn all_gene(pred_u: &[f32], x_p: &[f32], t_u: &[f32], t_p: &[f32], d_u: &[f32], nu: usize, k: usize) -> (f64, f64, f64) {
    let n = t_u.len() / nu;
    let rows: Vec<(usize, usize, usize, usize, f64)> = (0..n).into_par_iter().map(|i| {
        let cat = |a: &[f32], b: &[f32]| -> Vec<f64> { a[i * nu..(i + 1) * nu].iter().chain(b[i * k..(i + 1) * k].iter()).map(|&v| v as f64).collect() };
        let (v, t, dd) = (cat(pred_u, x_p), cat(t_u, t_p), cat(d_u, x_p));
        let nde = t.iter().filter(|x| x.abs() > DE_Z).count();
        let nc = dd.iter().filter(|x| x.abs() > DE_Z).count();
        let mut o: Vec<usize> = (0..v.len()).collect();
        o.sort_by(|&a, &b| v[b].abs().partial_cmp(&v[a].abs()).unwrap().then(a.cmp(&b)));
        let (mut tp, mut disc) = (0, 0);
        for &j in &o[..nc.min(v.len())] { let same = v[j].signum() == t[j].signum(); if !same { disc += 1; } if same && t[j].abs() > DE_Z { tp += 1; } }
        (nde, nc, tp, disc, spa_de(&v, &t))
    }).collect();
    let s = |f: &dyn Fn(&(usize, usize, usize, usize, f64)) -> usize| rows.iter().map(f).sum::<usize>() as f64;
    (s(&|r| r.2) / s(&|r| r.0), s(&|r| r.3) / s(&|r| r.1), median(rows.iter().map(|r| r.4).collect()))
}

pub fn run(a: &[String]) -> Result<()> {
    let t0 = std::time::Instant::now();
    let c = setup(a)?;
    let (d, s, u) = (&c.d, &c.s, &c.u);
    let (g, nu, k) = (d.genes.len(), u.len(), s.len());
    let out = req(a, "--out")?;
    let mut log = BufWriter::new(std::fs::File::create(format!("{out}.explore.test.log"))?);
    macro_rules! say { ($($t:tt)*) => {{ let s = format!($($t)*); println!("{s}"); writeln!(log, "{s}")?; }} }
    let gi: HashMap<&str, usize> = d.genes.iter().enumerate().map(|(i, x)| (x.as_str(), i)).collect();
    let upos: HashMap<usize, usize> = u.iter().enumerate().map(|(q, &j)| (j, q)).collect();
    let mut tss: HashMap<String, (String, f64)> = HashMap::new();
    if let Some(f) = arg(a, "--tss") { for (name, v) in read_tss(&f)? { tss.entry(name).or_insert(v); } }
    // genes outside the panel with a TSS within 2 kb of the target of each perturbation
    let near_n = |name: &str| -> Vec<usize> { match tss.get(name) { Some((chr, pos)) => u.iter().enumerate().filter(|(_, &j)| d.genes[j] != name && tss.get(&d.genes[j]).map(|(c2, p2)| c2 == chr && (p2 - pos).abs() < 2e3).unwrap_or(false)).map(|(q, _)| q).collect(), None => vec![] } };
    let own_n = |name: &str| -> Option<usize> { gi.get(name).and_then(|j| upos.get(j)).copied() };
    let near = |p: usize| near_n(&d.perts[p]);
    let own = |p: usize| own_n(&d.perts[p]);
    // test perturbations: the test fold of the reference screen, or every perturbation of another screen (--target)
    let target = match arg(a, "--target") { Some(t) => { let x = Data::open(&t)?; if x.genes != d.genes { bail!("target genes differ from the reference") } Some(x) } None => None };
    let dt: &Data = target.as_ref().unwrap_or(d);
    let tests: Vec<usize> = if target.is_some() { (0..dt.perts.len()).collect() } else { c.tests.clone() };
    let tnames: Vec<String> = tests.iter().map(|&p| dt.perts[p].clone()).collect();
    let inames: Vec<String> = c.iho.iter().map(|&p| d.perts[p].clone()).collect();

    // one training set -> all fitted pieces; evaluation rows -> predictions of every method
    let umi: HashMap<String, f64> = match arg(a, "--umi") { Some(f) => std::fs::read_to_string(f)?.lines().filter_map(|l| { let x: Vec<&str> = l.split('\t').collect(); x.get(1).and_then(|v| v.parse::<f64>().ok()).map(|v| (x[0].to_string(), v)) }).collect(), None => HashMap::new() };
    struct Fitted { w0: Vec<f32>, cr: Cross, pc: Vec<f32>, own_z: f32, near_z: f32, cuts: Vec<f64>, own_bin: Vec<f32> }
    let fit = |tr: &[usize]| -> Result<Fitted> {
        let c0 = sgram(&d.rows(0, tr), tr.len(), g);
        let w0 = ridge32(&c0, g, s, u, c.alpha)?;
        drop(c0);
        let cr = cross_stats(d, tr, s, u, None);
        let den32: Vec<f32> = cr.den.iter().map(|&v| v as f32).collect();
        let pc = top_eig32(&den32, k, 20, 30)?; // k x 20 leading components of the panel block
        // identity prior: median z of the target gene, and of genes within 2 kb, in the held-out half of the training knockdowns
        let z2 = d.rows(2, tr);
        let (mut zo, mut zn, mut ze) = (vec![], vec![], vec![]);
        for (b, &p) in tr.iter().enumerate() {
            if let Some(q) = own(p) { zo.push(z2[b * g + u[q]] as f64); if let Some(&f) = umi.get(&d.perts[p]) { ze.push((f, z2[b * g + u[q]] as f64)); } }
            for q in near(p) { zn.push(z2[b * g + u[q]] as f64); }
        }
        // deciles of the control-cell UMI fraction of the target gene (expression level), median z in each
        ze.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
        let (mut cuts, mut own_bin) = (vec![], vec![]);
        if ze.len() >= 100 { for ch in ze.chunks(ze.len().div_ceil(10)) { cuts.push(ch.last().unwrap().0); own_bin.push(median(ch.iter().map(|x| x.1).collect()) as f32); } }
        Ok(Fitted { w0, cr, pc, own_z: median(zo) as f32, near_z: if zn.is_empty() { 0.0 } else { median(zn) as f32 }, cuts, own_bin })
    };
    let ident = |pred: &mut [f32], names: &[String], f: &Fitted, with_near: bool, by_expr: bool| { for (b, p) in names.iter().enumerate() {
        if let Some(q) = own_n(p) { pred[b * nu + q] = match (by_expr && !f.cuts.is_empty(), umi.get(p)) { (true, Some(&fr)) => f.own_bin[f.cuts.iter().position(|&c| fr <= c).unwrap_or(f.cuts.len() - 1)], _ => f.own_z }; }
        if with_near { for q in near_n(p) { pred[b * nu + q] = f.near_z; } } } };
    let comps = |x: &[f32], f: &Fitted| -> Vec<f32> { sgemm(false, false, x.len() / k, 20, k, x, &f.pc) };
    let quad = |x: &[f32], f: &Fitted, sd: &[f32]| -> Vec<f32> { let t = comps(x, f); t.par_chunks(20).flat_map_iter(|r| { let mut o = Vec::with_capacity(210); for i in 0..20 { for j in i..20 { o.push(r[i] * r[j] / (sd[i] * sd[j])); } } o }).collect() };

    // ---------------- choices on the inner holdout (models from the inner training set)
    let fi = fit(&c.itr)?;
    let (xh, th) = (cols(&d.rows(1, &c.iho), g, s), cols(&d.rows(2, &c.iho), g, u));
    let nh = c.iho.len();
    let p0 = sgemm(false, true, nh, nu, k, &xh, &fi.w0);
    say!("inner holdout, mean correlation per perturbation: all-cell ridge {:.4}", mean_r(&p0, &th, nu));
    let lambdas = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0];
    let mut b1 = (f64::NEG_INFINITY, 0.0);
    for &l in &lambdas { let r = mean_r(&sgemm(false, true, nh, nu, k, &xh, &cross_w(&fi.cr, k, l)?), &th, nu); say!("  H1 cross-half ridge, lambda {l}: {r:.4}"); if r > b1.0 { b1 = (r, l); } }
    let w1 = cross_w(&fi.cr, k, b1.1)?;
    let p1 = sgemm(false, true, nh, nu, k, &xh, &w1);
    // H2: does the identity prior lower the squared error on the inner holdout
    let sq = |p: &[f32]| -> f64 { p.par_iter().zip(th.par_iter()).map(|(x, y)| ((*x - *y) as f64).powi(2)).sum::<f64>() };
    let (mut p_own, mut p_near) = (p1.clone(), p1.clone());
    ident(&mut p_own, &inames, &fi, false, false); ident(&mut p_near, &inames, &fi, true, false);
    let use_near = sq(&p_near) < sq(&p_own);
    let mut p_expr = p1.clone(); ident(&mut p_expr, &inames, &fi, use_near, true);
    let mut p_cst = p1.clone(); ident(&mut p_cst, &inames, &fi, use_near, false);
    let use_expr = sq(&p_expr) < sq(&p_cst);
    say!("  H2 identity prior: training median z of the target gene {:.2}, of genes within 2 kb {:.2}; squared error H1 {:.0}, + target {:.0}, + target and neighbours {:.0} -> neighbours {}; by expression decile {:.0} -> {} (medians {:?})", fi.own_z, fi.near_z, sq(&p1), sq(&p_own), sq(&p_near), if use_near { "kept" } else { "dropped" }, sq(&p_expr), if use_expr { "kept" } else { "dropped" }, fi.own_bin.iter().map(|v| (v * 10.0).round() / 10.0).collect::<Vec<_>>());
    // H3: local experts (cross-half ridge per cluster of training knockdowns, blended with the global statistics)
    let x_tr1 = cols(&d.rows(1, &c.itr), g, s);
    let mut t_tr = comps(&x_tr1, &fi); unit_rows(&mut t_tr, 20);
    let mut t_ho = comps(&xh, &fi); unit_rows(&mut t_ho, 20);
    let mut b3 = (b1.0, 0usize, 0.0f32);
    for kc in [4usize, 8] { for omega in [0.5f32, 0.8] {
        let cen = kmeans(&t_tr, 20, kc, 15);
        let lab_tr: Vec<usize> = t_tr.par_chunks(20).map(|r| nearest(r, &cen, 20)).collect();
        let lab_ho: Vec<usize> = t_ho.par_chunks(20).map(|r| nearest(r, &cen, 20)).collect();
        let mut p3 = p1.clone();
        for cl in 0..kc {
            let rows_ho: Vec<usize> = (0..nh).filter(|&b| lab_ho[b] == cl).collect();
            let n_cl = lab_tr.iter().filter(|&&l| l == cl).count();
            if rows_ho.is_empty() || n_cl < 50 { continue; }
            // sample weights: omega on every knockdown plus (1 - omega) concentrated on the cluster
            let wts: Vec<f32> = lab_tr.iter().map(|&l| omega + if l == cl { (1.0 - omega) * c.itr.len() as f32 / n_cl as f32 } else { 0.0 }).collect();
            let wl = cross_w(&cross_stats(d, &c.itr, s, u, Some(&wts)), k, b1.1)?;
            let xs: Vec<f32> = rows_ho.iter().flat_map(|&b| xh[b * k..(b + 1) * k].to_vec()).collect();
            let pl = sgemm(false, true, rows_ho.len(), nu, k, &xs, &wl);
            for (i, &b) in rows_ho.iter().enumerate() { p3[b * nu..(b + 1) * nu].copy_from_slice(&pl[i * nu..(i + 1) * nu]); }
        }
        let r = mean_r(&p3, &th, nu);
        say!("  H3 local experts, {kc} clusters, global share {omega}: {r:.4}");
        if r > b3.0 { b3 = (r, kc, omega); }
    } }
    // H4: quadratic terms of the 20 leading panel components, fitted to the cross-half residual of the training knockdowns
    let sd: Vec<f32> = { let t = comps(&x_tr1, &fi); (0..20).map(|i| (t.chunks(20).map(|r| r[i] * r[i]).sum::<f32>() / (t.len() / 20) as f32).sqrt().max(1e-6)).collect() };
    let quad_fit = |tr: &[usize], f: &Fitted, w: &[f32], mu: f64| -> Result<Vec<f32>> {
        let (x1, t2) = (cols(&d.rows(1, tr), g, s), cols(&d.rows(2, tr), g, u));
        let mut res = t2; let p = sgemm(false, true, tr.len(), nu, k, &x1, w);
        res.par_iter_mut().zip(&p).for_each(|(a, b)| *a -= b);
        let q = quad(&x1, f, &sd);
        let mut qq: Vec<f64> = sgemm(true, false, 210, 210, tr.len(), &q, &q).iter().map(|&v| v as f64).collect();
        let trc: f64 = (0..210).map(|i| qq[i * 211]).sum::<f64>() / 210.0;
        for i in 0..210 { qq[i * 211] += mu * trc; }
        let mut rhs: Vec<f64> = sgemm(true, false, nu, 210, tr.len(), &res, &q).iter().map(|&v| v as f64).collect();
        la::posv(qq, 210, &mut rhs, nu)?;
        Ok(rhs.iter().map(|&v| v as f32).collect())
    };
    let mut b4 = (b1.0, f64::NAN);
    let qh = quad(&xh, &fi, &sd);
    for mu in [0.1, 1.0, 10.0, 100.0] {
        let gm = quad_fit(&c.itr, &fi, &w1, mu)?;
        let r = mean_r(&add(&p1, &[(&sgemm(false, true, nh, nu, 210, &qh, &gm), 1.0)]), &th, nu);
        say!("  H4 quadratic terms, penalty {mu}: {r:.4}");
        if r > b4.0 { b4 = (r, mu); }
    }
    say!("chosen: H1 lambda {} ({:.4}); H3 {} ; H4 {}; resident {:.0} MB at {:.1}s", b1.1, b1.0, if b3.1 > 0 { format!("{} clusters, global share {} ({:.4})", b3.1, b3.2, b3.0) } else { "no gain on the inner holdout, not used".into() },
         if b4.1.is_finite() { format!("penalty {} ({:.4})", b4.1, b4.0) } else { "no gain on the inner holdout, not used".into() }, rss_mb(), t0.elapsed().as_secs_f64());
    drop((fi, x_tr1, t_tr, p0, p1, p_own, p_near, p_expr, p_cst));

    // ---------------- final models on all references, fold test set
    let ff = fit(&c.refs)?;
    let nt = tests.len();
    let (z1, z2) = (dt.rows(1, &tests), dt.rows(2, &tests));
    let (xt, tt, dd, tp_) = (cols(&z1, g, s), cols(&z2, g, u), cols(&z1, g, u), cols(&z2, g, s));
    // evaluation frame of the manuscript (all-cell Gram of the references)
    let (basis, den) = {
        let c0 = sgram(&d.rows(0, &c.refs), c.refs.len(), g);
        let v160 = top_eig32(&c0, g, 160, 30)?;
        let v10: Vec<f64> = (0..g).flat_map(|j| v160[j * 160..j * 160 + NPC].iter().map(|&v| v as f64).collect::<Vec<_>>()).collect();
        (basis_on(&v10, NPC, u), cols(&sgemm(false, true, nt, g, 160, &sgemm(false, false, nt, 160, g, &z1, &v160), &v160), g, u))
    };
    drop((z1, z2));
    let mut preds: Vec<(String, Vec<f32>, f64)> = vec![];
    let tm = std::time::Instant::now();
    let r0 = sgemm(false, true, nt, nu, k, &xt, &ff.w0);
    preds.push(("ridge (all cells), baseline".into(), r0, 1e3 * tm.elapsed().as_secs_f64()));
    preds.push(("direct (input half)".into(), dd.clone(), 0.0));
    preds.push(("whole half on 160 programs".into(), den.clone(), 0.0));
    let tm = std::time::Instant::now();
    let w1 = cross_w(&ff.cr, k, b1.1)?;
    let h1 = sgemm(false, true, nt, nu, k, &xt, &w1);
    preds.push(("H1 cross-half ridge".into(), h1.clone(), 1e3 * tm.elapsed().as_secs_f64()));
    let mut h2 = preds[0].1.clone(); ident(&mut h2, &tnames, &ff, false, false);
    preds.push(("H2a ridge + target gene".into(), h2, 0.0));
    let mut h2 = preds[0].1.clone(); ident(&mut h2, &tnames, &ff, use_near, false);
    preds.push(("H2b + neighbours within 2 kb".into(), h2, 0.0));
    let mut h2 = preds[0].1.clone(); ident(&mut h2, &tnames, &ff, use_near, use_expr);
    preds.push(("H2c + expression deciles".into(), h2, 0.0));
    let mut best_lin = h1.clone();
    if b3.1 > 0 {
        let tm = std::time::Instant::now();
        let x_tr1 = cols(&d.rows(1, &c.refs), g, s);
        let mut t_tr = comps(&x_tr1, &ff); unit_rows(&mut t_tr, 20);
        let mut t_te = comps(&xt, &ff); unit_rows(&mut t_te, 20);
        let cen = kmeans(&t_tr, 20, b3.1, 15);
        let lab_tr: Vec<usize> = t_tr.par_chunks(20).map(|r| nearest(r, &cen, 20)).collect();
        let lab_te: Vec<usize> = t_te.par_chunks(20).map(|r| nearest(r, &cen, 20)).collect();
        let mut h3 = h1.clone();
        for cl in 0..b3.1 {
            let rows: Vec<usize> = (0..nt).filter(|&b| lab_te[b] == cl).collect();
            let n_cl = lab_tr.iter().filter(|&&l| l == cl).count();
            if rows.is_empty() || n_cl < 50 { continue; }
            let wts: Vec<f32> = lab_tr.iter().map(|&l| b3.2 + if l == cl { (1.0 - b3.2) * c.refs.len() as f32 / n_cl as f32 } else { 0.0 }).collect();
            let wl = cross_w(&cross_stats(d, &c.refs, s, u, Some(&wts)), k, b1.1)?;
            let xs: Vec<f32> = rows.iter().flat_map(|&b| xt[b * k..(b + 1) * k].to_vec()).collect();
            let pl = sgemm(false, true, rows.len(), nu, k, &xs, &wl);
            for (i, &b) in rows.iter().enumerate() { h3[b * nu..(b + 1) * nu].copy_from_slice(&pl[i * nu..(i + 1) * nu]); }
        }
        preds.push(("H3 local experts".into(), h3.clone(), 1e3 * tm.elapsed().as_secs_f64()));
        best_lin = h3;
    }
    if b4.1.is_finite() {
        let tm = std::time::Instant::now();
        let gm = quad_fit(&c.refs, &ff, &w1, b4.1)?;
        let h4 = add(&h1, &[(&sgemm(false, true, nt, nu, 210, &quad(&xt, &ff, &sd), &gm), 1.0)]);
        preds.push(("H4 cross-half + quadratic".into(), h4, 1e3 * tm.elapsed().as_secs_f64()));
    }
    let mut h5 = best_lin; ident(&mut h5, &tnames, &ff, use_near, use_expr);
    preds.push(("H5 H1/H3 + identity".into(), h5, 0.0));

    say!("\ntest set{}: {nt} perturbations.", arg(a, "--target").map(|t| format!(" {t}")).unwrap_or_default()); say!(" Genes outside the panel (as in the manuscript), then a simplified all-gene readout (panel genes as measured).");
    say!("{:30} {:>6} {:>8} {:>7} {:>7} {:>7} {:>8} {:>9} | {:>7} {:>7} {:>8} | {:>8}", "method", "n_resp", "recovery", "spec_r", "recall", "oppos.", "rank_med", "rank_mean", "recall", "oppos.", "rank_med", "fit_ms");
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.explore.test.tsv"))?);
    writeln!(w, "method\tn_responding\trecovery\tspecific_r\tde_recall\topposite_sign\tsize_rank_median\tsize_rank_mean\tall_gene_recall\tall_gene_opposite\tall_gene_size_rank\tfit_ms")?;
    let mut scores = vec![];
    for (name, p, ms) in &preds {
        let sc = score(p, &tt, &dd, &den, nu, &basis);
        let ag = all_gene(p, &xt, &tt, &tp_, &dd, nu, k);
        say!("{name:30} {:>6} {:>8.4} {:>7.4} {:>7.4} {:>7.4} {:>8.4} {:>9.4} | {:>7.4} {:>7.4} {:>8.4} | {:>8.0}", sc.n_resp, sc.recovery, sc.spec, sc.recall, sc.opp, sc.spa, sc.spa_mean, ag.0, ag.1, ag.2, ms);
        writeln!(w, "{name}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.1}", sc.n_resp, sc.recovery, sc.spec, sc.recall, sc.opp, sc.spa, sc.spa_mean, ag.0, ag.1, ag.2, ms)?;
        scores.push(sc);
    }
    say!("\npaired against the baseline: size rank [n, mean difference, share higher, Wilcoxon P] | specific r of responders [n, median difference, share higher, P]");
    for (i, (name, _, _)) in preds.iter().enumerate().skip(3) {
        let pair = |x: &[f64], y: &[f64]| -> Vec<f64> { x.iter().zip(y).filter(|(p, q)| p.is_finite() && q.is_finite()).map(|(p, q)| p - q).collect() };
        let (d1, d2) = (pair(&scores[i].spa_v, &scores[0].spa_v), pair(&scores[i].spec_v, &scores[0].spec_v));
        say!("{name:30} size rank {} {:+.4} {:.3} {:.2e} | specific r {} {:+.4} {:.3} {:.2e}", d1.len(), d1.iter().sum::<f64>() / d1.len() as f64, d1.iter().filter(|v| **v > 0.0).count() as f64 / d1.len() as f64, upg::wilcoxon(&d1),
             d2.len(), median(d2.clone()), d2.iter().filter(|v| **v > 0.0).count() as f64 / d2.len() as f64, upg::wilcoxon(&d2));
    }
    say!("\nresident {:.0} MB, total {:.1}s", rss_mb(), t0.elapsed().as_secs_f64());
    Ok(())
}
