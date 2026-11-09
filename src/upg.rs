//! Step 9 upgrades: specific-subspace greedy (U1a), spectral imputation with empirical-Bayes scaling
//! (U1b), sparse residual correction on leave-one-out panel residuals (U3) and a leave-one-out
//! out-of-reference score (U4). All choices use the inner holdout of the reference only.
use crate::*;

/// Low-rank greedy MMSE. Signal covariance Sigma = U U' (U: g x r, row-major), noise nz_j,
/// gain of gene j = trace reduction / w_j. Residual covariance kept as U A U' (A: r x r).
pub fn greedy_lowrank(u: &[f64], g: usize, r: usize, k: usize, nz: &[f64], w: &[f64]) -> Vec<usize> {
    let gm = la::gemm(true, false, r, r, g, u, u);
    let mut a = vec![0.0; r * r];
    for i in 0..r { a[i * r + i] = 1.0; }
    let mut used = vec![false; g];
    let mut sel = Vec::with_capacity(k);
    for _ in 0..k {
        let b = la::gemm(false, false, g, r, r, u, &a);
        let h = la::gemm(false, false, g, r, r, &b, &gm);
        let best = (0..g).into_par_iter().filter(|&j| !used[j] && nz[j].is_finite() && w[j] > 0.0).map(|j| {
            let (bj, hj, uj) = (&b[j * r..(j + 1) * r], &h[j * r..(j + 1) * r], &u[j * r..(j + 1) * r]);
            let num: f64 = hj.iter().zip(bj).map(|(x, y)| x * y).sum();
            let den: f64 = bj.iter().zip(uj).map(|(x, y)| x * y).sum::<f64>() + nz[j];
            (num / den / w[j], j)
        }).reduce(|| (f64::NEG_INFINITY, usize::MAX), |x, y| if y.0 > x.0 || (y.0 == x.0 && y.1 < x.1) { y } else { x });
        let j = best.1;
        if j == usize::MAX { break; }
        used[j] = true;
        sel.push(j);
        let bj = b[j * r..(j + 1) * r].to_vec();
        let den: f64 = bj.iter().zip(&u[j * r..(j + 1) * r]).map(|(x, y)| x * y).sum::<f64>() + nz[j];
        a.par_chunks_mut(r).enumerate().for_each(|(p, row)| { for q in 0..r { row[q] -= bj[p] * bj[q] / den; } });
    }
    sel
}

/// Reference spectrum: top nk eigenvectors (g x nk), signal eigenvalues (C/n - s2), diag(C)/n.
pub struct Spec { pub v: Vec<f64>, pub lam: Vec<f64>, pub nk: usize, pub s2: f64, pub dc: Vec<f64>, pub kmp: usize }
pub fn spec_from(c: &[f64], g: usize, n: usize, s2: f64, nk: usize) -> Result<Spec> {
    let (v, ev) = la::top_eig(c, g, nk)?;
    let edge = s2 * (1.0 + (g as f64 / n as f64).sqrt()).powi(2);
    let kmp = ev.iter().filter(|&&e| e / n as f64 > edge).count();
    let lam = ev.iter().map(|e| e / n as f64 - s2).collect();
    let dc = (0..g).map(|j| c[j * g + j] / n as f64).collect();
    Ok(Spec { v, lam, nk, s2, dc, kmp })
}

/// Factor of the specific signal covariance: programs lo..hi, columns scaled by sqrt(lambda).
pub fn spec_factor(sp: &Spec, g: usize, lo: usize, hi: usize) -> (Vec<f64>, usize) {
    let r = hi - lo;
    let u = (0..g).flat_map(|j| (lo..hi).map(move |k| (j, k))).map(|(j, k)| sp.v[j * sp.nk + k] * sp.lam[k].max(0.0).sqrt()).collect();
    (u, r)
}

/// Observed panel covariance (S x S) and target covariance (U x S) under the low-rank plus diagonal model.
fn spec_cov(sp: &Spec, s: &[usize], u: &[usize], kk: usize, tau: f64) -> (Vec<f64>, Vec<f64>) {
    let phi: Vec<f64> = (0..kk).map(|k| { let l = sp.lam[k].max(0.0); if l > 0.0 { l * l / (l + tau) } else { 0.0 } }).collect();
    let fac = |rows: &[usize], w: &dyn Fn(usize) -> f64| -> Vec<f64> { rows.iter().flat_map(|&j| (0..kk).map(move |k| (j, k))).map(|(j, k)| sp.v[j * sp.nk + k] * w(k)).collect() };
    let fs = fac(s, &|k| phi[k].sqrt());
    let mut css = la::gemm(false, true, s.len(), s.len(), kk, &fs, &fs);
    for (i, &j) in s.iter().enumerate() {
        let expl: f64 = (0..kk).map(|k| sp.v[j * sp.nk + k].powi(2) * sp.lam[k].max(0.0)).sum();
        css[i * s.len() + i] += (sp.dc[j] - sp.s2 - expl).max(0.0) + sp.s2;
    }
    let fu = fac(u, &|k| phi[k]);
    let fsr = fac(s, &|_| 1.0);
    let cus = la::gemm(false, true, u.len(), s.len(), kk, &fu, &fsr);
    (css, cus)
}

/// Spectral BLUP weights W' (|U| x |S|) and the panel precision Theta (S x S).
pub fn spec_w(sp: &Spec, s: &[usize], u: &[usize], kk: usize, tau: f64) -> Result<(Vec<f64>, Vec<f64>)> {
    let k = s.len();
    let (css, mut cus) = spec_cov(sp, s, u, kk, tau);
    la::posv(css.clone(), k, &mut cus, u.len())?;
    let mut th = vec![0.0; k * k];
    for i in 0..k { th[i * k + i] = 1.0; }
    la::posv(css, k, &mut th, k)?;
    Ok((cus, th))
}

/// Ridge precision Theta = (C_SS + alpha tr I)^-1 (same scaling as `ridge`).
pub(crate) fn ridge_theta(c: &[f64], g: usize, s: &[usize], alpha: f64) -> Result<Vec<f64>> {
    let k = s.len();
    let mut a = vec![0.0; k * k];
    for (i, &si) in s.iter().enumerate() { for (j, &sj) in s.iter().enumerate() { a[i * k + j] = c[si * g + sj]; } }
    let tr: f64 = (0..k).map(|i| a[i * k + i]).sum::<f64>() / k as f64;
    for i in 0..k { a[i * k + i] += alpha * tr; }
    let mut th = vec![0.0; k * k];
    for i in 0..k { th[i * k + i] = 1.0; }
    la::posv(a, k, &mut th, k)?;
    Ok(th)
}

/// Leave-one-gene-out residuals of panel rows X (n x k): e = (X Theta) / diag(Theta).
pub fn loo_resid(x: &[f64], n: usize, k: usize, th: &[f64]) -> Vec<f64> {
    let mut e = la::gemm(false, false, n, k, k, x, th);
    e.par_chunks_mut(k).for_each(|row| { for g in 0..k { row[g] /= th[g * k + g]; } });
    e
}

/// Sparse residual correction (CSR over target genes): beta of target residual on the top-K
/// leave-one-out panel residuals with |correlation| >= t.
pub struct Csr { pub ptr: Vec<usize>, pub idx: Vec<usize>, pub val: Vec<f64> }
pub fn sparse_fit(e: &[f64], r: &[f64], n: usize, k: usize, nu: usize, kk: usize, t: f64) -> Result<Csr> {
    let ete = la::gemm(true, false, k, k, n, e, e);
    let etr = la::gemm(true, false, k, nu, n, e, r);
    let en: Vec<f64> = (0..k).map(|g| ete[g * k + g].sqrt()).collect();
    let rn: Vec<f64> = (0..nu).map(|j| (0..n).map(|i| r[i * nu + j].powi(2)).sum::<f64>().sqrt()).collect();
    let rows: Vec<Vec<(usize, f64)>> = (0..nu).into_par_iter().map(|j| {
        if kk == 0 || rn[j] == 0.0 { return vec![]; }
        let mut cand: Vec<(usize, f64)> = (0..k).filter(|&g| en[g] > 0.0).map(|g| (g, etr[g * nu + j] / (en[g] * rn[j]))).filter(|x| x.1.abs() >= t).collect();
        cand.sort_by(|a, b| b.1.abs().partial_cmp(&a.1.abs()).unwrap().then(a.0.cmp(&b.0)));
        cand.truncate(kk);
        let m = cand.len();
        if m == 0 { return vec![]; }
        let mut a = vec![0.0; m * m];
        for (p, &(gp, _)) in cand.iter().enumerate() { for (q, &(gq, _)) in cand.iter().enumerate() { a[p * m + q] = ete[gp * k + gq]; } }
        let tr: f64 = (0..m).map(|p| a[p * m + p]).sum::<f64>() / m as f64;
        for p in 0..m { a[p * m + p] += 1e-3 * tr; }
        let mut b: Vec<f64> = cand.iter().map(|&(gp, _)| etr[gp * nu + j]).collect();
        if la::posv(a, m, &mut b, 1).is_err() { return vec![]; }
        cand.iter().zip(b).map(|(&(gp, _), v)| (gp, v)).collect()
    }).collect();
    let mut csr = Csr { ptr: vec![0], idx: vec![], val: vec![] };
    for row in rows { for (g, v) in row { csr.idx.push(g); csr.val.push(v); } csr.ptr.push(csr.idx.len()); }
    Ok(csr)
}
/// pred (n x nu) += E (n x k) * B' where B is the CSR (nu rows over k columns)
pub fn sparse_apply(pred: &mut [f64], e: &[f64], n: usize, k: usize, nu: usize, b: &Csr) {
    pred.par_chunks_mut(nu).enumerate().for_each(|(i, row)| {
        let ei = &e[i * k..(i + 1) * k];
        for j in 0..nu { let mut s = 0.0; for p in b.ptr[j]..b.ptr[j + 1] { s += b.val[p] * ei[b.idx[p]]; } row[j] += s; }
    });
    let _ = n;
}

/// Empirical-Bayes per-gene scale from inner-holdout predictions (n x nu) and truth.
pub fn eb_scale(im: &[f64], tr: &[f64], n: usize, nu: usize) -> Vec<f64> {
    let st: Vec<(f64, f64)> = (0..nu).into_par_iter().map(|j| {
        let (mut sxx, mut sxy) = (0.0, 0.0);
        for i in 0..n { let (x, y) = (im[i * nu + j], tr[i * nu + j]); sxx += x * x; sxy += x * y; }
        if sxx <= 0.0 { return (1.0, f64::INFINITY); }
        let b = sxy / sxx;
        let rss: f64 = (0..n).map(|i| (tr[i * nu + j] - b * im[i * nu + j]).powi(2)).sum();
        (b, rss / (n as f64 - 1.0) / sxx)
    }).collect();
    let fin: Vec<&(f64, f64)> = st.iter().filter(|x| x.1.is_finite()).collect();
    let m = median(fin.iter().map(|x| x.0).collect());
    let mean_b = fin.iter().map(|x| x.0).sum::<f64>() / fin.len() as f64;
    let var_b = fin.iter().map(|x| (x.0 - mean_b).powi(2)).sum::<f64>() / fin.len() as f64;
    let mse = fin.iter().map(|x| x.1).sum::<f64>() / fin.len() as f64;
    let t2 = (var_b - mse).max(1e-6);
    st.iter().map(|&(b, se)| if se.is_finite() { (b / se + m / t2) / (1.0 / se + 1.0 / t2) } else { m }).collect()
}

fn ranks_mean(a: &[f64], b: &[f64]) -> Vec<f64> {
    let (ra, rb) = (ranks(a), ranks(b));
    ra.iter().zip(&rb).map(|(x, y)| (x + y) / 2.0).collect()
}

/// fold-level benchmark of baseline and upgrades; writes {out}.upg.tsv
#[allow(clippy::too_many_arguments)]
pub fn run(a: &[String], z: &Z, zt: &Z, fold: &[i64], group: &[i64], excl: &HashSet<(usize, i64)>, out: &str) -> Result<()> {
    let t0 = std::time::Instant::now();
    let g = z.genes.len();
    let transfer = !std::ptr::eq(z, zt);
    let runs: Vec<i64> = if transfer { vec![-9] } else { arg(a, "--folds-run").map(|s| s.split(',').map(|x| x.parse().unwrap()).collect()).unwrap_or(vec![0]) };
    let k = 500;
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.upg.tsv"))?);
    writeln!(w, "fold\tpanel\tmethod\tpert\tr\trperp\tden\tr_ceil\trperp_ceil\tn_de\tm_call\ttp\tdisc\tm_tp\tm_disc\tspa\tspa_dir\toor\tloo\tcomb")?;
    let mut wl = BufWriter::new(std::fs::File::create(format!("{out}.upg.log"))?);
    let key = |i: usize| -> u64 { let kk = if group[i] >= 0 { group[i] as u64 } else { 1_000_000 + i as u64 }; (kk.wrapping_mul(0x9E3779B97F4A7C15) >> 32) % 10 };
    for &f in &runs {
        let refs: Vec<usize> = (0..z.perts.len()).filter(|&i| transfer || (fold[i] != f && !excl.contains(&(i, f)))).collect();
        let tests: Vec<usize> = if transfer { (0..zt.perts.len()).collect() } else { (0..z.perts.len()).filter(|&i| fold[i] == f).collect() };
        let itr: Vec<usize> = refs.iter().copied().filter(|&i| key(i) != 0).collect();
        let iho: Vec<usize> = refs.iter().copied().filter(|&i| key(i) == 0).collect();
        let tg = std::time::Instant::now();
        let gr = grams(z, &refs, &itr)?;
        // noise variance of full-cell z from half differences of the reference
        let s2 = refs.par_iter().map(|&i| (0..g).map(|j| (z.at(1, i, j) - z.at(2, i, j)).powi(2)).sum::<f64>()).sum::<f64>() / (2.0 * (refs.len() * g) as f64);
        let sp = spec_from(&gr.c, g, refs.len(), s2, 500)?;
        let sp_in = spec_from(&gr.c_in, g, itr.len(), s2, 500)?;
        writeln!(wl, "fold {f}: refs {} itr {} iho {} tests {} s2 {s2:.4} kMP {} (inner {}) grams+eig {:.1}s", refs.len(), itr.len(), iho.len(), tests.len(), sp.kmp, sp_in.kmp, tg.elapsed().as_secs_f64())?;
        // panels
        let tp = std::time::Instant::now();
        let p_g = sim::greedy(&gr, refs.len(), g, k, &vec![1.0; g]);
        let t_g = tp.elapsed().as_secs_f64();
        let tp = std::time::Instant::now();
        let (uf, r) = spec_factor(&sp, g, NPC, sp.kmp);
        let p_s = greedy_lowrank(&uf, g, r, k, &vec![sp.s2; g], &vec![1.0; g]);
        writeln!(wl, "fold {f}: greedy full {t_g:.1}s, specific low-rank greedy (r = {r}) {:.2}s, overlap {}", tp.elapsed().as_secs_f64(), p_s.iter().filter(|j| p_g.contains(j)).count())?;
        let inp: HashSet<usize> = p_g.iter().chain(&p_s).copied().collect();
        let u: Vec<usize> = (0..g).filter(|j| !inp.contains(j)).collect();
        let nu = u.len();
        let basis = basis_on(&gr.v, NPC, &u);
        // denoised half-0 profiles (K = 160) on u, truth and direct
        let nk = *DEN_K.iter().max().unwrap();
        let h0full = zt.mat(1, &tests, &(0..g).collect::<Vec<_>>());
        let coef = la::gemm(false, false, tests.len(), nk, g, &h0full, &gr.vk);
        let ck: Vec<f64> = (0..tests.len()).flat_map(|i| coef[i * nk..i * nk + 160].to_vec()).collect();
        let vu: Vec<f64> = u.iter().flat_map(|&j| gr.vk[j * nk..j * nk + 160].to_vec()).collect();
        let den = la::gemm(false, true, tests.len(), nu, 160, &ck, &vu);
        let tr = zt.mat(2, &tests, &u);
        let dr = zt.mat(1, &tests, &u);
        let truth_ho = z.mat(2, &iho, &u);
        let nt = tests.len();
        let nh = iho.len();
        for (pname, s) in [("greedy", &p_g), ("specgreedy", &p_s)] {
            let tpn = std::time::Instant::now();
            let xh = z.mat(1, &iho, s);
            let xt = zt.mat(1, &tests, s);
            let mean_r = |im: &[f64], t: &[f64], n: usize| -> f64 { (0..n).map(|i| pearson(&im[i * nu..(i + 1) * nu], &t[i * nu..(i + 1) * nu])).filter(|x| x.is_finite()).sum::<f64>() / n as f64 };
            // baseline ridge
            let mut best = (f64::NEG_INFINITY, 0.0);
            for &al in &ALPHAS {
                let wt = ridge(&gr.c_in, g, s, &u, al)?;
                let im = la::gemm(false, true, nh, nu, k, &xh, &wt);
                let m = mean_r(&im, &truth_ho, nh);
                if m > best.0 { best = (m, al); }
            }
            let alpha = best.1;
            let wt_r = ridge(&gr.c, g, s, &u, alpha)?;
            let th_r = ridge_theta(&gr.c, g, s, alpha)?;
            let im_r = la::gemm(false, true, nt, nu, k, &xt, &wt_r);
            // spectral: choose (K, tau) on the inner holdout
            let kgrid = [sp_in.kmp, 300, 500];
            let mut bs = (f64::NEG_INFINITY, 0usize, 0.0);
            for &kk in &kgrid { for &tau in &[0.0, 0.5, 2.0, 8.0] {
                let (wt, _) = spec_w(&sp_in, s, &u, kk, tau)?;
                let im = la::gemm(false, true, nh, nu, k, &xh, &wt);
                let m = mean_r(&im, &truth_ho, nh);
                if m > bs.0 { bs = (m, kk, tau); }
            } }
            let (wt_in, th_in) = spec_w(&sp_in, s, &u, bs.1, bs.2)?;
            let im_ho = la::gemm(false, true, nh, nu, k, &xh, &wt_in);
            let ebs = eb_scale(&im_ho, &truth_ho, nh, nu);
            let kfin = if bs.1 == sp_in.kmp { sp.kmp } else { bs.1 };
            let (wt_s, th_s) = spec_w(&sp, s, &u, kfin, bs.2)?;
            let mut im_s = la::gemm(false, true, nt, nu, k, &xt, &wt_s);
            im_s.par_chunks_mut(nu).for_each(|row| for j in 0..nu { row[j] *= ebs[j]; });
            // stage 2: residual correction chosen on the inner holdout (size ranking among DE genes)
            let xi = z.mat(0, &itr, s);
            let yi = z.mat(0, &itr, &u);
            let mut ri = la::gemm(false, true, itr.len(), nu, k, &xi, &wt_in);
            ri.par_iter_mut().zip(&yi).for_each(|(p, y)| *p = y - *p);
            let ei = loo_resid(&xi, itr.len(), k, &th_in);
            let eh = loo_resid(&xh, nh, k, &th_in);
            let base_spa = median((0..nh).map(|i| spa_de(&im_ho[i * nu..(i + 1) * nu], &truth_ho[i * nu..(i + 1) * nu])).collect());
            let mut b2 = (base_spa, 0usize, 0.0);
            for &kk in &[3usize, 5, 10] { for &t in &[0.05, 0.1, 0.2] {
                let csr = sparse_fit(&ei, &ri, itr.len(), k, nu, kk, t)?;
                let mut p = im_ho.clone();
                sparse_apply(&mut p, &eh, nh, k, nu, &csr);
                let m = median((0..nh).map(|i| spa_de(&p[i * nu..(i + 1) * nu], &truth_ho[i * nu..(i + 1) * nu])).collect());
                if m > b2.0 { b2 = (m, kk, t); }
            } }
            let mut im_c = im_s.clone();
            if b2.1 > 0 {
                let xr = z.mat(0, &refs, s);
                let yr = z.mat(0, &refs, &u);
                let mut rr = la::gemm(false, true, refs.len(), nu, k, &xr, &wt_s);
                rr.par_iter_mut().zip(&yr).for_each(|(p, y)| *p = y - *p);
                let er = loo_resid(&xr, refs.len(), k, &th_s);
                let csr = sparse_fit(&er, &rr, refs.len(), k, nu, b2.1, b2.2)?;
                let et = loo_resid(&xt, nt, k, &th_s);
                let mut add = vec![0.0; nt * nu];
                sparse_apply(&mut add, &et, nt, k, nu, &csr);
                im_c.par_iter_mut().zip(&add).enumerate().for_each(|(ix, (p, a))| *p += a * ebs[ix % nu]);
                writeln!(wl, "  {pname}: stage 2 nnz {}", csr.val.len())?;
            }
            writeln!(wl, "  {pname}: ridge alpha {alpha} (inner r {:.4}); spectral K {} tau {} (inner r {:.4}); stage 2 K {} t {} (inner spa {:.4} vs {:.4}); {:.1}s",
                best.0, kfin, bs.2, bs.0, b2.1, b2.2, b2.0, base_spa, tpn.elapsed().as_secs_f64())?;
            // OOR: subspace (as eval_panel), leave-one-out fraction, rank-mean combination
            let mut css = vec![0.0; k * k];
            for (i, &si) in s.iter().enumerate() { for (j, &sj) in s.iter().enumerate() { css[i * k + j] = gr.c[si * g + sj]; } }
            let kq = 160.min(k / 2);
            let (pv, _) = la::top_eig(&css, k, kq)?;
            let proj = la::gemm(false, false, nt, kq, k, &xt, &pv);
            let oor: Vec<f64> = (0..nt).map(|i| { let nx: f64 = xt[i * k..(i + 1) * k].iter().map(|v| v * v).sum(); 1.0 - proj[i * kq..(i + 1) * kq].iter().map(|v| v * v).sum::<f64>() / nx }).collect();
            let et = loo_resid(&xt, nt, k, &th_r);
            let loo: Vec<f64> = (0..nt).map(|i| { let nx: f64 = xt[i * k..(i + 1) * k].iter().map(|v| v * v).sum(); et[i * k..(i + 1) * k].iter().map(|v| v * v).sum::<f64>() / nx }).collect();
            let comb = ranks_mean(&oor, &loo);
            for (mname, im) in [("ridge", &im_r), ("spectral", &im_s), ("spectral+sparse", &im_c)] {
                let rows: Vec<String> = (0..nt).into_par_iter().map(|i| {
                    let (yi, ti, di, dd) = (&im[i * nu..(i + 1) * nu], &tr[i * nu..(i + 1) * nu], &dr[i * nu..(i + 1) * nu], &den[i * nu..(i + 1) * nu]);
                    let rt = resid(ti, &basis);
                    let nde = ti.iter().filter(|x| x.abs() > DE_Z).count();
                    let mc = di.iter().filter(|x| x.abs() > DE_Z).count();
                    let matched = |v: &[f64]| -> (usize, usize) {
                        let mut o: Vec<usize> = (0..nu).collect();
                        o.sort_by(|&a, &b| v[b].abs().partial_cmp(&v[a].abs()).unwrap());
                        let (mut tp, mut disc) = (0, 0);
                        for &j in &o[..mc] { let same = v[j].signum() == ti[j].signum(); if !same { disc += 1; } if same && ti[j].abs() > DE_Z { tp += 1; } }
                        (tp, disc)
                    };
                    let (tp, disc) = matched(yi);
                    let (mtp, mdisc) = matched(di);
                    format!("{f}\t{pname}\t{mname}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{nde}\t{mc}\t{tp}\t{disc}\t{mtp}\t{mdisc}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                        zt.perts[tests[i]], pearson(yi, ti), pearson(&resid(yi, &basis), &rt), pearson(&resid(dd, &basis), &rt), pearson(di, ti), pearson(&resid(di, &basis), &rt),
                        spa_de(yi, ti), spa_de(di, ti), oor[i], loo[i], comb[i])
                }).collect();
                for row in rows { writeln!(w, "{row}")?; }
            }
        }
        writeln!(wl, "fold {f} done at {:.1}s", t0.elapsed().as_secs_f64())?;
        wl.flush()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn rnd(n: usize, seed: u64) -> Vec<f64> { let mut st = seed; (0..n).map(|_| (sm64(&mut st) >> 11) as f64 / (1u64 << 53) as f64 - 0.5).collect() }
    /// low-rank greedy == dense greedy on the materialised covariance
    #[test]
    fn lowrank_greedy_matches_dense() {
        let (g, r, k) = (120, 9, 25);
        let u = rnd(g * r, 7);
        let nz: Vec<f64> = rnd(g, 9).iter().map(|x| 0.2 + x.abs() * 3.0).collect();
        let w: Vec<f64> = rnd(g, 11).iter().map(|x| 0.5 + x.abs()).collect();
        let fast = greedy_lowrank(&u, g, r, k, &nz, &w);
        let mut rm = la::gemm(false, true, g, g, r, &u, &u);
        let mut used = vec![false; g];
        let mut dense = vec![];
        for _ in 0..k {
            let mut best = (f64::NEG_INFINITY, 0);
            for j in 0..g { if used[j] { continue; } let num: f64 = (0..g).map(|i| rm[i * g + j].powi(2)).sum(); let v = num / (rm[j * g + j] + nz[j]) / w[j]; if v > best.0 { best = (v, j); } }
            let j = best.1; used[j] = true; dense.push(j);
            let col: Vec<f64> = (0..g).map(|i| rm[i * g + j]).collect(); let d = col[j] + nz[j];
            for i in 0..g { for l in 0..g { rm[i * g + l] -= col[i] * col[l] / d; } }
        }
        assert_eq!(fast, dense);
    }
    /// closed-form leave-one-out residual == explicit regression of each gene on the others
    #[test]
    fn loo_matches_refit() {
        let (k, n) = (12, 40);
        let x = rnd(n * k, 3);
        let mut s = la::gemm(true, false, k, k, n, &x, &x);
        for i in 0..k { s[i * k + i] += 0.7; }
        let mut th = vec![0.0; k * k];
        for i in 0..k { th[i * k + i] = 1.0; }
        la::posv(s.clone(), k, &mut th, k).unwrap();
        let q = rnd(k, 5);
        let e = loo_resid(&q, 1, k, &th);
        for gsel in 0..k {
            let oth: Vec<usize> = (0..k).filter(|&j| j != gsel).collect();
            let m = oth.len();
            let a: Vec<f64> = oth.iter().flat_map(|&p| oth.iter().map(move |&r| (p, r))).map(|(p, r)| s[p * k + r]).collect();
            let mut b: Vec<f64> = oth.iter().map(|&p| s[p * k + gsel]).collect();
            la::posv(a, m, &mut b, 1).unwrap();
            let pred: f64 = oth.iter().zip(&b).map(|(&p, bb)| bb * q[p]).sum();
            assert!((e[gsel] - (q[gsel] - pred)).abs() < 1e-9, "gene {gsel}: {} vs {}", e[gsel], q[gsel] - pred);
        }
    }
    /// sparse CSR application == dense product with the same coefficients
    #[test]
    fn sparse_matches_dense() {
        let (n, k, nu) = (60, 15, 20);
        let e = rnd(n * k, 13);
        let r = rnd(n * nu, 17);
        let csr = sparse_fit(&e, &r, n, k, nu, 4, 0.05).unwrap();
        let mut dense = vec![0.0; nu * k];
        for j in 0..nu { for p in csr.ptr[j]..csr.ptr[j + 1] { dense[j * k + csr.idx[p]] = csr.val[p]; } }
        let mut pred = vec![0.0; n * nu];
        sparse_apply(&mut pred, &e, n, k, nu, &csr);
        let ref_ = la::gemm(false, true, n, nu, k, &e, &dense);
        for (a, b) in pred.iter().zip(&ref_) { assert!((a - b).abs() < 1e-12); }
        assert!(csr.val.len() <= 4 * nu);
    }
    /// spectral weights solve the normal equations Sigma_SS W' = Sigma_SU
    #[test]
    fn spectral_solves_normal_equations() {
        let (g, n, nk) = (40, 200, 10);
        let y = rnd(n * g, 21);
        let c = la::gram(&y, n, g);
        let sp = spec_from(&c, g, n, 0.01, nk).unwrap();
        let s: Vec<usize> = (0..12).collect();
        let u: Vec<usize> = (12..g).collect();
        let (wt, th) = spec_w(&sp, &s, &u, 8, 0.5).unwrap();
        let (css, cus) = spec_cov(&sp, &s, &u, 8, 0.5);
        let back = la::gemm(false, false, u.len(), s.len(), s.len(), &wt, &css);
        for (a, b) in back.iter().zip(&cus) { assert!((a - b).abs() < 1e-9); }
        let id = la::gemm(false, false, s.len(), s.len(), s.len(), &th, &css);
        for i in 0..s.len() { for j in 0..s.len() { assert!((id[i * s.len() + j] - if i == j { 1.0 } else { 0.0 }).abs() < 1e-9); } }
    }
}

/// Wilcoxon signed-rank test, normal approximation with tie and zero handling (Pratt: zeros dropped).
pub fn wilcoxon(d: &[f64]) -> f64 {
    let x: Vec<f64> = d.iter().copied().filter(|v| v.is_finite() && *v != 0.0).collect();
    let n = x.len() as f64;
    if n < 10.0 { return f64::NAN; }
    let r = ranks(&x.iter().map(|v| v.abs()).collect::<Vec<_>>());
    let wp: f64 = x.iter().zip(&r).filter(|(v, _)| **v > 0.0).map(|(_, rr)| rr).sum();
    let mut sorted: Vec<f64> = x.iter().map(|v| v.abs()).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut tie = 0.0;
    let mut i = 0;
    while i < sorted.len() { let mut j = i; while j + 1 < sorted.len() && sorted[j + 1] == sorted[i] { j += 1; } let t = (j - i + 1) as f64; tie += t * t * t - t; i = j + 1; }
    let mu = n * (n + 1.0) / 4.0;
    let sd = (n * (n + 1.0) * (2.0 * n + 1.0) / 24.0 - tie / 48.0).sqrt();
    2.0 * phi(-(((wp - mu).abs() - 0.5).max(0.0) / sd))
}

/// Summary of {pre}.upg.tsv: per panel x method, responding perturbations (r_ceil >= 0.2).
pub fn summary(a: &[String]) -> Result<()> {
    let pre = req(a, "--in")?;
    let txt = std::fs::read_to_string(format!("{pre}.upg.tsv"))?;
    let mut lines = txt.lines();
    let hdr: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| hdr.iter().position(|h| *h == n).unwrap();
    let rows: Vec<Vec<String>> = lines.map(|l| l.split('\t').map(String::from).collect()).collect();
    let f = |r: &Vec<String>, n: &str| -> f64 { r[col(n)].parse().unwrap_or(f64::NAN) };
    let mut combos: Vec<(String, String)> = vec![];
    for r in &rows { let c = (r[col("panel")].clone(), r[col("method")].clone()); if !combos.contains(&c) { combos.push(c); } }
    let key = |r: &Vec<String>| format!("{}:{}", r[col("fold")], r[col("pert")]);
    let base: HashMap<String, f64> = rows.iter().filter(|r| r[col("panel")] == "greedy" && r[col("method")] == "ridge").map(|r| (key(r), f(r, "rperp"))).collect();
    println!("panel       method           n    ratio   spec_r  DEsens  DEdisc  sizeRank | direct: sens disc sizeRank | vs greedy+ridge: diff wins p | OOR AUC: subspace loo comb (fails)");
    for (p, m) in &combos {
        let sel: Vec<&Vec<String>> = rows.iter().filter(|r| &r[col("panel")] == p && &r[col("method")] == m).collect();
        let resp: Vec<&&Vec<String>> = sel.iter().filter(|r| f(r, "r_ceil") >= 0.2).collect();
        let ratio = median(resp.iter().map(|r| f(r, "rperp") / f(r, "den")).collect());
        let spec = median(resp.iter().map(|r| f(r, "rperp")).collect());
        let s = |n: &str| sel.iter().map(|r| f(r, n)).sum::<f64>();
        let d: Vec<f64> = resp.iter().map(|r| f(r, "rperp") - base[&key(r)]).collect();
        let fails: Vec<bool> = resp.iter().map(|r| f(r, "rperp") / f(r, "rperp_ceil") < 0.5).collect();
        let au = |n: &str| auc(&resp.iter().map(|r| f(r, n)).collect::<Vec<_>>(), &fails).unwrap_or(f64::NAN);
        println!("{p:11} {m:16} {:4} {ratio:.4} {spec:.4}  {:.4}  {:.4}  {:.4}   | {:.4} {:.4} {:.4} | {:+.4} {:.2} {:.1e} | {:.3} {:.3} {:.3} ({})",
            resp.len(), s("tp") / s("n_de"), s("disc") / s("m_call"), median(sel.iter().map(|r| f(r, "spa")).collect()),
            s("m_tp") / s("n_de"), s("m_disc") / s("m_call"), median(sel.iter().map(|r| f(r, "spa_dir")).collect()),
            median(d.clone()), d.iter().filter(|x| **x > 0.0).count() as f64 / d.len() as f64, wilcoxon(&d), au("oor"), au("loo"), au("comb"), fails.iter().filter(|x| **x).count());
    }
    Ok(())
}
#[cfg(test)]
mod wtests {
    use super::*;
    /// exact small-sample check of the normal approximation against the permutation distribution
    #[test]
    fn wilcoxon_close_to_exact() {
        let d: Vec<f64> = (1..=20).map(|i| if i % 3 == 0 { -(i as f64) } else { i as f64 * 0.7 + 0.1 }).collect();
        let n = d.len();
        let r = ranks(&d.iter().map(|v| v.abs()).collect::<Vec<_>>());
        let w: f64 = d.iter().zip(&r).filter(|(v, _)| **v > 0.0).map(|(_, x)| x).sum();
        let mu = (n * (n + 1)) as f64 / 4.0;
        let mut extreme = 0u64;
        for mask in 0u32..(1 << n) { let s: f64 = (0..n).filter(|&i| mask >> i & 1 == 1).map(|i| r[i]).sum(); if (s - mu).abs() >= (w - mu).abs() - 1e-9 { extreme += 1; } }
        let exact = extreme as f64 / (1u64 << n) as f64;
        let approx = wilcoxon(&d);
        assert!((approx - exact).abs() < 0.01, "approx {approx} exact {exact}");
    }
}

/// Simulation summary by budget for a set of folds: median specific r and DE sensitivity per
/// panel x assay, and paired difference against targeted L1000.
pub fn simsum(a: &[String]) -> Result<()> {
    let pre = req(a, "--in")?;
    let folds: Vec<String> = req(a, "--folds-sel")?.split(',').map(String::from).collect();
    let txt = std::fs::read_to_string(format!("{pre}.sim.tsv"))?;
    let mut lines = txt.lines();
    let hdr: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| hdr.iter().position(|h| *h == n).unwrap();
    let rows: Vec<Vec<String>> = lines.map(|l| l.split('\t').map(String::from).collect()).filter(|r: &Vec<String>| folds.contains(&r[col("fold")])).collect();
    let f = |r: &Vec<String>, n: &str| -> f64 { r[col(n)].parse().unwrap_or(f64::NAN) };
    let mut qs: Vec<String> = vec![];
    for r in &rows { if !qs.contains(&r[col("q")]) { qs.push(r[col("q")].clone()); } }
    let lmname = rows.iter().map(|r| r[col("panel")].clone()).find(|p| p.starts_with("l1000")).unwrap();
    for q in &qs {
        let rq: Vec<&Vec<String>> = rows.iter().filter(|r| &r[col("q")] == q).collect();
        let dir: Vec<&&Vec<String>> = rq.iter().filter(|r| r[col("panel")] == lmname && r[col("assay")] == "extract").collect();
        println!("q {q}: whole transcriptome rp {:.4} sens {:.3} (n {})", median(dir.iter().map(|r| f(r, "rp_dir")).collect()),
            dir.iter().map(|r| f(r, "tp_dir")).sum::<f64>() / dir.iter().map(|r| f(r, "n_de")).sum::<f64>(), dir.len());
        let refm: HashMap<String, f64> = rq.iter().filter(|r| r[col("panel")] == lmname && r[col("assay")] == "targeted").map(|r| (format!("{}:{}", r[col("fold")], r[col("pert")]), f(r, "rp_imp"))).collect();
        let mut combos: Vec<(String, String)> = vec![];
        for r in &rq { let c = (r[col("panel")].clone(), r[col("assay")].clone()); if c.1 != "extract" && !combos.contains(&c) { combos.push(c); } }
        for (p, asy) in combos {
            let s: Vec<&&Vec<String>> = rq.iter().filter(|r| r[col("panel")] == p && r[col("assay")] == asy).collect();
            let d: Vec<f64> = s.iter().map(|r| f(r, "rp_imp") - refm[&format!("{}:{}", r[col("fold")], r[col("pert")])]).collect();
            println!("  {p:15} {asy:13} share {:.3} rate {:.3} rp {:.4} sens {:.3} | vs targeted L1000 {:+.4} wins {:.2} p {:.1e}",
                f(s[0], "share"), f(s[0], "rate"), median(s.iter().map(|r| f(r, "rp_imp")).collect()),
                s.iter().map(|r| f(r, "tp_imp")).sum::<f64>() / s.iter().map(|r| f(r, "n_de")).sum::<f64>(),
                median(d.clone()), d.iter().filter(|x| **x > 0.0).count() as f64 / d.len() as f64, wilcoxon(&d));
        }
    }
    Ok(())
}

/// Summary of a cv/transfer run ({pre}.perts.tsv, {pre}.genes.tsv): per panel, responding
/// perturbations (r_ceil >= 0.2): recovery vs denoised (K = 160), share >= 0.8, INFO slope, OOR AUC.
pub fn cvsum(a: &[String]) -> Result<()> {
    let pre = req(a, "--in")?;
    let txt = std::fs::read_to_string(format!("{pre}.perts.tsv"))?;
    let mut lines = txt.lines();
    let hdr: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| hdr.iter().position(|h| *h == n).unwrap();
    let rows: Vec<Vec<String>> = lines.map(|l| l.split('\t').map(String::from).collect()).collect();
    let f = |r: &Vec<String>, n: &str| -> f64 { r[col(n)].parse().unwrap_or(f64::NAN) };
    let mut panels: Vec<String> = vec![];
    for r in &rows { if !panels.contains(&r[col("panel")]) { panels.push(r[col("panel")].clone()); } }
    let gtxt = std::fs::read_to_string(format!("{pre}.genes.tsv"))?;
    println!("panel            n_resp  recovery  share>=0.8  spec_r  raw_half  r_imp  r_direct  INFO_slope  OOR_AUC (fails)");
    for p in &panels {
        let s: Vec<&Vec<String>> = rows.iter().filter(|r| &r[col("panel")] == p && f(r, "r_ceil") >= 0.2).collect();
        let ratio: Vec<f64> = s.iter().map(|r| f(r, "r_imp_perp") / f(r, "r_den160_perp")).collect();
        let fails: Vec<bool> = s.iter().map(|r| f(r, "r_imp_perp") / f(r, "r_ceil_perp") < 0.5).collect();
        let au = auc(&s.iter().map(|r| f(r, "oor")).collect::<Vec<_>>(), &fails).unwrap_or(f64::NAN);
        let (mut sx, mut sy, mut sxx, mut sxy, mut n) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for l in gtxt.lines().skip(1) {
            let g: Vec<&str> = l.split('\t').collect();
            if g[1] != p { continue; }
            let (x, y): (f64, f64) = (g[3].parse().unwrap_or(f64::NAN), g[4].parse().unwrap_or(f64::NAN));
            if !x.is_finite() || !y.is_finite() { continue; }
            sx += x; sy += y; sxx += x * x; sxy += x * y; n += 1.0;
        }
        let slope = (sxy - sx * sy / n) / (sxx - sx * sx / n);
        println!("{p:16} {:6}  {:.3}     {:.3}       {:.3}   {:.3}     {:.3}  {:.3}     {:.3}       {:.3} ({})", s.len(), median(ratio.clone()),
            ratio.iter().filter(|x| **x >= 0.8).count() as f64 / ratio.len() as f64, median(s.iter().map(|r| f(r, "r_imp_perp")).collect()),
            median(s.iter().map(|r| f(r, "r_ceil_perp")).collect()), median(s.iter().map(|r| f(r, "r_imp")).collect()), median(s.iter().map(|r| f(r, "r_ceil")).collect()),
            slope, au, fails.iter().filter(|x| **x).count());
    }
    Ok(())
}
