//! Count-level simulation at equal sequencing budget (step 6).
//! Binomial thinning of single-cell UMIs from a coexsc cache, pseudobulk CP10k per half, z as in
//! coexsc z (s2 from half differences across the scenario's perturbations). Assays: whole
//! transcriptome at budget q (direct), panel genes extracted from it, and a targeted panel with the
//! same UMI budget concentrated on panel genes (rate min(1, q / panel UMI share), library = panel UMIs).
use crate::cache::Cache;
use crate::*;

const MED_CHI2_1: f64 = 0.454936423119572;

/// sums of CP10k over cells for each (group, half): returns ng x 2 x ncol
#[allow(clippy::too_many_arguments)]
pub(crate) fn pseudobulk(c: &Cache, groups: &[usize], colmap: &[i32], ncol: usize, rate: &[f64], wt: bool, q: f64, seed: u64) -> Vec<f64> {
    // work items: (group slot, half, cells)
    let mut items: Vec<(usize, usize, Vec<u32>)> = vec![];
    for (s, &gid) in groups.iter().enumerate() {
        let cells = c.cells(gid);
        for h in 0..2 {
            let hc: Vec<u32> = cells.iter().enumerate().filter(|(i, _)| i % 2 == h).map(|(_, &x)| x).collect();
            for ch in hc.chunks(512) { items.push((s, h, ch.to_vec())); }
        }
    }
    let exact = wt && q >= 1.0;
    let parts: Vec<(usize, usize, Vec<f64>)> = items.into_par_iter().map(|(s, h, cells)| {
        let mut acc = vec![0.0f64; ncol];
        let mut buf: Vec<(usize, f64)> = vec![];
        for &cell in &cells {
            let (idx, cnt) = c.row(cell as usize);
            buf.clear();
            let scale: f64;
            if exact {
                scale = c.scale[cell as usize] as f64;
                for (&i, &n) in idx.iter().zip(cnt) { buf.push((i as usize, n as f64)); }
            } else {
                let mut st = seed ^ (cell as u64 + 1).wrapping_mul(0xD1B54A32D192ED03);
                let (mut lib, mut raw_cached) = (0.0, 0.0);
                for (&i, &n) in idx.iter().zip(cnt) {
                    let r = rate[i as usize];
                    raw_cached += n as f64;
                    if r <= 0.0 { continue; }
                    let t = if r >= 1.0 { n as u32 } else { (0..n).filter(|_| ((sm64(&mut st) >> 11) as f64) < r * (1u64 << 53) as f64).count() as u32 };
                    lib += t as f64;
                    if t > 0 { buf.push((i as usize, t as f64)); }
                }
                if wt {
                    // uncached genes contribute their expected thinned count to the library
                    let full = if c.scale[cell as usize] > 0.0 { 1e4 / c.scale[cell as usize] as f64 } else { 0.0 };
                    lib += q * (full - raw_cached).max(0.0);
                }
                scale = if lib > 0.0 { 1e4 / lib } else { 0.0 };
            }
            for &(i, t) in &buf { let k = colmap[i]; if k >= 0 { acc[k as usize] += t * scale; } }
        }
        (s, h, acc)
    }).collect();
    let mut out = vec![0.0; groups.len() * 2 * ncol];
    for (s, h, acc) in parts {
        let o = &mut out[(s * 2 + h) * ncol..(s * 2 + h + 1) * ncol];
        for (a, b) in o.iter_mut().zip(&acc) { *a += b; }
    }
    out
}

/// z for half 0 and half 1 of each query group (query = groups[1..], groups[0] = control)
pub(crate) fn zhalves(c: &Cache, groups: &[usize], sums: &[f64], ncol: usize) -> (Vec<f64>, Vec<f64>) {
    let n = |s: usize| c.groups[groups[s]].2;
    let nq = groups.len() - 1;
    let nc = n(0) as f64;
    let mc: Vec<f64> = (0..ncol).map(|j| (sums[j] + sums[ncol + j]) / nc).collect();
    let mean = |s: usize, h: usize, j: usize| { let (n0, n1) = ((n(s) + 1) / 2, n(s) / 2); sums[(s * 2 + h) * ncol + j] / [n0, n1][h] as f64 };
    let s2: Vec<f64> = (0..ncol).into_par_iter().map(|j| {
        let mut v: Vec<f64> = (1..=nq).filter(|&s| n(s) >= 2).map(|s| {
            let (n0, n1) = (((n(s) + 1) / 2) as f64, (n(s) / 2) as f64);
            (mean(s, 0, j) - mean(s, 1, j)).powi(2) / (1.0 / n0 + 1.0 / n1)
        }).collect();
        median(std::mem::take(&mut v)) / MED_CHI2_1
    }).collect();
    let mut z0 = vec![0.0; nq * ncol];
    let mut z1 = vec![0.0; nq * ncol];
    for s in 1..=nq {
        let (n0, n1) = (((n(s) + 1) / 2) as f64, (n(s) / 2) as f64);
        for j in 0..ncol {
            if s2[j] > 0.0 {
                z0[(s - 1) * ncol + j] = (mean(s, 0, j) - mc[j]) / (s2[j] * (1.0 / n0 + 1.0 / nc)).sqrt();
                z1[(s - 1) * ncol + j] = (mean(s, 1, j) - mc[j]) / (s2[j] * (1.0 / n1 + 1.0 / nc)).sqrt();
            }
        }
    }
    (z0, z1)
}

/// ridge with an extra diagonal on the panel block (heteroscedastic panel noise, Gram units)
fn ridge_diag(c: &[f64], g: usize, s: &[usize], u: &[usize], alpha: f64, extra: &[f64]) -> Result<Vec<f64>> {
    let k = s.len();
    let mut a = vec![0.0; k * k];
    for (i, &si) in s.iter().enumerate() { for (j, &sj) in s.iter().enumerate() { a[i * k + j] = c[si * g + sj]; } }
    let tr: f64 = (0..k).map(|i| a[i * k + i]).sum::<f64>() / k as f64;
    for i in 0..k { a[i * k + i] += alpha * tr + extra[i]; }
    let mut bt = vec![0.0; u.len() * k];
    bt.par_chunks_mut(k).zip(u).for_each(|(r, &ui)| { for (j, &sj) in s.iter().enumerate() { r[j] = c[ui * g + sj]; } });
    la::posv(a, k, &mut bt, u.len())?;
    Ok(bt)
}

/// alpha on inner holdout (as eval_panel) and final ridge weights (|u| x k)
pub(crate) fn fit(z: &Z, gr: &Grams, inner_ho: &[usize], p: &Panel, u: &[usize]) -> Result<(f64, Vec<f64>)> {
    let g = z.genes.len();
    let nu = u.len();
    let truth = z.mat(2, inner_ho, u);
    let mut best = (f64::NEG_INFINITY, ALPHAS[0]);
    for &al in &ALPHAS {
        let wt = ridge(&gr.c_in, g, &p.genes, u, al)?;
        let im = impute(z, inner_ho, &p.genes, &wt, nu, 1);
        let r: f64 = (0..inner_ho.len()).map(|i| pearson(&im[i * nu..(i + 1) * nu], &truth[i * nu..(i + 1) * nu])).filter(|x| x.is_finite()).sum::<f64>() / inner_ho.len() as f64;
        if r > best.0 { best = (r, al); }
    }
    Ok((best.1, ridge(&gr.c, g, &p.genes, u, best.1)?))
}

/// greedy forward MMSE selection: R = signal covariance, noise variance nz_j per gene
pub(crate) fn greedy(gr: &Grams, nref: usize, g: usize, k: usize, nz: &[f64]) -> Vec<usize> { greedy_cap(gr, nref, g, k, nz, &vec![0.0; g], f64::INFINITY, None).0 }
/// share-capped greedy: skip genes whose share would push the panel total above cap; returns (genes, objective)
fn greedy_cap(gr: &Grams, nref: usize, g: usize, k: usize, nz: &[f64], sh: &[f64], cap: f64, w: Option<&[f64]>) -> (Vec<usize>, f64) {
    let mut r: Vec<f64> = gr.c.iter().map(|v| v / nref as f64).collect();
    for j in 0..g { r[j * g + j] = (r[j * g + j] - 1.0).max(0.0); }
    let mut norm: Vec<f64> = r.par_chunks(g).map(|row| row.iter().map(|v| v * v).sum()).collect();
    let mut sel: Vec<usize> = vec![];
    let mut used = vec![false; g];
    let (mut tot, mut obj) = (0.0, 0.0);
    // candidates sorted by share, for the reserve of the remaining slots
    let mut bysh: Vec<usize> = (0..g).filter(|&j| nz[j].is_finite()).collect();
    bysh.sort_by(|&a, &b| sh[a].partial_cmp(&sh[b]).unwrap().then(a.cmp(&b)));
    for t in 0..k {
        let need = k - t - 1;
        let reserve: f64 = if cap.is_finite() { bysh.iter().filter(|&&j| !used[j]).take(need + 1).map(|&j| sh[j]).collect::<Vec<_>>().iter().sum() } else { 0.0 };
        let small: Vec<usize> = if cap.is_finite() { bysh.iter().copied().filter(|&j| !used[j]).take(need + 1).collect() } else { vec![] };
        let last = small.last().map(|&j| sh[j]).unwrap_or(0.0);
        let mut in_small = vec![false; g];
        for &j in &small { in_small[j] = true; }
        let best = (0..g).into_par_iter().filter(|&j| !used[j] && nz[j].is_finite() && {
                // sum of the `need` smallest unused shares excluding j
                let res = if !cap.is_finite() { 0.0 } else if in_small[j] { reserve - sh[j] } else { reserve - last };
                tot + sh[j] + res <= cap * (1.0 + 1e-12)
            })
            .map(|j| (norm[j] / (r[j * g + j] + nz[j]) / w.map_or(1.0, |w| w[j]), j))
            .reduce(|| (f64::NEG_INFINITY, usize::MAX), |a, b| if b.0 > a.0 || (b.0 == a.0 && b.1 < a.1) { b } else { a });
        let j = best.1;
        if j == usize::MAX { break; }
        used[j] = true;
        sel.push(j);
        tot += sh[j];
        obj += best.0;
        let col: Vec<f64> = r[j * g..(j + 1) * g].to_vec();
        let d = col[j] + nz[j];
        r.par_chunks_mut(g).zip(norm.par_iter_mut()).enumerate().for_each(|(i, (row, nm))| {
            let f = col[i] / d;
            let mut s = 0.0;
            for (v, c) in row.iter_mut().zip(&col) { *v -= f * c; s += *v * *v; }
            *nm = s;
        });
    }
    (sel, obj)
}

/// per-row metrics: r, specific r, matched DE (top n_de by |v|) true positives and sign discordance
pub(crate) fn metrics(v: &[f64], t: &[f64], basis: &[Vec<f64>]) -> (f64, f64, usize, usize, usize) {
    let nde = t.iter().filter(|x| x.abs() > DE_Z).count();
    let mut o: Vec<usize> = (0..v.len()).collect();
    o.sort_by(|&a, &b| v[b].abs().partial_cmp(&v[a].abs()).unwrap());
    let (mut tp, mut disc) = (0, 0);
    for &j in &o[..nde] { let same = v[j].signum() == t[j].signum(); if !same { disc += 1; } if same && t[j].abs() > DE_Z { tp += 1; } }
    (pearson(v, t), pearson(&resid(v, basis), &resid(t, basis)), nde, tp, disc)
}

#[allow(clippy::too_many_arguments)]
pub fn run(a: &[String], z: &Z, fold: &[i64], group: &[i64], excl: &HashSet<(usize, i64)>, l1000: &HashSet<String>, sizes: &[usize], out: &str) -> Result<()> {
    let t0 = std::time::Instant::now();
    let c = Cache::open(&req(a, "--cache")?)?;
    let qs: Vec<f64> = arg(a, "--budgets").map(|s| s.split(',').map(|x| x.parse().unwrap()).collect()).unwrap_or(vec![0.02, 0.05, 0.1, 0.25, 1.0]);
    let runs: Vec<i64> = arg(a, "--folds-run").map(|s| s.split(',').map(|x| x.parse().unwrap()).collect()).unwrap_or((0..5).collect());
    let g = z.genes.len();
    let cgi: HashMap<&str, usize> = c.genes.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let gmap: Vec<usize> = z.genes.iter().map(|s| cgi.get(s.as_str()).copied().ok_or_else(|| anyhow!("gene {s} not in cache"))).collect::<Result<_>>()?;
    let ncg = c.genes.len();
    let mut colmap_all = vec![-1i32; ncg];
    for (j, &cg) in gmap.iter().enumerate() { colmap_all[cg] = j as i32; }
    let pi: HashMap<&str, usize> = z.perts.iter().enumerate().map(|(i, p)| (p.as_str(), i)).collect();
    let ctrl = c.group("control").ok_or_else(|| anyhow!("no control"))?;
    // truth: full depth, whole transcriptome, all perturbations (s2 over all), half 1
    let allq: Vec<usize> = std::iter::once(ctrl).chain((0..c.groups.len()).filter(|&i| i != ctrl && c.groups[i].2 >= 2)).collect();
    let sums = pseudobulk(&c, &allq, &colmap_all, g, &vec![1.0; ncg], true, 1.0, 0);
    let (tz0, tz1) = zhalves(&c, &allq, &sums, g);
    {
        // dump truth z (layers half 0, half 1) for validation against coexsc z
        let mut w = BufWriter::new(std::fs::File::create(format!("{out}.truth.tsv"))?);
        writeln!(w, "pert\tgene\tz0\tz1")?;
        for (s, &gid) in allq.iter().enumerate().skip(1).step_by(97) {
            for j in (0..g).step_by(53) { writeln!(w, "{}\t{}\t{:.6}\t{:.6}", c.groups[gid].0, z.genes[j], tz0[(s - 1) * g + j], tz1[(s - 1) * g + j])?; }
        }
    }
    let trow: HashMap<usize, usize> = allq.iter().enumerate().skip(1).map(|(s, &gid)| (gid, s - 1)).collect();
    eprintln!("truth z: {} perturbations at {:.1}s", allq.len() - 1, t0.elapsed().as_secs_f64());
    // control UMI share of each cached gene (for targeted rates)
    let (mut gs, mut libsum) = (vec![0.0f64; ncg], 0.0);
    let (mut m1, mut mt, mut m2) = (vec![0.0f64; ncg], vec![0.0f64; ncg], vec![0.0f64; ncg]);
    for &cell in c.cells(ctrl) {
        let (idx, cnt) = c.row(cell as usize);
        let sc = c.scale[cell as usize] as f64;
        for (&i, &n) in idx.iter().zip(cnt) {
            let (i, n) = (i as usize, n as f64);
            gs[i] += n; m1[i] += n * sc; mt[i] += n * sc * sc; m2[i] += (n * sc).powi(2);
        }
        libsum += 1e4 / c.scale[cell as usize] as f64;
    }
    let ncc = c.cells(ctrl).len() as f64;
    // per universe gene: (cv2_bio, cv2_tech, cv2_tot)
    let cv: Vec<(f64, f64, f64)> = gmap.iter().map(|&cg| {
        let (a, t, b) = (m1[cg] / ncc, mt[cg] / ncc, m2[cg] / ncc);
        if a <= 0.0 { return (f64::NAN, f64::NAN, f64::NAN); }
        let tot = (b - a * a) / (a * a);
        let tech = t / (a * a);
        ((tot - tech).max(0.0), tech, tot)
    }).collect();
    let nzf = |r: f64| -> Vec<f64> { cv.iter().map(|&(bio, tech, tot)| if tot.is_finite() && tot > 0.0 { (bio + tech / r) / tot } else { f64::INFINITY }).collect() };
    let upg2 = a.contains(&"--upg2".to_string());
    let budgets: Vec<f64> = arg(a, "--design-budgets").map(|s| s.split(',').map(|x| x.parse().unwrap()).collect()).unwrap_or(vec![]);
    {
        let mut w = BufWriter::new(std::fs::File::create(format!("{out}.cv2.tsv"))?);
        writeln!(w, "gene\tcv2_bio\tcv2_tech\tcv2_tot\tshare")?;
        for (j, &(b, t, o)) in cv.iter().enumerate() { writeln!(w, "{}\t{b:.6e}\t{t:.6e}\t{o:.6e}\t{:.6e}", z.genes[j], gs[gmap[j]] / libsum)?; }
    }
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.sim.tsv"))?);
    let mut wden = BufWriter::new(std::fs::File::create(format!("{out}.den.tsv"))?);
    writeln!(wden, "fold\tq\tpert\tn_de\tr_den\trp_den\ttp_den\tdisc_den\trp_dir\ttp_dir\tdisc_dir")?;
    writeln!(w, "fold\tpanel\tassay\tq\tpert\tn_cells\tshare\trate\talpha\tn_de\tr_imp\trp_imp\ttp_imp\tdisc_imp\tr_dir\trp_dir\ttp_dir\tdisc_dir")?;
    let inner_split = |refs: &[usize]| -> (Vec<usize>, Vec<usize>) {
        let key = |i: usize| -> u64 { let k = if group[i] >= 0 { group[i] as u64 } else { 1_000_000 + i as u64 }; (k.wrapping_mul(0x9E3779B97F4A7C15) >> 32) % 10 };
        (refs.iter().copied().filter(|&i| key(i) != 0).collect(), refs.iter().copied().filter(|&i| key(i) == 0).collect())
    };
    for &f in &runs {
        let refs: Vec<usize> = (0..z.perts.len()).filter(|&i| fold[i] != f && !excl.contains(&(i, f))).collect();
        let (itr, iho) = inner_split(&refs);
        let gr = grams(z, &refs, &itr)?;
        let mut ps = panels(z, &refs, sizes, l1000, 1)?;
        if upg2 {
            let shg: Vec<f64> = (0..g).map(|j| gs[gmap[j]] / libsum).collect();
            let tg = std::time::Instant::now();
            ps.push(Panel { name: "greedyfull_500".into(), genes: greedy(&gr, refs.len(), g, 500, &vec![1.0; g]) });
            for &al in &[0.5f64, 1.0] {
                let wv: Vec<f64> = shg.iter().map(|s| (s + 1e-4).powf(al)).collect();
                ps.push(Panel { name: format!("cost{al}_500"), genes: greedy_cap(&gr, refs.len(), g, 500, &vec![1.0; g], &shg, f64::INFINITY, Some(&wv)).0 });
            }
            ps.push(Panel { name: "cap0.12_500".into(), genes: greedy_cap(&gr, refs.len(), g, 500, &vec![1.0; g], &shg, 0.12, None).0 });
            for p in &ps { eprintln!("fold {f}: panel {} share {:.3}", p.name, p.genes.iter().map(|&j| shg[j]).sum::<f64>()); }
            eprintln!("fold {f}: upg2 panels in {:.1}s", tg.elapsed().as_secs_f64());
        }
        if !budgets.is_empty() {
            let share = |gs_: &[usize]| gs_.iter().map(|&j| gs[gmap[j]]).sum::<f64>() / libsum;
            let k = 500;
            let tg = std::time::Instant::now();
            ps.push(Panel { name: format!("greedyfull_{k}"), genes: greedy(&gr, refs.len(), g, k, &nzf(1.0)) });
            let shg: Vec<f64> = (0..g).map(|j| gs[gmap[j]] / libsum).collect();
            for &q in &budgets {
                let mut best = (f64::NEG_INFINITY, vec![], 0.0);
                for &cap in &[0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.2, 0.3, 0.5] {
                    let r = (q / cap).min(1.0);
                    let (sel, obj) = greedy_cap(&gr, refs.len(), g, k, &nzf(r), &shg, cap, None);
                    eprintln!("fold {f} budget {q} cap {cap}: rate {r:.3} genes {} share {:.3} objective {obj:.2}", sel.len(), share(&sel));
                    if sel.len() == k && obj > best.0 { best = (obj, sel, cap); }
                }
                eprintln!("fold {f} budget {q}: chosen cap {}", best.2);
                ps.push(Panel { name: format!("budget{q}_{k}"), genes: best.1 });
            }
            eprintln!("fold {f}: greedy panels in {:.1}s", tg.elapsed().as_secs_f64());
        }
        let u = common_set(&ps, g, true).unwrap();
        let nu = u.len();
        let basis = basis_on(&gr.v, NPC, &u);
        // query: essential-screen perturbations whose target is in test fold f
        let qg: Vec<usize> = (0..c.groups.len()).filter(|&i| i != ctrl && c.groups[i].2 >= 4 && pi.get(c.groups[i].0.as_str()).map(|&p| fold[p] == f).unwrap_or(false)).collect();
        let grp: Vec<usize> = std::iter::once(ctrl).chain(qg.iter().copied()).collect();
        let nq = qg.len();
        let truth: Vec<f64> = (0..nq).flat_map(|s| { let r = trow[&qg[s]]; u.iter().map(move |&j| (r, j)) }).map(|(r, j)| tz1[r * g + j]).collect();
        let fits: Vec<(f64, Vec<f64>)> = ps.iter().map(|p| fit(z, &gr, &iho, p, &u)).collect::<Result<_>>()?;
        eprintln!("fold {f}: refs {} query {nq} common genes {nu} at {:.1}s", refs.len(), t0.elapsed().as_secs_f64());
        for (qi, &q) in qs.iter().enumerate() {
            let rate_wt = vec![q.min(1.0); ncg];
            let sw = pseudobulk(&c, &grp, &colmap_all, g, &rate_wt, true, q, (f as u64 + 11) * 1000 + qi as u64);
            let (zw, _) = zhalves(&c, &grp, &sw, g);
            let dir: Vec<f64> = (0..nq).flat_map(|s| u.iter().map(move |&j| (s, j))).map(|(s, j)| zw[s * g + j]).collect();
            let dm: Vec<(f64, f64, usize, usize, usize)> = (0..nq).into_par_iter().map(|s| metrics(&dir[s * nu..(s + 1) * nu], &truth[s * nu..(s + 1) * nu], &basis)).collect();
            // control: the thinned whole transcriptome projected on the 160 leading reference programs
            {
                let nk = *DEN_K.iter().max().unwrap();
                let v160: Vec<f64> = (0..g).flat_map(|j| gr.vk[j * nk..j * nk + 160].to_vec()).collect();
                let full: Vec<f64> = zw[..nq * g].to_vec();
                let den = la::gemm(false, true, nq, g, 160, &la::gemm(false, false, nq, 160, g, &full, &v160), &v160);
                for s in 0..nq {
                    let d: Vec<f64> = u.iter().map(|&j| den[s * g + j]).collect();
                    let m = metrics(&d, &truth[s * nu..(s + 1) * nu], &basis);
                    writeln!(wden, "{f}\t{q}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t{}\t{}", c.groups[qg[s]].0, m.2, m.0, m.1, m.3, m.4, dm[s].1, dm[s].3, dm[s].4)?;
                }
            }
            for (pk, p) in ps.iter().enumerate() {
                let k = p.genes.len();
                let share: f64 = p.genes.iter().map(|&j| gs[gmap[j]]).sum::<f64>() / libsum;
                let rate = (q / share).min(1.0);
                // targeted assay: only panel genes sequenced
                let mut rt = vec![0.0; ncg];
                let mut cm = vec![-1i32; ncg];
                for (kk, &j) in p.genes.iter().enumerate() { rt[gmap[j]] = rate; cm[gmap[j]] = kk as i32; }
                let st = pseudobulk(&c, &grp, &cm, k, &rt, false, q, (f as u64 + 11) * 1000 + 100 * (pk as u64 + 1) + qi as u64);
                let (zt0, _) = zhalves(&c, &grp, &st, k);
                let xe: Vec<f64> = (0..nq).flat_map(|s| p.genes.iter().map(move |&j| (s, j))).map(|(s, j)| zw[s * g + j]).collect();
                // heteroscedastic ridge: per-gene thinning noise, input rescaled by 1/sqrt(f_j)
                let inv_f: Vec<f64> = nzf(rate);
                let extra: Vec<f64> = p.genes.iter().map(|&j| refs.len() as f64 * (inv_f[j] - 1.0).max(0.0)).collect();
                let wt_h = ridge_diag(&gr.c, g, &p.genes, &u, fits[pk].0, &extra)?;
                let xh: Vec<f64> = (0..nq).flat_map(|s| p.genes.iter().enumerate().map(move |(kk, &j)| (s, kk, j))).map(|(s, kk, j)| zt0[s * k + kk] * inv_f[j].sqrt()).collect();
                for (assay, x, wsel) in [("extract", &xe, &fits[pk].1), ("targeted", &zt0, &fits[pk].1), ("targeted_het", &xh, &wt_h)] {
                    if assay == "targeted_het" && !upg2 { continue; }
                    let im = la::gemm(false, true, nq, nu, k, x, wsel);
                    let rows: Vec<String> = (0..nq).into_par_iter().map(|s| {
                        let m = metrics(&im[s * nu..(s + 1) * nu], &truth[s * nu..(s + 1) * nu], &basis);
                        let d = dm[s];
                        format!("{f}\t{}\t{assay}\t{q}\t{}\t{}\t{share:.5}\t{rate:.5}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}", p.name, c.groups[qg[s]].0, c.groups[qg[s]].2, fits[pk].0, m.2, m.0, m.1, m.3, m.4, d.0, d.1, d.3, d.4)
                    }).collect();
                    for r in rows { writeln!(w, "{r}")?; }
                }
            }
            eprintln!("fold {f} q {q} done at {:.1}s", t0.elapsed().as_secs_f64());
        }
    }
    Ok(())
}

/// validation helper: pseudobulk sums for named groups at budget q, whole transcriptome (wt) or
/// targeted on a gene list at a given rate; writes group, half, gene, sum
pub fn pbdump(a: &[String]) -> Result<()> {
    let c = Cache::open(&req(a, "--cache")?)?;
    let names: Vec<String> = req(a, "--groups")?.split(',').map(String::from).collect();
    let q: f64 = req(a, "--q")?.parse()?;
    let seed: u64 = req(a, "--seed")?.parse()?;
    let grp: Vec<usize> = names.iter().map(|n| c.group(n).ok_or_else(|| anyhow!("group {n}"))).collect::<Result<_>>()?;
    let ncg = c.genes.len();
    let (colmap, rate, wt, ncol): (Vec<i32>, Vec<f64>, bool, usize) = match arg(a, "--panel") {
        None => ((0..ncg as i32).collect(), vec![q.min(1.0); ncg], true, ncg),
        Some(pf) => {
            let gi: HashMap<&str, usize> = c.genes.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
            let r: f64 = req(a, "--rate")?.parse()?;
            let (mut cm, mut rt) = (vec![-1i32; ncg], vec![0.0; ncg]);
            let gl: Vec<usize> = std::fs::read_to_string(pf)?.lines().filter_map(|s| gi.get(s).copied()).collect();
            for (k, &j) in gl.iter().enumerate() { cm[j] = k as i32; rt[j] = r; }
            (cm, rt, false, gl.len())
        }
    };
    let s = pseudobulk(&c, &grp, &colmap, ncol, &rate, wt, q, seed);
    let mut w = BufWriter::new(std::fs::File::create(req(a, "--out")?)?);
    writeln!(w, "group\thalf\tcol\tsum")?;
    for (gi, n) in names.iter().enumerate() {
        for h in 0..2 { for j in 0..ncol { let v = s[(gi * 2 + h) * ncol + j]; if v != 0.0 { writeln!(w, "{n}\t{h}\t{j}\t{v:.9e}")?; } } }
    }
    Ok(())
}

/// validation helper: greedy selection on all perturbations of --z with noise from --noise (one per gene)
pub fn greedytest(a: &[String], z: &Z) -> Result<()> {
    let g = z.genes.len();
    let refs: Vec<usize> = (0..z.perts.len()).collect();
    let gr = grams(z, &refs, &refs)?;
    let nz: Vec<f64> = std::fs::read_to_string(req(a, "--noise")?)?.lines().map(|l| l.parse().unwrap()).collect();
    let k: usize = req(a, "--k")?.parse()?;
    let sh: Vec<f64> = match arg(a, "--share") { Some(f) => std::fs::read_to_string(f)?.lines().map(|l| l.parse().unwrap()).collect(), None => vec![0.0; g] };
    let cap: f64 = arg(a, "--cap").map(|v| v.parse().unwrap()).unwrap_or(f64::INFINITY);
    let (sel, obj) = greedy_cap(&gr, refs.len(), g, k, &nz, &sh, cap, None);
    eprintln!("objective {obj:.10e}");
    for j in sel { println!("{j}"); }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn rnd(n: usize, seed: u64) -> Vec<f64> { let mut st = seed; (0..n).map(|_| (sm64(&mut st) >> 11) as f64 / (1u64 << 53) as f64 - 0.5).collect() }
    fn small() -> (Grams, usize, usize) {
        let (n, g) = (80, 30);
        let y: Vec<f64> = rnd(n * g, 5).iter().map(|v| v * 6.0).collect();
        let c = la::gram(&y, n, g);
        (Grams { c: c.clone(), c_in: c, v: vec![], vk: vec![], ev: vec![] }, n, g)
    }
    #[test]
    fn ridge_diag_zero_equals_ridge() {
        let (gr, _, g) = small();
        let s: Vec<usize> = (0..8).collect();
        let u: Vec<usize> = (8..g).collect();
        let a = ridge(&gr.c, g, &s, &u, 0.1).unwrap();
        let b = ridge_diag(&gr.c, g, &s, &u, 0.1, &vec![0.0; s.len()]).unwrap();
        for (x, y) in a.iter().zip(&b) { assert!((x - y).abs() < 1e-12); }
    }
    #[test]
    fn unit_weights_equal_plain_greedy() {
        let (gr, n, g) = small();
        let plain = greedy(&gr, n, g, 10, &vec![1.0; g]);
        let w1 = greedy_cap(&gr, n, g, 10, &vec![1.0; g], &vec![0.0; g], f64::INFINITY, Some(&vec![1.0; g])).0;
        assert_eq!(plain, w1);
        let wv: Vec<f64> = (0..g).map(|j| 1.0 + j as f64).collect();
        let ww = greedy_cap(&gr, n, g, 10, &vec![1.0; g], &vec![0.0; g], f64::INFINITY, Some(&wv)).0;
        assert_ne!(plain, ww);
    }
}
