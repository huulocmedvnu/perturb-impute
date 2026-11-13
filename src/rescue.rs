//! Perturbations with few cells in an unevenly sampled screen (Additional file 1: Table S10).
//!   panelimp thin   --cache DIR --z PREFIX --out PREFIX [--cells 15,30,50]
//!   panelimp rescue --ref PREFIX --ref-folds F --ref-hubs H [--screen PREFIX --folds F --hubs H [--same-line]]
//!                   --corum FILE --genesets FILE --out PREFIX [--tss FILE --umi FILE]
//!                   [--thin PREFIX --thin-cells 15,30,50] [--min-test-in N] [--min-ref-cells N]
//! `thin` writes the z-scores of a cell cache on the genes of --z (all cells, input half, held-out half) and the
//! input half restricted to its first M cells. `rescue` tests every perturbation of the screen once, in the fold
//! that holds it, and compares readouts of its input half with the held-out half:
//!   direct measurement, the hybrid profile (500 panel genes measured, other genes imputed), the whole transcriptome
//!   projected on the reference programs, and the measured profile shrunk towards its projection.
//! External reference: the atlas (--ref), without the test fold when the screen is the atlas itself or shares its
//! cell line (--same-line). Internal reference: the perturbations of the screen itself in the other folds with at
//! least --min-ref-cells cells (default: the median of the screen).
use crate::bio::*;
use crate::cache::Cache;
use memmap2::Mmap;
use crate::causal::*;
use crate::explore::read_tss;
use crate::*;

const K_INTERNAL: [usize; 4] = [20, 50, 100, 160];
const STRATA: [(&str, usize, usize); 3] = [("<=50", 0, 50), ("51-150", 51, 150), (">150", 151, usize::MAX)];

pub(crate) struct Model { pub s: Vec<usize>, pub u: Vec<usize>, pub w: Vec<f32>, pub v: Vec<f32>, pub kp: usize, pub alpha: f64 }

fn inner(refs: &[usize], group: &[i64]) -> (Vec<usize>, Vec<usize>) {
    let key = |i: usize| -> u64 { let k = if group[i] >= 0 { group[i] as u64 } else { 1_000_000 + i as u64 }; (k.wrapping_mul(0x9E3779B97F4A7C15) >> 32) % 10 };
    (refs.iter().copied().filter(|&i| key(i) != 0).collect(), refs.iter().copied().filter(|&i| key(i) == 0).collect())
}
fn mean_r(p: &[f32], t: &[f32], w: usize) -> f64 {
    let v: Vec<f64> = p.par_chunks(w).zip(t.par_chunks(w)).map(|(x, y)| pearson(&x.iter().map(|&v| v as f64).collect::<Vec<_>>(), &y.iter().map(|&v| v as f64).collect::<Vec<_>>())).filter(|v| v.is_finite()).collect();
    v.iter().sum::<f64>() / v.len().max(1) as f64
}
fn project(z: &[f32], g: usize, v: &[f32], kp: usize) -> Vec<f32> { let n = z.len() / g; sgemm(false, true, n, g, kp, &sgemm(false, false, n, kp, g, z, v), v) }
fn first_cols(v: &[f32], g: usize, k: usize, kp: usize) -> Vec<f32> { (0..g).flat_map(|j| v[j * k..j * k + kp].to_vec()).collect() }

/// panel, ridge weights and programs from the reference perturbations `refs` (penalty, and the number of programs
/// when `kfix` is None, chosen on an inner holdout of the reference)
pub(crate) fn fit_model(d: &Data, refs: &[usize], group: &[i64], kfix: Option<usize>) -> Result<Model> {
    let g = d.genes.len();
    let (itr, iho) = inner(refs, group);
    // one Gram matrix in memory at a time: panel on all reference perturbations, penalty and programs on the inner
    // training set, then weights and programs on all reference perturbations
    let s = { let c = sgram(&d.rows(0, refs), refs.len(), g); greedy32(&c, refs.len(), g, 500) };
    let ins: HashSet<usize> = s.iter().copied().collect();
    let u: Vec<usize> = (0..g).filter(|j| !ins.contains(j)).collect();
    let (alpha, kp) = {
        let c_in = sgram(&d.rows(0, &itr), itr.len(), g);
        let (h1, h2) = (d.rows(1, &iho), d.rows(2, &iho));
        let (xh, th) = (cols(&h1, g, &s), cols(&h2, g, &u));
        let mut best = (f64::NEG_INFINITY, ALPHAS[0]);
        for &al in &ALPHAS {
            let r = mean_r(&sgemm(false, true, iho.len(), u.len(), s.len(), &xh, &ridge32(&c_in, g, &s, &u, al)?), &th, u.len());
            if r > best.0 { best = (r, al); }
        }
        let kp = match kfix {
            Some(k) => k,
            None => {
                let vin = top_eig32(&c_in, g, 160, 30)?;
                let mut bk = (f64::NEG_INFINITY, 160);
                for &k in &K_INTERNAL { let r = mean_r(&project(&h1, g, &first_cols(&vin, g, 160, k), k), &h2, g); if r > bk.0 { bk = (r, k); } }
                bk.1
            }
        };
        (best.1, kp)
    };
    let c = sgram(&d.rows(0, refs), refs.len(), g);
    let w = ridge32(&c, g, &s, &u, alpha)?;
    let v160 = top_eig32(&c, g, 160, 30)?;
    Ok(Model { v: first_cols(&v160, g, 160, kp), s, u, w, kp, alpha })
}
fn hybrid(z: &[f32], g: usize, m: &Model) -> Vec<f32> {
    let (nt, k, nu) = (z.len() / g, m.s.len(), m.u.len());
    let xp = cols(z, g, &m.s);
    let imp = sgemm(false, true, nt, nu, k, &xp, &m.w);
    let mut h = vec![0.0f32; nt * g];
    h.par_chunks_mut(g).enumerate().for_each(|(b, r)| { for (q, &j) in m.s.iter().enumerate() { r[j] = xp[b * k + q]; } for (q, &j) in m.u.iter().enumerate() { r[j] = imp[b * nu + q]; } });
    h
}
/// Measured profile shrunk towards its projection. The part of a z-score profile outside the reference programs has
/// unit noise variance, and its mean square minus one estimates the variance of the signal there. This signal variance
/// is taken proportional to the number of cells, with one constant c for the screen, and the measured part receives
/// the weight w = c n / (c n + 1) for a perturbation with n cells.
pub(crate) fn shrink(z: &[f32], proj: &[f32], g: usize, kp: usize, n_in: &[usize]) -> (Vec<f32>, f64) {
    let ms: Vec<f64> = z.par_chunks(g).zip(proj.par_chunks(g)).map(|(a, b)| a.iter().zip(b).map(|(x, y)| ((x - y) as f64).powi(2)).sum::<f64>() / (g - kp) as f64).collect();
    let c = (ms.iter().map(|m| m - 1.0).sum::<f64>() / n_in.iter().sum::<usize>() as f64).max(0.0);
    let mut o = proj.to_vec();
    o.par_chunks_mut(g).enumerate().for_each(|(i, r)| { let s2 = c * n_in[i] as f64; let w = (s2 / (s2 + 1.0)) as f32; for (x, y) in r.iter_mut().zip(&z[i * g..(i + 1) * g]) { *x += w * (y - *x); } });
    (o, c)
}

struct Row { fold: i64, pert: String, n_in: usize, readout: usize, r: f64, rp: f64, rp_own: f64, nde: usize, tp: usize, disc: usize }
/// results of one input (all cells of the input half, or its first m cells)
struct Acc { label: String, m: Option<usize>, thin: Option<Mmap>, rows: Vec<Row>, path: Vec<Vec<f64>>, pstr: Vec<usize>, pairs: Vec<Vec<(Vec<f64>, Vec<bool>)>> }

pub fn run(a: &[String]) -> Result<()> {
    let dref = Data::open(&req(a, "--ref")?)?;
    let screen = match arg(a, "--screen") { Some(p) => { let x = Data::open(&p)?; if x.genes != dref.genes { bail!("screen genes differ from the reference") } Some((x, p)) } None => None };
    let same_line = screen.is_none() || a.contains(&"--same-line".to_string());
    let (ds, spre): (&Data, String) = match &screen { Some((x, p)) => (x, p.clone()), None => (&dref, req(a, "--ref")?) };
    let (g, np) = (ds.genes.len(), ds.perts.len());
    let out = req(a, "--out")?;
    let read_folds = |d: &Data, ff: &str, hf: &str| -> Result<(Vec<i64>, Vec<i64>, HashSet<(usize, i64)>)> {
        let pi: HashMap<&str, usize> = d.perts.iter().enumerate().map(|(i, p)| (p.as_str(), i)).collect();
        let (mut fold, mut group) = (vec![-1i64; d.perts.len()], vec![-1i64; d.perts.len()]);
        for l in std::fs::read_to_string(ff)?.lines().skip(1) { let x: Vec<&str> = l.split('\t').collect(); if let Some(&i) = pi.get(x[0]) { group[i] = x[1].parse()?; fold[i] = x[2].parse()?; } }
        let mut excl = HashSet::new();
        for l in std::fs::read_to_string(hf)?.lines().skip(1) { let x: Vec<&str> = l.split('\t').collect(); if let Some(&i) = pi.get(x[0]) { excl.insert((i, x[1].parse()?)); } }
        Ok((fold, group, excl))
    };
    let (rfold, rgroup, rexcl) = read_folds(&dref, &req(a, "--ref-folds")?, &req(a, "--ref-hubs")?)?;
    let (sfold, sgroup, sexcl) = if screen.is_some() { read_folds(ds, &req(a, "--folds")?, &req(a, "--hubs")?)? } else { (rfold.clone(), rgroup.clone(), rexcl.clone()) };
    // cells of each perturbation of the screen: the input half holds (n + 1) / 2
    let ncell: Vec<usize> = std::fs::read_to_string(format!("{spre}.meta.tsv"))?.lines().skip(1).map(|l| l.split('\t').nth(1).unwrap().parse::<f64>().unwrap() as usize).collect();
    let min_in: usize = arg(a, "--min-test-in").map(|v| v.parse().unwrap()).unwrap_or(0);
    let min_ref: usize = match arg(a, "--min-ref-cells") { Some(v) => v.parse()?, None => { let mut v = ncell.clone(); v.sort_unstable(); v[v.len() / 2] } };
    let mut tss: HashMap<String, (String, f64)> = HashMap::new();
    if let Some(tf) = arg(a, "--tss") { for (name, v) in read_tss(&tf)? { tss.entry(name).or_insert(v); } }
    let umi = read_umi(arg(a, "--umi"))?;
    // CORUM complexes and gene sets
    let corum = folds::tsv_quoted(&std::fs::read_to_string(req(a, "--corum")?)?);
    let hdr: Vec<&str> = corum[0].iter().map(|x| x.as_str()).collect();
    let (co, cs) = (hdr.iter().position(|h| *h == "organism").unwrap(), hdr.iter().position(|h| *h == "subunits_gene_name").unwrap());
    let all_p: HashSet<&str> = ds.perts.iter().map(|x| x.as_str()).collect();
    let gi: HashMap<&str, usize> = ds.genes.iter().enumerate().map(|(i, x)| (x.as_str(), i)).collect();
    let sets: Vec<Vec<usize>> = std::fs::read_to_string(req(a, "--genesets")?)?.lines().filter_map(|l| { let m: Vec<usize> = l.split('\t').skip(1).filter_map(|n| gi.get(n.trim()).copied()).collect(); if m.len() >= 10 { Some(m) } else { None } }).collect();
    let ns = sets.len();
    let scores = |m: &[f32]| -> Vec<f64> { m.par_chunks(g).flat_map_iter(|r| { let mu = r.iter().map(|&v| v as f64).sum::<f64>() / g as f64; sets.iter().map(move |idx| idx.iter().map(|&j| r[j] as f64).sum::<f64>() / idx.len() as f64 - mu).collect::<Vec<_>>() }).collect() };

    let names = ["direct", "hybrid, external", "hybrid + prior, external", "projected, external", "shrunk, external", "hybrid, internal", "hybrid + prior, internal", "projected, internal", "shrunk, internal"];
    let nr = names.len();
    let new_acc = |label: String, m: Option<usize>, thin: Option<Mmap>| Acc { label, m, thin, rows: vec![], path: vec![vec![]; nr + 1], pstr: vec![], pairs: (0..STRATA.len() + 1).map(|_| (0..nr).map(|_| (vec![], vec![])).collect()).collect() };
    let mut accs: Vec<Acc> = vec![];
    if let Some(pre) = arg(a, "--thin") { for m in req(a, "--thin-cells")?.split(',') { accs.push(new_acc(m.to_string(), Some(m.parse()?), Some(unsafe { Mmap::map(&std::fs::File::open(format!("{pre}.thin{m}.f32"))?)? }))); } }
    accs.push(new_acc("all".into(), None, None));
    let stratum = |n: usize| STRATA.iter().position(|s| n >= s.1 && n <= s.2).unwrap();
    let whole = if same_line { None } else { Some(fit_model(&dref, &(0..dref.perts.len()).collect::<Vec<_>>(), &rgroup, Some(160))?) };
    let mut log = BufWriter::new(std::fs::File::create(format!("{out}.rescue.log"))?);
    for f in 0..5i64 {
        let tests: Vec<usize> = (0..np).filter(|&i| sfold[i] == f && ncell[i] >= 4 && (ncell[i] + 1) / 2 >= min_in).collect();
        if tests.is_empty() { continue; }
        let nt = tests.len();
        let ext_fold;
        let (ext, erefs): (&Model, Vec<usize>) = match &whole {
            Some(m) => (m, (0..dref.perts.len()).collect()),
            None => { let refs: Vec<usize> = (0..dref.perts.len()).filter(|&i| rfold[i] != f && !rexcl.contains(&(i, f))).collect(); ext_fold = fit_model(&dref, &refs, &rgroup, Some(160))?; (&ext_fold, refs) }
        };
        let irefs: Vec<usize> = (0..np).filter(|&i| sfold[i] != f && !sexcl.contains(&(i, f)) && ncell[i] >= min_ref).collect();
        let int = fit_model(ds, &irefs, &sgroup, None)?;
        let z2 = ds.rows(2, &tests);
        let tn: Vec<&str> = tests.iter().map(|&p| ds.perts[p].as_str()).collect();
        // specific response: the 10 leading programs of the external reference removed (rp), or those of the internal one (rp_own)
        let mk_basis = |m: &Model| -> (Vec<f32>, Vec<Vec<f64>>) { let v10 = first_cols(&m.v, g, m.kp, NPC); let b = (0..NPC).map(|c| (0..g).map(|j| v10[j * NPC + c] as f64).collect()).collect(); (v10, b) };
        let ((v10, basis), (_, basis_own)) = (mk_basis(ext), mk_basis(&int));
        let priors = [guide_prior(&dref, &erefs, &ext.s, &ext.u, &tn, &tss, &umi).0, guide_prior(ds, &irefs, &int.s, &int.u, &tn, &tss, &umi).0];
        // complexes of the tested knockdowns of this fold
        let pi: HashMap<&str, usize> = tn.iter().enumerate().map(|(i, p)| (*p, i)).collect();
        let mut member: Vec<Vec<usize>> = vec![vec![]; nt];
        for (ci, l) in corum.iter().skip(1).enumerate() {
            if l.len() <= cs || l[co] != "Human" { continue; }
            let mut nm: Vec<&str> = l[cs].split(';').map(|x| x.trim()).filter(|x| all_p.contains(x)).collect();
            nm.sort_unstable(); nm.dedup();
            if !(2..=50).contains(&nm.len()) { continue; }
            let here: Vec<usize> = nm.iter().filter_map(|n| pi.get(n).copied()).collect();
            if here.len() >= 2 { for &i in &here { member[i].push(ci); } }
        }
        let ann: Vec<usize> = (0..nt).filter(|&i| !member[i].is_empty()).collect();
        let na = ann.len();
        for acc in accs.iter_mut() {
            let z1: Vec<f32> = match &acc.thin { Some(m) => { let t: &[f32] = bytemuck::cast_slice(m); tests.iter().flat_map(|&p| t[p * g..(p + 1) * g].to_vec()).collect() } None => ds.rows(1, &tests) };
            let n_in: Vec<usize> = tests.iter().map(|&p| match acc.m { Some(m) => m.min((ncell[p] + 1) / 2), None => (ncell[p] + 1) / 2 }).collect();
            // one readout in memory at a time
            let models = [(ext, erefs.len()), (&int, irefs.len())];
            for ri in 0..nr {
                let rd: Vec<f32> = if ri == 0 { z1.clone() } else {
                    let (which, kind) = ((ri - 1) / 4, (ri - 1) % 4);
                    let (m, nref) = models[which];
                    match kind {
                        0 => hybrid(&z1, g, m),
                        1 => { let mut hp = hybrid(&z1, g, m); for &(b, j, v) in &priors[which] { hp[b * g + j] = v; } hp }
                        2 => project(&z1, g, &m.v, m.kp),
                        _ => {
                            let (sc, c) = shrink(&z1, &project(&z1, g, &m.v, m.kp), g, m.kp, &n_in);
                            writeln!(log, "fold {f}, input cells {}: reference {nref} perturbations, penalty {}, programs {}, signal variance per cell outside the programs {c:.5}", acc.label, m.alpha, m.kp)?;
                            sc
                        }
                    }
                };
                let m = &rd;
                let per: Vec<((f64, f64, usize, usize, usize), f64)> = (0..nt).into_par_iter().map(|i| { let f64s = |x: &[f32]| x[i * g..(i + 1) * g].iter().map(|&v| v as f64).collect::<Vec<_>>(); let (p, t) = (f64s(m), f64s(&z2)); (sim::metrics(&p, &t, &basis), pearson(&resid(&p, &basis_own), &resid(&t, &basis_own))) }).collect();
                for (i, (x, own)) in per.iter().enumerate() { acc.rows.push(Row { fold: f, pert: tn[i].to_string(), n_in: n_in[i], readout: ri, r: x.0, rp: x.1, rp_own: *own, nde: x.2, tp: x.3, disc: x.4 }); }
                acc.path[ri].extend(scores(m));
                // pairs of tested knockdowns of this fold, both in the same stratum of input cells
                let sub: Vec<f32> = ann.iter().flat_map(|&i| m[i * g..(i + 1) * g].to_vec()).collect();
                let r = row_corr(&deflate_rows(&sub, g, &v10), g);
                for x in 0..na { for y in x + 1..na {
                    let lab = member[ann[x]].iter().any(|c1| member[ann[y]].contains(c1));
                    let (sx, sy) = (stratum(n_in[ann[x]]), stratum(n_in[ann[y]]));
                    acc.pairs[STRATA.len()][ri].0.push(r[x * na + y] as f64); acc.pairs[STRATA.len()][ri].1.push(lab);
                    if sx == sy { acc.pairs[sx][ri].0.push(r[x * na + y] as f64); acc.pairs[sx][ri].1.push(lab); }
                } }
            }
            acc.path[nr].extend(scores(&z2));
            acc.pstr.extend(&n_in);
        }
        eprintln!("fold {f}: {nt} tested, external reference {} (programs {}), internal reference {} (programs {})", erefs.len(), ext.kp, irefs.len(), int.kp);
    }
    // per-perturbation table and summary
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.rescue.tsv"))?);
    writeln!(w, "input_cells\tfold\tpert\tn_in\treadout\tr\trp\trp_own\tn_de\ttp\tdisc")?;
    let mut ws = BufWriter::new(std::fs::File::create(format!("{out}.rescue.summary.tsv"))?);
    writeln!(ws, "input_cells\tstratum\treadout\tn\tmedian_n_in\tmedian_r\tmedian_rp\tmedian_rp_own\tde_recall\topposite\tshare_above_direct\twilcoxon_p\tpathway_r\tcomplex_auroc\tn_partner_pairs")?;
    for acc in &accs {
        for r in &acc.rows { writeln!(w, "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}", acc.label, r.fold, r.pert, r.n_in, names[r.readout], r.r, r.rp, r.rp_own, r.nde, r.tp, r.disc)?; }
        let ntot = acc.pstr.len();
        for si in 0..=STRATA.len() {
            let inside = |n: usize| si == STRATA.len() || stratum(n) == si;
            let label = if si == STRATA.len() { "all" } else { STRATA[si].0 };
            let dir: HashMap<(i64, &str), f64> = acc.rows.iter().filter(|r| r.readout == 0 && inside(r.n_in)).map(|r| ((r.fold, r.pert.as_str()), r.rp)).collect();
            if dir.len() < 20 || (acc.m.is_some() && si < STRATA.len()) { continue; }
            let keep: Vec<usize> = (0..ntot).filter(|&i| inside(acc.pstr[i])).collect();
            for ri in 0..nr {
                let x: Vec<&Row> = acc.rows.iter().filter(|r| r.readout == ri && inside(r.n_in)).collect();
                let (nde, tp, disc) = x.iter().fold((0usize, 0usize, 0usize), |s, r| (s.0 + r.nde, s.1 + r.tp, s.2 + r.disc));
                let d: Vec<f64> = x.iter().map(|r| r.rp - dir[&(r.fold, r.pert.as_str())]).filter(|v| v.is_finite()).collect();
                let pw = median((0..ns).map(|q| pearson(&keep.iter().map(|&i| acc.path[ri][i * ns + q]).collect::<Vec<_>>(), &keep.iter().map(|&i| acc.path[nr][i * ns + q]).collect::<Vec<_>>())).collect());
                let (sc, lab) = &acc.pairs[si][ri];
                writeln!(ws, "{}\t{label}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.3e}\t{:.4}\t{:.4}\t{}", acc.label, names[ri], x.len(), median(x.iter().map(|r| r.n_in as f64).collect()), median(x.iter().map(|r| r.r).collect()), median(x.iter().map(|r| r.rp).collect()), median(x.iter().map(|r| r.rp_own).collect()),
                    tp as f64 / nde as f64, disc as f64 / nde as f64, if ri == 0 { f64::NAN } else { d.iter().filter(|v| **v > 0.0).count() as f64 / d.len() as f64 }, if ri == 0 { f64::NAN } else { upg::wilcoxon(&d) },
                    pw, auc(sc, lab).unwrap_or(f64::NAN), lab.iter().filter(|l| **l).count())?;
            }
        }
    }
    eprintln!("peak memory {:.0} MB", rss_mb());
    Ok(())
}

/// z-scores of a cell cache on the genes of --z, and the input half restricted to its first M cells
pub fn thin(a: &[String]) -> Result<()> {
    let c = Cache::open(&req(a, "--cache")?)?;
    let genes: Vec<String> = std::fs::read_to_string(format!("{}.genes.tsv", req(a, "--z")?))?.lines().skip(1).map(|l| l.split('\t').next().unwrap().to_string()).collect();
    let out = req(a, "--out")?;
    let ms: Vec<usize> = arg(a, "--cells").map(|s| s.split(',').map(|x| x.parse().unwrap()).collect()).unwrap_or(vec![15, 30, 50]);
    let g = genes.len();
    let cgi: HashMap<&str, usize> = c.genes.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let mut map = vec![-1i32; c.genes.len()];
    for (j, s) in genes.iter().enumerate() { map[*cgi.get(s.as_str()).ok_or_else(|| anyhow!("gene {s} not in cache"))?] = j as i32; }
    let ctrl = c.group("control").ok_or_else(|| anyhow!("no control"))?;
    let allq: Vec<usize> = std::iter::once(ctrl).chain((0..c.groups.len()).filter(|&i| i != ctrl && c.groups[i].2 >= 4)).collect();
    let nq = allq.len() - 1;
    let sums = sim::pseudobulk(&c, &allq, &map, g, &vec![1.0; c.genes.len()], true, 1.0, 0);
    let (z0, z1) = sim::zhalves(&c, &allq, &sums, g);
    // all cells: mean of the perturbation against the control mean, with the null variance of the halves
    let n = |s: usize| c.groups[allq[s]].2 as f64;
    let zall: Vec<f64> = (0..nq * g).into_par_iter().map(|x| {
        let (s, j) = (x / g + 1, x % g);
        let (n0, n1) = (((n(s) + 1.0) / 2.0).floor(), (n(s) / 2.0).floor());
        let (m0, m1, mc) = (sums[(s * 2) * g + j] / n0, sums[(s * 2 + 1) * g + j] / n1, (sums[j] + sums[g + j]) / n(0));
        // the standard error of half 0 gives the null variance: z0 = (m0 - mc) / se0
        let se0 = if z0[x] != 0.0 { (m0 - mc) / z0[x] } else { return 0.0 };
        let s2 = se0 * se0 / (1.0 / n0 + 1.0 / n(0));
        ((m0 * n0 + m1 * n1) / n(s) - mc) / (s2 * (1.0 / n(s) + 1.0 / n(0))).sqrt()
    }).collect();
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.meta.tsv"))?);
    writeln!(w, "pert\tn")?;
    for s in 1..=nq { writeln!(w, "{}\t{}", c.groups[allq[s]].0, c.groups[allq[s]].2)?; }
    std::fs::write(format!("{out}.genes.tsv"), format!("gene\n{}\n", genes.join("\n")))?;
    let f32s = |v: &[f64]| -> Vec<f32> { v.iter().map(|&x| x as f32).collect() };
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.z.f32"))?);
    for layer in [&zall, &z0, &z1] { w.write_all(bytemuck::cast_slice(&f32s(layer)))?; }
    for &m in &ms { std::fs::write(format!("{out}.thin{m}.f32"), bytemuck::cast_slice::<f32, u8>(&f32s(&real::zsub(&c, &allq, &map, g, false, Some(m), 0))))?; }
    eprintln!("{nq} perturbations, {g} genes, input cells {ms:?}");
    Ok(())
}
