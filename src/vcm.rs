//! Measurement versus context transfer (virtual cell question). For a knockdown X in a target
//! screen, compare on the held-out half of the target cells:
//!   copy      zero-shot transfer: the profile of X measured in another atlas (no target data)
//!   panel     500 greedy MMSE genes measured in the target input half, other genes imputed by ridge
//!   panel+copy panel imputation plus beta times the gene-specific part of the copied profile
//!             (its residual after the same ridge map from its own panel genes); beta is fitted on
//!             a calibration fifth of the target knockdowns and reported on the other four fifths
//!   direct    the target input half measured for all genes (replicate ceiling)
//! Specific response: residual after the 10 leading programs of the target input half.
use crate::*;

struct Src { name: String, z: Z, idx: HashMap<String, usize>, rms: Vec<f64> }

fn load_src(name: &str, pre: &str) -> Result<Src> {
    let z = Z::load(pre)?;
    let g = z.genes.len();
    let idx = z.perts.iter().enumerate().map(|(i, p)| (p.clone(), i)).collect();
    let rms = (0..z.perts.len()).map(|i| ((0..g).map(|j| z.at(0, i, j).powi(2)).sum::<f64>() / g as f64).sqrt().max(1e-9)).collect();
    Ok(Src { name: name.into(), z, idx, rms })
}

/// DE metrics of a prediction p against truth t with nc calls: (n_de, tp, disc, bins[(n, tp)])
fn de_eval(p: &[f64], t: &[f64], nc: usize) -> (usize, usize, usize, [(usize, usize); 3]) {
    let mut o: Vec<usize> = (0..p.len()).collect();
    o.sort_by(|&a, &b| p[b].abs().partial_cmp(&p[a].abs()).unwrap().then(a.cmp(&b)));
    let bin = |z: f64| if z < 6.0 { 0 } else if z < 10.0 { 1 } else { 2 };
    let mut bins = [(0usize, 0usize); 3];
    let mut nde = 0;
    for j in 0..t.len() { if t[j].abs() > DE_Z { nde += 1; bins[bin(t[j].abs())].0 += 1; } }
    let (mut tp, mut disc) = (0, 0);
    for &j in &o[..nc.min(o.len())] {
        let same = p[j].signum() == t[j].signum();
        if !same { disc += 1; }
        if same && t[j].abs() > DE_Z { tp += 1; bins[bin(t[j].abs())].1 += 1; }
    }
    (nde, tp, disc, bins)
}

pub fn run(a: &[String]) -> Result<()> {
    let t0 = std::time::Instant::now();
    let (rname, rpre) = req(a, "--ref")?.split_once('=').map(|(k, v)| (k.to_string(), v.to_string())).unwrap();
    let rz = Z::load(&rpre)?;
    let g = rz.genes.len();
    let k: usize = arg(a, "--panel-size").map(|v| v.parse().unwrap()).unwrap_or(500);
    let srcs: Vec<Src> = req(a, "--sources")?.split(',').map(|kv| { let (n, p) = kv.split_once('=').unwrap(); load_src(n, p) }).collect::<Result<_>>()?;
    for s in &srcs { if s.z.genes != rz.genes { bail!("gene order differs for source {}", s.name) } }
    // reference model: greedy MMSE panel and ridge weights on the reference second-moment matrix
    let all_r: Vec<usize> = (0..rz.perts.len()).collect();
    let mut sig = la::gram(&rz.mat(0, &all_r, &(0..g).collect::<Vec<_>>()), all_r.len(), g);
    sig.par_iter_mut().for_each(|v| *v /= all_r.len() as f64);
    let s2 = all_r.par_iter().map(|&i| (0..g).map(|j| (rz.at(1, i, j) - rz.at(2, i, j)).powi(2)).sum::<f64>()).sum::<f64>() / (2.0 * (all_r.len() * g) as f64);
    let gr = Grams { c: sig.iter().map(|v| v / s2.max(1.0)).collect(), c_in: vec![], v: vec![], vk: vec![], ev: vec![] };
    let s = sim::greedy(&gr, 1, g, k, &vec![1.0; g]);
    let inp: HashSet<usize> = s.iter().copied().collect();
    let u: Vec<usize> = (0..g).filter(|j| !inp.contains(j)).collect();
    let nu = u.len();
    let w = ridge(&sig, g, &s, &u, 1.0)?; // nu x k
    eprintln!("vcm: reference {rname} {} perturbations, panel {k}, {:.1}s", all_r.len(), t0.elapsed().as_secs_f64());
    let mut out = BufWriter::new(std::fs::File::create(format!("{}.vcm.tsv", req(a, "--out")?))?);
    writeln!(out, "target\tpert\tsplit\tmethod\tscope\tr\trperp\trperp_ceil\tnde\tncall\ttp\tdisc\tspa\tn46\ttp46\tn610\ttp610\tn10\ttp10\tbeta")?;
    for kv in req(a, "--targets")?.split(',') {
        let (tname, tpre) = kv.split_once('=').unwrap();
        let tg = std::time::Instant::now();
        let zt = Z::load(tpre)?;
        if zt.genes != rz.genes { bail!("gene order differs for target {tname}") }
        let ss: Vec<&Src> = srcs.iter().filter(|x| x.name != tname).collect();
        // target frame: 10 leading programs of the input-half profiles (subspace iteration)
        let all_t: Vec<usize> = (0..zt.perts.len()).collect();
        let (v10, _) = la::top_eig_sub(&la::gram(&zt.mat(1, &all_t, &(0..g).collect::<Vec<_>>()), all_t.len(), g), g, NPC, 40)?;
        let bu = basis_on(&v10, NPC, &u);
        let ball = basis_on(&v10, NPC, &(0..g).collect::<Vec<_>>());
        // tests: target knockdowns also present in the reference atlas (for the K562 copy)
        let tests: Vec<usize> = all_t.iter().copied().filter(|&i| ss.iter().any(|x| x.name == rname && x.idx.contains_key(&zt.perts[i]))).collect();
        let cal = |i: usize| -> bool { let h = zt.perts[i].bytes().fold(0xcbf29ce484222325u64, |h, b| (h ^ b as u64).wrapping_mul(0x100000001b3)); h % 5 == 0 };
        // per test: panel imputation and copied profiles (single reference atlas and pooled over sources)
        struct P { imp: Vec<f64>, copy: Vec<Vec<f64>>, rc: Vec<Vec<f64>> }
        let copies = |i: usize| -> Vec<Vec<f64>> {
            let name = &zt.perts[i];
            let one = |x: &Src| -> Option<Vec<f64>> { x.idx.get(name).map(|&q| (0..g).map(|j| x.z.at(0, q, j) / x.rms[q]).collect()) };
            let single = ss.iter().find(|x| x.name == rname).and_then(|x| one(x)).unwrap();
            let avail: Vec<Vec<f64>> = ss.iter().filter_map(|x| one(x)).collect();
            let pooled: Vec<f64> = (0..g).map(|j| avail.iter().map(|v| v[j]).sum::<f64>() / avail.len() as f64).collect();
            vec![single, pooled]
        };
        let preds: Vec<P> = tests.par_iter().map(|&i| {
            let x: Vec<f64> = s.iter().map(|&j| zt.at(1, i, j)).collect();
            let imp: Vec<f64> = (0..nu).map(|q| (0..k).map(|r| w[q * k + r] * x[r]).sum()).collect();
            let copy = copies(i);
            let rc = copy.iter().map(|c| (0..nu).map(|q| c[u[q]] - (0..k).map(|r| w[q * k + r] * c[s[r]]).sum::<f64>()).collect()).collect();
            P { imp, copy, rc }
        }).collect();
        // beta for panel+copy on the calibration fifth (genes outside the panel, held-out half as truth)
        let beta: Vec<f64> = (0..2).map(|m| {
            let (mut num, mut den) = (0.0, 0.0);
            for (p, &i) in preds.iter().zip(&tests) { if !cal(i) { continue; } for q in 0..nu { let t = zt.at(2, i, u[q]); num += (t - p.imp[q]) * p.rc[m][q]; den += p.rc[m][q].powi(2); } }
            num / den.max(1e-12)
        }).collect();
        // two-coefficient least squares a * imputed + b * copied profile, on the calibration fifth
        let ab: Vec<(f64, f64)> = (0..2).map(|m| {
            let (mut s11, mut s12, mut s22, mut r1, mut r2) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for (p, &i) in preds.iter().zip(&tests) { if !cal(i) { continue; } for q in 0..nu { let (x1, x2, t) = (p.imp[q], p.copy[m][u[q]], zt.at(2, i, u[q])); s11 += x1 * x1; s12 += x1 * x2; s22 += x2 * x2; r1 += x1 * t; r2 += x2 * t; } }
            let det = s11 * s22 - s12 * s12;
            ((r1 * s22 - r2 * s12) / det, (r2 * s11 - r1 * s12) / det)
        }).collect();
        eprintln!("vcm {tname}: least squares (a, b) single {:?} pooled {:?}", ab[0], ab[1]);
        eprintln!("vcm {tname}: {} tests ({} calibration), sources {:?}, beta single {:.3} pooled {:.3}", tests.len(), tests.iter().filter(|&&i| cal(i)).count(), ss.iter().map(|x| x.name.clone()).collect::<Vec<_>>(), beta[0], beta[1]);
        let rows: Vec<String> = preds.par_iter().zip(&tests).flat_map(|(p, &i)| {
            let split = if cal(i) { "cal" } else { "eval" };
            let tu: Vec<f64> = u.iter().map(|&j| zt.at(2, i, j)).collect();
            let du: Vec<f64> = u.iter().map(|&j| zt.at(1, i, j)).collect();
            let ta: Vec<f64> = (0..g).map(|j| zt.at(2, i, j)).collect();
            let da: Vec<f64> = (0..g).map(|j| zt.at(1, i, j)).collect();
            let (rtu, rta) = (resid(&tu, &bu), resid(&ta, &ball));
            let (ceil_u, ceil_a) = (pearson(&resid(&du, &bu), &rtu), pearson(&resid(&da, &ball), &rta));
            let (ncu, nca) = (du.iter().filter(|v| v.abs() > DE_Z).count(), da.iter().filter(|v| v.abs() > DE_Z).count());
            let mut o = vec![];
            let mut push = |m: &str, scope: &str, y: &[f64], bt: f64| {
                let (t, rt, b, ceil, nc) = if scope == "out" { (&tu, &rtu, &bu, ceil_u, ncu) } else { (&ta, &rta, &ball, ceil_a, nca) };
                let (nde, tp, disc, bins) = de_eval(y, t, nc);
                o.push(format!("{tname}\t{}\t{split}\t{m}\t{scope}\t{:.5}\t{:.5}\t{ceil:.5}\t{nde}\t{nc}\t{tp}\t{disc}\t{:.5}\t{}\t{}\t{}\t{}\t{}\t{}\t{bt:.4}",
                    zt.perts[i], pearson(y, t), pearson(&resid(y, b), rt), spa_de(y, t), bins[0].0, bins[0].1, bins[1].0, bins[1].1, bins[2].0, bins[2].1));
            };
            // genes outside the panel
            push("direct", "out", &du, 0.0);
            push("panel", "out", &p.imp, 0.0);
            for (m, lab) in [(0, "copy_single"), (1, "copy_pooled")] {
                let cu: Vec<f64> = u.iter().map(|&j| p.copy[m][j]).collect();
                push(lab, "out", &cu, 0.0);
                let comb: Vec<f64> = (0..nu).map(|q| p.imp[q] + beta[m] * p.rc[m][q]).collect();
                push(&format!("panel+{lab}"), "out", &comb, beta[m]);
                let ls: Vec<f64> = (0..nu).map(|q| ab[m].0 * p.imp[q] + ab[m].1 * p.copy[m][u[q]]).collect();
                push(&format!("panel_ls_{lab}"), "out", &ls, ab[m].1);
            }
            // all genes: zero-shot copy, hybrid panel (measured panel genes + imputed), hybrid + copy
            push("direct", "all", &da, 0.0);
            let mut hyb = da.clone();
            for (q, &j) in u.iter().enumerate() { hyb[j] = p.imp[q]; }
            push("panel", "all", &hyb, 0.0);
            for (m, lab) in [(0, "copy_single"), (1, "copy_pooled")] {
                push(lab, "all", &p.copy[m], 0.0);
                let mut h = hyb.clone();
                for (q, &j) in u.iter().enumerate() { h[j] += beta[m] * p.rc[m][q]; }
                push(&format!("panel+{lab}"), "all", &h, beta[m]);
                let mut h2 = hyb.clone();
                for (q, &j) in u.iter().enumerate() { h2[j] = ab[m].0 * p.imp[q] + ab[m].1 * p.copy[m][j]; }
                push(&format!("panel_ls_{lab}"), "all", &h2, ab[m].1);
            }
            o
        }).collect();
        for r in rows { writeln!(out, "{r}")?; }
        out.flush()?;
        eprintln!("vcm {tname}: done in {:.1}s", tg.elapsed().as_secs_f64());
    }
    Ok(())
}

/// summary on the evaluation knockdowns: median specific correlation (specific responders,
/// replicate specific correlation >= 0.2), DE recall and opposite-sign rate, size ranking,
/// recall by strength, and paired comparisons of panel with copy and of panel+copy with panel
pub fn summary(a: &[String]) -> Result<()> {
    let txt = std::fs::read_to_string(format!("{}.vcm.tsv", req(a, "--in")?))?;
    let rows: Vec<Vec<&str>> = txt.lines().skip(1).map(|l| l.split('\t').collect()).collect();
    let f = |r: &Vec<&str>, i: usize| -> f64 { r[i].parse().unwrap_or(f64::NAN) };
    let mut targets: Vec<&str> = vec![];
    for r in &rows { if !targets.contains(&r[0]) { targets.push(r[0]); } }
    let methods = ["direct", "copy_single", "copy_pooled", "panel", "panel+copy_single", "panel+copy_pooled", "panel_ls_copy_pooled"];
    for scope in ["out", "all"] {
        println!("== scope {scope} (out: genes outside the panel, all: all genes)");
        println!("target   method              n    spec_r  recall  opp    rank   weak   mid    strong");
        for t in &targets {
            for m in methods {
                let s: Vec<&Vec<&str>> = rows.iter().filter(|r| r[0] == *t && r[2] == "eval" && r[3] == m && r[4] == scope).collect();
                if s.is_empty() { continue; }
                let resp: Vec<&&Vec<&str>> = s.iter().filter(|r| f(r, 7) >= 0.2).collect();
                let sum = |i: usize| s.iter().map(|r| f(r, i)).sum::<f64>();
                println!("{t:8} {m:19} {:4} {:.3}   {:.3}   {:.3}  {:.3}  {:.3}  {:.3}  {:.3}", resp.len(), median(resp.iter().map(|r| f(r, 6)).collect()),
                    sum(10) / sum(8), sum(11) / sum(9).max(1.0), median(s.iter().map(|r| f(r, 12)).collect()), sum(14) / sum(13), sum(16) / sum(15), sum(18) / sum(17));
            }
        }
    }
    println!("== paired specific correlation, genes outside the panel, specific responders (share A higher, Wilcoxon P)");
    for t in &targets {
        let get = |m: &str| -> HashMap<&str, (f64, f64)> { rows.iter().filter(|r| r[0] == *t && r[2] == "eval" && r[3] == m && r[4] == "out").map(|r| (r[1], (f(r, 6), f(r, 7)))).collect() };
        for (ma, mb) in [("panel", "copy_single"), ("panel", "copy_pooled"), ("panel+copy_pooled", "panel"), ("panel_ls_copy_pooled", "panel"), ("panel_ls_copy_pooled", "copy_pooled")] {
            let (x, y) = (get(ma), get(mb));
            let d: Vec<f64> = x.iter().filter(|(_, v)| v.1 >= 0.2).filter_map(|(p, v)| y.get(p).map(|w| v.0 - w.0)).filter(|v| v.is_finite()).collect();
            let hi = d.iter().filter(|v| **v > 0.0).count() as f64 / d.len().max(1) as f64;
            println!("{t:8} {ma:18} vs {mb:12} n {:4}  {hi:.2}  P {:.1e}", d.len(), crate::upg::wilcoxon(&d));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    /// subspace iteration matches the exact LAPACK eigenpairs on a spiked covariance
    #[test]
    fn subspace_matches_exact() {
        let (n, p) = (400, 3000);
        let mut st = 7u64;
        let mut x: Vec<f64> = (0..p * n).map(|_| (sm64(&mut st) >> 11) as f64 / (1u64 << 53) as f64 - 0.5).collect();
        for i in 0..p { for f in 0..12 { let a = ((sm64(&mut st) >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * (12 - f) as f64; for j in 0..n { x[i * n + j] += a * (((f * 31 + j * 7) % 13) as f64 - 6.0) / 6.0; } } }
        let c = la::gram(&x, p, n);
        let (ve, le) = la::top_eig(&c, n, 10).unwrap();
        let (vs, ls) = la::top_eig_sub(&c, n, 10, 40).unwrap();
        for j in 0..10 {
            assert!((le[j] - ls[j]).abs() / le[j] < 1e-9, "eigenvalue {j}: {} vs {}", le[j], ls[j]);
            let d: f64 = (0..n).map(|i| ve[i * 10 + j] * vs[i * 10 + j]).sum();
            assert!(d.abs() > 1.0 - 1e-8, "eigenvector {j}: |cos| {}", d.abs());
        }
    }
}
