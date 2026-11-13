//! Biological inference from imputed profiles (test perturbations of one cross-validation fold):
//!   panelimp bio --z PREFIX --folds F --hubs H --corum FILE --genesets FILE --out PREFIX [--fold 0] [--target PREFIX]
//! With --target, every knockdown of another screen is tested with the whole reference atlas (--fold -9).
//! A. Protein complexes: do two knockdowns of the same CORUM complex have similar profiles
//!    (area under the ROC curve over pairs, and share of knockdowns whose nearest neighbours are complex partners)?
//! B. Pathways: score of each gene set (mean z of its genes minus the mean of all genes) against the same score
//!    in the held-out half of the cells.
//! C. Strength: number of DE genes of each knockdown against the held-out half.
//! Readouts of the input half: whole transcriptome as measured, projected on 160 reference programs, the 500 panel
//! genes alone, and the hybrid profile (panel genes measured, other genes imputed by ridge regression).
use crate::causal::*;
use crate::explore::*;
use crate::*;

pub(crate) fn deflate_rows(m: &[f32], w: usize, b: &[f32]) -> Vec<f32> { let (n, k) = (m.len() / w, b.len() / w); let mut o = m.to_vec(); let back = sgemm(false, true, n, w, k, &sgemm(false, false, n, k, w, m, b), b); o.par_iter_mut().zip(&back).for_each(|(x, y)| *x -= y); o }
/// correlation matrix between the rows of m (n x w)
pub(crate) fn row_corr(m: &[f32], w: usize) -> Vec<f32> {
    let n = m.len() / w;
    let mut c = m.to_vec();
    c.par_chunks_mut(w).for_each(|r| { let mu = r.iter().sum::<f32>() / w as f32; for v in r.iter_mut() { *v -= mu; } let nr = r.iter().map(|v| v * v).sum::<f32>().sqrt().max(1e-12); for v in r.iter_mut() { *v /= nr; } });
    sgemm(false, true, n, n, w, &c, &c)
}

pub(crate) fn read_umi(f: Option<String>) -> Result<HashMap<String, f64>> {
    Ok(match f { Some(uf) => std::fs::read_to_string(uf)?.lines().filter_map(|l| { let x: Vec<&str> = l.split('\t').collect(); x.get(1).and_then(|v| v.parse::<f64>().ok()).map(|v| (x[0].to_string(), v)) }).collect(), None => HashMap::new() })
}
/// Guide-target prior: entries (row of `names`, gene, value) for the target gene of each knockdown (median z of the
/// reference knockdowns in its expression decile, held-out half) and for genes with a TSS within 2 kb (one median),
/// when these genes are outside the panel `s`. Also returns the two overall medians.
pub(crate) fn guide_prior(d: &Data, refs: &[usize], s: &[usize], u: &[usize], names: &[&str], tss: &HashMap<String, (String, f64)>, umi: &HashMap<String, f64>) -> (Vec<(usize, usize, f32)>, f32, f32) {
    let g = d.genes.len();
    let gix: HashMap<&str, usize> = d.genes.iter().enumerate().map(|(i, x)| (x.as_str(), i)).collect();
    let inp: HashSet<usize> = s.iter().copied().collect();
    let own = |name: &str| -> Option<usize> { gix.get(name).copied().filter(|j| !inp.contains(j)) };
    let near = |name: &str| -> Vec<usize> { match tss.get(name) { Some((chr, pos)) => u.iter().copied().filter(|&j| d.genes[j] != name && tss.get(&d.genes[j]).map(|(c2, p2)| c2 == chr && (p2 - pos).abs() < 2e3).unwrap_or(false)).collect(), None => vec![] } };
    let r2 = d.rows(2, refs);
    let (mut zo, mut zn, mut ze) = (vec![], vec![], vec![]);
    for (b, &p) in refs.iter().enumerate() {
        if let Some(j) = own(&d.perts[p]) { zo.push(r2[b * g + j] as f64); if let Some(&e) = umi.get(&d.perts[p]) { ze.push((e, r2[b * g + j] as f64)); } }
        for j in near(&d.perts[p]) { zn.push(r2[b * g + j] as f64); }
    }
    ze.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
    let (mut cuts, mut bins) = (vec![], vec![]);
    if ze.len() >= 100 { for ch in ze.chunks(ze.len().div_ceil(10)) { cuts.push(ch.last().unwrap().0); bins.push(median(ch.iter().map(|x| x.1).collect()) as f32); } }
    let (own_med, near_med) = (median(zo) as f32, if zn.is_empty() { 0.0 } else { median(zn) as f32 });
    let mut set = vec![];
    for (b, name) in names.iter().enumerate() {
        if let Some(j) = own(name) { set.push((b, j, match umi.get(*name) { Some(&e) if !cuts.is_empty() => bins[cuts.iter().position(|&cc| e <= cc).unwrap_or(cuts.len() - 1)], _ => own_med })); }
        for j in near(name) { set.push((b, j, near_med)); }
    }
    (set, own_med, near_med)
}

pub fn run(a: &[String]) -> Result<()> {
    let c = setup(a)?;
    let (d, s, u) = (&c.d, &c.s, &c.u);
    let (g, nu, k) = (d.genes.len(), u.len(), s.len());
    let out = req(a, "--out")?;
    let f: i64 = arg(a, "--fold").map(|v| v.parse().unwrap()).unwrap_or(0);
    let target = match arg(a, "--target") { Some(t) => { let x = Data::open(&t)?; if x.genes != d.genes { bail!("target genes differ from the reference") } Some(x) } None => None };
    let dt: &Data = target.as_ref().unwrap_or(d);
    let tests_v: Vec<usize> = if target.is_some() { (0..dt.perts.len()).collect() } else { c.tests.clone() };
    let mut w = BufWriter::new(std::fs::File::create(if target.is_some() { format!("{out}.bio.tsv") } else { format!("{out}.bio.fold{f}.tsv") })?);
    writeln!(w, "analysis\tframe\treadout\tmetric\tn\tvalue")?;
    macro_rules! rec { ($an:expr, $fr:expr, $ro:expr, $me:expr, $n:expr, $v:expr) => {{ let val: f64 = $v; println!("{:10} {:9} {:28} {:22} {:>7} {:.4}", $an, $fr, $ro, $me, $n, val); writeln!(w, "{}\t{}\t{}\t{}\t{}\t{:.6}", $an, $fr, $ro, $me, $n, val)?; }} }
    // reference model
    let c0 = sgram(&d.rows(0, &c.refs), c.refs.len(), g);
    let w0 = ridge32(&c0, g, s, u, c.alpha)?;
    let v160 = top_eig32(&c0, g, 160, 30)?;
    drop(c0);
    let v10: Vec<f32> = (0..g).flat_map(|j| v160[j * 160..j * 160 + NPC].to_vec()).collect();
    let v10p: Vec<f32> = { let v64: Vec<f64> = v10.iter().map(|&x| x as f64).collect(); let b = basis_on(&v64, NPC, s); (0..k).flat_map(|q| b.iter().map(|col| col[q] as f32).collect::<Vec<_>>()).collect() };
    let tests = &tests_v;
    let nt = tests.len();
    let (z0, z1, z2) = (dt.rows(0, tests), dt.rows(1, tests), dt.rows(2, tests));
    let xp = cols(&z1, g, s);
    let imp = sgemm(false, true, nt, nu, k, &xp, &w0);
    let mut hyb = vec![0.0f32; nt * g];
    hyb.par_chunks_mut(g).enumerate().for_each(|(b, r)| { for (q, &j) in s.iter().enumerate() { r[j] = xp[b * k + q]; } for (q, &j) in u.iter().enumerate() { r[j] = imp[b * nu + q]; } });
    let proj = sgemm(false, true, nt, g, 160, &sgemm(false, false, nt, 160, g, &z1, &v160), &v160);
    // hybrid profile with the guide-target prior: the target gene of each knockdown (by expression decile) and genes with a
    // TSS within 2 kb take the median z-score of the reference knockdowns (held-out half), when they are outside the panel
    let hybp = {
        let mut tss: HashMap<String, (String, f64)> = HashMap::new();
        if let Some(tf) = arg(a, "--tss") { for (name, v) in read_tss(&tf)? { tss.entry(name).or_insert(v); } }
        let umi = read_umi(arg(a, "--umi"))?;
        let names: Vec<&str> = tests.iter().map(|&p| dt.perts[p].as_str()).collect();
        let (set, own_med, near_med) = guide_prior(d, &c.refs, s, u, &names, &tss, &umi);
        let mut h = hyb.clone();
        for &(b, j, v) in &set { h[b * g + j] = v; }
        println!("guide-target prior: {} entries set (target median {own_med:.2}, neighbours {near_med:.2})", set.len());
        h
    };
    let readouts: Vec<(&str, &[f32], usize)> = vec![("whole transcriptome, all cells", &z0, g), ("whole transcriptome, half", &z1, g), ("projected on 160 programs", &proj, g), ("panel genes alone", &xp, k), ("hybrid (panel + imputed)", &hyb, g), ("hybrid + guide-target prior", &hybp, g)];

    // ---------------- A. complexes
    let pi: HashMap<&str, usize> = tests.iter().enumerate().map(|(i, &p)| (dt.perts[p].as_str(), i)).collect();
    let rows = folds::tsv_quoted(&std::fs::read_to_string(req(a, "--corum")?)?);
    let hdr: Vec<&str> = rows[0].iter().map(|x| x.as_str()).collect();
    let (co, cs) = (hdr.iter().position(|h| *h == "organism").unwrap(), hdr.iter().position(|h| *h == "subunits_gene_name").unwrap());
    let all_p: HashSet<&str> = dt.perts.iter().map(|x| x.as_str()).collect();
    let mut member: Vec<Vec<usize>> = vec![vec![]; nt];
    let mut ncx = 0;
    for (ci, l) in rows.iter().skip(1).enumerate() {
        if l.len() <= cs || l[co] != "Human" { continue; }
        let mut names: Vec<&str> = l[cs].split(';').map(|x| x.trim()).filter(|x| all_p.contains(x)).collect();
        names.sort_unstable(); names.dedup();
        if !(2..=50).contains(&names.len()) { continue; }
        let here: Vec<usize> = names.iter().filter_map(|n| pi.get(n).copied()).collect();
        if here.len() >= 2 { ncx += 1; for &i in &here { member[i].push(ci); } }
    }
    let ann: Vec<usize> = (0..nt).filter(|&i| !member[i].is_empty()).collect();
    let na = ann.len();
    let partner = |x: usize, y: usize| member[ann[x]].iter().any(|c1| member[ann[y]].contains(c1));
    let npos = (0..na).map(|x| (x + 1..na).filter(|&y| partner(x, y)).count()).sum::<usize>();
    println!("fold {f}: {nt} test knockdowns, {na} in {ncx} complexes with two or more members tested, {npos} partner pairs of {}", na * (na - 1) / 2);
    for (name, m, wd) in &readouts {
        let sub: Vec<f32> = ann.iter().flat_map(|&i| m[i * wd..(i + 1) * wd].to_vec()).collect();
        for (frame, mat) in [("all", sub.clone()), ("specific", deflate_rows(&sub, *wd, if *wd == g { &v10 } else { &v10p }))] {
            let r = row_corr(&mat, *wd);
            let (mut sc, mut lab) = (vec![], vec![]);
            for x in 0..na { for y in x + 1..na { sc.push(r[x * na + y] as f64); lab.push(partner(x, y)); } }
            rec!("complex", frame, name, "pair AUROC", sc.len(), auc(&sc, &lab).unwrap_or(f64::NAN));
            let (mut t1, mut t5, mut nn) = (0usize, 0usize, 0usize);
            for x in 0..na {
                let mut o: Vec<usize> = (0..na).filter(|&y| y != x).collect();
                o.sort_by(|&p, &q| r[x * na + q].partial_cmp(&r[x * na + p]).unwrap());
                nn += 1; if partner(x, o[0]) { t1 += 1; } if o[..5].iter().any(|&y| partner(x, y)) { t5 += 1; }
            }
            rec!("complex", frame, name, "nearest is partner", nn, t1 as f64 / nn as f64);
            rec!("complex", frame, name, "partner in nearest 5", nn, t5 as f64 / nn as f64);
        }
    }
    // ---------------- B. pathways (truth: held-out half)
    let gi: HashMap<&str, usize> = d.genes.iter().enumerate().map(|(i, x)| (x.as_str(), i)).collect();
    let sets: Vec<(String, Vec<usize>)> = std::fs::read_to_string(req(a, "--genesets")?)?.lines().filter_map(|l| { let x: Vec<&str> = l.split('\t').collect(); let m: Vec<usize> = x.iter().skip(1).filter_map(|n| gi.get(n.trim()).copied()).collect(); if m.len() >= 10 { Some((x[0].to_string(), m)) } else { None } }).collect();
    let ns = sets.len();
    let scores = |m: &[f32]| -> Vec<f64> { m.par_chunks(g).flat_map_iter(|r| { let mu = r.iter().map(|&v| v as f64).sum::<f64>() / g as f64; sets.iter().map(move |(_, idx)| idx.iter().map(|&j| r[j] as f64).sum::<f64>() / idx.len() as f64 - mu).collect::<Vec<_>>() }).collect() };
    let truth = scores(&z2);
    // responding knockdowns: halves correlate at 0.2 or more over all genes
    let resp: Vec<usize> = (0..nt).filter(|&i| pearson(&z1[i * g..(i + 1) * g].iter().map(|&v| v as f64).collect::<Vec<_>>(), &z2[i * g..(i + 1) * g].iter().map(|&v| v as f64).collect::<Vec<_>>()) >= 0.2).collect();
    println!("{ns} gene sets with 10 or more genes measured; {} responding knockdowns", resp.len());
    for (name, m) in [("whole transcriptome, half", &z1), ("projected on 160 programs", &proj), ("hybrid (panel + imputed)", &hyb), ("hybrid + guide-target prior", &hybp)] {
        let sc = scores(m);
        let per_set: Vec<f64> = (0..ns).map(|q| pearson(&(0..nt).map(|i| sc[i * ns + q]).collect::<Vec<_>>(), &(0..nt).map(|i| truth[i * ns + q]).collect::<Vec<_>>())).collect();
        rec!("pathway", "all", name, "median r per gene set", ns, median(per_set.clone()));
        let rs: Vec<f64> = (0..ns).map(|q| pearson(&resp.iter().map(|&i| sc[i * ns + q]).collect::<Vec<_>>(), &resp.iter().map(|&i| truth[i * ns + q]).collect::<Vec<_>>())).collect();
        rec!("pathway", "all", name, "median r, responders", ns, median(rs));
        // strongest gene set of each responding knockdown in the held-out half (score scaled by the square root of the set size)
        let (mut top1, mut top3, mut sign) = (0usize, 0usize, 0usize);
        for &i in &resp {
            let zs = |v: &[f64]| -> Vec<f64> { (0..ns).map(|q| v[i * ns + q] * (sets[q].1.len() as f64).sqrt()).collect() };
            let (t, r) = (zs(&truth), zs(&sc));
            let best = (0..ns).max_by(|&x, &y| t[x].abs().partial_cmp(&t[y].abs()).unwrap()).unwrap();
            let mut o: Vec<usize> = (0..ns).collect();
            o.sort_by(|&x, &y| r[y].abs().partial_cmp(&r[x].abs()).unwrap());
            if o[0] == best { top1 += 1; } if o[..3].contains(&best) { top3 += 1; } if r[best].signum() == t[best].signum() { sign += 1; }
        }
        rec!("pathway", "all", name, "top set identical", resp.len(), top1 as f64 / resp.len() as f64);
        rec!("pathway", "all", name, "top set in first 3", resp.len(), top3 as f64 / resp.len() as f64);
        rec!("pathway", "all", name, "top set same sign", resp.len(), sign as f64 / resp.len() as f64);
    }
    // ---------------- C. strength of each knockdown: DE genes in the held-out half against DE genes of the readout
    let nde = |m: &[f32]| -> Vec<f64> { m.par_chunks(g).map(|r| r.iter().filter(|v| v.abs() > DE_Z as f32).count() as f64).collect() };
    let t = nde(&z2);
    for (name, m) in [("whole transcriptome, half", &z1), ("projected on 160 programs", &proj), ("hybrid (panel + imputed)", &hyb), ("hybrid + guide-target prior", &hybp)] {
        rec!("strength", "all", name, "Spearman of DE counts", nt, spearman(&nde(m), &t));
        let en: Vec<f64> = m.par_chunks(g).map(|r| r.iter().map(|v| (*v as f64).powi(2)).sum::<f64>()).collect();
        rec!("strength", "all", name, "Spearman of squared norm", nt, spearman(&en, &t));
    }
    Ok(())
}
