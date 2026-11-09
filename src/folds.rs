//! Grouped 5-fold split without complex leakage (Rust port of scripts/folds.py).
//! Groups = union-find over specific-profile links (Pearson > 0.4 after removing the top 10
//! principal components of the centred profiles) and CORUM complexes (2-50 subunits) among
//! non-hub genes. Hubs (in > 3 complexes) are never tested and are excluded from the reference
//! of any fold holding a complex partner or a linked perturbation.
use crate::*;

fn find(par: &mut [usize], mut i: usize) -> usize {
    while par[i] != i { par[i] = par[par[i]]; i = par[i]; }
    i
}

pub fn run(a: &[String], z: &Z) -> Result<()> {
    let t0 = std::time::Instant::now();
    let (p, g) = (z.perts.len(), z.genes.len());
    let all_p: Vec<usize> = (0..p).collect();
    let mut x = z.mat(0, &all_p, &(0..g).collect::<Vec<_>>());
    // centre columns, top 10 right singular vectors = top eigenvectors of the centred Gram
    let mean: Vec<f64> = (0..g).map(|j| (0..p).map(|i| x[i * g + j]).sum::<f64>() / p as f64).collect();
    let mut xc = x.clone();
    xc.par_chunks_mut(g).for_each(|r| for j in 0..g { r[j] -= mean[j]; });
    let (v, _) = la::top_eig(&la::gram(&xc, p, g), g, 10)?;
    drop(xc);
    // remove the 10 programs from the uncentred profiles, then centre and scale each row
    let proj = la::gemm(false, false, p, 10, g, &x, &v);
    let back = la::gemm(false, true, p, g, 10, &proj, &v);
    x.par_chunks_mut(g).zip(back.par_chunks(g)).for_each(|(r, b)| {
        for j in 0..g { r[j] -= b[j]; }
        let m = r.iter().sum::<f64>() / g as f64;
        r.iter_mut().for_each(|v| *v -= m);
        let n = r.iter().map(|v| v * v).sum::<f64>().sqrt() + 1e-12;
        r.iter_mut().for_each(|v| *v /= n);
    });
    // links: correlation > 0.4, computed in row blocks
    let thr: f64 = arg(a, "--link").map(|v| v.parse().unwrap()).unwrap_or(0.4);
    let blk = 2048;
    let mut links: Vec<(usize, usize)> = vec![];
    for b0 in (0..p).step_by(blk) {
        let b1 = (b0 + blk).min(p);
        let c = la::gemm(false, true, b1 - b0, p, g, &x[b0 * g..b1 * g], &x);
        let found: Vec<(usize, usize)> = (b0..b1).into_par_iter().flat_map_iter(|i| { let row = &c[(i - b0) * p..(i - b0 + 1) * p]; (0..p).filter(move |&j| j != i && row[j] > thr).map(move |j| (i, j)).collect::<Vec<_>>() }).collect();
        links.extend(found);
    }
    eprintln!("folds: {} directed links > {thr} at {:.1}s", links.len(), t0.elapsed().as_secs_f64());
    // CORUM complexes (human, 2-50 subunits present)
    let idx: HashMap<&str, usize> = z.perts.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let txt = std::fs::read_to_string(req(a, "--corum")?)?;
    let rows = tsv_quoted(&txt);
    let mut lines = rows.iter();
    let hdr: Vec<&str> = lines.next().unwrap().iter().map(|s| s.as_str()).collect();
    let (co, cs) = (hdr.iter().position(|h| *h == "organism").unwrap(), hdr.iter().position(|h| *h == "subunits_gene_name").unwrap());
    let mut cx: Vec<Vec<usize>> = vec![];
    for l in lines {
        let f: Vec<&str> = l.iter().map(|s| s.as_str()).collect();
        if f.len() <= cs || f[co] != "Human" { continue; }
        let mut c: Vec<usize> = f[cs].split(';').map(|s| s.trim()).filter_map(|s| idx.get(s).copied()).collect();
        c.sort_unstable(); c.dedup();
        if (2..=50).contains(&c.len()) { cx.push(c); }
    }
    let mut deg = vec![0usize; p];
    for c in &cx { for &i in c { deg[i] += 1; } }
    let hub: Vec<bool> = deg.iter().map(|&d| d > 3).collect();
    let mut par: Vec<usize> = (0..p).collect();
    for &(i, j) in &links { if i < j && !hub[i] && !hub[j] { let (ri, rj) = (find(&mut par, i), find(&mut par, j)); par[ri] = rj; } }
    for c in &cx {
        let c: Vec<usize> = c.iter().copied().filter(|&i| !hub[i]).collect();
        for &b in c.iter().skip(1) { let (ra, rb) = (find(&mut par, c[0]), find(&mut par, b)); par[ra] = rb; }
    }
    let root: Vec<usize> = (0..p).map(|i| find(&mut par, i)).collect();
    let mut size: HashMap<usize, usize> = HashMap::new();
    for i in 0..p { if !hub[i] { *size.entry(root[i]).or_default() += 1; } }
    let mut order: Vec<(usize, usize)> = size.iter().map(|(&r, &n)| (r, n)).collect();
    order.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let mut fold = vec![-1i64; p];
    let mut load = [0usize; 5];
    let mut fold_of: HashMap<usize, i64> = HashMap::new();
    for &(r, n) in &order { let k = (0..5).min_by_key(|&k| (load[k], k)).unwrap(); fold_of.insert(r, k as i64); load[k] += n; }
    for i in 0..p { if !hub[i] { fold[i] = fold_of[&root[i]]; } }
    let mut excl: std::collections::BTreeSet<(usize, i64)> = Default::default();
    for c in &cx {
        let fs: HashSet<i64> = c.iter().filter(|&&i| !hub[i]).map(|&i| fold[i]).collect();
        for &h in c.iter().filter(|&&i| hub[i]) { for &k in &fs { excl.insert((h, k)); } }
    }
    for &(i, j) in &links { if hub[i] && !hub[j] { excl.insert((i, fold[j])); } }
    let out = req(a, "--out")?;
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.folds.tsv"))?);
    writeln!(w, "pert\tgroup\tfold")?;
    for i in 0..p { writeln!(w, "{}\t{}\t{}", z.perts[i], root[i], fold[i])?; }
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.hub_excl.tsv"))?);
    writeln!(w, "pert\tfold")?;
    for (h, k) in &excl { writeln!(w, "{}\t{k}", z.perts[*h])?; }
    eprintln!("folds: hubs {}, groups {}, largest {:?}, fold sizes {:?}, hub exclusions {}, {:.1}s",
        hub.iter().filter(|&&h| h).count(), order.len(), order.iter().take(4).map(|x| x.1).collect::<Vec<_>>(), load, excl.len(), t0.elapsed().as_secs_f64());
    Ok(())
}

/// TSV with double-quoted fields that may contain tabs, newlines and doubled quotes
pub fn tsv_quoted(s: &str) -> Vec<Vec<String>> {
    let (mut rows, mut row, mut f) = (vec![], vec![], String::new());
    let mut q = false;
    let mut it = s.chars().peekable();
    while let Some(c) = it.next() {
        match (q, c) {
            (true, '"') => { if it.peek() == Some(&'"') { f.push('"'); it.next(); } else { q = false; } }
            (true, _) => f.push(c),
            (false, '"') if f.is_empty() => q = true,
            (false, '\t') => row.push(std::mem::take(&mut f)),
            (false, '\n') => { row.push(std::mem::take(&mut f)); rows.push(std::mem::take(&mut row)); }
            (false, '\r') => {}
            _ => f.push(c),
        }
    }
    if !f.is_empty() || !row.is_empty() { row.push(f); rows.push(row); }
    rows
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn quoted_tsv() {
        let r = tsv_quoted("a\tb\tc\n1\t\"x\ty\nz\"\t\"say \"\"hi\"\"\"\n");
        assert_eq!(r.len(), 2);
        assert_eq!(r[1], vec!["1".to_string(), "x\ty\nz".into(), "say \"hi\"".into()]);
    }
}
