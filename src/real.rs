//! Measured targeted readout: the same K562 CRISPRi cells read as whole transcriptome and
//! after hybridisation enrichment for the L1000 landmarks (Replogle 2020, GSE146194 exp7 / exp11).
//! `mtx2cache` converts a 10x MatrixMarket stream to the sparse cell cache; `run` compares
//! readouts of half 0 with the whole-transcriptome half 1; `summary` prints medians and tests.
use crate::cache::Cache;
use crate::sim::{fit, greedy, metrics, pseudobulk, zhalves};
use crate::*;
use std::io::BufRead;

fn lines(p: &str) -> Result<Vec<String>> {
    Ok(std::fs::read_to_string(p).map_err(|e| anyhow!("{p}: {e}"))?.lines().map(String::from).collect())
}

/// stdin: MatrixMarket coordinate (gene, cell, count), sorted by cell. Keeps cells with a single
/// confident guide call whose barcode is also in --other-barcodes. Group = Ensembl id of the
/// target (control for Non-Targeting); members sorted by barcode.
pub fn mtx2cache(a: &[String]) -> Result<()> {
    let bc = lines(&req(a, "--barcodes")?)?;
    let other: HashSet<String> = lines(&req(a, "--other-barcodes")?)?.into_iter().collect();
    let feats: Vec<String> = lines(&req(a, "--features")?)?.iter().map(|l| l.split('\t').next().unwrap().to_string()).collect();
    if feats.len() > u16::MAX as usize + 1 { bail!("too many genes for u16 index") }
    let out = req(a, "--out")?;
    std::fs::create_dir_all(&out)?;
    let id = lines(&req(a, "--ident")?)?;
    let hdr: Vec<&str> = id[0].split(',').collect();
    let col = |n: &str| hdr.iter().position(|h| *h == n).ok_or_else(|| anyhow!("no column {n}"));
    let (cb, cg, ct, cgood, cn) = (col("cell_barcode")?, col("gene")?, col("target")?, col("good_coverage")?, col("number_of_cells")?);
    let mut grp_of: HashMap<&str, String> = HashMap::new();
    for l in &id[1..] {
        let f: Vec<&str> = l.split(',').collect();
        if f.len() != hdr.len() { bail!("unexpected comma in identities: {l}") }
        if f[cgood] != "True" || f[cn] != "1" || !other.contains(f[cb]) { continue; }
        let name = if f[cg] == "Non-Targeting" { "control".to_string() } else if f[ct].starts_with("ENSG") { f[ct].to_string() } else { f[cg].to_string() };
        grp_of.insert(f[cb], name);
    }
    // cell ids in column order
    let mut cid = vec![u32::MAX; bc.len()];
    let mut n = 0u32;
    let mut by_group: HashMap<&str, Vec<(&str, u32)>> = HashMap::new();
    for (c, b) in bc.iter().enumerate() {
        if let Some(gn) = grp_of.get(b.as_str()) { cid[c] = n; by_group.entry(gn.as_str()).or_default().push((b.as_str(), n)); n += 1; }
    }
    let ncell = n as usize;
    let (mut indptr, mut lib): (Vec<u64>, Vec<f32>) = (vec![0], vec![0.0; ncell]);
    let (mut idx, mut cnt): (Vec<u16>, Vec<u16>) = (vec![], vec![]);
    let stdin = std::io::stdin();
    let mut rd = std::io::BufReader::with_capacity(1 << 22, stdin.lock());
    let mut line = Vec::with_capacity(64);
    let mut dims_seen = false;
    let (mut cur, mut next_id) = (usize::MAX, 0usize);
    let (mut bi, mut bc_, mut bl): (Vec<u16>, Vec<u16>, f64) = (vec![], vec![], 0.0);
    let mut flush = |cur: usize, bi: &mut Vec<u16>, bc_: &mut Vec<u16>, bl: &mut f64, next_id: &mut usize, indptr: &mut Vec<u64>| -> Result<()> {
        if cur != usize::MAX && cid[cur] != u32::MAX {
            let k = cid[cur] as usize;
            if k < *next_id { bail!("matrix not sorted by cell") }
            while *next_id < k { indptr.push(idx.len() as u64); *next_id += 1; }
            idx.extend_from_slice(bi); cnt.extend_from_slice(bc_); lib[k] = *bl as f32;
            indptr.push(idx.len() as u64); *next_id += 1;
        }
        bi.clear(); bc_.clear(); *bl = 0.0;
        Ok(())
    };
    loop {
        line.clear();
        if rd.read_until(b'\n', &mut line)? == 0 { break; }
        if line[0] == b'%' { continue; }
        let mut it = line.split(|&b| b == b' ' || b == b'\n').filter(|s| !s.is_empty()).map(|s| s.iter().fold(0u64, |v, &d| v * 10 + (d - b'0') as u64));
        let (r, c, v) = (it.next().unwrap(), it.next().unwrap(), it.next().unwrap());
        if !dims_seen {
            dims_seen = true;
            if r as usize != feats.len() || c as usize != bc.len() { bail!("matrix {r} x {c}, features {} barcodes {}", feats.len(), bc.len()) }
            continue;
        }
        let c = c as usize - 1;
        if c != cur { flush(cur, &mut bi, &mut bc_, &mut bl, &mut next_id, &mut indptr)?; cur = c; }
        bi.push((r - 1) as u16); bc_.push(v.min(u16::MAX as u64) as u16); bl += v as f64;
    }
    if !dims_seen { bail!("empty matrix stream") }
    flush(cur, &mut bi, &mut bc_, &mut bl, &mut next_id, &mut indptr)?;
    drop(flush);
    while next_id < ncell { indptr.push(idx.len() as u64); next_id += 1; }
    let mut names: Vec<&str> = by_group.keys().copied().collect();
    names.sort();
    let (mut members, mut gl): (Vec<u32>, String) = (vec![], String::new());
    for gn in names {
        let v = by_group.get_mut(gn).unwrap();
        v.sort();
        gl.push_str(&format!("{gn}\t{}\t{}\n", members.len(), v.len()));
        members.extend(v.iter().map(|x| x.1));
    }
    std::fs::write(format!("{out}/genes.txt"), feats.join("\n") + "\n")?;
    std::fs::write(format!("{out}/groups.tsv"), gl)?;
    std::fs::write(format!("{out}/members.bin"), bytemuck::cast_slice::<u32, u8>(&members))?;
    std::fs::write(format!("{out}/indptr.bin"), bytemuck::cast_slice::<u64, u8>(&indptr))?;
    std::fs::write(format!("{out}/lib.bin"), bytemuck::cast_slice::<f32, u8>(&lib))?;
    std::fs::write(format!("{out}/idx.bin"), bytemuck::cast_slice::<u16, u8>(&idx))?;
    std::fs::write(format!("{out}/cnt.bin"), bytemuck::cast_slice::<u16, u8>(&cnt))?;
    eprintln!("cells {ncell} groups {} nnz {} mean UMIs {:.0}", by_group.len(), idx.len(), lib.iter().map(|&x| x as f64).sum::<f64>() / ncell as f64);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn run(a: &[String], z: &Z, fold: &[i64], group: &[i64], excl: &HashSet<(usize, i64)>, l1000: &HashSet<String>, out: &str) -> Result<()> {
    let t0 = std::time::Instant::now();
    let cw = Cache::open(&req(a, "--cache-wt")?)?;
    let ct = Cache::open(&req(a, "--cache-tg")?)?;
    if cw.groups.len() != ct.groups.len() || cw.groups.iter().zip(&ct.groups).any(|(x, y)| x.0 != y.0 || x.2 != y.2) { bail!("caches differ in groups") }
    let mut ens: HashMap<String, String> = HashMap::new();
    for l in lines(&req(a, "--gene-ids")?)? {
        let f: Vec<&str> = l.split('\t').collect();
        ens.entry(f[0].to_string()).or_insert(f[1].to_string());
    }
    let g = z.genes.len();
    let index = |c: &Cache| -> HashMap<String, usize> { c.genes.iter().enumerate().map(|(i, s)| (s.clone(), i)).collect() };
    let (iw, it) = (index(&cw), index(&ct));
    let gw: Vec<Option<usize>> = z.genes.iter().map(|s| ens.get(s).and_then(|e| iw.get(e).copied())).collect();
    let gt: Vec<Option<usize>> = z.genes.iter().map(|s| ens.get(s).and_then(|e| it.get(e).copied())).collect();
    let present: Vec<bool> = (0..g).map(|j| gw[j].is_some() && gt[j].is_some()).collect();
    eprintln!("universe genes {g}, matched by Ensembl id in both libraries {}", present.iter().filter(|&&p| p).count());
    // perturbation index of each cache group
    let pens: HashMap<&str, usize> = z.perts.iter().enumerate().filter_map(|(i, p)| ens.get(p).map(|e| (e.as_str(), i))).collect();
    let ctrl = cw.group("control").ok_or_else(|| anyhow!("no control"))?;
    let allq: Vec<usize> = std::iter::once(ctrl).chain((0..cw.groups.len()).filter(|&i| i != ctrl && cw.groups[i].2 >= 4)).collect();
    let pidx: Vec<Option<usize>> = allq.iter().map(|&gid| pens.get(cw.groups[gid].0.as_str()).copied()).collect();
    let (ncw, nct) = (cw.genes.len(), ct.genes.len());
    let mut colmap_w = vec![-1i32; ncw];
    for j in 0..g { if present[j] { colmap_w[gw[j].unwrap()] = j as i32; } }
    // L1000 panel present in both libraries
    let lm: Vec<usize> = (0..g).filter(|&j| present[j] && l1000.contains(&z.genes[j])).collect();
    let k = lm.len();
    // UMI budgets per cell
    let (mut lw, mut lt, mut pw, mut pt, mut nc) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut is_w, mut is_t) = (vec![false; ncw], vec![false; nct]);
    for &j in &lm { is_w[gw[j].unwrap()] = true; is_t[gt[j].unwrap()] = true; }
    for &gid in &allq {
        for &cell in cw.cells(gid) {
            let (i, n) = cw.row(cell as usize);
            for (&i, &n) in i.iter().zip(n) { lw += n as f64; if is_w[i as usize] { pw += n as f64; } }
            nc += 1.0;
        }
        for &cell in ct.cells(gid) {
            let (i, n) = ct.row(cell as usize);
            for (&i, &n) in i.iter().zip(n) { lt += n as f64; if is_t[i as usize] { pt += n as f64; } }
        }
    }
    let (lw, lt, pw, pt) = (lw / nc, lt / nc, pw / nc, pt / nc);
    let q_eq = lt / lw;
    let rate4 = pt / pw;
    eprintln!("cells {nc} groups {} L1000 genes {k}", allq.len() - 1);
    eprintln!("UMIs per cell: whole transcriptome {lw:.0} (panel genes {pw:.0}, {:.3}); enriched {lt:.0} (panel genes {pt:.0}, {:.3}); q_eq {q_eq:.4} panel UMI ratio {rate4:.4}", pw / lw, pt / lt);
    // readouts
    let sums = pseudobulk(&cw, &allq, &colmap_w, g, &vec![1.0; ncw], true, 1.0, 0);
    let (tz0, tz1) = zhalves(&cw, &allq, &sums, g);
    let targeted = |c: &Cache, gmap: &[Option<usize>], rate: f64, seed: u64| -> Vec<f64> {
        let nc_ = c.genes.len();
        let (mut rt, mut cm) = (vec![0.0; nc_], vec![-1i32; nc_]);
        for (kk, &j) in lm.iter().enumerate() { rt[gmap[j].unwrap()] = rate; cm[gmap[j].unwrap()] = kk as i32; }
        let s = pseudobulk(c, &allq, &cm, k, &rt, false, 1.0, seed);
        zhalves(c, &allq, &s, k).0
    };
    let z3 = targeted(&ct, &gt, 1.0, 0);
    let z4 = if rate4 < 1.0 { Some(targeted(&cw, &gw, rate4, 4401)) } else { None };
    let mut qs: Vec<(String, f64)> = [0.05, 0.1, 0.25, 0.5].iter().map(|&q| (format!("direct_q{q}"), q)).collect();
    if q_eq < 1.0 { qs.push(("direct_qeq".into(), q_eq)); }
    let zq: Vec<(String, Vec<f64>)> = qs.iter().enumerate().map(|(i, (nm, q))| {
        let s = pseudobulk(&cw, &allq, &colmap_w, g, &vec![*q; ncw], true, *q, 5500 + i as u64);
        (nm.clone(), zhalves(&cw, &allq, &s, g).0)
    }).collect();
    eprintln!("readouts at {:.1}s", t0.elapsed().as_secs_f64());
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.real.tsv"))?);
    writeln!(w, "fold\treadout\tpert\tn_cells\tn_de\tr\trp\ttp\tdisc")?;
    let inner_split = |refs: &[usize]| -> (Vec<usize>, Vec<usize>) {
        let key = |i: usize| -> u64 { let k = if group[i] >= 0 { group[i] as u64 } else { 1_000_000 + i as u64 }; (k.wrapping_mul(0x9E3779B97F4A7C15) >> 32) % 10 };
        (refs.iter().copied().filter(|&i| key(i) != 0).collect(), refs.iter().copied().filter(|&i| key(i) == 0).collect())
    };
    for f in 0..5i64 {
        let slots: Vec<usize> = (1..allq.len()).filter(|&s| pidx[s].map(|p| fold[p] == f).unwrap_or(false)).collect();
        if slots.is_empty() { continue; }
        let refs: Vec<usize> = (0..z.perts.len()).filter(|&i| fold[i] != f && !excl.contains(&(i, f))).collect();
        let (itr, iho) = inner_split(&refs);
        let mut gr = grams_c(z, &refs, &itr);
        let vk = la::top_eig_sub(&gr.c, g, 160, 40)?.0;
        gr.v = (0..g).flat_map(|i| vk[i * 160..i * 160 + NPC].to_vec()).collect();
        let denoise = |full: &[f64], n: usize| -> Vec<f64> { la::gemm(false, true, n, g, 160, &la::gemm(false, false, n, 160, g, full, &vk), &vk) };
        let gall: Vec<usize> = greedy(&gr, refs.len(), g, k, &vec![1.0; g]);
        let gp: Vec<usize> = gall[..500].iter().copied().filter(|&j| present[j]).collect();
        let gk: Vec<usize> = gall.iter().copied().filter(|&j| present[j]).collect();
        let ps = [Panel { name: "l1000".into(), genes: lm.clone() }, Panel { name: "greedy500".into(), genes: gp }, Panel { name: "greedy_lm".into(), genes: gk }];
        let used: HashSet<usize> = ps.iter().flat_map(|p| p.genes.iter().copied()).collect();
        let u: Vec<usize> = (0..g).filter(|j| present[*j] && !used.contains(j)).collect();
        let nu = u.len();
        let basis = basis_on(&gr.v, NPC, &u);
        let fits: Vec<(f64, Vec<f64>)> = ps.iter().map(|p| fit(z, &gr, &iho, p, &u)).collect::<Result<_>>()?;
        let nq = slots.len();
        let pick = |zz: &[f64], ncol: usize, cols: &[usize]| -> Vec<f64> { slots.iter().flat_map(|&s| cols.iter().map(move |&j| (s, j))).map(|(s, j)| zz[(s - 1) * ncol + j]).collect() };
        let truth = pick(&tz1, g, &u);
        let all_k: Vec<usize> = (0..k).collect();
        let mut reads: Vec<(String, Vec<f64>)> = vec![("R1_direct_full".into(), pick(&tz0, g, &u))];
        let all_g: Vec<usize> = (0..g).collect();
        { let d = denoise(&pick(&tz0, g, &all_g), nq); reads.push(("D160_denoised".into(), (0..nq).flat_map(|i| u.iter().map(move |&j| (i, j))).map(|(i, j)| d[i * g + j]).collect())); }
        reads.push(("R2_l1000_extract".into(), la::gemm(false, true, nq, nu, k, &pick(&tz0, g, &lm), &fits[0].1)));
        reads.push(("R3_l1000_enriched".into(), la::gemm(false, true, nq, nu, k, &pick(&z3, k, &all_k), &fits[0].1)));
        if let Some(z4) = &z4 { reads.push(("R4_l1000_thinned".into(), la::gemm(false, true, nq, nu, k, &pick(z4, k, &all_k), &fits[0].1))); }
        for (nm, zz) in &zq { reads.push((if nm == "direct_qeq" { "R5_direct_qeq".into() } else { nm.clone() }, pick(zz, g, &u))); }
        reads.push(("R6_greedy500_extract".into(), la::gemm(false, true, nq, nu, ps[1].genes.len(), &pick(&tz0, g, &ps[1].genes), &fits[1].1)));
        reads.push(("R7_greedyLM_extract".into(), la::gemm(false, true, nq, nu, ps[2].genes.len(), &pick(&tz0, g, &ps[2].genes), &fits[2].1)));
        // panel genes themselves: input half against the whole-transcriptome held-out half
        let (pt, pe, pn) = (pick(&tz1, g, &lm), pick(&tz0, g, &lm), pick(&z3, k, &all_k));
        for (nm, v) in [("P_extract", &pe), ("P_enriched", &pn)] {
            for (i, &s) in slots.iter().enumerate() {
                let m = metrics(&v[i * k..(i + 1) * k], &pt[i * k..(i + 1) * k], &[]);
                writeln!(w, "{f}\t{nm}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}", z.perts[pidx[s].unwrap()], cw.groups[allq[s]].2, m.2, m.0, m.1, m.3, m.4)?;
            }
        }
        // the same targets inside the atlas (input half -> held-out half of the atlas cells)
        {
            let rows: Vec<usize> = slots.iter().map(|&s| pidx[s].unwrap()).collect();
            let at = z.mat(2, &rows, &u);
            let mut ar: Vec<(&str, Vec<f64>)> = vec![("A1_direct_atlas", z.mat(1, &rows, &u))];
            ar.push(("A2_l1000_atlas", impute(z, &rows, &ps[0].genes, &fits[0].1, nu, 1)));
            ar.push(("A6_greedy500_atlas", impute(z, &rows, &ps[1].genes, &fits[1].1, nu, 1)));
            ar.push(("A7_greedyLM_atlas", impute(z, &rows, &ps[2].genes, &fits[2].1, nu, 1)));
            for (nm, v) in &ar {
                for (i, &s) in slots.iter().enumerate() {
                    let m = metrics(&v[i * nu..(i + 1) * nu], &at[i * nu..(i + 1) * nu], &basis);
                    writeln!(w, "{f}\t{nm}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}", z.perts[pidx[s].unwrap()], cw.groups[allq[s]].2, m.2, m.0, m.1, m.3, m.4)?;
                }
            }
        }
        for (nm, v) in &reads {
            for (i, &s) in slots.iter().enumerate() {
                let m = metrics(&v[i * nu..(i + 1) * nu], &truth[i * nu..(i + 1) * nu], &basis);
                writeln!(w, "{f}\t{nm}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}", z.perts[pidx[s].unwrap()], cw.groups[allq[s]].2, m.2, m.0, m.1, m.3, m.4)?;
            }
        }
        eprintln!("fold {f}: targets {nq} refs {} eval genes {nu} greedy present {} alpha {} / {} at {:.1}s", refs.len(), ps[1].genes.len(), fits[0].0, fits[1].0, t0.elapsed().as_secs_f64());
    }
    Ok(())
}

pub fn summary(a: &[String]) -> Result<()> {
    let txt = std::fs::read_to_string(format!("{}.real.tsv", req(a, "--in")?))?;
    let mut by: HashMap<String, HashMap<String, (f64, f64, f64, f64, f64)>> = HashMap::new();
    let mut order: Vec<String> = vec![];
    for l in txt.lines().skip(1) {
        let f: Vec<&str> = l.split('\t').collect();
        if !order.contains(&f[1].to_string()) { order.push(f[1].to_string()); }
        by.entry(f[1].to_string()).or_default().insert(f[2].to_string(), (f[6].parse()?, f[4].parse()?, f[7].parse()?, f[8].parse()?, f[5].parse()?));
    }
    println!("readout\tn\tmedian_r\tmedian_rp\tde_recall\topposite");
    for r in &order {
        let m = &by[r];
        let (nde, tp, disc) = m.values().fold((0.0, 0.0, 0.0), |s, v| (s.0 + v.1, s.1 + v.2, s.2 + v.3));
        println!("{r}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}", m.len(), median(m.values().map(|v| v.4).collect()), median(m.values().map(|v| v.0).collect()), tp / nde, disc / nde);
    }
    println!("\npaired\tn\tshare_a_higher\tmedian_diff\twilcoxon_p");
    for (x, y) in [("R3_l1000_enriched", "R2_l1000_extract"), ("R3_l1000_enriched", "R4_l1000_thinned"), ("R3_l1000_enriched", "R5_direct_qeq"), ("R3_l1000_enriched", "R1_direct_full"), ("R3_l1000_enriched", "D160_denoised"), ("R2_l1000_extract", "D160_denoised"), ("R6_greedy500_extract", "R2_l1000_extract"), ("R7_greedyLM_extract", "R2_l1000_extract"), ("P_enriched", "P_extract"), ("A7_greedyLM_atlas", "A2_l1000_atlas"), ("A6_greedy500_atlas", "A2_l1000_atlas")] {
        let (Some(mx), Some(my)) = (by.get(x), by.get(y)) else { continue };
        let d: Vec<f64> = mx.iter().filter_map(|(p, v)| my.get(p).map(|u| v.0 - u.0)).collect();
        println!("{x} vs {y}\t{}\t{:.3}\t{:.4}\t{:.3e}", d.len(), d.iter().filter(|v| **v > 0.0).count() as f64 / d.len() as f64, median(d.clone()), upg::wilcoxon(&d));
    }
    Ok(())
}

const MED_CHI2_1: f64 = 0.454936423119572;

/// z of half 0 restricted to its first m cells (all if None); control mean and per-cell null
/// variance from all cells, as in `zhalves`. `panel_norm`: cells normalised by the UMIs of the
/// mapped genes (targeted readout) instead of the whole library.
pub(crate) fn zsub(c: &Cache, allq: &[usize], map: &[i32], ncol: usize, panel_norm: bool, m: Option<usize>, rep: usize) -> Vec<f64> {
    let rate: Vec<f64> = map.iter().map(|&k| if panel_norm { if k >= 0 { 1.0 } else { 0.0 } } else { 1.0 }).collect();
    let sums = pseudobulk(c, allq, map, ncol, &rate, !panel_norm, 1.0, 0);
    let n = |s: usize| c.groups[allq[s]].2;
    let nq = allq.len() - 1;
    let nc = n(0) as f64;
    let mc: Vec<f64> = (0..ncol).map(|j| (sums[j] + sums[ncol + j]) / nc).collect();
    let mean = |s: usize, h: usize, j: usize| sums[(s * 2 + h) * ncol + j] / [(n(s) + 1) / 2, n(s) / 2][h] as f64;
    let s2: Vec<f64> = (0..ncol).into_par_iter().map(|j| {
        let v: Vec<f64> = (1..=nq).map(|s| { let (n0, n1) = (((n(s) + 1) / 2) as f64, (n(s) / 2) as f64); (mean(s, 0, j) - mean(s, 1, j)).powi(2) / (1.0 / n0 + 1.0 / n1) }).collect();
        median(v) / MED_CHI2_1
    }).collect();
    let mut z = vec![0.0; nq * ncol];
    z.par_chunks_mut(ncol).enumerate().for_each(|(i, zr)| {
        let cells: Vec<u32> = c.cells(allq[i + 1]).iter().step_by(2).copied().collect();
        let use_n = m.map(|m| m.min(cells.len())).unwrap_or(cells.len());
        let start = (rep * use_n).min(cells.len() - use_n);
        let mut acc = vec![0.0; ncol];
        for &cell in &cells[start..start + use_n] {
            let (idx, cnt) = c.row(cell as usize);
            let lib: f64 = if panel_norm { idx.iter().zip(cnt).filter(|(i, _)| map[**i as usize] >= 0).map(|(_, &n)| n as f64).sum() } else if c.scale[cell as usize] > 0.0 { 1e4 / c.scale[cell as usize] as f64 } else { 0.0 };
            if lib <= 0.0 { continue; }
            for (&i, &n) in idx.iter().zip(cnt) { let k = map[i as usize]; if k >= 0 { acc[k as usize] += n as f64 * 1e4 / lib; } }
        }
        for j in 0..ncol { if s2[j] > 0.0 { zr[j] = (acc[j] / use_n as f64 - mc[j]) / (s2[j] * (1.0 / use_n as f64 + 1.0 / nc)).sqrt(); } }
    });
    z
}

/// step 14b: readouts by number of input cells
#[allow(clippy::too_many_arguments)]
pub fn cells(a: &[String], z: &Z, fold: &[i64], group: &[i64], excl: &HashSet<(usize, i64)>, l1000: &HashSet<String>, out: &str) -> Result<()> {
    let cw = Cache::open(&req(a, "--cache-wt")?)?;
    let ct = Cache::open(&req(a, "--cache-tg")?)?;
    let mut ens: HashMap<String, String> = HashMap::new();
    for l in lines(&req(a, "--gene-ids")?)? { let f: Vec<&str> = l.split('\t').collect(); ens.entry(f[0].to_string()).or_insert(f[1].to_string()); }
    let g = z.genes.len();
    let index = |c: &Cache| -> HashMap<String, usize> { c.genes.iter().enumerate().map(|(i, s)| (s.clone(), i)).collect() };
    let (iw, it) = (index(&cw), index(&ct));
    let gw: Vec<Option<usize>> = z.genes.iter().map(|s| ens.get(s).and_then(|e| iw.get(e).copied())).collect();
    let gt: Vec<Option<usize>> = z.genes.iter().map(|s| ens.get(s).and_then(|e| it.get(e).copied())).collect();
    let present: Vec<bool> = (0..g).map(|j| gw[j].is_some() && gt[j].is_some()).collect();
    let pens: HashMap<&str, usize> = z.perts.iter().enumerate().filter_map(|(i, p)| ens.get(p).map(|e| (e.as_str(), i))).collect();
    let ctrl = cw.group("control").ok_or_else(|| anyhow!("no control"))?;
    let allq: Vec<usize> = std::iter::once(ctrl).chain((0..cw.groups.len()).filter(|&i| i != ctrl && cw.groups[i].2 >= 4)).collect();
    let pidx: Vec<Option<usize>> = allq.iter().map(|&gid| pens.get(cw.groups[gid].0.as_str()).copied()).collect();
    let mut map_w = vec![-1i32; cw.genes.len()];
    for j in 0..g { if present[j] { map_w[gw[j].unwrap()] = j as i32; } }
    let lm: Vec<usize> = (0..g).filter(|&j| present[j] && l1000.contains(&z.genes[j])).collect();
    let k = lm.len();
    let mut map_t = vec![-1i32; ct.genes.len()];
    for (kk, &j) in lm.iter().enumerate() { map_t[gt[j].unwrap()] = kk as i32; }
    // random panel of the same size among matched genes
    let mut rnd: Vec<usize> = (0..g).filter(|&j| present[j]).collect();
    let mut st = 0x5EEDu64;
    for i in (1..rnd.len()).rev() { let x = sm64(&mut st); rnd.swap(i, (x % (i as u64 + 1)) as usize); }
    rnd.truncate(k);
    let truth_all = { let s = pseudobulk(&cw, &allq, &map_w, g, &vec![1.0; cw.genes.len()], true, 1.0, 0); zhalves(&cw, &allq, &s, g).1 };
    let ms: Vec<Option<usize>> = vec![Some(25), Some(50), Some(100), Some(200), None];
    let mr: Vec<(Option<usize>, usize)> = ms.iter().flat_map(|&m| (0..if m.is_some() { 8 } else { 1 }).map(move |r| (m, r))).collect();
    let zs: Vec<(Vec<f64>, Vec<f64>)> = mr.iter().map(|&(m, r)| (zsub(&cw, &allq, &map_w, g, false, m, r), zsub(&ct, &allq, &map_t, k, true, m, r))).collect();
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.cells.tsv"))?);
    writeln!(w, "fold\treadout\tm\trep\tpert\tn_cells\tn_de\tr\trp\ttp\tdisc")?;
    let inner_split = |refs: &[usize]| -> (Vec<usize>, Vec<usize>) {
        let key = |i: usize| -> u64 { let k = if group[i] >= 0 { group[i] as u64 } else { 1_000_000 + i as u64 }; (k.wrapping_mul(0x9E3779B97F4A7C15) >> 32) % 10 };
        (refs.iter().copied().filter(|&i| key(i) != 0).collect(), refs.iter().copied().filter(|&i| key(i) == 0).collect())
    };
    for f in 0..5i64 {
        let slots: Vec<usize> = (1..allq.len()).filter(|&s| pidx[s].map(|p| fold[p] == f).unwrap_or(false)).collect();
        if slots.is_empty() { continue; }
        let refs: Vec<usize> = (0..z.perts.len()).filter(|&i| fold[i] != f && !excl.contains(&(i, f))).collect();
        let (itr, iho) = inner_split(&refs);
        let mut gr = grams_c(z, &refs, &itr);
        let vk = la::top_eig_sub(&gr.c, g, 160, 40)?.0;
        gr.v = (0..g).flat_map(|i| vk[i * 160..i * 160 + NPC].to_vec()).collect();
        let denoise = |full: &[f64], n: usize| -> Vec<f64> { la::gemm(false, true, n, g, 160, &la::gemm(false, false, n, 160, g, full, &vk), &vk) };
        let gk: Vec<usize> = greedy(&gr, refs.len(), g, k, &vec![1.0; g]).into_iter().filter(|&j| present[j]).collect();
        let ps = [Panel { name: "l1000".into(), genes: lm.clone() }, Panel { name: "greedy".into(), genes: gk }, Panel { name: "random".into(), genes: rnd.clone() }];
        let used: HashSet<usize> = ps.iter().flat_map(|p| p.genes.iter().copied()).collect();
        let u: Vec<usize> = (0..g).filter(|j| present[*j] && !used.contains(j)).collect();
        let nu = u.len();
        let basis = basis_on(&gr.v, NPC, &u);
        let fits: Vec<(f64, Vec<f64>)> = ps.iter().map(|p| fit(z, &gr, &iho, p, &u)).collect::<Result<_>>()?;
        let nq = slots.len();
        let pick = |zz: &[f64], ncol: usize, cols: &[usize]| -> Vec<f64> { slots.iter().flat_map(|&s| cols.iter().map(move |&j| (s, j))).map(|(s, j)| zz[(s - 1) * ncol + j]).collect() };
        let truth = pick(&truth_all, g, &u);
        let all_k: Vec<usize> = (0..k).collect();
        let all_g: Vec<usize> = (0..g).collect();
        for (mi, &(m, rep)) in mr.iter().enumerate() {
            let (zw, zt) = &zs[mi];
            let mut reads: Vec<(&str, Vec<f64>)> = vec![("direct", pick(zw, g, &u))];
            { let d = denoise(&pick(zw, g, &all_g), nq); reads.push(("denoised", (0..nq).flat_map(|i| u.iter().map(move |&j| (i, j))).map(|(i, j)| d[i * g + j]).collect())); }
            for (pk, nm) in [(0, "l1000_extract"), (1, "greedy_extract"), (2, "random_extract")] {
                reads.push((nm, la::gemm(false, true, nq, nu, ps[pk].genes.len(), &pick(zw, g, &ps[pk].genes), &fits[pk].1)));
            }
            reads.push(("l1000_enriched", la::gemm(false, true, nq, nu, k, &pick(zt, k, &all_k), &fits[0].1)));
            for (nm, v) in &reads {
                for (i, &s) in slots.iter().enumerate() {
                    if let Some(mm) = m { if (cw.groups[allq[s]].2 + 1) / 2 < (rep + 1) * mm { continue; } }
                    let mt = metrics(&v[i * nu..(i + 1) * nu], &truth[i * nu..(i + 1) * nu], &basis);
                    writeln!(w, "{f}\t{nm}\t{}\t{rep}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}", m.map(|x| x.to_string()).unwrap_or("all".into()), z.perts[pidx[s].unwrap()], cw.groups[allq[s]].2, mt.2, mt.0, mt.1, mt.3, mt.4)?;
                }
            }
        }
        eprintln!("fold {f}: targets {nq} eval genes {nu}");
    }
    Ok(())
}

pub fn cellsum(a: &[String]) -> Result<()> {
    let txt = std::fs::read_to_string(format!("{}.cells.tsv", req(a, "--in")?))?;
    let mut by: HashMap<(String, String), HashMap<String, (f64, f64, f64, f64)>> = HashMap::new();
    let (mut ro, mut mo): (Vec<String>, Vec<String>) = (vec![], vec![]);
    let mut dsc: HashMap<(String, String), f64> = HashMap::new();
    for l in txt.lines().skip(1) {
        let f: Vec<&str> = l.split('\t').collect();
        if !ro.contains(&f[1].to_string()) { ro.push(f[1].to_string()); }
        if !mo.contains(&f[2].to_string()) { mo.push(f[2].to_string()); }
        let e = by.entry((f[1].to_string(), f[2].to_string())).or_default().entry(f[4].to_string()).or_insert((0.0, 0.0, 0.0, 0.0));
        let n = e.3 + 1.0;
        // running mean of specific r over replicates; DE counts summed
        *e = ((e.0 * e.3 + f[8].parse::<f64>()?) / n, e.1 + f[6].parse::<f64>()?, e.2 + f[9].parse::<f64>()?, n);
        dsc.entry((f[1].to_string(), f[2].to_string())).and_modify(|v: &mut f64| *v += f[10].parse::<f64>().unwrap()).or_insert(f[10].parse::<f64>()?);
    }
    println!("m\treadout\tn\tmedian_rp\tde_recall\topposite");
    for m in &mo { for r in &ro {
        let x = &by[&(r.clone(), m.clone())];
        let (nde, tp) = x.values().fold((0.0, 0.0), |s, v| (s.0 + v.1, s.1 + v.2));
        println!("{m}\t{r}\t{}\t{:.4}\t{:.4}\t{:.4}", x.len(), median(x.values().map(|v| v.0).collect()), tp / nde, dsc[&(r.clone(), m.clone())] / nde);
    } }
    println!("\nm\tpaired\tshare_a_higher\tmedian_diff\twilcoxon_p");
    for m in &mo {
        for (p, q) in [("greedy_extract", "l1000_extract"), ("greedy_extract", "random_extract"), ("random_extract", "l1000_extract"), ("l1000_enriched", "l1000_extract"), ("l1000_enriched", "direct"), ("greedy_extract", "direct")] {
            let (mx, my) = (&by[&(p.to_string(), m.clone())], &by[&(q.to_string(), m.clone())]);
            let d: Vec<f64> = mx.iter().filter_map(|(k, v)| my.get(k).map(|u| v.0 - u.0)).collect();
            println!("{m}\t{p} vs {q}\t{:.3}\t{:.4}\t{:.3e}", d.iter().filter(|v| **v > 0.0).count() as f64 / d.len() as f64, median(d.clone()), upg::wilcoxon(&d));
        }
    }
    Ok(())
}
