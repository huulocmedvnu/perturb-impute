//! orion-ingest: streaming pseudobulk summariser for X-Atlas/Orion parquet shards.
//! Only non-zero (gene, count) entries are touched. Per cell: cp10k scale = 1e4 / total_counts,
//! scatter-add into per-perturbation accumulators for half 0 / half 1 (alternating within each
//! GEM batch) and for the high-sgRNA cells (top quartile of guide UMIs within perturbation x batch).
//!
//!   orion-ingest schema FILE.parquet
//!   orion-ingest scan   --shard FILE --meta gene_metadata.parquet       (gene detection rates)
//!   orion-ingest run    --cell HCT116 --shards list.txt --meta gene_metadata.parquet
//!                       --guides guide_calls.csv.gz --genes genes.txt --work DIR
//!   orion-ingest finish --work DIR --out PREFIX [--shared shared_genes.txt]
use anyhow::{anyhow, bail, Context, Result};
use arrow_array::{Array, Float64Array, Int64Array, ListArray, StringArray};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::time::Instant;

const MED_CHI2_1: f64 = 0.454936423119572;
const CONTROL: &str = "Non-Targeting";

fn arg(a: &[String], k: &str) -> Option<String> { a.iter().position(|x| x == k).and_then(|i| a.get(i + 1).cloned()) }
fn req(a: &[String], k: &str) -> Result<String> { arg(a, k).ok_or_else(|| anyhow!("missing {k}")) }

/// 2-bit packed 16-bp barcode (bits 0..32) + batch id (bits 32..48); None if not ACGT x16
fn pack(bc: &str, batch: u16) -> Option<u64> {
    let s = bc.as_bytes();
    if s.len() < 16 { return None; }
    let mut v = 0u64;
    for &c in &s[..16] { v = (v << 2) | match c { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => return None }; }
    Some(v | (batch as u64) << 32)
}
fn batch_of(bc: &str) -> Option<u16> { bc.rsplit("_Batch").next().and_then(|x| x.parse().ok()) }

/// token -> (ensembl, symbol)
fn gene_meta(path: &str) -> Result<Vec<(i64, String, String)>> {
    let r = ParquetRecordBatchReaderBuilder::try_new(std::fs::File::open(path)?)?.build()?;
    let mut out = vec![];
    for b in r {
        let b = b?;
        let e = b.column(0).as_any().downcast_ref::<StringArray>().unwrap();
        let n = b.column(1).as_any().downcast_ref::<StringArray>().unwrap();
        let t = b.column(2).as_any().downcast_ref::<Int64Array>().unwrap();
        for i in 0..b.num_rows() { out.push((t.value(i), e.value(i).to_string(), n.value(i).to_string())); }
    }
    Ok(out)
}

/// one shard decoded into flat sparse arrays (cells x tokens)
struct Shard { barcode: Vec<String>, target: Vec<String>, total: Vec<f64>, pass: Vec<bool>, off: Vec<usize>, tok: Vec<i64>, val: Vec<f64> }
fn read_shard(path: &str) -> Result<Shard> {
    let b = ParquetRecordBatchReaderBuilder::try_new(std::fs::File::open(path).with_context(|| path.to_string())?)?;
    let schema = b.schema().clone();
    let idx = |n: &str| schema.index_of(n).map_err(|e| anyhow!("{e}"));
    let (ci, cv, cb, ct, cn, cp) = (idx("gene_token_id")?, idx("gene_expression")?, idx("cell_barcode")?, idx("gene_target")?, idx("total_counts")?, idx("pass_guide_filter")?);
    let mut s = Shard { barcode: vec![], target: vec![], total: vec![], pass: vec![], off: vec![0], tok: vec![], val: vec![] };
    for rb in b.with_batch_size(8192).build()? {
        let rb = rb?;
        let li = rb.column(ci).as_any().downcast_ref::<ListArray>().unwrap();
        let lv = rb.column(cv).as_any().downcast_ref::<ListArray>().unwrap();
        let bc = rb.column(cb).as_any().downcast_ref::<StringArray>().unwrap();
        let tg = rb.column(ct).as_any().downcast_ref::<StringArray>().unwrap();
        let tc = rb.column(cn).as_any().downcast_ref::<arrow_array::Float64Array>().unwrap();
        let ps = rb.column(cp).as_any().downcast_ref::<Int64Array>().unwrap();
        let it = li.values().as_any().downcast_ref::<Int64Array>().unwrap();
        let vv = lv.values().as_any().downcast_ref::<Float64Array>().unwrap();
        let (io, vo) = (li.value_offsets(), lv.value_offsets());
        for r in 0..rb.num_rows() {
            let (a, z) = (io[r] as usize, io[r + 1] as usize);
            let (a2, z2) = (vo[r] as usize, vo[r + 1] as usize);
            if z - a != z2 - a2 { bail!("token and value lists differ in length") }
            s.tok.extend_from_slice(&it.values()[a..z]);
            s.val.extend_from_slice(&vv.values()[a2..z2]);
            s.off.push(s.tok.len());
            s.barcode.push(bc.value(r).to_string());
            s.target.push(tg.value(r).to_string());
            s.total.push(tc.value(r));
            s.pass.push(ps.value(r) == 1);
        }
    }
    Ok(s)
}

fn schema(path: &str) -> Result<()> {
    let b = ParquetRecordBatchReaderBuilder::try_new(std::fs::File::open(path)?)?;
    println!("rows {}", b.metadata().file_metadata().num_rows());
    for f in b.schema().fields() { println!("{}\t{:?}", f.name(), f.data_type()); }
    Ok(())
}

/// detection: mean count per cell for each token in one shard (to define the expressed gene set)
fn scan(a: &[String]) -> Result<()> {
    let t0 = Instant::now();
    let s = read_shard(&req(a, "--shard")?)?;
    let meta = gene_meta(&req(a, "--meta")?)?;
    let maxt = meta.iter().map(|m| m.0).max().unwrap() as usize;
    let mut sum = vec![0f64; maxt + 1];
    for (t, v) in s.tok.iter().zip(&s.val) { sum[*t as usize] += v; }
    let n = s.total.len() as f64;
    let mut w = BufWriter::new(std::fs::File::create(req(a, "--out")?)?);
    writeln!(w, "token\tensembl\tsymbol\tmean_count")?;
    for (t, e, g) in &meta { writeln!(w, "{t}\t{e}\t{g}\t{:.6}", sum[*t as usize] / n)?; }
    let nnz = s.tok.len();
    eprintln!("scan: {} cells, {nnz} non-zeros ({:.0} per cell), {:.1}s", s.total.len(), nnz as f64 / n, t0.elapsed().as_secs_f64());
    for th in [0.01, 0.05, 0.1, 0.2] { eprintln!("  tokens with mean count >= {th}: {}", sum.iter().filter(|&&x| x / n >= th).count()); }
    Ok(())
}

/// accumulator state (perturbation rows x genes), saved as checkpoints
struct Acc { perts: Vec<String>, g: usize, h: [Vec<f64>; 2], hh: [Vec<f32>; 2], n: Vec<[u32; 4]>, done: Vec<String> }
impl Acc {
    fn save(&self, path: &str) -> Result<()> {
        let tmp = format!("{path}.tmp");
        let mut w = BufWriter::new(std::fs::File::create(&tmp)?);
        writeln!(w, "{}\t{}", self.perts.len(), self.g)?;
        writeln!(w, "{}", self.perts.join("\t"))?;
        writeln!(w, "{}", self.done.join("\t"))?;
        for v in &self.h { w.write_all(bytemuck::cast_slice(v))?; }
        for v in &self.hh { w.write_all(bytemuck::cast_slice(v))?; }
        for x in &self.n { w.write_all(bytemuck::cast_slice(x))?; }
        drop(w);
        std::fs::rename(&tmp, path)?;
        Ok(())
    }
    fn load(path: &str) -> Result<Acc> {
        let mut r = BufReader::new(std::fs::File::open(path)?);
        let mut l = String::new();
        r.read_line(&mut l)?;
        let d: Vec<usize> = l.trim().split('\t').map(|x| x.parse().unwrap()).collect();
        let (p, g) = (d[0], d[1]);
        l.clear(); r.read_line(&mut l)?;
        let perts = l.trim_end_matches('\n').split('\t').map(String::from).collect();
        l.clear(); r.read_line(&mut l)?;
        let done = l.trim_end_matches('\n').split('\t').filter(|s| !s.is_empty()).map(String::from).collect();
        let mut rd = |n: usize, sz: usize| -> Result<Vec<u8>> { let mut b = vec![0u8; n * sz]; r.read_exact(&mut b)?; Ok(b) };
        let h0: Vec<f64> = bytemuck::cast_slice(&rd(p * g, 8)?).to_vec();
        let h1: Vec<f64> = bytemuck::cast_slice(&rd(p * g, 8)?).to_vec();
        let q0: Vec<f32> = bytemuck::cast_slice(&rd(p * g, 4)?).to_vec();
        let q1: Vec<f32> = bytemuck::cast_slice(&rd(p * g, 4)?).to_vec();
        let nn: Vec<u32> = bytemuck::cast_slice(&rd(p * 4, 4)?).to_vec();
        Ok(Acc { perts, g, h: [h0, h1], hh: [q0, q1], n: nn.chunks(4).map(|c| [c[0], c[1], c[2], c[3]]).collect(), done })
    }
}

fn run(a: &[String]) -> Result<()> {
    let t0 = Instant::now();
    let cell = req(a, "--cell")?;
    let work = req(a, "--work")?;
    std::fs::create_dir_all(&work)?;
    let shards: Vec<String> = std::fs::read_to_string(req(a, "--shards")?)?.lines().filter(|l| !l.is_empty()).map(String::from).collect();
    // genes: symbols in order; LUT token -> compact index
    let genes: Vec<String> = std::fs::read_to_string(req(a, "--genes")?)?.lines().map(String::from).collect();
    let gi: FxHashMap<&str, u32> = genes.iter().enumerate().map(|(i, s)| (s.as_str(), i as u32)).collect();
    let meta = gene_meta(&req(a, "--meta")?)?;
    let maxt = meta.iter().map(|m| m.0).max().unwrap() as usize;
    let mut lut = vec![u32::MAX; maxt + 1];
    for (t, _, sym) in &meta { if let Some(&j) = gi.get(sym.as_str()) { if lut[*t as usize] == u32::MAX { lut[*t as usize] = j; } } }
    let g = genes.len();
    eprintln!("{cell}: {g} genes, {} tokens mapped", lut.iter().filter(|&&x| x != u32::MAX).count());
    // guide calls: packed barcode -> (target, total sgRNA UMI); high flag = top quartile within target x batch
    let tg = Instant::now();
    let mut gd: FxHashMap<u64, (u32, u32)> = FxHashMap::default();
    let mut targets: Vec<String> = vec![];
    let mut tix: FxHashMap<String, u32> = FxHashMap::default();
    let rd = BufReader::new(flate2::read::MultiGzDecoder::new(std::fs::File::open(req(a, "--guides")?)?));
    for (i, l) in rd.lines().enumerate() {
        let l = l?;
        if i == 0 { continue; }
        let f: Vec<&str> = l.split(',').collect();
        if f.len() < 6 || f[5] != "True" { continue; }
        let Some(b) = batch_of(f[0]) else { continue };
        let Some(key) = pack(f[0], b) else { continue };
        let umi: u32 = f[3].split('|').map(|x| x.parse::<u32>().unwrap_or(0)).sum();
        let t = *tix.entry(f[4].to_string()).or_insert_with(|| { targets.push(f[4].to_string()); (targets.len() - 1) as u32 });
        gd.insert(key, (t, umi));
    }
    let mut groups: FxHashMap<(u32, u16), Vec<u32>> = FxHashMap::default();
    for (k, (t, u)) in &gd { groups.entry((*t, (k >> 32) as u16)).or_default().push(*u); }
    let mut cut: FxHashMap<(u32, u16), u32> = FxHashMap::default();
    for (k, mut v) in groups { v.sort_unstable(); let q = v[((v.len() as f64) * 0.75).floor() as usize % v.len()]; cut.insert(k, q); }
    eprintln!("guide calls: {} cells, {} targets, in {:.1}s", gd.len(), targets.len(), tg.elapsed().as_secs_f64());
    // accumulators (resume from checkpoint)
    let ck = format!("{work}/{cell}.ckpt");
    let mut acc = if std::path::Path::new(&ck).exists() { let x = Acc::load(&ck)?; eprintln!("resumed: {} shards done", x.done.len()); x }
        else { let p = targets.len(); Acc { perts: targets.clone(), g, h: [vec![0.0; p * g], vec![0.0; p * g]], hh: [vec![0.0; p * g], vec![0.0; p * g]], n: vec![[0; 4]; p], done: vec![] } };
    if acc.perts != targets { bail!("checkpoint perturbation list differs") }
    let p = acc.perts.len();
    let todo: Vec<String> = shards.into_iter().filter(|s| !acc.done.contains(s)).collect();
    // prefetch thread: download shard k+1 while shard k is processed (a local file is used as is)
    let (tx, rx) = std::sync::mpsc::sync_channel::<(String, Result<String>, f64)>(1);
    let list = todo.clone();
    let fw = work.clone();
    let prod = std::thread::spawn(move || {
        let fetch = |s: &str| -> Result<String> {
            let path = format!("{fw}/{s}");
            if std::path::Path::new(&path).exists() { return Ok(path); }
            for attempt in 0..5 {
                let st = std::process::Command::new("curl").args(["-sfL", "--retry", "3", "-o", &path, &format!("https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/resolve/main/data/{s}")]).status()?;
                if st.success() { return Ok(path); }
                eprintln!("download {s} failed (attempt {attempt}), retrying");
                std::thread::sleep(std::time::Duration::from_secs(10));
            }
            bail!("download {s} failed")
        };
        for s in list { let t = Instant::now(); let r = fetch(&s); if tx.send((s, r, t.elapsed().as_secs_f64())).is_err() { break; } }
    });
    let (mut ncell, mut nnz_tot, mut t_parse, mut t_acc, mut t_dl) = (0usize, 0usize, 0.0, 0.0, 0.0);
    for (k, (s, path, dl)) in rx.iter().enumerate() {
        let path = path?;
        t_dl += dl;
        let tp = Instant::now();
        let sh = read_shard(&path)?;
        std::fs::remove_file(&path)?;
        t_parse += tp.elapsed().as_secs_f64();
        let ta = Instant::now();
        let b: u16 = s.trim_end_matches(".parquet").rsplit("Batch").next().and_then(|x| x.parse().ok()).ok_or_else(|| anyhow!("no batch number in {s}"))?;
        // cells -> (perturbation row, half, high flag, high half); alternate within perturbation in this batch
        let mut cnt = vec![[0u32; 2]; p];
        let mut cells: Vec<Vec<(usize, u8)>> = vec![vec![]; p];
        let mut miss = 0usize;
        for i in 0..sh.total.len() {
            if !sh.pass[i] || sh.total[i] <= 0.0 { continue; }
            let key = match pack(&sh.barcode[i], b) { Some(x) => x, None => { miss += 1; continue } };
            let Some(&(t, umi)) = gd.get(&key) else { miss += 1; continue };
            let t = t as usize;
            if acc.perts[t] != sh.target[i] { miss += 1; continue; }
            let half = (cnt[t][0] % 2) as u8;
            cnt[t][0] += 1;
            let mut code = half;
            if umi >= *cut.get(&(t as u32, b)).unwrap_or(&u32::MAX) { code |= (((cnt[t][1] % 2) as u8) << 1) | 4; cnt[t][1] += 1; }
            cells[t].push((i, code));
        }
        let g_ = g;
        let (h0, h1) = acc.h.split_at_mut(1);
        let (q0, q1) = acc.hh.split_at_mut(1);
        h0[0].par_chunks_mut(g_).zip(h1[0].par_chunks_mut(g_)).zip(q0[0].par_chunks_mut(g_).zip(q1[0].par_chunks_mut(g_))).zip(acc.n.par_iter_mut()).enumerate()
            .for_each(|(t, (((a0, a1), (b0, b1)), nn))| {
                for &(i, code) in &cells[t] {
                    let sc = 1e4 / sh.total[i];
                    let dst = if code & 1 == 0 { &mut *a0 } else { &mut *a1 };
                    let (lo, hi) = (sh.off[i], sh.off[i + 1]);
                    for e in lo..hi { let j = lut[sh.tok[e] as usize]; if j != u32::MAX { dst[j as usize] += sh.val[e] * sc; } }
                    nn[(code & 1) as usize] += 1;
                    if code & 4 != 0 {
                        let dq = if code & 2 == 0 { &mut *b0 } else { &mut *b1 };
                        for e in lo..hi { let j = lut[sh.tok[e] as usize]; if j != u32::MAX { dq[j as usize] += (sh.val[e] * sc) as f32; } }
                        nn[2 + ((code >> 1) & 1) as usize] += 1;
                    }
                }
            });
        t_acc += ta.elapsed().as_secs_f64();
        ncell += sh.total.len();
        nnz_tot += sh.tok.len();
        acc.done.push(s.clone());
        eprintln!("[{}/{}] {s}: {} cells ({} unmatched), {:.0} nnz/cell; download {dl:.1}s, parse {:.1}s, accumulate {:.2}s; elapsed {:.0}s",
            k + 1, todo.len(), sh.total.len(), miss, sh.tok.len() as f64 / sh.total.len() as f64, tp.elapsed().as_secs_f64() - ta.elapsed().as_secs_f64(), ta.elapsed().as_secs_f64(), t0.elapsed().as_secs_f64());
        if acc.done.len() % 10 == 0 { let tc = Instant::now(); acc.save(&ck)?; eprintln!("  checkpoint saved in {:.1}s", tc.elapsed().as_secs_f64()); }
    }
    prod.join().map_err(|_| anyhow!("download thread panicked"))?;
    acc.save(&ck)?;
    eprintln!("{cell} done: {ncell} cells, {nnz_tot} non-zeros; download {t_dl:.0}s (overlapped), parse {t_parse:.0}s, accumulate {t_acc:.0}s, wall {:.0}s", t0.elapsed().as_secs_f64());
    Ok(())
}

/// means per half, null variance from half differences, z for all / half 0 / half 1 (main and high-sgRNA)
fn finish(a: &[String]) -> Result<()> {
    let work = req(a, "--work")?;
    let cell = req(a, "--cell")?;
    let acc = Acc::load(&format!("{work}/{cell}.ckpt"))?;
    let genes: Vec<String> = std::fs::read_to_string(req(a, "--genes")?)?.lines().map(String::from).collect();
    let g = acc.g;
    if genes.len() != g { bail!("gene list does not match the checkpoint") }
    let out = req(a, "--out")?;
    let ci = acc.perts.iter().position(|x| x == CONTROL).ok_or_else(|| anyhow!("no control"))?;
    for (tag, hi) in [("", false), ("_high_sgrna", true)] {
        let n = |t: usize, h: usize| -> f64 { acc.n[t][if hi { 2 + h } else { h }] as f64 };
        let s = |t: usize, h: usize, j: usize| -> f64 { if hi { acc.hh[h][t * g + j] as f64 } else { acc.h[h][t * g + j] } };
        // control mean over all control cells (main control cells for both)
        let ncc = n(ci, 0) + n(ci, 1);
        let ncm = acc.n[ci][0] as f64 + acc.n[ci][1] as f64;
        let mc: Vec<f64> = (0..g).map(|j| (acc.h[0][ci * g + j] + acc.h[1][ci * g + j]) / ncm).collect();
        let _ = ncc;
        let perts: Vec<usize> = (0..acc.perts.len()).filter(|&t| t != ci && n(t, 0) >= 2.0 && n(t, 1) >= 2.0).collect();
        let s2: Vec<f64> = (0..g).into_par_iter().map(|j| {
            let mut v: Vec<f64> = perts.iter().map(|&t| { let (a0, a1) = (n(t, 0), n(t, 1)); (s(t, 0, j) / a0 - s(t, 1, j) / a1).powi(2) / (1.0 / a0 + 1.0 / a1) }).collect();
            v.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let m = v.len();
            (if m % 2 == 1 { v[m / 2] } else { 0.5 * (v[m / 2 - 1] + v[m / 2]) }) / MED_CHI2_1
        }).collect();
        let np = perts.len();
        let mut z = vec![0f32; 3 * np * g];
        for layer in 0..3 {
            z[layer * np * g..(layer + 1) * np * g].par_chunks_mut(g).zip(&perts).for_each(|(row, &t)| {
                let (a0, a1) = (n(t, 0), n(t, 1));
                let (m, nn) = match layer { 0 => (None, a0 + a1), 1 => (Some(0), a0), _ => (Some(1), a1) };
                for j in 0..g {
                    let mean = match m { None => (s(t, 0, j) + s(t, 1, j)) / nn, Some(h) => s(t, h, j) / nn };
                    row[j] = if s2[j] > 0.0 { ((mean - mc[j]) / (s2[j] * (1.0 / nn + 1.0 / ncm)).sqrt()) as f32 } else { 0.0 };
                }
            });
        }
        let pre = format!("{out}{tag}");
        std::fs::write(format!("{pre}.z.f32"), bytemuck::cast_slice(&z))?;
        let mut w = BufWriter::new(std::fs::File::create(format!("{pre}.meta.tsv"))?);
        writeln!(w, "pert\tn\tn_h0\tn_h1")?;
        for &t in &perts { writeln!(w, "{}\t{}\t{}\t{}", acc.perts[t], n(t, 0) + n(t, 1), n(t, 0), n(t, 1))?; }
        let mut w = BufWriter::new(std::fs::File::create(format!("{pre}.genes.tsv"))?);
        writeln!(w, "gene\ts2\tctrl_mean")?;
        for j in 0..g { writeln!(w, "{}\t{}\t{}", genes[j], s2[j], mc[j])?; }
        eprintln!("{cell}{tag}: {np} perturbations x {g} genes; control cells {ncm}; median cells per perturbation {}",
            { let mut v: Vec<f64> = perts.iter().map(|&t| n(t, 0) + n(t, 1)).collect(); v.sort_by(|x, y| x.partial_cmp(y).unwrap()); v[v.len() / 2] });
    }
    Ok(())
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    match a.get(1).map(|s| s.as_str()) {
        Some("schema") => schema(&a[2]),
        Some("scan") => scan(&a),
        Some("run") => run(&a),
        Some("finish") => finish(&a),
        _ => bail!("usage: orion-ingest schema|scan|run|finish"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn barcode_packing_is_unique_and_rejects_n() {
        let a = pack("AAACCAAAGACATGTT-HCT116_Batch1", 1).unwrap();
        let b = pack("AAACCAAAGACATGTT-HCT116_Batch2", 2).unwrap();
        let c = pack("AAACCAAAGACATGTA-HCT116_Batch1", 1).unwrap();
        assert!(a != b && a != c);
        assert_eq!(a & 0xFFFF_FFFF, b & 0xFFFF_FFFF);
        assert!(pack("AAACCAAAGACATGNT-x", 1).is_none());
        assert_eq!(batch_of("AAACCAAAGACATGTT-HCT116_Batch17"), Some(17));
    }
}
