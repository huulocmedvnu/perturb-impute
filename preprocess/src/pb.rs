//! Pseudobulk means per group: mean CP10k of every gene, for all cells and
//! for two interleaved halves (even / odd member position).
//! The cache is read once, front to back, in chunks of cells: each chunk is
//! prefetched sequentially, then all groups add the cells that fall in it in
//! parallel (one group per task, so accumulators are never shared).
//! Output: <out>.f32 (3 x n_groups x n_genes, f32: all, h0, h1),
//!         <out>.groups.tsv (name, n_cells), <out>.genes.txt
use crate::cache::Cache;
use anyhow::Result;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

const CHUNK_ENTRIES: usize = 256 << 20; // ~1 GB of idx+cnt per chunk

pub fn run(dir: &str, out: &str) -> Result<()> {
    let t = Instant::now();
    let c = Cache::open(dir)?;
    let ng = c.genes.len();
    let n = c.groups.len();
    let nc = c.n_cells();
    // chunk boundaries in cell rows
    let mut bounds = vec![0usize];
    let mut start = 0;
    while start < nc {
        let (a, _) = c.entry_range(start, start);
        let mut end = start + 1;
        // grow by blocks of 4096 cells until the chunk holds CHUNK_ENTRIES entries
        while end < nc {
            let nxt = (end + 4096).min(nc);
            if c.entry_range(start, nxt).1 - a > CHUNK_ENTRIES {
                break;
            }
            end = nxt;
        }
        bounds.push(end);
        start = end;
    }
    // per group: member cells with their half, already sorted by row
    let mut acc: Vec<[Vec<f64>; 2]> = (0..n).map(|_| [vec![0f64; ng], vec![0f64; ng]]).collect();
    let mut cur: Vec<usize> = vec![0; n]; // position in each group's member list
    for w in bounds.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        let (ea, eb) = c.entry_range(lo, hi);
        c.prefetch(ea, eb);
        acc.par_iter_mut().zip(cur.par_iter_mut()).enumerate().for_each(|(g, (a, k))| {
            let cells = c.cells(g);
            while *k < cells.len() && (cells[*k] as usize) < hi {
                let cell = cells[*k] as usize;
                let h = *k & 1;
                let s = c.scale[cell] as f64;
                let (ii, vv) = c.row(cell);
                for (&i, &v) in ii.iter().zip(vv) {
                    a[h][i as usize] += v as f64 * s;
                }
                *k += 1;
            }
        });
    }
    let mut w = BufWriter::new(File::create(format!("{out}.f32"))?);
    let halves: Vec<[usize; 2]> = (0..n).map(|g| {
        let k = c.cells(g).len();
        [(k + 1) / 2, k / 2]
    }).collect();
    for part in 0..3 {
        for g in 0..n {
            let [n0, n1] = halves[g];
            let row: Vec<f32> = (0..ng)
                .map(|j| match part {
                    0 => ((acc[g][0][j] + acc[g][1][j]) / (n0 + n1).max(1) as f64) as f32,
                    1 => (acc[g][0][j] / n0.max(1) as f64) as f32,
                    _ => (acc[g][1][j] / n1.max(1) as f64) as f32,
                })
                .collect();
            w.write_all(bytemuck::cast_slice(&row))?;
        }
    }
    w.flush()?;
    let mut gw = BufWriter::new(File::create(format!("{out}.groups.tsv"))?);
    for (name, _, k) in &c.groups {
        writeln!(gw, "{name}\t{k}")?;
    }
    std::fs::write(format!("{out}.genes.txt"), c.genes.join("\n") + "\n")?;
    eprintln!("pb: {n} groups x {ng} genes, {} chunks, {:.1}s", bounds.len() - 1, t.elapsed().as_secs_f64());
    Ok(())
}
