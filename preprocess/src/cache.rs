//! Read-only view of a sparse cell cache (genes.txt, groups.tsv, members.bin u32,
//! indptr.bin u64, lib.bin f32, idx.bin u16, cnt.bin u16).
use anyhow::{Context, Result};
use memmap2::Mmap;
use std::fs::File;
use std::path::Path;

pub struct Cache {
    pub genes: Vec<String>,
    pub groups: Vec<(String, usize, usize)>,
    members: Mmap,
    indptr: Mmap,
    idx: Mmap,
    cnt: Mmap,
    pub scale: Vec<f32>,
}

fn map(p: &Path) -> Result<Mmap> {
    let f = File::open(p).with_context(|| format!("open {}", p.display()))?;
    Ok(unsafe { Mmap::map(&f)? })
}

impl Cache {
    pub fn open(dir: &str) -> Result<Self> {
        let d = Path::new(dir);
        let genes = std::fs::read_to_string(d.join("genes.txt"))?.lines().map(String::from).collect();
        let groups = std::fs::read_to_string(d.join("groups.tsv"))?
            .lines()
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f[0].to_string(), f[1].parse().unwrap(), f[2].parse().unwrap())
            })
            .collect();
        let libm = map(&d.join("lib.bin"))?;
        let lib: &[f32] = bytemuck::cast_slice(&libm);
        let scale = lib.iter().map(|&l| if l > 0.0 { 1e4 / l } else { 0.0 }).collect();
        Ok(Cache { genes, groups, members: map(&d.join("members.bin"))?, indptr: map(&d.join("indptr.bin"))?, idx: map(&d.join("idx.bin"))?, cnt: map(&d.join("cnt.bin"))?, scale })
    }
    pub fn cells(&self, g: usize) -> &[u32] {
        let m: &[u32] = bytemuck::cast_slice(&self.members);
        let (_, o, n) = self.groups[g];
        &m[o..o + n]
    }
    pub fn n_cells(&self) -> usize {
        self.indptr.len() / 8 - 1
    }
    /// byte range offsets (entries) of cells [a, b)
    pub fn entry_range(&self, a: usize, b: usize) -> (usize, usize) {
        let ip: &[u64] = bytemuck::cast_slice(&self.indptr);
        (ip[a] as usize, ip[b] as usize)
    }
    /// ask the kernel to read entries [a, b) of idx and cnt ahead, sequentially
    pub fn prefetch(&self, a: usize, b: usize) {
        use memmap2::Advice;
        let _ = self.idx.advise_range(Advice::WillNeed, a * 2, (b - a) * 2);
        let _ = self.cnt.advise_range(Advice::WillNeed, a * 2, (b - a) * 2);
    }
    #[inline]
    pub fn row(&self, c: usize) -> (&[u16], &[u16]) {
        let ip: &[u64] = bytemuck::cast_slice(&self.indptr);
        let (a, b) = (ip[c] as usize, ip[c + 1] as usize);
        let idx: &[u16] = bytemuck::cast_slice(&self.idx);
        let cnt: &[u16] = bytemuck::cast_slice(&self.cnt);
        (&idx[a..b], &cnt[a..b])
    }
}
