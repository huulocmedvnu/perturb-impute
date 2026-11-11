//! Additional analysis (Additional file 1: Table S1, Note S1): does a sparse gene-gene graph
//! built inside the training fold improve on ridge imputation from a 500-gene panel?
//!
//!   panelimp causal --z PREFIX --folds folds.tsv --hubs hub_excl.tsv --out PREFIX [--fold 0] [--metal]
//!
//! Final estimate for the genes outside the panel: ridge + s * (graph correction), where the
//! correction propagates residuals (after removal of R global factors) over
//!   dir   directed edges from knockdowns whose target gene is measured (network deconvolution of
//!         per-unit total effects, top-k incoming edges per gene),
//!   cond  sparse conditional-dependency edges (regression form of the regularised inverse
//!         covariance of the residual profiles, top-k per gene),
//!   omp   a sparse combination of training knockdowns that explains the panel residual.
//! Every hyperparameter and the scale s are chosen on the inner holdout of the training fold.
//! All matrices are f32; graphs are CSR; z-scores are memory-mapped.
use crate::*;
use memmap2::Mmap;

extern "C" {
    fn cblas_sgemm(o: i32, ta: i32, tb: i32, m: i32, n: i32, k: i32, al: f32, a: *const f32, lda: i32, b: *const f32, ldb: i32, be: f32, c: *mut f32, ldc: i32);
    fn cblas_ssyrk(o: i32, up: i32, t: i32, n: i32, k: i32, al: f32, a: *const f32, lda: i32, be: f32, c: *mut f32, ldc: i32);
    fn spotrf_(uplo: *const u8, n: *const i32, a: *mut f32, lda: *const i32, info: *mut i32);
    fn spotri_(uplo: *const u8, n: *const i32, a: *mut f32, lda: *const i32, info: *mut i32);
    fn sgesv_(n: *const i32, nrhs: *const i32, a: *mut f32, lda: *const i32, ipiv: *mut i32, b: *mut f32, ldb: *const i32, info: *mut i32);
}
pub(crate) fn sgemm(ta: bool, tb: bool, m: usize, n: usize, k: usize, a: &[f32], b: &[f32]) -> Vec<f32> {
    let mut c = vec![0.0f32; m * n];
    let (lda, ldb) = (if ta { m } else { k }, if tb { k } else { n });
    unsafe { cblas_sgemm(101, if ta { 112 } else { 111 }, if tb { 112 } else { 111 }, m as i32, n as i32, k as i32, 1.0, a.as_ptr(), lda as i32, b.as_ptr(), ldb as i32, 0.0, c.as_mut_ptr(), n as i32) };
    c
}
/// A^T A (full symmetric), A row-major r x c
pub(crate) fn sgram(a: &[f32], r: usize, c: usize) -> Vec<f32> {
    let mut g = vec![0.0f32; c * c];
    unsafe { cblas_ssyrk(101, 121, 112, c as i32, r as i32, 1.0, a.as_ptr(), c as i32, 0.0, g.as_mut_ptr(), c as i32) };
    for i in 0..c { for j in 0..i { g[i * c + j] = g[j * c + i]; } }
    g
}
/// leading k eigenvectors (g x k row-major) by block subspace iteration in f32 with a Rayleigh-Ritz step
pub(crate) fn top_eig32(c: &[f32], g: usize, k: usize, iters: usize) -> Result<Vec<f32>> {
    let b = (k + 16).min(g);
    let mut st = 0x9E3779B97F4A7C15u64;
    let mut q: Vec<f32> = (0..g * b).map(|_| { st ^= st << 13; st ^= st >> 7; st ^= st << 17; ((st >> 11) as f64 / (1u64 << 53) as f64 - 0.5) as f32 }).collect();
    let orth = |y: &mut Vec<f32>| {
        for j in 0..b {
            for _ in 0..2 { for i in 0..j { let d: f64 = (0..g).map(|r| y[r * b + i] as f64 * y[r * b + j] as f64).sum(); for r in 0..g { y[r * b + j] -= (d * y[r * b + i] as f64) as f32; } } }
            let nn = (0..g).map(|r| (y[r * b + j] as f64).powi(2)).sum::<f64>().sqrt().max(1e-30);
            for r in 0..g { y[r * b + j] = (y[r * b + j] as f64 / nn) as f32; }
        }
    };
    orth(&mut q);
    for _ in 0..iters { q = sgemm(false, false, g, b, g, c, &q); orth(&mut q); }
    let sq = sgemm(false, false, g, b, g, c, &q);
    let small: Vec<f64> = sgemm(true, false, b, b, g, &q, &sq).iter().map(|&v| v as f64).collect();
    let (u, _) = la::top_eig(&small, b, k)?;
    let u32v: Vec<f32> = u.iter().map(|&v| v as f32).collect();
    Ok(sgemm(false, false, g, k, b, &q, &u32v))
}
/// ridge weights (|u| x k) from the Gram matrix, penalty alpha x mean diagonal of the panel block
pub(crate) fn ridge32(c: &[f32], g: usize, s: &[usize], u: &[usize], alpha: f64) -> Result<Vec<f32>> {
    let k = s.len();
    let mut a = vec![0.0f64; k * k];
    for (i, &si) in s.iter().enumerate() { for (j, &sj) in s.iter().enumerate() { a[i * k + j] = c[si * g + sj] as f64; } }
    let tr: f64 = (0..k).map(|i| a[i * k + i]).sum::<f64>() / k as f64;
    for i in 0..k { a[i * k + i] += alpha * tr; }
    let mut bt = vec![0.0f64; u.len() * k];
    bt.par_chunks_mut(k).zip(u).for_each(|(r, &ui)| { for (j, &sj) in s.iter().enumerate() { r[j] = c[ui * g + sj] as f64; } });
    la::posv(a, k, &mut bt, u.len())?;
    Ok(bt.iter().map(|&v| v as f32).collect())
}
/// greedy MMSE selection with unit noise, as the pipeline (sim::greedy), in f32
pub(crate) fn greedy32(c: &[f32], nref: usize, g: usize, k: usize) -> Vec<usize> {
    let mut r: Vec<f32> = c.iter().map(|v| v / nref as f32).collect();
    for j in 0..g { r[j * g + j] = (r[j * g + j] - 1.0).max(0.0); }
    let mut norm: Vec<f64> = r.par_chunks(g).map(|row| row.iter().map(|&v| (v as f64) * (v as f64)).sum()).collect();
    let (mut sel, mut used) = (vec![], vec![false; g]);
    for _ in 0..k {
        let best = (0..g).into_par_iter().filter(|&j| !used[j]).map(|j| (norm[j] / (r[j * g + j] as f64 + 1.0), j))
            .reduce(|| (f64::NEG_INFINITY, usize::MAX), |a, b| if b.0 > a.0 || (b.0 == a.0 && b.1 < a.1) { b } else { a });
        let j = best.1;
        used[j] = true; sel.push(j);
        let col: Vec<f32> = r[j * g..(j + 1) * g].to_vec();
        let d = col[j] + 1.0;
        r.par_chunks_mut(g).zip(norm.par_iter_mut()).enumerate().for_each(|(i, (row, nm))| {
            let f = col[i] / d; let mut s = 0.0f64;
            for (v, cc) in row.iter_mut().zip(&col) { *v -= f * cc; s += (*v as f64) * (*v as f64); }
            *nm = s;
        });
    }
    sel
}

/// sparse graph in CSR: row = target gene, entries = (source gene, weight)
pub struct Csr { pub indptr: Vec<u32>, pub idx: Vec<u32>, pub val: Vec<f32> }
impl Csr {
    fn from_rows(rows: Vec<Vec<(u32, f32)>>) -> Csr {
        let mut c = Csr { indptr: vec![0], idx: vec![], val: vec![] };
        for r in rows { for (i, v) in r { c.idx.push(i); c.val.push(v); } c.indptr.push(c.idx.len() as u32); }
        c
    }
    fn nnz(&self) -> usize { self.idx.len() }
}
/// one propagation step on the CPU: out[j, b] = gamma * sum_e val[e] * cur[idx[e], b] + (1 - gamma) * r0[j, b], clamped rows copied from r0
fn spmm_cpu(w: &Csr, cur: &[f32], r0: &[f32], clamp: &[bool], nb: usize, gamma: f32) -> Vec<f32> {
    let g = clamp.len();
    let mut out = vec![0.0f32; g * nb];
    out.par_chunks_mut(nb).enumerate().for_each(|(j, o)| {
        if clamp[j] { o.copy_from_slice(&r0[j * nb..(j + 1) * nb]); return; }
        for e in w.indptr[j] as usize..w.indptr[j + 1] as usize {
            let (i, v) = (w.idx[e] as usize, w.val[e] * gamma);
            for (x, c) in o.iter_mut().zip(&cur[i * nb..(i + 1) * nb]) { *x += v * c; }
        }
        if gamma < 1.0 { for (x, c) in o.iter_mut().zip(&r0[j * nb..(j + 1) * nb]) { *x += (1.0 - gamma) * c; } }
    });
    out
}

#[cfg(feature = "metal")]
mod gpu {
    use super::Csr;
    use metal::*;
    const SRC: &str = r#"
#include <metal_stdlib>
using namespace metal;
struct P { uint g; uint nb; float gamma; uint pad; };
// one thread per (gene j, perturbation b): CSR row j of the graph times column b of the state
kernel void spmm(device const uint* indptr [[buffer(0)]], device const uint* idx [[buffer(1)]], device const float* val [[buffer(2)]],
                 device const float* cur [[buffer(3)]], device const float* r0 [[buffer(4)]], device const uchar* clamp [[buffer(5)]],
                 device float* out [[buffer(6)]], constant P& p [[buffer(7)]], uint2 gid [[thread_position_in_grid]]) {
    uint b = gid.x, j = gid.y;
    if (j >= p.g || b >= p.nb) return;
    uint o = j * p.nb + b;
    if (clamp[j]) { out[o] = r0[o]; return; }
    float acc = 0.0f;
    for (uint e = indptr[j]; e < indptr[j + 1]; e++) acc += val[e] * cur[idx[e] * p.nb + b];
    out[o] = p.gamma * acc + (1.0f - p.gamma) * r0[o];
}
"#;
    pub struct Gpu { dev: Device, pso: ComputePipelineState, q: CommandQueue }
    impl Gpu {
        pub fn new() -> anyhow::Result<Gpu> {
            let dev = Device::system_default().ok_or_else(|| anyhow::anyhow!("no Metal device"))?;
            let lib = dev.new_library_with_source(SRC, &CompileOptions::new()).map_err(|e| anyhow::anyhow!(e))?;
            let pso = dev.new_compute_pipeline_state_with_function(&*lib.get_function("spmm", None).map_err(|e| anyhow::anyhow!(e))?).map_err(|e| anyhow::anyhow!(e))?;
            let q = dev.new_command_queue();
            Ok(Gpu { dev, pso, q })
        }
        fn buf<T: Copy>(&self, v: &[T]) -> Buffer { self.dev.new_buffer_with_data(v.as_ptr() as *const _, (std::mem::size_of_val(v).max(4)) as u64, MTLResourceOptions::StorageModeShared) }
        /// `hops` propagation steps on shared (unified memory) buffers; returns the final state
        pub fn propagate(&self, w: &Csr, start: &[f32], r0: &[f32], clamp: &[bool], nb: usize, gamma: f32, hops: usize) -> Vec<f32> {
            let g = clamp.len();
            let cl: Vec<u8> = clamp.iter().map(|&c| c as u8).collect();
            let (bp, bi, bv, br, bc) = (self.buf(&w.indptr), self.buf(&w.idx), self.buf(&w.val), self.buf(r0), self.buf(&cl));
            let mut a = self.buf(start);
            let mut b = self.dev.new_buffer((g * nb * 4) as u64, MTLResourceOptions::StorageModeShared);
            let p = [g as u32, nb as u32, gamma.to_bits(), 0u32];
            for _ in 0..hops {
                let cb = self.q.new_command_buffer();
                let enc = cb.new_compute_command_encoder();
                enc.set_compute_pipeline_state(&self.pso);
                for (i, x) in [&bp, &bi, &bv, &a, &br, &bc, &b].iter().enumerate() { enc.set_buffer(i as u64, Some(x), 0); }
                enc.set_bytes(7, 16, p.as_ptr() as *const _);
                let (tw, th) = (self.pso.thread_execution_width(), self.pso.max_total_threads_per_threadgroup() / self.pso.thread_execution_width());
                enc.dispatch_threads(MTLSize::new(nb as u64, g as u64, 1), MTLSize::new(tw, th, 1));
                enc.end_encoding();
                cb.commit();
                cb.wait_until_completed();
                std::mem::swap(&mut a, &mut b);
            }
            unsafe { std::slice::from_raw_parts(a.contents() as *const f32, g * nb) }.to_vec()
        }
    }
}

#[cfg(not(feature = "metal"))]
mod gpu {
    use super::Csr;
    /// placeholder when the crate is built without the `metal` feature
    pub struct Gpu;
    impl Gpu {
        pub fn new() -> anyhow::Result<Gpu> { anyhow::bail!("built without the metal feature (cargo build --release --features metal)") }
        pub fn propagate(&self, _w: &Csr, _start: &[f32], _r0: &[f32], _clamp: &[bool], _nb: usize, _gamma: f32, _hops: usize) -> Vec<f32> { unreachable!() }
    }
}

pub(crate) struct Data { pub(crate) perts: Vec<String>, pub(crate) genes: Vec<String>, map: Mmap }
impl Data {
    pub(crate) fn open(pre: &str) -> Result<Data> {
        let col0 = |f: String| -> Result<Vec<String>> { Ok(std::fs::read_to_string(&f)?.lines().skip(1).map(|l| l.split('\t').next().unwrap().to_string()).collect()) };
        let (perts, genes) = (col0(format!("{pre}.meta.tsv"))?, col0(format!("{pre}.genes.tsv"))?);
        let map = unsafe { Mmap::map(&std::fs::File::open(format!("{pre}.z.f32"))?)? };
        if map.len() != 3 * perts.len() * genes.len() * 4 { bail!("z size") }
        Ok(Data { perts, genes, map })
    }
    pub(crate) fn z(&self) -> &[f32] { bytemuck::cast_slice(&self.map) }
    /// rows x all genes of one layer (0 all cells, 1 input half, 2 held-out half)
    pub(crate) fn rows(&self, layer: usize, rows: &[usize]) -> Vec<f32> {
        let (g, n, z) = (self.genes.len(), self.perts.len(), self.z());
        let mut m = vec![0.0f32; rows.len() * g];
        m.par_chunks_mut(g).zip(rows).for_each(|(r, &p)| r.copy_from_slice(&z[(layer * n + p) * g..(layer * n + p + 1) * g]));
        m
    }
}
pub(crate) fn cols(m: &[f32], g: usize, c: &[usize]) -> Vec<f32> { m.par_chunks(g).flat_map_iter(|r| c.iter().map(move |&j| r[j])).collect() }

/// everything fitted on one training set (inner training or all references)
struct Fit { w: Vec<f32>, v: Vec<f32>, e: Vec<f32>, n: usize, tr: Vec<usize> }
/// residual training profiles after removal of the r leading factors of the Gram matrix
fn fit_train(d: &Data, tr: &[usize], s: &[usize], u: &[usize], alpha: f64, r: usize) -> Result<(Fit, Vec<f32>)> {
    let g = d.genes.len();
    let y = d.rows(0, tr);
    let c = sgram(&y, tr.len(), g);
    let w = ridge32(&c, g, s, u, alpha)?;
    let v = top_eig32(&c, g, r, 30)?;
    let mut e = y;
    for blk in e.chunks_mut(512 * g) {
        let nb = blk.len() / g;
        let back = sgemm(false, true, nb, g, r, &sgemm(false, false, nb, r, g, blk, &v), &v);
        blk.iter_mut().zip(&back).for_each(|(a, b)| *a -= b);
    }
    Ok((Fit { w, v, e, n: tr.len(), tr: tr.to_vec() }, c))
}
/// conditional-dependency graph: regression form b_ji = -Theta_ji / Theta_jj of the regularised inverse
/// covariance of the residual profiles, top-k per gene by absolute partial correlation
fn graph_cond(f: &Fit, g: usize, lambda: f32, ks: &[usize]) -> Result<Vec<Csr>> {
    let mut a = sgram(&f.e, f.n, g);
    let md: f32 = (0..g).map(|j| a[j * g + j]).sum::<f32>() / g as f32;
    for j in 0..g { a[j * g + j] += lambda * md; }
    let (n32, mut info) = (g as i32, 0);
    unsafe { spotrf_(b"L".as_ptr(), &n32, a.as_mut_ptr(), &n32, &mut info) };
    if info != 0 { bail!("spotrf {info}") }
    unsafe { spotri_(b"L".as_ptr(), &n32, a.as_mut_ptr(), &n32, &mut info) };
    if info != 0 { bail!("spotri {info}") }
    // column-major lower triangle = row-major upper triangle
    let th = |i: usize, j: usize| if i <= j { a[i * g + j] } else { a[j * g + i] };
    let kmax = *ks.iter().max().unwrap();
    let top: Vec<Vec<(u32, f32, f32)>> = (0..g).into_par_iter().map(|j| {
        let tjj = th(j, j);
        let mut c: Vec<(u32, f32, f32)> = (0..g).filter(|&i| i != j).map(|i| { let t = th(i, j); (i as u32, -t / tjj, (t / (th(i, i) * tjj).sqrt()).abs()) }).collect();
        c.select_nth_unstable_by(kmax, |x, y| y.2.partial_cmp(&x.2).unwrap());
        c.truncate(kmax);
        c.sort_by(|x, y| y.2.partial_cmp(&x.2).unwrap());
        c
    }).collect();
    Ok(ks.iter().map(|&k| Csr::from_rows(top.iter().map(|r| r[..k].iter().map(|x| (x.0, x.1)).collect()).collect())).collect())
}
/// directed graph from knockdowns whose target gene is measured: per-unit total effects
/// M0[i, j] = E[i, j] / z[i, i] (zero for j = i), direct effects D = (I + beta M0_SS)^-1 M0, top-k incoming per gene
fn graph_dir(d: &Data, f: &Fit, g: usize, beta: f32, ks: &[usize], tau: f32) -> Result<(Vec<Csr>, usize)> {
    let gi: HashMap<&str, usize> = d.genes.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let z0 = d.z();
    // regulators: training knockdowns with a measured target gene and a strong on-target effect (z < -6, all cells)
    let reg: Vec<(usize, usize, f32)> = f.tr.iter().enumerate().filter_map(|(row, &p)| gi.get(d.perts[p].as_str()).map(|&gj| (row, gj, z0[p * g + gj]))).filter(|x| x.2 < -6.0).collect();
    let ns = reg.len();
    // rhs[j * ns + i] = M0[i, j]: column-major for LAPACK, and the incoming edges of gene j are contiguous
    let mut rhs = vec![0.0f32; g * ns];
    rhs.par_chunks_mut(ns).enumerate().for_each(|(j, col)| { for (i, &(row, gj, zii)) in reg.iter().enumerate() { if gj != j { col[i] = f.e[row * g + j] / zii; } } });
    if beta > 0.0 {
        let mut a = vec![0.0f32; ns * ns]; // column-major (I + beta M0_SS)
        a.par_chunks_mut(ns).enumerate().for_each(|(c, col)| { let gj = reg[c].1; for i in 0..ns { col[i] = beta * rhs[gj * ns + i]; } col[c] += 1.0; });
        let (n32, nrhs, mut info) = (ns as i32, g as i32, 0);
        let mut ipiv = vec![0i32; ns];
        unsafe { sgesv_(&n32, &nrhs, a.as_mut_ptr(), &n32, ipiv.as_mut_ptr(), rhs.as_mut_ptr(), &n32, &mut info) };
        if info != 0 { bail!("sgesv {info}") }
    }
    let kmax = *ks.iter().max().unwrap();
    let top: Vec<Vec<(u32, f32)>> = (0..g).into_par_iter().map(|j| {
        // candidate edges: observed residual effect of the knockdown on gene j beyond tau
        let mut c: Vec<(u32, f32)> = (0..ns).filter(|&i| reg[i].1 != j && f.e[reg[i].0 * g + j].abs() > tau).map(|i| (reg[i].1 as u32, rhs[j * ns + i])).filter(|x| x.1.is_finite()).collect();
        c.sort_by(|x, y| y.1.abs().partial_cmp(&x.1.abs()).unwrap());
        c.truncate(kmax);
        c
    }).collect();
    Ok((ks.iter().map(|&k| Csr::from_rows(top.iter().map(|r| r[..k.min(r.len())].to_vec()).collect())).collect(), ns))
}

/// state of a batch of profiles for the graph: residuals r0 (gene-major, g x nb) with panel genes observed
/// and the other genes at their ridge estimate, both after removal of the global factors
struct Batch { nb: usize, ridge: Vec<f32>, r0: Vec<f32> }
fn batch(d: &Data, f: &Fit, rows: &[usize], s: &[usize], u: &[usize]) -> Batch {
    let (g, nb, k, nu) = (d.genes.len(), rows.len(), s.len(), u.len());
    let xp = cols(&d.rows(1, rows), g, s);
    let ridge = sgemm(false, true, nb, nu, k, &xp, &f.w);
    let mut full = vec![0.0f32; nb * g];
    full.par_chunks_mut(g).enumerate().for_each(|(b, r)| { for (q, &j) in s.iter().enumerate() { r[j] = xp[b * k + q]; } for (q, &j) in u.iter().enumerate() { r[j] = ridge[b * nu + q]; } });
    let rr = f.v.len() / g;
    let glob = sgemm(false, true, nb, g, rr, &sgemm(false, false, nb, rr, g, &full, &f.v), &f.v);
    let mut r0 = vec![0.0f32; g * nb];
    r0.par_chunks_mut(nb).enumerate().for_each(|(j, o)| for b in 0..nb { o[b] = full[b * g + j] - glob[b * g + j]; });
    Batch { nb, ridge, r0 }
}
/// graph correction on the genes outside the panel after `hops` steps (nb x nu, row-major)
#[allow(clippy::too_many_arguments)]
fn propagate(w: &Csr, bt: &Batch, clamp: &[bool], u: &[usize], gamma: f32, hops: usize, init_ridge: bool, gpu: Option<&gpu::Gpu>) -> Vec<f32> {
    let (nb, g) = (bt.nb, clamp.len());
    let mut start = bt.r0.clone();
    if !init_ridge { start.par_chunks_mut(nb).enumerate().for_each(|(j, o)| if !clamp[j] { o.fill(0.0) }); }
    // the constant term of the damped iteration is the starting state
    let cur = match gpu { Some(gp) => gp.propagate(w, &start, &start, clamp, nb, gamma, hops), None => { let mut c = start.clone(); for _ in 0..hops { c = spmm_cpu(w, &c, &start, clamp, nb, gamma); } c } };
    let _ = g;
    let mut delta = vec![0.0f32; nb * u.len()];
    delta.par_chunks_mut(u.len()).enumerate().for_each(|(b, o)| for (q, &j) in u.iter().enumerate() { o[q] = cur[j * nb + b] - start[j * nb + b]; });
    delta
}
/// sparse combination of at most m training knockdowns (orthogonal matching pursuit on their panel residuals);
/// returns the corrections on the genes outside the panel for every m in ms
fn omp(f: &Fit, bt: &Batch, g: usize, s: &[usize], u: &[usize], ms: &[usize]) -> Vec<Vec<f32>> {
    let (k, nu, n, nb) = (s.len(), u.len(), f.n, bt.nb);
    let dict = cols(&f.e, g, s); // n x k, residual panel profiles of the training knockdowns
    let nrm: Vec<f32> = dict.par_chunks(k).map(|r| r.iter().map(|v| v * v).sum::<f32>().sqrt().max(1e-6)).collect();
    let mmax = *ms.iter().max().unwrap();
    let per: Vec<Vec<Vec<f32>>> = (0..nb).into_par_iter().map(|b| {
        let y: Vec<f32> = s.iter().map(|&j| bt.r0[j * nb + b]).collect();
        let (mut res, mut sel): (Vec<f32>, Vec<usize>) = (y.clone(), vec![]);
        let mut outs = vec![];
        for m in 1..=mmax {
            let best = (0..n).filter(|i| !sel.contains(i)).map(|i| ((dict[i * k..(i + 1) * k].iter().zip(&res).map(|(a, c)| a * c).sum::<f32>() / nrm[i]).abs(), i))
                .fold((f32::NEG_INFINITY, 0usize), |a, c| if c.0 > a.0 { c } else { a });
            sel.push(best.1);
            // least squares on the selected atoms (normal equations, small)
            let q = sel.len();
            let mut a = vec![0.0f64; q * q]; let mut rhs = vec![0.0f64; q];
            for (x, &i) in sel.iter().enumerate() {
                rhs[x] = dict[i * k..(i + 1) * k].iter().zip(&y).map(|(p, c)| (*p as f64) * (*c as f64)).sum();
                for (z, &i2) in sel.iter().enumerate() { a[x * q + z] = dict[i * k..(i + 1) * k].iter().zip(&dict[i2 * k..(i2 + 1) * k]).map(|(p, c)| (*p as f64) * (*c as f64)).sum(); }
                a[x * q + x] *= 1.0 + 1e-6;
            }
            if la::posv(a, q, &mut rhs, 1).is_err() { sel.pop(); break; }
            let coef = rhs;
            res = y.clone();
            for (x, &i) in sel.iter().enumerate() { for (r, p) in res.iter_mut().zip(&dict[i * k..(i + 1) * k]) { *r -= (coef[x] as f32) * p; } }
            if ms.contains(&m) {
                let mut o = vec![0.0f32; nu];
                for (x, &i) in sel.iter().enumerate() { for (q2, &j) in u.iter().enumerate() { o[q2] += (coef[x] as f32) * f.e[i * g + j]; } }
                outs.push(o);
            }
        }
        while outs.len() < ms.len() { outs.push(outs.last().cloned().unwrap_or(vec![0.0; nu])); }
        outs
    }).collect();
    (0..ms.len()).map(|q| per.iter().flat_map(|p| p[q].iter().copied()).collect()).collect()
}

pub(crate) struct Score { pub(crate) recovery: f64, pub(crate) spec: f64, pub(crate) recall: f64, pub(crate) opp: f64, pub(crate) spa: f64, pub(crate) n_resp: usize, pub(crate) spa_mean: f64, pub(crate) spa_v: Vec<f64>, pub(crate) spec_v: Vec<f64> }
/// the metrics of the manuscript for one prediction of the genes outside the panel
pub(crate) fn score(pred: &[f32], truth: &[f32], direct: &[f32], den: &[f32], nu: usize, basis: &[Vec<f64>]) -> Score {
    let nt = truth.len() / nu;
    let rows: Vec<(f64, f64, f64, usize, usize, usize, usize, f64)> = (0..nt).into_par_iter().map(|i| {
        let f = |m: &[f32]| -> Vec<f64> { m[i * nu..(i + 1) * nu].iter().map(|&v| v as f64).collect() };
        let (p, t, d, dn) = (f(pred), f(truth), f(direct), f(den));
        let rt = resid(&t, basis);
        let (sp, sd) = (pearson(&resid(&p, basis), &rt), pearson(&resid(&dn, basis), &rt));
        let nde = t.iter().filter(|v| v.abs() > DE_Z).count();
        let nc = d.iter().filter(|v| v.abs() > DE_Z).count();
        let mut o: Vec<usize> = (0..nu).collect();
        o.sort_by(|&a, &b| p[b].abs().partial_cmp(&p[a].abs()).unwrap().then(a.cmp(&b)));
        let (mut tp, mut disc) = (0, 0);
        for &j in &o[..nc.min(nu)] { let same = p[j].signum() == t[j].signum(); if !same { disc += 1; } if same && t[j].abs() > DE_Z { tp += 1; } }
        (pearson(&d, &t), sp, sd, nde, nc, tp, disc, spa_de(&p, &t))
    }).collect();
    let resp: Vec<&(f64, f64, f64, usize, usize, usize, usize, f64)> = rows.iter().filter(|r| r.0 >= 0.2).collect();
    let s = |f: &dyn Fn(&(f64, f64, f64, usize, usize, usize, usize, f64)) -> usize| rows.iter().map(f).sum::<usize>() as f64;
    Score { recovery: median(resp.iter().map(|r| r.1 / r.2).collect()), spec: median(resp.iter().map(|r| r.1).collect()), recall: s(&|r| r.5) / s(&|r| r.3), opp: s(&|r| r.6) / s(&|r| r.4),
            spa: median(rows.iter().map(|r| r.7).collect()), n_resp: resp.len(), spa_mean: { let v: Vec<f64> = rows.iter().map(|r| r.7).filter(|x| x.is_finite()).collect(); v.iter().sum::<f64>() / v.len() as f64 },
            spa_v: rows.iter().map(|r| r.7).collect(), spec_v: rows.iter().map(|r| if r.0 >= 0.2 { r.1 } else { f64::NAN }).collect() }
}
pub(crate) fn add(base: &[f32], terms: &[(&[f32], f32)]) -> Vec<f32> { let mut o = base.to_vec(); for (d, s) in terms { o.par_iter_mut().zip(d.par_iter()).for_each(|(a, b)| *a += s * b); } o }
/// least-squares scale of a correction against the ridge error (through the origin)
pub(crate) fn scale(delta: &[f32], err: &[f32]) -> f32 { let (mut xy, mut xx) = (0.0f64, 0.0f64); for (d, e) in delta.iter().zip(err) { xy += (*d as f64) * (*e as f64); xx += (*d as f64) * (*d as f64); } if xx > 0.0 { (xy / xx) as f32 } else { 0.0 } }
/// mean over perturbations of the rank correlation of absolute changes among DE genes
pub(crate) fn spa_med(pred: &[f32], truth: &[f32], nu: usize) -> f64 {
    let mean = |v: Vec<f64>| { let x: Vec<f64> = v.into_iter().filter(|a| a.is_finite()).collect(); x.iter().sum::<f64>() / x.len().max(1) as f64 };
    mean(pred.par_chunks(nu).zip(truth.par_chunks(nu)).map(|(p, t)| spa_de(&p.iter().map(|&v| v as f64).collect::<Vec<_>>(), &t.iter().map(|&v| v as f64).collect::<Vec<_>>())).collect())
}

/// resident memory of this process in MB (ps)
pub(crate) fn rss_mb() -> f64 {
    std::process::Command::new("ps").args(["-o", "rss=", "-p", &std::process::id().to_string()]).output().ok().and_then(|o| String::from_utf8_lossy(&o.stdout).trim().parse::<f64>().ok()).map(|k| k / 1024.0).unwrap_or(f64::NAN)
}

pub fn run(a: &[String]) -> Result<()> {
    let t0 = std::time::Instant::now();
    let d = Data::open(&req(a, "--z")?)?;
    let out = req(a, "--out")?;
    let f: i64 = arg(a, "--fold").map(|v| v.parse().unwrap()).unwrap_or(0);
    let (g, np) = (d.genes.len(), d.perts.len());
    let pi: HashMap<&str, usize> = d.perts.iter().enumerate().map(|(i, p)| (p.as_str(), i)).collect();
    let (mut fold, mut group) = (vec![-1i64; np], vec![-1i64; np]);
    for l in std::fs::read_to_string(req(a, "--folds")?)?.lines().skip(1) { let x: Vec<&str> = l.split('\t').collect(); if let Some(&i) = pi.get(x[0]) { group[i] = x[1].parse()?; fold[i] = x[2].parse()?; } }
    let mut excl: HashSet<(usize, i64)> = HashSet::new();
    for l in std::fs::read_to_string(req(a, "--hubs")?)?.lines().skip(1) { let x: Vec<&str> = l.split('\t').collect(); if let Some(&i) = pi.get(x[0]) { excl.insert((i, x[1].parse()?)); } }
    let refs: Vec<usize> = (0..np).filter(|&i| fold[i] != f && !excl.contains(&(i, f))).collect();
    let tests: Vec<usize> = (0..np).filter(|&i| fold[i] == f).collect();
    let key = |i: usize| -> u64 { let k = if group[i] >= 0 { group[i] as u64 } else { 1_000_000 + i as u64 }; (k.wrapping_mul(0x9E3779B97F4A7C15) >> 32) % 10 };
    let (itr, iho): (Vec<usize>, Vec<usize>) = (refs.iter().copied().filter(|&i| key(i) != 0).collect(), refs.iter().copied().filter(|&i| key(i) == 0).collect());
    let gpu = if a.contains(&"--metal".to_string()) { Some(gpu::Gpu::new()?) } else { None };
    let mut log = BufWriter::new(std::fs::File::create(format!("{out}.causal.fold{f}.log"))?);
    macro_rules! say { ($($t:tt)*) => {{ let s = format!($($t)*); eprintln!("{s}"); writeln!(log, "{s}")?; }} }
    say!("fold {f}: references {} (inner training {}, inner holdout {}), tests {}, genes {g}", refs.len(), itr.len(), iho.len(), tests.len());

    // ---- panel (greedy MMSE on all references) and ridge penalty (inner holdout), as the pipeline
    let (s, u, alpha) = {
        let y = d.rows(0, &refs);
        let c = sgram(&y, refs.len(), g);
        drop(y);
        let s = greedy32(&c, refs.len(), g, 500);
        drop(c);
        let ins: HashSet<usize> = s.iter().copied().collect();
        let u: Vec<usize> = (0..g).filter(|j| !ins.contains(j)).collect();
        let c_in = sgram(&d.rows(0, &itr), itr.len(), g);
        let (xh, th) = (cols(&d.rows(1, &iho), g, &s), cols(&d.rows(2, &iho), g, &u));
        let mut best = (f64::NEG_INFINITY, ALPHAS[0]);
        for &al in &ALPHAS {
            let p = sgemm(false, true, iho.len(), u.len(), 500, &xh, &ridge32(&c_in, g, &s, &u, al)?);
            let r: f64 = p.par_chunks(u.len()).zip(th.par_chunks(u.len())).map(|(x, t)| pearson(&x.iter().map(|&v| v as f64).collect::<Vec<_>>(), &t.iter().map(|&v| v as f64).collect::<Vec<_>>())).filter(|v| v.is_finite()).sum::<f64>() / iho.len() as f64;
            if r > best.0 { best = (r, al); }
        }
        (s, u, best.1)
    };
    let nu = u.len();
    let mut clamp = vec![false; g];
    for &j in &s { clamp[j] = true; }
    say!("panel 500 genes, ridge alpha {alpha} at {:.1}s, resident {:.0} MB", t0.elapsed().as_secs_f64(), rss_mb());

    // ---- tuning on the inner holdout: graphs built on the inner training set only
    let (rs, ks, lambdas, betas, gammas, hops_g, ms) = ([20usize, 50], [16usize, 32, 64], [0.1f32, 1.0], [0.0f32, 0.5, 1.0], [0.5f32, 1.0], [1usize, 2, 3], [1usize, 2, 4, 8]);
    let th = cols(&d.rows(2, &iho), g, &u);
    #[derive(Clone, Debug)]
    struct Cfg { r: usize, k: usize, par: f32, gamma: f32, hops: usize, init_ridge: bool, s: f32, spa: f64 }
    let none = Cfg { r: 0, k: 0, par: 0.0, gamma: 0.0, hops: 0, init_ridge: false, s: 0.0, spa: f64::NEG_INFINITY };
    let (mut b_dir, mut b_cond, mut b_omp) = (none.clone(), none.clone(), none.clone());
    let mut base_spa;
    let mut t_spmm = (0.0f64, 0usize);
    for &r in &rs {
        let (fi, _c) = fit_train(&d, &itr, &s, &u, alpha, r)?;
        drop(_c);
        let bt = batch(&d, &fi, &iho, &s, &u);
        let err: Vec<f32> = th.iter().zip(&bt.ridge).map(|(t, p)| t - p).collect();
        base_spa = spa_med(&bt.ridge, &th, nu);
        let mut try_graph = |w: &Csr, best: &mut Cfg, k: usize, par: f32| {
            for &gamma in &gammas { for &hops in &hops_g { for init_ridge in [false, true] {
                let ts = std::time::Instant::now();
                let dl = propagate(w, &bt, &clamp, &u, gamma, hops, init_ridge, gpu.as_ref());
                t_spmm.0 += ts.elapsed().as_secs_f64(); t_spmm.1 += hops;
                let sc = scale(&dl, &err);
                let spa = spa_med(&add(&bt.ridge, &[(&dl, sc)]), &th, nu);
                if spa > best.spa { *best = Cfg { r, k, par, gamma, hops, init_ridge, s: sc, spa }; }
            } } }
        };
        for &lam in &lambdas { let gs = graph_cond(&fi, g, lam, &ks)?; for (w, &k) in gs.iter().zip(&ks) { try_graph(w, &mut b_cond, k, lam); } }
        for &beta in &betas { let (gs, ns) = graph_dir(&d, &fi, g, beta, &ks, 3.0)?; for (w, &k) in gs.iter().zip(&ks) { try_graph(w, &mut b_dir, k, beta); }
            if beta == 0.0 { say!("R {r}: {ns} regulators with a measured target; directed edges at k 64: {}", gs[2].nnz()); } }
        for (dl, &m) in omp(&fi, &bt, g, &s, &u, &ms).iter().zip(&ms) {
            let sc = scale(dl, &err);
            let spa = spa_med(&add(&bt.ridge, &[(dl, sc)]), &th, nu);
            if spa > b_omp.spa { b_omp = Cfg { r, k: m, par: 0.0, gamma: 0.0, hops: 0, init_ridge: false, s: sc, spa }; }
        }
        say!("R {r} tuned at {:.1}s, resident {:.0} MB (inner-holdout mean size rank: ridge {base_spa:.4}, best dir {:.4}, cond {:.4}, omp {:.4})", t0.elapsed().as_secs_f64(), rss_mb(), b_dir.spa, b_cond.spa, b_omp.spa);
    }
    if let Some(gp) = gpu.as_ref() {
        let (fi, _c) = fit_train(&d, &itr, &s, &u, alpha, 20)?;
        drop(_c);
        let bt = batch(&d, &fi, &iho, &s, &u);
        let w = &graph_cond(&fi, g, 0.1, &[64])?[0];
        let t1 = std::time::Instant::now(); let mut c = bt.r0.clone(); for _ in 0..3 { c = spmm_cpu(w, &c, &bt.r0, &clamp, bt.nb, 0.5); } let t_cpu = t1.elapsed().as_secs_f64();
        let _warm = gp.propagate(w, &bt.r0, &bt.r0, &clamp, bt.nb, 0.5, 1);
        let t2 = std::time::Instant::now(); let m = gp.propagate(w, &bt.r0, &bt.r0, &clamp, bt.nb, 0.5, 3); let t_gpu = t2.elapsed().as_secs_f64();
        let dmax = c.iter().zip(&m).map(|(x, y)| (x - y).abs()).fold(0.0f32, f32::max);
        say!("backend check, 3 steps, {} edges x {} profiles: CPU {:.1} ms, Metal {:.1} ms (with buffer setup), largest difference {dmax:.2e}", w.nnz(), bt.nb, 1e3 * t_cpu, 1e3 * t_gpu);
    }
    say!("chosen dir  {:?}", b_dir); say!("chosen cond {:?}", b_cond); say!("chosen omp  {:?}", b_omp);
    say!("propagation: {:.2} ms per step on the inner holdout ({} profiles x {g} genes), backend {}", 1e3 * t_spmm.0 / t_spmm.1.max(1) as f64, iho.len(), if cfg!(feature = "metal") && a.contains(&"--metal".to_string()) { "Metal" } else { "CPU" });
    // combined correction: least squares of the ridge error on the three corrections (inner holdout)
    let comb = {
        let mut dl: Vec<Vec<f32>> = vec![];
        let mut errs: Vec<f32> = vec![];
        for (which, cf) in [(0, &b_dir), (1, &b_cond), (2, &b_omp)] {
            let (fi, _c) = fit_train(&d, &itr, &s, &u, alpha, cf.r)?;
            drop(_c);
            let bt = batch(&d, &fi, &iho, &s, &u);
            if errs.is_empty() { errs = th.iter().zip(&bt.ridge).map(|(t, p)| t - p).collect(); }
            dl.push(match which {
                0 => propagate(&graph_dir(&d, &fi, g, cf.par, &[cf.k], 3.0)?.0[0], &bt, &clamp, &u, cf.gamma, cf.hops, cf.init_ridge, gpu.as_ref()),
                1 => propagate(&graph_cond(&fi, g, cf.par, &[cf.k])?[0], &bt, &clamp, &u, cf.gamma, cf.hops, cf.init_ridge, gpu.as_ref()),
                _ => omp(&fi, &bt, g, &s, &u, &[cf.k]).remove(0),
            });
        }
        let mut m = vec![0.0f64; 9]; let mut rhs = vec![0.0f64; 3];
        for x in 0..3 { rhs[x] = dl[x].iter().zip(&errs).map(|(p, e)| (*p as f64) * (*e as f64)).sum(); for y in 0..3 { m[x * 3 + y] = dl[x].iter().zip(&dl[y]).map(|(p, q)| (*p as f64) * (*q as f64)).sum(); } m[x * 3 + x] *= 1.0 + 1e-9; }
        la::posv(m, 3, &mut rhs, 1)?;
        say!("combined weights (dir, cond, omp): {:.4} {:.4} {:.4}", rhs[0], rhs[1], rhs[2]);
        rhs
    };

    // ---- final models on all references, evaluation on the test fold
    let (tt, dd) = (cols(&d.rows(2, &tests), g, &u), cols(&d.rows(1, &tests), g, &u));
    let mut res: Vec<(String, Score, f64)> = vec![];
    let mut deltas: Vec<Vec<f32>> = vec![];
    let mut ridge_pred: Vec<f32> = vec![];
    let (mut basis, mut den): (Vec<Vec<f64>>, Vec<f32>) = (vec![], vec![]);
    let mut cache: Option<(usize, Fit, Batch)> = None;
    for (which, cf) in [(0, &b_dir), (1, &b_cond), (2, &b_omp)] {
        if cache.as_ref().map(|c| c.0) != Some(cf.r) {
            drop(cache.take()); // release the previous fit before the next one is built
            let (fi, c) = fit_train(&d, &refs, &s, &u, alpha, cf.r)?;
            if basis.is_empty() {
                // evaluation frame of the manuscript: 10 leading programs removed, 160 programs for the denoised profile
                let v160 = top_eig32(&c, g, 160, 30)?;
                let v10: Vec<f64> = (0..g).flat_map(|j| v160[j * 160..j * 160 + NPC].iter().map(|&v| v as f64).collect::<Vec<_>>()).collect();
                basis = basis_on(&v10, NPC, &u);
                let x = d.rows(1, &tests);
                den = cols(&sgemm(false, true, tests.len(), g, 160, &sgemm(false, false, tests.len(), 160, g, &x, &v160), &v160), g, &u);
            }
            drop(c);
            let bt = batch(&d, &fi, &tests, &s, &u);
            cache = Some((cf.r, fi, bt));
        }
        let (_, fi, bt) = cache.as_ref().unwrap();
        if ridge_pred.is_empty() { ridge_pred = bt.ridge.clone(); res.push(("ridge".into(), score(&ridge_pred, &tt, &dd, &den, nu, &basis), 0.0)); }
        let ts = std::time::Instant::now();
        let dl = match which {
            0 => propagate(&graph_dir(&d, fi, g, cf.par, &[cf.k], 3.0)?.0[0], bt, &clamp, &u, cf.gamma, cf.hops, cf.init_ridge, gpu.as_ref()),
            1 => propagate(&graph_cond(fi, g, cf.par, &[cf.k])?[0], bt, &clamp, &u, cf.gamma, cf.hops, cf.init_ridge, gpu.as_ref()),
            _ => omp(fi, bt, g, &s, &u, &[cf.k]).remove(0),
        };
        let ms_ = 1e3 * ts.elapsed().as_secs_f64();
        res.push((["ridge + dir", "ridge + cond", "ridge + omp"][which].into(), score(&add(&ridge_pred, &[(&dl, cf.s)]), &tt, &dd, &den, nu, &basis), ms_));
        deltas.push(dl);
    }
    res.push(("ridge + dir + cond + omp".into(), score(&add(&ridge_pred, &[(&deltas[0], comb[0] as f32), (&deltas[1], comb[1] as f32), (&deltas[2], comb[2] as f32)]), &tt, &dd, &den, nu, &basis), 0.0));
    res.insert(1, ("direct (input half)".into(), score(&dd, &tt, &dd, &den, nu, &basis), 0.0));
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.causal.fold{f}.tsv"))?);
    writeln!(w, "method\tn_responding\trecovery\tspecific_r\tde_recall\topposite_sign\tsize_rank\tgraph_ms")?;
    say!("\n{:28} {:>6} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}", "method (fold test set)", "n_resp", "recovery", "spec_r", "DE_recall", "opposite", "rank_med", "rank_mean", "graph_ms");
    for (name, sc, ms_) in &res {
        say!("{name:28} {:>6} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {:>9.1}", sc.n_resp, sc.recovery, sc.spec, sc.recall, sc.opp, sc.spa, sc.spa_mean, ms_);
        writeln!(w, "{name}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.2}", sc.n_resp, sc.recovery, sc.spec, sc.recall, sc.opp, sc.spa, ms_)?;
    }
    say!("\npaired against ridge (test perturbations): size rank [n, mean difference, share higher, Wilcoxon P] | specific r of responders [n, median difference, share higher, P]");
    let base = &res[0].1;
    for (name, sc, _) in res.iter().skip(2) {
        let pair = |a: &[f64], b: &[f64]| -> Vec<f64> { a.iter().zip(b).filter(|(x, y)| x.is_finite() && y.is_finite()).map(|(x, y)| x - y).collect() };
        let (d1, d2) = (pair(&sc.spa_v, &base.spa_v), pair(&sc.spec_v, &base.spec_v));
        say!("{name:28} size rank {} {:+.4} {:.3} {:.2e} | specific r {} {:+.4} {:.3} {:.2e}", d1.len(), d1.iter().sum::<f64>() / d1.len() as f64, d1.iter().filter(|v| **v > 0.0).count() as f64 / d1.len() as f64, upg::wilcoxon(&d1),
             d2.len(), median(d2.clone()), d2.iter().filter(|v| **v > 0.0).count() as f64 / d2.len() as f64, upg::wilcoxon(&d2));
    }
    say!("resident at the end {:.0} MB", rss_mb());
    say!("total {:.1}s", t0.elapsed().as_secs_f64());
    Ok(())
}
