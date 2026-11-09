//! Thin wrappers over Apple Accelerate (CBLAS + LAPACK, LP64). Matrices are row-major Vec<f64>.
#![allow(clippy::too_many_arguments)]
extern "C" {
    fn cblas_dgemm(o: i32, ta: i32, tb: i32, m: i32, n: i32, k: i32, al: f64, a: *const f64, lda: i32, b: *const f64, ldb: i32, be: f64, c: *mut f64, ldc: i32);
    fn cblas_dsyrk(o: i32, up: i32, t: i32, n: i32, k: i32, al: f64, a: *const f64, lda: i32, be: f64, c: *mut f64, ldc: i32);
    fn dposv_(uplo: *const u8, n: *const i32, nrhs: *const i32, a: *mut f64, lda: *const i32, b: *mut f64, ldb: *const i32, info: *mut i32);
    fn dgeqp3_(m: *const i32, n: *const i32, a: *mut f64, lda: *const i32, jpvt: *mut i32, tau: *mut f64, work: *mut f64, lwork: *const i32, info: *mut i32);
    fn dsyevr_(jobz: *const u8, range: *const u8, uplo: *const u8, n: *const i32, a: *mut f64, lda: *const i32, vl: *const f64, vu: *const f64, il: *const i32, iu: *const i32, abstol: *const f64, m: *mut i32, w: *mut f64, z: *mut f64, ldz: *const i32, isuppz: *mut i32, work: *mut f64, lwork: *const i32, iwork: *mut i32, liwork: *const i32, info: *mut i32);
}
const RM: i32 = 101;
const NT: i32 = 111;
const TR: i32 = 112;
const UP: i32 = 121;

/// C (m x n) = op(A) op(B); A is (m x k) or (k x m) if ta; B is (k x n) or (n x k) if tb.
pub fn gemm(ta: bool, tb: bool, m: usize, n: usize, k: usize, a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut c = vec![0.0; m * n];
    let lda = if ta { m } else { k };
    let ldb = if tb { k } else { n };
    unsafe { cblas_dgemm(RM, if ta { TR } else { NT }, if tb { TR } else { NT }, m as i32, n as i32, k as i32, 1.0, a.as_ptr(), lda as i32, b.as_ptr(), ldb as i32, 0.0, c.as_mut_ptr(), n as i32) };
    c
}

/// A^T A for row-major A (r x c), full symmetric c x c.
pub fn gram(a: &[f64], r: usize, c: usize) -> Vec<f64> {
    let mut g = vec![0.0; c * c];
    unsafe { cblas_dsyrk(RM, UP, TR, c as i32, r as i32, 1.0, a.as_ptr(), c as i32, 0.0, g.as_mut_ptr(), c as i32) };
    for i in 0..c {
        for j in 0..i {
            g[i * c + j] = g[j * c + i];
        }
    }
    g
}

/// Solve (A) X = B for SPD A (n x n); bt is B^T row-major (nrhs x n), overwritten by X^T.
pub fn posv(mut a: Vec<f64>, n: usize, bt: &mut [f64], nrhs: usize) -> anyhow::Result<()> {
    let mut info = 0;
    unsafe { dposv_(b"L".as_ptr(), &(n as i32), &(nrhs as i32), a.as_mut_ptr(), &(n as i32), bt.as_mut_ptr(), &(n as i32), &mut info) };
    if info != 0 { anyhow::bail!("dposv info {info}") }
    Ok(())
}

/// Column-pivoted QR of row-major A (r x c); returns the pivot order (0-based column indices).
pub fn qr_pivots(a: &[f64], r: usize, c: usize) -> anyhow::Result<Vec<usize>> {
    // column-major copy
    let mut cm = vec![0.0; r * c];
    for i in 0..r {
        for j in 0..c {
            cm[j * r + i] = a[i * c + j];
        }
    }
    let (m, n) = (r as i32, c as i32);
    let mut jpvt = vec![0i32; c];
    let mut tau = vec![0.0; r.min(c)];
    let mut info = 0;
    let mut q = 0.0;
    unsafe { dgeqp3_(&m, &n, cm.as_mut_ptr(), &m, jpvt.as_mut_ptr(), tau.as_mut_ptr(), &mut q, &-1, &mut info) };
    let lw = q as i32;
    let mut work = vec![0.0; lw as usize];
    unsafe { dgeqp3_(&m, &n, cm.as_mut_ptr(), &m, jpvt.as_mut_ptr(), tau.as_mut_ptr(), work.as_mut_ptr(), &lw, &mut info) };
    if info != 0 { anyhow::bail!("dgeqp3 info {info}") }
    Ok(jpvt.iter().map(|&v| (v - 1) as usize).collect())
}

/// Top-k eigenvectors of symmetric S (n x n); returns row-major n x k, largest eigenvalue first,
/// and the eigenvalues.
pub fn top_eig(s: &[f64], n: usize, k: usize) -> anyhow::Result<(Vec<f64>, Vec<f64>)> {
    let mut a = s.to_vec();
    let (nn, il, iu) = (n as i32, (n - k + 1) as i32, n as i32);
    let mut m = 0;
    let mut w = vec![0.0; n];
    let mut z = vec![0.0; n * k];
    let mut isuppz = vec![0i32; 2 * k];
    let mut info = 0;
    let (mut wq, mut iq) = (0.0, 0i32);
    unsafe { dsyevr_(b"V".as_ptr(), b"I".as_ptr(), b"L".as_ptr(), &nn, a.as_mut_ptr(), &nn, &0.0, &0.0, &il, &iu, &0.0, &mut m, w.as_mut_ptr(), z.as_mut_ptr(), &nn, isuppz.as_mut_ptr(), &mut wq, &-1, &mut iq, &-1, &mut info) };
    let (lw, liw) = (wq as i32, iq);
    let mut work = vec![0.0; lw as usize];
    let mut iwork = vec![0i32; liw as usize];
    unsafe { dsyevr_(b"V".as_ptr(), b"I".as_ptr(), b"L".as_ptr(), &nn, a.as_mut_ptr(), &nn, &0.0, &0.0, &il, &iu, &0.0, &mut m, w.as_mut_ptr(), z.as_mut_ptr(), &nn, isuppz.as_mut_ptr(), work.as_mut_ptr(), &lw, iwork.as_mut_ptr(), &liw, &mut info) };
    if info != 0 || m as usize != k { anyhow::bail!("dsyevr info {info} m {m}") }
    // z column-major n x k, ascending eigenvalues -> row-major n x k descending
    let mut out = vec![0.0; n * k];
    for j in 0..k {
        for i in 0..n {
            out[i * k + j] = z[(k - 1 - j) * n + i];
        }
    }
    let vals = (0..k).map(|j| w[k - 1 - j]).collect();
    Ok((out, vals))
}

/// Leading k eigenpairs of a symmetric positive semi-definite n x n matrix by block subspace
/// iteration (all work in dgemm, which Accelerate runs on the AMX units across cores), with a
/// Rayleigh-Ritz step on the (k + 16)-dimensional subspace. Same output layout as `top_eig`.
pub fn top_eig_sub(s: &[f64], n: usize, k: usize, iters: usize) -> anyhow::Result<(Vec<f64>, Vec<f64>)> {
    let b = (k + 16).min(n);
    let mut st = 0x9E3779B97F4A7C15u64;
    let mut q: Vec<f64> = (0..n * b).map(|_| { st ^= st << 13; st ^= st >> 7; st ^= st << 17; (st >> 11) as f64 / (1u64 << 53) as f64 - 0.5 }).collect();
    let orth = |y: &mut Vec<f64>| {
        // modified Gram-Schmidt on the b columns of the row-major n x b matrix
        for j in 0..b {
            for i in 0..j { let d: f64 = (0..n).map(|r| y[r * b + i] * y[r * b + j]).sum(); for r in 0..n { y[r * b + j] -= d * y[r * b + i]; } }
            let nn = (0..n).map(|r| y[r * b + j].powi(2)).sum::<f64>().sqrt().max(1e-300);
            for r in 0..n { y[r * b + j] /= nn; }
        }
    };
    orth(&mut q);
    for _ in 0..iters { q = gemm(false, false, n, b, n, s, &q); orth(&mut q); }
    let sq = gemm(false, false, n, b, n, s, &q);
    let small = gemm(true, false, b, b, n, &q, &sq);
    let (u, vals) = top_eig(&small, b, k)?;
    Ok((gemm(false, false, n, k, b, &q, &u), vals))
}
