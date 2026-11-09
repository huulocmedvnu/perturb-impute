//! Metal GPU check: Gram matrix C = X' X (n x g, f32) with a tiled compute shader on shared
//! (zero-copy) buffers, compared with Accelerate sgemm (f32) and dsyrk (f64) on the AMX units.
use metal::*;
use std::time::Instant;
extern "C" {
    fn cblas_sgemm(o: i32, ta: i32, tb: i32, m: i32, n: i32, k: i32, al: f32, a: *const f32, lda: i32, b: *const f32, ldb: i32, be: f32, c: *mut f32, ldc: i32);
    fn cblas_dsyrk(o: i32, up: i32, t: i32, n: i32, k: i32, al: f64, a: *const f64, lda: i32, be: f64, c: *mut f64, ldc: i32);
}
const SRC: &str = r#"
#include <metal_stdlib>
using namespace metal;
#define T 32
// C[i,j] = sum_r X[r,i] * X[r,j]; X row-major n x g
kernel void gram(device const float* X [[buffer(0)]], device float* C [[buffer(1)]], constant uint2& dims [[buffer(2)]],
                 uint2 gid [[thread_position_in_grid]], uint2 lid [[thread_position_in_threadgroup]], uint2 tg [[threadgroup_position_in_grid]]) {
    uint n = dims.x, g = dims.y;
    threadgroup float A[T][T];
    threadgroup float B[T][T];
    uint i = tg.y * T + lid.y, j = tg.x * T + lid.x;
    float acc = 0.0f;
    for (uint r0 = 0; r0 < n; r0 += T) {
        uint ra = r0 + lid.x, rb = r0 + lid.y;
        A[lid.y][lid.x] = (ra < n && (tg.y * T + lid.y) < g) ? X[ra * g + tg.y * T + lid.y] : 0.0f;
        B[lid.y][lid.x] = (rb < n && j < g) ? X[rb * g + j] : 0.0f;
        threadgroup_barrier(mem_flags::mem_threadgroup);
        for (uint k = 0; k < T; k++) acc += A[lid.y][k] * B[k][lid.x];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (i < g && j < g) C[i * g + j] = acc;
}
"#;
fn main() -> anyhow::Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let n: usize = a.get(1).map(|v| v.parse().unwrap()).unwrap_or(18294);
    let g: usize = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(11948);
    let mut st = 7u64;
    let mut rnd = || { st ^= st << 13; st ^= st >> 7; st ^= st << 17; (st >> 40) as f32 / (1u64 << 24) as f32 - 0.5 };
    let dev = Device::system_default().ok_or_else(|| anyhow::anyhow!("no Metal device"))?;
    println!("device: {} (unified memory {})", dev.name(), dev.has_unified_memory());
    // shared buffers: the CPU writes X in place, the GPU reads it without a copy
    let bx = dev.new_buffer((n * g * 4) as u64, MTLResourceOptions::StorageModeShared);
    let bc = dev.new_buffer((g * g * 4) as u64, MTLResourceOptions::StorageModeShared);
    let x: &mut [f32] = unsafe { std::slice::from_raw_parts_mut(bx.contents() as *mut f32, n * g) };
    for v in x.iter_mut() { *v = rnd(); }
    let lib = dev.new_library_with_source(SRC, &CompileOptions::new()).map_err(|e| anyhow::anyhow!(e))?;
    let pso = dev.new_compute_pipeline_state_with_function(&*lib.get_function("gram", None).map_err(|e| anyhow::anyhow!(e))?).map_err(|e| anyhow::anyhow!(e))?;
    let q = dev.new_command_queue();
    let dims = [n as u32, g as u32];
    let run = || {
        let cb = q.new_command_buffer();
        let enc = cb.new_compute_command_encoder();
        enc.set_compute_pipeline_state(&pso);
        enc.set_buffer(0, Some(&bx), 0);
        enc.set_buffer(1, Some(&bc), 0);
        enc.set_bytes(2, 8, dims.as_ptr() as *const _);
        let groups = MTLSize::new(g.div_ceil(32) as u64, g.div_ceil(32) as u64, 1);
        enc.dispatch_thread_groups(groups, MTLSize::new(32, 32, 1));
        enc.end_encoding();
        cb.commit();
        cb.wait_until_completed();
    };
    run();
    let t = Instant::now(); run(); let t_gpu = t.elapsed().as_secs_f64();
    let cg: &[f32] = unsafe { std::slice::from_raw_parts(bc.contents() as *const f32, g * g) };
    let flops = 2.0 * n as f64 * g as f64 * g as f64;
    println!("Metal tiled f32 gram: {:.2} s ({:.2} TFLOP/s)", t_gpu, flops / t_gpu / 1e12);
    let mut cs = vec![0f32; g * g];
    let t = Instant::now();
    unsafe { cblas_sgemm(101, 112, 111, g as i32, g as i32, n as i32, 1.0, x.as_ptr(), g as i32, x.as_ptr(), g as i32, 0.0, cs.as_mut_ptr(), g as i32) };
    let t_s = t.elapsed().as_secs_f64();
    println!("Accelerate sgemm f32: {:.2} s ({:.2} TFLOP/s)", t_s, flops / t_s / 1e12);
    let xd: Vec<f64> = x.iter().map(|&v| v as f64).collect();
    let mut cd = vec![0f64; g * g];
    let t = Instant::now();
    unsafe { cblas_dsyrk(101, 121, 112, g as i32, n as i32, 1.0, xd.as_ptr(), g as i32, 0.0, cd.as_mut_ptr(), g as i32) };
    let t_d = t.elapsed().as_secs_f64();
    println!("Accelerate dsyrk f64 (upper triangle): {:.2} s ({:.2} TFLOP/s effective)", t_d, flops / 2.0 / t_d / 1e12);
    let mut mg = 0f64; let mut ms = 0f64; let mut mx = 0f64;
    for i in 0..g { for j in i..g { let r = cd[i * g + j]; mx = mx.max(r.abs()); mg = mg.max((cg[i * g + j] as f64 - r).abs()); ms = ms.max((cs[i * g + j] as f64 - r).abs()); } }
    println!("max abs error vs f64: Metal {:.2e}, sgemm {:.2e} (largest entry {:.1})", mg, ms, mx);
    Ok(())
}
