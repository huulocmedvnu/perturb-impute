// Linear algebra backend: Apple Accelerate on macOS, OpenBLAS (CBLAS and LAPACK) on other systems.
fn main() {
    let os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    if os == "macos" {
        println!("cargo:rustc-link-lib=framework=Accelerate");
    } else {
        println!("cargo:rustc-link-lib=openblas");
    }
}
