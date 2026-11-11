# Environment

- Rust: version 1.99.0, pinned in `rust-toolchain.toml` at the repository root (rustup installs it on first build).
- Linear algebra: Apple Accelerate on macOS (part of the system). On Linux, OpenBLAS with CBLAS and LAPACK
  (`apt-get install libopenblas-dev` on Debian or Ubuntu).
- Python 3.12 with `pip install -r environment/requirements.txt` for the verification and figure scripts.
- Figure export (kaleido) needs the system library `libexpat1`, which minimal Linux images lack (`apt-get install libexpat1`).
- `rsvg-convert` (librsvg) is needed only for the schematic of Fig. 1.
- The analyses of the manuscript were run on a MacBook Pro (Apple M3 Pro, 18 GB) with macOS and Accelerate.
