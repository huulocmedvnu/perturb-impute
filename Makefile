# perturb-impute: build, test, verify the manuscript numbers and regenerate the figures.
# Requirements: Rust (see rust-toolchain.toml), Python 3 with environment/requirements.txt,
# and on Linux the OpenBLAS development package (libopenblas-dev).
PYTHON ?= python3
BIN := target/release/perturb-impute

.PHONY: help build test verify figures all clean

help:  ## list the targets
	@grep -E '^[a-z]+:.*##' $(MAKEFILE_LIST) | awk -F':.*## ' '{printf "  make %-8s %s\n", $$1, $$2}'

build:  ## compile perturb-impute and the analysis binary in release mode
	cargo build --release

test: build  ## unit tests, end-to-end command-line test and a printed run on the toy data in tests/data
	cargo test --release
	sh tests/smoke.sh

verify:  ## recompute the numbers of the manuscript from source_data/ and print a PASS/FAIL table
	$(PYTHON) scripts/verify_manuscript_numbers.py

figures:  ## regenerate Figs. 2 to 8 and Fig. S1 (and Fig. 1 if rsvg-convert is installed) into output/
	$(PYTHON) scripts/make_figures.py output
	@if command -v rsvg-convert >/dev/null 2>&1; then $(PYTHON) scripts/make_flowchart.py output; else echo "rsvg-convert not found: Fig. 1 (schematic) skipped"; fi

all: build test verify figures  ## build, test, verify and figures in sequence

clean:  ## remove build products and regenerated figures
	cargo clean
	rm -rf output
