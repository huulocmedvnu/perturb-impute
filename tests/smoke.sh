#!/bin/sh
# Command-line workflow on the toy data (tests/data): two panels, two models, imputation, projection.
set -e
cd "$(dirname "$0")/.."
P=target/release/perturb-impute; D=tests/data; O=output/smoke
rm -rf "$O"; mkdir -p "$O"
$P select --ref $D/toy_reference.tsv --panel-size 40 --umi-fraction $D/toy_umi_fraction.tsv --out $O/panel_standard.tsv
$P select --ref $D/toy_reference.tsv --panel-size 40 --mode low-budget --umi-fraction $D/toy_umi_fraction.tsv --out $O/panel_low_budget.tsv
$P fit-reference --ref $D/toy_reference.tsv --panel $O/panel_standard.tsv --groups $D/toy_groups.tsv --out $O/model_standard.bin
$P fit-reference --ref $D/toy_reference.tsv --panel $O/panel_low_budget.tsv --groups $D/toy_groups.tsv --noise $D/toy_noise.tsv --budget 0.02 --out $O/model_low_budget.bin
$P impute --model $O/model_standard.bin --input $D/toy_screen_input.tsv --out-dir $O/imputed
$P project --ref $D/toy_reference.tsv --input $D/toy_screen_input.tsv --programs 8 --out $O/projected.tsv
# internal reference: the perturbations of the reference with 150 cells or more, and profiles shrunk by number of cells
I="--ref $D/toy_reference.tsv --min-cells 150 --ref-cells $D/toy_reference_cells.tsv"
$P select $I --panel-size 40 --out $O/panel_internal.tsv
$P fit-reference $I --panel $O/panel_internal.tsv --out $O/model_internal.bin
$P impute --model $O/model_internal.bin --input $D/toy_screen_input.tsv --out-dir $O/imputed_internal
$P project $I --input $D/toy_screen_input.tsv --programs 8 --shrink-by-cells --cells $D/toy_screen_cells.tsv --out $O/shrunk.tsv
for f in imputed_internal/imputed_zscores.tsv shrunk.tsv shrunk.tsv.weights.tsv imputed/imputed_zscores.tsv imputed/imputed_lfsr.tsv imputed/gene_info.tsv imputed/perturbation_qc.tsv projected.tsv; do test -s "$O/$f" || { echo "missing $O/$f"; exit 1; }; done
echo "smoke test passed: outputs in $O"
