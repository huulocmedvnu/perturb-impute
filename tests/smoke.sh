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
for f in imputed/imputed_zscores.tsv imputed/imputed_lfsr.tsv imputed/gene_info.tsv imputed/perturbation_qc.tsv projected.tsv; do test -s "$O/$f" || { echo "missing $O/$f"; exit 1; }; done
echo "smoke test passed: outputs in $O"
