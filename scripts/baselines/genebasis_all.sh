#!/bin/bash
# Runs geneBasisR gene_search(n_genes_total = 680) for folds 0-4 and -9 on 3,000 subsampled reference profiles.
# Usage: genebasis_all.sh [PARALLEL=2] [FOLDS...]
cd "$(dirname "$0")/../.."/ext/baselines
P=${1:-2}; shift; FOLDS=${@:-"0 1 2 3 4 -9"}
python3 ../../scripts/baselines/genebasis_export.py --sub 3000
mkdir -p genebasis_parts logs
for f in $FOLDS; do echo $f; done | xargs -P $P -I{} bash -c \
  '/opt/homebrew/bin/Rscript ../../scripts/baselines/genebasis_run.R {} 680 genebasis_parts/fold{}.tsv 2>&1 | while IFS= read -r l; do echo "$(date +%H:%M:%S) $l"; done > logs/genebasis_fold{}.log'
python3 ../../scripts/baselines/genebasis_merge.py
