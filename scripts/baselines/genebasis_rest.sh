#!/bin/bash
# Continuation queue used in the actual run: fold -9 immediately, folds 3 and 4 once folds 1 / 2 finish; then merge.
cd ~/Downloads/panelimp/ext/baselines
run() { /opt/homebrew/bin/Rscript ../../scripts/baselines/genebasis_run.R $1 680 genebasis_parts/fold$1.tsv 2>&1 | while IFS= read -r l; do echo "$(date +%H:%M:%S) $l"; done > logs/genebasis_fold$1.log; }
run -9 &
( until grep -q "gene_search done" logs/genebasis_fold1.log; do sleep 15; done; run 3 ) &
( until grep -q "gene_search done" logs/genebasis_fold2.log; do sleep 15; done; run 4 ) &
wait
python3 ../../scripts/baselines/genebasis_merge.py
