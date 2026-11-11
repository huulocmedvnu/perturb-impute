#!/bin/sh
# step 14: real targeted readout (Replogle 2020 exp7 / exp11); run from the repository root
A="--z out/gwps --l1000 ext/l1000_landmarks.txt --folds out/folds.tsv --hubs out/hub_excl.tsv --cache-wt out/real/wt --cache-tg out/real/tg --gene-ids models/inputs/gene_ids_gencode_v44.tsv --out out/real/r"
target/release/panelimp real $A 2> out/real/real.log
target/release/panelimp realsum --in out/real/r > out/real/real.summary.tsv
target/release/panelimp realcells $A 2> out/real/cells.log
target/release/panelimp cellsum --in out/real/r > out/real/cells.summary.tsv
