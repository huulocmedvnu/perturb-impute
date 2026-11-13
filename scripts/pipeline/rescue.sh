#!/bin/sh
# Perturbations with few cells (Additional file 1: Table S10): external and internal references, shrinkage,
# and the input half of the K562 essential-gene screen restricted to 15, 30 and 50 cells.
# Inputs: z-score prefixes in out/, the cell cache of the K562 essential-gene screen (CACHE), ext/ and models/inputs/.
set -e
cd "$(dirname "$0")/../.."
P=target/release/panelimp; O=out/rescue; CACHE=${CACHE:-../coexsc/data/k562ess}
B="--corum ext/corum_human.txt --genesets ext/genesets/hallmark_2020.txt --tss models/inputs/tss_grch38_refseq_select.tsv --umi models/inputs/control_umi_fraction_k562.tsv"
K="--ref-folds out/folds.tsv --ref-hubs out/hub_excl.tsv"
mkdir -p $O
for s in rpe1 hepg2 jurkat; do
  $P folds --z out/$s --corum ext/corum_human.txt --l1000 ext/l1000_landmarks.txt --out $O/$s
done
$P rescue --ref out/gwps $K $B --out $O/k562
$P rescue --ref out/gwps $K --screen out/rpe1 --folds $O/rpe1.folds.tsv --hubs $O/rpe1.hub_excl.tsv $B --out $O/rpe1
for s in hepg2 jurkat; do
  $P rescue --ref out/gwpsHJ $K --screen out/$s --folds $O/$s.folds.tsv --hubs $O/$s.hub_excl.tsv $B --out $O/$s
done
$P thin --cache $CACHE --z out/gwpsE --out $O/k562ess --cells 15,30,50
$P rescue --ref out/gwpsE $K --screen $O/k562ess --folds out/folds.tsv --hubs out/hub_excl.tsv --same-line \
  --thin $O/k562ess --thin-cells 15,30,50 --min-test-in 50 $B --out $O/k562ess
