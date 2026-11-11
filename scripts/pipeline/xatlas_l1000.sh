#!/bin/sh
cd "$(dirname "$0")/../.."
B=target/release/panelimp; W=orion/work; L="--l1000 ext/l1000_landmarks.txt"
for c in hct116 hek293t; do
  /usr/bin/time -l $B cv --z $W/$c --folds $W/${c}_fx.folds.tsv --hubs $W/${c}_fx.hub_excl.tsv $L --out $W/${c}_lm --sizes lm --seeds 1 > $W/${c}_lm.log 2>&1
done
echo done > $W/l1000.done
