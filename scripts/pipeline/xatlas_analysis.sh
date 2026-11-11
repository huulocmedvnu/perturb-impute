#!/bin/sh
# Experiments 1 and 2 after ingest (runs sequentially to stay within memory)
cd "$(dirname "$0")/../.."
B=target/release/panelimp; W=orion/work; L="--l1000 ext/l1000_landmarks.txt"
for c in hct116 hek293t; do $B mp --z $W/$c $L --out $W/${c}_mp > $W/${c}_mp.log 2>&1; done
for c in hct116 hek293t; do
  /usr/bin/time -l $B folds --z $W/$c --corum ext/corum_human.txt $L --out $W/${c}_fx > $W/${c}_folds.log 2>&1
  /usr/bin/time -l $B cv --z $W/$c --folds $W/${c}_fx.folds.tsv --hubs $W/${c}_fx.hub_excl.tsv $L --out $W/${c}_cv --sizes 500 --seeds 1 --baselines > $W/${c}_cv.log 2>&1
  /usr/bin/time -l $B cv --z $W/${c}_high_sgrna --folds $W/${c}_fx.folds.tsv --hubs $W/${c}_fx.hub_excl.tsv $L --out $W/${c}_high_cv --sizes 500 --seeds 0 --baselines > $W/${c}_high_cv.log 2>&1
  for f in 0 1 2 3 4; do $B de --z $W/$c --folds $W/${c}_fx.folds.tsv --hubs $W/${c}_fx.hub_excl.tsv $L --out $W/${c}_de_f$f --fold $f > $W/${c}_de_f$f.log 2>&1; done
  head -1 $W/${c}_de_f0.dex.tsv > $W/${c}_de_all.dex.tsv; for f in 0 1 2 3 4; do tail -n +2 $W/${c}_de_f$f.dex.tsv >> $W/${c}_de_all.dex.tsv; done
done
# common gene set for the multi-atlas experiment
tail -n +2 out/gwps.genes.tsv | cut -f1 > $W/g0.txt
for f in $W/hct116.genes.tsv $W/hek293t.genes.tsv out/hepg2.genes.tsv; do tail -n +2 $f | cut -f1 > $W/g1.txt; awk 'NR==FNR{a[$1]=1;next} ($1 in a)' $W/g1.txt $W/g0.txt > $W/g2.txt; mv $W/g2.txt $W/g0.txt; done
mv $W/g0.txt $W/common_genes.txt
for x in "out/gwps m_k562" "$W/hct116 m_hct116" "$W/hek293t m_hek293t" "out/rpe1 m_rpe1" "out/hepg2 m_hepg2" "out/jurkat m_jurkat"; do set -- $x; $B subset --z $1 --genes $W/common_genes.txt --out $W/$2 >> $W/subset.log 2>&1; done
/usr/bin/time -l $B multi --refs K562=$W/m_k562,HCT116=$W/m_hct116,HEK293T=$W/m_hek293t --targets RPE1=$W/m_rpe1,HepG2=$W/m_hepg2,Jurkat=$W/m_jurkat,K562=$W/m_k562,HCT116=$W/m_hct116,HEK293T=$W/m_hek293t --out $W/multi > $W/multi.log 2>&1
echo done > $W/analysis.done
