#!/bin/sh
# X-Atlas/Orion derived models (CC BY-NC-SA 4.0), built after the L1000 runs
cd "$(dirname "$0")/../.."
W=orion/work; D=deliverables_xatlas; P=target/release/perturb-impute
while [ ! -f $W/l1000.done ]; do sleep 30; done
for c in hct116 hek293t; do
  cut -f1,2 $W/${c}_fx.folds.tsv | tail -n +2 > $W/${c}_groups.tsv
  /usr/bin/time $P select --ref $W/$c --panel-size 500 --mode standard --out $D/panels/${c}_500_standard.tsv > $W/${c}_select.log 2>&1
  /usr/bin/time $P fit-reference --ref $W/$c --panel $D/panels/${c}_500_standard.tsv --groups $W/${c}_groups.tsv --out $D/models/${c}_500_standard.bin > $W/${c}_fit.log 2>&1
done
# pan-atlas reference on the common genes, equal weight per atlas
target/release/panelimp concat --in K562=$W/m_k562,HCT116=$W/m_hct116,HEK293T=$W/m_hek293t --out $W/pan3 > $W/pan3.log 2>&1
( cut -f1,2 out/folds.tsv 2>/dev/null | tail -n +2 | sed 's/^/K562:/; s/\t/\tK562_/'
  sed 's/^/HCT116:/; s/\t/\tHCT116_/' $W/hct116_groups.tsv
  sed 's/^/HEK293T:/; s/\t/\tHEK293T_/' $W/hek293t_groups.tsv ) > $W/pan3_groups.tsv
/usr/bin/time $P select --ref $W/pan3 --panel-size 500 --mode standard --out $D/panels/pan_atlas_500_standard.tsv > $W/pan3_select.log 2>&1
/usr/bin/time $P fit-reference --ref $W/pan3 --panel $D/panels/pan_atlas_500_standard.tsv --groups $W/pan3_groups.tsv --out $D/models/pan_atlas_500_standard.bin > $W/pan3_fit.log 2>&1
echo done > $W/models.done
