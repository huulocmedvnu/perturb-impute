#!/bin/sh
cd ~/Downloads/panelimp/orion
while pgrep -f "orion-ingest run --cell HCT116" >/dev/null; do sleep 30; done
../target/release/orion-ingest finish --work work --cell HCT116 --genes meta/HCT116_genes.txt --out work/hct116 > work/HCT116_finish.log 2>&1
../target/release/orion-ingest run --cell HEK293T --shards meta/HEK293T_shards.txt --meta meta/gene_metadata.parquet --guides meta/HEK293T_guide_calls.csv.gz --genes meta/HEK293T_genes.txt --work work > work/HEK293T.log 2>&1
../target/release/orion-ingest finish --work work --cell HEK293T --genes meta/HEK293T_genes.txt --out work/hek293t > work/HEK293T_finish.log 2>&1
echo done > work/chain.done
