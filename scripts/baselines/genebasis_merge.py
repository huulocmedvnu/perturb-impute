"""Merge geneBasisR outputs (rank, gene; greedy order) into genebasis.tsv: first 200 / 500 / 680 by rank."""
import os, sys, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import OUT, ROOT
univ = set(pd.read_csv(f"{ROOT}/out/gwps.genes.tsv", sep="\t")["gene"].astype(str))
rows = []
for f in [0, 1, 2, 3, 4, -9]:
    r = pd.read_csv(f"{OUT}/genebasis_parts/fold{f}.tsv", sep="\t", dtype={"gene": str})
    assert list(r["rank"]) == list(range(1, len(r) + 1)), "rank column is not 1..n in order"
    assert len(r) == 680 and r.gene.is_unique and set(r.gene) <= univ, f
    for k in [200, 500, 680]:
        rows += [(f, k, g) for g in r.sort_values("rank").gene.iloc[:k]]
d = pd.DataFrame(rows, columns=["fold", "size", "gene"])
print(d.groupby(["fold", "size"]).size())
d.to_csv(f"{OUT}/genebasis.tsv", sep="\t", index=False)
