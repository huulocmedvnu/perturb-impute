"""Collect the biological-inference analysis (panelimp bio) into source_data/bio/summary.tsv:
five K562 folds and three screens of other cell lines."""
import pandas as pd
O = "out/bio/"; rows = []
for i in range(5): rows.append(pd.read_csv(f"{O}k562.bio.fold{i}.tsv", sep="\t").assign(dataset=f"K562 fold {i}"))
for n, c in (("RPE1", "rpe1"), ("HepG2", "hepg2"), ("Jurkat", "jurkat")): rows.append(pd.read_csv(f"{O}{c}.bio.tsv", sep="\t").assign(dataset=n))
b = pd.concat(rows); b = b[["dataset", "analysis", "frame", "readout", "metric", "n", "value"]]
b.to_csv("source_data/bio/summary.tsv", sep="\t", index=False)
print(len(b), "rows")
