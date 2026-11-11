"""Collect the comparison of alternative imputation models and the guide-target prior into source_data/alt
(Additional file 1: Table S1, Note S1). Inputs: outputs of `panelimp explore --stage test` (five K562 folds and
three transfer screens), of `panelimp explore` (diagnostic, fold 0) and of `panelimp causal` (gene graphs)."""
import os, re, pandas as pd
O = "out/causal/"; T = "source_data/alt/"; os.makedirs(T, exist_ok=True)
sets = [(f"K562 fold {i}", f"k562_f{i}", f"explore_f{i}.txt") for i in range(5)] + [(n, f"tr_{c}", f"explore_{c}.txt") for n, c in (("RPE1", "rpe1"), ("HepG2", "hepg2"), ("Jurkat", "jurkat"))]
bench, paired = [], []
for name, pre, log in sets:
    bench.append(pd.read_csv(f"{O}{pre}.explore.test.tsv", sep="\t").assign(dataset=name))
    for l in open(O + log):
        m = re.match(r"(.{30}) size rank (\d+) (\S+) (\S+) (\S+) \| specific r (\d+) (\S+) (\S+) (\S+)", l)
        if m: paired.append((name, m.group(1).strip(), int(m.group(2)), float(m.group(3)), float(m.group(4)), float(m.group(5)), int(m.group(6)), float(m.group(7)), float(m.group(8)), float(m.group(9))))
b = pd.concat(bench); b = b[["dataset"] + [c for c in b.columns if c != "dataset"]]
b.to_csv(T + "benchmark.tsv", sep="\t", index=False)
pd.DataFrame(paired, columns=["dataset", "method", "n_size_rank", "size_rank_mean_difference", "size_rank_share_higher", "size_rank_p", "n_responding", "specific_r_median_difference", "specific_r_share_higher", "specific_r_p"]).to_csv(T + "paired_against_ridge.tsv", sep="\t", index=False)
g = pd.concat([pd.read_csv(f"{O}k562.causal.fold{i}.tsv", sep="\t").assign(dataset=f"K562 fold {i}") for i in range(5)])
mean_rank = []
for i in range(5):
    for l in open(f"{O}run{i}.txt"):
        t = l.split()
        if l.startswith(("ridge", "direct")) and "size rank" not in l and len(t) >= 9: mean_rank.append((f"K562 fold {i}", l[:28].strip(), float(t[-2])))
g = g.merge(pd.DataFrame(mean_rank, columns=["dataset", "method", "size_rank_mean"]), on=["dataset", "method"]).rename(columns={"size_rank": "size_rank_median"})
g[["dataset", "method", "n_responding", "recovery", "specific_r", "de_recall", "opposite_sign", "size_rank_median", "size_rank_mean"]].to_csv(T + "gene_graphs.tsv", sep="\t", index=False)
open(T + "diagnostic_fold0.txt", "w").write(open(O + "k562.explore.diag.log").read())
print(b.groupby("dataset").size().to_dict(), len(paired), "paired rows", len(g), "graph rows")
