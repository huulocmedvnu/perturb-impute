"""Summary: budget-aware panels vs designed / L1000 / random, targeted assay."""
import sys, pandas as pd, numpy as np
from scipy.stats import wilcoxon
d = pd.read_csv(sys.argv[1] + ".sim.tsv", sep="\t"); d = d[d.assay == "targeted"]
comp = ["designed_500", "greedyfull_500", "l1000_679", "random1_500"]
for q, g in d.groupby("q"):
    print(f"q {q}: direct WT rp {g[g.panel == 'designed_500'].rp_dir.median():.3f}")
    names = comp + sorted(p for p in g.panel.unique() if p.startswith("budget"))
    for p in names:
        h = g[g.panel == p]
        print(f"  {p:16s} share {h.share.iloc[0]:.3f} rate {h.rate.iloc[0]:.3f} rp {h.rp_imp.median():.4f} r {h.r_imp.median():.3f} sens {h.tp_imp.sum() / h.n_de.sum():.3f}")
    own = f"budget{q}_500"
    if own in set(g.panel):
        b = g[g.panel == own].set_index(["fold", "pert"]).rp_imp
        for p in comp:
            x = (b - g[g.panel == p].set_index(["fold", "pert"]).rp_imp).dropna()
            print(f"    {own} - {p}: median {x.median():+.4f} wins {(x > 0).mean():.2f} p {wilcoxon(x).pvalue:.1e}")
