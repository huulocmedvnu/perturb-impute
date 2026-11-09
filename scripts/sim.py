"""Step 6 summary: targeted panel vs whole transcriptome at equal UMI budget."""
import sys, pandas as pd, numpy as np
from scipy.stats import wilcoxon
d = pd.read_csv(sys.argv[1] + ".sim.tsv", sep="\t")
dirs = d[(d.panel == "designed_500") & (d.assay == "extract")]
print("direct whole transcriptome (median r, specific r, matched DE sens):")
for q, g in dirs.groupby("q"):
    print(f"  q {q:<5} r {g.r_dir.median():.3f} rp {g.rp_dir.median():.3f} sens {g.tp_dir.sum() / max(g.n_de.sum(), 1):.3f} disc {g.disc_dir.sum() / max(g.n_de.sum(), 1):.3f}")
for p in ["designed_200", "designed_500", "designed_679", "random1_500", "random1_679", "l1000_679"]:
    for a in ["extract", "targeted"]:
        g0 = d[(d.panel == p) & (d.assay == a)]
        if not len(g0): continue
        print(f"{p} {a}: share {g0.share.iloc[0]:.3f}")
        for q, g in g0.groupby("q"):
            x = (g.rp_imp - g.rp_dir).dropna()
            print(f"  q {q:<5} rate {g.rate.iloc[0]:.3f} r {g.r_imp.median():.3f} rp {g.rp_imp.median():.3f} sens {g.tp_imp.sum() / max(g.n_de.sum(), 1):.3f} disc {g.disc_imp.sum() / max(g.n_de.sum(), 1):.3f}"
                  f" | rp - direct {x.median():+.3f} wins {(x > 0).mean():.2f} p {wilcoxon(x).pvalue:.1e}")
