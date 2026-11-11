"""Summary: ridge vs knn vs hybrid (size ranking among DE, matched calls, r)."""
import sys, pandas as pd, numpy as np
from scipy.stats import wilcoxon
d = pd.read_csv(sys.argv[1] + ".knn.tsv", sep="\t")
for p in ["designed_500", "l1000_680", "random1_680"]:
    g = d[d.panel == p]
    if not len(g): continue
    print(p, "params knn", sorted(set(zip(g.k_knn, g.p_knn))), "hyb", sorted(set(zip(g.k_hyb, g.p_hyb))))
    nde, nm = g.n_de.sum(), g.m_call.sum()
    for m in ["ridge", "knn", "hyb"]:
        print(f"  {m:5s} spa {g[m+'_spa'].median():.3f} sp {g[m+'_sp'].median():.3f} r {g[m+'_r'].median():.3f} rperp {g[m+'_rperp'].median():.3f} sens {g[m+'_ktp'].sum()/nde:.3f} disc {g[m+'_kdisc'].sum()/nm:.3f}")
    for m in ["knn", "hyb"]:
        x = (g[m + "_spa"] - g.ridge_spa).dropna(); y = (g[m + "_r"] - g.ridge_r).dropna()
        print(f"  {m} - ridge: spa median diff {x.median():+.3f} wins {(x>0).mean():.2f} p {wilcoxon(x).pvalue:.1e} | r diff {y.median():+.3f}")
