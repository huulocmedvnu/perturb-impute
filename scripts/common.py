"""Common-gene, common-perturbation comparison of panels (paired by perturbation)."""
import sys, numpy as np, pandas as pd
from scipy.stats import wilcoxon
R = pd.read_csv(sys.argv[1] + ".perts.tsv", sep="\t")
dk = [c for c in R.columns if c.startswith("r_den")]
R["key"] = R.fold.astype(str) + ":" + R.pert
base = R[R.panel == R.panel.iloc[0]].set_index("key")
resp = base.index[base.r_ceil >= 0.2]
best = max(dk, key=lambda c: base.loc[resp, c].median())
print(f"responding (same set for all panels): {len(resp)}; raw full {base.loc[resp, 'r_ceil_perp'].median():.3f}; "
      f"denoised full ({best}) {base.loc[resp, best].median():.3f}")
W = R.pivot(index="key", columns="panel", values="r_imp_perp").loc[resp]
for c in W.columns:
    print(f"  {c:13s} specific r {W[c].median():.3f}  / denoised full {np.median(W[c] / base.loc[resp, best]):.3f}  "
          f"share >= 0.8x {np.mean(W[c] >= 0.8 * base.loc[resp, best]):.2f}  AUC DE {R[R.panel == c].set_index('key').loc[resp, 'auc_imp'].median():.3f}")
lm = [c for c in W.columns if c.startswith("l1000")][0]
for c in W.columns:
    if c.startswith("designed"):
        d = W[c] - W[lm]; print(f"  {c} - {lm}: median diff {d.median():+.3f}, wins {np.mean(d > 0):.2f}, Wilcoxon p {wilcoxon(d).pvalue:.1e}")
