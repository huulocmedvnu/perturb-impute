"""Summary: designed vs baselines at each size (common genes, responding perturbations, paired)."""
import sys, re, numpy as np, pandas as pd
from scipy.stats import wilcoxon
R = pd.read_csv(sys.argv[1] + ".perts.tsv", sep="\t"); D = pd.read_csv(sys.argv[1] + ".de.tsv", sep="\t")
R["key"] = R.fold.astype(str) + ":" + R.pert; D["key"] = D.fold.astype(str) + ":" + D.pert
base = R[R.panel == R.panel.iloc[0]].set_index("key"); resp = base.index[base.r_ceil >= 0.2]
den = base.loc[resp, "r_den160_perp"]
W = R.pivot(index="key", columns="panel", values="r_imp_perp").loc[resp]
print(f"responding perturbations {len(resp)}, denoised full (K=160) specific r {den.median():.3f}")
sizes = sorted({int(re.search(r"_(\d+)$", c).group(1)) for c in W.columns})
for k in sizes:
    cols = [c for c in W.columns if c.endswith(f"_{k}")]
    dz = f"designed_{k}"
    print(f"size {k}:")
    for c in sorted(cols, key=lambda c: -np.median(W[c] / den)):
        g = D[D.panel == c].sum(numeric_only=True)
        line = f"  {c:16s} ratio {np.median(W[c] / den):.3f}  matched DE sens {g.k_tp / g.n_de:.3f}"
        if c != dz and dz in W:
            d = (W[dz] - W[c]).dropna(); line += f" | designed - this {d.median():+.4f} wins {np.mean(d > 0):.2f} p {wilcoxon(d).pvalue:.1e}"
        print(line)
