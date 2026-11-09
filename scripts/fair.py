"""Panel imputation vs full transcriptome denoised on the same reference programs (top-K, best K
chosen on test = optimistic for full sequencing)."""
import sys, numpy as np, pandas as pd
R = pd.read_csv(sys.argv[1] + ".perts.tsv", sep="\t"); R = R[R.r_ceil >= 0.2]
dk = [c for c in R.columns if c.startswith("r_den")]
for pan in ["designed_100", "designed_200", "designed_500", "designed_680", "l1000_680", "random1_500"]:
    t = R[R.panel == pan]
    if t.empty: continue
    med = {c: t[c].median() for c in dk}; best = max(med, key=med.get)
    print(f"{pan:13s} n {len(t)} | specific r: panel {t.r_imp_perp.median():.3f}, raw full {t.r_ceil_perp.median():.3f}, "
          f"denoised full best {best} {med[best]:.3f} | panel/denoised median ratio {np.median(t.r_imp_perp / t[best]):.3f}, "
          f"share >= 0.8x {np.mean(t.r_imp_perp >= 0.8 * t[best]):.2f}")
