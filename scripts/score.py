"""Score K1-K3 from panelimp output (criteria fixed before the analysis)."""
import sys, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
pre = sys.argv[1] if len(sys.argv) > 1 else "./out/cv"
R = pd.read_csv(pre + ".perts.tsv", sep="\t"); Gt = pd.read_csv(pre + ".genes.tsv", sep="\t")
R["kind"] = R.panel.str.extract(r"^(designed|random|l1000)")[0]; R["size"] = R.panel.str.extract(r"_(\d+)$")[0].astype(int)
R["ratio_perp"] = R.r_imp_perp / R.r_ceil_perp; R["ratio"] = R.r_imp / R.r_ceil
resp = R.r_ceil >= 0.2
print(f"test rows {len(R)}; responding perturbations (r_ceil >= 0.2): {R[resp & (R.panel == R.panel.iloc[0])].pert.nunique()}")
S = R[resp].groupby(["kind", "size"]).agg(total=("ratio", "median"), specific=("ratio_perp", "median"),
        r_imp_perp=("r_imp_perp", "median"), r_ceil_perp=("r_ceil_perp", "median"), auc_imp=("auc_imp", "median"), auc_ceil=("auc_ceil", "median")).round(3)
print(S.to_string())
d500 = S.loc[("designed", 500), "specific"]; dlm = S.loc[("designed", S.loc["l1000"].index[0]), "specific"]; lm = S.loc["l1000"].specific.iloc[0]
print(f"K1 designed 500 specific ratio {d500:.3f} (need >= 0.6): {'PASS' if d500 >= 0.6 else 'FAIL'}")
print(f"K2 designed {S.loc['l1000'].index[0]} {dlm:.3f} vs L1000 {lm:.3f} (diff {dlm - lm:+.3f}); designed 500 vs L1000 {d500 - lm:+.3f}: "
      f"{'PASS' if dlm - lm >= 0.1 or d500 >= lm else 'FAIL'}")
g = Gt[Gt.panel == "designed_500"].dropna()
slope = np.polyfit(g["info"], g["real"], 1)[0]
r5 = R[resp & (R.panel == "designed_500")].dropna(subset=["ratio_perp", "oor"])
fail = (r5.ratio_perp < 0.5).astype(int)
a = roc_auc_score(fail, r5.oor) if 0 < fail.sum() < len(fail) else np.nan
print(f"K3 INFO calibration slope {slope:.3f} (need 0.8-1.2); OOR AUC for failure {a:.3f} (need >= 0.7; failures {fail.sum()}/{len(fail)})")
