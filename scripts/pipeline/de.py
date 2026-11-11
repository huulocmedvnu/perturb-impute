"""Summarize DE call / direction / relative-size preservation (step 3)."""
import sys, pandas as pd, numpy as np
d = pd.read_csv(sys.argv[1] + ".de.tsv", sep="\t")
for p in ["designed_200", "designed_500", "designed_680", "random1_680", "l1000_680"]:
    g = d[d.panel == p]
    if not len(g): continue
    s = g.sum(numeric_only=True); nde = s.n_de
    print(f"{p}: perts {len(g)}, truth DE {int(nde)}")
    for k, lab in [("m", "direct"), ("r", "raw imp"), ("c", "calib")]:
        nc, tp, di = s[k + "_call"], s[k + "_tp"], s[k + "_disc"]
        ex = f" exp {s.c_exp_disc / max(nc, 1):.3f}" if k == "c" else ""
        print(f"  {lab:8s} calls {int(nc):7d} sens {tp / nde:.3f} prec {tp / max(nc, 1):.3f} disc {di / max(nc, 1):.3f}{ex}"
              f" slope {s[k + '_sxy'] / s[k + '_sxx']:.3f} sp {g['sp_' + k].median():.3f} spa {g['spa_' + k].median():.3f} (n {g['sp_' + k].notna().sum()})")
