"""Collect the per-perturbation result tables behind the figures and the manuscript numbers into
source_data/tables (gzip, floats written with 5 significant digits, correlations with 6 decimals as in the pipeline output). Run from the repository root after the
analysis pipeline has written its outputs to out/."""
import os, pandas as pd
O = "out/"; T = "source_data/tables/"; os.makedirs(T, exist_ok=True)
def put(df, name, fmt="%.5g"): df.to_csv(T + name + ".gz", sep="\t", index=False, float_format=fmt, compression={"method": "gzip", "mtime": 0}); print(f"{name:28s} {len(df):8d} rows")
# selector comparison (common evaluation genes): responding perturbations only
for pre in ["cvg", "trg", "tr_hepg2", "tr_jurkat"]:
    R = pd.read_csv(O + pre + ".perts.tsv", sep="\t"); b = R[R.panel == R.panel.iloc[0]]
    keys = set(zip(b.fold[b.r_ceil >= 0.2], b.pert[b.r_ceil >= 0.2]))
    R = R[[k in keys for k in zip(R.fold, R.pert)]]
    put(R[["fold", "panel", "pert", "r_ceil", "r_imp", "r_ceil_perp", "r_imp_perp", "r_den160_perp"]], pre + ".perts.tsv", "%.6f")
# main characterisation of the greedy MMSE 500 panel (all genes outside the panel)
for pre in ["cvm", "trm_rpe1", "trm_hepg2", "trm_jurkat"]:
    R = pd.read_csv(O + pre + ".perts.tsv", sep="\t"); R = R[R.panel == "greedy_500"]
    put(R[["fold", "panel", "pert", "r_ceil", "r_imp", "r_ceil_perp", "r_imp_perp", "auc_ceil", "auc_imp", "oor", "r_den160_perp"]], pre + ".perts.tsv", "%.6f")
    D = pd.read_csv(O + pre + ".de.tsv", sep="\t"); put(D[D.panel == "greedy_500"], pre + ".de.tsv")
    G = pd.read_csv(O + pre + ".genes.tsv", sep="\t"); put(G[G.panel == "greedy_500"], pre + ".genes.tsv")
for f in ["de_k562_all", "de_rpe1", "de_hepg2", "de_jurkat"]: put(pd.read_csv(O + f + ".dex.tsv", sep="\t"), f + ".dex.tsv")
S = pd.read_csv(O + "simd.sim.tsv", sep="\t")
put(S[["fold", "panel", "assay", "q", "pert", "share", "n_de", "rp_imp", "tp_imp", "disc_imp", "rp_dir", "tp_dir", "disc_dir"]], "simd.sim.tsv")
put(pd.read_csv(O + "simd.den.tsv", sep="\t"), "simd.den.tsv")
B = pd.read_csv(O + "simb.sim.tsv", sep="\t"); put(B[B.assay == "targeted"][["fold", "panel", "assay", "q", "pert", "rp_imp"]], "simb.sim.tsv")
put(pd.read_csv(O + "simb.cv2.tsv", sep="\t"), "simb.cv2.tsv")
K = pd.read_csv(O + "cvk.knn.tsv", sep="\t"); put(K[K.panel == "designed_500"][["fold", "panel", "pert", "ridge_spa", "knn_spa", "hyb_spa"]], "cvk.knn.tsv")
put(pd.read_csv(O + "real/r.real.tsv", sep="\t"), "paired.real.tsv"); put(pd.read_csv(O + "real/r.cells.tsv", sep="\t"), "paired.cells.tsv")
