"""All manuscript numbers from the final outputs (greedy MMSE = main selector) -> docs/manuscript/numbers.json."""
import os
import json, numpy as np, pandas as pd
from scipy.stats import wilcoxon, spearmanr
from sklearn.metrics import roc_auc_score
O = "./out/"; N = {}
CT = {"K562": "cvg", "RPE1": "trg", "HepG2": "tr_hepg2", "Jurkat": "tr_jurkat"}
SEL = ["greedy", "hvgresp", "designed", "genebasis", "hvg", "module", "random1"]
def load(pre):
    R = pd.read_csv(O + pre + ".perts.tsv", sep="\t"); R["key"] = R.fold.astype(str) + ":" + R.pert
    D = pd.read_csv(O + pre + ".de.tsv", sep="\t"); D["key"] = D.fold.astype(str) + ":" + D.pert
    b = R[R.panel == R.panel.iloc[0]].set_index("key"); resp = b.index[b.r_ceil >= 0.2]
    return R, D, b, resp
for ct, pre in CT.items():
    R, D, b, resp = load(pre); den = b.loc[resp, "r_den160_perp"]
    W = R.pivot(index="key", columns="panel", values="r_imp_perp").loc[resp]
    lm = [c for c in W if c.startswith("l1000")][0]; klm = lm.split("_")[1]
    n = {"n_perts": int(b.shape[0]), "n_resp": len(resp), "den160": float(den.median()), "raw_full": float(b.loc[resp, "r_ceil_perp"].median()),
         "n_common_genes": None, "l1000_size": int(klm)}
    for k in sorted({int(c.split("_")[-1]) for c in W}):
        for s in SEL + ["l1000"]:
            c = f"{s}_{k}"
            if c in W: n[f"ratio_{c}"] = float((W[c] / den).median()); n[f"spec_{c}"] = float(W[c].median())
            if c in W and s != "greedy" and f"greedy_{k}" in W:
                d = (W[f"greedy_{k}"] - W[c]).dropna(); n[f"diff_greedy_{c}"] = [float(d.median()), float((d > 0).mean()), float(wilcoxon(d).pvalue)]
    n["share_ge08_greedy_500"] = float((W["greedy_500"] >= 0.8 * den).mean())
    # absolute per-perturbation Pearson r on all unmeasured genes, responding perturbations
    A = R.pivot(index="key", columns="panel", values="r_imp").loc[resp]; n["r_greedy_500"] = float(A["greedy_500"].median()); n["r_direct"] = float(b.loc[resp, "r_ceil"].median())
    # OOR: failure = specific ratio to raw half-0 measurement < 0.5
    g = R[(R.panel == "greedy_500")].set_index("key").loc[resp].dropna(subset=["oor", "r_imp_perp", "r_ceil_perp"])
    fail = (g.r_imp_perp / g.r_den160_perp < 0.5).astype(int)
    n["oor_median"] = float(g.oor.median()); n["oor_fail_n"] = int(fail.sum()); n["oor_n"] = len(fail)
    n["oor_auc"] = float(roc_auc_score(fail, g.oor)) if 0 < fail.sum() < len(fail) else None
    n["oor_spearman_ratio"] = list(map(float, spearmanr(g.oor, g.r_imp_perp / g.r_den160_perp)))
    # DE preservation for greedy_500 (all test perturbations)
    s = D[D.panel == "greedy_500"].sum(numeric_only=True)
    n["de"] = {"truth_de": int(s.n_de), "direct_sens": float(s.m_tp / s.n_de), "direct_disc": float(s.m_disc / s.m_call), "matched_sens": float(s.k_tp / s.n_de),
               "matched_disc": float(s.k_disc / s.m_call), "slope_cal": float(s.c_sxy / s.c_sxx), "slope_raw": float(s.r_sxy / s.r_sxx),
               "lfsr_obs_exp": float(s.c_disc / s.c_exp_disc), "lfsr_calls": int(s.c_call),
               "sp_signed_imp": float(D[D.panel == "greedy_500"].sp_c.median()), "sp_signed_dir": float(D[D.panel == "greedy_500"].sp_m.median()),
               "sp_abs_imp": float(D[D.panel == "greedy_500"].spa_r.median()), "sp_abs_dir": float(D[D.panel == "greedy_500"].spa_m.median()),
               "n_sp": int(D[D.panel == "greedy_500"].sp_m.notna().sum())}
    for o in ["l1000_" + klm, "designed_500", "random1_500"]:
        so = D[D.panel == o].sum(numeric_only=True); n["de"][f"matched_sens_{o}"] = float(so.k_tp / so.n_de)
    # INFO calibration (greedy_500)
    G = pd.read_csv(O + pre + ".genes.tsv", sep="\t"); G = G[G.panel == "greedy_500"].dropna()
    n["info_slope"] = float(np.polyfit(G["info"], G["real"], 1)[0]); n["info_r"] = float(np.corrcoef(G["info"], G["real"])[0, 1])
    n["info_median"] = float(G["info"].median()); n["real_median"] = float(G["real"].median()); n["n_common_genes"] = int(G.gene.nunique())
    N[ct] = n
# main characterisation of greedy_500 on all unmeasured genes (non-common runs)
CTM = {"K562": "cvm", "RPE1": "trm_rpe1", "HepG2": "trm_hepg2", "Jurkat": "trm_jurkat"}
for ct, pre in CTM.items():
    R, D, b, resp = load(pre); g = R[R.panel == "greedy_500"].set_index("key"); resp = g.index[g.r_ceil >= 0.2]; g = g.loc[resp]
    m = {"n_resp": len(resp), "ratio_den": float((g.r_imp_perp / g.r_den160_perp).median()), "spec_imp": float(g.r_imp_perp.median()), "den160": float(g.r_den160_perp.median()),
         "raw_half": float(g.r_ceil_perp.median()), "r_imp": float(g.r_imp.median()), "r_direct": float(g.r_ceil.median()),
         "share_ge08": float((g.r_imp_perp >= 0.8 * g.r_den160_perp).mean()), "auc_imp": float(g.auc_imp.median()), "auc_dir": float(g.auc_ceil.median())}
    gg = g.dropna(subset=["oor", "r_imp_perp", "r_ceil_perp"]); fail = (gg.r_imp_perp / gg.r_ceil_perp < 0.5).astype(int)
    m["oor_median"] = float(gg.oor.median()); m["oor_fail_n"] = int(fail.sum()); m["oor_auc"] = float(roc_auc_score(fail, gg.oor)) if 0 < fail.sum() < len(fail) else None
    m["oor_spearman"] = list(map(float, spearmanr(gg.oor, gg.r_imp_perp / gg.r_den160_perp)))
    s_ = D[D.panel == "greedy_500"]; t = s_.sum(numeric_only=True)
    m["de"] = {"truth_de": int(t.n_de), "direct_sens": float(t.m_tp / t.n_de), "direct_disc": float(t.m_disc / t.m_call), "matched_sens": float(t.k_tp / t.n_de),
               "matched_disc": float(t.k_disc / t.m_call), "slope_cal": float(t.c_sxy / t.c_sxx), "slope_raw": float(t.r_sxy / t.r_sxx),
               "lfsr_obs_exp": float(t.c_disc / t.c_exp_disc), "sp_signed_imp": float(s_.sp_c.median()), "sp_signed_dir": float(s_.sp_m.median()),
               "sp_abs_imp": float(s_.spa_r.median()), "sp_abs_dir": float(s_.spa_m.median()), "n_sp": int(s_.sp_m.notna().sum())}
    G = pd.read_csv(O + pre + ".genes.tsv", sep="\t"); G = G[G.panel == "greedy_500"].dropna()
    m["info_slope"] = float(np.polyfit(G["info"], G["real"], 1)[0]); m["info_r"] = float(np.corrcoef(G["info"], G["real"])[0, 1]); m["n_genes_eval"] = int(G.gene.nunique())
    q = pd.qcut(G["info"], 5, labels=False); m["info_bins"] = [[float(G["info"][q == i].median()), float(G["real"][q == i].median())] for i in range(5)]
    N.setdefault("main", {})[ct] = m
    print("MAIN", ct, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in m.items() if k != "de" and k != "info_bins"}); print("   DE", {k: round(v, 3) for k, v in m["de"].items()})
# K562 OOR by cell type summary
N["oor_by_ct"] = {ct: N[ct]["oor_median"] for ct in CT}
# equal-budget simulation: out/simu (fold 0 = development for the UMI-cost panel, report folds 1-4)
d = pd.read_csv(O + "simu.sim.tsv", sep="\t"); d = d[d.fold.isin([1, 2, 3, 4])]; S = {}
for q, g in d.groupby("q"):
    w = g[(g.panel == "l1000_679") & (g.assay == "extract")]
    S[str(q)] = {"direct_rp": float(w.rp_dir.median()), "direct_r": float(w.r_dir.median()), "direct_sens": float(w.tp_dir.sum() / w.n_de.sum())}
    ref = g[(g.panel == "l1000_679") & (g.assay == "targeted")].set_index(["fold", "pert"]).rp_imp
    for p, a in [("greedyfull_500", "targeted"), ("greedyfull_500", "targeted_het"), ("cost0.5_500", "targeted_het"), ("cost0.5_500", "targeted"),
                 ("l1000_679", "targeted"), ("random1_500", "targeted"), ("designed_500", "targeted"), ("cost1_500", "targeted_het"), ("cap0.12_500", "targeted_het"), ("greedyfull_500", "extract")]:
        h = g[(g.panel == p) & (g.assay == a)]
        x = (h.rp_imp - h.rp_dir).dropna(); y = (h.set_index(["fold", "pert"]).rp_imp - ref).dropna()
        S[str(q)][f"{p}|{a}"] = {"rp": float(h.rp_imp.median()), "r": float(h.r_imp.median()), "sens": float(h.tp_imp.sum() / h.n_de.sum()), "share": float(h.share.median()),
                                 "rate": float(h.rate.median()), "vs_direct": [float(x.median()), float((x > 0).mean()), float(wilcoxon(x).pvalue)],
                                 "vs_l1000": [float(y.median()), float((y > 0).mean()), float(wilcoxon(y).pvalue) if (y != 0).sum() > 10 else None]}
S["n_query"] = int(d.pert.nunique()); S["n_cells_total"] = 310385; N["sim"] = S
# budget-aware MMSE panels (out/simb, all folds), Additional file Fig. S1
d = pd.read_csv(O + "simb.sim.tsv", sep="\t"); SB = {}
for q, g in d.groupby("q"):
    SB[str(q)] = {p: float(g[(g.panel == p) & (g.assay == "targeted")].rp_imp.median()) for p in g.panel.unique()}
N["sim_budget"] = SB
# DE call variants and hybrid readout (out/de_*.dex.tsv; K562 = five folds)
H = {}
for ct, f in [("K562", "de_k562_all"), ("RPE1", "de_rpe1"), ("HepG2", "de_hepg2"), ("Jurkat", "de_jurkat")]:
    d = pd.read_csv(O + f + ".dex.tsv", sep="\t"); H[ct] = {}
    for m, g in d.groupby("method"):
        H[ct][m] = {"recall": float(g.tp.sum() / g.n_de.sum()), "opp": float(g.disc.sum() / g.ncall.sum()), "rank": float(g.spa.median()), "n_rank": int(g.spa.notna().sum())}
N["hybrid"] = H
# cv2 decomposition
cv = pd.read_csv(O + "simb.cv2.tsv", sep="\t").dropna(); N["cv2_tech_dominated"] = float((cv.cv2_tech > cv.cv2_bio).mean()); N["cv2_tech_share_median"] = float((cv.cv2_tech / cv.cv2_tot).median())
# effect-size ranking: projection control (first fold, scripts/pipeline/oracle_rank.py) and neighbour copying (cvk)
N["rank_oracle"] = {"direct": 0.311, "K20": 0.200, "K160": 0.209, "K500": 0.252, "K2000": 0.286, "ridge500": 0.191}
k = pd.read_csv(O + "cvk.knn.tsv", sep="\t"); k = k[k.panel == "designed_500"]
N["rank_knn"] = {"ridge": float(k.ridge_spa.median()), "knn": float(k.knn_spa.median()), "hyb": float(k.hyb_spa.median()), "p_hyb": float(wilcoxon((k.hyb_spa - k.ridge_spa).dropna()).pvalue)}
# noise tails
m = pd.read_csv(O + "gwps.meta.tsv", sep="\t"); gg = pd.read_csv(O + "gwps.genes.tsv", sep="\t")
Z = np.fromfile(O + "gwps.z.f32", np.float32).reshape(3, len(m), len(gg)); dd = (Z[1] - Z[2]) / np.sqrt(2)
N["noise"] = {"sd": float(dd.std()), "tail4": float((np.abs(dd) > 4).mean()), "normal_tail4": 6.334e-5, "n_genes_universe": len(gg), "n_perts_gwps": len(m)}
NJ = "./source_data/numbers.json"
N["mp"] = json.load(open(NJ)).get("mp", {}) if os.path.exists(NJ) else {}
json.dump(N, open(NJ, "w"), indent=1)
for ct in CT: print(ct, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in N[ct].items() if k.startswith(("n_", "ratio_greedy_500", "ratio_hvgresp_500", "ratio_designed_500", "ratio_l1000", "oor", "info", "r_", "share"))})
