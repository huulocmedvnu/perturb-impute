"""Recompute the key numbers of the manuscript from source_data/ and compare them with the values
printed in the text. Each check states its source: `tables` (recomputed from per-perturbation tables),
`summary` (read from a summary written by the Rust tool) or `json` (numbers.json, for values whose
inputs are too large to ship: program count and noise tails need the full atlas).
Exit status 1 if any check fails.  Usage: python3 scripts/verify_manuscript_numbers.py"""
import json, os, re, sys
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

SD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "source_data"); T = os.path.join(SD, "tables")
tab = lambda name: pd.read_csv(os.path.join(T, name + ".gz"), sep="\t")
N = json.load(open(os.path.join(SD, "numbers.json")))
rows = []

def check(section, claim, expected, value, source="tables"):
    """expected is the string printed in the manuscript ('0.89', '79%', '763', 'P < 0.001', 'P = 0.40')."""
    if expected.startswith("P <"):
        ok = value < float(expected.split("<")[1]); shown = f"{value:.1e}"
    elif expected.startswith("P ="):
        d = len(expected.split(".")[1]) if "." in expected else 0
        ok = round(value, d) == float(expected.split("=")[1]); shown = f"{value:.4f}"
    elif expected.endswith("%"):
        e = expected[:-1]; d = len(e.split(".")[1]) if "." in e else 0
        ok = abs(round(100 * value, d) - float(e)) < 1e-9 or abs(100 * value - float(e)) <= 0.5 * 10 ** -d + 1e-6; shown = f"{100 * value:.2f}%"
    elif "." in expected:
        d = len(expected.split(".")[1]); ok = abs(value - float(expected)) <= 0.5 * 10 ** -d + 1e-6; shown = f"{value:.4f}"
    else:
        ok = int(round(value)) == int(expected.replace(",", "")); shown = f"{int(round(value)):,}"
    rows.append(("PASS" if ok else "FAIL", section, claim, expected, shown, source))

def in_range(section, claim, lo, hi, values, source="tables"):
    d = len(lo.rstrip("%").split(".")[1]) if "." in lo else 0; pct = lo.endswith("%"); f = 100 if pct else 1
    r = [round(f * v, d) for v in values]; ok = min(r) == float(lo.rstrip("%")) and max(r) == float(hi.rstrip("%"))
    rows.append(("PASS" if ok else "FAIL", section, claim, f"{lo} to {hi}", ", ".join(f"{f * v:.3f}" for v in values), source))

# ---------------- imputation from the greedy MMSE 500 panel (Fig. 3)
def main_run(pre):
    g = tab(pre + ".perts.tsv"); g = g[(g.panel == "greedy_500") & (g.r_ceil >= 0.2)]; return g
k = main_run("cvm"); S = "Imputation, K562"
check(S, "responding test perturbations", "763", len(k))
check(S, "median recovery", "0.89", (k.r_imp_perp / k.r_den160_perp).median())
check(S, "perturbations with recovery of 0.8 or more", "79%", (k.r_imp_perp >= 0.8 * k.r_den160_perp).mean())
check(S, "total correlation, imputed", "0.45", k.r_imp.median()); check(S, "total correlation, measured input half", "0.30", k.r_ceil.median())
in_range("Imputation, other cell lines", "median recovery RPE1, HepG2, Jurkat", "0.75", "0.76",
         [(g.r_imp_perp / g.r_den160_perp).median() for g in map(main_run, ["trm_rpe1", "trm_hepg2", "trm_jurkat"])])
G = tab("cvm.genes.tsv").dropna()
check(S, "INFO slope", "0.85", np.polyfit(G["info"], G["real"], 1)[0]); check(S, "INFO correlation", "0.81", np.corrcoef(G["info"], G["real"])[0, 1]); check(S, "genes evaluated", "6,680", G.gene.nunique())
kk = k.dropna(subset=["oor", "r_imp_perp", "r_ceil_perp"]); fail = kk.r_imp_perp / kk.r_ceil_perp < 0.5
check(S, "median out-of-reference score", "0.14", kk.oor.median()); check(S, "perturbations below half of the measured halves", "9", fail.sum())
r = kk.oor.rank(); n1, n0 = fail.sum(), (~fail).sum(); check(S, "out-of-reference AUC", "0.90", (r[fail].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))

# ---------------- panel selectors (Fig. 2)
def selectors(pre):
    R = tab(pre + ".perts.tsv"); R["key"] = R.fold.astype(str) + ":" + R.pert
    den = R[R.panel == R.panel.iloc[0]].set_index("key").r_den160_perp
    return R.pivot(index="key", columns="panel", values="r_imp_perp"), den
W, den = selectors("cvg"); S = "Selectors, K562"
check(S, "responding test perturbations", "515", len(W))
for name, exp in [("greedy_500", "0.87"), ("hvgresp_500", "0.86"), ("designed_500", "0.85"), ("genebasis_500", "0.83"), ("hvg_500", "0.79"), ("random1_500", "0.60"),
                  ("designed_200", "0.77"), ("greedy_200", "0.76"), ("l1000_680", "0.69")]:
    check(S, f"recovery {name}", exp, (W[name] / den).median())
check(S, "recovery random panel of 680 genes", "0.67", (W["random1_680"] / den).median())
d = (W["greedy_680"] - W["l1000_680"]).dropna(); check(S, "greedy 680 above L1000 680", "91%", (d > 0).mean()); check(S, "greedy 680 against L1000 680", "P < 0.001", wilcoxon(d).pvalue)
lm, gr = [], []
for pre in ["trg", "tr_hepg2", "tr_jurkat"]:
    W2, den2 = selectors(pre); gr.append((W2["greedy_500"] / den2).median()); lm.append((W2[[c for c in W2 if c.startswith("l1000")][0]] / den2).median())
in_range("Selectors, other cell lines", "L1000 recovery", "0.57", "0.61", lm); in_range("Selectors, other cell lines", "greedy MMSE recovery", "0.73", "0.73", gr)

# ---------------- differential expression (Fig. 4)
D = tab("de_k562_all.dex.tsv"); S = "DE, K562"
rec = lambda m: D[D.method == m].tp.sum() / D[D.method == m].n_de.sum(); opp = lambda m: D[D.method == m].disc.sum() / D[D.method == m].ncall.sum()
check(S, "DE pairs outside the panel", "120,209", D[D.method == "direct_out"].n_de.sum())
check(S, "imputed genes, DE found", "21%", rec("baseline")); check(S, "direct outside panel, DE found", "21%", rec("direct_out"))
check(S, "imputed genes, opposite sign", "14%", opp("baseline")); check(S, "direct outside panel, opposite sign", "21%", opp("direct_out"))
check(S, "hybrid profile, DE found", "36%", rec("hybrid_all")); check(S, "direct all genes, DE found", "35%", rec("direct_all"))
check(S, "hybrid profile, opposite sign", "9%", opp("hybrid_all")); check(S, "direct all genes, opposite sign", "15%", opp("direct_all"))
med = lambda m: D[D.method == m].spa.median()
check(S, "size ranking, imputed", "0.10", med("baseline")); check(S, "size ranking, two measured halves", "0.20", med("direct_out"))
check(S, "size ranking, hybrid", "0.24", med("hybrid_all")); check(S, "size ranking, direct all genes", "0.34", med("direct_all"))
E = tab("cvm.de.tsv").sum(numeric_only=True)
check(S, "slope of measured on imputed", "0.80", E.r_sxy / E.r_sxx); check(S, "slope after calibration", "0.98", E.c_sxy / E.c_sxx); check(S, "observed over expected sign errors", "1.09", E.c_disc / E.c_exp_disc)
oth = [tab(f"de_{c}.dex.tsv") for c in ("rpe1", "hepg2", "jurkat")]
in_range("DE, other cell lines", "hybrid profile, DE found", "37%", "47%", [x[x.method == "hybrid_all"].tp.sum() / x[x.method == "hybrid_all"].n_de.sum() for x in oth])
in_range("DE, other cell lines", "direct all genes, DE found", "42%", "57%", [x[x.method == "direct_all"].tp.sum() / x[x.method == "direct_all"].n_de.sum() for x in oth])

# ---------------- simulation at equal sequencing budget (Fig. 5), folds 1 to 4
s = tab("simd.sim.tsv"); s = s[s.fold > 0]; dn = tab("simd.den.tsv"); dn = dn[dn.fold > 0]; S = "Simulation"
check(S, "test perturbations", "1,217", s.pert.nunique())
def sim(q, panel, assay): return s[(s.q == q) & (s.panel == panel) & (s.assay == assay)].set_index("pert")
def den(q): return dn[dn.q == q].set_index("pert")
for q, e_panel, e_den in [(0.02, "0.036", "0.001"), (0.05, "0.052", "0.026"), (0.1, "0.060", "0.059"), (0.25, "0.072", "0.075"), (1.0, "0.080", "0.092")]:
    h = sim(q, "greedyfull_500", "targeted"); g = den(q)
    check(S, f"targeted greedy MMSE, specific r at {q:g}", e_panel, h.rp_imp.median()); check(S, f"projected whole transcriptome, specific r at {q:g}", e_den, g.rp_den.median())
h, g = sim(0.02, "greedyfull_500", "targeted"), den(0.02); x = h.rp_imp - g.rp_den.reindex(h.index)
check(S, "greedy above projected whole transcriptome at 2%", "80%", (x > 0).mean()); check(S, "same comparison", "P < 0.001", wilcoxon(x).pvalue)
h5, g5 = sim(0.05, "greedyfull_500", "targeted"), den(0.05); check(S, "greedy above projected whole transcriptome at 5%", "72%", (h5.rp_imp - g5.rp_den.reindex(h5.index) > 0).mean())
in_range(S, "greedy above measured whole transcriptome, 2% to 25%", "77%", "81%", [((v := sim(q, "greedyfull_500", "targeted")).rp_imp > v.rp_dir).mean() for q in (0.02, 0.05, 0.1, 0.25)])
check(S, "DE found at 2%, targeted greedy", "43%", h.tp_imp.sum() / h.n_de.sum()); check(S, "DE found at 2%, measured whole transcriptome", "11%", g.tp_dir.sum() / g.n_de.sum()); check(S, "DE found at 2%, projected whole transcriptome", "9%", g.tp_den.sum() / g.n_de.sum())
c, l = sim(0.02, "cost0.5_500", "targeted_het"), sim(0.02, "l1000_679", "targeted"); y = c.rp_imp - l.rp_imp.reindex(c.index)
check(S, "UMI-cost panel at 2%", "0.054", c.rp_imp.median()); check(S, "targeted L1000 at 2%", "0.046", l.rp_imp.median()); check(S, "UMI-cost above L1000 at 2%", "62%", (y > 0).mean()); check(S, "same comparison", "P < 0.001", wilcoxon(y).pvalue)
check(S, "DE found at 2%, UMI-cost panel", "47%", c.tp_imp.sum() / c.n_de.sum()); check(S, "DE found at 2%, targeted L1000", "43%", l.tp_imp.sum() / l.n_de.sum())
c5, l5 = sim(0.05, "cost0.5_500", "targeted_het"), sim(0.05, "l1000_679", "targeted"); check(S, "UMI-cost above L1000 at 5%", "63%", (c5.rp_imp - l5.rp_imp.reindex(c5.index) > 0).mean())
check(S, "UMI-cost above projected whole transcriptome at 2%", "85%", (c.rp_imp - g.rp_den.reindex(c.index) > 0).mean()); check(S, "UMI-cost above projected whole transcriptome at 5%", "81%", (c5.rp_imp - g5.rp_den.reindex(c5.index) > 0).mean())
c10, g10 = sim(0.1, "cost0.5_500", "targeted_het"), den(0.1); check(S, "UMI-cost against projected whole transcriptome at 10%", "P = 0.09", wilcoxon(c10.rp_imp - g10.rp_den.reindex(c10.index)).pvalue)
check(S, "UMI share of the greedy MMSE panel", "52%", sim(1.0, "greedyfull_500", "targeted").share.median()); check(S, "UMI share of the L1000 landmarks", "6%", sim(1.0, "l1000_679", "targeted").share.median())
cv = tab("simb.cv2.tsv").dropna(); check(S, "counting share of gene-level variance, median gene", "94%", (cv.cv2_tech / cv.cv2_tot).median())

# ---------------- paired whole-transcriptome and L1000-enriched libraries (Fig. 6)
R = tab("paired.real.tsv"); S = "Paired libraries"
P = R.pivot(index="pert", columns="readout", values="rp"); rr = lambda a: P[a].median()
check(S, "testable knockdowns", "74", len(P)); check(S, "median input cells per knockdown", "447", np.ceil(R[R.readout == "R1_direct_full"].n_cells.median() / 2)); check(S, "median cells per knockdown", "893", R[R.readout == "R1_direct_full"].n_cells.median())
check(S, "enriched library, specific r", "0.117", rr("R3_l1000_enriched")); check(S, "extracted landmark genes, specific r", "0.124", rr("R2_l1000_extract"))
check(S, "enriched against extracted", "P = 0.40", wilcoxon(P["R3_l1000_enriched"] - P["R2_l1000_extract"]).pvalue)
check(S, "recovery, enriched", "0.69", rr("R3_l1000_enriched") / rr("D160_denoised")); check(S, "recovery, extracted", "0.72", rr("R2_l1000_extract") / rr("D160_denoised"))
check(S, "whole transcriptome at 52% of UMIs", "0.128", rr("R5_direct_qeq")); x = P["R5_direct_qeq"] - P["R3_l1000_enriched"]
check(S, "whole transcriptome above enriched at equal UMIs", "76%", (x > 0).mean()); check(S, "same comparison", "P < 0.001", wilcoxon(x).pvalue)
check(S, "whole transcriptome at 25% of UMIs", "0.097", rr("direct_q0.25")); check(S, "whole transcriptome at 10% of UMIs", "0.065", rr("direct_q0.1"))
de = lambda a, col: R[R.readout == a][col].sum() / R[R.readout == a].n_de.sum()
check(S, "enriched, DE found", "43%", de("R3_l1000_enriched", "tp")); check(S, "enriched, opposite sign", "8.9%", de("R3_l1000_enriched", "disc"))
check(S, "extracted, DE found", "45%", de("R2_l1000_extract", "tp")); check(S, "extracted, opposite sign", "6.1%", de("R2_l1000_extract", "disc"))
C = tab("paired.cells.tsv"); C["m"] = C["m"].astype(str)
M = C.groupby(["readout", "m", "pert"]).rp.mean().unstack("readout")
def share(m, a, b): d = (M.loc[m, a] - M.loc[m, b]).dropna(); return (d > 0).mean(), wilcoxon(d).pvalue
for m, e, pe in [("25", "73%", "P < 0.001"), ("50", "66%", "P < 0.001"), ("100", "62%", "P = 0.002")]:
    sh, pv = share(m, "l1000_enriched", "l1000_extract"); check(S, f"enriched above extracted, {m} cells", e, sh); check(S, "same comparison", pe, pv)
for m, e, pe in [("25", "65%", "P = 0.002"), ("50", "58%", "P = 0.010")]:
    sh, pv = share(m, "greedy_extract", "l1000_extract"); check(S, f"greedy MMSE above L1000, {m} cells", e, sh); check(S, "same comparison", pe, pv)
dc = lambda m, a, col: C[(C.m == m) & (C.readout == a)][col].sum() / C[(C.m == m) & (C.readout == a)].n_de.sum()
check(S, "direct measurement, DE found with 25 cells", "27%", dc("25", "direct", "tp")); check(S, "projected whole transcriptome, DE found with 25 cells", "46%", dc("25", "denoised", "tp"))
check(S, "direct measurement, DE found with 50 cells", "32%", dc("50", "direct", "tp")); check(S, "projected whole transcriptome, DE found with 50 cells", "48%", dc("50", "denoised", "tp"))
in_range(S, "projected whole transcriptome, specific r with 25 to 200 cells", "0.068", "0.133", [M.loc[m, "denoised"].median() for m in ("25", "50", "100", "200")])
in_range(S, "greedy MMSE genes, specific r with 25 to 200 cells", "0.047", "0.088", [M.loc[m, "greedy_extract"].median() for m in ("25", "50", "100", "200")])
log = open(os.path.join(SD, "real", "real.log")).read(); u = [float(v) for v in re.findall(r"whole transcriptome (\d+) \(panel genes (\d+), [\d.]+\); enriched (\d+) \(panel genes (\d+)", log)[0]]
check(S, "UMIs per cell, whole transcriptome", "10,793", u[0], "summary"); check(S, "landmark UMIs, whole transcriptome", "598", u[1], "summary")
check(S, "UMIs per cell, enriched", "5,664", u[2], "summary"); check(S, "landmark UMIs, enriched", "4,541", u[3], "summary"); check(S, "ratio of landmark UMIs", "7.6", u[3] / u[1], "summary")

# ---------------- independent atlases and transferred profiles (Figs. 7 and 8), summaries of the Rust tool
def xrec(f, panel):
    for line in open(os.path.join(SD, "xatlas", f)):
        t = line.split()
        if t and t[0] == panel: return float(t[2])
S = "Independent atlases"
check(S, "greedy MMSE recovery, HCT116", "0.88", xrec("hct116_cv.txt", "greedy_500"), "summary"); check(S, "greedy MMSE recovery, HEK293T", "0.73", xrec("hek293t_cv.txt", "greedy_500"), "summary")
def xde(f, method):
    for line in open(os.path.join(SD, "xatlas", f)):
        t = line.split()
        if t and t[0] == method: return float(t[1])
check(S, "hybrid profile, DE found, HCT116", "7.9%", xde("hct116_de.txt", "hybrid_all"), "summary"); check(S, "direct measurement, DE found, HCT116", "10.7%", xde("hct116_de.txt", "direct_all"), "summary")
check(S, "hybrid profile, DE found, HEK293T", "7.0%", xde("hek293t_de.txt", "hybrid_all"), "summary"); check(S, "direct measurement, DE found, HEK293T", "9.3%", xde("hek293t_de.txt", "direct_all"), "summary")
check(S, "random panel recovery, HCT116", "0.57", xrec("hct116_cv.txt", "random1_500"), "summary"); check(S, "random panel recovery, HEK293T", "0.40", xrec("hek293t_cv.txt", "random1_500"), "summary")
V = {}; scope = None
for line in open(os.path.join(SD, "vcm", "vcmsum.txt")):
    if line.startswith("== scope"): scope = line.split()[2]; continue
    t = line.split()
    if scope and len(t) > 5 and t[0] in ("RPE1", "HepG2", "Jurkat", "HCT116", "HEK293T"):
        try: V[(scope, t[0], t[1])] = [float(v) for v in t[3:6]]
        except ValueError: pass
S = "Transferred profiles"
check(S, "hybrid profile in RPE1, DE found", "49%", V[("all", "RPE1", "panel")][1], "summary"); check(S, "hybrid profile in RPE1, opposite sign", "3.7%", V[("all", "RPE1", "panel")][2], "summary")
check(S, "pooled transferred profile in RPE1, DE found", "29%", V[("all", "RPE1", "copy_pooled")][1], "summary"); check(S, "pooled transferred profile in RPE1, opposite sign", "18%", V[("all", "RPE1", "copy_pooled")][2], "summary")
in_range(S, "pooled transferred profile, specific r outside panel", "0.19", "0.23", [V[("out", c, "copy_pooled")][0] for c in ("RPE1", "HepG2", "Jurkat", "HCT116", "HEK293T")], "summary")
in_range(S, "panel imputation, specific r outside panel", "0.11", "0.14", [V[("out", c, "panel")][0] for c in ("RPE1", "HepG2", "Jurkat", "HCT116", "HEK293T")], "summary")

# ---------------- alternative imputation models and the guide-target prior (Table S1, Note S1), summaries of the Rust tool
A = pd.read_csv(os.path.join(SD, "alt", "benchmark.tsv"), sep="\t"); PA = pd.read_csv(os.path.join(SD, "alt", "paired_against_ridge.tsv"), sep="\t"); GG = pd.read_csv(os.path.join(SD, "alt", "gene_graphs.tsv"), sep="\t")
K = A[A.dataset.str.startswith("K562")].groupby("method").mean(numeric_only=True); S = "Effect size ranking"
RID, PRI = "ridge (all cells), baseline", "H2c + expression deciles"
check(S, "ridge, mean size rank, five folds", "0.073", K.loc[RID, "size_rank_mean"], "summary"); check(S, "two measured halves, mean size rank", "0.198", K.loc["direct (input half)", "size_rank_mean"], "summary")
check(S, "whole half on 160 programs, mean size rank", "0.095", K.loc["whole half on 160 programs", "size_rank_mean"], "summary")
alt = [K.loc[m, "size_rank_mean"] for m in ("H1 cross-half ridge", "H3 local experts", "H4 cross-half + quadratic")] + list(GG[GG.method.str.startswith("ridge +")].groupby("method").size_rank_mean.mean())
in_range(S, "alternative models, mean size rank", "0.070", "0.075", alt, "summary")
check(S, "guide-target prior, mean size rank", "0.139", K.loc[PRI, "size_rank_mean"], "summary")
check(S, "part of the difference closed by the prior", "0.066", K.loc[PRI, "size_rank_mean"] - K.loc[RID, "size_rank_mean"], "summary"); check(S, "difference between two measured halves and ridge", "0.125", K.loc["direct (input half)", "size_rank_mean"] - K.loc[RID, "size_rank_mean"], "summary")
check(S, "recovery, ridge (mean of fold medians)", "0.883", K.loc[RID, "recovery"], "summary"); check(S, "recovery, guide-target prior", "0.928", K.loc[PRI, "recovery"], "summary")
check(S, "DE found outside the panel, ridge", "0.206", K.loc[RID, "de_recall"], "summary"); check(S, "DE found outside the panel, guide-target prior", "0.229", K.loc[PRI, "de_recall"], "summary")
check(S, "opposite sign, ridge", "0.141", K.loc[RID, "opposite_sign"], "summary"); check(S, "opposite sign, guide-target prior", "0.126", K.loc[PRI, "opposite_sign"], "summary")
hp = PA[PA.method.str.startswith("H2c")].set_index("dataset")
check(S, "guide-target prior against ridge, largest P of the five folds", "P < 1e-17", hp[hp.index.str.startswith("K562")].size_rank_p.max(), "summary")
for ds, e0, e1 in [("RPE1", "0.095", "0.126"), ("HepG2", "0.055", "0.077"), ("Jurkat", "0.108", "0.150")]:
    x = A[A.dataset == ds].set_index("method"); check(S, f"{ds}, mean size rank, ridge", e0, x.loc[RID, "size_rank_mean"], "summary"); check(S, f"{ds}, mean size rank, guide-target prior", e1, x.loc[PRI, "size_rank_mean"], "summary")
check(S, "cross-half ridge, recovery", "0.820", K.loc["H1 cross-half ridge", "recovery"], "summary")
dg = open(os.path.join(SD, "alt", "diagnostic_fold0.txt")).read(); num = lambda pat: float(re.search(pat, dg).group(1))
check(S, "reproducible share of the variance of a half", "8.3%", num(r"all programs:.*= ([\d.]+)% of the total") / 100, "summary"); check(S, "reproducible share after removal of 10 programs", "3.5%", num(r"specific.*= ([\d.]+)% of the total") / 100, "summary")
check(S, "knockdowns with the target gene outside the panel", "511", num(r"outside the panel for (\d+) of"), "summary"); check(S, "inner holdout knockdowns", "787", num(r"inner holdout (\d+),"), "summary")
check(S, "median z of the target gene", "-3.8", num(r"median measured z (-[\d.]+)"), "summary"); check(S, "median z of genes within 2 kb", "-2.05", num(r"0 to 2 kb from the target: \d+ pairs, median z (-[\d.]+)"), "summary")
check(S, "leading eigenvalue of the noise covariance", "108.6", num(r"leading noise eigenvalues ([\d.]+)"), "summary")

# ---------------- inferences on complexes, pathways and strength (Fig. 5, Table S9), summaries of the Rust tool
BI = pd.read_csv(os.path.join(SD, "bio", "summary.tsv"), sep="\t"); BI["cell"] = BI.dataset.str.replace(r" fold \d", "", regex=True)
BIm = BI.groupby(["cell", "analysis", "frame", "readout", "metric"]).value.mean(); S = "Biological inference"
bi = lambda cell, an, ro, me: float(BIm[(cell, an, "all", ro, me)])
HY, WT, PA, AL = "hybrid + guide-target prior", "whole transcriptome, half", "panel genes alone", "whole transcriptome, all cells"
kn = BI[(BI.cell == "K562") & (BI.frame == "all") & (BI.readout == HY)]
check(S, "K562 knockdowns with a complex partner tested", "1,878", kn[kn.metric == "nearest is partner"].n.sum(), "summary"); check(S, "K562 responding knockdowns", "1,227", kn[kn.metric == "top set identical"].n.sum(), "summary")
for ro, e_auc, e_nn in [(HY, "0.715", "24.2%"), (WT, "0.718", "22.3%"), (PA, "0.720", "24.9%"), (AL, "0.730", "25.2%")]:
    check(S, f"complex AUC, {ro}", e_auc, bi("K562", "complex", ro, "pair AUROC"), "summary"); check(S, f"nearest knockdown is a partner, {ro}", e_nn, bi("K562", "complex", ro, "nearest is partner"), "summary")
for cell, e_h, e_w in [("RPE1", "0.740", "0.701"), ("HepG2", "0.755", "0.740"), ("Jurkat", "0.749", "0.738")]:
    check(S, f"complex AUC, hybrid, {cell}", e_h, bi(cell, "complex", HY, "pair AUROC"), "summary"); check(S, f"complex AUC, whole transcriptome, {cell}", e_w, bi(cell, "complex", WT, "pair AUROC"), "summary")
for me, e_h, e_w in [("median r per gene set", "0.453", "0.357"), ("median r, responders", "0.737", "0.698"), ("top set identical", "53.1%", "52.0%"), ("top set same sign", "97.4%", "97.7%")]:
    check(S, f"pathway {me}, hybrid", e_h, bi("K562", "pathway", HY, me), "summary"); check(S, f"pathway {me}, whole transcriptome", e_w, bi("K562", "pathway", WT, me), "summary")
for cell, e_h, e_w in [("RPE1", "0.737", "0.755"), ("HepG2", "0.641", "0.678"), ("Jurkat", "0.625", "0.596")]:
    check(S, f"pathway r, hybrid, {cell}", e_h, bi(cell, "pathway", HY, "median r per gene set"), "summary"); check(S, f"pathway r, whole transcriptome, {cell}", e_w, bi(cell, "pathway", WT, "median r per gene set"), "summary")
check(S, "DE count rank correlation, hybrid", "0.604", bi("K562", "strength", HY, "Spearman of DE counts"), "summary"); check(S, "DE count rank correlation, whole transcriptome", "0.609", bi("K562", "strength", WT, "Spearman of DE counts"), "summary")
in_range(S, "DE count rank correlation, hybrid, other cell lines", "0.798", "0.904", [bi(c, "strength", HY, "Spearman of DE counts") for c in ("RPE1", "HepG2", "Jurkat")], "summary")
in_range(S, "DE count rank correlation, whole transcriptome, other cell lines", "0.804", "0.910", [bi(c, "strength", WT, "Spearman of DE counts") for c in ("RPE1", "HepG2", "Jurkat")], "summary")
dmax = max(abs(float(BIm[(c, an, "all", HY, me)]) - float(BIm[(c, an, "all", "hybrid (panel + imputed)", me)])) for (c, an, fr, ro, me) in BIm.index if ro == HY and fr == "all")
check(S, "largest difference between hybrid with and without the prior", "0.011", dmax, "summary")

# ---------------- perturbations with few cells and internal references (Table S10)
S = "Few cells"; RS = pd.read_csv(os.path.join(SD, "rescue", "summary.tsv"), sep="\t", dtype={"input_cells": str})
DIR, HYE, PRE, PRI, SHI = "direct", "hybrid, external", "projected, external", "projected, internal", "shrunk, internal"
def rs(screen, level, stratum, ro, col): return float(RS[(RS.screen == screen) & (RS.input_cells == level) & (RS.stratum == stratum) & (RS.readout == ro)][col].iloc[0])
def low(c):
    """per-perturbation table of one screen, knockdowns with 50 or fewer input cells"""
    t = tab(f"rescue.{c}.tsv"); return t[(t.input_cells.astype(str) == "all") & (t.n_in <= 50)]
def pooled(t, ro, col): x = t[t.readout == ro]; return x[col].sum() / x.n_de.sum()
E15 = tab("rescue.k562ess.tsv"); E15["input_cells"] = E15.input_cells.astype(str); e15 = E15[E15.input_cells == "15"]
check(S, "knockdowns of the essential-gene screen tested", "1,055", e15[e15.readout == DIR].shape[0])
for ro, e_tp, e_op in [(DIR, "36%", "14%"), (HYE, "53%", "6%"), (PRE, "56%", "3%"), (PRI, "58%", "3%")]:
    check(S, f"15 cells, DE genes found, {ro}", e_tp, pooled(e15, ro, "tp")); check(S, f"15 cells, opposite-sign calls, {ro}", e_op, pooled(e15, ro, "disc"))
for ro, e_pw, e_cx in [(DIR, "0.570", "0.848"), (PRE, "0.726", "0.869"), (PRI, "0.728", "0.891")]:
    check(S, f"15 cells, pathway r, {ro}", e_pw, rs("K562 essential", "15", "all", ro, "pathway_r"), "summary"); check(S, f"15 cells, complex AUC, {ro}", e_cx, rs("K562 essential", "15", "all", ro, "complex_auroc"), "summary")
k50 = low("k562")
check(S, "K562 knockdowns with 50 or fewer input cells", "1,678", k50[k50.readout == DIR].shape[0])
check(S, "K562, DE genes found, direct", "24%", pooled(k50, DIR, "tp")); check(S, "K562, DE genes found, projected on the atlas", "30%", pooled(k50, PRE, "tp"))
check(S, "K562, opposite-sign calls, direct", "22%", pooled(k50, DIR, "disc")); check(S, "K562, opposite-sign calls, projected on the atlas", "10%", pooled(k50, PRE, "disc"))
o50 = [low(c) for c in ("rpe1", "hepg2", "jurkat")]
rows.append(("PASS" if (min(len(t[t.readout == DIR]) for t in o50), max(len(t[t.readout == DIR]) for t in o50)) == (1229, 1796) else "FAIL", S, "knockdowns with 50 or fewer input cells, other screens", "1,229 to 1,796", ", ".join(str(len(t[t.readout == DIR])) for t in o50), "tables"))
in_range(S, "DE genes found, projected on the atlas, other screens", "31%", "45%", [pooled(t, PRE, "tp") for t in o50]); in_range(S, "DE genes found, direct, other screens", "33%", "53%", [pooled(t, DIR, "tp") for t in o50])
in_range(S, "DE genes found, projected on the internal reference, other screens", "36%", "56%", [pooled(t, PRI, "tp") for t in o50])
in_range(S, "opposite-sign calls, direct, other screens", "5.2%", "14.5%", [pooled(t, DIR, "disc") for t in o50]); in_range(S, "opposite-sign calls, internal reference, other screens", "1.6%", "4.9%", [pooled(t, PRI, "disc") for t in o50])
def paired(t, ro):
    w = t.pivot(index=["fold", "pert"], columns="readout", values="rp_own"); d = (w[ro] - w[DIR]).dropna(); return (d > 0).mean(), wilcoxon(d).pvalue
in_range(S, "knockdowns with a higher specific correlation, internal projection, other screens", "75%", "82%", [paired(t, PRI)[0] for t in o50])
check(S, "largest P value, internal projection against direct", "P < 0.001", max(paired(t, PRI)[1] for t in o50))
in_range(S, "knockdowns with a higher specific correlation, shrunk profile, 15 cells to all cells", "86%", "89%", [paired(E15[E15.input_cells == l], SHI)[0] for l in ("15", "30", "50", "all")])
eall = E15[E15.input_cells == "all"]
for ro, e in [(SHI, "0.200"), (DIR, "0.116"), (PRI, "0.189")]: check(S, f"all cells, median specific correlation, {ro}", e, eall[eall.readout == ro].rp_own.median())
RF = pd.read_csv(os.path.join(SD, "rescue", "reference.tsv"), sep="\t"); RI = RF[RF.reference == "internal"]
rng = lambda sel: (int(sel.n_reference.min()), int(sel.n_reference.max()))
for nm, sel, e in [("K562", RI[RI.screen == "K562"], (3847, 3974)), ("other screens", RI[RI.screen.isin(["RPE1", "HepG2", "Jurkat"])], (854, 977)), ("essential-gene screen", RI[RI.screen == "K562 essential"], (712, 826))]:
    rows.append(("PASS" if rng(sel) == e else "FAIL", S, f"internal reference knockdowns, {nm}", f"{e[0]:,} to {e[1]:,}", f"{rng(sel)[0]}, {rng(sel)[1]}", "summary"))
rows.append(("PASS" if set(RI.programs) == {100, 160} else "FAIL", S, "programs of the internal reference", "100 or 160", ", ".join(map(str, sorted(set(RI.programs)))), "summary"))
# Table S10 of the supplement agrees with the summary, row by row
TS = os.path.join(SD, "rescue", "table_s10.md")
if os.path.exists(TS):
    n_ok = n_all = 0
    for l in open(TS):
        f = [x.strip() for x in l.strip().strip("|").split("|")]
        if len(f) != 9 or not f[1].isdigit(): continue
        lev, strat, scr = (f[0].lower(), "all", "K562 essential") if f[0] in ("15", "30", "50", "All") else ("all", "<=50", f[0])
        ro = {"Direct measurement": "direct", "Hybrid profile": "hybrid", "Hybrid profile with guide-target prior": "hybrid + prior", "Projected whole transcriptome": "projected", "Shrunk whole transcriptome": "shrunk"}[f[3]]
        ro = ro if ro == "direct" else ro + (", external" if f[2] == "Atlas" else ", internal")
        vals = [rs(scr, lev, strat, ro, c) for c in ("median_rp_own", "de_recall", "opposite", "pathway_r", "complex_auroc")]
        n_all += 1; n_ok += all(abs(float(a) - b) <= 5.1e-4 for a, b in zip(f[4:], vals))
    rows.append(("PASS" if n_ok == n_all == 72 else "FAIL", S, "rows of Table S10 equal to the summary", "72", str(n_ok), "summary"))

# ---------------- values that need the full atlas
S = "Reference atlas"
check(S, "programs above the noise level, K562", "172", N["mp"]["n_programs"], "json"); check(S, "genes analysed", "7,094", N["noise"]["n_genes_universe"], "json"); check(S, "knockdowns", "9,866", N["noise"]["n_perts_gwps"], "json")
check(S, "excess of differences beyond 4 SD over the normal expectation (eight times)", "8", N["noise"]["tail4"] / N["noise"]["normal_tail4"], "json")
check(S, "size ranking, all genes projected on 160 programs", "0.21", N["rank_oracle"]["K160"], "json"); check(S, "size ranking, imputation from 500 genes", "0.19", N["rank_oracle"]["ridge500"], "json")

w = [max(len(str(r[i])) for r in rows) for i in range(6)]
print("  ".join(h.ljust(w[i]) for i, h in enumerate(["", "section", "claim", "manuscript", "recomputed", "source"])))
for r in rows: print("  ".join(str(v).ljust(w[i]) for i, v in enumerate((f"[{r[0]}]",) + r[1:])))
nf = sum(r[0] == "FAIL" for r in rows)
print(f"\n{len(rows) - nf} of {len(rows)} checks passed" + (f", {nf} FAILED" if nf else ""))
sys.exit(1 if nf else 0)
