"""Figures 2 to 8 and Additional file Fig. S1 from source_data/ (tables, summaries and numbers.json).
Usage: python3 scripts/make_figures.py [OUTPUT_DIR]   (default: output/)
Style: plotly_white, default Plotly colours, Arial, black text, no gridlines, lowercase panel letters."""
import os, re, json, numpy as np, pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "source_data"); O = os.path.join(D, "tables") + "/"
F = (sys.argv[1] if len(sys.argv) > 1 else os.path.join(D, "..", "output")).rstrip("/") + "/"; os.makedirs(F, exist_ok=True)
N = json.load(open(os.path.join(D, "numbers.json")))
C = ["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3", "#FF6692"]; BLACK = "#111111"; GREY = "#9a9a9a"
SEL = [("greedy", "Greedy MMSE", C[0]), ("hvgresp", "Response variance", C[2]), ("designed", "Pivoted QR", C[3]), ("genebasis", "geneBasis", C[4]),
       ("hvg", "Highly variable (control cells)", C[5]), ("module", "Program modules", C[6]), ("random1", "Random", "#FECB52")]
CTS = ["K562", "RPE1", "HepG2", "Jurkat"]

def style(fig, w, h, l=60, b=60, t=40, r=20):
    fig.update_layout(template="plotly_white", width=w, height=h, font=dict(family="Arial", size=13, color=BLACK),
                      margin=dict(l=l, r=r, t=t, b=b), legend=dict(font=dict(color=BLACK, size=12)))
    fig.update_xaxes(showgrid=False, zeroline=False, showline=True, linecolor=BLACK, ticks="outside", tickfont=dict(color=BLACK), title_font=dict(color=BLACK))
    fig.update_yaxes(showgrid=False, zeroline=False, showline=True, linecolor=BLACK, ticks="outside", tickfont=dict(color=BLACK), title_font=dict(color=BLACK))
    for a in fig.layout.annotations:
        if a.text and len(a.text) > 2 and a.text[1:3] == "  ": a.font = dict(family="Arial", size=14, color=BLACK); a.x = a.x; a.xanchor = "center"
    return fig
def save(fig, name): fig.write_image(F + name + ".png", scale=4); fig.write_image(F + name + ".pdf")
def resp_pairs(pre, panels):
    R = pd.read_csv(O + pre + ".perts.tsv.gz", sep="\t"); R["key"] = R.fold.astype(str) + ":" + R.pert
    b = R[R.panel == R.panel.iloc[0]].set_index("key"); resp = b.index[b.r_ceil >= 0.2]
    W = R.pivot(index="key", columns="panel", values="r_imp_perp").loc[resp]; return W[panels], b.loc[resp]

# ---------------- Figure 2: panel selectors
fig = make_subplots(1, 3, subplot_titles=["a  K562, by panel size", "b  500 genes, by cell line", "c  K562, 680 genes"], horizontal_spacing=0.1, column_widths=[0.34, 0.36, 0.30])
k = N["K562"]
for s, lab, col in SEL:
    xs = [x for x in (200, 500, 680) if f"ratio_{s}_{x}" in k]
    fig.add_trace(go.Scatter(x=xs, y=[k[f"ratio_{s}_{x}"] for x in xs], mode="lines+markers", name=lab, line=dict(color=col, width=2), marker=dict(size=8, color=col), legendgroup=s), 1, 1)
fig.add_trace(go.Scatter(x=[680], y=[k["ratio_l1000_680"]], mode="markers", name="L1000 landmarks", marker=dict(size=11, symbol="diamond", color=C[1]), legendgroup="l1000"), 1, 1)
for i, (s, lab, col) in enumerate(SEL):
    fig.add_trace(go.Scatter(x=[j + (i - 3) * 0.1 for j in range(4)], y=[N[ct].get(f"ratio_{s}_500", np.nan) for ct in CTS], mode="markers",
                             marker=dict(size=9, color=col, symbol="diamond"), showlegend=False, legendgroup=s), 1, 2)
fig.update_xaxes(tickvals=[200, 500, 680], title_text="Genes in panel", range=[150, 730], row=1, col=1)
fig.update_xaxes(tickvals=list(range(4)), ticktext=CTS, range=[-0.6, 3.6], row=1, col=2)
fig.update_yaxes(title_text="Specific response recovered<br>(fraction of denoised full profile)", range=[0.35, 0.95], row=1, col=1)
fig.update_yaxes(range=[0.35, 0.95], row=1, col=2)
W, b = resp_pairs("cvg", ["greedy_680", "l1000_680"])
fig.add_trace(go.Scatter(x=W["l1000_680"], y=W["greedy_680"], mode="markers", marker=dict(size=4, color=C[0], opacity=0.5), showlegend=False), 1, 3)
fig.add_trace(go.Scatter(x=[-0.1, 0.7], y=[-0.1, 0.7], mode="lines", line=dict(color=GREY, dash="dash", width=1), showlegend=False), 1, 3)
fig.update_xaxes(title_text="L1000 landmarks, specific r", range=[-0.1, 0.7], row=1, col=3); fig.update_yaxes(title_text="Greedy MMSE, specific r", range=[-0.1, 0.7], row=1, col=3)
style(fig, 1300, 470, b=70); fig.update_layout(legend=dict(orientation="h", y=-0.22, x=0.0)); save(fig, "Figure2")

# ---------------- Figure 3: greedy 500 imputation, INFO and OOR
M = N["main"]; PRE = {"K562": "cvm", "RPE1": "trm_rpe1", "HepG2": "trm_hepg2", "Jurkat": "trm_jurkat"}
fig = make_subplots(2, 3, specs=[[{}, {}, {}], [{}, {"colspan": 2}, None]], subplot_titles=["a  Recovery per perturbation", "b  INFO calibration", "c  Out-of-reference score", "d  K562", "e  MBTPS2 knockdown, K562 test perturbation"], horizontal_spacing=0.08, vertical_spacing=0.16)
POS = {1: (1, 1), 2: (1, 2), 3: (1, 3), 4: (2, 1)}
for j, ct in enumerate(CTS):
    R = pd.read_csv(O + PRE[ct] + ".perts.tsv.gz", sep="\t"); g = R[(R.panel == "greedy_500") & (R.r_ceil >= 0.2)]
    fig.add_trace(go.Box(y=g.r_imp_perp / g.r_den160_perp, name=ct, marker=dict(color=C[j], size=3), line=dict(color=C[j]), boxpoints=False, showlegend=False), *POS[1])
    fig.add_trace(go.Box(y=g.oor, name=ct, marker=dict(color=C[j], size=3), line=dict(color=C[j]), boxpoints=False, showlegend=False), *POS[3])
    bi = np.array(M[ct]["info_bins"])
    fig.add_trace(go.Scatter(x=bi[:, 0], y=bi[:, 1], mode="lines+markers", name=ct, line=dict(color=C[j]), marker=dict(size=7, color=C[j])), *POS[2])
    if ct == "K562":
        gg = g.dropna(subset=["oor"]); fail = gg.r_imp_perp / gg.r_ceil_perp < 0.5
        fig.add_trace(go.Scatter(x=gg.oor[~fail], y=(gg.r_imp_perp / gg.r_den160_perp)[~fail], mode="markers", marker=dict(size=4, color=C[0], opacity=0.5), showlegend=False), *POS[4])
        fig.add_trace(go.Scatter(x=gg.oor[fail], y=(gg.r_imp_perp / gg.r_den160_perp)[fail], mode="markers", marker=dict(size=7, color=C[1], symbol="x"), showlegend=False), *POS[4])
fig.add_trace(go.Scatter(x=[0, 0.6], y=[0, 0.6], mode="lines", line=dict(color=GREY, dash="dash", width=1), showlegend=False), *POS[2])
fig.update_yaxes(title_text="Specific response recovered", range=[-0.2, 1.6], row=POS[1][0], col=POS[1][1])
fig.update_xaxes(title_text="INFO (predicted)", range=[0, 0.6], row=POS[2][0], col=POS[2][1]); fig.update_yaxes(title_text="Observed per-gene r", range=[0, 0.6], row=POS[2][0], col=POS[2][1])
fig.update_yaxes(title_text="Out-of-reference score", range=[0, 0.7], row=POS[3][0], col=POS[3][1])
fig.update_xaxes(title_text="Out-of-reference score", row=POS[4][0], col=POS[4][1]); fig.update_yaxes(title_text="Specific response recovered", range=[-0.2, 1.6], row=POS[4][0], col=POS[4][1])
E = pd.read_csv(os.path.join(D, "example", "mbtps2_fold1.tsv"), sep="\t")
SREBP = ["HMGCS1", "HMGCR", "LDLR", "INSIG1", "SQLE", "FDFT1", "IDI1", "MVD", "DHCR7", "SC5D", "CYP51A1", "MSMO1", "NSDHL", "FDPS", "ACAT2", "TMEM97", "EBP", "SCD", "FASN"]
o = E[~E.gene.isin(SREBP)]
fig.add_trace(go.Scatter(x=o.imputed, y=o.half1, mode="markers", marker=dict(size=3, color=GREY, opacity=0.35), name="Other genes"), 2, 2)
for flag, lab, col, sym in ((1, "SREBP targets in the panel (measured)", C[0], "circle"), (0, "SREBP targets outside the panel (imputed)", C[1], "diamond")):
    t = E[E.gene.isin(SREBP) & (E.in_panel == flag)]
    lab_genes = {"ACAT2", "SQLE", "INSIG1", "LDLR", "HMGCR", "FDFT1", "NSDHL"}
    fig.add_trace(go.Scatter(x=t.imputed, y=t.half1, mode="markers+text", text=[g if g in lab_genes else "" for g in t.gene], textposition="middle right", textfont=dict(size=9, color=col), marker=dict(size=8, color=col, symbol=sym), name=lab), 2, 2)
fig.add_vline(x=0, line=dict(color=GREY, width=1), row=2, col=2); fig.add_hline(y=0, line=dict(color=GREY, width=1), row=2, col=2)
fig.update_xaxes(title_text="Input half: measured (panel) or imputed (other genes), z", range=[-12, 12], row=2, col=2); fig.update_yaxes(range=[-12, 12], row=2, col=2); fig.update_yaxes(title_text="Held-out half, measured z", row=2, col=2)
style(fig, 1300, 860); fig.update_layout(legend=dict(orientation="v", yanchor="top", y=0.42, x=0.36, font=dict(size=11), bgcolor="rgba(255,255,255,0.8)")); save(fig, "Figure3")

# ---------------- Figure 4: DE presence, direction and size, out of panel and hybrid readout
H = N["hybrid"]
fig = make_subplots(2, 3, specs=[[{}, {}, {}], [{}, {}, None]], vertical_spacing=0.2, horizontal_spacing=0.09,
                    subplot_titles=["a  DE genes recovered", "b  Opposite sign", "c  Size ranking among DE genes", "d  Effect size calibration", "e  Size ranking, K562 fold 0"])
x = list(range(4))
bars = [("baseline", "Imputed, genes outside the panel", C[0], -0.3), ("direct_out", "Measured, genes outside the panel", C[4], -0.1),
        ("hybrid_all", "Hybrid, panel measured and other genes imputed", C[2], 0.1), ("direct_all", "Measured, all genes", C[1], 0.3)]
for key, lab, col, off in bars:
    for c_, met in ((1, "recall"), (2, "opp"), (3, "rank")):
        fig.add_trace(go.Bar(x=[i + off for i in x], y=[H[ct][key][met] for ct in CTS], width=0.19, marker_color=col, name=lab, showlegend=c_ == 1, legend="legend", marker_line=dict(color="white", width=1)), 1, c_)
for c_ in (1, 2, 3): fig.update_xaxes(tickvals=x, ticktext=CTS, row=1, col=c_)
fig.update_yaxes(title_text="Fraction of DE genes,<br>correct sign", range=[0, 0.62], row=1, col=1)
fig.update_yaxes(title_text="Fraction of calls with opposite sign", range=[0, 0.25], row=1, col=2)
fig.update_yaxes(title_text="Rank correlation of |change|", range=[0, 0.5], row=1, col=3)
for key, lab, col, off in [("slope_raw", "Before calibration", C[3], -0.1), ("slope_cal", "After calibration", C[5], 0.1)]:
    fig.add_trace(go.Scatter(x=[i + off for i in x], y=[M[ct]["de"][key] for ct in CTS], mode="markers", marker=dict(size=10, symbol="diamond", color=col), name=lab, legend="legend2"), 2, 1)
fig.add_hline(y=1, line=dict(color=GREY, dash="dash", width=1), row=2, col=1)
fig.update_xaxes(tickvals=x, ticktext=CTS, row=2, col=1); fig.update_yaxes(title_text="Slope, measured on imputed", range=[0.6, 1.4], row=2, col=1)
o = N["rank_oracle"]
fig.add_trace(go.Scatter(x=[20, 160, 500, 2000], y=[o["K20"], o["K160"], o["K500"], o["K2000"]], mode="lines+markers", line=dict(color=C[2]), marker=dict(size=8, color=C[2]), name="All genes projected on K programs", legend="legend3"), 2, 2)
fig.add_trace(go.Scatter(x=[20, 2000], y=[o["direct"]] * 2, mode="lines", line=dict(color=C[1], dash="dash"), name="Measured, all genes", legend="legend3"), 2, 2)
fig.add_trace(go.Scatter(x=[20, 2000], y=[o["ridge500"]] * 2, mode="lines", line=dict(color=C[0], dash="dot"), name="Imputed from 500 genes (pivoted QR)", legend="legend3"), 2, 2)
fig.update_xaxes(type="log", title_text="Programs (K)", tickvals=[20, 160, 500, 2000], row=2, col=2)
fig.update_yaxes(title_text="Rank correlation of |change|", range=[0, 0.35], row=2, col=2)
style(fig, 1300, 900)
y1 = fig.layout.yaxis.domain; xd, xe = fig.layout.xaxis4.domain, fig.layout.xaxis5.domain; y2 = fig.layout.yaxis4.domain
L4 = dict(yanchor="top", xanchor="left", font=dict(color=BLACK, size=12), bgcolor="rgba(255,255,255,0.85)")
fig.update_layout(barmode="overlay",
    legend=dict(L4, orientation="h", x=0.0, y=y1[0] - 0.075),
    legend2=dict(L4, orientation="v", x=xd[0] + 0.12, y=y2[0] + 0.16),
    legend3=dict(L4, orientation="v", x=xe[0] + 0.02, y=y2[0] + 0.16))
save(fig, "Figure4")

# ---------------- Figure 5: inferences on complexes and pathways (summary written by the Rust tool, five K562 folds averaged)
B = pd.read_csv(os.path.join(D, "bio", "summary.tsv"), sep="\t"); B["cell"] = B.dataset.str.replace(r" fold \d", "", regex=True)
BM = B.groupby(["cell", "analysis", "frame", "readout", "metric"]).value.mean()
bv = lambda cell, an, ro, me: float(BM[(cell, an, "all", ro, me)])
RO = [("whole transcriptome, all cells", "Whole transcriptome, all cells", "#4d4d4d"), ("whole transcriptome, half", "Whole transcriptome, measured", GREY),
      ("panel genes alone", "Panel genes alone", C[4]), ("hybrid + guide-target prior", "Hybrid profile with guide-target prior", C[2])]
fig = make_subplots(1, 3, subplot_titles=["a  Complex partners", "b  Pathway scores", "c  Strongest pathway"], horizontal_spacing=0.09)
for key, lab, col in RO:
    fig.add_trace(go.Bar(x=CTS, y=[bv(c, "complex", key, "nearest is partner") for c in CTS], name=lab, marker_color=col, legendgroup=lab), 1, 1)
    if key in ("whole transcriptome, half", "hybrid + guide-target prior"):
        fig.add_trace(go.Bar(x=CTS, y=[bv(c, "pathway", key, "median r per gene set") for c in CTS], name=lab, marker_color=col, legendgroup=lab, showlegend=False), 1, 2)
        fig.add_trace(go.Bar(x=CTS, y=[bv(c, "pathway", key, "top set identical") for c in CTS], name=lab, marker_color=col, legendgroup=lab, showlegend=False), 1, 3)
fig.update_yaxes(title_text="Nearest knockdown in the same complex", range=[0, 0.5], row=1, col=1)
fig.update_yaxes(title_text="Median r with the held-out half", range=[0, 0.85], row=1, col=2)
fig.update_yaxes(title_text="Strongest gene set identical", range=[0, 0.9], row=1, col=3)
style(fig, 1100, 430); fig.update_layout(barmode="group", legend=dict(orientation="h", y=-0.14, x=0.0))
fig.for_each_trace(lambda t: t.update(marker_line=dict(color="white", width=1)))
save(fig, "Figure5")

# ---------------- Figure 6: equal sequencing budget (folds 1-4)
S = N["sim"]; Q = ["0.02", "0.05", "0.1", "0.25", "1.0"]; qx = [float(q) * 100 for q in Q]
fig = make_subplots(1, 2, subplot_titles=["a  Specific response", "b  DE genes recovered"], horizontal_spacing=0.12)
lines = [("direct", "Whole transcriptome, measured", GREY, "dash"), ("greedyfull_500|targeted", "Targeted greedy MMSE panel", C[0], "solid"),
         ("cost0.5_500|targeted_het", "Targeted UMI-cost panel, noise-weighted imputation", C[3], "solid"),
         ("l1000_679|targeted", "Targeted L1000 landmarks", C[1], "solid"), ("random1_500|targeted", "Targeted random panel", C[2], "solid"),
         ("greedyfull_500|extract", "Greedy MMSE genes within whole transcriptome", C[0], "dot")]
for key, lab, col, dash in lines:
    if key == "direct": y1 = [S[q]["direct_rp"] for q in Q]; y2 = [S[q]["direct_sens"] for q in Q]
    else: y1 = [S[q][key]["rp"] for q in Q]; y2 = [S[q][key]["sens"] for q in Q]
    for c, y in ((1, y1), (2, y2)):
        fig.add_trace(go.Scatter(x=qx, y=y, mode="lines+markers", name=lab, line=dict(color=col, dash=dash, width=2), marker=dict(size=7, color=col), showlegend=c == 1), 1, c)
DEN = pd.read_csv(os.path.join(D, "sim", "denoised.tsv"), sep="\t"); DEN = DEN[DEN.readout == "denoised"].set_index("q")
for c, col_ in ((1, "specific_r"), (2, "de_found")):
    fig.add_trace(go.Scatter(x=qx, y=[DEN.loc[float(q), col_] for q in Q], mode="lines+markers", name="Whole transcriptome, projected on reference programs", line=dict(color=BLACK, dash="dash", width=2), marker=dict(size=7, color=BLACK), showlegend=c == 1), 1, c)
for c in (1, 2): fig.update_xaxes(type="log", title_text="Sequencing budget (% of reads)", tickvals=qx, ticktext=["2", "5", "10", "25", "100"], row=1, col=c)
fig.update_yaxes(title_text="Median specific r", range=[-0.005, 0.1], row=1, col=1); fig.update_yaxes(title_text="Fraction of DE genes, correct sign", range=[0, 0.6], row=1, col=2)
style(fig, 1000, 500); fig.update_layout(legend=dict(orientation="h", y=-0.25, x=0.0)); save(fig, "Figure6")

# ---------------- Figure 7: paired whole-transcriptome and L1000-enriched libraries (summaries written by the Rust tool)
RS = os.path.join(D, "real") + "/"
rs = {}; 
for l in open(RS + "real.summary.tsv"):
    t = l.rstrip("\n").split("\t")
    if len(t) == 6 and t[0] != "readout": rs[t[0]] = dict(rp=float(t[3]), sens=float(t[4]), opp=float(t[5]))
cs = {}
for l in open(RS + "cells.summary.tsv"):
    t = l.rstrip("\n").split("\t")
    if len(t) == 6 and t[0] != "m": cs[(t[0], t[1])] = dict(rp=float(t[3]), sens=float(t[4]), opp=float(t[5]))
umi = [float(x) for x in re.findall(r"whole transcriptome (\d+) \(panel genes (\d+), [\d.]+\); enriched (\d+) \(panel genes (\d+)", open(RS + "real.log").read())[0]]
fig = make_subplots(2, 2, subplot_titles=["a  UMIs per cell", "b  Equal number of UMIs", "c  Specific response by input cells", "d  DE genes found by input cells"], horizontal_spacing=0.14, vertical_spacing=0.2)
fig.add_trace(go.Bar(x=["Whole transcriptome", "L1000-enriched"], y=[umi[1], umi[3]], name="L1000 landmark genes", marker_color=C[1], legend="legend"), 1, 1)
fig.add_trace(go.Bar(x=["Whole transcriptome", "L1000-enriched"], y=[umi[0] - umi[1], umi[2] - umi[3]], name="Other genes", marker_color=GREY, legend="legend"), 1, 1)
fig.update_yaxes(title_text="Mean UMIs per cell", row=1, col=1)
qeq = umi[2] / umi[0] * 100
xq = [5, 10, 25, 50, 100]; yq = [rs["direct_q0.05"]["rp"], rs["direct_q0.1"]["rp"], rs["direct_q0.25"]["rp"], rs["direct_q0.5"]["rp"], rs["R1_direct_full"]["rp"]]
fig.add_trace(go.Scatter(x=xq, y=yq, mode="lines+markers", name="Whole transcriptome, measured", line=dict(color=GREY, dash="dash", width=2), marker=dict(size=7, color=GREY), legend="legend2"), 1, 2)
fig.add_trace(go.Scatter(x=[qeq], y=[rs["R3_l1000_enriched"]["rp"]], mode="markers", name="L1000-enriched library, imputed", marker=dict(size=13, symbol="diamond", color=C[1]), legend="legend2"), 1, 2)
fig.add_trace(go.Scatter(x=[100], y=[rs["R2_l1000_extract"]["rp"]], mode="markers", name="L1000 genes of whole transcriptome, imputed", marker=dict(size=12, symbol="diamond-open", color=C[1], line=dict(width=2)), legend="legend2"), 1, 2)
fig.update_xaxes(type="log", title_text="UMIs per cell (% of whole transcriptome)", tickvals=[5, 10, 25, 50, 100], row=1, col=2); fig.update_yaxes(title_text="Median specific r", range=[0, 0.2], row=1, col=2)
MS = ["25", "50", "100", "200", "all"]; mx = [25, 50, 100, 200, 447]
RL = [("direct", "Whole transcriptome, measured", GREY, "dash"), ("denoised", "Whole transcriptome, projected on reference programs", BLACK, "dash"), ("l1000_enriched", "L1000-enriched library, imputed", C[1], "solid"), ("l1000_extract", "L1000 genes of whole transcriptome, imputed", C[1], "dot"),
      ("greedy_extract", "Greedy MMSE genes of whole transcriptome, imputed", C[0], "dot"), ("random_extract", "Random genes of whole transcriptome, imputed", "#FECB52", "dot")]
for key, lab, col, dash in RL:
    for c, f in ((1, "rp"), (2, "sens")):
        fig.add_trace(go.Scatter(x=mx, y=[cs[(m, key)][f] for m in MS], mode="lines+markers", name=lab, line=dict(color=col, dash=dash, width=2), marker=dict(size=7, color=col), showlegend=c == 1, legend="legend3"), 2, c)
for c in (1, 2): fig.update_xaxes(type="log", title_text="Input cells per perturbation", tickvals=mx, ticktext=["25", "50", "100", "200", "all"], row=2, col=c)
fig.update_yaxes(title_text="Median specific r", range=[0, 0.2], row=2, col=1); fig.update_yaxes(title_text="Fraction of DE genes, correct sign", range=[0, 0.65], row=2, col=2)
style(fig, 1000, 900)
fig.update_layout(barmode="stack", legend=dict(orientation="h", x=0.0, y=0.56, xanchor="left", yanchor="top"), legend2=dict(orientation="v", x=0.58, y=0.995, xanchor="left", yanchor="top", bgcolor="rgba(255,255,255,0)"),
                  legend3=dict(orientation="h", x=0.0, y=-0.1, xanchor="left", yanchor="top"))
save(fig, "Figure7")

# ---------------- Additional file Fig. S1: budget-aware MMSE panels (out/simb, all folds)
SB = N["sim_budget"]
fig = make_subplots(1, 1)
for key, lab, col in [("greedyfull_500", "Greedy MMSE (budget not modelled)", C[0]), ("designed_500", "Pivoted QR", C[3]), ("l1000_679", "L1000 landmarks", C[1]), ("random1_500", "Random", C[2])]:
    fig.add_trace(go.Scatter(x=qx[:4], y=[SB[q][key] for q in Q[:4]], mode="lines+markers", name=lab, line=dict(color=col), marker=dict(color=col, size=7)))
fig.add_trace(go.Scatter(x=qx[:4], y=[SB[q][f"budget{q}_500"] for q in Q[:4]], mode="markers", name="Panel designed for that budget by depth-aware MMSE", marker=dict(size=11, symbol="diamond", color=BLACK)))
fig.update_xaxes(type="log", title_text="Sequencing budget (% of reads)", tickvals=qx[:4], ticktext=["2", "5", "10", "25"]); fig.update_yaxes(title_text="Median specific r", range=[0, 0.1])
style(fig, 760, 450); save(fig, "FigureS1")
# ---------------- Figure 8: X-Atlas/Orion atlases and pooled references (summaries written by the Rust tool)
X = os.path.join(D, "xatlas") + "/"
def table(f):
    L = [l.split() for l in open(X + f) if l.strip()]; return L[0], L[1:]
def cv(f):
    _, rows = table(f); return {r[0]: float(r[2]) for r in rows}
def de(f):
    _, rows = table(f); return {r[0]: (float(r[1]), float(r[2])) for r in rows}
fig = make_subplots(2, 3, specs=[[{"colspan": 1}, {}, {}], [{"colspan": 3}, None, None]], vertical_spacing=0.27, horizontal_spacing=0.09,
                    subplot_titles=["a  Recovery by panel selector", "b  DE genes recovered", "   Opposite sign", "c  Pooled and single-atlas references"])
ATL = ["K562", "HCT116", "HEK293T"]
k = N["K562"]
rec = {"K562": {s: k[f"ratio_{s}_500"] for s, _, _ in SEL if f"ratio_{s}_500" in k}}
rec["K562"]["l1000"] = k["ratio_l1000_680"]; rec["K562"]["random_lm"] = k["ratio_random1_680"]
for a, f in (("HCT116", "hct116"), ("HEK293T", "hek293t")):
    c = cv(f + "_cv.txt"); lm = cv(f + "_lm.txt")
    rec[a] = {s: c.get(f"{s}_500", np.nan) for s, _, _ in SEL}
    rec[a]["l1000"] = [v for kk, v in lm.items() if kk.startswith("l1000")][0]; rec[a]["random_lm"] = [v for kk, v in lm.items() if kk.startswith("random1")][0]
# Plotly colours, selectors keep their Fig. 2 colours
show = [("greedy", "Greedy MMSE", C[0]), ("hvgresp", "Response variance", C[2]), ("designed", "Pivoted QR", C[3]), ("module", "Program modules", C[6]),
        ("random1", "Random, 500 genes", "#FECB52"), ("l1000", "L1000 landmarks", C[1]), ("random_lm", "Random, L1000 size", "#B6E880")]
w = 0.11
for i, (s_, lab, col) in enumerate(show):
    fig.add_trace(go.Bar(x=[j + (i - 3) * w for j in range(3)], y=[rec[a].get(s_, np.nan) for a in ATL], width=w, marker_color=col, name=lab, legend="legend"), 1, 1)
fig.update_xaxes(tickvals=list(range(3)), ticktext=ATL, row=1, col=1); fig.update_yaxes(title_text="Specific response recovered", range=[0, 1], row=1, col=1)
H = N["hybrid"]; dd = {"K562": {kk: (H["K562"][kk]["recall"], H["K562"][kk]["opp"]) for kk in ("hybrid_all", "direct_all")}}
for a, f in (("HCT116", "hct116"), ("HEK293T", "hek293t")): dd[a] = de(f + "_de.txt")
for key, lab, col, off in (("hybrid_all", "Hybrid, panel measured and other genes imputed", C[2], -0.18), ("direct_all", "Measured, all genes", C[0], 0.18)):
    for c_, m in ((2, 0), (3, 1)):
        fig.add_trace(go.Bar(x=[j + off for j in range(3)], y=[dd[a][key][m] for a in ATL], width=0.34, marker_color=col, name=lab, showlegend=c_ == 2, legend="legend2"), 1, c_)
for c_ in (2, 3): fig.update_xaxes(tickvals=list(range(3)), ticktext=ATL, row=1, col=c_)
fig.update_yaxes(title_text="Fraction of DE genes, correct sign", range=[0, 0.4], row=1, col=2); fig.update_yaxes(title_text="Fraction of calls, opposite sign", range=[0, 0.35], row=1, col=3)
_, rows = table("multi_tf.txt"); _, trows = table("multi_tf_test.txt")
TG = ["RPE1", "HepG2", "Jurkat", "K562", "HCT116", "HEK293T"]
mr = {(r[0], r[2]): float(r[4]) for r in rows if r[1] == "greedy_pan"}
share = {}
for t in TG:
    singles = [(mr[(t, a)], a) for a in ATL if (t, a) in mr]; best = max(singles)[1]
    share[t] = [float(r[4]) for r in trows if r[0] == t and r[1] == best][0]
refs = [("K562", "K562 atlas", C[0]), ("HCT116", "HCT116 atlas", C[4]), ("HEK293T", "HEK293T atlas", C[5]), ("pan_equal", "Pooled across atlases", C[1])]
for i, (rf, lab, col) in enumerate(refs):
    xs = [j + (i - 1.5) * 0.2 for j in range(len(TG))]
    fig.add_trace(go.Bar(x=xs, y=[mr.get((t, rf), np.nan) for t in TG], width=0.19, marker_color=col, name=lab, legend="legend3",
                         text=[f"{share[t]:.0%}" if rf == "pan_equal" else "" for t in TG], textposition="outside", textfont=dict(size=11, color=BLACK)), 2, 1)
fig.add_vline(x=2.5, line=dict(color=GREY, dash="dash", width=1), row=2, col=1)
fig.update_xaxes(tickvals=list(range(len(TG))), ticktext=["RPE1", "HepG2", "Jurkat", "K562 (left out)", "HCT116 (left out)", "HEK293T (left out)"], row=2, col=1)
fig.update_yaxes(title_text="Specific response recovered<br>(target-derived programs)", range=[0, 0.7], row=2, col=1)
style(fig, 1300, 1000, b=40)
xa, xb, xc = fig.layout.xaxis.domain, fig.layout.xaxis2.domain, fig.layout.xaxis4.domain; ya, yc = fig.layout.yaxis.domain, fig.layout.yaxis4.domain
LG = dict(orientation="v", yanchor="top", xanchor="left", font=dict(color=BLACK, size=12), tracegroupgap=0)
fig.update_layout(barmode="overlay",
    legend=dict(LG, x=xa[0], y=ya[0] - 0.07),
    legend2=dict(LG, x=xb[0], y=ya[0] - 0.07),
    legend3=dict(LG, x=xc[0] + 0.02, y=yc[1]))
fig.for_each_trace(lambda t: t.update(marker_line=dict(color="white", width=1.5)) if t.type == "bar" else None)
save(fig, "Figure8")
# ---------------- Figure 9: measured panel versus profiles transferred from other cell lines (vcm summary, Rust)
V = {}; scope = None
for l in open(os.path.join(D, "vcm", "vcmsum.txt")):
    if l.startswith("== scope"): scope = l.split()[2]; continue
    if l.startswith("==") : scope = None; continue
    t = l.split()
    if scope and len(t) == 10 and t[0] != "target": V[(scope, t[0], t[1])] = [float(x) for x in t[3:]]
TG7 = ["RPE1", "HepG2", "Jurkat", "HCT116", "HEK293T"]
M7 = [("copy_single", "K562 profile of the same knockdown", C[5]), ("copy_pooled", "Pooled profiles of the same knockdown, other cell lines", C[4]),
      ("panel", "500 measured genes, others imputed", C[0]), ("panel+copy_pooled", "Measured panel combined with pooled profiles", C[2]), ("direct", "Direct measurement, replicate", C[3])]
fig = make_subplots(1, 3, subplot_titles=["a  Specific response, genes outside the panel", "b  DE genes found, all genes", "c  Opposite sign, all genes"], horizontal_spacing=0.08)
for i, (m, lab, col) in enumerate(M7):
    xs = [j + (i - 2) * 0.16 for j in range(len(TG7))]
    for c_, sc, k_ in ((1, "out", 0), (2, "all", 1), (3, "all", 2)):
        fig.add_trace(go.Bar(x=xs, y=[V[(sc, t, m)][k_] for t in TG7], width=0.15, marker_color=col, name=lab, showlegend=c_ == 1), 1, c_)
for c_ in (1, 2, 3): fig.update_xaxes(tickvals=list(range(len(TG7))), ticktext=TG7, row=1, col=c_)
fig.update_yaxes(title_text="Median specific r", range=[0, 0.35], row=1, col=1)
fig.update_yaxes(title_text="Fraction of DE genes, correct sign", range=[0, 0.65], row=1, col=2)
fig.update_yaxes(title_text="Fraction of calls, opposite sign", range=[0, 0.36], row=1, col=3)
style(fig, 1400, 480, b=90); fig.update_layout(barmode="overlay", legend=dict(orientation="h", y=-0.2, x=0.0))
fig.for_each_trace(lambda t: t.update(marker_line=dict(color="white", width=1)) if t.type == "bar" else None); save(fig, "Figure9")
print("figures written")
