"""Simulation control: thinned whole transcriptome projected on 160 reference programs (out/simd.den.tsv),
compared with targeted panels (out/simd.sim.tsv), folds 1-4. Writes docs/manuscript/source_sim/denoised.tsv."""
import os, pandas as pd
from scipy.stats import wilcoxon
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
d = pd.read_csv(f"{R}/out/simd.den.tsv", sep="\t"); d = d[d.fold > 0]
s = pd.read_csv(f"{R}/out/simd.sim.tsv", sep="\t"); s = s[s.fold > 0]
rows = []
for q, g in d.groupby("q"):
    rows.append((q, "direct", g.rp_dir.median(), g.tp_dir.sum() / g.n_de.sum(), g.disc_dir.sum() / g.n_de.sum(), "", ""))
    rows.append((q, "denoised", g.rp_den.median(), g.tp_den.sum() / g.n_de.sum(), g.disc_den.sum() / g.n_de.sum(), "", ""))
    for p, a in [("greedyfull_500", "targeted"), ("cost0.5_500", "targeted_het"), ("l1000_679", "targeted"), ("random1_500", "targeted")]:
        h = s[(s.q == q) & (s.panel == p) & (s.assay == a)].set_index("pert"); x = h.rp_imp - g.set_index("pert").rp_den.reindex(h.index)
        rows.append((q, f"{p}|{a}", h.rp_imp.median(), h.tp_imp.sum() / h.n_de.sum(), h.disc_imp.sum() / h.n_de.sum(), (x > 0).mean(), wilcoxon(x).pvalue))
pd.DataFrame(rows, columns=["q", "readout", "specific_r", "de_found", "opposite", "share_above_denoised", "p_vs_denoised"]).to_csv(f"{R}/docs/manuscript/source_sim/denoised.tsv", sep="\t", index=False, float_format="%.6g")
print(d.pert.nunique(), "perturbations")
