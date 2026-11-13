"""Collect the analysis of perturbations with few cells (panelimp rescue, scripts/pipeline/rescue.sh) into
source_data/: per-perturbation tables (tables/rescue.<screen>.tsv.gz), the summary by stratum of input cells
(rescue/summary.tsv), the reference models of each fold (rescue/reference.tsv) and Table S10 (rescue/table_s10.md)."""
import re, pandas as pd
O = "out/rescue/"; SD = "source_data/"
SCREENS = [("K562", "k562"), ("RPE1", "rpe1"), ("HepG2", "hepg2"), ("Jurkat", "jurkat"), ("K562 essential", "k562ess")]
summ, refs = [], []
for name, c in SCREENS:
    t = pd.read_csv(f"{O}{c}.rescue.tsv", sep="\t")
    t.to_csv(f"{SD}tables/rescue.{c}.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    summ.append(pd.read_csv(f"{O}{c}.rescue.summary.tsv", sep="\t").assign(screen=name))
    kind = ["external", "internal"]; i = 0
    for l in open(f"{O}{c}.rescue.log"):
        m = re.match(r"fold (\d), input cells (\S+): reference (\d+) perturbations, penalty (\S+), programs (\d+), signal variance per cell outside the programs (\S+)", l)
        refs.append((name, int(m[1]), m[2], kind[i % 2], int(m[3]), float(m[4]), int(m[5]), float(m[6]))); i += 1
S = pd.concat(summ); S = S[["screen"] + [c for c in S.columns if c != "screen"]]
S.to_csv(f"{SD}rescue/summary.tsv", sep="\t", index=False)
pd.DataFrame(refs, columns=["screen", "fold", "input_cells", "reference", "n_reference", "penalty", "programs", "signal_variance_per_cell"]).to_csv(f"{SD}rescue/reference.tsv", sep="\t", index=False)

LAB = {"direct": "Direct measurement", "hybrid, external": "Hybrid profile", "hybrid + prior, external": "Hybrid profile with guide-target prior", "projected, external": "Projected whole transcriptome",
       "shrunk, external": "Shrunk whole transcriptome", "hybrid, internal": "Hybrid profile", "hybrid + prior, internal": "Hybrid profile with guide-target prior",
       "projected, internal": "Projected whole transcriptome", "shrunk, internal": "Shrunk whole transcriptome"}
def rows(sel, first):
    out = []
    for _, r in sel.iterrows():
        ref = "None" if r.readout == "direct" else ("Atlas" if r.readout.endswith("external") else "Internal")
        out.append(f"| {first} | {int(r.n)} | {ref} | {LAB[r.readout]} | {r.median_rp_own:.3f} | {r.de_recall:.3f} | {r.opposite:.3f} | {r.pathway_r:.3f} | {r.complex_auroc:.3f} |")
    return out
head = ["| {} | n | Reference | Readout | Specific r | DE recall | Opposite sign | Pathway r | Complex AUC |", "|---------|---------|--------|---------------------------------|--------|-------|--------|--------|--------|"]
md = ["### a. Input half restricted to 15, 30 and 50 cells, K562 essential-gene screen", "", head[0].format("Input cells"), head[1]]
e = S[S.screen == "K562 essential"]
for lev in ["15", "30", "50", "all"]:
    md += rows(e[(e.input_cells.astype(str) == lev) & (e.stratum == "all")], "All" if lev == "all" else lev)
md += ["", "### b. Knockdowns with 50 or fewer input cells in four screens", "", head[0].format("Cell line"), head[1]]
for name, _ in SCREENS[:4]:
    x = S[(S.screen == name) & (S.input_cells.astype(str) == "all")]
    md += rows(x[x.stratum == "<=50"], name)
open(f"{SD}rescue/table_s10.md", "w").write("\n".join(md) + "\n")
print(len(S), "summary rows;", len(md), "table lines")
