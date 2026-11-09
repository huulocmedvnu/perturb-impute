"""Rust control-cell CV^2 decomposition vs numpy."""
import numpy as np, pandas as pd
C = "../coexsc/data/k562ess/"
genes = open(C + "genes.txt").read().split(); gr = pd.read_csv(C + "groups.tsv", sep="\t", header=None, names=["g", "o", "n"])
mem = np.memmap(C + "members.bin", np.uint32, "r"); ip = np.memmap(C + "indptr.bin", np.uint64, "r")
idx = np.memmap(C + "idx.bin", np.uint16, "r"); cnt = np.memmap(C + "cnt.bin", np.uint16, "r"); lib = np.memmap(C + "lib.bin", np.float32, "r")
r = gr[gr.g == "control"].iloc[0]; cells = mem[r.o:r.o + r.n]; G = len(genes)
m1 = np.zeros(G); mt = np.zeros(G); m2 = np.zeros(G); gs = np.zeros(G); L = 0.0
for c in cells:
    a, b = int(ip[c]), int(ip[c + 1]); ii = idx[a:b].astype(int); n = cnt[a:b].astype(float); sc = float(np.float32(1e4) / lib[c])
    np.add.at(m1, ii, n * sc); np.add.at(mt, ii, n * sc * sc); np.add.at(m2, ii, (n * sc) ** 2); np.add.at(gs, ii, n); L += 1e4 / sc
N = len(cells); a, t, b = m1 / N, mt / N, m2 / N
R = pd.read_csv("./out/simb.cv2.tsv", sep="\t"); gi = {g: i for i, g in enumerate(genes)}; j = R.gene.map(gi).values
with np.errstate(all="ignore"):
    tot = (b - a * a) / a ** 2; tech = t / a ** 2; bio = np.maximum(tot - tech, 0)
for k, v in [("cv2_tot", tot), ("cv2_tech", tech), ("cv2_bio", bio), ("share", gs / L)]:
    x, y = R[k].values, v[j]; ok = np.isfinite(y)
    print(k, "max rel diff", np.nanmax(np.abs(x[ok] - y[ok]) / np.maximum(np.abs(y[ok]), 1e-12)))
print("genes bio-dominated (bio > tech):", (bio[j] > tech[j]).mean().round(3), " median tech share of cv2:", np.nanmedian(tech[j] / tot[j]).round(3))
