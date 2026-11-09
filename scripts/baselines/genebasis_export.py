"""Export reference z matrices (layer 0) per fold for geneBasisR, as float32 binary genes x perts (column-major
for R: we write the (n_ref, G) C-order matrix, which R reads as a G x n_ref matrix). Also writes row/col names.
Usage: python3 genebasis_export.py [--sub N]   (N = optional row subsample, seed 0; default none)"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load, ref_matrix, OUT
sub = int(sys.argv[sys.argv.index("--sub") + 1]) if "--sub" in sys.argv else 0
perts, genes, z, fold, excl = load()
d = f"{OUT}/genebasis_in"; os.makedirs(d, exist_ok=True)
open(f"{d}/genes.txt", "w").write("\n".join(genes) + "\n")
for f in [0, 1, 2, 3, 4, -9]:
    idx, X = ref_matrix(z, fold, excl, f)
    if sub and len(idx) > sub:
        keep = np.sort(np.random.default_rng(0).choice(len(idx), sub, replace=False)); idx, X = idx[keep], X[keep]
    X.astype("<f4").tofile(f"{d}/fold{f}.f32")
    open(f"{d}/fold{f}.perts.txt", "w").write("\n".join(perts[i] for i in idx) + "\n")
    print("fold", f, X.shape)
