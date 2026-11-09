"""Shared loader for baseline panel selection: reference sets identical to src/main.rs (line ~486).
ref(f) = {i : fold[i] != f and (i, f) not in hub_excl}; transfer (f = -9) = all perturbations."""
import os, numpy as np, pandas as pd
ROOT = os.path.expanduser("~/Downloads/panelimp")
OUT = os.path.join(ROOT, "ext/baselines")

def load():
    perts = pd.read_csv(f"{ROOT}/out/gwps.meta.tsv", sep="\t")["pert"].astype(str).tolist()
    genes = pd.read_csv(f"{ROOT}/out/gwps.genes.tsv", sep="\t")["gene"].astype(str).tolist()
    P, G = len(perts), len(genes)
    z = np.memmap(f"{ROOT}/out/gwps.z.f32", dtype=np.float32, mode="r", shape=(3, P, G))
    pi = {p: i for i, p in enumerate(perts)}
    fold = np.full(P, -1, dtype=np.int64)
    fd = pd.read_csv(f"{ROOT}/out/folds.tsv", sep="\t")
    for p, f in zip(fd.pert.astype(str), fd.fold):
        if p in pi: fold[pi[p]] = int(f)
    excl = set()
    hx = pd.read_csv(f"{ROOT}/out/hub_excl.tsv", sep="\t")
    for p, f in zip(hx.pert.astype(str), hx.fold):
        if p in pi: excl.add((pi[p], int(f)))
    return perts, genes, z, fold, excl

def ref_idx(fold, excl, f):
    if f == -9: return np.arange(len(fold))
    return np.array([i for i in range(len(fold)) if fold[i] != f and (i, f) not in excl])

def ref_matrix(z, fold, excl, f):
    idx = ref_idx(fold, excl, f)
    return idx, np.ascontiguousarray(z[0][idx], dtype=np.float32)

if __name__ == "__main__":
    perts, genes, z, fold, excl = load()
    print("P", len(perts), "G", len(genes), "fold counts", np.unique(fold, return_counts=True), "excl", len(excl))
    for f in [0, 1, 2, 3, 4, -9]:
        idx = ref_idx(fold, excl, f); print("fold", f, "n_ref", len(idx))
    X = z[0]; print("layer0 nan", np.isnan(X[:200]).sum(), "mean", float(X[:200].mean()), "sd", float(X[:200].std()))
