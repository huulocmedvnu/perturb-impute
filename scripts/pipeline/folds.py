"""Grouped 5-fold split without complex leakage.
Groups = union-find over specific-profile links (corr > 0.4 after removing top 10 PCs) and CORUM
complexes (<= 50 subunits) among non-hub genes. Hubs (in > 3 complexes) are never test perturbations;
in fold f they are dropped from the reference if they share a complex with a fold-f perturbation."""
import collections, numpy as np, pandas as pd
O = "./"
m = pd.read_csv(O + "out/gwps.meta.tsv", sep="\t"); G = len(pd.read_csv(O + "out/gwps.genes.tsv", sep="\t"))
P = len(m); z = np.fromfile(O + "out/gwps.z.f32", np.float32).reshape(3, P, G)[0].astype(np.float64)
_, _, Vt = np.linalg.svd(z - z.mean(0), full_matrices=False); V = Vt[:10].T
zs = z - (z @ V) @ V.T
zc = zs - zs.mean(1, keepdims=True); zc /= np.linalg.norm(zc, axis=1, keepdims=True) + 1e-12
C = zc @ zc.T; np.fill_diagonal(C, 0)
idx = {p: i for i, p in enumerate(m.pert)}
cor = pd.read_csv(O + "ext/corum_human.txt", sep="\t", usecols=["organism", "subunits_gene_name"])
cx = [[idx[g] for g in set(x.strip() for x in s.split(";")) if g in idx] for s in cor[cor.organism == "Human"].subunits_gene_name.astype(str)]
cx = [c for c in cx if 2 <= len(c) <= 50]
deg = collections.Counter(i for c in cx for i in c); hub = np.array([deg[i] > 3 for i in range(P)])
par = list(range(P))
def f(i):
    while par[i] != i: par[i] = par[par[i]]; i = par[i]
    return i
for a, b in zip(*np.where(np.triu(C) > 0.4)):
    if not hub[a] and not hub[b]: par[f(a)] = f(b)
for c in cx:
    c = [i for i in c if not hub[i]]
    for a in c[1:]: par[f(c[0])] = f(a)
root = np.array([f(i) for i in range(P)])
sz = pd.Series(root[~hub]).value_counts()
fold = np.full(P, -1); load = np.zeros(5)
for r in sz.index:
    k = int(np.argmin(load)); fold[(root == r) & ~hub] = k; load[k] += sz[r]
print(f"hubs {hub.sum()}; groups {len(sz)}; largest {sz.iloc[:4].tolist()}; fold sizes {load.astype(int).tolist()}")
# hub exclusion per fold: complex partners in fold, or strong specific-profile link to a fold member
rows = set()
for c in cx:
    hs = [i for i in c if hub[i]]; fs = {fold[i] for i in c if not hub[i]}
    for h in hs:
        for k in fs: rows.add((m.pert[h], k))
for a, b in zip(*np.where(C > 0.4)):
    if hub[a] and not hub[b]: rows.add((m.pert[a], fold[b]))
pd.DataFrame({"pert": m.pert, "group": root, "fold": fold}).to_csv(O + "out/folds.tsv", sep="\t", index=False)
pd.DataFrame(sorted(rows), columns=["pert", "fold"]).to_csv(O + "out/hub_excl.tsv", sep="\t", index=False)
print("hub exclusions", len(rows), "per fold", pd.DataFrame(sorted(rows), columns=["p", "k"]).k.value_counts().sort_index().tolist())
