"""scanpy HVG baseline: control (non-targeting) K562 essential-screen cells, raw UMI counts,
genes restricted to the 7,094 panelimp universe (gwps.genes.tsv), then
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=k) separately for each k.
Output ext/baselines/hvg.tsv (fold = all). Usage: venv/bin/python hvg_scanpy.py"""
import os, sys, time, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc, skmisc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT, OUT
D = os.path.expanduser("~/Downloads/coexsc/data/k562ess")
SIZES = [200, 500, 680]

cache_genes = [l.rstrip("\n") for l in open(f"{D}/genes.txt")]
grp = pd.read_csv(f"{D}/groups.tsv", sep="\t", header=None, names=["name", "start", "count"])
g = grp[grp.name == "control"].iloc[0]
members = np.memmap(f"{D}/members.bin", dtype=np.uint32, mode="r")
cells = np.asarray(members[g.start:g.start + g["count"]], dtype=np.int64)
indptr = np.memmap(f"{D}/indptr.bin", dtype=np.uint64, mode="r")
idx = np.memmap(f"{D}/idx.bin", dtype=np.uint16, mode="r")
cnt = np.memmap(f"{D}/cnt.bin", dtype=np.uint16, mode="r")
print(f"control cells {len(cells)}; cache genes {len(cache_genes)}; total cells {len(indptr) - 1}")
rows_ptr, cols, vals = [0], [], []
for c in cells:
    a, b = int(indptr[c]), int(indptr[c + 1])
    cols.append(np.asarray(idx[a:b], dtype=np.int32)); vals.append(np.asarray(cnt[a:b], dtype=np.float32))
    rows_ptr.append(rows_ptr[-1] + (b - a))
X = sp.csr_matrix((np.concatenate(vals), np.concatenate(cols), np.array(rows_ptr)), shape=(len(cells), len(cache_genes)))
X.sum_duplicates()

universe = pd.read_csv(f"{ROOT}/out/gwps.genes.tsv", sep="\t")["gene"].astype(str).tolist()
assert len(set(cache_genes)) == len(cache_genes), "duplicate symbols in cache"
gi = {s: i for i, s in enumerate(cache_genes)}
found = [s for s in universe if s in gi]
print(f"universe genes {len(universe)}; found in single-cell cache {len(found)}; missing {len(universe) - len(found)}")
adata = ad.AnnData(X[:, [gi[s] for s in found]].tocsr(), obs=pd.DataFrame(index=[f"c{c}" for c in cells]),
                   var=pd.DataFrame(index=found))
print(adata, "\nUMI/cell median", float(np.median(np.asarray(adata.X.sum(1)).ravel())))
print("scanpy", sc.__version__, "anndata", ad.__version__, "scikit-misc", skmisc.__version__)
rows, ranks = [], {}
for k in SIZES:
    a = adata.copy(); t0 = time.time()
    sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=k)
    sel = a.var.index[a.var.highly_variable].tolist()
    assert len(sel) == k
    sel = sorted(sel, key=lambda s: a.var.loc[s, "highly_variable_rank"])   # rank order
    ranks[k] = sel
    print(f"k={k}: {len(sel)} genes in {time.time() - t0:.2f}s; top5 {sel[:5]}")
    rows += [f"all\t{k}\t{s}" for s in sel]
for a_, b_ in [(200, 500), (500, 680), (200, 680)]:
    print(f"nested {a_} in {b_}: {set(ranks[a_]) <= set(ranks[b_])}")
os.makedirs(OUT, exist_ok=True)
open(f"{OUT}/hvg.tsv", "w").write("fold\tsize\tgene\n" + "\n".join(rows) + "\n")
print("wrote", f"{OUT}/hvg.tsv")
