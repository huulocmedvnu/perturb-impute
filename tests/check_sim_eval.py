"""Rust sim evaluation (fold 0, q = 1, deterministic) vs an independent numpy/scipy implementation."""
import numpy as np, pandas as pd, scipy.linalg as sl, scipy.sparse as sp, time
O = "./out/"; C = "../coexsc/data/k562ess/"; t0 = time.time()
m = pd.read_csv(O + "gwpsE.meta.tsv", sep="\t"); gn = pd.read_csv(O + "gwpsE.genes.tsv", sep="\t").gene.values
P, G = len(m), len(gn); Z = np.fromfile(O + "gwpsE.z.f32", np.float32).reshape(3, P, G).astype(np.float64)
fo = pd.read_csv(O + "folds.tsv", sep="\t").set_index("pert"); hub = set(zip(*pd.read_csv(O + "hub_excl.tsv", sep="\t").values.T))
fold = fo.fold.reindex(m.pert).fillna(-1).values.astype(int); group = fo.group.reindex(m.pert).fillna(-1).values.astype(int)
f = 0; refs = np.array([i for i in range(P) if fold[i] != f and (m.pert[i], f) not in hub])
key = lambda i: ((int(group[i] if group[i] >= 0 else 1_000_000 + i) * 0x9E3779B97F4A7C15) % 2**64 >> 32) % 10
itr = np.array([i for i in refs if key(i) != 0]); iho = np.array([i for i in refs if key(i) == 0])
Y = Z[0][refs]; Cm = Y.T @ Y; Ci = Z[0][itr].T @ Z[0][itr]
w, V = np.linalg.eigh(Cm); V = V[:, ::-1][:, :10]
_, _, piv = sl.qr(Y, pivoting=True, mode="economic")
def perm(seed, g):
    st = (seed * 0x9E3779B97F4A7C15) % 2**64; idx = list(range(g)); M = 2**64
    for i in range(g - 1, 0, -1):
        st = (st + 0x9E3779B97F4A7C15) % M; x = st
        x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) % M; x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) % M; x ^= x >> 31
        j = x % (i + 1); idx[i], idx[j] = idx[j], idx[i]
    return idx
pan = {"designed_500": list(piv[:500]), "random1_500": perm(1, G)[:500]}
used = set().union(*map(set, pan.values())); u = np.array([j for j in range(G) if j not in used])
Q, _ = np.linalg.qr(V[u]); res = lambda x: x - Q @ (Q.T @ x)
print("setup", round(time.time() - t0, 1), "common", len(u))
def ridge(Cx, s, al):
    A = Cx[np.ix_(s, s)].copy(); A += al * np.trace(A) / len(s) * np.eye(len(s)); return np.linalg.solve(A, Cx[np.ix_(s, u)])
# cache pseudobulk at full depth
genes = open(C + "genes.txt").read().split(); cgi = {g: i for i, g in enumerate(genes)}
gr = pd.read_csv(C + "groups.tsv", sep="\t", header=None, names=["g", "o", "n"])
mem = np.memmap(C + "members.bin", np.uint32, "r"); ip = np.memmap(C + "indptr.bin", np.uint64, "r")
idx = np.memmap(C + "idx.bin", np.uint16, "r"); cnt = np.memmap(C + "cnt.bin", np.uint16, "r"); lib = np.memmap(C + "lib.bin", np.float32, "r")
pf = dict(zip(m.pert, fold))
q = gr[(gr.g != "control") & (gr.n >= 4) & gr.g.map(lambda x: pf.get(x, -99) == f)].reset_index(drop=True)
grp = pd.concat([gr[gr.g == "control"], q]).reset_index(drop=True)
cols = np.array([cgi[g] for g in gn]); cmap = -np.ones(len(genes), int); cmap[cols] = np.arange(G)
def sums(colmap, ncol, targeted):
    S = np.zeros((len(grp), 2, ncol))
    for s, r in grp.iterrows():
        cells = mem[r.o:r.o + r.n]
        for h in (0, 1):
            for cell in cells[h::2]:
                a, b = int(ip[cell]), int(ip[cell + 1]); ii = idx[a:b].astype(int); cc = cnt[a:b].astype(float); k = colmap[ii]
                if targeted:
                    L = cc[k >= 0].sum(); sc = 1e4 / L if L > 0 else 0.0
                else: sc = float(np.float32(1e4) / lib[cell])
                np.add.at(S[s, h], k[k >= 0], cc[k >= 0] * sc)
    return S
def zh(S):
    n = grp.n.values; nc = n[0]; mc = S[0].sum(0) / nc; n0, n1 = (n + 1) // 2, n // 2
    m0 = S[1:, 0] / n0[1:, None]; m1 = S[1:, 1] / n1[1:, None]
    s2 = np.median((m0 - m1) ** 2 / (1 / n0[1:, None] + 1 / n1[1:, None]), 0) / 0.454936423119572
    with np.errstate(all="ignore"): z0 = np.where(s2 > 0, (m0 - mc) / np.sqrt(s2 * (1 / n0[1:, None] + 1 / nc)), 0)
    return z0
Sw = sums(cmap, G, False); zw = zh(Sw); print("wt pseudobulk", round(time.time() - t0, 1))
# truth half 1 = coexsc ess.z layer 2 (validated == Rust truth to 4e-6)
em = pd.read_csv("../coexsc/out/ess.meta.tsv", sep="\t"); eg = pd.read_csv("../coexsc/out/ess.genes.tsv", sep="\t")
EZ = np.fromfile("../coexsc/out/ess.z.f32", np.float32).reshape(3, len(em), len(eg)); epi = {p: i for i, p in enumerate(em.pert)}; egi = {g: i for i, g in enumerate(eg.gene)}
truth = np.zeros((len(q), len(u)))
for j, gj in enumerate(u):
    if gn[gj] in egi: truth[:, j] = EZ[2][[epi[x] for x in q.g], egi[gn[gj]]]
R = pd.read_csv(O + "simv.sim.tsv", sep="\t"); worst = {}
def met(v, t):
    nde = int((np.abs(t) > 4).sum()); o = np.argsort(-np.abs(v), kind="stable")[:nde]; same = np.sign(v[o]) == np.sign(t[o])
    return np.corrcoef(v, t)[0, 1], np.corrcoef(res(v), res(t))[0, 1], nde, int((same & (np.abs(t[o]) > 4)).sum()), int((~same).sum())
dir_m = [met(zw[s, u], truth[s]) for s in range(len(q))]
for name, s in pan.items():
    best = None
    for al in [1e-3, 1e-2, 1e-1, 1, 10]:
        W = ridge(Ci, s, al); im = Z[1][np.ix_(iho, s)] @ W; tr = Z[2][np.ix_(iho, u)]
        r = np.sum([x for x in [np.corrcoef(im[i], tr[i])[0, 1] for i in range(len(iho))] if np.isfinite(x)]) / len(iho)
        if best is None or r > best[0]: best = (r, al)
    W = ridge(Cm, s, best[1])
    cm2 = -np.ones(len(genes), int); cm2[cols[s]] = np.arange(len(s)); zt = zh(sums(cm2, len(s), True))
    for assay, X in [("extract", zw[:, s]), ("targeted", zt)]:
        rr = R[(R.panel == name) & (R.assay == assay)].set_index("pert"); im = X @ W
        worst["alpha"] = max(worst.get("alpha", 0), abs(rr.alpha.iloc[0] - best[1]))
        for si, pn in enumerate(q.g):
            a = met(im[si], truth[si]); d = dir_m[si]; row = rr.loc[pn]
            for k, x, y in [("r_imp", a[0], row.r_imp), ("rp_imp", a[1], row.rp_imp), ("r_dir", d[0], row.r_dir), ("rp_dir", d[1], row.rp_dir)]:
                if np.isfinite(x): worst[k] = max(worst.get(k, 0), abs(x - y))
            for k, x, y in [("n_de", a[2], row.n_de), ("tp_imp", a[3], row.tp_imp), ("disc_imp", a[4], row.disc_imp), ("tp_dir", d[3], row.tp_dir), ("disc_dir", d[4], row.disc_dir)]:
                worst[k + "_mismatch"] = worst.get(k + "_mismatch", 0) + int(x != y)
print("rows", len(q), "max abs diff / count mismatches Rust vs numpy:", {k: (f"{v:.2e}" if isinstance(v, float) else v) for k, v in worst.items()}, round(time.time() - t0, 1), "s")
