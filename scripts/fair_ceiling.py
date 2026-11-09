"""Fair comparison: full-transcriptome half-0 profile denoised by projection on the top-K reference
programs (best K on test = optimistic ceiling) vs designed-500 panel imputation, on the same
non-panel genes and responding perturbations."""
import numpy as np, pandas as pd, scipy.linalg as sl
O = "./out/"
m = pd.read_csv(O + "gwps.meta.tsv", sep="\t"); G = len(pd.read_csv(O + "gwps.genes.tsv", sep="\t")); P = len(m)
Z = np.fromfile(O + "gwps.z.f32", np.float32).reshape(3, P, G)
fo = pd.read_csv(O + "folds.tsv", sep="\t"); hub = pd.read_csv(O + "hub_excl.tsv", sep="\t"); ex = set(zip(hub.pert, hub.fold))
R = pd.read_csv(O + "cv.perts.tsv", sep="\t"); R = R[R.panel == "designed_500"].set_index("pert")
fold = fo.set_index("pert").fold.reindex(m.pert).values; pidx = {p: i for i, p in enumerate(m.pert)}
Ks = [20, 50, 100, 160, 300, 500]; out = []
def pc(a, b):
    a = a - a.mean(1, keepdims=True); b = b - b.mean(1, keepdims=True)
    return (a * b).sum(1) / np.sqrt((a * a).sum(1) * (b * b).sum(1))
for f in range(5):
    refs = np.array([i for i in range(P) if fold[i] != f and (m.pert[i], f) not in ex]); tests = np.where(fold == f)[0]
    Yr = Z[0][refs].astype(np.float64)
    w, V = np.linalg.eigh(Yr.T @ Yr); V = V[:, ::-1]
    _, _, piv = sl.qr(Yr, pivoting=True, mode="economic"); S = set(piv[:500]); U = np.array([j for j in range(G) if j not in S])
    Q, _ = np.linalg.qr(V[U, :10]); res = lambda X: X - (X @ Q) @ Q.T
    h0 = Z[1][tests].astype(np.float64); h1 = Z[2][tests][:, U].astype(np.float64); t = res(h1)
    row = {"pert": m.pert.values[tests]}
    row["raw"] = pc(res(h0[:, U]), t)
    for K in Ks:
        VK = V[:, :K]; d = (h0 @ VK) @ VK.T
        row[f"K{K}"] = pc(res(d[:, U]), t)
    out.append(pd.DataFrame(row)); print("fold", f, flush=True)
D = pd.concat(out).set_index("pert").join(R[["r_ceil", "r_imp_perp", "r_ceil_perp"]])
D = D[D.r_ceil >= 0.2]
print(f"responding n={len(D)}; check raw == r_ceil_perp: max diff {np.abs(D.raw - D.r_ceil_perp).max():.2e}")
print("median specific r:", {"raw full": round(D.raw.median(), 3), **{k: round(D[k].median(), 3) for k in [f'K{K}' for K in Ks]}, "panel500": round(D.r_imp_perp.median(), 3)})
best = max([f"K{K}" for K in Ks], key=lambda k: D[k].median())
print(f"best denoised full = {best}; panel500 / denoised-full median ratio {np.median(D.r_imp_perp / D[best]):.3f}; "
      f"fraction of perts where panel >= 0.8 x denoised full {np.mean(D.r_imp_perp >= 0.8 * D[best]):.2f}")
D.to_csv(O + "fair_ceiling.tsv", sep="\t")
