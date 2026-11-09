"""Rust native baselines (hvgresp, module k-means) vs numpy on the small subset made by check_small.py."""
import subprocess, numpy as np, pandas as pd
O = "./out/"; S = O + "small"
subprocess.run(["./target/release/panelimp", "cv", "--z", S, "--folds", O + "folds.tsv", "--hubs", O + "hub_excl.tsv",
                "--l1000", "./ext/l1000_landmarks.txt", "--out", O + "small_base", "--sizes", "30,50", "--folds-run", "0,3", "--baselines", "--seeds", "1"], check=True, capture_output=True)
m = pd.read_csv(S + ".meta.tsv", sep="\t"); gn = pd.read_csv(S + ".genes.tsv", sep="\t").gene.values; P, G = len(m), len(gn)
z = np.fromfile(S + ".z.f32", np.float32).reshape(3, P, G).astype(np.float64)
fo = pd.read_csv(O + "folds.tsv", sep="\t").set_index("pert"); hub = set(zip(*pd.read_csv(O + "hub_excl.tsv", sep="\t").values.T))
fold = fo.fold.reindex(m.pert).fillna(-1).values
RG = pd.read_csv(O + "small_base.genes.tsv", sep="\t"); M64 = 2**64
def sm(st):
    st = (st + 0x9E3779B97F4A7C15) % M64; x = st
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) % M64; x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) % M64
    return st, x ^ (x >> 31)
def kmeans_panel(X, k):
    g, d = X.shape; xn = (X ** 2).sum(1); st = (1 * 0x9E3779B97F4A7C15) % M64
    st, r = sm(st); first = r % g; cen = [X[first]]; dmin = ((X - X[first]) ** 2).sum(1)
    for _ in range(1, k):
        st, r = sm(st); t = (r >> 11) / 2.0**53 * dmin.sum(); pick = g - 1
        for j, v in enumerate(dmin):
            if t < v: pick = j; break
            t -= v
        cen.append(X[pick]); dmin = np.minimum(dmin, ((X - X[pick]) ** 2).sum(1))
    cen = np.array(cen); asg = None
    for _ in range(300):
        new = np.argmin(xn[:, None] - 2 * X @ cen.T + (cen ** 2).sum(1)[None], 1)
        if asg is not None and (new == asg).all(): break
        asg = new
        for i in range(k):
            if (asg == i).any(): cen[i] = X[asg == i].mean(0)
    out = []
    for i in range(k):
        mem = [j for j in np.where(asg == i)[0] if j not in out]
        if mem: out.append(mem[int(np.argmin(((X[mem] - cen[i]) ** 2).sum(1)))])
    return set(out)
for f in [0, 3]:
    refs = np.array([i for i in range(P) if fold[i] != f and (m.pert[i], f) not in hub]); Y = z[0][refs]
    w, V = np.linalg.eigh(Y.T @ Y); w, V = w[::-1], V[:, ::-1]
    X = V[:, :160] * np.sqrt(np.maximum(w[:160], 0)); var = Y.var(0, ddof=1)
    for k in [30, 50]:
        exp = {f"hvgresp_{k}": set(np.lexsort((np.arange(G), -var))[:k]), f"module_{k}": kmeans_panel(X, k)}
        for name, e in exp.items():
            u = set(RG[(RG.fold == f) & (RG.panel == name)].gene); rust = {j for j in range(G) if gn[j] not in u}
            print(f"fold {f} {name}: rust {len(rust)} numpy {len(e)} identical {rust == e} overlap {len(rust & e)}")
