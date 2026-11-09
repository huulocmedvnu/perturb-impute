"""Rust panelimp vs an independent numpy/scipy implementation on a subset (600 genes, 1,500 perts)."""
import os, subprocess, numpy as np, pandas as pd, scipy.linalg as sl
O = "./out/"; S = O + "small"
m = pd.read_csv(O + "gwps.meta.tsv", sep="\t"); gn = pd.read_csv(O + "gwps.genes.tsv", sep="\t")
P, G = len(m), len(gn); Z = np.fromfile(O + "gwps.z.f32", np.float32).reshape(3, P, G)
fo = pd.read_csv(O + "folds.tsv", sep="\t"); hub = pd.read_csv(O + "hub_excl.tsv", sep="\t")
rng = np.random.default_rng(0); pi = np.sort(rng.choice(P, 1500, replace=False)); gi = np.sort(rng.choice(G, 600, replace=False))
L = set(open("./ext/l1000_landmarks.txt").read().split())
Z[:, pi][:, :, gi].astype(np.float32).tofile(S + ".z.f32")
m.iloc[pi].to_csv(S + ".meta.tsv", sep="\t", index=False); gn.iloc[gi].to_csv(S + ".genes.tsv", sep="\t", index=False)
nlm = sum(g in L for g in gn.gene.iloc[gi]); print("landmarks in subset", nlm)
subprocess.run(["./target/release/panelimp", "cv", "--z", S, "--folds", O + "folds.tsv", "--hubs", O + "hub_excl.tsv",
                "--l1000", "./ext/l1000_landmarks.txt", "--out", O + "small_rust", "--sizes", f"50,{nlm}", "--folds-run", "0,3"] + (["--common"] if os.environ.get("COMMON") else []), check=True)
R = pd.read_csv(O + "small_rust.perts.tsv", sep="\t"); RD = pd.read_csv(O + "small_rust.de.tsv", sep="\t"); RK = pd.read_csv(O + "small_rust.knn.tsv", sep="\t"); RG = pd.read_csv(O + "small_rust.genes.tsv", sep="\t")
z = Z[:, pi][:, :, gi].astype(np.float64); perts = m.pert.values[pi]; genes = gn.gene.values[gi]
fmap = fo.set_index("pert"); fold = fmap.fold.reindex(perts).values; group = fmap.group.reindex(perts).values
ex = set(zip(hub.pert, hub.fold))
def key(i):
    k = group[i] if group[i] >= 0 else 1_000_000 + i
    return ((int(k) * 0x9E3779B97F4A7C15) % 2**64 >> 32) % 10
def splitmix_perm(seed, g):
    st = (seed * 0x9E3779B97F4A7C15) % 2**64; idx = list(range(g)); M64 = 2**64
    for i in range(g - 1, 0, -1):
        st = (st + 0x9E3779B97F4A7C15) % M64; x = st
        x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) % M64; x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) % M64; x ^= x >> 31
        j = x % (i + 1); idx[i], idx[j] = idx[j], idx[i]
    return idx
def pear(a, b): return np.corrcoef(a, b)[0, 1]
def auc(s, l):
    if l.sum() < 5 or (~l).sum() < 5: return np.nan
    from scipy.stats import rankdata
    r = rankdata(s); n1 = l.sum(); return (r[l].sum() - n1 * (n1 + 1) / 2) / (n1 * (~l).sum())
def ridge(refs, s, u, al):
    X = z[0][np.ix_(refs, s)]; Y = z[0][np.ix_(refs, u)]; A = X.T @ X
    A += al * np.trace(A) / len(s) * np.eye(len(s)); return np.linalg.solve(A, X.T @ Y)
worst = {}
def upd(k, v): worst[k] = max(worst.get(k, 0), v)
for f in [0, 3]:
    refs = np.array([i for i in range(len(perts)) if fold[i] != f and (perts[i], f) not in ex]); tests = np.where(fold == f)[0]
    itr = np.array([i for i in refs if key(i) != 0]); iho = np.array([i for i in refs if key(i) == 0])
    Yr = z[0][refs]; w, Vv = np.linalg.eigh(Yr.T @ Yr); V = Vv[:, ::-1][:, :10]
    _, _, piv = sl.qr(Yr, pivoting=True, mode="economic")
    lm = [j for j in range(len(genes)) if genes[j] in L]
    pan = {}
    for k in [50, nlm]:
        pan[f"designed_{k}"] = list(piv[:k])
        for sd in (1, 2, 3): pan[f"random{sd}_{k}"] = splitmix_perm(sd, len(genes))[:k]
        if k == len(lm): pan[f"l1000_{k}"] = lm
    E = [j for j in range(len(genes)) if j not in set().union(*map(set, pan.values()))]
    for name, s in pan.items():
        u = E if os.environ.get("COMMON") else [j for j in range(len(genes)) if j not in set(s)]
        best = None
        for al in [1e-3, 1e-2, 1e-1, 1, 10]:
            W = ridge(itr, s, u, al); im = z[1][np.ix_(iho, s)] @ W; tr = z[2][np.ix_(iho, u)]
            r = np.nanmean([pear(im[i], tr[i]) for i in range(len(iho))]) * 1.0
            r = np.sum([x for x in [pear(im[i], tr[i]) for i in range(len(iho))] if np.isfinite(x)]) / len(iho)
            if best is None or r > best[0]: best = (r, al, im, tr, W)
        _, al, imh, trh, Win = best
        info = np.array([pear(imh[:, j], trh[:, j]) for j in range(len(u))])
        W = ridge(refs, s, u, al); im = z[1][np.ix_(tests, s)] @ W; tr = z[2][np.ix_(tests, u)]; dr = z[1][np.ix_(tests, u)]
        real = np.array([pear(im[:, j], tr[:, j]) for j in range(len(u))])
        rg = RG[(RG.fold == f) & (RG.panel == name)]
        upd("info", np.nanmax(np.abs(rg["info"].values - info))); upd("real", np.nanmax(np.abs(rg["real"].values - real)))
        Q, _ = np.linalg.qr(V[u]); res = lambda x: x - Q @ (Q.T @ x)
        Xr = z[0][np.ix_(refs, s)]; ew, ev = np.linalg.eigh(Xr.T @ Xr); kk = min(160, len(s) // 2); PB = ev[:, ::-1][:, :kk]
        rr = R[(R.fold == f) & (R.panel == name)].set_index("pert")
        upd("alpha", abs(rr.alpha.iloc[0] - al))
        # DE preservation: calibration on inner holdout, then calls / signs / ranks on tests
        from scipy.stats import norm, spearmanr
        bb = (imh * trh).sum(0) / (imh ** 2).sum(0); ee = np.maximum(((trh - bb * imh) ** 2).mean(0) - 1, 0.05)
        dd = RD[(RD.fold == f) & (RD.panel == name)].set_index("pert")
        for i, t in enumerate(tests):
            y, ti, di = im[i], tr[i], dr[i]; c = bb * y; de = np.abs(ti) > 4
            def calls(v, sel):
                same = np.sign(v) == np.sign(ti); return sel.sum(), (sel & same & de).sum(), (sel & ~same).sum()
            cs = norm.cdf(-np.abs(c) / np.sqrt(ee)) < 0.01
            vals = dict(zip(["m_call", "m_tp", "m_disc"], calls(di, np.abs(di) > 4)))
            vals.update(zip(["r_call", "r_tp", "r_disc"], calls(y, np.abs(y) > 4)))
            vals.update(zip(["c_call", "c_tp", "c_disc"], calls(c, cs)))
            top = np.zeros(len(c), bool); top[np.argsort(-np.abs(c), kind="stable")[:int((np.abs(di) > 4).sum())]] = True
            _, vals["k_tp"], vals["k_disc"] = calls(c, top)
            vals["n_de"] = de.sum(); vals["c_exp_disc"] = norm.cdf(-np.abs(c[cs]) / np.sqrt(ee[cs] + 1)).sum()
            for k3, v in [("m", di), ("r", y), ("c", c)]:
                vals[k3 + "_sxy"] = v @ ti; vals[k3 + "_sxx"] = v @ v
                if de.sum() >= 5:
                    vals["sp_" + k3] = spearmanr(v[de], ti[de])[0]; vals["spa_" + k3] = spearmanr(np.abs(v[de]), np.abs(ti[de]))[0]
            row = dd.loc[perts[t]]
            for k2, v2 in vals.items():
                rel = abs(row[k2] - v2) / max(1.0, abs(v2)); upd("de_" + k2.split("_")[0] + ("_int" if "call" in k2 or "tp" in k2 or k2.startswith("k_") or "disc" in k2 and "exp" not in k2 else ""), rel)
        for i, t in enumerate(tests):
            x = z[1][t, s]; oor = 1 - ((PB.T @ x) ** 2).sum() / (x @ x)
            vals = dict(r_ceil=pear(dr[i], tr[i]), r_imp=pear(im[i], tr[i]), r_ceil_perp=pear(res(dr[i]), res(tr[i])), r_imp_perp=pear(res(im[i]), res(tr[i])),
                        auc_ceil=auc(np.abs(dr[i]), np.abs(tr[i]) > 4), auc_imp=auc(np.abs(im[i]), np.abs(tr[i]) > 4), oor=oor)
            row = rr.loc[perts[t]]
            for k2, v2 in vals.items():
                if np.isfinite(v2): upd(k2, abs(row[k2] - v2))
                elif np.isfinite(row[k2]): upd(k2 + "_nan_mismatch", 1)
        # step 4: knn / hybrid
        def nbrs(xq, xp):
            st = lambda x: (lambda c: c / np.maximum(np.linalg.norm(c, axis=1, keepdims=True), 1e-300))(x - x.mean(1, keepdims=True))
            C = st(xq) @ st(xp).T; km = min(100, len(xp)); out = []
            for q in range(len(xq)):
                o = np.lexsort((np.arange(len(xp)), -C[q]))[:km]; pn = (xp[o] ** 2).sum(1)
                out.append((o, C[q, o], np.where(pn > 0, xp[o] @ xq[q] / pn, 0)))
            return out
        def kpred(nb, base, pool, k, pw):
            out = np.zeros((len(nb), pool.shape[1])) if base is None else base.copy()
            for q, (o, r, sc) in enumerate(nb):
                w = np.maximum(r[:k], 0) ** pw
                if w.sum() > 0: out[q] += (w * sc[:k]) @ pool[o[:k]] / w.sum()
            return out
        def spa(v, t):
            d = np.abs(t) > 4
            return spearmanr(np.abs(v[d]), np.abs(t[d]))[0] if d.sum() >= 5 else np.nan
        def pools(rows, Wx):
            xp = z[0][np.ix_(rows, s)]; yp = z[0][np.ix_(rows, u)]; return xp, yp, yp - xp @ Wx
        xp, yp, rp = pools(itr, Win); nb = nbrs(z[1][np.ix_(iho, s)], xp); bk = {0: (-np.inf,), 1: (-np.inf,)}
        for kk in [5, 10, 20, 50, 100]:
            for pw in [1.0, 2.0, 4.0]:
                for md in (0, 1):
                    pr = kpred(nb, None, yp, kk, pw) if md == 0 else kpred(nb, imh, rp, kk, pw)
                    ob = np.nanmedian([spa(pr[i], trh[i]) for i in range(len(iho))])
                    if ob > bk[md][0]: bk[md] = (ob, kk, pw)
        xp, yp, rp = pools(refs, W); nb = nbrs(z[1][np.ix_(tests, s)], xp)
        pk = kpred(nb, None, yp, bk[0][1], bk[0][2]); ph = kpred(nb, im, rp, bk[1][1], bk[1][2])
        kd = RK[(RK.fold == f) & (RK.panel == name)].set_index("pert")
        upd("knn_param", abs(kd.k_knn.iloc[0] - bk[0][1]) + abs(kd.p_knn.iloc[0] - bk[0][2]) + abs(kd.k_hyb.iloc[0] - bk[1][1]) + abs(kd.p_hyb.iloc[0] - bk[1][2]))
        upd("knn_obj", abs(kd.obj_knn.iloc[0] - bk[0][0]) + abs(kd.obj_hyb.iloc[0] - bk[1][0]))
        for i, t in enumerate(tests):
            row = kd.loc[perts[t]]; ti = tr[i]; nm = int((np.abs(dr[i]) > 4).sum()); dmask = np.abs(ti) > 4
            for mn, v in [("ridge", im[i]), ("knn", pk[i]), ("hyb", ph[i])]:
                o = np.argsort(-np.abs(v), kind="stable")[:nm]; same = np.sign(v[o]) == np.sign(ti[o])
                vals = {"r": pear(v, ti), "rperp": pear(res(v), res(ti)), "spa": spa(v, ti),
                        "sp": spearmanr(v[dmask], ti[dmask])[0] if dmask.sum() >= 5 else np.nan, "ktp": (same & dmask[o]).sum(), "kdisc": (~same).sum()}
                for k2, v2 in vals.items():
                    rv = row[f"{mn}_{k2}"]
                    if np.isfinite(v2): upd("knn_" + k2, abs(rv - v2))
                    elif np.isfinite(rv): upd("knn_" + k2 + "_nan", 1)
print("max abs diff Rust vs numpy:", {k: f"{v:.2e}" for k, v in worst.items()})
