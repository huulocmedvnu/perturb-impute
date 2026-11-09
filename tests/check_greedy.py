"""Rust greedy MMSE panel selection vs numpy on the small subset."""
import os; os.makedirs("/tmp/panelimp_check", exist_ok=True)
import subprocess, numpy as np, pandas as pd
O = "./out/"; S = O + "small"; T = "/tmp/panelimp_check/"
m = pd.read_csv(S + ".meta.tsv", sep="\t"); G = len(pd.read_csv(S + ".genes.tsv", sep="\t"))
z = np.fromfile(S + ".z.f32", np.float32).reshape(3, len(m), G).astype(np.float64)[0]
rng = np.random.default_rng(1)
sh = rng.exponential(1, G); sh /= sh.sum() * 2; np.savetxt(T + "sh.txt", sh)
for case, nz, cap in [("unit", np.ones(G), np.inf), ("random", rng.uniform(0.2, 20, G), np.inf), ("capped", rng.uniform(0.2, 20, G), 0.05), ("tight", rng.uniform(0.2, 20, G), 0.03)]:
    np.savetxt(T + "nz.txt", nz)
    out = subprocess.run(["./target/release/panelimp", "greedytest", "--z", S, "--l1000", "./ext/l1000_landmarks.txt", "--out", T + "gt", "--noise", T + "nz.txt", "--k", "80", "--share", T + "sh.txt", "--cap", str(cap if np.isfinite(cap) else 1e300)], check=True, capture_output=True, text=True)
    robj = float(out.stderr.split()[-1]); out = out.stdout.split()
    rust = [int(x) for x in out]
    R = z.T @ z / len(z); R[np.diag_indices(G)] = np.maximum(np.diag(R) - 1, 0); sel = []; tot = 0.0; obj = 0.0
    for _ in range(80):
        gain = (R ** 2).sum(0) / (np.diag(R) + nz); gain[sel] = -np.inf; need = 80 - len(sel) - 1
        if np.isfinite(cap):
            un = np.array([j for j in np.lexsort((np.arange(G), sh)) if j not in sel])[:need + 1]; reserve = sh[un].sum(); last = sh[un[-1]] if len(un) else 0
            res = np.where(np.isin(np.arange(G), un), reserve - sh, reserve - last); gain[tot + sh + res > cap * (1 + 1e-12)] = -np.inf
        if not np.isfinite(gain.max()): break
        j = int(np.argmax(gain)); sel.append(j); tot += sh[j]; obj += gain[j]
        col = R[:, j].copy(); R -= np.outer(col, col) / (col[j] + nz[j])
    print(case, len(sel), "genes, identical order:", rust == sel, "objective rel diff", abs(robj - obj) / obj, "first diff at", next((i for i, (a, b) in enumerate(zip(rust, sel)) if a != b), None))
