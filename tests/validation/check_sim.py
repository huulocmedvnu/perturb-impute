"""Rust sim (pbdump thinning + truth z) vs an independent numpy implementation."""
import os; os.makedirs("/tmp/panelimp_check", exist_ok=True)
import subprocess, numpy as np, pandas as pd
C = "../coexsc/data/k562ess/"; B = "./target/release/panelimp"
T = "/tmp/panelimp_check/"
genes = open(C + "genes.txt").read().split(); gr = pd.read_csv(C + "groups.tsv", sep="\t", header=None, names=["g", "o", "n"])
mem = np.memmap(C + "members.bin", np.uint32, "r"); ip = np.memmap(C + "indptr.bin", np.uint64, "r")
idx = np.memmap(C + "idx.bin", np.uint16, "r"); cnt = np.memmap(C + "cnt.bin", np.uint16, "r"); lib = np.memmap(C + "lib.bin", np.float32, "r")
M = 2**64
def sm(st):
    st = (st + 0x9E3779B97F4A7C15) % M; x = st
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) % M; x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) % M
    return st, x ^ (x >> 31)
def pb(names, q, seed, panel=None, rate=None):
    out = {}
    pidx = {genes.index(g): k for k, g in enumerate(panel)} if panel else None
    for n in names:
        r = gr[gr.g == n].iloc[0]; cells = mem[r.o:r.o + r.n]
        for h in (0, 1):
            acc = {}
            for cell in cells[h::2]:
                a, b = int(ip[cell]), int(ip[cell + 1]); ii, cc = idx[a:b], cnt[a:b]
                scale32 = np.float32(1e4) / lib[cell]
                if panel is None and q >= 1:
                    sc = float(scale32); vals = [(int(i), float(c)) for i, c in zip(ii, cc)]
                else:
                    st = (seed ^ (((int(cell) + 1) * 0xD1B54A32D192ED03) % M)) % M; L = 0.0; raw = 0.0; vals = []
                    for i, c in zip(ii, cc):
                        raw += float(c); rr = (q if panel is None else (rate if int(i) in pidx else 0.0))
                        if rr <= 0: continue
                        if rr >= 1: t = int(c)
                        else:
                            t = 0
                            for _ in range(int(c)):
                                st, x = sm(st); t += float(x >> 11) < rr * 2.0**53
                        L += t
                        if t > 0: vals.append((int(i), float(t)))
                    if panel is None: L += q * max(1e4 / float(scale32) - raw, 0.0)
                    sc = 1e4 / L if L > 0 else 0.0
                for i, t in vals:
                    k = i if panel is None else pidx.get(i)
                    if k is not None: acc[k] = acc.get(k, 0.0) + t * sc
            for k, v in acc.items(): out[(n, h, k)] = v
    return out
small = gr[(gr.g != "control") & (gr.n.between(8, 14))].g.iloc[:3].tolist(); print("groups", small)
panel = genes[100:400:3]; open(T + "panel.txt", "w").write("\n".join(panel))
worst = 0
for q, pan, rate in [(1.0, None, None), (0.25, None, None), (0.05, None, None), (0.1, panel, 0.6)]:
    cmd = [B, "pbdump", "--cache", C, "--groups", ",".join(small), "--q", str(q), "--seed", "7", "--out", T + "pb.tsv"] + (["--panel", T + "panel.txt", "--rate", str(rate)] if pan else [])
    subprocess.run(cmd, check=True)
    R = pd.read_csv(T + "pb.tsv", sep="\t"); rd = {(a, b, c): d for a, b, c, d in R.itertuples(index=False)}
    P = pb(small, q, 7, pan, rate)
    keys = set(rd) | set(P); d = max(abs(rd.get(k, 0) - P.get(k, 0)) / max(1, abs(P.get(k, 0))) for k in keys)
    print(f"q {q} panel {pan is not None}: entries {len(keys)} max rel diff {d:.2e}"); worst = max(worst, d)
# truth z vs coexsc ess.z (same cache, same z definition)
m = pd.read_csv("../coexsc/out/ess.meta.tsv", sep="\t"); g = pd.read_csv("../coexsc/out/ess.genes.tsv", sep="\t")
Z = np.fromfile("../coexsc/out/ess.z.f32", np.float32).reshape(3, len(m), len(g))
tr = pd.read_csv("./out/simv.truth.tsv", sep="\t")
pi = {p: i for i, p in enumerate(m.pert)}; gi = {x: i for i, x in enumerate(g.gene)}
tr = tr[tr.pert.isin(pi) & tr.gene.isin(gi)]
e0 = np.abs(tr.z0.values - Z[1][tr.pert.map(pi), tr.gene.map(gi)]); e1 = np.abs(tr.z1.values - Z[2][tr.pert.map(pi), tr.gene.map(gi)])
print(f"truth z vs coexsc ess.z: {len(tr)} pairs, max abs diff z0 {e0.max():.2e} z1 {e1.max():.2e}")
