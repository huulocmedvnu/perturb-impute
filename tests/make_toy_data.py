"""Synthetic toy data for the command-line tests (tests/cli.rs): a reference atlas of 240 perturbations
and 300 genes whose responses follow 8 programs, and a screen of 40 new perturbations plus one profile
outside the reference programs. Regenerate with: python3 tests/make_toy_data.py"""
import os, numpy as np
rng = np.random.default_rng(2026); D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
G, K, NREF, NNEW = 300, 8, 240, 40
genes = [f"G{j:03d}" for j in range(G)]
L = rng.normal(size=(G, K)) * rng.gamma(2.0, 0.5, size=(G, 1))            # gene loadings on programs
def profiles(n): return (rng.standard_t(4, size=(n, K)) * 1.5) @ L.T      # noise-free responses
def write(name, rows, M):
    with open(os.path.join(D, name), "w") as f:
        f.write("perturbation\t" + "\t".join(genes) + "\n")
        for r, v in zip(rows, M): f.write(r + "\t" + "\t".join(f"{x:.3f}" for x in v) + "\n")
ref = profiles(NREF); write("toy_reference.tsv", [f"REF{i:03d}" for i in range(NREF)], ref + rng.normal(size=ref.shape))
new = np.vstack([profiles(NNEW), rng.normal(size=(1, G)) * 6.0])            # last row: no program structure
names = [f"NEW{i:02d}" for i in range(NNEW)] + ["OUT_OF_REFERENCE"]
write("toy_screen_input.tsv", names, new + rng.normal(size=new.shape))      # measured half used as input
write("toy_screen_truth.tsv", names, new + rng.normal(size=new.shape))      # independent half used as truth
frac = rng.dirichlet(np.full(G, 0.3)); cv2_bio = rng.gamma(2.0, 0.05, G); cv2_tech = 1e-3 / frac
with open(os.path.join(D, "toy_umi_fraction.tsv"), "w") as f: f.writelines(f"{g}\t{u:.6e}\n" for g, u in zip(genes, frac))
with open(os.path.join(D, "toy_noise.tsv"), "w") as f:
    f.write("gene\tcv2_bio\tcv2_tech\tcv2_tot\tshare\n"); f.writelines(f"{g}\t{b:.6e}\t{t:.6e}\t{b + t:.6e}\t{u:.6e}\n" for g, b, t, u in zip(genes, cv2_bio, cv2_tech, frac))
with open(os.path.join(D, "toy_groups.tsv"), "w") as f:
    f.write("pert\tgroup\n"); f.writelines(f"REF{i:03d}\t{i // 3}\n" for i in range(NREF))   # complexes of three perturbations
