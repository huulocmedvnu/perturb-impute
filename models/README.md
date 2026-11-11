# Pre-trained reference models

| File | Reference | Panel | License |
|---|---|---|---|
| `k562_500_standard.bin` | K562 genome-wide atlas, 9,866 perturbations, 7,094 genes | greedy MMSE, 500 genes | MIT |
| `k562_500_lowbudget.bin` | same | UMI-cost greedy, 500 genes, noise weights for 2% of the reads | MIT |
| `xatlas/hct116_500_standard.bin` | X-Atlas/Orion HCT116 | greedy MMSE, 500 genes | CC BY-NC-SA 4.0 |
| `xatlas/hek293t_500_standard.bin` | X-Atlas/Orion HEK293T | greedy MMSE, 500 genes | CC BY-NC-SA 4.0 |
| `xatlas/pan_atlas_500_standard.bin` | K562, HCT116 and HEK293T pooled | greedy MMSE, 500 genes | CC BY-NC-SA 4.0 |

The models are 26 to 50 MB each. They are part of the Git repository but not of the journal
archive (Additional file 2), whose files are limited to 20 MB. Download them from
https://github.com/huulocmedvnu/perturb-impute/tree/main/models or rebuild them:

```sh
perturb-impute fit-reference --ref <k562 z-score prefix> --panel panels/standard_500.tsv \
  --groups models/inputs/k562_groups.tsv --out models/k562_500_standard.bin
perturb-impute fit-reference --ref <k562 z-score prefix> --panel panels/low-budget_500.tsv \
  --groups models/inputs/k562_groups.tsv --noise models/inputs/control_cv2_k562.tsv --budget 0.02 \
  --out models/k562_500_lowbudget.bin
```

`inputs/` holds what these commands need besides the z-scores: control-cell UMI fractions and squared
coefficients of variation (K562 essential-gene screen), Ensembl identifiers and perturbation groups
(CORUM complexes and correlated profiles). `scripts/pipeline/xatlas_models.sh` builds the X-Atlas/Orion models.
