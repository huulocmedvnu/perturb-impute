# perturb-impute

Gene panels and transcriptome-wide imputation for CRISPR perturbation screens.

A targeted readout (hybridisation capture, targeted Perturb-seq, probe panels, MERFISH or Perturb-FISH) measures a few hundred genes per cell. `perturb-impute` chooses those genes on a reference perturbation atlas, fits an imputation model on the atlas, and returns the transcriptome-wide response of every perturbation in a new targeted screen, with a reliability score per gene (INFO), a local false sign rate per imputed change (lfsr) and an out-of-reference score per perturbation (OOR).

Written in Rust. Linear algebra runs through Apple Accelerate (macOS). No Python or R is needed.

## Install

```sh
cargo build --release          # binary: target/release/perturb-impute
```

Tested with Rust 1.99 on macOS (Accelerate). The build takes about 30 seconds.

## Quick start

Design a probe panel (wet lab), from the K562 genome-wide reference:

```sh
perturb-impute select --ref k562_gwps --panel-size 500 --mode standard \
  --umi-fraction control_umi_fraction_k562.tsv --gene-ids gene_ids.tsv --out panel.tsv
```

Ready-made panels are in `deliverables/Table_S2_Panels.tsv`, so this step is only needed for another reference or panel size.

Impute a targeted screen after sequencing, with a pre-trained model:

```sh
perturb-impute impute --model deliverables/models/k562_500_standard.bin \
  --input my_screen_panel_zscores.tsv --out-dir imputed/
```

Fit your own model on a reference atlas:

```sh
perturb-impute fit-reference --ref k562_gwps --panel panel.tsv --groups k562_groups.tsv --out model.bin
```

## Which panel

| Sequencing budget of the targeted assay (reads per cell relative to a whole-transcriptome experiment) | Panel | Model |
|---|---|---|
| 10% or more | `standard_500` (greedy minimum mean squared error) | `k562_500_standard.bin` |
| 2% to 10% | `low-budget_500` (greedy with a penalty on the UMI share of each gene, exponent 0.5) | `k562_500_lowbudget.bin` (noise weights set for 2% of the reads) |

The standard panel holds 52% of the UMIs of a K562 control cell and the low-budget panel 6%. At 2% of the reads, the low-budget panel recovered more of the specific response than the L1000 landmarks in a read-thinning simulation (manuscript, Fig. 5). For another budget, refit with `--noise control_cv2_k562.tsv --budget 0.05`.

## Input formats

- **Reference (`--ref`)**: either a z-score prefix (`PREFIX.meta.tsv` with perturbation names in the first column, `PREFIX.genes.tsv` with gene symbols in the first column, `PREFIX.z.f32` with float32 values of shape layers x perturbations x genes, layers = all cells, half 0, half 1), or a TSV/CSV matrix with perturbations in rows, genes in columns and names in the first column. With the three-layer format, INFO and calibration use half 0 to predict half 1 of held-out perturbations. With a single layer they are in-sample on the held-out perturbations.
- **Screen to impute (`--input`)**: TSV/CSV, one row per perturbation, one column per panel gene (gene symbols in the header). Values are z-scores of the pseudobulk difference from control cells, computed as in the reference (see Methods of the manuscript). Missing panel genes are allowed: the weights are solved on the genes present. `deliverables/e2e_test/example_panel_input.tsv` shows the format.
- **`--groups`**: perturbation and group (for example CORUM complex clusters). The 10% inner holdout is drawn by group, so that members of one complex are not split between training and holdout. Without groups, INFO is higher (median 0.194 against 0.171 for the K562 standard model) because related perturbations sit on both sides.

## Outputs of `impute`

| File | Content |
|---|---|
| `imputed_zscores.tsv` | perturbations x imputed genes, calibrated z-scores |
| `imputed_lfsr.tsv` | local false sign rate of each imputed change |
| `perturbation_qc.tsv` | out-of-reference score, panel norm, panel genes used |
| `gene_info.tsv` | INFO, calibration slope and error variance of each imputed gene |

## Filtering and calling DE genes

- **Genes**: keep imputed genes with INFO of 0.3 or more for gene-level analyses. In the K562 standard model, 804 of 6,594 imputed genes pass, and their median observed correlation in RPE1 cells was 0.46, against 0.38 for INFO of 0.2 or more (2,494 genes) and 0.28 below 0.3. Analyses of programs, clusters of perturbations and signatures can use all genes.
- **DE calls**: call an imputed gene DE when its lfsr is below 0.01, which bounds the probability that the sign of the change is wrong. In K562 cross-validation, the observed sign errors were 1.09 times the expected number. Outside the reference cell type, imputed changes were 13% to 18% smaller than measured changes, so effect sizes are conservative there.
- **Perturbations**: the OOR score rises when a whole screen differs from the reference (median 0.14 in K562, 0.23 to 0.25 in RPE1, HepG2 and Jurkat cells). It does not identify single failing perturbations outside the reference cell type.
- **Effect size ranking**: the ranking of DE genes by the size of their change is half as reliable as direct measurement. Add genes whose effect size matters to the panel.

## Deliverables (`deliverables/`)

- `Table_S2_Panels.tsv`: standard and low-budget panels of 200 and 500 genes from the K562 genome-wide atlas (9,866 perturbations, 7,094 genes), in selection order, with Ensembl gene identifiers (GENCODE v44, where available), control-cell UMI fractions and the variance explained under the reference model.
- `Table_S3_Gene_INFO_Scores.tsv`: for all 7,094 genes and both 500-gene panels, the K562 INFO score, calibration slope and error variance from the grouped inner holdout, and the observed per-gene correlation in RPE1, HepG2 and Jurkat cells (half 0 imputed, half 1 measured, perturbations whose halves correlate at 0.2 or more). For HepG2 and Jurkat cells, imputation used the panel genes measured there (450 of 500 standard and 392 of 500 low-budget genes).
- `models/k562_500_standard.bin`, `models/k562_500_lowbudget.bin`: pre-trained models (28.6 MB each).
- `inputs/`: control-cell UMI fractions and squared coefficients of variation (K562 essential-gene screen), gene identifiers and perturbation groups.

## Speed (MacBook Pro, M3 Pro, 18 GB)

| Step | Time |
|---|---|
| `select`, 500 genes from 9,866 perturbations x 7,094 genes | 4.8 s |
| `fit-reference` | 3.3 s |
| `impute`, 2,393 perturbations | 1.3 s (0.11 s computation) |

## Validation

`cargo test --release` runs the unit tests (model round trip, weights on a subset of panel genes, and checks of the analysis code against direct computations). `select` reproduces the greedy selection of the analysis pipeline gene for gene, and `validate` reproduces its per-gene correlations in RPE1 cells exactly.

## Reproducing the analyses

The analysis pipeline is the `panelimp` binary (`cargo build --release`, commands in `src/main.rs`) with helper scripts in `scripts/`:

- Pseudobulk z-scores: `dosecurve cache <screen.h5ad> <cache_dir>` (https://github.com/huulocmedvnu/dosecurve) writes a cell cache, `preprocess pb <cache_dir> <pb>` computes mean CP10k per perturbation over all cells and two interleaved halves, and `preprocess z --pb <pb> --genes <genes.txt> --out <z>` computes z-scores against non-targeting control cells (crate `preprocess/`).
- `panelimp folds` builds the complex-aware cross-validation folds from the z-scores and CORUM.
- `panelimp cv`, `panelimp de` and `panelimp sim` run panel selection with imputation, DE calls and the read-thinning simulation. `scripts/ms_numbers.py` collects the summary numbers.
- `orion-ingest` and `orion/scripts/*.sh` process the X-Atlas/Orion HCT116 and HEK293T atlases. `panelimp multi` and `panelimp vcm` run the pooled-reference and cross-cell-line transfer analyses.
- `scripts/baselines/` contains the highly variable gene and geneBasis baselines (see `ext/baselines/README.md`).
- `tests/` compares the Rust code with independent numpy implementations on small subsets.

Raw data: K562 and RPE1 screens from Figshare (https://doi.org/10.25452/figshare.plus.20029387), HepG2 and Jurkat screens from GEO (GSE264667), X-Atlas/Orion from Hugging Face (Xaira-Therapeutics/X-Atlas-Orion), CORUM complexes and the L1000 landmark list (GSE92742).

## License

Code: MIT (see `LICENSE`). The K562-derived panels and models in `deliverables/` are released under the same license. Panels and models derived from X-Atlas/Orion in `deliverables_xatlas/` follow the CC BY-NC-SA 4.0 license of these data (see `deliverables_xatlas/README.md`).

## Citation

Tran THL, Nguyen TV. Imputing transcriptome-wide CRISPR perturbation responses from small gene panels. Manuscript in preparation.
