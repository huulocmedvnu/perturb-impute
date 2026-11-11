# perturb-impute

Gene panels and transcriptome-wide imputation for CRISPR perturbation screens.

## Overview

A targeted readout (hybridisation capture, targeted Perturb-seq, probe panels, MERFISH or Perturb-FISH) measures a few hundred genes per cell. `perturb-impute` chooses those genes on a reference perturbation atlas by greedy minimum mean squared error (MMSE) selection, with a UMI-cost variant for low sequencing budgets, and imputes all other genes by ridge regression fitted on the atlas. Each imputed gene carries an INFO score, a calibration slope and a local false sign rate, and each perturbation an out-of-reference score. An optional guide-target prior sets the target gene of each knockdown, and genes with a transcription start site within 2 kb, to the values observed for reference knockdowns. For whole-transcriptome screens, `perturb-impute project` projects profiles on the leading programs of the reference atlas, which removes counting noise with the same reference information.

The repository holds the command-line tool, the analysis pipeline of the manuscript "Imputing transcriptome-wide CRISPR perturbation responses from small gene panels", the released panels and models, and the source data behind every figure.

## Quick start for reviewers

```sh
pip install -r environment/requirements.txt    # numpy, pandas, scipy, plotly, kaleido
make build      # compiles the Rust binaries (about 40 s on first build)
make test       # unit tests and the command-line workflow on toy data (under 10 s after the build)
make verify     # recomputes 162 numbers of the manuscript from source_data/ (a few seconds)
make figures    # regenerates Figs. 2 to 8 and Fig. S1 from source_data/ into output/ (about 10 s)
make all        # the four steps in sequence
```

No atlas has to be downloaded for these steps. Requirements are listed in `environment/README.md`: Rust 1.99.0 (pinned in `rust-toolchain.toml`), Python 3.12, and on Linux the OpenBLAS development package (`libopenblas-dev`) and, for the figures, `libexpat1`. On macOS the linear algebra runs through Apple Accelerate, which is part of the system.

`make verify` prints one line per number with the value in the manuscript, the recomputed value and its source. Most checks recompute medians, proportions and tests from the per-perturbation tables in `source_data/tables/`. Checks marked `summary` read summaries written by the Rust tool, and checks marked `json` read `source_data/numbers.json` for values that need the full atlas (number of programs, noise tails).

## Repository layout

| Path | Content |
|---|---|
| `src/bin/perturb_impute.rs` | command-line tool `perturb-impute` |
| `src/` | analysis pipeline `panelimp` (cross-validation, DE calls, simulation, paired libraries, pooled references, transferred profiles, and the alternative models of Table S1 in `causal.rs`, `explore.rs`, `explore_test.rs`) and `la.rs` (BLAS and LAPACK bindings) |
| `preprocess/`, `orion-ingest/` | pseudobulk z-scores from cell caches, and streaming ingest of the X-Atlas/Orion Parquet files |
| `panels/` | released panels and per-gene tables (Additional files 3 to 5 of the manuscript) |
| `models/` | pre-trained models and the inputs needed to rebuild them (`models/README.md`) |
| `source_data/` | per-perturbation tables (`tables/`), summaries (`real/`, `sim/`, `xatlas/`, `vcm/`, `alt/`, `example/`) and `numbers.json` behind the figures and tables |
| `scripts/` | `verify_manuscript_numbers.py`, `make_figures.py`, `make_flowchart.py`, `make_source_data.py`, `ms_numbers.py`, and the drivers of the full pipeline in `scripts/pipeline/` and `scripts/baselines/` |
| `tests/` | `cli.rs` (end-to-end test), `data/` (toy data, 0.7 MB), `smoke.sh`, and `validation/` (comparisons of the Rust code with independent numpy implementations, which need the full outputs) |
| `environment/` | pinned Python packages and notes on the toolchain |
| `ext/` | CORUM complexes, L1000 landmark list, gene universe, and the baseline panels (highly variable genes, geneBasis) |

## Using `perturb-impute`

Inputs are plain matrices: one row per perturbation, one column per gene, names in the first column, values as z-scores of the pseudobulk difference from control cells (Methods of the manuscript). A reference can also be given as a binary z-score prefix (`PREFIX.meta.tsv`, `PREFIX.genes.tsv`, `PREFIX.z.f32`).

1. Select a panel on a reference atlas. The standard panel is for 10% or more of the reads of a whole-transcriptome experiment, the low-budget panel for less.

```sh
perturb-impute select --ref reference.tsv --panel-size 500 --mode standard --out panel_standard.tsv
perturb-impute select --ref reference.tsv --panel-size 500 --mode low-budget \
  --umi-fraction models/inputs/control_umi_fraction_k562.tsv --out panel_low_budget.tsv
```

The panels of the manuscript are ready in `panels/` (`standard_500.tsv`, `low-budget_500.tsv`, 200-gene versions, and `Table_S2_Panels.tsv` with all four).

2. Fit a reference model. `--groups` keeps the members of a complex together when INFO scores and calibration are estimated.

```sh
perturb-impute fit-reference --ref reference.tsv --panel panel_standard.tsv \
  --groups models/inputs/k562_groups.tsv --out model.bin
# low budget: noise weights for the planned share of reads
perturb-impute fit-reference --ref reference.tsv --panel panel_low_budget.tsv --groups models/inputs/k562_groups.tsv \
  --noise models/inputs/control_cv2_k562.tsv --budget 0.02 --out model_low_budget.bin
```

3. Impute a targeted screen and keep the reliable genes.

```sh
perturb-impute impute --model models/k562_500_standard.bin --input my_screen_panel_zscores.tsv --out-dir imputed/
# genes with an INFO score of 0.3 or more
awk -F'\t' 'NR > 1 && $2 >= 0.3 {print $1}' imputed/gene_info.tsv > imputed/genes_info_0.3.txt
```

| Output | Content |
|---|---|
| `imputed_zscores.tsv` | perturbations by imputed genes, calibrated z-scores |
| `imputed_lfsr.tsv` | local false sign rate of each imputed change |
| `gene_info.tsv` | INFO score, calibration slope and error variance of each imputed gene |
| `perturbation_qc.tsv` | out-of-reference score, panel norm, number of panel genes used |

Missing panel genes are allowed: the weights are solved again on the genes present.

4. Project whole-transcriptome profiles on the reference programs.

```sh
perturb-impute project --ref reference.tsv --input my_screen_all_genes.tsv --programs 160 --out projected.tsv
```

5. Use the guide annotation (optional). The released K562 models hold the guide-target prior. With `--guide-target`, the imputed value of the target gene of each knockdown is replaced by the median z-score of that gene among reference knockdowns of the same expression decile, and the values of genes with a TSS within 2 kb by the median of such neighbours. Every other value is the ridge imputation. Row names of the input are read as target gene symbols, or `--target-map` gives them.

```sh
perturb-impute impute --model models/k562_500_standard.bin --input my_screen_panel_zscores.tsv \
  --guide-target --out-dir imputed/                 # also writes prior_applied.tsv
# a model with the prior, from your own reference
perturb-impute fit-reference --ref reference.tsv --panel panel_standard.tsv --groups models/inputs/k562_groups.tsv \
  --tss-annotation models/inputs/tss_grch38_refseq_select.tsv --expression models/inputs/control_umi_fraction_k562.tsv --out model.bin
```

The prior values are z-scores at the number of cells of a reference half (K562 models: a median of 89 cells per knockdown and half). A z-score grows with the square root of the number of cells, and `--prior-scale` multiplies the prior for screens with more or fewer cells. In K562 cross-validation the prior raised the mean rank correlation of effect sizes among DE genes outside the panel from 0.073 to 0.139 and the DE genes found from 0.206 to 0.229 (Additional file 1: Table S1).

Guidance from the manuscript: use a panel with imputation at 5% or less of the reads of a whole-transcriptome experiment or for readouts limited to a few hundred genes. With more reads, sequence the whole transcriptome and project. Use imputed genes with INFO of 0.3 or more for gene-level analyses (804 of 6,594 genes for the standard K562 panel), and all genes for analyses of programs. The out-of-reference score indicates when a screen differs from the reference (median 0.14 in K562 cells, 0.23 to 0.25 in other cell lines). It does not identify single failing perturbations outside the reference cell line. Effect sizes outside the reference cell line were 13% to 18% smaller than measured ones.

## Manuscript to code and source data

All paths are relative to the repository root. "pipeline" means the `panelimp` binary (`src/main.rs`), with the subcommand given.

| Item | Source data | Generated by | Manuscript section |
|---|---|---|---|
| Fig. 1 | `source_data/numbers.json`, `real/cells.summary.tsv`, `xatlas/*_cv.txt`, `vcm/vcmsum.txt` | `scripts/make_flowchart.py` | Study design |
| Fig. 2a, b | `numbers.json` (from `tables/cvg`, `trg`, `tr_hepg2`, `tr_jurkat` `.perts.tsv.gz`) | pipeline `cv --common --baselines`, `transfer`; `scripts/make_figures.py` | Panels chosen on perturbation responses |
| Fig. 2c | `tables/cvg.perts.tsv.gz` | pipeline `cv`; `make_figures.py` | same |
| Fig. 3a, c, d | `tables/cvm`, `trm_rpe1`, `trm_hepg2`, `trm_jurkat` `.perts.tsv.gz` | pipeline `cv`, `transfer`; `make_figures.py` | Imputation from 500 genes |
| Fig. 3b | `numbers.json` (from `tables/*.genes.tsv.gz`) | same | same |
| Fig. 3e | `example/mbtps2_fold1.tsv` | pipeline `example`; `make_figures.py` | same |
| Fig. 4a to d | `numbers.json` (from `tables/de_*.dex.tsv.gz`, `tables/*.de.tsv.gz`) | pipeline `upg`, `de`; `make_figures.py` | Measured and imputed genes combined |
| Fig. 4e | `numbers.json` (`rank_oracle`) | `scripts/pipeline/oracle_rank.py` | same |
| Fig. 5a, b | `numbers.json` (`sim`, from `tables/simd.sim.tsv.gz`), `sim/denoised.tsv` (from `tables/simd.den.tsv.gz`) | pipeline `sim --upg2`; `scripts/sim_denoised.py`; `make_figures.py` | Read-thinning simulation |
| Fig. 6a to d | `real/real.log`, `real/real.summary.tsv`, `real/cells.summary.tsv` (from `tables/paired.real.tsv.gz`, `paired.cells.tsv.gz`) | pipeline `mtx2cache`, `real`, `realcells`, `realsum`, `cellsum` (`scripts/pipeline/real.sh`); `make_figures.py` | Enriched library of the L1000 landmarks |
| Fig. 7a, b | `xatlas/hct116_*.txt`, `xatlas/hek293t_*.txt`, `numbers.json` | `orion-ingest`, pipeline `cv`, `upg`, `cvsum` (`scripts/pipeline/xatlas_analysis.sh`, `xatlas_l1000.sh`); `make_figures.py` | Two independent atlases |
| Fig. 7c | `xatlas/multi_panframe.txt` | pipeline `multi`, `multisum` (`scripts/pipeline/xatlas_multi_tf.sh`); `make_figures.py` | Pooled reference |
| Fig. 8a to c | `vcm/vcmsum.txt` | pipeline `vcm`, `vcmsum`; `make_figures.py` | Transferred profiles |
| Fig. S1 | `numbers.json` (`sim_budget`, from `tables/simb.sim.tsv.gz`) | pipeline `sim --design-budgets`; `make_figures.py` | Read-thinning simulation |
| Table 1 | none (recommendations, with the evidence cited in the table) | | Discussion |
| Table S1 | `numbers.json` (`rank_oracle`, `rank_knn`), `tables/cvk.knn.tsv.gz` | pipeline `cv` (neighbour copying), `upg`; `scripts/pipeline/oracle_rank.py` | Effect size ranking |
| Table S1 (continued), Note S1 | `alt/benchmark.tsv`, `alt/paired_against_ridge.tsv`, `alt/gene_graphs.tsv`, `alt/diagnostic_fold0.txt` | pipeline `explore` (diagnostic), `explore --stage test` (per fold, and per screen with `--target`), `causal`; `scripts/make_alt_source.py` | Effect size ranking |
| Table S2 (Additional file 3) | `panels/Table_S2_Panels.tsv` | `perturb-impute select` | Data availability |
| Table S3 (Additional file 4) | `panels/Table_S3_Gene_INFO_Scores.tsv` | `perturb-impute fit-reference`, `validate` | Discussion |
| Table S4 | `tables/de_k562_all`, `de_rpe1`, `de_hepg2`, `de_jurkat` `.dex.tsv.gz` | pipeline `upg` | DE calls |
| Table S5 | `xatlas/*_high.txt`, `xatlas/multi_tf.txt`, `xatlas/multi_tf_test.txt` | pipeline `cv`, `multi` | Independent atlases |
| Table S6 | not in `source_data/` (per-stratum counts of the DE runs) | pipeline `de --panel-size N [--same-cells]` | DE calls |
| Table S7 | `vcm/vcmsum.txt` | pipeline `vcm`, `vcmsum` | Transferred profiles |
| Table S8 | `real/real.summary.tsv`, `real/cells.summary.tsv` | pipeline `realsum`, `cellsum` | Enriched library |
| Additional file 5 | `panels/xatlas/*.tsv` | `perturb-impute select` (`scripts/pipeline/xatlas_models.sh`) | Data availability |

`scripts/make_source_data.py` writes `source_data/tables/` from the pipeline outputs, and `scripts/ms_numbers.py` writes `numbers.json`. File names in `tables/` are those of the pipeline outputs. Columns: `r_ceil` and `r_imp` are total correlations of the measured input half and of the imputed profile with the held-out half, `r_ceil_perp` and `r_imp_perp` the specific correlations, `r_den160_perp` the specific correlation of the profile projected on 160 programs, `oor` the out-of-reference score, `n_de`, `tp` and `disc` the numbers of DE genes, of calls with the correct sign and of calls with the opposite sign.

## Full pipeline from the public atlases

The steps below repeat the analyses from raw data. They need several tens of GB of free disk space (the working files of X-Atlas/Orion alone take 25 GB) and 18 GB of memory. On a MacBook Pro (M3 Pro, 18 GB), the simulation takes about 6 minutes, the paired libraries 4 minutes, and the ingest of the two X-Atlas/Orion atlases about 70 minutes plus the download.

| Data | Source |
|---|---|
| K562 genome-wide, K562 essential and RPE1 screens (Replogle et al. 2022) | Figshare, https://doi.org/10.25452/figshare.plus.20029387 |
| HepG2 and Jurkat screens (Nadig et al. 2025) | GEO GSE264667 |
| Paired whole-transcriptome and L1000-enriched libraries (Replogle et al. 2020, samples exp7 and exp11) | GEO GSE146194 |
| X-Atlas/Orion HCT116 and HEK293T | Hugging Face, `Xaira-Therapeutics/X-Atlas-Orion` (CC BY-NC-SA 4.0) |
| CORUM complexes | https://mips.helmholtz-muenchen.de/corum (copy in `ext/corum_human.txt`) |
| L1000 landmark genes | GEO GSE92742 gene annotation (list in `ext/l1000_landmarks.txt`) |

1. Pseudobulk z-scores. `dosecurve cache <screen.h5ad> <cache_dir>` (https://github.com/huulocmedvnu/dosecurve) writes a sparse cell cache. `preprocess pb <cache_dir> <pb>` computes mean counts per 10,000 UMIs per perturbation over all cells and two alternating halves, and `preprocess z --pb <pb> --genes ext/universe.txt --out out/gwps` the z-scores against control cells.
2. Folds. `panelimp folds` builds the cross-validation folds from the z-scores and CORUM (`out/folds.tsv`, `out/hub_excl.tsv`).
3. Selectors, imputation and DE. `panelimp cv`, `transfer`, `de` and `upg`, with the baseline panels of `scripts/baselines/` (`ext/baselines/README.md`).
4. Simulation. `panelimp sim --upg2` on the cell cache of the K562 essential-gene screen, then `scripts/sim_denoised.py`.
5. Paired libraries. `scripts/pipeline/real.sh` after the two matrices of GSE146194 are placed in `ext/replogle2020/`.
6. X-Atlas/Orion. `orion-ingest` and `scripts/pipeline/xatlas_*.sh`, then `panelimp multi` and `panelimp vcm`.
7. Alternative models and the guide-target prior. `panelimp explore`, `panelimp explore --stage test` for each fold and with `--target` for the other screens, and `panelimp causal` (optional Metal backend with `--features metal` on macOS), then `scripts/make_alt_source.py`.
8. Numbers and figures. `scripts/ms_numbers.py`, `scripts/make_source_data.py`, then `make verify figures`.

The usage of every pipeline subcommand is documented at the top of `src/main.rs`.

## License and citation

Code, the K562 models and the panels chosen on the K562 atlas: MIT (`LICENSE`). Panels and models derived from X-Atlas/Orion (`panels/xatlas/`, `models/xatlas/`): CC BY-NC-SA 4.0, the license of these data (`models/xatlas/README.md`).

Tran THL, Nguyen TV. Imputing transcriptome-wide CRISPR perturbation responses from small gene panels. Manuscript submitted. See `CITATION.cff`.
