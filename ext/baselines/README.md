# Gene-panel baselines (highly variable genes and geneBasis)

Each baseline was run with its own official package. Outputs are tab-separated with the header
`fold	size	gene`. Each (fold, size) pair has exactly `size` unique genes from the 7,094-gene universe.

- `hvg.tsv`: fold `all`, sizes 200, 500 and 680
- `genebasis.tsv`: folds 0 to 4 and -9 (all reference perturbations), sizes 200, 500 and 680

## Shared inputs (`scripts/baselines/common.py`)

The z matrix is layer 0 (all cells) of the K562 genome-wide profiles, 9,866 perturbations by 7,094 genes.
The reference set for fold f contains the perturbations with fold other than f that are not excluded as
hubs for f, as in `src/main.rs`. Fold -9 uses all 9,866 perturbations.

## Highly variable genes (`scripts/baselines/hvg_scanpy.py`)

- Cells: the 10,691 non-targeting control cells of the K562 essential-gene screen, raw UMI counts.
- Genes: the universe genes found in the single-cell data (7,079 of 7,094).
- Call: `sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=k)` for k = 200, 500 and 680,
  other arguments at their defaults. Genes are written in `highly_variable_rank` order. The panels are nested.
- Environment: Python 3.12, scanpy 1.12.5, anndata 0.13.4, scikit-misc 0.5.3.

## geneBasis (`genebasis_export.py`, `genebasis_run.R`, `genebasis_all.sh`, `genebasis_rest.sh`, `genebasis_merge.py`)

- Package: geneBasisR from github MarioniLab/geneBasisR (commit 1e81fa70ab4457466706314e13fcd3223770495e),
  R 4.4.1, Bioconductor 3.20. Install with `scripts/baselines/install_genebasis.sh`.
- Input: the reference perturbation profiles as observations and the 7,094 genes as features, with z stored
  in the `logcounts` assay that `gene_search` reads. geneBasis uses this assay only for PCA, kNN graphs and
  Minkowski distances, so it does not assume positive values.
- Call: `gene_search(sce, n_genes_total = 680, nPC.selection = NULL)` with package defaults otherwise
  (n.neigh = 5, p.minkowski = 3, nPC.all = 50). The 200 and 500 panels are the first 200 and 500 genes by rank.
- Notes on the run:
  1. BiocNeighbors 1.22.0 and bluster 1.14.0 (Bioconductor 3.19) are loaded from a private library, because
     `gene_search` fails with BiocNeighbors 2.0.1 at its second step, also on the package example.
  2. Reference rows are subsampled to 3,000 per fold (`genebasis_export.py --sub 3000`, numpy default_rng(0)),
     because a run on all reference perturbations took more than 1.5 h per fold. All genes stay candidates.
- Runtime on 3,000 rows: 4,684 to 6,249 s per fold, single-threaded R.

## Checks

- Every (fold, size) has exactly `size` unique genes from the universe.
- The geneBasis 200 and 500 panels are rank prefixes of the 680 panel, and the HVG panels are nested.
- geneBasis 500-gene panels of different folds share 211 to 236 genes. geneBasis 500 and HVG 500 share 180 to 197 genes per fold.
