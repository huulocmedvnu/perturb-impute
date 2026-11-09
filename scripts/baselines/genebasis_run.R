# geneBasis (Missarova et al. 2021, Genome Biol), official package geneBasisR (github MarioniLab/geneBasisR).
# "Cells" = reference perturbation profiles, features = 7,094 genes; the z matrix is stored in the
# 'logcounts' assay (the assay gene_search reads; z is not a log count, see README).
# Usage: Rscript genebasis_run.R FOLD N_GENES_TOTAL OUTFILE [NPC_SELECTION]   (NPC_SELECTION omitted = package default NULL)
# geneBasisR calls BiocNeighbors::queryKNN on a 1-column matrix that R drops to a vector when the selection
# holds one gene; BiocNeighbors 2.x (Bioc 3.20) rejects vectors, so the package's own example fails.
# We therefore load BiocNeighbors 1.22.0 + bluster 1.14.0 (Bioc 3.19, same R 4.4) from a private library.
.libPaths(c(file.path(Sys.getenv("HOME"), "Downloads/panelimp/ext/baselines/rlib_bn1"), .libPaths()))
suppressPackageStartupMessages({library(SingleCellExperiment); library(geneBasisR)})
cat("versions: geneBasisR", format(packageVersion("geneBasisR")), "BiocNeighbors", format(packageVersion("BiocNeighbors")),
    "bluster", format(packageVersion("bluster")), "SingleCellExperiment", format(packageVersion("SingleCellExperiment")), R.version.string, "\n")
a <- commandArgs(TRUE); f <- a[1]; n <- as.integer(a[2]); out <- a[3]
npc <- if (length(a) >= 4) as.integer(a[4]) else NULL
d <- file.path(Sys.getenv("HOME"), "Downloads/panelimp/ext/baselines/genebasis_in")
genes <- readLines(file.path(d, "genes.txt")); perts <- readLines(file.path(d, paste0("fold", f, ".perts.txt")))
con <- file(file.path(d, paste0("fold", f, ".f32")), "rb")
x <- readBin(con, "numeric", n = length(genes) * length(perts), size = 4); close(con)
m <- matrix(x, nrow = length(genes), ncol = length(perts)); rm(x)
rownames(m) <- genes; colnames(m) <- paste0("p", seq_along(perts), "_", perts)
sce <- SingleCellExperiment(assays = list(logcounts = m)); rm(m)
cat("fold", f, ": sce", nrow(sce), "genes x", ncol(sce), "profiles; n_genes_total =", n, "; nPC.selection =", if (is.null(npc)) "NULL" else npc, "\n")
t0 <- Sys.time()
res <- gene_search(sce, n_genes_total = n, nPC.selection = npc, verbose = TRUE)   # other arguments at package defaults
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("fold %s: gene_search done in %.1f s, %d genes\n", f, el, nrow(res)))
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(sprintf("%.1f", el), paste0(out, ".time"))
