#!/bin/bash
# Installs geneBasisR (MarioniLab GitHub; the method is not on Bioconductor 3.20) into the user R library.
# Toolchain workarounds for this Mac (macOS 27 SDK in CommandLineTools is unreadable by ld-1267;
# Homebrew gfortran 14 maps darwin27 to an invalid -mmacosx-version-min=18.0):
export SDKROOT=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk
export MACOSX_DEPLOYMENT_TARGET=26.0
/opt/homebrew/bin/Rscript -e '
options(repos = BiocManager::repositories())
if (!requireNamespace("Deriv", quietly = TRUE)) tryCatch(install.packages("Deriv"), error = function(e) NULL)
if (!requireNamespace("Deriv", quietly = TRUE)) remotes::install_version("Deriv", version = "4.1.6", upgrade = "never")
remotes::install_github("MarioniLab/geneBasisR", ref = "1e81fa70ab4457466706314e13fcd3223770495e", upgrade = "never", dependencies = TRUE, Ncpus = 3)
library(geneBasisR); cat("geneBasisR", as.character(packageVersion("geneBasisR")), "\n")'
# geneBasisR fails with BiocNeighbors 2.x (1-gene selection -> vector -> queryKNN signature error, also in the
# package's own example). Private library with the Bioc 3.19 versions, used only by genebasis_run.R:
mkdir -p ~/Downloads/panelimp/ext/baselines/rlib_bn1
/opt/homebrew/bin/Rscript -e '
lib <- path.expand("~/Downloads/panelimp/ext/baselines/rlib_bn1")
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("RcppHNSW", quietly = TRUE)) install.packages("RcppHNSW")
install.packages("https://bioconductor.org/packages/3.19/bioc/src/contrib/BiocNeighbors_1.22.0.tar.gz", repos = NULL, type = "source", lib = lib)
install.packages("https://bioconductor.org/packages/3.19/bioc/src/contrib/bluster_1.14.0.tar.gz", repos = NULL, type = "source", lib = lib)'
