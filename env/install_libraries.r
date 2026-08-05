# Set CRAN repos
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Install other CRAN packages
cran_packages <- c("foreach", "doParallel", "ggplot2", "cowplot", "tidyverse")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

# Install BiocManager
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(version = "3.22", ask = FALSE)

# Install Bioconductor packages
bioc_packages <- c("SummarizedExperiment", 
                   "SingleCellExperiment", "LoomExperiment",
                   "MungeSumstats", "rhdf5", "hdf5array",
                   "IRanges", "S4Vectors", "S4Arrays",
                   "DelayedArray", "XVector", "SNPlocs.Hsapiens.dbSNP155.GRCh38")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

# Install remotes for GitHub packages
if (!require("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

github_packages <- list(
    "scverse/anndataR" = "anndataR",
    "satijalab/seurat-object" = "SeuratObject",
    "satijalab/seurat" = "Seurat",
    "mojaveazure/seurat-disk" = "SeuratDisk",
    "NathanSkene/EWCE" = "EWCE",
    "neurogenomics/orthogene" = "orthogene",
    "neurogenomics/MAGMA_Celltyping" = "MAGMA.Celltyping",
    "bschilder/scKirby" = "scKirby"
)

for (repo in names(github_packages)) {
  pkg <- github_packages[[repo]]
  if (!require(pkg, character.only = TRUE)) {
    tryCatch({
      cat(sprintf("Installing %s...\n", pkg))
      remotes::install_github(repo, dependencies = TRUE, build_vignettes = FALSE)
    }, error = function(e) {
      cat(sprintf("Warning: Failed to install %s\n", pkg))
    })
  }
}

cat("\n=== Installation Complete ===\n")