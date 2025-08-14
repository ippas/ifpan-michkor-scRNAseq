# Install required packages (if not already installed)
# install.packages("remotes")  # For installing packages from GitHub
# remotes::install_github("satijalab/seurat", ref = "release/5.0.0")  # Install Seurat v5
# install.packages("openxlsx", repos = "https://cloud.r-project.org")  # Excel writing support
# 
# BiocManager::install("SingleR")
# BiocManager::install("celldex") 
# BiocManager::install("TENxPBMCData")
# 
# install.packages("patchwork")
# install.packages("enrichR")


# CRAN
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("enrichR", quietly = TRUE)) install.packages("enrichR")

# GitHub
if (!requireNamespace("Seurat", quietly = TRUE)) remotes::install_github("satijalab/seurat", ref = "release/5.0.0")

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("SingleR", quietly = TRUE)) BiocManager::install("SingleR")
if (!requireNamespace("celldex", quietly = TRUE)) BiocManager::install("celldex")
if (!requireNamespace("TENxPBMCData", quietly = TRUE)) BiocManager::install("TENxPBMCData")


# Load required libraries
library(Seurat)        # Main package for single-cell RNA-seq analysis
library(Matrix)        # Sparse matrix support
library(data.table)    # Fast reading/writing of data (especially CSV)
library(tibble)        # Enhanced data frames (tibbles)
library(dplyr)         # Data manipulation
library(tidyr)         # Data reshaping (e.g. separate, pivot)
library(purrr)         # Functional programming tools (e.g. map)
library(openxlsx)      # Writing Excel files (e.g. with filters and styling)
library(SingleR)
library(celldex)
library(patchwork)
library(ggplot2)
library(stringr)
library(enrichR)
library(magrittr)
library(scales)
library(SingleR) 
library(celldex) 
library(SummarizedExperiment) 
library(plyr)
library(SingleCellExperiment)
library(scuttle)
library(zellkonverter)
library(HDF5Array)
library(glue)


# ##############################################################################
# ---- load functions ----
# ##############################################################################
lapply(list.files("preprocessing/src/seurat_utils", "\\.R$", full.names = TRUE), source)

lapply(list.files("preprocessing/src/visualization_utils/umap/", "\\.R$", full.names = TRUE), source)

lapply(list.files("preprocessing/src/io/", "\\.R$", full.names = TRUE), source)

lapply(list.files("preprocessing/src/enrichr_utils/", "\\.R$", full.names = TRUE), source)

lapply(list.files("preprocessing/src/singler_annotation_utils/", "\\.R$", full.names = TRUE), source)

lapply(list.files("preprocessing/src/visualization_utils/barplot/", "\\.R$", full.names = TRUE), source)


# ##############################################################################
# ---- load Rdata ----
# ##############################################################################
# load("data/pilotDexCortCtrl/r_outputs/merged_pilotDexCortCtrl.RData")

load("data/pilotDexCortCtrl/r_outputs/integrated_pilotDexCortCtrl.RData")

# load("data/pilotDexCortCtrl/r_outputs/integrated_pilotDexCortCtrl_backupFirstIntegration.RData")

