#!/bin/bash

# ========================================================================================
# === PCA
# ========================================================================================

# Perform principal components analysis
plink --noweb --bfile results/all.cleaned.LDpruned --pca --out results/all.cleaned.LDpruned

# Make plot of PCA
scripts/plot_pca.R

exit