#!/bin/bash

# ========================================================================================
# === MDS
# ========================================================================================

# Make genome file based on pruned file
plink --noweb --bfile results/all.cleaned.LDpruned --genome --out results/all.cleaned.LDpruned
mv results/all.cleaned.LDpruned.log reports/plink_genome_file.log

# Perform multidimensional scaling of similarities (proportions of alleles IBS)
plink --noweb --bfile results/all.cleaned.LDpruned --read-genome results/all.cleaned.LDpruned.genome --cluster --mds-plot 2 --out results/all.cleaned.LDpruned
mv results/all.cleaned.LDpruned.log reports/plink_mds.log

# Generate distance matrix (1 - IBS)
plink --noweb --bfile results/all.cleaned.LDpruned --distance 1-ibs flat-missing --out results/all.cleaned.LDpruned
mv results/all.cleaned.LDpruned.log reports/plink_dist_matrix.log

# The below line is for backwards compatibility with PLINK 1.07
# plink --noweb --bfile results/all.cleaned.LDpruned --read-genome results/all.cleaned.LDpruned.genome --cluster --distance-matrix --out results/all.cleaned.LDpruned

# Make MDS plot of IBS (depends on MDS file and individual info input file)
scripts/plot_mds.R

exit