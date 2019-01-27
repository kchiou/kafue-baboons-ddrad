#!/bin/bash

# Identify SNPs in LD
plink --noweb --bfile results/all.cleaned --indep-pairwise 50 5 0.5 --out results/all.cleaned
mv results/all.cleaned.log reports/plink_LD_pruning_part1.log

# Prune SNPS in LD
plink --noweb --bfile results/all.cleaned --exclude results/all.cleaned.prune.out --make-bed --out results/all.cleaned.LDpruned
mv results/all.cleaned.LDpruned.log reports/plink_LD_pruning_part2.log

exit