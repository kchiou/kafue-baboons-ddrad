#!/bin/bash

# Remove repetitive regions (autosomal markers only)
plink --noweb --file data/kafue.pass.snp.noX --exclude data/excluded_repetitive_regions.plink.txt --range  --make-bed --out results/all.norep --missing
mv results/all.norep.log reports/plink_repetitive_region_exclusion.log

# ========================================================================================
# === Make all.cleaned dataset
# ========================================================================================

# Remove SNPs with too much missing data
plink --noweb --bfile results/all.norep --geno 0.2 --make-bed --out results/all.cleaned --missing
mv results/all.cleaned.log reports/plink_missing_exclusion_main1.log

plink --noweb --bfile results/all.cleaned --mind 0.8 --make-bed --out results/all.cleaned --missing
mv results/all.cleaned.log reports/plink_missing_exclusion_main2.log

# Load file to get final genotyping rate
plink --noweb --bfile results/all.cleaned --make-bed --out results/all.cleaned
mv results/all.cleaned.log reports/plink_info_full.log

# ========================================================================================
# === Make strict.cleaned dataset
# ========================================================================================

# Filter stringently to create smaller dataset
plink --noweb --bfile results/all.cleaned --geno 0.1 --recode --out results/strict.cleaned.tmp
mv results/strict.cleaned.tmp.log reports/plink_missing_exclusion_strict_part1.log

plink --noweb --file results/strict.cleaned.tmp --maf 0.05 --recode --out results/strict.cleaned --missing
mv results/strict.cleaned.log reports/plink_missing_exclusion_strict_part2.log

# Convert smaller dataset SNP file to binary format
plink --noweb --file results/strict.cleaned --make-bed --out results/strict.cleaned
mv results/strict.cleaned.log reports/plink_info_strict.log

exit