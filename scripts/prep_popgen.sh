#!/bin/bash

# ========================================================================================
# === Subset data
# ========================================================================================

# Determine list of alleles deemed to be the "reference" in plink
# List of reference alleles depends on pruned binary PED file
plink --noweb --bfile results/all.cleaned --recode lgen-ref --out results/all.cleaned
mv results/all.cleaned.log reports/plink_ref_alleles.log
cut -d' ' -f1-2 results/all.cleaned.ref > results/all.cleaned.just.ref

# Get list of pure animals
scripts/find_pure_animals.R

# Individual IDs in the bedfile annoyingly have .PE tacked on to the end
cat results/pure_kind_list.txt | sed -re 's/([A-Za-z0-9_]+)/\1.PE/g' > results/pure_kind_list.tmp
cat results/pure_chac_list.txt | sed -re 's/([A-Za-z0-9_]+)/\1.PE/g' > results/pure_chac_list.tmp
cat results/hybrid_list.txt | sed -re 's/([A-Za-z0-9_]+)/\1.PE/g' > results/hybrid_list.tmp

# Rename .tmp files
mv results/pure_kind_list.tmp results/pure_kind_list.txt
mv results/pure_chac_list.tmp results/pure_chac_list.txt
mv results/hybrid_list.tmp results/hybrid_list.txt

# Subset pure Kindas (forcing ref allele to be same)
plink --noweb --bfile results/all.cleaned --keep results/pure_kind_list.txt --reference-allele results/all.cleaned.just.ref --make-bed --out results/pure.kinda
mv results/pure.kinda.log reports/plink_subset_kind.log

# Subset pure chacmas (forcing ref allele to be same)
plink --noweb --bfile results/all.cleaned --keep results/pure_chac_list.txt --reference-allele results/all.cleaned.just.ref --make-bed --out results/pure.chacma
mv results/pure.chacma.log reports/plink_subset_kind.log

# Subset hybrids (forcing ref allele to be same)
plink --noweb --bfile results/all.cleaned --keep results/hybrid_list.txt --reference-allele results/all.cleaned.just.ref --make-bed --out results/hybrids
mv results/hybrids.log reports/plink_subset_hybr.log

# ========================================================================================
# === Convert data
# ========================================================================================

# Make pure Kinda VCF
pseq pure.kinda.pseq new-project
pseq pure.kinda.pseq load-plink --file results/pure.kinda --id pure.kinda
pseq pure.kinda.pseq write-vcf > results/pure.kinda.vcf
rm -r pure.kinda.pseq*

# Make pure chacma VCF
pseq pure.chacma.pseq new-project
pseq pure.chacma.pseq load-plink --file results/pure.chacma --id pure.chacma
pseq pure.chacma.pseq write-vcf > results/pure.chacma.vcf
rm -r pure.chacma.pseq*

# Make hybrids VCF
pseq hybrids.pseq new-project
pseq hybrids.pseq load-plink --file results/hybrids --id hybrids
pseq hybrids.pseq write-vcf > results/hybrids.vcf
rm -r hybrids.pseq*

# Sort and compress Kinda VCF
vcf-sort results/pure.kinda.vcf > pure.kinda.sort
bgzip -c pure.kinda.sort > results/pure.kinda.vcf.gz
rm pure.kinda.sort
tabix -p vcf results/pure.kinda.vcf.gz

# Sort and compress chacma VCF
vcf-sort results/pure.chacma.vcf > pure.chacma.sort
bgzip -c pure.chacma.sort > results/pure.chacma.vcf.gz
rm pure.chacma.sort
tabix -p vcf results/pure.chacma.vcf.gz

# Sort and compress hybrids VCF
vcf-sort results/hybrids.vcf > hybrids.sort
bgzip -c hybrids.sort > results/hybrids.vcf.gz
rm hybrids.sort
tabix -p vcf results/hybrids.vcf.gz

# Make VCF file containing all pure animals
vcf-merge results/pure.kinda.vcf.gz results/pure.chacma.vcf.gz > results/all.pure.vcf

exit