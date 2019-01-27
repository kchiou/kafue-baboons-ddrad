#!/bin/bash

# -------------------------------------------------------------------------------------- #
# --- Calculate pi, nucleotide diversity, and H, observed heterozygosity, among others
# -------------------------------------------------------------------------------------- #

#	--freq              per-site frequency              *.frq
#	--counts            per-site counts                 *.frq.count
#	--depth             mean depth per individual       *.idepth            (fails)
#	--site-depth        depth for each site (summed)    *.ldepth            (fails)
#	--site-mean-depth   depth for each site (mean)      *.ldepth.mean       (fails)
#	--geno-depth        depth for each genotype         *.gdepth            (fails)
#	--site-quality      per-site SNP quality            *.lqual             (fails)
#	--het               inbreeding coefficient, F       *.het
#	--hardy             p-value for HWE                 *.hwe
#	--missing-indv      missingness per ind             *.imiss
#	--missing-site      missingness per site            *.lmiss
#	--geno-r2           squard corr coef btwn genotypes *.geno.ld           (takes awhile)
#	--SNPdensity <int>  SNP # and density for bin size  *.snpden
#	--TsTv <int>        Ts / Tv ratio for bin size      *.TsTv, out.TsTv.summary
#	--singletons        location of singletons          *.singletons
#	--site-pi           nuc. diversity per-site         *.sites.pi
#	--window-pi <int>   nuc. diversity in windows       *.windowed.pi
#	--TajimaD <int>     TajD in bins of size <int>      *.Tajima.D

# Calculate stats for Kindas
vcftools --vcf results/pure.kinda.vcf --out results/pure.kinda --freq
vcftools --vcf results/pure.kinda.vcf --out results/pure.kinda --counts
vcftools --vcf results/pure.kinda.vcf --out results/pure.kinda --het
vcftools --vcf results/pure.kinda.vcf --out results/pure.kinda --hardy
vcftools --vcf results/pure.kinda.vcf --out results/pure.kinda --missing-indv
vcftools --vcf results/pure.kinda.vcf --out results/pure.kinda --missing-site
vcftools --vcf results/pure.kinda.vcf --out results/pure.kinda --SNPdensity 1000000
vcftools --vcf results/pure.kinda.vcf --out results/pure.kinda --TsTv 1000000
vcftools --vcf results/pure.kinda.vcf --out results/pure.kinda --singletons
vcftools --vcf results/pure.kinda.vcf --out results/pure.kinda --site-pi
vcftools --vcf results/pure.kinda.vcf --out results/pure.kinda --TajimaD 300000

# Calculate stats for chacmas
vcftools --vcf results/pure.chacma.vcf --out results/pure.chacma --freq
vcftools --vcf results/pure.chacma.vcf --out results/pure.chacma --counts
vcftools --vcf results/pure.chacma.vcf --out results/pure.chacma --het
vcftools --vcf results/pure.chacma.vcf --out results/pure.chacma --hardy
vcftools --vcf results/pure.chacma.vcf --out results/pure.chacma --missing-indv
vcftools --vcf results/pure.chacma.vcf --out results/pure.chacma --missing-site
vcftools --vcf results/pure.chacma.vcf --out results/pure.chacma --SNPdensity 1000000
vcftools --vcf results/pure.chacma.vcf --out results/pure.chacma --TsTv 1000000
vcftools --vcf results/pure.chacma.vcf --out results/pure.chacma --singletons
vcftools --vcf results/pure.chacma.vcf --out results/pure.chacma --site-pi
vcftools --vcf results/pure.chacma.vcf --out results/pure.chacma --TajimaD 300000

# Calculate stats for hybrids
vcftools --vcf results/hybrids.vcf --out results/hybrids --freq
vcftools --vcf results/hybrids.vcf --out results/hybrids --counts
vcftools --vcf results/hybrids.vcf --out results/hybrids --het
vcftools --vcf results/hybrids.vcf --out results/hybrids --hardy
vcftools --vcf results/hybrids.vcf --out results/hybrids --missing-indv
vcftools --vcf results/hybrids.vcf --out results/hybrids --missing-site
vcftools --vcf results/hybrids.vcf --out results/hybrids --SNPdensity 1000000
vcftools --vcf results/hybrids.vcf --out results/hybrids --TsTv 1000000
vcftools --vcf results/hybrids.vcf --out results/hybrids --singletons
vcftools --vcf results/hybrids.vcf --out results/hybrids --site-pi
vcftools --vcf results/hybrids.vcf --out results/hybrids --TajimaD 300000

# Calculate Fst
# Need lists of animals in each population
# This generates results/all.pure.weir.fst
# Based on Weir and Cockerham 1984
vcftools --vcf results/all.pure.vcf --out results/all.pure --weir-fst-pop results/pure_kind_list.txt --weir-fst-pop results/pure_chac_list.txt

exit