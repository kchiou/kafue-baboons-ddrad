#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

vcf_path = 'results/all.cleaned.LDpruned.vcf'
ind_info = read.csv('data/individual_info.csv')

vcf_header = system(paste('grep -nr "#CHROM"',vcf_path,'| cut -d ":" -f 1 | sed -n 1p'),intern=TRUE)

# Pull out VCF genotype matrix header info
header_info = system(paste0('cat ',vcf_path,' | sed -n ',vcf_header,'p'),intern=TRUE)

header_info = strsplit(header_info,'\t')[[1]]

vcf.individuals = header_info[(which(header_info %in% 'FORMAT')+1):length(header_info)]

# This corrects a bug in which individual names are concatenated back to back earlier in the pipeline (may be due to underscores?)
vcf.individuals = gsub('\\.PE.*','',vcf.individuals)

# Subset ind_info to retain only individuals in the VCF file
ind.info = ind_info[ind_info$Individual.ID %in% vcf.individuals,]

vcf.coords = ind.info[,c('Longitude','Latitude')][match(vcf.individuals,ind.info$Individual.ID),]

write.table(vcf.coords,file='data/all.cleaned.LDpruned.coord',col.names=FALSE,row.names=FALSE,sep=' ')