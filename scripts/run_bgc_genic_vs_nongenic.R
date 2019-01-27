#!/usr/bin/env Rscript

snps.bgc = read.csv('results/snps_bgc.csv')

gene.snp.hits = read.csv('results/gene.snp.hits.csv')

gene.snp.hits$id = paste(gene.snp.hits$chr,gene.snp.hits$pos,sep=':')

# Check if snp is in a gene

snps.bgc$gene = snps.bgc$name %in% gene.snp.hits$id
snps.bgc$gene = factor(snps.bgc$gene,levels=c('TRUE','FALSE'))

print(wilcox.test(beta~gene,snps.bgc,alternative='greater'))