#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

fst = read.table('results/all.pure.weir.fst',header=TRUE)

# Read in gene CDS, chromosome numbers, and gene positions
genes = read.csv('results/baboon_genes_info.csv')
full.genes = genes

# No scaffolds, or MT/X
genes = genes[!genes$chr %in% c('MT','X'),]
genes$chr = paste0('chr',genes$chr)

names(fst) = c('chr','pos','fst')
fst$id = paste(fst$chr,fst$pos,sep=':')
fst = fst[,c('id','chr','pos','fst')]

# Get rid of SNPs that are not variable in these two pure samples
fst.real = fst[!is.nan(fst$fst),]

library(IRanges)
library(GenomicRanges)

isnps = with(fst.real,GRanges(chr,IRanges(pos,width=1,names=id),'*'))
igenes = with(genes,GRanges(chr,IRanges(sta.50,sto.50,names=nam),'+'))

olaps = findOverlaps(isnps, igenes)

genes$nam = factor(genes$nam,levels=genes$nam)

snp.hits = fst.real$id[queryHits(olaps)]
fst.hits = fst.real$fst[queryHits(olaps)]
gene.hits = genes$nam[subjectHits(olaps)]

# This gives the Fst of real genes
genes$fst = as.numeric(tapply(fst.hits,gene.hits,mean))

# Create a representative distribution of genome-wide FST values
fst.genomewide = fst$fst[!is.nan(fst$fst)]

gene.snp.overlaps = droplevels(data.frame(snp=snp.hits,gene=gene.hits))
write.csv(gene.snp.overlaps,file='results/gene_snp_overlaps.csv',row.names=FALSE,quote=FALSE)

# Now create a fake table with all the Gene/Fst overlaps
fst.gene.overlaps = droplevels(data.frame(gene=gene.hits,fst=fst.hits))

# Number of Fst values to simulate (total number of unique SNP-gene combinations)
n.fst = nrow(fst.gene.overlaps)

# Make cluster
library(parallel)
cl = makeCluster(detectCores()-1)

# Export variables to cluster
clusterExport(cl,c('fst.genomewide','n.fst','fst.gene.overlaps'))

# Start with 10,000 simulations in which Fst is sampled with replacement from genome-wide values
# Genes are then reassigned a simulated Fst calculated as the average of all SNPs in the gene+100bp region

fst.perm.results = parSapply(cl,1:10000,function(i) {
	fst.sim = sample(fst.genomewide,n.fst,replace=TRUE)
	as.numeric(tapply(fst.sim,fst.gene.overlaps$gene,mean))
})

# Stop cluster
stopCluster(cl)

# Unparallelized version
# fst.perm.results = replicate(10000,{
# 	fst.sim = sample(fst.genomewide,n.fst,replace=TRUE)
# 	as.numeric(tapply(fst.sim,fst.gene.overlaps$gene,mean))
# })

# Name the rows for easier merging
rownames(fst.perm.results) = names(table(fst.gene.overlaps$gene))

# Keep only non-NA genes and the name and Fst columns
real.fst = genes[!is.na(genes$fst),c('nam','fst')]
rownames(real.fst) = real.fst$nam

# Quick check that order of genes to be joined is identical
if (!identical(names(table(fst.gene.overlaps$gene)),as.character(real.fst$nam))) stop('Check that gene orders are compatible')

# P-value = proportion of permutation results that are as or more extreme than the real value
p.values = rowMeans(fst.perm.results >= real.fst$fst)

# Record p-values and number of simulations
real.fst$pval = p.values
real.fst$nreps = 10000

save(real.fst,file=paste0('checkpoints/fst.pval.perm.1e',log(10000,base=10),'.RData'))

nreps = 10000

# Explore further if p-value is close to alpha (explore further for alpha=0.05
while (sum(real.fst$pval >= (0.05 - 500/nreps) & real.fst$pval <= (0.05 + 500/nreps))) {

	this.round = as.character(real.fst$nam[which((real.fst$pval >= (0.05 - 500/nreps) & real.fst$pval <= (0.05 + 500/nreps))),drop=TRUE])
	this.round.overlaps = droplevels(fst.gene.overlaps[fst.gene.overlaps$gene %in% this.round,])
	n.fst = nrow(this.round.overlaps)
	
	cl = makeCluster(detectCores()-1)
	clusterExport(cl,c('fst.genomewide','n.fst','this.round.overlaps'))
	
	new.perm.results = parSapply(cl,1:(nreps*9),function(i) {
		fst.sim = sample(fst.genomewide,n.fst,replace=FALSE)
		as.numeric(tapply(fst.sim,this.round.overlaps$gene,mean))
	})

	stopCluster(cl)

	if (!identical(names(table(this.round.overlaps$gene)),as.character(this.round))) stop('Check that gene orders are compatible')

	# Proportion of permutation results that are more extreme than the real value
	p.values = rowMeans(new.perm.results >= real.fst[this.round,]$fst)

	# Calculate new weighted average of p-value (90% newest chain, 10% previous chain)
	real.fst$pval[match(this.round,rownames(real.fst))] = p.values * 0.9 + real.fst[this.round,]$pval * 0.1
	real.fst$nreps[match(this.round,rownames(real.fst))] = nreps * 10

	save(real.fst,file=paste0('checkpoints/fst.pval.perm.1e',log(nreps*10,base=10),'.RData'))
	
	nreps = nreps * 10

}

write.csv(real.fst,file='results/scan_selection_fst.csv',row.names=FALSE)
