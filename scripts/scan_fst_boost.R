#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

arguments = commandArgs(trailingOnly=TRUE)

# First argument gives the progress (10^n reps)
n = arguments[1]

# Second argument gives the random number
i = arguments[2]

nreps.completed = paste0('1e',n)

load(paste0('checkpoints/pval.perm.',nreps.completed,'.RData'))

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

fst.hits = fst.real$fst[queryHits(olaps)]
gene.hits = genes$nam[subjectHits(olaps)]

# This gives the Fst of real genes
genes$fst = as.numeric(tapply(fst.hits,gene.hits,mean))

# Create a representative distribution of genome-wide FST values
fst.genomewide = fst$fst[!is.nan(fst$fst)]

# Now create a fake table with all the Gene/Fst overlaps
fst.gene.overlaps = droplevels(data.frame(gene=gene.hits,fst=fst.hits))

nreps = as.numeric(nreps.completed)

# Figure out which genes have p-values that are not yet resolved
unfinished = real.fst[(real.fst$pval >= (0.05 - 500/nreps) & real.fst$pval <= (0.05 + 500/nreps)) | (real.fst$pval >= (0.01 - 500/nreps) & real.fst$pval <= (0.01 + 500/nreps)),]
unfinished$nam = as.character(unfinished$nam)

# Start with 10,000 simulations in which Fst is sampled with replacement from genome-wide values
# Genes are then reassigned a simulated Fst calculated as the average of all SNPs in the gene+100bp region

fst.gene.overlaps.subset = droplevels(fst.gene.overlaps[fst.gene.overlaps$gene %in% unfinished$nam,])
n.fst = nrow(fst.gene.overlaps.subset)

library(parallel)
cl = makeCluster(detectCores()-1)

# Export variables to cluster
clusterExport(cl,c('fst.genomewide','n.fst','fst.gene.overlaps.subset','nreps'))

fst.perm.results = parSapply(cl,1:nreps,function(i) {
	fst.sim = sample(fst.genomewide,n.fst,replace=FALSE)
	as.numeric(tapply(fst.sim,fst.gene.overlaps.subset$gene,mean))
})

# Stop cluster
stopCluster(cl)

out = data.frame(i=i,unfinished)

out$pval = rowMeans(matrix(fst.perm.results,nrow=nrow(unfinished)) >= unfinished$fst)

write.table(out,col.names=FALSE,row.names=FALSE,file=paste0('checkpoints/p.booster.',i,'.txt'))
