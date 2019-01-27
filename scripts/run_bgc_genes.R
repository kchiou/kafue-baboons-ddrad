#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

library(ggplot2)
library(reshape2)
library(abind)
library(coda)
library(parallel)

# Bring in BGC SNP metadata
snps.bgc = read.csv('results/snps_bgc.csv')

snps.bgc = within(snps.bgc,{
	chr = as.numeric(gsub('chr([0-9]+):([0-9]+)','\\1',name))
	pos = as.numeric(gsub('chr([0-9]+):([0-9]+)','\\2',name))
	chr.sort = formatC(chr,digits=2,width=2,flag=0)
	pos.sort = formatC(pos,digits=9,width=9,flag=0)
	id.sort = paste('chr',chr.sort,':',pos.sort,sep='')
})

snps.ids = data.frame(id=sort(snps.bgc$id.sort))

bgc.match = snps.bgc[snps.bgc$id.sort %in% snps.ids$id,c('id.sort','alpha','beta','alpha.excess','beta.excess','alpha.outlier','beta.outlier','crazy.a','crazy.b')]
bgc.match = bgc.match[order(bgc.match$id.sort),]

snps = data.frame(snps.ids,bgc.match)

if (identical(bgc.match$id.sort,snps.ids$id)) bgc.match$id.sort = NULL else stop('SNP IDs do not match')

snps$chr = paste0('chr',as.numeric(gsub('chr([0-9]+?):.+','\\1',snps$id)))
snps$pos = as.numeric(gsub('chr[0-9]+?:(.+)','\\1',snps$id))

# Read in gene CDS, chromosome numbers, and gene positions
genes = read.csv('results/baboon_genes_info.csv')
full.genes = genes

# Ensure that there are no unplaced scaffolds, or MT/X
genes = genes[!genes$chr %in% c('MT','X'),]
genes$chr = paste0('chr',genes$chr)

library(IRanges)
library(GenomicRanges)

# Fix the line below
isnps = with(snps,GRanges(chr,IRanges(pos,width=1,names=id),'*'))
igenes = with(genes,GRanges(chr,IRanges(sta.50,sto.50,names=nam),'+'))

olaps = findOverlaps(isnps, igenes)

# Match SNPs with their nearest genes
nnbrs = nearest(isnps, igenes)

# Throw an error if chromosomes of matching SNPs/genes do not match
if (!all(snps$chr == genes$chr[nnbrs])) stop('Nearest neighbors found on non-matching chromosomes')

snps$nearest.gene = genes$nam[nnbrs]
snps$pos.1 = genes$sta[nnbrs]
snps$pos.2 = genes$sto[nnbrs]

# Calculate the distance to the nearest genes
snps$nearest.distance = do.call(c,lapply(1:nrow(snps),function(i) {
	pos = snps$pos[i]
	pos1 = snps$pos.1[i]
	pos2 = snps$pos.2[i]
	if (pos >= pos1 && pos <= pos2) {
		0
	} else {
		min(abs(pos - pos1),abs(pos - pos2))
	}
}))

genes$nam = factor(genes$nam,levels=genes$nam)

# Sanity check. We must ensure that the BGC SNP order is equivalent to the SNPs
bgc.snps = read.delim('results/bgc_snp_info.txt')

if (!identical(paste0(bgc.snps$chr,':',bgc.snps$pos),paste0(snps$chr,':',snps$pos)))
	stop('SNPs do not match!')

snp.ids = snps$id[queryHits(olaps)]
snp.indices = match(snp.ids,snps$id)
chr.ids = snps$chr[queryHits(olaps)]
pos.ids = snps$pos[queryHits(olaps)]

alpha.hits = snps$alpha[queryHits(olaps)]
beta.hits = snps$beta[queryHits(olaps)]
gene.hits = genes$nam[subjectHits(olaps)]

gene.snp.hits = data.frame(snp=snp.indices,chr=chr.ids,pos=pos.ids,gene=gene.hits)

# List of SNP ID and gene ID matches. The "snp" column gives the line number in the BGC posterior parameter files
write.csv(gene.snp.hits,'results/gene.snp.hits.csv',row.names=FALSE,quote=FALSE)

unique.genes = unique(gene.snp.hits$gene)

# For each unique gene, read in the posterior sample for the matching SNPs
# and recalculate a mean, median, and HPD interval

# Accomplish this using sed. sed cannot handle too many line number arguments, but this
# shouldn't be a problem as we will not approach that constraint. See vcf_to_bgc.R for
# a case where that constraint was overcome by calculating line ranges.

# PATHs to alpha and beta parameter chains
path.to.alpha = 'results/bgc_stat_a0'
path.to.beta = 'results/bgc_stat_b0'

cl = makeCluster(detectCores()-1)
clusterExport(cl,c('unique.genes','gene.snp.hits','path.to.alpha','path.to.beta'))

clusterEvalQ(cl, {
	library(coda)
})

new.mcmc = parLapply(cl,unique.genes,function(x) {

	# SNP indices
	i = gene.snp.hits$snp[gene.snp.hits$gene %in% x]
	cmd1 = paste0('sed -n -e ',paste(i,collapse='p -e '),'p ',path.to.alpha)
	cmd2 = paste0('sed -n -e ',paste(i,collapse='p -e '),'p ',path.to.beta)

	# Read in the sample estimates for alpha and beta from the posterior sample files
	a = do.call(rbind,lapply(strsplit(do.call(c,strsplit(system(cmd1,intern=TRUE),'\n')),','),as.numeric))
	b = do.call(rbind,lapply(strsplit(do.call(c,strsplit(system(cmd2,intern=TRUE),'\n')),','),as.numeric))

	# Next summarize the sample estimates for each gene such that each gene has 1 estimate per sample
	# This is the new MCMC sample for the GENE rather than the SNPs
	a = colMeans(a)
	b = colMeans(b)

	# Posterior point estimates can now be recalculated
	a.mean = mean(a)
	b.mean = mean(b)
	a.median = median(a)
	b.median = median(b)

	# Highest posterior density. Do not use this
	hpd.a = HPDinterval(mcmc(a))
	hpd.b = HPDinterval(mcmc(b))

	# Equal tail probability. This is used in Gompert & Buerkle (2012)
	etp.a = quantile(a, c(0.025, 0.975))
	etp.b = quantile(b, c(0.025, 0.975))

	p.a.gt = sum(a > 0) / length(a)
	p.a.lt = sum(a < 0) / length(a)
	p.b.gt = sum(b > 0) / length(b)
	p.b.lt = sum(b < 0) / length(b)

	a.lb = as.numeric(etp.a['2.5%'])
	a.ub = as.numeric(etp.a['97.5%'])
	b.lb = as.numeric(etp.b['2.5%'])
	b.ub = as.numeric(etp.b['97.5%'])

	data.frame(nam=x,bgc.a=a.mean,bgc.alb=a.lb,bgc.aub=a.ub,bgc.b=b.mean,bgc.blb=b.lb,bgc.bub=b.ub,a.pos=p.a.gt,a.neg=p.a.lt,b.pos=p.b.gt,b.neg=p.b.lt)
})

stopCluster(cl)

new.mcmc = do.call(rbind,new.mcmc)

write.table(new.mcmc,file='results/bayesianetp_bgc.tsv',sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)
