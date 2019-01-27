#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

library(ggplot2)

chr.lengths = read.delim('data/papAnu2_chr_lengths.txt',header=FALSE,col.names=c('chr','len'),row.names=NULL)
chr.lengths = chr.lengths[order(chr.lengths$chr),]
chr.lengths$chr = paste0('chr',chr.lengths$chr)

# No need for X chromosome
chr.lengths = chr.lengths[chr.lengths$chr %in% 1:20,]

next.lengths = c(0,chr.lengths$len[1:(length(chr.lengths$len)-1)])

# Create a vector of starting lengths). For each chromosome, the corresponding
# value in the vector + 1 will mark the unique index of position 1 on that chromosome

for (i in 2:length(next.lengths)) {
	next.lengths[i] = next.lengths[i] + next.lengths[i-1]
}

names(next.lengths) = chr.lengths$chr

chr.lengths$sta = next.lengths
chr.lengths$sto = chr.lengths$len + next.lengths

chr.lengths$color = c('#ffffff','#cccccc')

# Get all genes
genes = read.csv('results/baboon_genes_info.csv')

# No scaffolds, or MT/X
genes = genes[!genes$chr %in% c('MT','X'),]
genes$chr = factor(paste0('chr',genes$chr),levels=paste0('chr',1:20))

genes$mid = round((genes$sta.50 + genes$sto.50) / 2)

# This will give the position scaled across all chromosomes of the genome
genes$pos = genes$mid + next.lengths[genes$chr]

genes = genes[order(genes$pos),]

# make a data frame just for merging the chr and pos columns
gene.pos = genes[,c('nam','chr','pos')]

# ========================================================================================
# === Fst scan results
# ========================================================================================

fst.scan = read.csv('results/scan_selection_fst.csv')
fst.scan$sig = factor(fst.scan$pval <= 0.05,levels=c('TRUE','FALSE'))

fst.scan = merge(fst.scan,gene.pos,by='nam',all.x=TRUE,all.y=FALSE,sort=FALSE)

p = ggplot() +
	geom_rect(data=chr.lengths,aes(xmin=sta,xmax=sto,ymin=0,ymax=1,fill=color)) +
	scale_fill_manual(values=c('#ffffff','#cccccc')[(1:length(unique(fst.scan$chr))+1) %% 2 + 1],guide=FALSE) +
	geom_point(data=fst.scan,aes(pos,fst,color=sig),size=0.6) + theme_bw() +
	scale_color_manual(values=c('TRUE'='#ff0000','FALSE'='#000000')[(1:length(unique(fst.scan$chr))+1) %% 2 + 1],guide=FALSE) +
	theme_classic() +
	theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
	ylab(expression(italic(F[ST]))) + xlab('Position')
ggsave(p,file='results/gene_scans_fst.pdf',width=6,height=1,useDingbats=FALSE)

# ========================================================================================
# === BGC scan results
# ========================================================================================

bgc.scan = read.delim('results/bayesianetp_bgc.tsv')

# Sort bgc.scan by the gene name

bgc.scan = bgc.scan[order(bgc.scan$nam),]

bgc.pars = list()
for (i in 1:2) {
	parameter = bgc.scan[,ceiling(((1:7)-1)/3) %in% i]
	names(parameter) = c('bgc','lb','ub')
	parameter$par = c('alpha','beta')[i]
	bgc.pars[[c('alpha','beta')[i]]] = parameter
}
bgc.scan = data.frame(nam=bgc.scan$nam,do.call(rbind,bgc.pars),row.names=NULL)

bgc.scan$sig = factor(bgc.scan$lb > 0 | bgc.scan$ub < 0,levels=c('TRUE','FALSE'))
bgc.scan = merge(bgc.scan,gene.pos,by='nam',all.x=TRUE,all.y=FALSE,sort=FALSE)

p = ggplot() +
	geom_rect(data=chr.lengths,aes(xmin=sta,xmax=sto,ymin=min(bgc.scan$bgc),ymax=max(bgc.scan$bgc),fill=color)) +
	scale_fill_manual(values=c('#ffffff','#cccccc')[(1:length(unique(fst.scan$chr))+1) %% 2 + 1],guide=FALSE) +
	geom_point(data=bgc.scan,aes(pos,bgc,color=sig),size=0.6) + theme_bw() +
	scale_color_manual(values=c('TRUE'='#ff0000','FALSE'='#000000')[(1:length(unique(bgc.scan$chr))+1) %% 2 + 1],guide=FALSE) +
	theme_classic() +
	theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),strip.text.x = element_blank()) + 
	ylab('value') + xlab('Position') + facet_wrap(~par,ncol=1,scale='free_y')
ggsave(p,file='results/gene_scans_bgc.pdf',width=6,height=3,useDingbats=FALSE)

# Merge all measures together
bgc.a = bgc.scan[bgc.scan$par %in% 'alpha',c('nam','bgc','lb','ub','sig')]
bgc.b = bgc.scan[bgc.scan$par %in% 'beta',c('nam','bgc','lb','ub','sig')]

names(bgc.a)[2:ncol(bgc.a)] = paste0(names(bgc.a)[2:ncol(bgc.a)],'.a')
names(bgc.b)[2:ncol(bgc.b)] = paste0(names(bgc.b)[2:ncol(bgc.b)],'.b')

bgc = data.frame(bgc.a,bgc.b)

# Get rid of redundant names (proven to be redundant by the previous line)
bgc = bgc[,-grep('nam\\.[0-9]',names(bgc))]

gene.fst = fst.scan[,c('nam','fst','pval','nreps','sig')]
names(gene.fst)[2:ncol(gene.fst)] = paste0(names(gene.fst)[2:ncol(gene.fst)],'.fst')

full.scans = merge(bgc,gene.fst,by='nam',all.x=TRUE,all.y=TRUE,sort=FALSE)

# Quick test to see if beta will be more extreme for extremely diverged loci (two-sided)
beta.by.fst = wilcox.test(bgc.b~sig.fst,data=full.scans) # significant as expected
print(beta.by.fst)

# Quick test to see if alpha will be more extreme for extremely diverged loci (two-sided)
alpha.by.fst = wilcox.test(bgc.a~sig.fst,data=full.scans) # significant as expected
print(alpha.by.fst)

write.csv(full.scans,file='results/full_scans.csv',row.names=FALSE,quote=FALSE)


