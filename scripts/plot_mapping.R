#!/usr/bin/env Rscript

library(reshape2)
library(ggplot2)

options(stringsAsFactors=FALSE)

# Read in individual and mapping info
ind.info = read.csv('data/individual_info.csv')

# Read in a dataset that includes only individuals passing filters
mds = read.table('results/mds.txt')
ind.info$Passed.filter = ind.info$Individual.ID %in% mds$FID

# Subset feces
feces = ind.info[ind.info$Sample.type %in% 'feces',]

# Create inverse datasets of host and nonhost (pre- and post-enrichment)
feces.mapped = data.frame(id=feces$ID,Host=feces$Perc.mapped,Nonhost=1-feces$Perc.mapped,Set='Fraction reads mapping to genome')
feces.sample = data.frame(id=feces$ID,Host=feces$Sample.Endogenousity,Nonhost=1-feces$Sample.Endogenousity,Set='Fraction host DNA in fecal extraction')

# Combine datasets
feces.host = rbind(feces.mapped,feces.sample)

# Order by decreasing sequencing efficiency (i.e., reads mapping percentage)
feces.host$id = factor(feces.host$id,levels=feces$ID[order(feces$Perc.mapped,decreasing=TRUE)])

ggplot(melt(feces.host),aes(id,value,fill=variable)) +
	geom_bar(stat='identity',position='stack') +
	scale_fill_manual(values=c('#31a354','#e5f5e0'),guide=FALSE) +
	scale_x_discrete(name=NULL) +
	scale_y_continuous(name=NULL) +
	facet_wrap(~Set,ncol=1) +
	theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5))
ggsave('results/reads_mapped_rate.pdf',width=7,height=7)

all.mapped = data.frame(id=ind.info$ID,Type=ind.info$Sample.type,Reads=ind.info$Reads.mapped,Passed=ind.info$Passed.filter)
all.mapped$Type[!all.mapped$Type %in% 'feces'] = 'blood'
all.mapped$id = factor(all.mapped$id,levels=all.mapped$id[order(all.mapped$Type,all.mapped$Reads,decreasing=TRUE)])

ggplot(all.mapped,aes(id,Reads,fill=Type,alpha=Passed)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Type',values=c('red','darkgoldenrod'),guide=FALSE) +
	scale_x_discrete(name=NULL) +
	scale_alpha_manual(values=c(0.5,1)) +
	coord_cartesian(ylim=c(0,2e6)) +
	facet_wrap(~Type,ncol=1,scales='free_x') +
	theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5,size=5),legend.position='none') +
	ggtitle('Reads mapped to genome')
ggsave('results/reads_mapped_count_clipped.pdf',width=7,height=7)

ggplot(all.mapped,aes(id,Reads,fill=Type,alpha=Passed)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Type',values=c('red','darkgoldenrod'),guide=FALSE) +
	scale_x_discrete(name=NULL) +
	scale_alpha_manual(values=c(0.5,1)) +
	facet_wrap(~Type,ncol=1,scales='free_x') +
	theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5,size=5),legend.position='none') +
	ggtitle('Reads mapped to genome')
ggsave('results/reads_mapped_count.pdf',width=7,height=7)

ggplot(feces,aes(Sample.Endogenousity,Perc.mapped)) +
	geom_point() +
	geom_smooth(method=lm,se=FALSE) +
	xlab('Fraction host DNA in fecal extraction') +
	ylab('Fraction reads mapping to genome') +
	ggtitle('Enrichment rate')
ggsave('results/enrichment_rate.pdf',width=7,height=7)

ggplot(feces,aes(Sample.Endogenousity,Perc.mapped)) +
	geom_point() +
	geom_smooth(method=lm,se=FALSE) +
	xlim(c(0,1)) +
	ylim(c(0,1)) +
	xlab('Fraction host DNA in fecal extraction') +
	ylab('Fraction reads mapping to genome') +
	ggtitle('Enrichment rate')
ggsave('results/enrichment_rate_scaled.pdf',width=7,height=7)

