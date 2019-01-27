#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# Generate individual_info.csv from the file on Google Drive
# Sequenced Kafue Animals
# https://docs.google.com/spreadsheets/d/19jdEEj9rRTSlXLAzoiAlZ0szq0LGEvA6pVb5FU8D-Q8

library(ggplot2)

ind_info = read.csv("data/individual_info.csv")

pca = read.table('results/all.cleaned.LDpruned.eigenvec', header = FALSE)

names(pca) = c('FID','IID',paste0('PC',1:(ncol(pca)-2)))

pca$FID = gsub('\\.PE','',pca$FID)
pca$IID = gsub('\\.PE','',pca$IID)

group  = as.character(pca$FID)
taxon  = as.character(pca$FID)
ind_id = as.character(pca$FID)

for (i in 1:dim(ind_info)[1]) {

	this.ind = ind_info$Individual.ID[i]
	this.ind.info = ind_info[ind_info$Individual.ID == this.ind,]

	group  = replace(group,  group==this.ind, this.ind.info$Group)
	taxon  = replace(taxon,  taxon==this.ind, this.ind.info$Taxon)
	ind_id = replace(ind_id, ind_id==this.ind, gsub('_','-',this.ind.info$Individual.ID))

}

pca$group = group
pca$taxon = factor(taxon)
pca$ind_id = ind_id

# Trim outliers
sd.threshold = 3     # keep records within 3 standard deviations

key = c("Chunga HQ", "Mwengwa Rp", "Namiyezhi R", "Malala Cmp", "Musa Br", "Top Musa", "N Nkala Rd", "Ngoma AS", "Dendro P", "Nanzhila Pl", "Choma", "L Zambezi NP")
names(key) = LETTERS[1:12]
key = key[names(key) %in% pca$group]

pca = pca[pca$PC1 > mean(pca$PC1) - sd.threshold * sd(pca$PC1) & pca$PC1 < mean(pca$PC1) + sd.threshold * sd(pca$PC1),]
pca = pca[pca$PC2 > mean(pca$PC2) - sd.threshold * sd(pca$PC2) & pca$PC2 < mean(pca$PC2) + sd.threshold * sd(pca$PC2),]

ggplot(pca,aes(PC1,PC2,shape=group,color=taxon)) +
	geom_point(cex=3) +
	scale_shape_manual(name='Group',values=sort(unique(pca$group)),labels=key) +
	scale_color_manual(name='Taxon',values=c('#e41a1c','#4daf4a','#377eb8')) +
	theme_bw() + theme(legend.position='bottom',panel.grid=element_blank()) +
	xlab('Eigenvector 1') + ylab('Eigenvector 2')
ggsave(filename='results/pca_plot.pdf',width=6,height=6,useDingbats=FALSE)
