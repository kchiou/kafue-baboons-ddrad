#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

library(ggplot2)

# Generate individual_info.csv from the file on Google Drive
# Sequenced Kafue Animals
# https://docs.google.com/spreadsheets/d/19jdEEj9rRTSlXLAzoiAlZ0szq0LGEvA6pVb5FU8D-Q8

ind_info = read.csv("data/individual_info.csv")

mds = read.table('results/all.cleaned.LDpruned.mds', header = TRUE)

mds$FID = gsub('\\.PE','',mds$FID)
mds$IID = gsub('\\.PE','',mds$IID)

group  = as.character(mds$FID)
taxon  = as.character(mds$FID)
ind_id = as.character(mds$FID)

for (i in 1:dim(ind_info)[1]) {

	this.ind = ind_info$Individual.ID[i]
	this.ind.info = ind_info[ind_info$Individual.ID == this.ind,]

	group  = replace(group,  group==this.ind, this.ind.info$Group)
	taxon  = replace(taxon,  taxon==this.ind, this.ind.info$Taxon)
	ind_id = replace(ind_id, ind_id==this.ind, gsub('_','-',this.ind.info$Individual.ID))

}

mds$group = group
mds$taxon = factor(taxon)
mds$ind_id = ind_id

key = c("Chunga HQ", "Mwengwa Rp", "Namiyezhi R", "Malala Cmp", "Musa Br", "Top Musa", "N Nkala Rd", "Ngoma Airstr", "Dendro P", "Nanzhila Pl", "Choma", "L Zambezi NP")
names(key) = LETTERS[1:12]
key = key[names(key) %in% mds$group]

ggplot(mds,aes(C1,C2,shape=group,color=taxon)) +
	geom_point(cex=3) +
	scale_shape_manual(name='Group',values=sort(unique(mds$group)),labels=key) +
	scale_color_manual(name='Taxon',values=c('#e41a1c','#4daf4a','#377eb8')) +
	theme_bw() + theme(legend.position='bottom',panel.grid=element_blank()) +
	xlab('Dimension 1') + ylab('Dimension 2') + ggtitle('IBS MDS plot')
ggsave(filename='results/ibs_mds_plot.pdf',width=7,height=7)


# Export mds now that info on individuals has been added.
write.table(mds, file="results/mds.txt")