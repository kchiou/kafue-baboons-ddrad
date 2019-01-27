#!/usr/bin/env Rscript

library(ggplot2)

options(stringsAsFactors = FALSE)

mds = read.table('results/mds.txt', header=TRUE)
adm = read.table('results/all.cleaned.LDpruned.2.Q')

het.hyb = read.table('results/hybrids.het', header=TRUE)
het.kin = read.table('results/pure.kinda.het', header=TRUE)
het.cha = read.table('results/pure.chacma.het', header=TRUE)

# Swap ADM1 and ADM2 if they are not the order I arbitrarily like them in.
if (adm[which(mds$ind_id == 'BZ11-001'),]$V1 > 0.5) {
	adm = data.frame(V1=adm$V2, V2=adm$V1)
}

names(adm) = paste0('ADM', 1:2)
mds.adm = cbind(mds, adm)

# Combine pures and hybrids
het = rbind(het.kin,het.cha,het.hyb)

het$INDV = gsub('\\.PE','',het$INDV)

adm.het = merge(mds.adm,het,by.x='FID',by.y='INDV')

# Arbitrarily order by ind_id
adm.het = adm.het[order(adm.het$ind_id),]

ggplot(adm.het,aes(ADM1,F,col=group)) + geom_point() + xlab('ADMIXTURE Inferred Ancestry') + ylab('F') + scale_color_discrete(name='Group')
ggsave('results/inbreeding_by_ancestry.pdf',width=7,height=7)